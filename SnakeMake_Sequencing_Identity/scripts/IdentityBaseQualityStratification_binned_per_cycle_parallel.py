#!/usr/bin/env python3

import os
import sys
import pysam
import re
import gzip
import multiprocessing as mp
from collections import defaultdict
import pandas as pd
import numpy as np

############################################
# FUNZIONE PER CARICARE POSIZIONI VARIANTI DAL VCF
############################################
def load_variant_positions(vcf_path):
    print("Inizio caricamento varianti dal VCF:", vcf_path)
    snv_dict = defaultdict(set)
    indel_dict = defaultdict(set)

    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            chrom = rec.contig
            pos = rec.pos - 1  # 0-based
            ref = rec.ref
            for alt in rec.alts:
                if len(ref) == 1 and len(alt) == 1:
                    snv_dict[chrom].add(pos)
                elif len(ref) > len(alt):
                    del_len = len(ref) - len(alt)
                    indel_dict[chrom].add((pos, pos + del_len, 'D'))
                elif len(ref) < len(alt):
                    ins_len = len(alt) - len(ref)
                    indel_dict[chrom].add((pos, pos + ins_len, 'I'))
    print("Fine caricamento varianti.")
    return snv_dict, indel_dict

############################################
# FUNZIONE PER OTTENERE POSIZIONI MISMATCH DA TAG MD
############################################
def get_mismatch_positions(read):
    md = read.get_tag('MD') if read.has_tag('MD') else ''
    if not md:
        return set()

    mismatches = set()
    pos = read.reference_start  # 0-based

    for match in re.finditer(r"(\d+)|([A-Z]|\^[A-Z]+)", md):
        if match.group(1):
            pos += int(match.group(1))
        elif match.group(2).startswith('^'):
            pos += len(match.group(2)) - 1
        else:
            mismatches.add(pos)
            pos += 1

    return mismatches


############################################
# FUNZIONE PER CALCOLO CICLO DI SEQUENZIAMENTO
############################################
def get_cycle(read, query_pos):
    if query_pos is None:
        return None

    H_start = 0
    if read.cigartuples:
        if read.is_reverse:
            if read.cigartuples[-1][0] == 5:
                H_start = read.cigartuples[-1][1]
        else:
            if read.cigartuples[0][0] == 5:
                H_start = read.cigartuples[0][1]

    read_len = len(read.query_sequence)
    if read.is_reverse:
        cycle = (read_len + H_start) - query_pos
    else:
        cycle = query_pos + 1 + H_start
    return cycle

############################################
# DEFINIZIONE STRUTTURA DATI DINAMICA E FUNZIONE DI CONVERSIONE
############################################
def per_quality_dict():
    """Inizializza e ritorna la struttura dati dinamica per qualità e ciclo."""
    per_quality = defaultdict(lambda: defaultdict(lambda: {
        'match': 0,
        'match_correct': 0,
        'match_error': 0,
        'mismatch': 0,
        'mismatch_correct': 0,
        'mismatch_error': 0,
        'insertion': 0,
        'insertion_correct': 0,
        'insertion_error': 0
    }))
    return per_quality

def dictify(d):
    """
    Converte ricorsivamente un defaultdict annidato in un dict normale.
    Utile per rendere serializzabili i risultati per multiprocessing.
    """
    if isinstance(d, defaultdict):
        d = {k: dictify(v) for k, v in d.items()}
    return d

############################################
# FUNZIONE PER PROCESSARE UN SINGOLO CROMOSOMA
############################################
def process_chromosome(args):
    bam_path, chrom, snv_dict, indel_dict, mode = args
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Inizializza struttura dati
    per_quality = per_quality_dict()

    processed = 0
    bad_cycle_count = 0

    for read in bam.fetch(chrom):
        if mode == "read1" and read.is_read2:
            continue
        if mode == "read2" and read.is_read1:
            continue
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        processed += 1
        mismatch_positions = get_mismatch_positions(read)
        last_ref_pos = None

        for query_pos, ref_pos, cigar_op in read.get_aligned_pairs(matches_only=False, with_cigar=True):
            cycle = get_cycle(read, query_pos)
            if cycle is None:
                continue
            if cycle < 1 or cycle > len(read.query_sequence):
                bad_cycle_count += 1
                continue

            # Delezione (ignorata)
            if query_pos is None and ref_pos is not None:
                last_ref_pos = ref_pos
                continue

            # Inserzione
            elif query_pos is not None and ref_pos is None and cigar_op == 1:
                q = read.query_qualities[query_pos]
                per_quality[q][cycle]['insertion'] += 1
                insertion_site = last_ref_pos
                if insertion_site is not None:
                    match_found = any(start <= insertion_site <= end and t=='I' for start,end,t in indel_dict[chrom])
                else:
                    match_found = False
                if match_found:
                    per_quality[q][cycle]['insertion_correct'] += 1
                else:
                    per_quality[q][cycle]['insertion_error'] += 1

            # Match / Mismatch
            elif query_pos is not None and ref_pos is not None:
                q = read.query_qualities[query_pos]
                if ref_pos in mismatch_positions:
                    per_quality[q][cycle]['mismatch'] += 1
                    if ref_pos in snv_dict[chrom]:
                        per_quality[q][cycle]['mismatch_correct'] += 1
                    else:
                        per_quality[q][cycle]['mismatch_error'] += 1
                else:
                    per_quality[q][cycle]['match'] += 1
                    if ref_pos in snv_dict[chrom]:
                        per_quality[q][cycle]['match_error'] += 1
                    else:
                        per_quality[q][cycle]['match_correct'] += 1
                last_ref_pos = ref_pos

    bam.close()
    per_quality_dictified = dictify(per_quality)
    print(f"✔ Cromosoma {chrom} completato. Ignorati {bad_cycle_count} cicli fuori range.")
    return per_quality_dictified

############################################
# FUNZIONE PER FONDO DI DIZIONARI
############################################
def merge_per_quality_dicts(dicts):
    merged = per_quality_dict()
    for d in dicts:
        for q, cycles in d.items():
            for cycle, vals in cycles.items():
                for key in vals:
                    merged[q][cycle][key] += vals[key]
    return merged

############################################
# FUNZIONE PER SALVARE RISULTATI
############################################
def save_results(per_quality, output_tot, output_binned):
    records = []
    for q, cycles in per_quality.items():
        for cycle, vals in cycles.items():
            total = vals['match'] + vals['mismatch']
            total_all = total + vals['insertion']
            records.append({
                'BaseQuality': q,
                'Cycle': cycle,
                **vals,
                'identity': 1 - (vals['mismatch']/total) if total>0 else np.nan,
                'identity_filtered': 1 - (vals['mismatch_error']/total) if total>0 else np.nan,
                'identity_with_ins': 1 - ((vals['mismatch']+vals['insertion'])/total_all) if total_all>0 else np.nan,
                'identity_with_ins_filtered': 1 - ((vals['mismatch_error']+vals['insertion_error'])/total_all) if total_all>0 else np.nan
            })
    df = pd.DataFrame(records).sort_values(by=['BaseQuality','Cycle'])
    df.to_csv(output_tot, sep="\t", index=False, float_format="%.6f", na_rep="NaN")
    extract_binned_identity_by_cycle_and_quality(df, output_binned)

############################################
# FUNZIONE PER CALCOLO IDENTITÀ BINNED PER CICLO E QUALITÀ
############################################
def extract_binned_identity_by_cycle_and_quality(df, output_binned):
    bins = [0,20,30,51]
    labels = ['0-19','20-29','30+']
    df['BQ_bin'] = pd.cut(df['BaseQuality'], bins=bins, labels=labels, right=False)
    grouped = df.groupby(['BQ_bin','Cycle'], observed=False).agg({
        'match':'sum','match_correct':'sum','match_error':'sum',
        'mismatch':'sum','mismatch_correct':'sum','mismatch_error':'sum',
        'insertion':'sum','insertion_correct':'sum','insertion_error':'sum'
    }).reset_index()
    total = grouped['match']+grouped['mismatch']
    total_all = total+grouped['insertion']
    grouped['identity'] = 1 - (grouped['mismatch']/total)
    grouped['identity_filtered'] = 1 - (grouped['mismatch_error']/total)
    grouped['identity_with_ins'] = 1 - ((grouped['mismatch']+grouped['insertion'])/total_all)
    grouped['identity_with_ins_filtered'] = 1 - ((grouped['mismatch_error']+grouped['insertion_error'])/total_all)
    grouped = grouped.sort_values(by=['BQ_bin','Cycle'])
    grouped.to_csv(output_binned, sep="\t", index=False, float_format="%.6f", na_rep="NaN")

############################################
# FUNZIONE PER LANCIARE LA PARALLELIZZAZIONE
############################################
def calculate_identity_parallel(bam_path, vcf_path, output_tot, output_binned, mode="all", threads=4):
    snv_dict, indel_dict = load_variant_positions(vcf_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    chroms = list(bam.references)
    bam.close()

    args_list = [(bam_path, chrom, snv_dict, indel_dict, mode) for chrom in chroms]

    with mp.Pool(threads) as pool:
        results = pool.map(process_chromosome, args_list)

    merged = merge_per_quality_dicts(results)
    save_results(merged, output_tot, output_binned)

############################################
# MAIN
############################################
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Uso: python script.py <input.bam> <input.vcf> <output_prefix> [mode] [threads]")
        sys.exit(1)

    bam = sys.argv[1]
    vcf = sys.argv[2]
    output_tot = sys.argv[3]
    output_binned = sys.argv[4]
    mode = sys.argv[5] if len(sys.argv) > 5 else "all"
    threads = int(sys.argv[6]) if len(sys.argv) > 6 else 1

    if not os.path.exists(bam) or not os.path.exists(vcf):
        print("Errore: file BAM o VCF non trovato.")
        sys.exit(1)

    print(f"File BAM: {bam}")
    print(f"File VCF: {vcf}")
    print(f"Output prefix: {output_tot}, {output_binned}")
    print(f"Modalità: {mode}")
    print(f"Threads: {threads}")

    calculate_identity_parallel(bam, vcf, output_tot, output_binned, mode, threads)
    print("✔ FINE calcolo identità per qualità e ciclo")