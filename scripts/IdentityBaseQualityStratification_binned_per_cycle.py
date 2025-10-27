# Script: Calcolo dell'identità di sequenziamento stratificata per la qualità della base e per ciclo di sequenziamento. Le qualità delle basi sono binnate in quattro diversi intervalli: [0-2], [3-17], [18-29], [30+].

# Questo script analizza un file BAM/SAM e un file VCF di varianti note per
# calcolare l'identità delle basi rispetto al genoma di riferimento.
# L'identità viene stratificata per livello di qualità delle basi (Phred score) assegnato dal sequenziatore e per ciclo di sequenziamento (numero della base all'interno della read).
# Lo script ignora le delezioni in quanto quest'ultime non presentano basi sequenziate e la rispettiva qualità.
# I match vengono considerati sempre corretti, sia se ricadono in varianti che se ricadono fuori da varianti (questo per via delle varianti in eterozigosi).


# Input
# -----
#   - BAM/SAM (allineamento con tag MD)
#   - VCF     (varianti note, SNV e indel)
#
# Output
# ------
#   - <output>.gz       : tabella dei conteggi per qualità
#   - <output>.gz.stats : tabella con statistiche di identità

# Uso
# ---
#   python script.py <input.bam> <input.vcf> <output.gz>

#esempio di esecuzione:
# nohup python -u ../../../scripts/IdentityBaseQualityStratification_binned_per_cycle.py NA12878_KAPA_HCR.bam ../../../HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz IDENTITY_NA12878_KAPA_full_basequal_strat > IDENTITY_NA12878_KAPA_full_basequal_strat_cycle.log &

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import sys
import pysam
import re
from collections import defaultdict
from tqdm import tqdm

############################################
# FUNZIONE PER CARICARE POSIZIONI VARIANTI DAL VCF
############################################

def load_variant_positions(vcf_path):
    """Carica le posizioni delle varianti SNV e INDEL dal file VCF."""
    print("Inizio caricamento varianti dal VCF:", vcf_path)
    snv_dict = defaultdict(set)
    indel_dict = defaultdict(set)

    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            chrom = rec.contig
            pos = rec.pos  # 1-based
            ref = rec.ref
            alts = rec.alts

            for alt in alts:
                if len(ref) == 1 and len(alt) == 1:
                    snv_dict[chrom].add(pos - 1)  # 0-based
                elif len(ref) > len(alt):  # Deletion
                    del_len = len(ref) - len(alt)
                    del_start = pos - 1
                    del_end = pos - 1 + del_len
                    indel_dict[chrom].add((del_start, del_end, 'D'))
                elif len(ref) < len(alt):  # Insertion
                    ins_len = len(alt) - len(ref)
                    ins_start = pos - 1
                    ins_end = pos - 1 + ins_len
                    indel_dict[chrom].add((ins_start, ins_end, 'I'))
    print("Fine caricamento varianti.")
    return snv_dict, indel_dict

############################################
# FUNZIONE PER OTTENERE POSIZIONI MISMATCH DA TAG MD NELLA READ
############################################
def get_mismatch_positions(read):
    """Ritorna le posizioni genomiche dei mismatch da un record BAM usando il tag MD."""
    md = read.get_tag('MD') if read.has_tag('MD') else ''
    if not md:
        return []
    mismatch_positions = []
    pos = read.reference_start
    for match in re.finditer(r"(\d+)|([A-Z]|\^[A-Z]+)", md):
        if match.group(1):
            pos += int(match.group(1))
        elif match.group(2):
            if match.group(2).startswith('^'):
                pos += len(match.group(2)) - 1
            else:
                mismatch_positions.append(pos)
                pos += 1
    return mismatch_positions

############################################
# FUNZIONE PER CALCOLO CICLO REALE (calcola il ciclo di sequenziamento a partire dalla posizione della base nella read, tenendo conto di soft/hard clip iniziali e dell'orientamento della read)
############################################
# def get_cycle(read, query_pos):
#     """Restituisce il ciclo di sequenziamento corretto per l'orientamento della read."""
#     return len(read.query_sequence) - query_pos if read.is_reverse else query_pos + 1

def get_cycle(read, query_pos):
    """
    Restituisce il ciclo di sequenziamento reale (1-based),
    considerando le hardclip iniziali nel senso del SEQUENZIAMENTO.
    Le softclip sono già incluse in query_sequence, quindi non conteggiate.
    """
    if query_pos is None:
        return None

    H_start = 0
    if read.cigartuples:
        if read.is_reverse:
            # Hardclip all'inizio del sequenziamento = alla FINE della CIGAR
            if read.cigartuples[-1][0] == 5:
                H_start = read.cigartuples[-1][1]
        else:
            # Hardclip all'inizio del sequenziamento = all'INIZIO della CIGAR
            if read.cigartuples[0][0] == 5:
                H_start = read.cigartuples[0][1]

    read_len = len(read.query_sequence)

    if read.is_reverse:
        # Le reverse vengono lette dal fondo → inverto il ciclo
        cycle = (read_len + H_start) - query_pos
    else:
        cycle = query_pos + 1 + H_start
    # print(f"[DEBUG] Read: {read.query_name}, cigar {read.cigarstring}, is_reverse: {read.is_reverse}, query_pos: {query_pos}, H_start: {H_start}, read_len: {read_len}, cycle: {cycle}")
    return cycle

############################################
# FUNZIONE PRINCIPALE
############################################

def calculate_identity_by_cycle_and_quality(bam_path, vcf_path, output_prefix):
    snv_dict, indel_dict = load_variant_positions(vcf_path)
    bam = pysam.AlignmentFile(bam_path, "rb")

    print()
    print("Inizio analisi BAM:", bam_path)

    
    # Conta il numero totale di reads nel BAM (richiede .bai)
    try:
        print("Conta il numero totale di reads nel BAM...")
        total_reads = bam.count()
        print(total_reads, " reads totali nel BAM.")
    except ValueError:
        print("Attenzione: BAM non indicizzato, impossibile stimare il numero totale di reads.")
        total_reads = None
    bam.close()
    
    # Inizializza barra di progresso
    if total_reads:
        pbar = tqdm(total=total_reads, unit=" reads", ncols=100, file=sys.stdout)
    else:
        pbar = tqdm(unit=" reads", ncols=100, file=sys.stdout)
    update_every = 1000000  # aggiorna la barra ogni 1M di reads
    processed = 0

    # Struttura dati: quality → cycle → conteggi
    per_quality = defaultdict(lambda: defaultdict(lambda: {
        'match': 0, 'match_correct': 0, 'match_error': 0,
        'mismatch': 0, 'mismatch_correct': 0, 'mismatch_error': 0,
        'insertion': 0, 'insertion_correct': 0, 'insertion_error': 0
    }))
    
    bam = pysam.AlignmentFile(bam_path, "rb")
    bad_cycle_count = 0
    warned = False

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        processed += 1
        if processed % update_every == 0:
            pbar.update(update_every)
            print(f"[DEBUG] Reads totali viste: {processed}")

        chrom = bam.get_reference_name(read.reference_id)
        mismatch_positions = set(get_mismatch_positions(read))
        read_len = len(read.query_sequence)

        last_ref_pos = None
        for query_pos, ref_pos, cigar_op in read.get_aligned_pairs(matches_only=False, with_cigar=True):

            # Calcolo del ciclo di sequenziamento
            cycle = get_cycle(read, query_pos)
            if cycle is None:
                continue
            if cycle < 1 or cycle > read_len:
                bad_cycle_count += 1
                print(f"[ATTENZIONE] Rilevati {bad_cycle_count} cicli fuori range in {bam_path}, nella read {read.query_name}. Verranno ignorati. CIGAR: {read.cigarstring}")
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
                    match_found = any(
                        (start <= insertion_site <= end and typ == 'I')
                        for (start, end, typ) in indel_dict[chrom]
                    )
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

    # Aggiorna eventuali ultime read non multiple di 1M
    if processed % update_every != 0:
        pbar.update(processed % update_every)
    bam.close()
    print("\nAnalisi BAM completata. Processate", processed, "reads.")
    print(f"Ignorate {bad_cycle_count} basi con ciclo fuori range (> read length).")
    print("Inizio calcolo dell'identità e salvataggio dei risultati...")

    # Salva file per qualità e ciclo
    records = []
    for q, cycles in per_quality.items():
        for cycle, vals in cycles.items():
            total = vals['match'] + vals['mismatch']
            total_all = vals['match'] + vals['mismatch'] + vals['insertion']
            identity = 1 - (vals['mismatch'] / total) if total > 0 else np.nan
            identity_filtered = 1 - (vals['mismatch_error'] / total) if total > 0 else np.nan
            identity_with_ins = 1 - ((vals['mismatch'] + vals['insertion']) / total_all) if total_all > 0 else np.nan
            identity_with_ins_filtered = 1 - ((vals['mismatch_error'] + vals['insertion_error']) / total_all) if total_all > 0 else np.nan

            records.append({
                'BaseQuality': q,
                'Cycle': cycle,
                **vals,
                'identity': identity,
                'identity_filtered': identity_filtered,
                'identity_with_ins': identity_with_ins,
                'identity_with_ins_filtered': identity_with_ins_filtered
            })

    df = pd.DataFrame(records).sort_values(by=['BaseQuality','Cycle'])
    df.to_csv(output_prefix + ".by_quality_cycle.tsv", sep="\t", index=False, float_format="%.6f", na_rep="NaN")

    # Calcola anche la versione per bin
    extract_binned_identity_by_cycle_and_quality(df, output_prefix)


############################################
# FUNZIONE PER BIN DI QUALITÀ
############################################

def extract_binned_identity_by_cycle_and_quality(df, output_prefix):
    bins = [0, 3, 18, 30, 51]
    labels = ['0-2', '3-17', '18-29', '30+']
    df['BQ_bin'] = pd.cut(df['BaseQuality'], bins=bins, labels=labels, right=False)

    grouped = df.groupby(['BQ_bin','Cycle'], observed=False).agg({
        'match': 'sum',
        'match_correct': 'sum',
        'match_error': 'sum',
        'mismatch': 'sum',
        'mismatch_correct': 'sum',
        'mismatch_error': 'sum',
        'insertion': 'sum',
        'insertion_correct': 'sum',
        'insertion_error': 'sum'
    }).reset_index()

    total = grouped['match'] + grouped['mismatch']
    total_all = grouped['match'] + grouped['mismatch'] + grouped['insertion']
    grouped['identity'] = 1 - (grouped['mismatch'] / total)
    grouped['identity_filtered'] = 1 - (grouped['mismatch_error'] / total)
    grouped['identity_with_ins'] = 1 - ((grouped['mismatch'] + grouped['insertion']) / total_all)
    grouped['identity_with_ins_filtered'] = 1 - ((grouped['mismatch_error'] + grouped['insertion_error']) / total_all)

    grouped = grouped.sort_values(by=['BQ_bin','Cycle'])
    grouped.to_csv(output_prefix + ".by_bin_cycle.tsv", sep="\t", index=False, float_format="%.6f", na_rep="NaN")


############################################
# MAIN
############################################

if __name__ == "__main__":
    print("INIZIO calcolo identità stratificata per qualità e ciclo")

    bam = sys.argv[1]
    vcf = sys.argv[2]
    output = sys.argv[3]

    if not os.path.exists(bam):
        print("Errore: file BAM non trovato:", bam)
        sys.exit(1)
    if not os.path.exists(vcf):
        print("Errore: file VCF non trovato:", vcf)
        sys.exit(1)

    calculate_identity_by_cycle_and_quality(bam, vcf, output)

    print("FINE calcolo identità per qualità e ciclo")
    print("File prodotti:")
    print(" -", output + ".by_quality_cycle.tsv")
    print(" -", output + ".by_bin_cycle.tsv")