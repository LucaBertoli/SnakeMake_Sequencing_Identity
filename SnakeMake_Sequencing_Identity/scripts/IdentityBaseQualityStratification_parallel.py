# Script: Calcolo dell'identità di sequenziamento stratificata per la qualità della base

# Questo script analizza un file BAM/SAM e un file VCF di varianti note per
# calcolare l'identità delle basi rispetto al genoma di riferimento.
# L'identità viene stratificata per livello di qualità delle basi (Phred score) assegnato dal sequenziatore.
# Lo script ignora le delezioni in quanto quest'ultime non presentano basi sequenziate e la rispettiva qualità.
# I match vengono considerati sempre corretti, sia se ricadono in varianti che se ricadono fuori da varianti (questo per via delle varianti in eterozigosi).

# Lo script:
#   1. Carica le posizioni delle varianti dal VCF (SNV e indel). Usare gold set se possibile.
#   2. Per ciascuna read del BAM:
#      - Identifica match, mismatch e inserzioni utilizzando i tag MD e CIGAR.
#      - Conta quanti eventi coincidono con varianti note (corretti) o non note (errori).
#      - Registra i conteggi in base alla qualità delle basi.
#   3. Produce una tabella tab-delimited (compressa in gzip) con i conteggi
#      per ogni livello di qualità da 0 a 50.
#   4. Calcola le metriche di identità:
#        - identity:              1 - (mismatch / (match + mismatch))
#        - identity_filtered:     1 - (mismatch_error / (match + mismatch))
#        - identity_with_ins:     1 - ((mismatch + insertion) / (match + mismatch + insertion))
#        - identity_with_ins_filtered:  1 - ((mismatch_error + insertion_error) / (match + mismatch + insertion))

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


import pysam
import pandas as pd
import os
import sys
import gzip
import re
import multiprocessing as mp
from collections import defaultdict

############################################
# FUNZIONE PER CARICARE POSIZIONI VARIANTI DAL VCF
############################################
def load_variant_positions(vcf_path):
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
                    del_start = pos
                    del_end = pos + del_len
                    indel_dict[chrom].add((del_start, del_end, 'D'))
                elif len(ref) < len(alt):
                    ins_len = len(alt) - len(ref)
                    ins_start = pos
                    ins_end = pos + ins_len
                    indel_dict[chrom].add((ins_start, ins_end, 'I'))
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
# FUNZIONE PER PROCESSARE UN SINGOLO CROMOSOMA
############################################
def process_chromosome(args):
    bam_path, chrom, snv_dict, indel_dict = args
    bam = pysam.AlignmentFile(bam_path, "rb")

    # base_quality_dictionary locale
    bq_dict = {str(i): {'match':0,'match_correct':0,'match_error':0,
                        'mismatch':0,'mismatch_correct':0,'mismatch_error':0,
                        'insertion':0,'insertion_correct':0,'insertion_error':0}
               for i in range(51)}

    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 60:
            continue

        mismatch_positions = get_mismatch_positions(read)
        last_ref_pos = None

        for query_pos, ref_pos, cigar_op in read.get_aligned_pairs(matches_only=False, with_cigar=True):
            if query_pos is None and ref_pos is not None:
                last_ref_pos = ref_pos
                continue

            elif query_pos is not None and ref_pos is None and cigar_op == 1:
                q = read.query_qualities[query_pos]
                bq_dict[str(q)]['insertion'] += 1
                insertion_site = last_ref_pos
                if insertion_site is not None and any((start <= insertion_site <= end and t=='I') for start,end,t in indel_dict[read.reference_name]):
                    bq_dict[str(q)]['insertion_correct'] += 1
                else:
                    bq_dict[str(q)]['insertion_error'] += 1

            elif query_pos is not None and ref_pos is not None:
                q = read.query_qualities[query_pos]
                if ref_pos in mismatch_positions:
                    bq_dict[str(q)]['mismatch'] += 1
                    if ref_pos in snv_dict[read.reference_name]:
                        bq_dict[str(q)]['mismatch_correct'] += 1
                    else:
                        bq_dict[str(q)]['mismatch_error'] += 1
                else:
                    bq_dict[str(q)]['match'] += 1
                    if ref_pos in snv_dict[read.reference_name]:
                        bq_dict[str(q)]['match_error'] += 1
                    else:
                        bq_dict[str(q)]['match_correct'] += 1
                last_ref_pos = ref_pos

    bam.close()
    print(f"[PID {os.getpid()}] ✔ Finished chromosome {chrom}", flush=True)
    return bq_dict

############################################
# FUNZIONE PER UNIRE I DIZIONARI DEI PROCESSI PARALLELI
############################################
def merge_bq_dicts(dicts):
    merged = {str(i): {'match':0,'match_correct':0,'match_error':0,
                       'mismatch':0,'mismatch_correct':0,'mismatch_error':0,
                       'insertion':0,'insertion_correct':0,'insertion_error':0}
              for i in range(51)}
    for d in dicts:
        for k in merged:
            for field in merged[k]:
                merged[k][field] += d[k][field]
    return merged

############################################
# FUNZIONE PER LANCIARE LA PARALLELIZZAZIONE E CALCOLARE L'IDENTITÀ STRATIFICATA
############################################
def calculate_stratified_identity_parallel(bam_path, vcf_path, output_path, threads=4):
    snv_dict, indel_dict = load_variant_positions(vcf_path)

    bam = pysam.AlignmentFile(bam_path, "rb")
    chroms = list(bam.references)
    bam.close()

    args = [(bam_path, c, snv_dict, indel_dict) for c in chroms]

    with mp.Pool(threads) as pool:
        results = pool.map(process_chromosome, args)

    merged = merge_bq_dicts(results)

    with gzip.open(output_path, 'wt') as out:
        out.write("BaseQuality\tMatch\tMatch_correct\tMatch_error\tMismatch\tMismatch_correct\tMismatch_error\tInsertion\tInsertion_correct\tInsertion_error\tIdentity\tIdentity_filtered\tIdentity_with_ins\tIdentity_with_ins_filtered\n")
        for q in sorted(merged.keys(), key=int):
            vals = merged[q]
            vals['identity'] = 1 - (vals['mismatch'] / (vals['match'] + vals['mismatch'])) if (vals['match'] + vals['mismatch']) > 0 else float('nan')
            vals['identity_filtered'] = 1 - ((vals['mismatch_error']) / (vals['match'] + vals['mismatch'])) if (vals['match'] + vals['mismatch']) > 0 else float('nan')
            vals['identity_with_ins'] = 1 - ((vals['mismatch'] + vals['insertion']) / (vals['match'] + vals['mismatch'] + vals['insertion'])) if (vals['match'] + vals['mismatch'] + vals['insertion']) > 0 else float('nan')
            vals['identity_with_ins_filtered'] = 1 - ((vals['mismatch_error'] + vals["insertion_error"]) / (vals['match'] + vals['mismatch'] + vals['insertion'])) if (vals['match'] + vals['mismatch'] + vals['insertion']) > 0 else float('nan')
            out.write(f"{q}\t{vals['match']}\t{vals['match_correct']}\t{vals['match_error']}\t{vals['mismatch']}\t{vals['mismatch_correct']}\t{vals['mismatch_error']}\t{vals['insertion']}\t{vals['insertion_correct']}\t{vals['insertion_error']}\t{vals['identity']:.6f}\t{vals['identity_filtered']:.6f}\t{vals['identity_with_ins']:.6f}\t{vals['identity_with_ins_filtered']:.6f}\n")

    print(f"✔ Done → Saved: {output_path}")

############################################
# FUNZIONE PER STRATIFICARE LE STATISTICHE DI IDENTITÀ IN BIN DI QUALITÀ
############################################
def extract_binned_identity_statistics(output_tot, output_binned):
    df = pd.read_csv(output_tot, sep="\t")

    bins = [0, 20, 30, 51]
    labels = ['0-19', '20-29', '30+']
    df['BQ_bin'] = pd.cut(df['BaseQuality'], bins=bins, labels=labels, right=False)

    binned_stats = df.groupby('BQ_bin', observed=False).agg({
        'Match': 'sum',
        'Match_correct': 'sum',
        'Match_error': 'sum',
        'Mismatch': 'sum',
        'Mismatch_correct': 'sum',
        'Mismatch_error': 'sum',
        'Insertion': 'sum',
        'Insertion_correct': 'sum',
        'Insertion_error': 'sum'
    }).reset_index()

    binned_stats['identity'] = 1 - (binned_stats['Mismatch'] / (binned_stats['Match'] + binned_stats['Mismatch']))
    binned_stats['identity_filtered'] = 1 - ((binned_stats['Mismatch_error']) / (binned_stats['Match'] + binned_stats['Mismatch']))
    binned_stats['identity_with_ins'] = 1 - ((binned_stats['Mismatch'] + binned_stats['Insertion']) / (binned_stats['Match'] + binned_stats['Mismatch'] + binned_stats['Insertion']))
    binned_stats['identity_with_ins_filtered'] = 1 - ((binned_stats['Mismatch_error'] + binned_stats["Insertion_error"]) / (binned_stats['Match'] + binned_stats['Mismatch'] + binned_stats['Insertion']))

    binned_stats.to_csv(output_binned, sep="\t", index=False, float_format="%.6f", na_rep="NaN")

############################################
# MAIN
############################################
if __name__ == "__main__":
    print("INIZIO del calcolo della identità stratificando i dati per qualità delle basi")

    bam=sys.argv[1]
    VCF=sys.argv[2]
    output_tot=sys.argv[3]
    output_binned=sys.argv[4]
    threads=int(sys.argv[5]) if len(sys.argv) > 5 else 1

    if not os.path.exists(bam):
        print("Il file pileup non esiste:", bam)
        sys.exit(1)

    if not os.path.exists(VCF):
        print("Il file VCF non esiste:", VCF)
        sys.exit(1)

    calculate_stratified_identity_parallel(bam, VCF, output_tot, threads)
    extract_binned_identity_statistics(output_tot, output_binned)