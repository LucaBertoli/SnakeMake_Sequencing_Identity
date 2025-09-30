#!/usr/bin/env python3

"""
Calcola l'identità delle reads allineate in un BAM, escludendo le posizioni con varianti note da VCF (SNV/INDEL).
Filtra le reads con più di N mismatch/indel e confronta:
- Phred stimato da identità
- Phred medio da sequencer

Genera output TSV e plot confronto distribuzioni Phred.

Uso:
    python IdentityRevelations.py input.bam input.vcf.gz output.tsv.gz nr_mismatches [all|read1|read2]

    esempio:
        python IdentityRevelations.py input.bam input.vcf.gz output.tsv.gz 5nohup python IdentityRevelations_select_mismatch.py 20250617_WGS_2_1/start_sorted_HCR.bam HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz IDENTITY_run_workflow_20250617_WGS_2_1_GIAB_more_equal_3_mismatch_plotting_2.tsv.gz 3 > IDENTITY_run_workflow_20250617_WGS_2_1_GIAB_more_equal_3_mismatch_plotting_2.log &
"""

import pysam
import re
import gzip
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

#catica le le posizioni delle varianti SNV e INDEL da un file VCF (da utilizzare per esempio il VCF di GIAB.)
def load_variant_positions(vcf_path):
    #definisce due dizionari per le posizioni delle varianti: uno per gli SNV e uno per gli INDEL.
    snv_dict = defaultdict(set)
    indel_dict = defaultdict(set)

    #apre il file VCF in modalità lettura.
    with pysam.VariantFile(vcf_path) as vcf:
        #scorre le varianti nel file VCF e salva le posizioni ed alleli
        for rec in vcf:
            chrom = rec.contig
            pos = rec.pos  # 1-based
            ref = rec.ref
            alts = rec.alts

            #per ogni allele della variante, controlla se è un SNV o un INDEL.
            # Se è un SNV, aggiunge la posizione al dizionario degli SNV.
            # Se è un INDEL, calcola la lunghezza e le posizioni di inserzione o delezione e le aggiunge al dizionario degli INDEL.
            for alt in alts:
                if len(ref) == 1 and len(alt) == 1:
                    snv_dict[chrom].add(pos)
                elif len(ref) > len(alt):  # Deletion
                    del_len = len(ref) - len(alt)
                    del_start = pos
                    del_end = pos + del_len
                    indel_dict[chrom].add((del_start, del_end, 'D'))
                elif len(ref) < len(alt):  # Insertion
                    ins_len = len(alt) - len(ref)
                    ins_start = pos
                    ins_end = pos + ins_len
                    indel_dict[chrom].add((ins_start, ins_end, 'I'))
    return snv_dict, indel_dict

# Estrae le posizioni dei mismatch da un read, basandosi sul tag MD.
# Il tag MD contiene informazioni sui mismatch rispetto al riferimento.
def get_mismatch_positions(read):
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

# Converte la percentuale di identità in un punteggio Phred.
# La percentuale di identità è un valore tra 0 e 1, dove 1.0 rappresenta il 100% di identità.
# Il punteggio Phred è calcolato come -10 * log10(1 - percentuale di identità).
# Se la percentuale di identità è maggiore di 0.999, viene limitata a 0.999 per evitare errori di logaritmo.
#la percentuale di identità a 0.999 equivale a Phred 30.
def percent_identity_to_phred(identity_percent):
    identity = min(float(identity_percent), 0.999)
    error_prob = 1.0 - identity
    if error_prob <= 0:
        return 60.0
    return round(-10 * np.log10(error_prob), 2)

#funzione centrale dello script.
def compute_identity(bam_path, vcf_path, output_path, mismatches, read_mode="all"):

    # Carica le posizioni delle varianti SNV e INDEL dal VCF
    snv_dict, indel_dict = load_variant_positions(vcf_path)
    # Apri il file BAM in modalità lettura
    bam = pysam.AlignmentFile(bam_path, "rb")

    #inizializza due file di output: uno per i dati di identità e uno per i dati Phred calcolati sulla qualità media delle basi delle read.
    with gzip.open(output_path, 'wt') as out, gzip.open(output_path + ".phred.tsv.gz", 'wt') as out_phred:
        out.write(
            "read\tidentity\tidentity_filtered\tidentity_with_indels\tidentity_filtered_with_indels\t"
            "aligned_length_query\taligned_length_total\tmismatches_total\tmismatches_filtered\t"
            "insertions\tins_bases\tins_filtered\tins_filtered_bases\tdeletions\t"
            "del_bases\tdel_filtered\tdel_filtered_bases\n"
        )

        #esclude le read non allineate, secondarie o supplementari.
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            #nel caso in cui si voglia analizzare solo le read di tipo read1 o read2, si filtra in base al flag is_read1 o is_read2.
            if read_mode == "read1" and not read.is_read1:
                continue
            if read_mode == "read2" and not read.is_read2:
                continue

            # Calcola la lunghezza allineata della query e se è zero skippa la read.
            aligned_length_query = read.query_alignment_length
            if aligned_length_query == 0:
                continue

            # Estrae il cromosoma di riferimento e le posizioni dei mismatch.
            chrom = bam.get_reference_name(read.reference_id)
            #per la read, estrae i mismatch utilizzando il tag MD.
            mismatch_positions = get_mismatch_positions(read)
            #filtra i mismatch della read, controllando se i mismatch nella read sono presenti nel dizionario delle varianti SNV (creato a partire dal VCF).
            mismatch_filtered = [p for p in mismatch_positions if (p + 1) not in snv_dict[chrom]]

            #Estrae la posizione di start della read ed inizializza liste per inserzioni e delezioni.
            ref_pos = read.reference_start
            ins_positions = []
            del_intervals = []
            inserted_bases = 0
            deleted_bases = 0

            #Scorre la CIGAR della read per estrarre le indel (non filtrate). Salva la lunghezza delle indel (non filtrate), appende gli eventi alle relative liste.
            for op, length in read.cigartuples:
                if op in (0, 7, 8):  # match or mismatch
                    ref_pos += length
                elif op == 1:  # insertion
                    ins_positions.append((ref_pos, ref_pos + length))
                    inserted_bases += length
                elif op == 2:  # deletion
                    del_intervals.append((ref_pos, ref_pos + length, 'D'))
                    deleted_bases += length
                    ref_pos += length

            aligned_length_total = aligned_length_query

            #filtra le indel trovate, confrontandole con le indel presenti nel dizionario delle varianti (creato a partire dal VCF).
            ins_filtered = [i for i in ins_positions if (i[0], i[1], 'I') not in indel_dict[chrom]]
            del_filtered = [d for d in del_intervals if (d[0], d[1], 'D') not in indel_dict[chrom]]

            # Calcola il numero totale di errori (mismatch + indel) e verifica se supera il numero massimo di mismatch consentiti.
            total_errors = (len(mismatch_filtered) + sum(i[1] - i[0] for i in ins_filtered) + sum(d[1] - d[0] for d in del_filtered))
            if total_errors > mismatches:
                #calcola le identità considerando i mismatch e le indel filtrate e non filtrate.
                identity = 1 - (len(mismatch_positions) / aligned_length_total)
                identity_filtered = 1 - (len(mismatch_filtered) / aligned_length_total)
                identity_with_indels = 1 - ((len(mismatch_positions) + inserted_bases + deleted_bases) / aligned_length_total)
                identity_filtered_with_indels = 1 - ((len(mismatch_filtered) +
                                                      sum(i[1] - i[0] for i in ins_filtered) +
                                                      sum(d[1] - d[0] for d in del_filtered)) / aligned_length_total)

                # Scrive i risultati nel file di output.
                out.write(
                    f"{read.query_name}\t"
                    f"{identity:.6f}\t"
                    f"{identity_filtered:.6f}\t"
                    f"{identity_with_indels:.6f}\t"
                    f"{identity_filtered_with_indels:.6f}\t"
                    f"{aligned_length_query}\t"
                    f"{aligned_length_total}\t"
                    f"{len(mismatch_positions)}\t{len(mismatch_filtered)}\t"
                    f"{len(ins_positions)}\t{inserted_bases}\t{len(ins_filtered)}\t{sum(i[1] - i[0] for i in ins_filtered)}\t"
                    f"{len(del_intervals)}\t{deleted_bases}\t{len(del_filtered)}\t{sum(d[1] - d[0] for d in del_filtered)}\n"
                )

                # Calcola il Phred stimato basato sull'identità.
                phred_identity_based = percent_identity_to_phred(identity_filtered_with_indels)
                #estrae le qualità delle basi della read
                quals = read.query_alignment_qualities
                #calcola il Phred medio basato sulla qualità delle basi della read e stampa i risultati nel file in output
                if quals:
                    phred_basequality_based = sum(quals) / len(quals)
                    out_phred.write(
                        f"{read.query_name}\t"
                        f"{phred_identity_based:.2f}\t"
                        f"{phred_basequality_based:.2f}\n"
                    )

    print(f"Output scritto in {output_path}")
    print("Elaborazione completata.")

    # --- Plotting ---
    #carica il file contenente le qualità (Phred) calcolate sulla base dell'identità e della qualità media delle basi.
    df_phred = pd.read_csv(output_path + ".phred.tsv.gz", sep='\t', names=['read', 'estimated', 'measured'])

    #definisce i bin dell'istogramma per i punteggi Phred.
    bins = np.arange(0, 61, 1)
    #genera i due istogrammi per i punteggi Phred stimati e misurati.
    counts_estimated, _ = np.histogram(df_phred['estimated'], bins=bins)
    counts_measured, _ = np.histogram(df_phred['measured'], bins=bins)

    # Crea un DataFrame per l'istogramma e salva in CSV
    df_hist = pd.DataFrame({
        'bin_start': bins[:-1],
        'bin_end': bins[1:],
        'count_estimated': counts_estimated,
        'count_measured': counts_measured
    })
    df_hist.to_csv(f'{output_path}.hist.csv', index=False)

    # Calcola i centri dei bin e raggruppa per bin
    df_hist['bin_center'] = (df_hist['bin_start'] + df_hist['bin_end']) / 2
    df_hist['bin_group'] = np.round(df_hist['bin_center']).astype(int)

    # Raggruppa per bin e somma i conteggi
    df_plot = df_hist.groupby('bin_group').agg({
        'count_estimated': 'sum',
        'count_measured': 'sum'
    }).reset_index().sort_values('bin_group')

    # Plotta l'istogramma sovrapposto
    plt.figure(figsize=(12, 6))
    width = 0.4
    plt.bar(df_plot['bin_group'] - width/2, df_plot['count_estimated'], width=width, label='Estimated by Alignment', color='orange')
    plt.bar(df_plot['bin_group'] + width/2, df_plot['count_measured'], width=width, label='Measured by Sequencer', color='steelblue')
    plt.xlabel('Phred Score')
    plt.ylabel('Read Count')
    plt.title(f'Distribution of Estimated vs Measured Phred Scores ({output_path})')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_path}.png', dpi=300)
    plt.close()

    #Plotta il subplot con i due istogrammi affiancati.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    axes[0].bar(df_plot['bin_group'], df_plot['count_estimated'], color='orange', width=0.8)
    axes[0].set_title('Estimated by Alignment')
    axes[0].set_xlabel('Phred Score')
    axes[0].set_ylabel('Read Count')
    axes[1].bar(df_plot['bin_group'], df_plot['count_measured'], color='steelblue', width=0.8)
    axes[1].set_title('Measured by Sequencer')
    axes[1].set_xlabel('Phred Score')
    plt.suptitle(f'Distribution of Phred Scores ({output_path})')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'{output_path}_subplot.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        print("Uso: python IdentityRevelations.py input.bam input.vcf.gz output.tsv.gz nr_mismatches [all|read1|read2]")
        sys.exit(1)

    bam_file = sys.argv[1]
    vcf_file = sys.argv[2]
    output_file = sys.argv[3]
    mismatches = int(sys.argv[4])
    read_mode = sys.argv[5] if len(sys.argv) == 6 else "all"

    if read_mode not in ("all", "read1", "read2"):
        print("Errore: il parametro facoltativo deve essere 'all', 'read1' o 'read2'")
        sys.exit(1)

    compute_identity(bam_file, vcf_file, output_file, mismatches, read_mode)
