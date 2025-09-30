#!/usr/bin/env python3

"""
Questo script calcola l'identità di ciascuna read allineata in un file BAM, escludendo le basi sovrapposte a varianti note (da un VCF)
e tenendo conto anche di inserzioni e delezioni nel calcolo dell'identità estesa.

Input:
- BAM con reads allineate (con tag MD)
- VCF con varianti germinali note (SNV e INDEL)

Output:
- TSV gzippato con per-read le metriche di identità e conteggi di varianti.

Formato colonne output (TSV):
- read: nome della read
- identity: identità considerando solo mismatch
- identity_filtered: identità escludendo mismatch in posizioni variant SNV
- identity_with_indels: identità con mismatch + indel
- identity_filtered_with_indels: identità filtrata da SNV/indel
- aligned_length_query: lunghezza dell'allineamento sulla read (escludendo indel)
- aligned_length_total: lunghezza dell'allineamento inclusi indel
- mismatches_total: numero totale di mismatch (da tag MD)
- mismatches_filtered: mismatch non sovrapposti a SNV
- insertions: numero di inserzioni (da CIGAR)
- ins_filtered: inserzioni non presenti come INDEL
- deletions: numero di delezioni (da CIGAR)
- del_filtered: delezioni non presenti come INDEL
"""

import pysam
import re
import gzip
import sys
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
    pos = read.reference_start  # 0-based

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

#funzione centrale dello script.
def compute_identity(bam_path, vcf_path, output_path, read_mode="all"):
    # Carica le posizioni delle varianti SNV e INDEL dal VCF
    snv_dict, indel_dict = load_variant_positions(vcf_path)
    # Apri il file BAM in modalità lettura
    bam = pysam.AlignmentFile(bam_path, "rb")

    #inizializza il file di output per i dati di identità.
    with gzip.open(output_path, 'wt') as out:
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

            #calcola le identità considerando i mismatch e le indel filtrate e non filtrate.
            identity = 1 - (len(mismatch_positions) / aligned_length_total)
            identity_filtered = 1 - (len(mismatch_filtered) / aligned_length_total)
            identity_with_indels = 1 - ((len(mismatch_positions) + inserted_bases + deleted_bases) / aligned_length_total)
            identity_filtered_with_indels = 1 - ((len(mismatch_filtered) + sum(i[1] - i[0] for i in ins_filtered) + sum(d[1] - d[0] for d in del_filtered)) / aligned_length_total)

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

    print(f"Output scritto in {output_path}")
    print("Elaborazione completata.")

if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        print("Uso: python IdentityRevelations.py input.bam input.vcf.gz output.tsv.gz [all|read1|read2]")
        sys.exit(1)

    bam_file = sys.argv[1]
    vcf_file = sys.argv[2]
    output_file = sys.argv[3]
    read_mode = sys.argv[4] if len(sys.argv) == 5 else "all"

    if read_mode not in ("all", "read1", "read2"):
        print("Errore: il parametro facoltativo deve essere 'all', 'read1' o 'read2'")
        sys.exit(1)

    compute_identity(bam_file, vcf_file, output_file, read_mode)
