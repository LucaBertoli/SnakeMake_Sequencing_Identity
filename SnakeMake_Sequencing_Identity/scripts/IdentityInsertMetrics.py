#!/usr/bin/env python3
import pysam
import re
import math
import numpy as np
import gzip
from collections import defaultdict
import argparse

# -----------------------------
# Funzioni di supporto
# -----------------------------

def load_variant_positions(vcf_path):
    """Carica posizioni SNV e INDEL dal VCF (come nel tuo script originale)."""
    snv_dict = defaultdict(set)
    indel_dict = defaultdict(set)
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            chrom, pos, ref, alts = rec.contig, rec.pos, rec.ref, rec.alts
            for alt in alts:
                if len(ref) == 1 and len(alt) == 1:
                    snv_dict[chrom].add(pos)
                elif len(ref) > len(alt):  # deletion
                    indel_dict[chrom].add((pos, pos + len(ref) - len(alt), 'D'))
                elif len(ref) < len(alt):  # insertion
                    indel_dict[chrom].add((pos, pos + len(alt) - len(ref), 'I'))
    return snv_dict, indel_dict


def get_mismatch_positions(read):
    """Estrae posizioni dei mismatch dal tag MD."""
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


def phred_log_mean(quals):
    """Media logaritmica dei Q-score (come -10*log10(mean(10^-Q/10)))."""
    if not quals:
        return np.nan
    probs = [10 ** (-q / 10) for q in quals]
    return -10 * math.log10(np.mean(probs))


def compute_identity_filtered_with_indels(read, snv_dict, indel_dict, bam):
    """Calcola identity_filtered_with_indels come nel tuo script originale."""
    chrom = bam.get_reference_name(read.reference_id)
    mismatch_positions = get_mismatch_positions(read)
    mismatch_filtered = [p for p in mismatch_positions if (p + 1) not in snv_dict[chrom]]

    # parse CIGAR per inserzioni e delezioni
    ref_pos = read.reference_start
    ins_positions, del_intervals = [], []
    inserted_bases = deleted_bases = 0
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

    aligned_length_total = read.query_alignment_length
    if aligned_length_total == 0:
        return np.nan

    ins_filtered = [i for i in ins_positions if (i[0], i[1], 'I') not in indel_dict[chrom]]
    del_filtered = [d for d in del_intervals if (d[0], d[1], 'D') not in indel_dict[chrom]]

    num_indel_filtered_bases = (
        sum(i[1] - i[0] for i in ins_filtered)
        + sum(d[1] - d[0] for d in del_filtered)
    )

    identity_filtered_with_indels = 1 - (
        (len(mismatch_filtered) + num_indel_filtered_bases) / aligned_length_total
    )
    return identity_filtered_with_indels


# -----------------------------
# Funzione principale
# -----------------------------

def collect_insert_identity_allreads(bam_path, vcf_path, output_path):
    snv_dict, indel_dict = load_variant_positions(vcf_path)
    bam = pysam.AlignmentFile(bam_path, "rb")

    with gzip.open(output_path, "wt") as out:
        out.write("read_name\tread_type\tinsert_size\tmean_qual_log\tidentity_filtered_with_indels\n")

        for read in bam:
            # considera solo allineamenti primari e mappati
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 60:
                continue

            # determina tipo di read
            if read.is_paired:
                read_type = "R1" if read.is_read1 else "R2"
            else:
                read_type = "S"

            # calcolo insert size
            insert_size = (
                abs(read.template_length)
                if read.is_paired and read.template_length != 0
                else read.query_length
            )

            # qualità media log
            mean_q = phred_log_mean(read.query_qualities)

            # identità filtrata con indel
            idfwi = compute_identity_filtered_with_indels(read, snv_dict, indel_dict, bam)

            out.write(f"{read.query_name}\t{read_type}\t{insert_size}\t{mean_q:.3f}\t{idfwi:.6f}\n")

    bam.close()
    print(f"Output scritto in {output_path}")


# -----------------------------
# Entry point
# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calcola insert size, qualità media e identità filtrata con indel per tutte le read mappate (R1/R2/S)"
    )
    parser.add_argument("-b", "--bam", required=True, help="BAM di input (con tag MD)")
    parser.add_argument("-v", "--vcf", required=True, help="VCF germinale (BGZipped + index)")
    parser.add_argument("-o", "--output", required=True, help="TSV.gz di output")
    args = parser.parse_args()

    collect_insert_identity_allreads(args.bam, args.vcf, args.output)

