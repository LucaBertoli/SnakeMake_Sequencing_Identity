#!/usr/bin/env python3

"""
Extract per-read quality score distributions stratified by error count.

For each read in Q30+, outputs:
- read_id
- num_errors: total errors in the read (0, 1, 2, 3, >3)
- pct_Q0_19: % of bases with Q < 20
- pct_Q20_29: % of bases with 20 <= Q < 30
- pct_Q30plus: % of bases with Q >= 30

This TSV is used to generate violin plots of quality score distributions.
"""

import pysam
import re
import sys
from collections import defaultdict
import math
import gzip
import multiprocessing as mp
import os


def load_variant_positions(vcf_path):
    """Load SNV and INDEL positions from VCF (to exclude known variants)."""
    snv_dict = defaultdict(set)
    indel_dict = defaultdict(set)
    
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            chrom = rec.contig
            pos = rec.pos
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
    
    return snv_dict, indel_dict


def get_mismatch_positions(read, snv_dict, chrom):
    """Extract mismatch positions from MD tag, excluding known SNVs."""
    mismatches = set()
    ref_pos = read.reference_start
    md = read.get_tag('MD') if read.has_tag('MD') else ''
    
    for match in re.finditer(r"(\d+)|([A-Z]|\^[A-Z]+)", md):
        if match.group(1):
            ref_pos += int(match.group(1))
        else:
            if match.group(2).startswith('^'):
                ref_pos += len(match.group(2)) - 1
            else:
                if (ref_pos + 1) not in snv_dict[chrom]:
                    mismatches.add(ref_pos)
                ref_pos += 1
    
    return mismatches


def process_bam_chrom(args):
    """
    Process BAM file and generate per-read quality distribution TSV.
    """

    chrom, bam_path, snv_dict, indel_dict, output_folder = args
    print(f'Processing BAM for read quality distribution for chromosome {chrom}...')
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # Output file
    output_path = f"{output_folder}/read_quality_dist_{chrom}.tsv.gz"
    os.makedirs(output_folder, exist_ok=True)
    with gzip.open(output_path, 'wt') as out:
        # Header
        out.write("read_id\tnum_errors\terror_class\tq0_19\tq20_29\tq30plus\tpct_Q0_19\tpct_Q20_29\tpct_Q30plus\tread_length\tmean_Q\n")
        
        read_count = 0
        processed_reads = 0
        

        for read in bam.fetch(chrom):
            read_count += 1
            
            # Filter: unmapped, secondary, supplementary, low MQ
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 60:
                continue
                
            processed_reads += 1
            chrom_ref = bam.get_reference_name(read.reference_id)
            quals = read.query_qualities
            read_length = read.query_length
            
            if read_length == 0 or not quals:
                continue
            
            # Get mismatch positions
            mismatches = get_mismatch_positions(read, snv_dict, chrom_ref)
            
            # Parse INDELs
            ref_pos = read.reference_start
            query_pos = 0
            ins_positions = []
            del_intervals = []
            
            if read.cigartuples:
                for op, length in read.cigartuples:
                    if op in (0, 7, 8):  # M, =, X
                        ref_pos += length
                        query_pos += length
                    elif op == 1:  # I
                        ins_positions.append((ref_pos, ref_pos + length))
                        query_pos += length
                    elif op == 2:  # D
                        del_intervals.append((ref_pos, ref_pos + length, 'D'))
                        ref_pos += length
            
            # Filter insertions/deletions to exclude known variants
            ins_filtered = [i for i in ins_positions if (i[0], i[1], 'I') not in indel_dict[chrom_ref]]
            del_filtered = [d for d in del_intervals if (d[0], d[1], 'D') not in indel_dict[chrom_ref]]
            
            # Count deletion bases
            deleted_bases_filtered = sum((end - start) for start, end, _ in del_filtered)
            
            # Count insertion errors
            inserted_errors = 0
            ref_pos = read.reference_start
            query_pos = 0
            
            if read.cigartuples:
                for op, length in read.cigartuples:
                    if op in (0, 7, 8):
                        ref_pos += length
                        query_pos += length
                    elif op == 1:
                        is_error_ins = (ref_pos, ref_pos + length) in ins_filtered
                        if is_error_ins:
                            inserted_errors += length
                        query_pos += length
                    elif op == 2:
                        ref_pos += length
            
            # Total errors
            total_errors = len(mismatches) + deleted_bases_filtered + inserted_errors
            
            # Classify error count
            if total_errors == 0:
                error_class = "0_errors"
            elif total_errors == 1:
                error_class = "1_error"
            elif total_errors == 2:
                error_class = "2_errors"
            elif total_errors == 3:
                error_class = "3_errors"
            else:
                error_class = "gt3_errors"
            
            # Count bases by quality bin
            q0_19 = sum(1 for q in quals if q < 20)
            q20_29 = sum(1 for q in quals if 20 <= q < 30)
            q30plus = sum(1 for q in quals if q >= 30)
            
            # Percentages
            pct_q0_19 = (q0_19 / read_length) * 100 if read_length > 0 else 0
            pct_q20_29 = (q20_29 / read_length) * 100 if read_length > 0 else 0
            pct_q30plus = (q30plus / read_length) * 100 if read_length > 0 else 0
            
            # Mean Q
            mean_q = -10 * math.log10(sum(10 ** (-q / 10) for q in quals) / len(quals)) if quals else 0
            
            # Write row
            out.write(
                f"{read.query_name}\t{total_errors}\t{error_class}\t"
                f"{q0_19}\t{q20_29}\t{q30plus}\t"
                f"{pct_q0_19:.4f}\t{pct_q20_29:.4f}\t{pct_q30plus:.4f}\t"
                f"{read_length}\t{mean_q:.4f}\n"
            )
        
        bam.close()
        print(f"✔ Processed {chrom}, output saved: {output_path}")

def compute_parallel(bam_path, vcf_path, output_folder, threads):
    snv_dict, indel_dict = load_variant_positions(vcf_path)

    bam = pysam.AlignmentFile(bam_path, "rb")
    chroms = list(bam.references)
    bam.close()

    args = [(chrom, bam_path, snv_dict, indel_dict, output_folder) for chrom in chroms]
    with mp.Pool(threads) as pool:
        pool.map(process_bam_chrom, args)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("USAGE: python read_quality_dist.py input.bam input.vcf.gz output_folder threads")
        sys.exit(1)
    
    bam_path = sys.argv[1]
    vcf_path = sys.argv[2]
    output_folder = sys.argv[3]
    threads = int(sys.argv[4])
    
    print("Starting read quality distribution computation...")
    print(f"Input BAM: {bam_path}")
    print(f"Input VCF: {vcf_path}")
    print(f"Output folder: {output_folder}")
    print(f"Threads: {threads}")
    compute_parallel(bam_path, vcf_path, output_folder, threads)