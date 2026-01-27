#!/usr/bin/env python3
import pysam
import re
import sys
import gzip
import multiprocessing as mp
from collections import defaultdict
import os


def load_variant_positions(vcf_path):
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


def process_chromosome(args):
    bam_path, chrom, error_classes, snv_dict, indel_dict = args
    bam = pysam.AlignmentFile(bam_path, "rb")

    result = {
        err: {
            "correct": {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0, "del": 0},
            "error": {
                "Q0_19": 0, "Q20_29": 0, "Q30plus": 0, "del": 0,
                "SNV":   {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0},
                "INDEL": {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0},
                "BED": []
            }
        } for err in error_classes
    }

    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 60:
            continue

        chrom_ref = bam.get_reference_name(read.reference_id)
        quals = read.query_qualities

        mismatches = get_mismatch_positions(read, snv_dict, chrom_ref)

        # ---------------- INDEL parsing ----------------
        ref_pos = read.reference_start
        query_pos = 0

        ins_positions = []
        del_intervals = []

        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                ref_pos += length
                query_pos += length
            elif op == 1:
                ins_positions.append((ref_pos, ref_pos + length))
                query_pos += length
            elif op == 2:
                del_intervals.append((ref_pos, ref_pos + length, 'D'))
                ref_pos += length

        ins_filtered = [i for i in ins_positions if (i[0], i[1], 'I') not in indel_dict[chrom_ref]]
        del_filtered = [d for d in del_intervals if (d[0], d[1], 'D') not in indel_dict[chrom_ref]]

        # ---------------- COUNTING DELETION ERROR BASES ----------------
        del_filtered_bases = set()
        deleted_bases_filtered = 0

        for start, end, _ in del_filtered:
            for pos in range(start, end):
                del_filtered_bases.add(pos)
            deleted_bases_filtered += (end - start)

        # ---------------- SAVING INSERTION ERROR AND CORRECT BASES (BASE QUALITY STRATIFIED) ----------------
        inserted_errors_by_qbin = {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0}
        inserted_correct_by_qbin = {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0}

        ref_pos = read.reference_start
        query_pos = 0

        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                ref_pos += length
                query_pos += length
            elif op == 1:
                is_error_ins = (ref_pos, ref_pos + length) in ins_filtered
                for i in range(query_pos, query_pos + length):
                    q = quals[i]
                    qbin = "Q0_19" if q < 20 else "Q20_29" if q < 30 else "Q30plus"
                    if is_error_ins:
                        inserted_errors_by_qbin[qbin] += 1
                    else:
                        inserted_correct_by_qbin[qbin] += 1
                query_pos += length
            elif op == 2:
                ref_pos += length

        # ---------------- STRATIFICATION OF THE READS BY ERROR NUMBER ----------------
        total_errors = len(mismatches) + len(del_filtered_bases) + sum(inserted_errors_by_qbin.values())

        if total_errors == 0:
            err_key = "0_errors"
        elif total_errors == 1:
            err_key = "1_error"
        elif total_errors == 2:
            err_key = "2_errors"
        elif total_errors == 3:
            err_key = "3_errors"
        else:
            err_key = "gt3_errors"

        aligned_pairs = read.get_aligned_pairs(matches_only=False)
        
        # ---------------- POPULATION OF ERROR AND CORRECT BASES RESULT DATA STRUCTURE ----------------
        #counting mismatch errors and correct bases 
        for qpos, rpos in aligned_pairs:
            if qpos is None or rpos is None:
                continue
            q = quals[qpos]
            qbin = "Q0_19" if q < 20 else "Q20_29" if q < 30 else "Q30plus"
            if rpos in mismatches:
                result[err_key]["error"]["SNV"][qbin] += 1
                result[err_key]["error"][qbin] += 1
                result[err_key]["error"]["BED"].append((chrom_ref, rpos, rpos + 1))
            else:
                result[err_key]["correct"][qbin] += 1

        #counting insertion errors and correct bases
        for qbin in inserted_correct_by_qbin:
            result[err_key]["correct"][qbin] += inserted_correct_by_qbin[qbin]

        for qbin in inserted_errors_by_qbin:
            result[err_key]["error"]["INDEL"][qbin] += inserted_errors_by_qbin[qbin]
            result[err_key]["error"][qbin] += inserted_errors_by_qbin[qbin]

        #deletions counting deletion errors
        result[err_key]["error"]["del"] += deleted_bases_filtered

    bam.close()
    print(f"[OK] {chrom}")
    return result


def merge_results(all_results, error_classes):
    merged = {
        err: {
            "correct": {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0, "del": 0},
            "error": {
                "Q0_19": 0, "Q20_29": 0, "Q30plus": 0, "del": 0,
                "SNV":   {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0},
                "INDEL": {"Q0_19": 0, "Q20_29": 0, "Q30plus": 0},
                "BED": []
            }
        } for err in error_classes
    }

    for partial in all_results:
        for e in partial:
            for key in merged[e]["correct"]:
                merged[e]["correct"][key] += partial[e]["correct"][key]
            for key in ("Q0_19", "Q20_29", "Q30plus", "del"):
                merged[e]["error"][key] += partial[e]["error"][key]
            for qbin in ("Q0_19", "Q20_29", "Q30plus"):
                merged[e]["error"]["SNV"][qbin] += partial[e]["error"]["SNV"][qbin]
                merged[e]["error"]["INDEL"][qbin] += partial[e]["error"]["INDEL"][qbin]
            merged[e]["error"]["BED"].extend(partial[e]["error"]["BED"])
    return merged


def compute_parallel(bam_path, vcf_path, output_path, threads):
    snv_dict, indel_dict = load_variant_positions(vcf_path)

    bam = pysam.AlignmentFile(bam_path, "rb")
    chroms = list(bam.references)
    bam.close()

    error_classes = ["0_errors", "1_error", "2_errors", "3_errors", "gt3_errors"]
    args = [(bam_path, c, error_classes, snv_dict, indel_dict) for c in chroms]

    with mp.Pool(threads) as pool:
        results = pool.map(process_chromosome, args)

    final = merge_results(results, error_classes)

    with gzip.open(output_path, "wt") as out:
        out.write(
            "error_class\tcorrect_or_error\tQ0_19\tQ20_29\tQ30plus\t"
            "SNV_Q0_19\tSNV_Q20_29\tSNV_Q30plus\t"
            "INDEL_Q0_19\tINDEL_Q20_29\tINDEL_Q30plus\tdeleted_bases\n"
        )

        for ec in error_classes:
            row_corr = final[ec]["correct"]
            out.write(
                f"{ec}\tcorrect\t{row_corr['Q0_19']}\t{row_corr['Q20_29']}\t{row_corr['Q30plus']}\t"
                f"0\t0\t0\t0\t0\t0\t{row_corr['del']}\n"
            )

            row_err = final[ec]["error"]
            snv = row_err["SNV"]
            indel = row_err["INDEL"]
            out.write(
                f"{ec}\terror\t{row_err['Q0_19']}\t{row_err['Q20_29']}\t{row_err['Q30plus']}\t"
                f"{snv['Q0_19']}\t{snv['Q20_29']}\t{snv['Q30plus']}\t"
                f"{indel['Q0_19']}\t{indel['Q20_29']}\t{indel['Q30plus']}\t"
                f"{row_err['del']}\n"
            )
    print(f"✔ Done → Saved: {output_path}")
    
    # --- BED: scrittura delle coordinate SNV per 1,2,3,>3 errori ---
    prefix = output_path
    if prefix.endswith(".tsv.gz"):
        prefix = prefix[:-7]
    elif prefix.endswith(".gz"):
        prefix = prefix[:-3]

    bed_dir = f"{prefix}_BED"
    os.makedirs(bed_dir, exist_ok=True)

    for ec in ["1_error", "2_errors", "3_errors", "gt3_errors"]:
        bed_path = os.path.join(bed_dir, f"{ec}.bed")
        with open(bed_path, "w") as bed:
            for chrom, start, end in final[ec]["error"]["BED"]:
                bed.write(f"{chrom}\t{start}\t{end}\n")
        print(f"✔ BED scritto → {bed_path}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("USAGE: python script.py input.bam input.vcf.gz output.tsv.gz threads")
        sys.exit(1)

    compute_parallel(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
