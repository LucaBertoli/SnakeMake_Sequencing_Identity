#!/usr/bin/env python3
import pysam
import re
import sys
import gzip
import multiprocessing as mp
from collections import defaultdict
import os


def load_variant_positions(vcf_path):
    """Carica dizionari SNV e INDEL da VCF (stessa logica di IdentityRevelations)."""
    snv_dict = defaultdict(set)
    indel_dict = defaultdict(set)

    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            chrom = rec.contig
            pos = rec.pos  # come nello script 3
            ref = rec.ref
            for alt in rec.alts:
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


def get_mismatch_positions(read, snv_dict, chrom):
    """Ritorna solo mismatch NON in VCF (stessa logica SNV di prima, coerente con script 3)."""
    mismatches = set()
    ref_pos = read.reference_start  # 0-based
    md = read.get_tag('MD') if read.has_tag('MD') else ''

    for match in re.finditer(r"(\d+)|([A-Z]|\^[A-Z]+)", md):
        if match.group(1):
            ref_pos += int(match.group(1))
        else:
            if match.group(2).startswith('^'):  # deletion nel MD → skip
                ref_pos += len(match.group(2)) - 1
            else:
                # controllo SNV come mismatch_filtered nello script 3
                if (ref_pos + 1) not in snv_dict[chrom]:
                    mismatches.add(ref_pos)
                ref_pos += 1

    return mismatches


def process_chromosome(args):
    bam_path, chrom, error_classes, snv_dict, indel_dict = args
    bam = pysam.AlignmentFile(bam_path, "rb")

    # struttura risultato con split SNV/INDEL e lista BED
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
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < 60: ##aggiunto filtraggio per mapping quality, test post consegna tesi
            continue

        chrom_ref = bam.get_reference_name(read.reference_id)

        # SNV filtrate come prima (logica già coerente con script 3)
        mismatches = get_mismatch_positions(read, snv_dict, chrom_ref)

        # --- INDEL: logica di filtraggio come nello script di identità (script 3) ---
        ref_pos = read.reference_start  # 0-based
        ins_positions = []   # liste di eventi (start, end)
        del_intervals = []   # liste di eventi (start, end, 'D')

        for op, length in read.cigartuples:
            if op in (0, 7, 8):  # match/mismatch
                ref_pos += length

            elif op == 1:  # insertion
                # come nello script 3: si usa ref_pos "tal quale"
                ins_positions.append((ref_pos, ref_pos + length))
                # ref_pos NON avanza

            elif op == 2:  # deletion
                del_intervals.append((ref_pos, ref_pos + length, 'D'))
                ref_pos += length

        # filtraggio eventi confrontando direttamente con indel_dict, come script 3
        ins_filtered = [i for i in ins_positions
                        if (i[0], i[1], 'I') not in indel_dict[chrom_ref]]
        del_filtered = [d for d in del_intervals
                        if (d[0], d[1], 'D') not in indel_dict[chrom_ref]]

        # basi degli indel filtrati (puntiformi sul ref, 0-based)
        indel_filtered_bases = set()
        for start, end in ins_filtered:
            for pos in range(start, end):
                indel_filtered_bases.add(pos)

        deleted_bases_filtered = 0
        for start, end, _ in del_filtered:
            for pos in range(start, end):
                indel_filtered_bases.add(pos)
            deleted_bases_filtered += (end - start)

        # totale errori = SNV filtrate + basi degli indel filtrati
        total_errors = len(mismatches) + len(indel_filtered_bases)

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

        quals = read.query_qualities
        aligned_pairs = read.get_aligned_pairs(matches_only=False)

        for qpos, rpos in aligned_pairs:
            if qpos is None or rpos is None:
                continue

            q = quals[qpos]
            qbin = "Q0_19" if q < 20 else "Q20_29" if q < 30 else "Q30plus"

            is_snv   = (rpos in mismatches)
            is_indel = (rpos in indel_filtered_bases)

            if is_snv:
                result[err_key]["error"]["SNV"][qbin] += 1
                result[err_key]["error"][qbin] += 1
                result[err_key]["error"]["BED"].append((chrom_ref, rpos, rpos + 1))

            elif is_indel:
                result[err_key]["error"]["INDEL"][qbin] += 1
                result[err_key]["error"][qbin] += 1
                result[err_key]["error"]["BED"].append((chrom_ref, rpos, rpos + 1))

            else:
                result[err_key]["correct"][qbin] += 1

        # nel TSV "deleted_bases" mettiamo solo delezioni filtrate (non-GIAB)
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
            # correct
            for key in merged[e]["correct"]:
                merged[e]["correct"][key] += partial[e]["correct"][key]

            # error totals
            for key in ("Q0_19", "Q20_29", "Q30plus", "del"):
                merged[e]["error"][key] += partial[e]["error"][key]

            # split SNV/INDEL
            for qbin in ("Q0_19", "Q20_29", "Q30plus"):
                merged[e]["error"]["SNV"][qbin] += partial[e]["error"]["SNV"][qbin]
                merged[e]["error"]["INDEL"][qbin] += partial[e]["error"]["INDEL"][qbin]

            # BED
            merged[e]["error"]["BED"].extend(partial[e]["error"]["BED"])

    return merged


def compute_parallel(bam_path, vcf_path, output_path, threads):
    snv_dict, indel_dict = load_variant_positions(vcf_path)

    bam = pysam.AlignmentFile(bam_path, "rb")
    chroms = list(bam.references)
    bam.close()

    error_classes = ["0_errors", "1_error", "2_errors", "3_errors", "gt3_errors"]
    args = [(bam_path, c, error_classes, snv_dict, indel_dict) for c in chroms]

    print(f">>> Running on {len(chroms)} chroms using {threads} processes")

    with mp.Pool(threads) as pool:
        results = pool.map(process_chromosome, args)

    final = merge_results(results, error_classes)

    # --- TSV come prima ---
    with gzip.open(output_path, "wt") as out:
        out.write("error_class\tcorrect_or_error\tQ0_19\tQ20_29\tQ30plus\t"
                  "SNV_Q0_19\tSNV_Q20_29\tSNV_Q30plus\t"
                  "INDEL_Q0_19\tINDEL_Q20_29\tINDEL_Q30plus\tdeleted_bases\n")

        for ec in error_classes:
            row_corr = final[ec]["correct"]
            out.write(f"{ec}\tcorrect\t"
                      f"{row_corr['Q0_19']}\t{row_corr['Q20_29']}\t{row_corr['Q30plus']}\t"
                      f"0\t0\t0\t0\t0\t0\t{row_corr['del']}\n")

            row_err = final[ec]["error"]
            snv = row_err["SNV"]
            indel = row_err["INDEL"]
            out.write(f"{ec}\terror\t"
                      f"{row_err['Q0_19']}\t{row_err['Q20_29']}\t{row_err['Q30plus']}\t"
                      f"{snv['Q0_19']}\t{snv['Q20_29']}\t{snv['Q30plus']}\t"
                      f"{indel['Q0_19']}\t{indel['Q20_29']}\t{indel['Q30plus']}\t"
                      f"{row_err['del']}\n")

    print(f"✔ Done → Saved: {output_path}")

    # --- BED: stessa base, ma salviamo solo per 1,2,3,>3 errori ---
    # prefix per la cartella BED
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
