# Lo script ScatterBam.py è stato creato per calcolare le qualità medie delle read in un file BAM e per separare le read in due file BAM distinti in base alla loro qualità media. 
# Le qualità considere sono le qualità in phread score delle basi allineate (esclude le softclipped).

import pysam
import sys
import gzip

def process_bam(input_bam, mode):
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    write_bams = mode in ("BOTH", "BAM")
    write_tsv = mode in ("BOTH", "QUAL")

    #condizione per aprire in modalità scrittura dei file bam (da utilzizare se si vogliono generare bam con soltanto read di una certa qualità)
    if write_bams:
        bam_out_30_up = pysam.AlignmentFile(input_bam + "_reads_30_up.bam", "wb", template=bam_in)
        # bam_out_30_35 = pysam.AlignmentFile(input_bam + "_reads_30_35.bam", "wb", template=bam_in)
        # bam_out_35_up = pysam.AlignmentFile(input_bam + "_reads_35_up.bam", "wb", template=bam_in)

    if write_tsv:
        out_tsv = gzip.open(input_bam + "_mean_quality.tsv.gz", "wt")
        out_tsv.write("read_name\tmean_quality\n")

    for read in bam_in.fetch(until_eof=True):
        quals = read.query_alignment_qualities
        if not quals:
            continue

        mean_qual = sum(quals) / len(quals)

        if write_tsv:
            out_tsv.write(f"{read.query_name}\t{mean_qual:.2f}\n")

        #condizione per scrivere uno o piu bam con read di una certa qualità
        if write_bams:
            if 30 <= mean_qual:
                bam_out_30_up.write(read)
            # if 30 <= mean_qual < 35:
            #     bam_out_30_35.write(read)
            # elif mean_qual >= 35:
            #     bam_out_35_up.write(read)

    bam_in.close()

    if write_bams:
        bam_out_30_35.close()
        bam_out_35_up.close()
    if write_tsv:
        out_tsv.close()

def main():
    if len(sys.argv) != 3:
        print(f"Uso: python ScatterBam.py input.bam BAM|QUAL|BOTH")
        sys.exit(1)

    input_bam = sys.argv[1]
    mode = sys.argv[2].upper()

    print(f"Input BAM: {input_bam}")
    print(f"Mode: {mode}")

    process_bam(input_bam, mode)

if __name__ == "__main__":
    main()