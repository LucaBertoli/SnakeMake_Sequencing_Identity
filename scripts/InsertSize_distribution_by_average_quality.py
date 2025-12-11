# Lo script è stato creato per estrarre un istogramma contenente la distribuzione di insert size stratificata per qualità media della read. 
# Le qualità considerate sono le qualità in phread score delle basi allineate (esclude le softclipped).
# La qualità media di ciascuna read viene calcolata come logaritmo della media delle probabilità di errore delle basi allineate.

import pysam
import sys
import gzip
import math

def process_bam(input_bam, mode):
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    write_bams = mode in ("BOTH", "BAM")
    write_tsv = mode in ("BOTH", "QUAL")

    #condizione per aprire in modalità scrittura dei file bam (da utilzizare se si vogliono generare bam con soltanto read di una certa qualità)
    if write_bams:
        base_name = input_bam[:-4] if input_bam.endswith(".bam") else input_bam
        bam_out_30_up = pysam.AlignmentFile(base_name + "_reads_30_up_LogOfMeans.bam", "wb", template=bam_in)

    if write_tsv:
        out_tsv = gzip.open(input_bam + "_mean_quality_LogOfMeans.tsv.gz", "wt")
        out_tsv.write("read_name\tmean_quality\n")

    for read in bam_in.fetch(until_eof=True):
        quals = read.query_alignment_qualities
        if not quals:
            continue

        #calcolo la qualità media di ciascuna read come logaritmo della media delle probabilità di errore
        mean_qual = -10 * math.log10(sum(10 ** (-q / 10) for q in quals) / len(quals))

        if write_tsv:
            out_tsv.write(f"{read.query_name}\t{mean_qual:.2f}\n")

        #condizione per scrivere uno o piu bam con read di una certa qualità
        if write_bams:
            if 30 <= mean_qual:
                bam_out_30_up.write(read)

    bam_in.close()

    if write_bams:
        bam_out_30_up.close()
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

    print("Elaborazione completata.")

if __name__ == "__main__":
    main()