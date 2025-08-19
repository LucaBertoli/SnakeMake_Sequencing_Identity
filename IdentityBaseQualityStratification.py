import pandas as pd
import numpy as np
import os
import sys
import pysam
import gzip
import re
from collections import defaultdict

####definizione di un dizionario con le posizioni varianti nel VCF in input
def load_variant_positions(vcf_path):
    #definisce due dizionari per le posizioni delle varianti: uno per gli SNV e uno per gli INDEL.
    print("inizio caricamento delle posizioni delle varianti dal VCF:", vcf_path)
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

            # per ogni allele della variante, controlla se è un SNV o un INDEL.
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
    
    print("fine caricamento delle posizioni delle varianti dal VCF:", vcf_path)
    return snv_dict, indel_dict

####estrazione mismatch da tag MD
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

def calculate_stratified_identity(bam, vcf, output):

    snv_dict, indel_dict = load_variant_positions(vcf)

    base_quality_dictionary = {
            str(i): {'match': 0, 'match_correct' : 0, 'match_error': 0, 'mismatch': 0, 'mismatch_correct': 0, 'mismatch_error': 0, 'insertion': 0, 'insertion_correct': 0, 'insertion_error': 0}
            for i in range(0, 51)
        }

    bam = pysam.AlignmentFile(bam, "rb")

    with gzip.open(output, 'wt') as out:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            chrom = bam.get_reference_name(read.reference_id)
            read_start = read.reference_start
            read_alignment_length = read.query_alignment_length
            mismatch_positions = get_mismatch_positions(read)
            aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=False, with_cigar=True)
            mismatch_positions = set(get_mismatch_positions(read))

            print(read.query_qualities)
            print(aligned_pairs)

            for query_pos, ref_pos, cigar_op in aligned_pairs:
                if cigar_op == 4:  # softclip → ignoriamo
                    continue
                elif query_pos is None and ref_pos is not None: # deletion → non ha la base nella read → saltiamo
                    continue
                elif query_pos is not None and ref_pos is None and cigar_op == 1:  # insertion
                    base_quality_dictionary[str(read.query_qualities[query_pos])]['insertion'] += 1

                    # calcolo della posizione di inserzione e la confronto con le posizioni di inserzioni nel VCF.
                    insertion_pos = read.reference_start + query_pos
                    print(insertion_pos)
                    match_found = any(
                        insertion_pos-1 >= start and insertion_pos-1 <= end and typ == 'I'
                        for (start, end, typ) in indel_dict[chrom]
                    )
                    if match_found:
                        base_quality_dictionary[str(read.query_qualities[query_pos])]['insertion_correct'] += 1
                    else:
                        base_quality_dictionary[str(read.query_qualities[query_pos])]['insertion_error'] += 1

                elif query_pos is not None and ref_pos is not None: # match or mismatch
                    q = read.query_qualities[query_pos]
                    if ref_pos in mismatch_positions:                                   # mismatch
                        base_quality_dictionary[str(q)]['mismatch'] += 1
                        if ref_pos not in snv_dict[chrom]:                                  #non nel vcf
                            base_quality_dictionary[str(q)]['mismatch_error'] += 1
                        else:                                                               # nel vcf
                            base_quality_dictionary[str(q)]['mismatch_correct'] += 1
                    else:                                                               # match
                        base_quality_dictionary[str(q)]['match'] += 1
                        if ref_pos not in snv_dict[chrom]:                                  # non nel vcf
                            base_quality_dictionary[str(q)]['match_correct'] += 1
                        else:                                                               # nel vcf
                            base_quality_dictionary[str(q)]['match_error'] += 1


        out.write("BaseQuality\tMatch\tMatch_correct\tMatch_error\tMismatch\tMismatch_correct\tMismatch_error\tInsertion\tInsertion_correct\tInsertion_error\n")
        for q in sorted(base_quality_dictionary.keys(), key=lambda x: int(x)):
            values = base_quality_dictionary[q]
            out.write(f"{q}\t{values['match']}\t{values['match_correct']}\t{values['match_error']}\t{values['mismatch']}\t{values['mismatch_correct']}\t{values['mismatch_error']}\t{values['insertion']}\t{values['insertion_correct']}\t{values['insertion_error']}\n")

    bam.close()
    out.close()

if __name__ == "__main__":
    print("INIZIO del calcolo della identità stratificando i dati per qualità delle basi")

    bam=sys.argv[1]
    VCF=sys.argv[2]
    output=sys.argv[3]

    if not os.path.exists(bam):
        print("Il file pileup non esiste:", bam)
        sys.exit(1)

    if not os.path.exists(VCF):
        print("Il file VCF non esiste:", VCF)
        sys.exit(1)

    print("BAM:", bam)
    print("VCF:", VCF)
    print("Output:", output)

    calculate_stratified_identity(bam, VCF, output)

    print("FINE calcolo dell'identità...")