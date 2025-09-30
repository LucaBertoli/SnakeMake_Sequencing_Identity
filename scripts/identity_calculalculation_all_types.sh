#!/bin/bash

FASTA=$1
BAM=$2
VCF=$3
TARGET=$4
OUTPUT=$5

if [[ -n $FASTA && -n $BAM && -n $VCF && -n $TARGET ]]; then
  PILEUP_FILE=$(dirname $(realpath $BAM))/$(basename $BAM .bam).${OUTPUT}.all.types.pileup
  # REPORT_FILE=$(dirname $(realpath $BAM))/$(basename $BAM .bam).all.types.pileup.report
  REPORT_FILE=${OUTPUT}.identity.all.types.pileup.report

  NEW_TARGET=$(dirname $(realpath $BAM))/$(basename $BAM .bam).all.types.pileup.SNP.target.bed

  echo "Writing output to: $REPORT_FILE"
  echo "FASTA: $FASTA" > $REPORT_FILE
  echo "BAM: $BAM" >> $REPORT_FILE
  echo "VCF: $VCF" >> $REPORT_FILE
  echo "Target: $TARGET" >> $REPORT_FILE

  bedtools subtract -a $TARGET -b $VCF > $NEW_TARGET

 if [[ ! -f $PILEUP_FILE ]] ; then
  echo "pileup not present"
  samtools mpileup -A -l $NEW_TARGET -f $FASTA -q 0 -Q 0 $BAM | bgzip > $PILEUP_FILE
 else
  echo "pileup present"
 fi

  # Calcolo delle basi totali
  TOTAL_BASES=$(zcat $PILEUP_FILE | awk -F "\t" '{SUM+=$4} END {print SUM}')
  echo "Total Bases: $TOTAL_BASES" >> $REPORT_FILE

  # Calcolo delle basi mismatch
  TOTAL_MM=$(zcat $PILEUP_FILE | awk -F "\t" '$3!="N" {print $5}' | grep -oe A -oe a -oe C -oe c -oe G -oe g -oe T -oe t | wc -l)
  echo "TOTAL mismatches: $TOTAL_MM" >> $REPORT_FILE

  # Calcolo del mismatch rate e dell'identità
  MISMATCH_RATE=$(awk "BEGIN {print ($TOTAL_MM / $TOTAL_BASES)*100}")
  IDENTITY=$(awk "BEGIN {print 100 - $MISMATCH_RATE}")

  echo "Mismatch rate: $MISMATCH_RATE%" >> $REPORT_FILE
  echo "Identity: $IDENTITY%" >> $REPORT_FILE

  # # Conteggio delle singole basi
  # for BASE in A a C c G g T t; do
  #   COUNT=$(zcat $PILEUP_FILE | awk -F "\t" '$3!="N" {print $5}' | grep -o $BASE | wc -l)
  #   echo "$BASE: $COUNT" >> $REPORT_FILE
  # done
else
  echo "missing input file"
  echo "bash identity_calculalculation.sh <FASTA> <BAM> <VCF> <TARGET>"
fi