#!/bin/bash

R1=$1;
R2=$2;
thread=$3;

fastp="/home/tools/fastp/fastp";
# FASTA="/home/db/hg38_bwa2/hg38_chr1-22-XYM/Homo_sapiens_assembly38_noalt.fasta";
FASTA="/home/db/hg38_bwa2/hg38_noalt_nodecoy/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"; # per test rianalisi

echo "Start trimming at $(date)"; 

$fastp -i $R1 -I $R2 -o trimmed1.fastq.gz -O trimmed2.fastq.gz --failed_out failed.fastq.gz --overrepresentation_analysis -h adapters_removal_report.html -w $thread -g --disable_quality_filtering --disable_length_filtering
perl /home/script/mapBWA2.pl $FASTA trimmed1.fastq.gz trimmed2.fastq.gz $thread

echo "Ended trimming at $(date)";