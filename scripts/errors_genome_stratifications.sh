##bash /home/comparazione_kit/calcolo_identità_sequenze_Sequenziatori/scripts/errors_genome_stratifications.sh 1_error.bed 2_errors.bed 3_errors.bed gt3_errors.bed 
##script per stratificare il numero di errori totali in base alle regioni complesse del genoma umano secondo le stratificazioni GIAB (tandem repeats, homopolymers, GC content extreme, low mappability, segmental duplications)

file_1_errors=$1
file_2_errors=$2
file_3_errors=$3
file_gt3_errors=$4

stratifications_dir="/home/db/hg38_bwa2/genome_stratification/v3.5/GRCh38@all"

tandem_repeats=$stratifications_dir"/LowComplexity/GRCh38_AllTandemRepeats.bed.gz"
homopolymers=$stratifications_dir"/LowComplexity/GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz"
tandem_repeats_and_homopolymers=$stratifications_dir"/LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz"
gclt25orgt65=$stratifications_dir"/GCcontent/GRCh38_gclt25orgt65_slop50.bed.gz"
gclt30orgt55=$stratifications_dir"/GCcontent/GRCh38_gclt30orgt55_slop50.bed.gz"
low_map=$stratifications_dir"/Mappability/GRCh38_lowmappabilityall.bed.gz"
seg_dup=$stratifications_dir"/SegmentalDuplications/GRCh38_segdups.bed.gz"

total_bp_errors=$(cat ${file_1_errors} ${file_2_errors} ${file_3_errors} ${file_gt3_errors} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_tandem_repeats=$(cat ${file_1_errors} ${file_2_errors} ${file_3_errors} ${file_gt3_errors} | bedtools intersect -a - -b "$tandem_repeats" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_homopolymers=$(cat ${file_1_errors} ${file_2_errors} ${file_3_errors} ${file_gt3_errors} | bedtools intersect -a - -b "$homopolymers" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_tandem_repeats_and_homopolymers=$(cat ${file_1_errors} ${file_2_errors} ${file_3_errors} ${file_gt3_errors} | bedtools intersect -a - -b "$tandem_repeats_and_homopolymers" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_gclt25orgt65=$(cat "${file_1_errors}" "${file_2_errors}" "${file_3_errors}" "${file_gt3_errors}" | bedtools intersect -a - -b "$gclt25orgt65" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_gclt30orgt55=$(cat "${file_1_errors}" "${file_2_errors}" "${file_3_errors}" "${file_gt3_errors}" | bedtools intersect -a - -b "$gclt30orgt55" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_low_map=$(cat "${file_1_errors}" "${file_2_errors}" "${file_3_errors}" "${file_gt3_errors}" | bedtools intersect -a - -b "$low_map" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_seg_dup=$(cat ${file_1_errors} ${file_2_errors} ${file_3_errors} ${file_gt3_errors} | bedtools intersect -a - -b "$seg_dup" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')

echo -e "total_errors\ttandem_repeats\thomopolymers\ttandem_repeats_and_homopolymers\tgclt25orgt65\tgclt30orgt55\tlow_map\tseg_dup"
echo -e "$total_bp_errors\t$count_tandem_repeats\t$count_homopolymers\t$count_tandem_repeats_and_homopolymers\t$count_gclt25orgt65\t$count_gclt30orgt55\t$count_low_map\t$count_seg_dup"
echo ""
echo "Only reads with 1 error"
total_bp_errors=$(cat ${file_1_errors} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_tandem_repeats=$(cat ${file_1_errors} | bedtools intersect -a - -b "$tandem_repeats" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_homopolymers=$(cat ${file_1_errors} | bedtools intersect -a - -b "$homopolymers" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_tandem_repeats_and_homopolymers=$(cat ${file_1_errors} | bedtools intersect -a - -b "$tandem_repeats_and_homopolymers" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_gclt25orgt65=$(cat "${file_1_errors}" | bedtools intersect -a - -b "$gclt25orgt65" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_gclt30orgt55=$(cat "${file_1_errors}" | bedtools intersect -a - -b "$gclt30orgt55" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_low_map=$(cat "${file_1_errors}" | bedtools intersect -a - -b "$low_map" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
count_seg_dup=$(cat ${file_1_errors} | bedtools intersect -a - -b "$seg_dup" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')

echo -e "total_errors\ttandem_repeats\thomopolymers\ttandem_repeats_and_homopolymers\tgclt25orgt65\tgclt30orgt55\tlow_map\tseg_dup"
echo -e "$total_bp_errors\t$count_tandem_repeats\t$count_homopolymers\t$count_tandem_repeats_and_homopolymers\t$count_gclt25orgt65\t$count_gclt30orgt55\t$count_low_map\t$count_seg_dup"