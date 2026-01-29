#nohup bash /home/comparazione_kit/calcolo_identità_sequenze_Sequenziatori/scripts/count_errors_by_type.sh IDENTITY_NA12878_KAPA_reads_30_up_LogOfMeans_all.tsv.gz > IDENTITY_NA12878_KAPA_reads_30_up_LogOfMeans_all.tsv.gz.error.stats &
#CALCOLO NEL NUMERO DI ERRORI FILTRATI, DIVISI PER TIPOLOGIA DI ERRORE

error_file=$1

zcat $error_file | awk '
BEGIN {
    printf "total_bases\ttotal_filtered_errors\ttotal_filtered_error_bases\tfiltered_mismatches\tfiltered_insertions\tfiltered_insertions_bases\tfiltered_deletions\tfiltered_deletions_bases\n"
}
NR>1 {
    bases += $7
    mm_filt += $9
    ins_filt += $12
    ins_filt_bases += $13
    del_filt += $16
    del_filt_bases += $17
}
END {
    printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
           bases,
           mm_filt + ins_filt + del_filt,
           mm_filt + ins_filt_bases + del_filt_bases,
           mm_filt,
           ins_filt,
           ins_filt_bases,
           del_filt,
           del_filt_bases
}
'
echo "reads with 1 error only:"
zcat $error_file | awk '
BEGIN {
    printf "total_bases\ttotal_filtered_errors\ttotal_filtered_error_bases\tfiltered_mismatches\tfiltered_insertions\tfiltered_insertions_bases\tfiltered_deletions\tfiltered_deletions_bases\n"
}
NR>1 && ($9 + $13 + $17) == 1{
    bases += $7
    mm_filt += $9
    ins_filt += $12
    ins_filt_bases += $13
    del_filt += $16
    del_filt_bases += $17
}
END {
    printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
           bases,
           mm_filt + ins_filt + del_filt,
           mm_filt + ins_filt_bases + del_filt_bases,
           mm_filt,
           ins_filt,
           ins_filt_bases,
           del_filt,
           del_filt_bases
}
'