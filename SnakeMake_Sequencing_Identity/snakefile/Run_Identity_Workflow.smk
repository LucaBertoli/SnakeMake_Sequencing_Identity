configfile:"config.yaml"

OUTPUT_NAME = config["output_name"]  # assegna a variabile

#separa read1 e read2 se richiesto, altrimenti una sola analisi
READS=[config["read_type"]] # può essere "all", "read1" e "read2"

rule all:
	input:
		expand("{output_name}/{output_name}_HCR.bam", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz.stats", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_HCR.bam_mean_quality_LogOfMeans.tsv.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/{output_name}_{read_type}_phred_hist_subplot.png", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_{read_type}_phred_hist.csv", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_full_{read_type}.done", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz.stats", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam_mean_quality_LogOfMeans.tsv.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/{output_name}_{read_type}_reads_30_up_LogOfMeans_phred_hist_subplot.png", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_{read_type}_reads_30_up_LogOfMeans_phred_hist.csv", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.done", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification.tsv.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification.tsv.gz.stats", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_cycle.tsv", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_bin_cycle.tsv", output_name=OUTPUT_NAME)
		expand("{output_name}/IDENTITY_{output_name}_full_all_insert_size_strat.tsv.gz", output_name=OUTPUT_NAME)
		expand("{output_name}/IDENTITY_{output_name}_full_all_insert_size_strat_10M_per_bin.png", output_name=OUTPUT_NAME)
		expand("{output_name}/IDENTITY_{output_name}_full_all_insert_size_strat_5M_per_bin_R1_R2.png", output_name=OUTPUT_NAME)
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}_error_quality_distribution.tsv.gz", output_name=OUTPUT_NAME, read_type=READS)


#intersection of raw bam against GIAB high confidence regions with samtools
import: "rules/intersect_bam.smk"
#identity calculation with IdentityRevelations.py
import: "rules/identity_calculation.smk"
#statistics on identity calculation results with IdentityRevelations_stats.py
import: "rules/identity_stats.smk"
#compute the average quality of the reads in the bam file using Scatter_BAM_LogOfMeans.py
import: "rules/scatter_bam.smk"
#plotting of identity results with Identity_plot.py
import: "rules/identity_plot.smk"
#identity calculation with samtools mpileup using identity_calculation_all_types.sh
import: "rules/identity_all_types.smk"
#filter reads, selecting only those with average Qscore >=30 using Scatter_BAM_LogOfMeans.py
import: "rules/filter_bam.smk"
#identity calculation on Q30+ reads using IdentityRevelations.py
import: "rules/identity_calculation_reads_qual_30.smk"
#statistics on identity calculation results on Q30+ reads with IdentityRevelations_stats.py
import: "rules/identity_stats_reads_qual_30.smk"
#compute average quality of Q30+ reads in the bam file using Scatter_BAM_LogOfMeans.py
import: "rules/scatter_bam_reads_qual_30.smk"
#plotting of identity results on Q30+ reads with Identity_plot.py
import: "rules/identity_plot_reads_qual_30.smk"
#identity calculation with samtools mpileup on Q30+ reads using identity_calculation_all_types.sh
import: "rules/identity_all_types_reads_qual_30.smk"
#compute identity stratified by base quality using IdentityBaseQualityStratification.py
import: "rules/identity_base_quality_stratification.smk"
#compute statistics on identity stratified by base quality using IdentityBaseQualityStratification_stats.py
import: "rules/identity_base_quality_stratification_stats.smk"
#compute the identity stratified by base quality and cycle using IdentityBaseQualityStratification_binned_per_cycle.py
import: "rules/IdentityBaseQualityStratification_per_cycle.smk"
#compute per-read insert size metrics using IdentityInsertMetrics.py
import: "rules/identity_insert_metrics.smk"
#creates plots for Qscore and identity stratified by insert size using IdentityInsertMetrics_plot_subplot_noQscore.py and IdentityInsertMetrics_plot_subplot_noQscore_R1_R2.py
import: "rules/identity_insert_metrics_plot.smk"
#compute quality distribution of the errors in Q30 reads using IdentityRevelations_error_by_base_quality_distribution.py
import: "rules/identity_error_by_base_quality_distribution_qual_Q30.smk"