configfile:"config.yaml"

OUTPUT_NAME = config["output_name"]  # assegna a variabile

#separa read1 e read2 se richiesto, altrimenti una sola analisi
# READS=[config["read_type"]] # può essere "all", "read1" e "read2"
READS = ["all", "read1", "read2"]

rule all:
	input:
		######## ANALYSIS ON ALL READS ########
		expand("{output_name}/{output_name}_HCR.bam", output_name=OUTPUT_NAME),
		expand("{output_name}/{output_name}_HCR.bam_mean_quality_LogOfMeans.tsv.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz.stats", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_{read_type}_phred_hist_subplot.png", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_{read_type}_phred_hist.csv", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}.all.types.pileup.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/{output_name}.all.types.pileup.target.bed", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_full_all.identity.all.types.pileup.report", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification.tsv.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification_binned.tsv.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_{read_type}_cycle.tsv", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_{read_type}_bin_cycle.tsv", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}.tsv.gz", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}_10M_per_bin.png", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}_5M_per_bin_R1_R2.png", output_name=OUTPUT_NAME, read_type=READS),
		######## ANALYSIS ON READS WITH AVERAGE QSCORE >=30 ########
		expand("{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz.stats", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam_mean_quality_LogOfMeans.tsv.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/{output_name}_{read_type}_reads_30_up_LogOfMeans_phred_hist_subplot.png", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_{read_type}_reads_30_up_LogOfMeans_phred_hist.csv", output_name=OUTPUT_NAME, read_type=READS),
		expand("{output_name}/{output_name}_reads_30_up_LogOfMeans.all.types.pileup.gz", output_name=OUTPUT_NAME),
		expand("{output_name}/{output_name}_reads_30_up_LogOfMeans.all.types.pileup.target.bed", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_full_all.identity.all.types.pileup.report", output_name=OUTPUT_NAME),
		expand("{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_error_quality_distribution.tsv.gz", output_name=OUTPUT_NAME, read_type=READS)

######## ANALYSIS ON ALL READS ########
#intersection of raw bam against GIAB high confidence regions with samtools
include: config["rules_folder"] + "intersect_bam.smk",
#compute the average quality of the reads in the bam file using Scatter_BAM_LogOfMeans.py
include: config["rules_folder"] + "scatter_bam.smk",
#identity calculation with IdentityRevelations.py
include: config["rules_folder"] + "identity_calculation.smk",
#statistics on identity calculation results with IdentityRevelations_stats.py
include: config["rules_folder"] + "identity_stats.smk",
#plotting of identity results with Identity_plot.py
include: config["rules_folder"] + "identity_plot.smk",
#identity calculation with samtools mpileup using identity_calculation_all_types.sh
include: config["rules_folder"] + "identity_all_types.smk",
#compute identity stratified by base quality using IdentityBaseQualityStratification.py
include: config["rules_folder"] + "identity_base_quality_stratification.smk",
#compute the identity stratified by base quality and cycle using IdentityBaseQualityStratification_binned_per_cycle.py
include: config["rules_folder"] + "identity_base_quality_stratification_per_cycle.smk",
#compute per-read insert size metrics using IdentityInsertMetrics.py
include: config["rules_folder"] + "identity_insert_metrics.smk",
#creates plots for Qscore and identity stratified by insert size using IdentityInsertMetrics_plot_subplot_noQscore.py and IdentityInsertMetrics_plot_subplot_noQscore_R1_R2.py
include: config["rules_folder"] + "identity_insert_metrics_plot.smk",

######## ANALYSIS ON READS WITH AVERAGE QSCORE >=30 ########
#filter reads, selecting only those with average Qscore >=30 using Scatter_BAM_LogOfMeans.py
include: config["rules_folder"] + "filter_bam.smk",
#identity calculation on Q30+ reads using IdentityRevelations.py
include: config["rules_folder"] + "identity_calculation_reads_qual_30.smk",
#statistics on identity calculation results on Q30+ reads with IdentityRevelations_stats.py
include: config["rules_folder"] + "identity_stats_reads_qual_30.smk",
#compute average quality of Q30+ reads in the bam file using Scatter_BAM_LogOfMeans.py
include: config["rules_folder"] + "scatter_bam_reads_qual_30.smk",
#plotting of identity results on Q30+ reads with Identity_plot.py
include: config["rules_folder"] + "identity_plot_reads_qual_30.smk",
#identity calculation with samtools mpileup on Q30+ reads using identity_calculation_all_types.sh
include: config["rules_folder"] + "identity_all_types_reads_qual_30.smk",
#compute quality distribution of the errors in Q30 reads using IdentityRevelations_error_by_base_quality_distribution.py
include: config["rules_folder"] + "identity_error_by_base_quality_distribution_qual_Q30.smk",