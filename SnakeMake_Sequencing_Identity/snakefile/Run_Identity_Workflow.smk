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

##intersezione con bedtools intersect del BAM con un BED file di regioni ad alta confidenza
rule intersect_bam:
	input:
		bam_temporaneo=config["bam_path"],
		bed=config["benchmark_bed"]
	output:
		bam="{output_name}/{output_name}_HCR.bam",
		bai="{output_name}/{output_name}_HCR.bam.bai"
	shell:
		"""
		bedtools intersect -abam {input.bam_temporaneo} -b {input.bed} -f 1 -u > {output.bam}
		samtools index {output.bam}
		"""
##output file tsv.gz contenente i mismatch, indel e valori di identità di ogni read (escludendo secondary, supplementary e unmapped fragments)
rule identity_calculation:
	input:
		bam="{output_name}/{output_name}_HCR.bam",
		vcf=config["vcf"]
	output:
		tsv="{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz"
	params:
		script_folder=config["script_folder"],
		read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
	shell:
		"""
		python  {params.script_folder}/IdentityRevelations.py {input.bam} {input.vcf} {output.tsv} {params.read_param}
		"""

##calcolo di alcune statistiche partendo dal file ottenuto tsv.gz generato nello step precedente 
rule identity_stats:
	input:
		tsv = "{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz"
	output:
		stats = "{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz.stats"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/IdentityRevelations_stats.py {input.tsv} > {output.stats}
		"""

#calcola le qualità medie delle read in un file BAM.
#Le qualità considerate sono le qualità in phread score delle basi allineate (esclude softclipped)
rule scatter_bam:
	input:
		bam="{output_name}/{output_name}_HCR.bam"
	output:
		bam_qual="{output_name}/{output_name}_HCR.bam_mean_quality_LogOfMeans.tsv.gz"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/Scatter_BAM_LogOfMeans.py {input.bam} QUAL
		"""

##lo script prende in input il tsv.gz con le identità (generato nella rule identity_calculation) e il tsv.gz con i phread score del BAM (generato da scatter_bam)
##converte le identità in Q-score con la formula Q-score=-10*log(1-identity)
##infine viene creato un plot con i Q-score assegnati dal sequenziatore e a quelli calcolati sulla base dell'allineamento
rule identity_plot:
	input:
		tsv="{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz",
		qual="{output_name}/{output_name}_HCR.bam_mean_quality_LogOfMeans.tsv.gz"
	output:
		png="{output_name}/{output_name}_{read_type}_phred_hist_subplot.png",
		csv="{output_name}/{output_name}_{read_type}_phred_hist.csv"
	wildcard_constraints:
		read_type="read1|read2|all"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/Identity_plot.py {input.tsv} {input.qual} {wildcards.output_name}/{wildcards.output_name}_{wildcards.read_type}_phred_hist
		"""

##prende in input BAM, VCF, fasta e le BED con regioni target
##genera un report calcolando l'identità dei dati allineati con samttols mpileup 
##per far questo calcola il numero totale di basi allineate che non sono in posizioni varianti (in vcf) e conta il numero di basi mismatch/indel
rule identity_all_types:
	input:
		fasta=config["FASTA"],
		bam="{output_name}/{output_name}_HCR.bam",
		vcf=config["vcf"],
		bed=config["BED"]
	output:
		flag="{output_name}/IDENTITY_{output_name}_full_{read_type}.done"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		bash {params.script_folder}/identity_calculation_all_types.sh {input.fasta} {input.bam} {input.vcf} {input.bed} {wildcards.output_name}/IDENTITY_{wildcards.output_name}_full_{wildcards.read_type}
		touch {output.flag}
		"""

######
#filtraggio delle reads, si selezionano quelle con qualità >=30
rule filter_bam:
	input:
		bam="{output_name}/{output_name}_HCR.bam"
	output:
		bam_filtered="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/Scatter_BAM_LogOfMeans.py {input.bam} BAM
		samtools index {output.bam_filtered}
		"""

##output file tsv.gz contenente i mismatch, indel e valori di identità di ogni read (escludendo secondary, supplementary e unmapped fragments)
rule identity_calculation_reads_qual_30:
	input:
		bam="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam",
		vcf=config["vcf"]
	output:
		tsv="{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz"
	params:
		script_folder=config["script_folder"],
		read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
	shell:
		"""
		python  {params.script_folder}/IdentityRevelations.py {input.bam} {input.vcf} {output.tsv} {params.read_param}
		"""
##calcolo di alcune statistiche partendo dal file ottenuto tsv.gz generato nello step precedente
rule identity_stats_reads_qual_30:
	input:
		tsv = "{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz"
	output:
		stats = "{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz.stats"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/IdentityRevelations_stats.py {input.tsv} > {output.stats}
		"""
##calcola le qualità medie delle read in un file BAM.
##Le qualità considerate sono le qualità in phread score delle basi allineate (esclude softclipped)
rule scatter_bam_reads_qual_30:
	input:
		bam="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam"
	output:
		bam_qual="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam_mean_quality_LogOfMeans.tsv.gz"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/Scatter_BAM_LogOfMeans.py {input.bam} QUAL
		"""
##lo script prende in input il tsv.gz con le identità (generato nella rule identity_calculation) e il tsv.gz con i phread score del BAM (generato da scatter_bam)
##converte le identità in Q-score con la formula Q-score=-10*log(1-identity)
##infine viene creato un plot con i Q-score assegnati dal sequenziatore e a quelli calcolati sulla base dell'allineamento
rule identity_plot_reads_qual_30:
	input:
		tsv="{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz",
		qual="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam_mean_quality_LogOfMeans.tsv.gz"
	output:
		png="{output_name}/{output_name}_{read_type}_reads_30_up_LogOfMeans_phred_hist_subplot.png",
		csv="{output_name}/{output_name}_{read_type}_reads_30_up_LogOfMeans_phred_hist.csv"
	wildcard_constraints:
		read_type="read1|read2|all"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/Identity_plot.py {input.tsv} {input.qual} {wildcards.output_name}/{wildcards.output_name}_{wildcards.read_type}_reads_30_up_LogOfMeans_phred_hist
		"""
##prende in input BAM, VCF, fasta e le BED con regioni target
##genera un report calcolando l'identità dei dati allineati con samttols mpileup
##per far questo calcola il numero totale di basi allineate che non sono in posizioni varianti (in vcf) e conta il numero di basi mismatch/indel
rule identity_all_types_reads_qual_30:
	input:
		fasta=config["FASTA"],
		bam="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam",
		vcf=config["vcf"],
		bed=config["BED"]
	output:
		flag="{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.done"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		bash {params.script_folder}/identity_calculation_all_types.sh {input.fasta} {input.bam} {input.vcf} {input.bed} {wildcards.output_name}/IDENTITY_{wildcards.output_name}_reads_30_up_LogOfMeans_{wildcards.read_type}
		touch {output.flag}
		"""


rule IdentityBaseQualityStratification:
	input:
		bam="{output_name}/{output_name}_HCR.bam",
		vcf=config["vcf"],
	output:
		tsv="{output_name}/IDENTITY_{output_name}_base_quality_stratification.tsv.gz",
		tsv_stats="{output_name}/IDENTITY_{output_name}_base_quality_stratification.tsv.gz.stats"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python {params.script_folder}/IdentityBaseQualityStratification.py {input.bam} {input.vcf} {output.tsv}
		"""

rule IdentityBaseQualityStratification_per_cycle:
	input:
		bam="{output_name}/{output_name}_HCR.bam",
		vcf=config["vcf"],
	output:
		tsv1="{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_cycle.tsv",
		tsv2="{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_bin_cycle.tsv"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		python -u {params.script_folder}/IdentityBaseQualityStratification_binned_per_cycle.py {input.bam} {input.vcf} {output.tsv1} {output.tsv2}
		"""
