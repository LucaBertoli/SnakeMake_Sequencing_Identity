#SnakeMake

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
