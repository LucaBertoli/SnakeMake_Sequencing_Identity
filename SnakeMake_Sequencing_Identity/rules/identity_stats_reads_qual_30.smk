#SnakeMake

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
