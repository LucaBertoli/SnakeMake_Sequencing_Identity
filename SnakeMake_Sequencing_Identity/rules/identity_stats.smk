#SnakeMake

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
