#SnakeMake

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
