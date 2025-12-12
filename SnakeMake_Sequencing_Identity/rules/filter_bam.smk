#SnakeMake

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
