#SnakeMake

##output file tsv.gz contenente i mismatch, indel e valori di identità di ogni read (escludendo secondary, supplementary e unmapped fragments)
rule identity_calculation_error_qual_distribution_qual_Q30:
	input:
		bam="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam",
		vcf=config["vcf"]
	output:
		tsv="{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_error_quality_distribution.tsv.gz"
	threads:
		config["threads"]
	params:
		script_folder=config["script_folder"],
		threads=config["threads"]
	shell:
		"""
		python {params.script_folder}/IdentityRevelations_error_by_base_quality_distribution.py {input.bam} {input.vcf} {output.tsv} {params.threads}
		"""