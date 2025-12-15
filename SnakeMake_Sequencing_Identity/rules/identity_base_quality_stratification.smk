#SnakeMake

rule IdentityBaseQualityStratification:
	input:
		bam="{output_name}/{output_name}_HCR.bam",
		vcf=config["vcf"],
	output:
		tsv="{output_name}/IDENTITY_{output_name}_base_quality_stratification.tsv.gz",
		tsv_stats="{output_name}/IDENTITY_{output_name}_base_quality_stratification_binned.tsv.gz"
	params:
		script_folder=config["script_folder"],
		threads = config["threads_quality_stratification"]
	shell:
		"""
		# python {params.script_folder}/IdentityBaseQualityStratification.py {input.bam} {input.vcf} {output.tsv}
		python {params.script_folder}/IdentityBaseQualityStratification_parallel.py {input.bam} {input.vcf} {output.tsv} {params.threads}
		"""
