#SnakeMake

rule IdentityBaseQualityStratification:
	input:
		bam="{output_name}/{output_name}_HCR.bam",
		vcf=config["vcf"],
	output:
		tsv="{output_name}/IDENTITY_{output_name}_base_quality_stratification.tsv.gz",
		tsv_binned="{output_name}/IDENTITY_{output_name}_base_quality_stratification_binned.tsv.gz"
	threads:
		config["threads"]
	params:
		script_folder=config["script_folder"],
		threads = config["threads"]
	shell:
		"""
		python {params.script_folder}/IdentityBaseQualityStratification_parallel.py {input.bam} {input.vcf} {output.tsv} {output.tsv_binned} {params.threads}
		"""
