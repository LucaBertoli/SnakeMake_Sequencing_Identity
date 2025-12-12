#SnakeMake

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

