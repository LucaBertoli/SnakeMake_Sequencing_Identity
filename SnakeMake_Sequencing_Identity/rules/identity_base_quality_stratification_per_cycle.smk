#SnakeMake

rule IdentityBaseQualityStratification_per_cycle:
	input:
		bam="{output_name}/{output_name}_HCR.bam",
		vcf=config["vcf"],
	output:
		tsv1="{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_{read_type}_cycle.tsv",
		tsv2="{output_name}/IDENTITY_{output_name}_base_quality_stratification_by_quality_{read_type}_bin_cycle.tsv"
	threads:
		config["threads"]
	params:
		script_folder=config["script_folder"],
		threads=config["threads"],
		read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
	wildcard_constraints:
		read_type = "all|read1|read2"
	shell:
		"""
		python -u {params.script_folder}/IdentityBaseQualityStratification_binned_per_cycle_parallel.py {input.bam} {input.vcf} {output.tsv1} {output.tsv2} {params.read_param} {params.threads}
		"""

