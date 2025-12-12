#SnakeMake

##prende in input BAM, VCF, fasta e le BED con regioni target
##genera un report calcolando l'identità dei dati allineati con samtools mpileup
##per far questo calcola il numero totale di basi allineate che non sono in posizioni varianti (in vcf) e conta il numero di basi mismatch/indel
rule identity_all_types_reads_qual_30:
	input:
		fasta=config["FASTA"],
		bam="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam",
		vcf=config["vcf"],
		bed=config["BED"]
	output:
		flag="{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.done"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		bash {params.script_folder}/identity_calculation_all_types.sh {input.fasta} {input.bam} {input.vcf} {input.bed} {wildcards.output_name}/IDENTITY_{wildcards.output_name}_reads_30_up_LogOfMeans_{wildcards.read_type}
		touch {output.flag}
		"""

