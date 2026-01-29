#SnakeMake

##prende in input BAM, VCF, fasta e le BED con regioni target
##genera un report calcolando l'identità dei dati allineati con samtools mpileup 
##per far questo calcola il numero totale di basi allineate che non sono in posizioni varianti (in vcf) e conta il numero di basi mismatch/indel
rule mpileup_qual_30:
	input:
		bam="{output_name}/{output_name}_HCR_reads_30_up_LogOfMeans.bam",
		vcf=config["vcf"],
		bed=config["BED"],
		fasta=config["FASTA"]
	output:
		pileup=temp("{output_name}/{output_name}_reads_30_up_LogOfMeans.all.types.pileup.gz"),
		new_target=temp("{output_name}/{output_name}_reads_30_up_LogOfMeans.all.types.pileup.target.bed")
	shell:
		"""
		bedtools subtract -a {input.bed} -b {input.vcf} > {output.new_target}
		samtools mpileup -A -l {output.new_target} -f {input.fasta} -q 0 -Q 0 {input.bam} | bgzip > {output.pileup}
		"""


rule identity_all_types_qual_30:
	input:
		pileup="{output_name}/{output_name}_reads_30_up_LogOfMeans.all.types.pileup.gz"
	output:
		out_report="{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_full_all.identity.all.types.pileup.report"
	params:
		script_folder=config["script_folder"]
	shell:
		"""
		# Calcolo delle basi totali
		TOTAL_BASES=$(zcat {input.pileup} | awk -F "\t" '{{SUM+=$4}} END {{print SUM}}')
		echo "Total Bases: $TOTAL_BASES" > {output.out_report}
		# Calcolo delle basi mismatch e indel
		TOTAL_MM=$(zcat {input.pileup} | awk -F "\t" '$3!="N" {{print $5}}' | grep -oe A -oe a -oe C -oe c -oe G -oe g -oe T -oe t | wc -l)
		echo "TOTAL mismatches: $TOTAL_MM" >> {output.out_report}
		# Calcolo dell'error rate
		MISMATCH_RATE=$(awk "BEGIN {{print ($TOTAL_MM / $TOTAL_BASES)*100}}")
		echo "Mismatch rate: $MISMATCH_RATE%" >> {output.out_report}
		# Calcolo dell'identità
		IDENTITY=$(awk "BEGIN {{print 100 - $MISMATCH_RATE}}")
		echo "Identity: $IDENTITY%" >> {output.out_report}
		"""

