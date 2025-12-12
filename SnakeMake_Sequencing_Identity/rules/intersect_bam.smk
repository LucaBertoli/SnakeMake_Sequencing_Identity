#SnakeMake

##intersezione con bedtools intersect del BAM con un BED file di regioni ad alta confidenza
rule intersect_bam:
	input:
		bam_temporaneo=config["bam_path"],
		bed=config["benchmark_bed"]
	output:
		bam="{output_name}/{output_name}_HCR.bam",
		bai="{output_name}/{output_name}_HCR.bam.bai"
	shell:
		"""
		bedtools intersect -abam {input.bam_temporaneo} -b {input.bed} -f 1 -u > {output.bam}
		samtools index {output.bam}
		"""
