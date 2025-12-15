#SnakeMake

##output file tsv.gz contenente i mismatch, indel e valori di identità di ogni read (escludendo secondary, supplementary e unmapped fragments)
rule identity_calculation:
    input:
        bam="{output_name}/{output_name}_HCR.bam",
        vcf=config["vcf"]
    output:
        tsv="{output_name}/IDENTITY_{output_name}_full_{read_type}.tsv.gz"
    params:
        script_folder=config["script_folder"],
        read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
    wildcard_constraints:
        read_type = "all|read1|read2"
    shell:
        """
        python {params.script_folder}/IdentityRevelations.py {input.bam} {input.vcf} {output.tsv} {params.read_param}
        """

