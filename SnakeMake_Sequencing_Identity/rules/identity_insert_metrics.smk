#SnakeMake

#the scripts computes per-read insert size filtered identity and average sequencer-assigned Qscore. The output is a tsv.gz file with the following columns:
#READ_ID: read identifier
#ORIENTATION: read orientation (R1, R2, S)
#INSERT_SIZE: insert size of the fragment corresponding to the read
#AVERAGE_QSCORE: average sequencer-assigned Qscore for the read
#IDENTITY: filtered alignment identity of the read 

rule identity_insert_metrics:
    input:
        bam="{output_name}/{output_name}_HCR.bam",
        vcf=config["vcf"]
    output:
        tsv="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}.tsv.gz"
    params:
        script_folder=config["script_folder"],
        read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
    shell:
        """
        python {params.script_folder}/IdentityInsertMetrics.py -b {input.bam} -o {output.tsv} -v {input.vcf}
        """