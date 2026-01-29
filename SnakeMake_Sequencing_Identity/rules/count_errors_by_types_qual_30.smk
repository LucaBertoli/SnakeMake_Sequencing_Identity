#SnakeMake

##calcolo di alcune statistiche partendo dal file ottenuto tsv.gz generato nello step precedente 
rule identity_error_stats_qual_30:
    input:
        stats = "{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz"
    output:
        error_stats = "{output_name}/IDENTITY_{output_name}_reads_30_up_LogOfMeans_{read_type}.tsv.gz.error.stats"
    params:
        script_folder=config["script_folder"]
    wildcard_constraints:
        read_type="read1|read2|all"
    shell:
        r"""
        zcat {input.stats} | awk '
        BEGIN {{
            OFS="\t";
            print "total_bases",
                "total_filtered_errors",
                "total_filtered_error_bases",
                "filtered_mismatches",
                "filtered_insertions",
                "filtered_insertions_bases",
                "filtered_deletions",
                "filtered_deletions_bases"
        }}
        NR>1 {{
            bases += $7
            mm_filt += $9
            ins_filt += $12
            ins_filt_bases += $13
            del_filt += $16
            del_filt_bases += $17
        }}
        END {{
            print bases,
                  mm_filt + ins_filt + del_filt,
                  mm_filt + ins_filt_bases + del_filt_bases,
                  mm_filt,
                  ins_filt,
                  ins_filt_bases,
                  del_filt,
                  del_filt_bases
        }}' > {output.error_stats}
        """