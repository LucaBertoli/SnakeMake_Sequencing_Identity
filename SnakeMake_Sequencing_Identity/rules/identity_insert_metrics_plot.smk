#SnakeMake

#The scripts creates violin plots with the per-read alignment identity and sequencer Qscore distributions stratified by insert size.
#First, the plot is created for the entire dataset, selecting 10M reads (5M fragments).
#Second, an analogous subplot is created separating Read1 and Read2 (5M reads each)

rule identity_insert_metrics_plot:
    input:
        tsv="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}.tsv.gz"
    output:
        plot="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}_10M_per_bin.png"
    params:
        script_folder=config["script_folder"],
        read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
    wildcard_constraints:
        read_type = "all|read1|read2"
    shell:
        """
        python {params.script_folder}/IdentityInsertMetrics_plot_subplot_noQscore.py -i {input.tsv} -o {output.plot} --y-margin 0.05 --max-reads-per-bin 10000000
        """


rule identity_insert_metrics_plot_R1_R2:
    input:
        tsv="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}.tsv.gz"
    output:
        plot="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}_5M_per_bin_R1_R2.png"
    params:
        script_folder=config["script_folder"],
        read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
    wildcard_constraints:
        read_type = "all|read1|read2"
    shell:
        """
        python {params.script_folder}/IdentityInsertMetrics_plot_subplot_noQscore_R1_R2.py -i {input.tsv} -o {output.plot} --y-margin 0.05 --max-reads-per-bin 10000000
        """

rule identity_insert_metrics_plot_GlobalErrorRate:
    input:
        tsv="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}.tsv.gz"
    output:
        plot="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}_10M_per_bin_GlobalErrorRate.png"
    params:
        script_folder=config["script_folder"],
        read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
    wildcard_constraints:
        read_type = "all|read1|read2"
    shell:
        """
        python {params.script_folder}/IdentityInsertMetrics_plot_subplot_noQscore_GlobalErrorRate.py -i {input.tsv} -o {output.plot} --y-margin 0.05 --max-reads-per-bin 10000000
        """

rule identity_insert_metrics_plot_GlobalErrorRate_R1_R2:
    input:
        tsv="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}.tsv.gz"
    output:
        plot="{output_name}/IDENTITY_{output_name}_insert_size_strat_{read_type}_5M_per_bin_GlobalErrorRate_R1_R2.png"
    params:
        script_folder=config["script_folder"],
        read_param=lambda wildcards: "" if wildcards.read_type == "all" else wildcards.read_type
    wildcard_constraints:
        read_type = "all|read1|read2"
    shell:
        """
        python {params.script_folder}/IdentityInsertMetrics_plot_subplot_noQscore_GlobalErrorRate_R1_R2.py -i {input.tsv} -o {output.plot} --y-margin 0.05 --max-reads-per-bin 5000000
        """