#SnakeMake

#The scripts creates violin plots with the per-read alignment identity and sequencer Qscore distributions stratified by insert size.
#First, the plot is created for the entire dataset, selecting 10M reads (5M fragments).
#Second, an analogous subplot is created separating Read1 and Read2 (5M reads each)

rule identity_insert_metrics:
    input:
        tsv="{output_name}/IDENTITY_{output_name}_full_all_insert_size_strat.tsv.gz"
    output:
        plot="{output_name}/IDENTITY_{output_name}_full_all_insert_size_strat_10M_per_bin.png"
    params:
		script_folder=config["script_folder"],
    shell:
        """
        python {params.script_folder}/IdentityInsertMetrics_plot_subplot_noQscore.py -i {input.tsv} -o {output.plot} --y-margin 0.05 --max-reads-per-bin 10000000
        """


rule identity_insert_metrics_R1_R2:
    input:
        tsv="{output_name}/IDENTITY_{output_name}_full_all_insert_size_strat.tsv.gz"
    output:
        plot="{output_name}/IDENTITY_{output_name}_full_all_insert_size_strat_5M_per_bin_R1_R2.png"
    params:
		script_folder=config["script_folder"],
    shell:
        """
        python {params.script_folder}/IdentityInsertMetrics_plot_subplot_noQscore_R1_R2.py -i {input.tsv} -o {output.plot} --y-margin 0.05 --max-reads-per-bin 10000000
        """