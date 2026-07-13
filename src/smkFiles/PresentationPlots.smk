rule plot_go_neighbor_sums:
    input:
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        summary_stats_plots = "work_folder/data/presentation_plots/go_{aspect}_depth_{depth}_cutoff_{cutoff}_neighbor_sums.png"
    script:
        "../pyScripts/presentation_plots/go_bait_count_sums_distribution.py"

rule plot_go_annotated_neighbors:
    input:
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        summary_stats_plots = "work_folder/data/presentation_plots/go_{aspect}_depth_{depth}_cutoff_{cutoff}_annotated_neighbors.png"
    script:
        "../pyScripts/presentation_plots/go_annotated_neighbors_distribution.py"

rule plot_hdo_annotated_neighbors:
    input:
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        summary_stats_plots = "work_folder/data/presentation_plots/hdo_depth_{depth}_cutoff_{cutoff}_annotated_neighbors.png"
    script:
        "../pyScripts/presentation_plots/hdo_annotated_neighbors.py"

rule plot_hdo_neighbor_sums:
    input:
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        summary_stats_plots = "work_folder/data/presentation_plots/hdo_depth_{depth}_cutoff_{cutoff}_neighbor_sums.png"
    script:
        "../pyScripts/presentation_plots/hdo_bait_count_sums.py"

rule create_GO_scatterplot:
    input:
        wait = "work_folder/data/GO/cutoff/done_files/complete_en_coefficients_done.txt",
        all_coefficients = "work_folder/data/dendrograms/GO/all_coefficients/{aspect}_all_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        scatterplot_png = "work_folder/data/dendrograms/GO/visualization/scatterplot/{aspect}_scatterplot_cutoff_{cutoff}.png"
    script: 
        "../pyScripts/presentation_plots/create_GO_scatterplot.py"

rule create_HDO_scatterplot:
    input:
        wait = "work_folder/data/HDO/cutoff/done_files/complete_en_coefficients_done.txt",
        all_coefficients = "work_folder/data/dendrograms/HDO/all_coefficients/all_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output: 
        scatterplot_png = "work_folder/data/dendrograms/HDO/visualization/scatterplot/scatterplot_cutoff_{cutoff}.png"
    script: 
        "../pyScripts/presentation_plots/create_HDO_scatterplot.py"