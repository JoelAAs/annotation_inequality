ASPECTS = ['BP', 'MF', 'CC']

checkpoint add_GO_annotations_to_bait_prey_publications_network:
    input: 
        bp_network = "work_folder/data/network/raw_networks/bait_prey_publications_network.pkl",
        complete_annotations = "work_folder/data/GO/{aspect}_annotations_per_gene.csv"
    output: 
        bp_network_with_attributes = "work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl"
    script: 
        "../pyScripts/network/add_GO_annotations_to_raw_network.py"

rule find_bait_count_sums_of_neighbors_for_each_GO_annotation_at_certain_depth_and_cutoff:
    input: 
        bp_network_with_attributes = "work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl",
        annotations = "work_folder/data/GO/cutoff/feature_matrices/single_depth/{aspect}_feature_matrix_depth_{depth}_cutoff_{cutoff}.csv"
    output: 
        neighbors_bait_count_sums_df = "work_folder/data/network/GO/all_annotations/{aspect}_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.csv",
        neighbors_bait_count_sums_pickle = "work_folder/data/network/GO/all_annotations/{aspect}_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    script: 
        "../pyScripts/network/compute_GO_neighbors_bait_count_sums.py"

rule find_bait_count_sums_of_neighbors_for_top_GO_annotation_coefficients_at_certain_depth_and_cutoff:
    input: 
        bp_network_with_attributes = "work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        neighbors_bait_count_sums_df = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.csv",
        neighbors_bait_count_sums_pickle = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    script: 
        "../pyScripts/network/compute_top_coefficients_GO_neighbors_bait_count_sums.py"

rule plot_bait_count_sums_for_each_GO_annotation_at_certain_depth_and_cutoff:
    input:
        neighbors_bait_count_sums_pickle = "work_folder/data/network/GO/all_annotations/{aspect}_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        plots_pdf = "work_folder/data/network/GO/plots/all_annotations/{aspect}_depth_{depth}_cutoff_{cutoff}_all_distributions.pdf"
    script:
        "../pyScripts/plotting/plot_GO_neighbors_bait_count_sums.py"

rule plot_bait_count_sums_of_neighbors_for_top_GO_annotation_coefficients_at_certain_depth_and_cutoff:
    input:
        neighbors_bait_count_sums_pickle = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        plots_pdf = "work_folder/data/network/GO/plots/top_coefficients/{aspect}_top_coefficients_depth_{depth}_cutoff_{cutoff}_all_distributions.pdf"
    script:
        "../pyScripts/plotting/plot_GO_neighbors_bait_count_sums.py"

rule build_GO_baseline_for_top_coefficients_at_certain_depth_and_cutoff:
    input:
        bp_network_with_attributes = "work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        baseline_bait_count_sums_df = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.csv",
        baseline_bait_count_sums_pickle = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    script:
        "../pyScripts/network/create_top_coefficients_GO_baseline.py"

rule plot_GO_baseline_for_top_coefficients_at_certain_depth_and_cutoff:
    input: 
        baseline_bait_count_sums_pickle = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        plots_pdf = "work_folder/data/network/GO/plots/baseline/top_coefficients/baseline_{aspect}_top_coefficients_depth_{depth}_cutoff_{cutoff}_all_distributions.pdf"
    script: 
        "../pyScripts/plotting/plot_baseline_GO_neighbors_bait_count_sums.py"

rule compute_correlation_between_bait_count_sums_of_annotated_GO_genes_and_baseline_for_top_coefficients:
    input: 
        observed = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv"
    script:
        "../pyScripts/network/compute_correlation_between_GO_bait_count_sums_of_annotated_genes_and_baseline.py"

rule visualize_GO_summary_stats:
    input:
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        summary_stats_table = "work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_summary_stats_table_depth_{depth}_cutoff_{cutoff}.png",
        fold_enrichment_plot = "work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_fold_enrichment_plot_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/visualize_GO_summary_stats.py"

rule visualize_correlation_between_bait_count_sums_of_annotated_GO_genes_and_baseline_for_top_coefficients:
    input:
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        summary_stats_plots = "work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/plot_correlation_between_GO_bait_count_sums_of_annotated_genes_and_baseline.py"

rule find_annotated_genes_among_neighbors_for_top_GO_annotation_coefficients_at_certain_depth_and_cutoff:
    input: 
        bp_network_with_attributes = "work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        neighbors_with_annotation_df = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.csv",
        neighbors_with_annotation_pickle = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    script: 
        "../pyScripts/network/compute_top_coefficients_GO_neighbors_annotated_genes.py"

rule plot_annotated_neighbors_for_top_GO_annotation_coefficients_at_certain_depth_and_cutoff:
    input:
        neighbors_with_annotation_pickle = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        plots_pdf = "work_folder/data/network/GO/plots/top_coefficients/{aspect}_top_coefficients_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pdf"
    script:
        "../pyScripts/plotting/plot_GO_annotated_neighbors.py"

rule build_baseline_for_annotated_genes_among_neighbors_for_top_GO_annotation_coefficients_at_certain_depth_and_cutoff:
    input:
        bp_network_with_attributes = "work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output:
        baseline_neighbors_with_annotation_df = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.csv",
        baseline_neighbors_with_annotation_pickle = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    script:
        "../pyScripts/network/create_GO_annotated_neighbor_genes_baseline.py"

rule plot_GO_baseline_for_annotated_coefficients_at_certain_depth_and_cutoff:
    input: 
        baseline_neighbors_with_annotation_pickle = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        plots_pdf = "work_folder/data/network/GO/plots/baseline/top_coefficients/baseline_{aspect}_top_coefficients_annotated_genes_depth_{depth}_cutoff_{cutoff}.pdf"
    script: 
        "../pyScripts/plotting/plot_baseline_GO_annotated_neighbors.py"

rule compute_correlation_between_annotated_neighbors_of_annotated_GO_genes_and_baseline_for_top_coefficients:
    input: 
        observed = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv"
    script:
        "../pyScripts/network/compute_correlation_between_GO_annotated_neighbors_of_annotated_genes_and_baseline.py"

rule visualize_GO_annotated_neighbors_summary_stats:
    input:
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        summary_stats_table = "work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_summary_stats_table_depth_{depth}_cutoff_{cutoff}.png",
        fold_enrichment_plot = "work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_fold_enrichment_plot_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/visualize_GO_annotated_neighbors_summary_stats.py"

rule visualize_correlation_between_annotated_neighbors_of_annotated_GO_genes_and_baseline_for_top_coefficients:
    input:
        summary_stats = "work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        summary_stats_plots = "work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/plot_correlation_between_GO_annotated_neighbors_of_annotated_genes_and_baseline.py"