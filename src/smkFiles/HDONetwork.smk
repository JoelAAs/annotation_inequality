checkpoint add_HDO_annotations_to_bait_prey_publications_network:
    input: 
        bp_network = "work_folder/data/network/raw_networks/bait_prey_publications_network.pkl",
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv"
    output: 
        bp_network_with_attributes = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl"
    script: 
        "../pyScripts/network/add_HDO_annotations_to_raw_network.py"

rule find_bait_count_sums_of_neighbors_for_each_HDO_annotation_at_certain_depth_and_cutoff:
    input: 
        bp_network_with_attributes = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        annotations = "work_folder/data/HDO/cutoff/feature_matrices/single_depth/depth_{depth}_feature_matrix_cutoff_{cutoff}.csv"
    output: 
        neighbors_bait_count_sums_df = "work_folder/data/network/HDO/HDO_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.csv",
        neighbors_bait_count_sums_pickle = "work_folder/data/network/HDO/HDO_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    script: 
        "../pyScripts/network/compute_HDO_neighbors_bait_count_sums.py"

rule find_bait_count_sums_of_neighbors_for_top_HDO_annotation_coefficients_at_certain_depth_and_cutoff:
    input: 
        bp_network_with_attributes = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        neighbors_bait_count_sums_df = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.csv",
        neighbors_bait_count_sums_pickle = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    script: 
        "../pyScripts/network/compute_top_coefficients_HDO_neighbors_bait_count_sums.py"

rule visualize_target_genes_neighborhoods:
    input: 
        observed = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        bp_network_with_attributes = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl"
    output: 
        interactive_network = "work_folder/data/network/HDO/interactive_network_target_genes_depth_{depth}_cutoff_{cutoff}.html"
    script: 
        "../pyScripts/network/visualize_target_genes_neighborhoods.py"

rule plot_bait_count_sums_for_each_HDO_annotation_at_certain_depth_and_cutoff:
    input:
        neighbors_bait_count_sums_pickle = "work_folder/data/network/HDO/HDO_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        plots_pdf = "work_folder/data/network/HDO/plots/depth_{depth}_cutoff_{cutoff}_all_distributions.pdf"
    script:
        "../pyScripts/plotting/plot_HDO_neighbors_bait_count_sums.py"

rule plot_bait_count_sums_of_neighbors_for_top_HDO_annotation_coefficients_at_certain_depth_and_cutoff:
    input:
        neighbors_bait_count_sums_pickle = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        plots_pdf = "work_folder/data/network/HDO/plots/top_coefficients_depth_{depth}_cutoff_{cutoff}.pdf"
    script:
        "../pyScripts/plotting/plot_HDO_neighbors_bait_count_sums.py"

rule build_HDO_baseline_for_top_coefficients_at_certain_depth_and_cutoff:
    input:
        bp_network_with_attributes = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        baseline_bait_count_sums_df = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.csv",
        baseline_bait_count_sums_pickle = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    script:
        "../pyScripts/network/create_top_coefficients_HDO_baseline.py"

rule plot_HDO_baseline_for_top_coefficients_at_certain_depth_and_cutoff:
    input: 
        baseline_bait_count_sums_pickle = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        plots_pdf = "work_folder/data/network/HDO/plots/baseline/baseline_top_coefficients_depth_{depth}_cutoff_{cutoff}.pdf"
    script: 
        "../pyScripts/plotting/plot_baseline_HDO_neighbors_bait_count_sums.py"

rule compute_correlation_between_bait_count_sums_of_annotated_HDO_genes_and_baseline_for_top_coefficients:
    input: 
        observed = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv"
    script:
        "../pyScripts/network/compute_correlation_between_HDO_bait_count_sums_of_annotated_genes_and_baseline.py"

rule visualize_HDO_neighbor_bait_count_sums_summary_stats:
    input:
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        summary_stats_table = "work_folder/data/network/HDO/plots/comparison/top_coefficients_summary_stats_table_depth_{depth}_cutoff_{cutoff}.png",
        fold_enrichment_plot = "work_folder/data/network/HDO/plots/comparison/top_coefficients_fold_enrichment_plot_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/visualize_HDO_summary_stats.py"

rule visualize_correlation_between_bait_count_sums_of_annotated_HDO_genes_and_baseline_for_top_coefficients:
    input:
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        summary_stats_plots = "work_folder/data/network/HDO/plots/comparison/top_coefficients_neighbor_bait_count_sums_comparison_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/plot_correlation_between_HDO_bait_count_sums_of_annotated_genes_and_baseline.py"

rule find_annotated_genes_among_neighbors_for_top_HDO_annotation_coefficients_at_certain_depth_and_cutoff:
    input: 
        bp_network_with_attributes = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        neighbors_with_annotation_df = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.csv",
        neighbors_with_annotation_pickle = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    script: 
        "../pyScripts/network/compute_top_coefficients_HDO_neighbors_annotated_genes.py"

rule plot_annotated_neighbors_for_top_HDO_annotation_coefficients_at_certain_depth_and_cutoff:
    input:
        neighbors_with_annotation_pickle = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        plots_pdf = "work_folder/data/network/HDO/plots/top_coefficients_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pdf"
    script:
        "../pyScripts/plotting/plot_HDO_annotated_neighbors.py"

rule build_baseline_for_annotated_genes_among_neighbors_for_top_HDO_annotation_coefficients_at_certain_depth_and_cutoff:
    input:
        bp_network_with_attributes = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output:
        baseline_neighbors_with_annotation_df = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.csv",
        baseline_neighbors_with_annotation_pickle = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    script:
        "../pyScripts/network/create_HDO_annotated_neighbor_genes_baseline.py"

rule plot_HDO_baseline_for_annotated_coefficients_at_certain_depth_and_cutoff:
    input: 
        baseline_neighbors_with_annotation_pickle = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        plots_pdf = "work_folder/data/network/HDO/plots/baseline/baseline_top_coefficients_annotated_genes_depth_{depth}_cutoff_{cutoff}.pdf"
    script: 
        "../pyScripts/plotting/plot_baseline_HDO_annotated_neighbors.py"

rule compute_correlation_between_annotated_neighbors_of_annotated_HDO_genes_and_baseline_for_top_coefficients:
    input: 
        observed = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv"
    script:
        "../pyScripts/network/compute_correlation_between_HDO_annotated_neighbors_of_annotated_genes_and_baseline.py"

rule visualize_HDO_annotated_neighbors_summary_stats:
    input:
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        summary_stats_table = "work_folder/data/network/HDO/plots/comparison/top_coefficients_annotated_neighbors_summary_stats_table_depth_{depth}_cutoff_{cutoff}.png",
        fold_enrichment_plot = "work_folder/data/network/HDO/plots/comparison/top_coefficients_annotated_neighbors_fold_enrichment_plot_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/visualize_HDO_annotated_neighbors_summary_stats.py"

rule visualize_correlation_between_annotated_neighbors_of_annotated_HDO_genes_and_baseline_for_top_coefficients:
    input:
        summary_stats = "work_folder/data/network/HDO/comparison/top_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.csv",
        observed = "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        baseline = "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_{depth}_cutoff_{cutoff}.pkl",
        ontology = "work_folder/data/HDO/doid.obo"
    output:
        summary_stats_plots = "work_folder/data/network/HDO/plots/comparison/top_coefficients_annotated_neighbors_comparison_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "../pyScripts/plotting/plot_correlation_between_HDO_annotated_neighbors_of_annotated_genes_and_baseline.py"