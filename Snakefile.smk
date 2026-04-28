include: "src/BaitUsage.smk"
include: "src/GetAnnotationData.smk"
include: "src/DfJoin.smk"
include: "src/PrintCorrelationValuesPlots.smk"
include: "src/ComputeCorrelationValues.smk"
include: "src/FeatureMatrix.smk"
include: "src/AddAncestors.smk"
include: "src/ElasticNet.smk"
include: "src/AnnotationDepth.smk"
include: "src/AnnotVSGeneCount.smk"
include: "src/HDOCutoffs.smk"
include: "src/GOCutoffs.smk"
include: "src/DendrogramOfCoefficients.smk"
include: "src/AdjustedRSquared.smk"
include: "src/CoefficientsPlotting.smk"
include: "src/smkFiles/RawNetwork.smk"
include: "src/smkFiles/HDONetwork.smk"
include: "src/smkFiles/GONetwork.smk"
include: "src/smkFiles/GOAnnotationDates.smk"

ASPECTS = ["BP", "CC", "MF"]

rule all:
    input:
        # --- INTACT SECTION ---
        "work_folder/data/intact/bait_prey_publications_no_gene_complete_dates.pq",
        "work_folder/data/intact/bait_prey_publications_complete_dates.pq",

        # --- ANNOTATIONS VS GENE COUNTS ---
        expand("work_folder/data/plots/GO_plots/{aspect}_annotations_vs_gene_counts_distrib.png",
               aspect = ASPECTS),
        expand("work_folder/data/GO/annotations_vs_gene_counts/{aspect}_annotations_gene_counts.csv",
               aspect = ASPECTS),

        # --- GO COMPLETE FEATURE MATRIX SECTION ---
        expand("work_folder/data/GO/{aspect}_annotations_per_gene.csv", aspect=ASPECTS),
        expand("work_folder/data/GO/{aspect}_all_annotations.csv", aspect=ASPECTS),
        "work_folder/data/GO/max_depths_file.csv",
        expand("work_folder/data/GO/annotations_counts/{aspect}_annotations_counts.csv", aspect = ASPECTS),
        expand("work_folder/data/GO/feature_matrices/complete/{aspect}_complete_feature_matrix.csv", aspect = ASPECTS),

        # --- GO SINGLE DEPTH FEATURE MATRIX SECTION --- 
        get_all_single_depth_matrices,      

        # --- GO COMPLETE EN SECTION ---
        expand("work_folder/data/ElasticNet/GO/EN_coefficients/complete/{aspect}_complete_elastic_net_coefficients.csv", aspect = ASPECTS),
        expand("work_folder/data/ElasticNet/GO/plots/Top/complete/complete_{aspect}_top_coefficients.png", aspect = ASPECTS),
        expand("work_folder/data/ElasticNet/GO/plots/Distribution/complete/complete_{aspect}_coefficients_distribution.png", aspect = ASPECTS),
        expand("work_folder/data/GO/annotations_per_depth_{aspect}.csv", aspect = ASPECTS),
        expand("work_folder/data/GO/genes_per_depth_{aspect}.csv", aspect = ASPECTS),
        expand("work_folder/data/ElasticNet/GO/plots/{aspect}_annotations_per_depth.png", aspect = ASPECTS),
        expand("work_folder/data/ElasticNet/GO/plots/{aspect}_genes_per_depth.png", aspect = ASPECTS),

        # --- GO SINGLE DEPTH EN SECTION ---
        get_all_en_coefficients,
        get_all_en_plots,

        # --- COMPLETE GO CUTOFF SECTION ---
        expand("work_folder/data/GO/cutoff/feature_matrices/complete/complete_{aspect}_feature_matrix_cutoff_{cutoff}.csv",
               aspect = ASPECTS, cutoff = CUTOFFS),
        expand("work_folder/data/GO/cutoff/cutoff_files/complete/complete_{aspect}_cutoff_{cutoff}_file.csv",
               aspect = ASPECTS, cutoff = CUTOFFS),
        expand("work_folder/data/GO/cutoff/cutoff_files/complete/complete_{aspect}_cutoff_file.csv", 
               aspect = ASPECTS),
        "work_folder/data/GO/cutoff/done_files/complete_matrices_done.txt",
        expand("work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/complete/complete_{aspect}_elastic_net_coefficients_cutoff_{cutoff}.csv",
               aspect = ASPECTS, cutoff = CUTOFFS),
        "work_folder/data/GO/cutoff/done_files/complete_en_coefficients_done.txt",
        expand("work_folder/data/ElasticNet/GO_cutoff/plots/Top/complete/complete_{aspect}_top_coefficients_cutoff_{cutoff}.png", 
               aspect = ASPECTS, cutoff = CUTOFFS),
        expand("work_folder/data/ElasticNet/GO_cutoff/plots/Distribution/complete/complete_{aspect}_coefficients_distribution_cutoff_{cutoff}.png",
               aspect = ASPECTS, cutoff = CUTOFFS),
        expand("work_folder/data/ElasticNet/GO_cutoff/plots/go_ids_lost/complete/complete_{aspect}_GO_ids_lost.png",
               aspect = ASPECTS),

        # --- SINGLE DEPTH GO CUTOFF SECTION ---
        get_all_GO_depth_cutoff_paths,
        get_all_GO_depth_summary_files,
        "work_folder/data/GO/cutoff/done_files/single_depth_matrices_done.txt",
        get_all_GO_single_depth_en_coefficients,
        "work_folder/data/GO/cutoff/done_files/single_depth_en_coefficients_done.txt",
        get_all_GO_single_depth_adj_r2_files,
        get_all_GO_single_depth_en_plots,
        get_all_GO_lost_ids_plots,
        get_all_GO_top_abs_plots,
        expand("work_folder/data/GO/cutoff/adj_r2_plots/{aspect}_cutoff_{cutoff}_adj_r2.png",
               aspect = ASPECTS, cutoff = CUTOFFS),

        # --- COMPLETE HDO CUTOFF SECTION ---
        expand("work_folder/data/HDO/cutoff/feature_matrices/complete/complete_feature_matrix_cutoff_{cutoff}.csv", cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete/complete_elastic_net_coefficients_cutoff_{cutoff}.csv", cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Top/complete/complete_top_coefficients_cutoff_{cutoff}.png", cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/complete/complete_coefficients_distribution_cutoff_{cutoff}.png", cutoff=CUTOFFS),
        "work_folder/data/ElasticNet/HDO_cutoff/plots/doids_lost/complete_HDO_doids_lost.png",

        # --- SINGLE DEPTH HDO CUTOFF SECTION ---
        expand("work_folder/data/HDO/cutoff/feature_matrices/single_depth/depth_{depth}_feature_matrix_cutoff_{cutoff}.csv", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Top/single_depth/depth_{depth}_top_coefficients_cutoff_{cutoff}.png", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/single_depth/depth_{depth}_coefficients_distribution_cutoff_{cutoff}.png", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/doids_lost/depth_{depth}_HDO_doids_lost.png", 
               depth=HDO_DEPTHS_WITH_ANCESTORS),
        expand("work_folder/data/HDO/cutoff/adj_r2_files/full/cutoff_{cutoff}_adj_r2_file.csv",
               cutoff = CUTOFFS),
        expand("work_folder/data/HDO/cutoff/adj_r2_plots/cutoff_{cutoff}_adj_r2.png",
               cutoff = CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Top_abs/top_abs_depth_{depth}_cutoff_{cutoff}.png",
               depth = HDO_DEPTHS_WITH_ANCESTORS, cutoff = CUTOFFS),

        # --- HDO DENDROGRAM OF COEFFICIENTS SECTION ---
        expand("work_folder/data/dendrograms/HDO/all_coefficients/all_coefficients_cutoff_{cutoff}.csv",
               cutoff = CUTOFFS),
        expand("work_folder/data/dendrograms/HDO/visualization/dendrogram/dendrogram_cutoff_{cutoff}.pdf",
               cutoff = CUTOFFS),
        expand("work_folder/data/dendrograms/HDO/visualization/treemap/treemap_cutoff_{cutoff}.pdf",
               cutoff = CUTOFFS),

        # --- GO DENDROGRAM OF COEFFICIENTS SECTION ---
        expand("work_folder/data/dendrograms/GO/all_coefficients/{aspect}_all_coefficients_cutoff_{cutoff}.csv",
               aspect = ASPECTS, cutoff = CUTOFFS),
        expand("work_folder/data/dendrograms/GO/visualization/dendrogram/{aspect}_dendrogram_cutoff_{cutoff}.pdf",
               aspect = ASPECTS, cutoff = CUTOFFS),
        expand("work_folder/data/dendrograms/GO/visualization/treemap/{aspect}_treemap_cutoff_{cutoff}.html",
               aspect = ASPECTS, cutoff = CUTOFFS),

        # --- RAW NETWORK CREATION SECTION ---
         "work_folder/data/network/raw_networks/raw_network_degree_frequencies.csv",
         "work_folder/data/network/raw_networks/bait_prey_publications_network.pkl",
         "work_folder/data/network/raw_networks/raw_network_degree_frequencies.png",
         "work_folder/data/network/raw_networks/network_genes.txt",
         "work_folder/data/network/raw_networks/network_genes_translated.txt",

        # --- HDO ANNOTATIONS NETWORK SECTION ---
         "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
         "work_folder/data/network/HDO/HDO_neighbors_bait_count_sums_depth_5_cutoff_20.csv",
         "work_folder/data/network/HDO/HDO_neighbors_bait_count_sums_depth_5_cutoff_20.pkl",
         "work_folder/data/network/HDO/plots/depth_5_cutoff_20_all_distributions.pdf",
         "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.csv",
         "work_folder/data/network/HDO/HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.pkl",
         "work_folder/data/network/HDO/plots/top_coefficients_depth_5_cutoff_20.pdf",
         "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.csv",
         "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.pkl",
         "work_folder/data/network/HDO/plots/baseline/baseline_top_coefficients_depth_5_cutoff_20.pdf",
         "work_folder/data/network/HDO/comparison/top_coefficients_neighbor_bait_count_sums_comparison_depth_5_cutoff_20.csv",
         "work_folder/data/network/HDO/plots/comparison/top_coefficients_neighbor_bait_count_sums_comparison_depth_5_cutoff_20.png",
         "work_folder/data/network/HDO/plots/comparison/top_coefficients_summary_stats_table_depth_5_cutoff_20.png",
         "work_folder/data/network/HDO/plots/comparison/top_coefficients_fold_enrichment_plot_depth_5_cutoff_20.png",
         "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.csv",
         "work_folder/data/network/HDO/HDO_top_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.pkl",
         "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.csv",
         "work_folder/data/network/HDO/baseline/baseline_HDO_top_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.pkl",
         "work_folder/data/network/HDO/plots/top_coefficients_annotated_neighbors_depth_5_cutoff_20.pdf",
         "work_folder/data/network/HDO/plots/baseline/baseline_top_coefficients_annotated_genes_depth_5_cutoff_20.pdf",
         "work_folder/data/network/HDO/comparison/top_coefficients_annotated_neighbors_comparison_depth_5_cutoff_20.csv",
         "work_folder/data/network/HDO/plots/comparison/top_coefficients_annotated_neighbors_summary_stats_table_depth_5_cutoff_20.png",
         "work_folder/data/network/HDO/plots/comparison/top_coefficients_annotated_neighbors_fold_enrichment_plot_depth_5_cutoff_20.png",
         "work_folder/data/network/HDO/plots/comparison/top_coefficients_annotated_neighbors_comparison_depth_5_cutoff_20.png",

        # --- GO ANNOTATIONS NETWORK SECTION ---
        expand("work_folder/data/network/GO/{aspect}_bait_prey_publications_network.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/all_annotations/{aspect}_neighbors_bait_count_sums_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/all_annotations/{aspect}_neighbors_bait_count_sums_depth_5_cutoff_20.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/all_annotations/{aspect}_depth_5_cutoff_20_all_distributions.pdf",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/top_coefficients/{aspect}_top_coefficients_depth_5_cutoff_20_all_distributions.pdf",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_neighbors_bait_count_sums_depth_5_cutoff_20.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/baseline/top_coefficients/baseline_{aspect}_top_coefficients_depth_5_cutoff_20_all_distributions.pdf",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_neighbor_bait_count_sums_comparison_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_summary_stats_table_depth_5_cutoff_20.png",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_fold_enrichment_plot_depth_5_cutoff_20.png",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_neighbor_bait_count_sums_comparison_depth_5_cutoff_20.png",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/top_coefficients/top_{aspect}_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/baseline/top_coefficients/baseline_{aspect}_top_coefficients_nodes_annotated_neighbors_depth_5_cutoff_20.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/baseline/top_coefficients/baseline_{aspect}_top_coefficients_annotated_genes_depth_5_cutoff_20.pdf",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_comparison_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_summary_stats_table_depth_5_cutoff_20.png",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_fold_enrichment_plot_depth_5_cutoff_20.png",
               aspect = ASPECTS),
        expand("work_folder/data/network/GO/plots/comparison/top_coefficients/top_{aspect}_coefficients_annotated_neighbors_comparison_depth_5_cutoff_20.png",
               aspect = ASPECTS),

        # --- GO ANNOTATIONS DATES SECTION ---
        "work_folder/data/dates/GO/downloads_complete.txt",
        get_all_parsed_snapshots,
        "work_folder/data/dates/GO/GO_first_annotation_dates.csv",
        expand("work_folder/data/dates/GO/networks_with_dates/{aspect}_network_with_dates.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/networks_with_dates/{aspect}_network_with_dates_complete.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_5_cutoff_20.csv",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_5_cutoff_20.pkl",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/network_statistics/{aspect}_depth_5_cutoff_20_networks_statistics.csv",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/plots/{aspect}_depth_5_cutoff_20_total_plots.png",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/plots/{aspect}_depth_5_cutoff_20_annotated_plots.png",
               aspect = ASPECTS),
        expand("work_folder/data/dates/GO/correlation/{aspect}_depth_5_cutoff_20_temporal_correlation_results.csv",
               aspect = ASPECTS),