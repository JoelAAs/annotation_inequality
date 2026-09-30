rule generate_GO_master_temporal_matrices_exact_degree_match:
    input: 
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl",
        top_annot_df = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        matrix_dir = directory("work_folder/data/dates/GO/ed_master_matrices/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}")        
    threads: 10
    script: 
        "../pyScripts/dates/GO/generate_GO_master_temporal_matrices_exact_degree_match.py"

rule compute_GO_mean_distance_of_future_annotation_and_random_genes_to_already_annotated_genes_using_probabilities_exact_degree_match:
    input:
        matrix_dir = "work_folder/data/dates/GO/ed_master_matrices/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}"
    output:
        mean_matrix = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    threads: 2
    script:
        "../pyScripts/dates/GO/compute_GO_mean_distance_of_future_annotation_and_random_genes_to_already_annotated_genes_using_probabilities.py"

rule compute_GO_true_annotated_genes_quantiles_exact_degree_match:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    output:
        quantile_file = "work_folder/data/dates/GO/ed_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet"
    threads: 2
    script:
        "../pyScripts/dates/GO/compute_GO_true_annotated_genes_quantiles.py"

rule plot_GO_true_vs_permutations_predictive_power_over_time_exact_degree_match:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    output:
        plot_file = "work_folder/data/dates/GO/plots/ed_predictive_power/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_true_vs_perm.png"
    script:
        "../pyScripts/plotting/plot_GO_true_vs_permutations_predictive_power_over_time.py"

rule plot_GO_true_annotated_genes_quantiles_over_time_exact_degree_match:
    input:
        quantile_file = "work_folder/data/dates/GO/ed_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet"
    output:
        plot_file = "work_folder/data/dates/GO/plots/ed_quantiles/quantiles_over_time/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles_over_time.png"
    script:
        "../pyScripts/plotting/plot_GO_true_annotated_genes_quantiles_over_time.py"

rule extract_GO_raw_tta_exact_degree_match:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv"
    output:
        raw_tta_file = "work_folder/data/dates/GO/ed_raw_tta/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_raw_tta.parquet"
    script:
        "../pyScripts/dates/GO/extract_GO_raw_tta.py"

rule plot_GO_raw_tta_distributions_exact_degree_match:
    input:
        raw_tta_file = "work_folder/data/dates/GO/ed_raw_tta/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_raw_tta.parquet"
    output:
        plot_file = "work_folder/data/dates/GO/plots/ed_raw_tta/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_raw_tta_distribution.png"
    script:
        "../pyScripts/plotting/plot_GO_raw_tta_distributions.py"

rule plot_GO_network_growth:
    input:
        network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl"
    output:
        plot = "work_folder/data/dates/GO/plots/{aspect}_network_growth_over_time.png"
    script:
        "../pyScripts/plotting/plot_GO_network_growth.py"

rule plot_GO_mean_adjacency_vs_time_to_annotation_exact_degree_match:
    input:
        mean_adj_dir = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}",
        annot_dates_dir = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}"
    output:
        plot_file = "work_folder/data/dates/GO/plots/ed_tta_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/mean_adjacency_vs_time_to_annotation.png"
    script:
        "../pyScripts/plotting/plot_GO_mean_adjacency_vs_time_to_annotation.py"

rule compute_GO_time_to_annotation_and_mean_adjacency_correlation_exact_degree_match:
    input:
        mean_adj_dir = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}",
        annot_dates_dir = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}"
    output:
        global_stats = "work_folder/data/dates/GO/stats/ed_tta_and_mean_adj_corr/{aspect}_depth_{depth}_cutoff_{cutoff}_spearman_correlation.csv"
    script:
        "../pyScripts/dates/GO/compute_GO_time_to_annotation_and_mean_adjacency_correlation.py"   

rule plot_GO_time_to_annotation_and_mean_adjacency_correlation_exact_degree_match:
    input:
        global_stats = "work_folder/data/dates/GO/stats/ed_tta_and_mean_adj_corr/{aspect}_depth_{depth}_cutoff_{cutoff}_spearman_correlation.csv"
    output:
        waterfall_plot = "work_folder/data/dates/GO/plots/ed_tta_and_mean_adj_corr/{aspect}_depth_{depth}_cutoff_{cutoff}_waterfall_plot.png"
    script:
        "../pyScripts/plotting/plot_GO_time_to_annotation_and_mean_adjacency_correlation_exact_degree_match.py"

rule plot_GO_mean_adjacency_distribution_for_various_thresholds_exact_degree_match:
    input:
        mean_adj_dir = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}",
        annot_dates_dir = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}"
    output:
        plot_file = "work_folder/data/dates/GO/plots/ed_tta_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/mean_adjacency_distribution_multiple_thresholds.pdf"
    script:
        "../pyScripts/plotting/plot_GO_mean_adjacency_distribution_for_various_thresholds.py"

rule plot_GO_mean_adjacency_vs_time_to_annotation_binned:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv"
    output:
        plot_file = "work_folder/data/dates/GO/plots/binned_tta_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacency_vs_time_to_annotation_binned.png"
    script:
        "../pyScripts/plotting/plot_GO_mean_adjacency_vs_time_to_annotation_binned.py"

rule plot_GO_quantiles_vs_time_to_annotation_binned:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv",
        network_file = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl"
    output:
        plot_quantiles = "work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned.png",
        plot_degrees = "work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_degree_vs_time_to_annotation_binned.png"
    script:
        "../pyScripts/plotting/plot_GO_quantiles_vs_time_to_annotation_binned.py"

rule plot_GO_binned_super_predictors_and_density:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv"
    output:
        plot_super_predictors = "work_folder/data/dates/GO/plots/distribution_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_super_predictors_over_time.png",
        plot_density_2d = "work_folder/data/dates/GO/plots/distribution_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_density_2d_quantiles.png"
    script:
        "../pyScripts/plotting/plot_GO_binned_super_predictors_and_density.py"

rule plot_GO_quantiles_vs_time_to_annotation_binned_sliding_window:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv"
    output:
        plot_quantiles = "work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_sliding_window_tta.png"
    script:
        "../pyScripts/plotting/plot_GO_quantiles_vs_time_to_annotation_binned_sliding_window.py"

rule plot_GO_quantiles_vs_time_to_annotation_binned_by_tertiles:
    input:
        mean_adj_file = "work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/GO/first_annotation_dates/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv",
        network_file = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl"
    output:
        plot_quantiles = "work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned_by_tertiles.png",
        plot_degrees = "work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_degree_vs_time_to_annotation_binned_by_tertiles.png"
    script:
        "../pyScripts/plotting/plot_GO_quantiles_vs_time_to_annotation_binned_by_tertiles.py"        

def get_all_GO_master_temporal_matrices_exact_degree_match(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                # Force the checkpoint to finish for this specific combination
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                # Read the newly created pickle file
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # Add the required term folders to our master target list
                target_folders.extend(
                    expand("work_folder/data/dates/GO/ed_master_matrices/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_folders

def get_all_GO_mean_adjacencies_exact_degree_match(wildcards):
    target_files = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                # Force the checkpoint to finish for this specific combination
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                # Read the newly created pickle file
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # Add the required MEAN ADJACENCY parquet files to our master target list
                target_files.extend(
                    expand("work_folder/data/dates/GO/ed_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_true_annotated_genes_quantiles_exact_degree_match(wildcards):
    target_files = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                # Force the checkpoint to finish for this specific combination
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                # Read the newly created pickle file
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # Add the required MEAN ADJACENCY parquet files to our master target list
                target_files.extend(
                    expand("work_folder/data/dates/GO/ed_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_true_vs_permutations_predictive_power_over_time_exact_degree_match(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/ed_predictive_power/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_true_vs_perm.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_true_annotated_genes_quantiles_over_time_exact_degree_match(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/ed_quantiles/quantiles_over_time/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles_over_time.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_raw_tta_exact_degree_match(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/ed_raw_tta/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_raw_tta.parquet",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_raw_tta_distributions_exact_degree_match(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/ed_raw_tta/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_raw_tta_distribution.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_mean_adjacency_vs_time_to_annotation_binned(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/binned_tta_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacency_vs_time_to_annotation_binned.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_quantile_vs_time_to_annotation_binned(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_binned_super_predictors_and_density(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/distribution_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_super_predictors_over_time.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_quantiles_vs_time_to_annotation_binned_sliding_window(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_sliding_window_tta.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_quantile_vs_time_to_annotation_binned_by_tertiles(wildcards):
    target_files = []
    
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                    aspect=a, depth=d, cutoff=c
                ).output.nodes_with_top_5_annotations_pickle
                
                top_annot_df = pd.read_pickle(annot_file)
                my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                
                # ALL terms go straight into the best_predictions folder
                target_files.extend(
                    expand("work_folder/data/dates/GO/plots/binned_tta_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned_by_tertiles.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files