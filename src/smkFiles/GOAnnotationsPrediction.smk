import pandas as pd

TEMPORAL_MATRICES_ASPECTS = ["BP"]
TEMPORAL_MATRICES_DEPTHS = [5]
TEMPORAL_MATRICES_CUTOFFS = [20]

rule generate_GO_network_degree_bins:
    input: 
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl"
    output: 
        bins_json = "work_folder/data/dates/GO/bins/{aspect}_degree_bins.json"
    script:
        "../pyScripts/dates/GO/generate_GO_network_degree_bins.py"

# Auxiliary function to get all the terms to build the various folders for the master temporal matrices to be called in rule all
def get_all_GO_master_matrices(wildcards):
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
                    expand("work_folder/data/dates/GO/probability_master_matrices/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_folders

rule generate_GO_master_temporal_matrices:
    input: 
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl",
        top_annot_df = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_{depth}_cutoff_{cutoff}.pkl",
        bins_json = "work_folder/data/dates/GO/bins/{aspect}_degree_bins.json"
    output: 
        matrix_dir = directory("work_folder/data/dates/GO/master_matrices/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}")        
    threads: 10
    script: 
        "../pyScripts/dates/GO/generate_GO_master_temporal_matrices.py"

rule generate_GO_master_temporal_matrices_using_probabilities:
    input:
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl",
        top_annot_df = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_{aspect}_annotations_depth_{depth}_cutoff_{cutoff}.pkl",
        bins_json = "work_folder/data/dates/GO/bins/{aspect}_degree_bins.json"
    output:
        matrix_dir = directory("work_folder/data/dates/GO/probability_master_matrices/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}")
    threads: 10
    script: 
        "../pyScripts/dates/GO/generate_GO_master_temporal_matrices_using_probabilities.py"

# Auxiliary function for rule all to get all the temporal matrices for one term and compute the mean adjacency values for each date and for each permutation
def get_all_GO_mean_adjacencies(wildcards):
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
                    expand("work_folder/data/dates/GO/probability_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

rule compute_GO_mean_distance_of_future_annotation_and_random_genes_to_already_annotated_genes_using_probabilities:
    input:
        matrix_dir = "work_folder/data/dates/GO/probability_master_matrices/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}"
    output:
        mean_matrix = "work_folder/data/dates/GO/probability_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    threads: 2
    script:
        "../pyScripts/dates/GO/compute_GO_mean_distance_of_future_annotation_and_random_genes_to_already_annotated_genes_using_probabilities.py"

def get_all_GO_true_annotated_genes_quantiles(wildcards):
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
                    expand("work_folder/data/dates/GO/quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

rule compute_GO_true_annotated_genes_quantiles:
    input:
        mean_adj_file = "work_folder/data/dates/GO/probability_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    output:
        quantile_file = "work_folder/data/dates/GO/quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet"
    threads: 2
    script:
        "../pyScripts/dates/GO/compute_GO_true_annotated_genes_quantiles.py"

def get_all_GO_best_quantile_plots(wildcards):
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
                    expand("work_folder/data/dates/GO/plots/quantiles/best_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_best_quantile.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

rule plot_GO_best_quantiles:
    input:
        mean_adj_file = "work_folder/data/dates/GO/probability_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    output:
        plot_file = "work_folder/data/dates/GO/plots/quantiles/best_quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_best_quantile.png"
    threads: 1
    script:
        "../pyScripts/plotting/plot_GO_best_quantiles.py"

def get_all_GO_true_annotated_genes_quantiles_over_time(wildcards):
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
                    expand("work_folder/data/dates/GO/plots/quantiles_over_time/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles_over_time.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

rule plot_GO_true_annotated_genes_quantiles_over_time:
    input:
        quantile_file = "work_folder/data/dates/GO/quantiles/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet"
    output:
        plot_file = "work_folder/data/dates/GO/plots/quantiles_over_time/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_quantiles_over_time.png"
    script:
        "../pyScripts/plotting/plot_GO_true_annotated_genes_quantiles_over_time.py"

def get_all_GO_true_vs_permutations_predictive_power_over_time(wildcards):
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
                    expand("work_folder/data/dates/GO/plots/predictive_power/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_true_vs_perm.png",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

rule plot_GO_true_vs_permutations_predictive_power_over_time:
    input:
        mean_adj_file = "work_folder/data/dates/GO/probability_mean_adjacencies/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    output:
        plot_file = "work_folder/data/dates/GO/plots/predictive_power/{aspect}_depth_{depth}_cutoff_{cutoff}/{term}_true_vs_perm.png"
    script:
        "../pyScripts/plotting/plot_GO_true_vs_permutations_predictive_power_over_time.py"