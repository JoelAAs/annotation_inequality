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

# Get all the terms to build the various folders for the master temporal matrices to be called in rule all
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

1014899

# rule compute_GO_distance_of_future_annotation_and_random_genes_to_already_annotated_genes:

# rule compute_GO_true_annotated_genes_quantiles:

# (maybe plot max and minimum quantile for each term but I don't know if it is necessary)

# rule plot_GO_true_annotated_genes_quantiles_over_time: