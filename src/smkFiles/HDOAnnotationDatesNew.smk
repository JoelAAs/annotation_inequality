rule download_alliance_of_genome_resources_tsv_data:
    output:
        "work_folder/data/dates/HDO_new/downloads/DISEASE-ALLIANCE_TSV_HUMAN.tsv"
    params:
        url="https://www.alliancegenome.org/download/DISEASE-ALLIANCE_TSV_HUMAN.tsv.gz"
    message:
        "--> [JOB] Downloading and extracting Alliance of Genome Resources data: DISEASE-ALLIANCE_TSV_HUMAN"
    shell:
        "wget -qO- {params.url} | gunzip -c > {output}"

rule download_hgnc_mapping_new:
    output:
        hgnc_file = "work_folder/data/dates/HDO_new/downloads/hgnc_complete_set.txt"
    shell:
        """
        curl -fL --retry 3 -o {output.hgnc_file} \
        https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
        """

rule find_gene_annotation_associations_in_common_with_alliance_of_genome_resources:
    input:
        agr_file = "work_folder/data/dates/HDO_new/downloads/DISEASE-ALLIANCE_TSV_HUMAN.tsv",
        annot_file = "work_folder/data/HDO/new_annotations_per_gene_with_ancestors.csv",
        hgnc_file = "work_folder/data/dates/HDO_new/downloads/hgnc_complete_set.txt"
    output:
        gene_annot_associations = "work_folder/data/dates/HDO_new/references/gene_annot_associations.parquet"
    script:
        "../pyScripts/dates/HDO_new/find_gene_annotation_associations_in_common_with_alliance_of_genome_resources.py"

rule translate_agr_and_rgd_to_pmid:
    input:
        associations = "work_folder/data/dates/HDO_new/references/gene_annot_associations.parquet"
    output:
        translated_refs = "work_folder/data/dates/HDO_new/references/translated_references.parquet"
    script:
        "../pyScripts/dates/HDO_new/translate_agr_and_rgd_to_pmid.py"

rule fetch_pubmed_dates:
    input:
        associations = "work_folder/data/dates/HDO_new/references/gene_annot_associations.parquet"
    output:
        pmid_dates = "work_folder/data/dates/HDO_new/references/pmid_publication_dates.parquet"
    script:
        "../pyScripts/dates/HDO_new/fetch_pubmed_dates.py"

rule merge_pubmed_ids_to_associations:
    input:
        associations = "work_folder/data/dates/HDO_new/references/gene_annot_associations.parquet",
        pmid_dates = "work_folder/data/dates/HDO_new/references/pmid_publication_dates.parquet"
    output:
        full_dataset = "work_folder/data/dates/HDO_new/associated/full_annotated_network.parquet",
        only_associated = "work_folder/data/dates/HDO_new/associated/time_resolved_network.parquet"
    script:
        "../pyScripts/dates/HDO_new/merge_pubmed_ids_to_associations.py"

rule inject_dates_to_HDO_new_network:
    input:
        network = "work_folder/data/network/HDO/HDO_bait_prey_publications_network.pkl",
        only_associated = "work_folder/data/dates/HDO_new/associated/time_resolved_network.parquet"
    output:
        network_with_dates = "work_folder/data/dates/HDO_new/networks/network_with_found_dates.pkl"
    script:
        "../pyScripts/dates/HDO_new/inject_dates_to_HDO_new_network.py"

# Add both missing dates and interaction dates
rule assign_dates_to_HDO_new_missing_ones:
    input:
        network_with_dates = "work_folder/data/dates/HDO_new/networks/network_with_found_dates.pkl",
        ontology = "work_folder/data/HDO/doid.obo",
        bp_publications = "work_folder/data/intact/bait_prey_publications_complete_dates.pq"
    output:
        final_network = "work_folder/data/dates/HDO_new/networks/final_network.pkl"
    script:
        "../pyScripts/dates/HDO_new/assign_dates_to_HDO_new_missing_ones.py"

checkpoint find_HDO_new_nodes_with_top_5_annotations:
    input:
        final_network = "work_folder/data/dates/HDO_new/networks/final_network.pkl",
        coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        nodes_with_top_5_annotations = "work_folder/data/dates/HDO_new/top_5/nodes_with_top_5_annot_depth_{depth}_cutoff_{cutoff}.pkl"
    script: 
        "../pyScripts/dates/HDO_new/find_HDO_new_nodes_with_top_5_annotations.py"

rule generate_HDO_master_temporal_matrices_exact_degree_match:
    input: 
        final_network = "work_folder/data/dates/HDO_new/networks/final_network.pkl",
        top_annot_df = "work_folder/data/dates/HDO_new/top_5/nodes_with_top_5_annot_depth_{depth}_cutoff_{cutoff}.pkl"
    output: 
        matrix_dir = directory("work_folder/data/dates/HDO_new/ed_master_matrices/depth_{depth}_cutoff_{cutoff}/{term}")        
    threads: 30
    script: 
        "../pyScripts/dates/HDO_new/generate_HDO_master_temporal_matrices_exact_degree_match.py"

rule compute_HDO_mean_distance_of_future_annotation_and_random_genes_to_already_annotated_genes_using_probabilities_exact_degree_match:
    input:
        matrix_dir = "work_folder/data/dates/HDO_new/ed_master_matrices/depth_{depth}_cutoff_{cutoff}/{term}"
    output:
        mean_matrix = "work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    threads: 1
    script:
        "../pyScripts/dates/HDO_new/compute_HDO_mean_distance_of_future_annotation_and_random_genes_to_already_annotated_genes_using_probabilities.py"

rule compute_HDO_true_annotated_genes_quantiles_exact_degree_match:
    input:
        mean_adj_file = "work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    output:
        quantile_file = "work_folder/data/dates/HDO_new/ed_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet"
    threads: 1
    script:
        "../pyScripts/dates/HDO_new/compute_HDO_true_annotated_genes_quantiles.py"

rule plot_HDO_true_vs_permutations_predictive_power_over_time_exact_degree_match:
    input:
        mean_adj_file = "work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet"
    output:
        plot_file = "work_folder/data/dates/HDO_new/plots/ed_predictive_power/depth_{depth}_cutoff_{cutoff}/{term}_true_vs_perm.png"
    script:
        "../pyScripts/dates/HDO_new/plots/plot_HDO_true_vs_permutations_predictive_power_over_time_exact_degree_match.py"

rule plot_HDO_true_annotated_genes_quantiles_over_time_exact_degree_match:
    input:
        quantile_file = "work_folder/data/dates/HDO_new/ed_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet"
    output:
        plot_file = "work_folder/data/dates/HDO_new/plots/ed_quantiles/quantiles_over_time/depth_{depth}_cutoff_{cutoff}/{term}_quantiles_over_time.png"
    script:
        "../pyScripts/dates/HDO_new/plots/plot_HDO_true_annotated_genes_quantiles_over_time_exact_degree_match.py"

rule retrieve_first_annotation_dates_for_top_5_HDO_new_annotations:
    input:
        final_network = "work_folder/data/dates/HDO_new/networks/final_network.pkl",
        nodes_with_top_5_annotations = "work_folder/data/dates/HDO_new/top_5/nodes_with_top_5_annot_depth_{depth}_cutoff_{cutoff}.pkl"
    output:
        first_annotation_dates_dir = directory("work_folder/data/dates/HDO_new/first_annotation_dates/depth_{depth}_cutoff_{cutoff}")
    threads: 5
    script:
        "../pyScripts/dates/HDO_new/retrieve_first_annotation_dates_for_top_5_HDO_new_annotations.py"

rule plot_HDO_mean_adjacency_vs_time_to_annotation_exact_degree_match:
    input:
        mean_adj_dir = "work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}",
        annot_dates_dir = "work_folder/data/dates/HDO_new/first_annotation_dates/depth_{depth}_cutoff_{cutoff}"
    output:
        plot_file = "work_folder/data/dates/HDO_new/plots/ed_tta_mean_adjacencies/depth_{depth}_cutoff_{cutoff}/mean_adjacency_vs_time_to_annotation.png"
    script:
        "../pyScripts/dates/HDO_new/plots/plot_HDO_mean_adjacency_vs_time_to_annotation_exact_degree_match.py"

rule compute_HDO_time_to_annotation_and_mean_adjacency_correlation_exact_degree_match:
    input:
        mean_adj_dir = "work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}",
        annot_dates_dir = "work_folder/data/dates/HDO_new/first_annotation_dates/depth_{depth}_cutoff_{cutoff}"
    output:
        global_stats = "work_folder/data/dates/HDO_new/stats/ed_tta_and_mean_adj_corr/depth_{depth}_cutoff_{cutoff}_spearman_correlation.csv"
    script:
        "../pyScripts/dates/HDO_new/compute_HDO_time_to_annotation_and_mean_adjacency_correlation_exact_degree_match.py"

rule plot_HDO_time_to_annotation_and_mean_adjacency_correlation_exact_degree_match:
    input:
        global_stats = "work_folder/data/dates/HDO_new/stats/ed_tta_and_mean_adj_corr/depth_{depth}_cutoff_{cutoff}_spearman_correlation.csv"
    output:
        waterfall_plot = "work_folder/data/dates/HDO_new/plots/ed_tta_and_mean_adj_corr/depth_{depth}_cutoff_{cutoff}_waterfall_plot.png"
    script:
        "../pyScripts/dates/HDO_new/plots/plot_HDO_time_to_annotation_and_mean_adjacency_correlation_exact_degree_match.py"

rule plot_HDO_quantiles_vs_time_to_annotation_binned:
    input:
        mean_adj_file = "work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/HDO_new/first_annotation_dates/depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv",
        network_file = "work_folder/data/dates/HDO_new/networks/final_network.pkl"
    output:
        plot_quantiles = "work_folder/data/dates/HDO_new/plots/binned_tta_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned.png",
        plot_degrees = "work_folder/data/dates/HDO_new/plots/binned_tta_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_degree_vs_time_to_annotation_binned.png"
    script:
        "../pyScripts/dates/HDO_new/plots/plot_HDO_quantiles_vs_time_to_annotation_binned.py"

rule plot_HDO_quantiles_vs_time_to_annotation_binned_by_tertiles:
    input:
        mean_adj_file = "work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
        annot_dates_file = "work_folder/data/dates/HDO_new/first_annotation_dates/depth_{depth}_cutoff_{cutoff}/{term}_first_annotation_dates.csv",
        network_file = "work_folder/data/dates/HDO_new/networks/final_network.pkl"
    output:
        plot_quantiles = "work_folder/data/dates/HDO_new/plots/binned_tta_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned_by_tertiles.png",
        plot_degrees = "work_folder/data/dates/HDO_new/plots/binned_tta_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_degree_vs_time_to_annotation_binned_by_tertiles.png"
    script:
        "../pyScripts/dates/HDO_new/plots/plot_HDO_quantiles_vs_time_to_annotation_binned_by_tertiles.py"

def get_all_HDO_master_temporal_matrices_exact_degree_match(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above (assuming these exist in your Snakefile)
    for d in TEMPORAL_MATRICES_DEPTHS:
        for c in TEMPORAL_MATRICES_CUTOFFS:
            
            # Force the checkpoint to finish for this specific combination.
            annot_file = checkpoints.find_HDO_new_nodes_with_top_5_annotations.get(
                depth=d, cutoff=c
            ).output.nodes_with_top_5_annotations
            
            # Read the newly created pickle file
            top_annot_df = pd.read_pickle(annot_file)
            
            # Extract the HDO terms and replace colons with underscores for safe directory names
            my_terms = [str(term).replace(":", "_") for term in top_annot_df['HDO_doid'].unique()]
            
            # Add the required term folders to our master target list (Removed 'aspect')
            target_folders.extend(
                expand("work_folder/data/dates/HDO_new/ed_master_matrices/depth_{depth}_cutoff_{cutoff}/{term}",
                       depth=d, cutoff=c, term=my_terms)
            )
            
    return target_folders

def get_all_HDO_mean_adjacencies_exact_degree_match(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above (assuming these exist in your Snakefile)
    for d in TEMPORAL_MATRICES_DEPTHS:
        for c in TEMPORAL_MATRICES_CUTOFFS:
            
            # Force the checkpoint to finish for this specific combination.
            annot_file = checkpoints.find_HDO_new_nodes_with_top_5_annotations.get(
                depth=d, cutoff=c
            ).output.nodes_with_top_5_annotations
            
            # Read the newly created pickle file
            top_annot_df = pd.read_pickle(annot_file)
            
            # Extract the HDO terms and replace colons with underscores for safe directory names
            my_terms = [str(term).replace(":", "_") for term in top_annot_df['HDO_doid'].unique()]
            
            # Add the required term folders to our master target list (Removed 'aspect')
            target_folders.extend(
                expand("work_folder/data/dates/HDO_new/ed_mean_adjacencies/depth_{depth}_cutoff_{cutoff}/{term}_mean_adjacencies.parquet",
                       depth=d, cutoff=c, term=my_terms)
            )
            
    return target_folders

def get_all_HDO_true_annotated_genes_quantiles_exact_degree_match(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above (assuming these exist in your Snakefile)
    for d in TEMPORAL_MATRICES_DEPTHS:
        for c in TEMPORAL_MATRICES_CUTOFFS:
            
            # Force the checkpoint to finish for this specific combination.
            annot_file = checkpoints.find_HDO_new_nodes_with_top_5_annotations.get(
                depth=d, cutoff=c
            ).output.nodes_with_top_5_annotations
            
            # Read the newly created pickle file
            top_annot_df = pd.read_pickle(annot_file)
            
            # Extract the HDO terms and replace colons with underscores for safe directory names
            my_terms = [str(term).replace(":", "_") for term in top_annot_df['HDO_doid'].unique()]
            
            # Add the required term folders to our master target list (Removed 'aspect')
            target_folders.extend(
                expand("work_folder/data/dates/HDO_new/ed_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_quantiles.parquet",
                       depth=d, cutoff=c, term=my_terms)
            )
            
    return target_folders

def get_all_HDO_true_vs_permutations_predictive_power_over_time_exact_degree_match(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above (assuming these exist in your Snakefile)
    for d in TEMPORAL_MATRICES_DEPTHS:
        for c in TEMPORAL_MATRICES_CUTOFFS:
            
            # Force the checkpoint to finish for this specific combination.
            annot_file = checkpoints.find_HDO_new_nodes_with_top_5_annotations.get(
                depth=d, cutoff=c
            ).output.nodes_with_top_5_annotations
            
            # Read the newly created pickle file
            top_annot_df = pd.read_pickle(annot_file)
            
            # Extract the HDO terms and replace colons with underscores for safe directory names
            my_terms = [str(term).replace(":", "_") for term in top_annot_df['HDO_doid'].unique()]
            
            # Add the required term folders to our master target list (Removed 'aspect')
            target_folders.extend(
                expand("work_folder/data/dates/HDO_new/plots/ed_predictive_power/depth_{depth}_cutoff_{cutoff}/{term}_true_vs_perm.png",
                       depth=d, cutoff=c, term=my_terms)
            )
            
    return target_folders

def get_all_HDO_true_annotated_genes_quantiles_over_time_exact_degree_match(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above (assuming these exist in your Snakefile)
    for d in TEMPORAL_MATRICES_DEPTHS:
        for c in TEMPORAL_MATRICES_CUTOFFS:
            
            # Force the checkpoint to finish for this specific combination.
            annot_file = checkpoints.find_HDO_new_nodes_with_top_5_annotations.get(
                depth=d, cutoff=c
            ).output.nodes_with_top_5_annotations
            
            # Read the newly created pickle file
            top_annot_df = pd.read_pickle(annot_file)
            
            # Extract the HDO terms and replace colons with underscores for safe directory names
            my_terms = [str(term).replace(":", "_") for term in top_annot_df['HDO_doid'].unique()]
            
            # Add the required term folders to our master target list (Removed 'aspect')
            target_folders.extend(
                expand("work_folder/data/dates/HDO_new/plots/ed_quantiles/quantiles_over_time/depth_{depth}_cutoff_{cutoff}/{term}_quantiles_over_time.png",
                       depth=d, cutoff=c, term=my_terms)
            )
            
    return target_folders

def get_all_HDO_quantiles_vs_time_to_annotation_binned(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above (assuming these exist in your Snakefile)
    for d in TEMPORAL_MATRICES_DEPTHS:
        for c in TEMPORAL_MATRICES_CUTOFFS:
            
            # Force the checkpoint to finish for this specific combination.
            annot_file = checkpoints.find_HDO_new_nodes_with_top_5_annotations.get(
                depth=d, cutoff=c
            ).output.nodes_with_top_5_annotations
            
            # Read the newly created pickle file
            top_annot_df = pd.read_pickle(annot_file)
            
            # Extract the HDO terms and replace colons with underscores for safe directory names
            my_terms = [str(term).replace(":", "_") for term in top_annot_df['HDO_doid'].unique()]
            
            # Add the required term folders to our master target list (Removed 'aspect')
            target_folders.extend(
                expand("work_folder/data/dates/HDO_new/plots/binned_tta_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned.png",
                       depth=d, cutoff=c, term=my_terms)
            )
            
    return target_folders

def get_all_HDO_quantiles_vs_time_to_annotation_binned_by_tertiles(wildcards):
    target_folders = []
    
    # Loop through the parameters defined above (assuming these exist in your Snakefile)
    for d in TEMPORAL_MATRICES_DEPTHS:
        for c in TEMPORAL_MATRICES_CUTOFFS:
            
            # Force the checkpoint to finish for this specific combination.
            annot_file = checkpoints.find_HDO_new_nodes_with_top_5_annotations.get(
                depth=d, cutoff=c
            ).output.nodes_with_top_5_annotations
            
            # Read the newly created pickle file
            top_annot_df = pd.read_pickle(annot_file)
            
            # Extract the HDO terms and replace colons with underscores for safe directory names
            my_terms = [str(term).replace(":", "_") for term in top_annot_df['HDO_doid'].unique()]
            
            # Add the required term folders to our master target list (Removed 'aspect')
            target_folders.extend(
                expand("work_folder/data/dates/HDO_new/plots/binned_tta_quantiles/depth_{depth}_cutoff_{cutoff}/{term}_quantile_vs_time_to_annotation_binned_by_tertiles.png",
                       depth=d, cutoff=c, term=my_terms)
            )
            
    return target_folders