import pandas as pd

rule compute_GO_full_prediction_matrix_without_permutations:
    input:
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl"
    output:
        full_matrix = "work_folder/data/prediction/GO/{aspect}/matrices/without/full_matrix.parquet"
    threads: 10
    script:
        "../pyScripts/prediction/GO/compute_GO_full_prediction_matrix_without_permutations.py"

rule extract_GO_y_for_each_term:
    input:
        final_network = "work_folder/data/dates/GO/networks_with_dates/{aspect}_final_network.pkl"
    output:
        y_file = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/y_files/{term}_y_file.parquet"
    script:
        "../pyScripts/prediction/GO/extract_GO_y_for_each_term.py"

rule optimize_GO_xgboost_params:
    input:
        full_matrix = "work_folder/data/prediction/GO/{aspect}/matrices/{wo}/full_matrix.parquet",
        y_file = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/y_files/{term}_y_file.parquet"
    output:
        saved_model = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/models/{wo}/{term}_model.joblib",
        best_params = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/parameters/{wo}/{term}_best_params.json",
        tuning_log = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/logs/{wo}/{term}_tuning.log",
        rest_ids = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/splits/{wo}/{term}_rest_ids.csv",
        test_ids = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/splits/{wo}/{term}_test_ids.csv"
    threads: 5
    script:
        "../pyScripts/prediction/GO/optimize_GO_xgboost_params.py"

rule plot_GO_xgboost_tuning:
    input:
        tuning_log = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/logs/{wo}/{term}_tuning.log"
    output:
        plot = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/plots/{wo}/{term}_tuning.png",
        history = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/history/{wo}/{term}_tuning_history.csv"
    params:
        n_initial = 10
    script:
        "../pyScripts/prediction/GO/plot_tuning.py"

rule fit_GO_xgboost_model:
    input:
        full_matrix = "work_folder/data/prediction/GO/{aspect}/matrices/{wo}/full_matrix.parquet",
        y_file = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/y_files/{term}_y_file.parquet",
        best_params = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/parameters/{wo}/{term}_best_params.json",
        rest_ids = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/splits/{wo}/{term}_rest_ids.csv",
        test_ids = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/splits/{wo}/{term}_test_ids.csv"
    output:
        saved_model = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/final_models/{wo}/{term}_model.joblib",
        test_metrics = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/{wo}/{term}_metrics.json",
        test_predictions = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/{wo}/{term}_predictions.parquet",
        shap_values = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/shap_values/{wo}/{term}_shap_values.npy",
        shap_expected_value = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/shap_values/{wo}/{term}_shap_expected_value.npy",
        x_test_features  = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/x_test_features/{wo}/{term}_X_test_features.parquet"
    threads: 5
    script:
        "../pyScripts/prediction/GO/fit_GO_xgboost_model.py"

rule plot_model_results:
    input:
        shap_values = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/shap_values/{wo}/{term}_shap_values.npy",
        x_test_features = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/x_test_features/{wo}/{term}_X_test_features.parquet",
        test_predictions = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/{wo}/{term}_predictions.parquet",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        shap_summary_plot = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/plots/shap_summary/{wo}/{term}_shap_summary.pdf",
        performance_curves = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/plots/performance_curves/{wo}/{term}_curves.pdf",
        shap_top_genes = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/shap_top_genes/{wo}/{term}_shap_top_genes.csv"
    threads: 1
    script:
        "../pyScripts/prediction/GO/plot_model_results.py"

rule extract_GO_novel_candidates:
    input:
        shap_values = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/shap_values/{wo}/{term}_shap_values.npy",
        x_test_features = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/x_test_features/{wo}/{term}_X_test_features.parquet",
        shap_expected_value = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/shap_values/{wo}/{term}_shap_expected_value.npy",
        test_predictions = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/{wo}/{term}_predictions.parquet"
    output:
        novel_candidates_csv = "work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/novel_candidates/{wo}/{term}_novel_candidates.csv",
        waterfall_dir        = directory("work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/plots/waterfalls/{wo}/{term}/")
    threads: 1
    script:
        "../pyScripts/prediction/GO/extract_novel_candidates.py"

def get_all_GO_y_files(wildcards):
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
                    expand("work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/y_files/{term}_y_file.parquet",
                           aspect=a, depth=d, cutoff=c, term=my_terms)
                )
                
    return target_files

def get_all_GO_optimized_params(wildcards):
    target_files = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                for w in WITH_OR_WITHOUT:
                    # Force the checkpoint to finish for this specific combination
                    annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                        aspect=a, depth=d, cutoff=c
                    ).output.nodes_with_top_5_annotations_pickle
                    
                    # Read the newly created pickle file
                    top_annot_df = pd.read_pickle(annot_file)
                    my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                    
                    # Add the required MEAN ADJACENCY parquet files to our master target list
                    target_files.extend(
                        expand("work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/parameters/{wo}/{term}_best_params.json",
                            aspect=a, depth=d, cutoff=c, term=my_terms, wo=w)
                    )
                
    return target_files

def get_all_GO_optimized_params_plots(wildcards):
    target_files = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                for w in WITH_OR_WITHOUT:
                    # Force the checkpoint to finish for this specific combination
                    annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                        aspect=a, depth=d, cutoff=c
                    ).output.nodes_with_top_5_annotations_pickle
                    
                    # Read the newly created pickle file
                    top_annot_df = pd.read_pickle(annot_file)
                    my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                    
                    # Add the required MEAN ADJACENCY parquet files to our master target list
                    target_files.extend(
                        expand("work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/plots/{wo}/{term}_tuning.png",
                            aspect=a, depth=d, cutoff=c, term=my_terms, wo=w)
                    )
                
    return target_files

def get_all_GO_xgboost_models(wildcards):
    target_files = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                for w in WITH_OR_WITHOUT:
                    # Force the checkpoint to finish for this specific combination
                    annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                        aspect=a, depth=d, cutoff=c
                    ).output.nodes_with_top_5_annotations_pickle
                    
                    # Read the newly created pickle file
                    top_annot_df = pd.read_pickle(annot_file)
                    my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                    
                    # Add the required MEAN ADJACENCY parquet files to our master target list
                    target_files.extend(
                        expand("work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/shap_values/{wo}/{term}_shap_values.npy",
                            aspect=a, depth=d, cutoff=c, term=my_terms, wo=w)
                    )
                
    return target_files

def get_all_GO_model_results_plots(wildcards):
    target_files = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                for w in WITH_OR_WITHOUT:
                    # Force the checkpoint to finish for this specific combination
                    annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                        aspect=a, depth=d, cutoff=c
                    ).output.nodes_with_top_5_annotations_pickle
                    
                    # Read the newly created pickle file
                    top_annot_df = pd.read_pickle(annot_file)
                    my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                    
                    # Add the required MEAN ADJACENCY parquet files to our master target list
                    target_files.extend(
                        expand("work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/plots/shap_summary/{wo}/{term}_shap_summary.pdf",
                            aspect=a, depth=d, cutoff=c, term=my_terms, wo=w)
                    )
                
    return target_files

def get_all_GO_novel_candidates(wildcards):
    target_files = []
    
    # Loop through the parameters defined above
    for a in TEMPORAL_MATRICES_ASPECTS:
        for d in TEMPORAL_MATRICES_DEPTHS:
            for c in TEMPORAL_MATRICES_CUTOFFS:
                for w in WITH_OR_WITHOUT:
                    # Force the checkpoint to finish for this specific combination
                    annot_file = checkpoints.find_GO_nodes_with_top_5_annotations.get(
                        aspect=a, depth=d, cutoff=c
                    ).output.nodes_with_top_5_annotations_pickle
                    
                    # Read the newly created pickle file
                    top_annot_df = pd.read_pickle(annot_file)
                    my_terms = [str(term).replace(":", "_") for term in top_annot_df['GO_id'].unique()]
                    
                    # Add the required MEAN ADJACENCY parquet files to our master target list
                    target_files.extend(
                        expand("work_folder/data/prediction/GO/{aspect}/{depth}/{cutoff}/evaluation/novel_candidates/{wo}/{term}_novel_candidates.csv",
                            aspect=a, depth=d, cutoff=c, term=my_terms, wo=w)
                    )
                
    return target_files