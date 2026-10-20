rule calculate_BP_HDO_jaccard:
    input:
        go_file = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_BP_annotations_depth_5_cutoff_20.pkl",
        hdo_file = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_5_cutoff_20.pkl",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        jaccard_tsv = "work_folder/data/dates/Covariates/BP_HDO_jaccard_overlap.tsv",
        jaccard_plot = "work_folder/data/dates/Covariates/plots/check_en_model/BP_HDO_jaccard_heatmap.png"
    script:
        "../pyScripts/dates/Covariates/calculate_BP_HDO_jaccard.py"

rule create_check_elastic_net_model:
    input:
        go_file = "work_folder/data/dates/GO/top_5_annotations/nodes_with_top_5_BP_annotations_depth_5_cutoff_20.pkl",
        hdo_file = "work_folder/data/dates/HDO/top_5_annotations/nodes_with_top_5_HDO_annotations_depth_5_cutoff_20.pkl",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        confusion_matrix = "work_folder/data/dates/Covariates/check_en_model/confusion_matrix.parquet",
        en_coefficients = "work_folder/data/dates/Covariates/check_en_model/en_coefficients.csv"
    script:
        "../pyScripts/dates/Covariates/create_check_elastic_net_model.py"

rule visualize_check_elastic_net_results:
    input:
        confusion_matrix = "work_folder/data/dates/Covariates/check_en_model/confusion_matrix.parquet",
        en_coefficients = "work_folder/data/dates/Covariates/check_en_model/en_coefficients.csv"
    output:
        matrix_plot = "work_folder/data/dates/Covariates/plots/check_en_model/confusion_matrix_plot.png",
        coef_plot = "work_folder/data/dates/Covariates/plots/check_en_model/en_coefficients_plot.png"
    script:
        "../pyScripts/plotting/visualize_check_elastic_net_results.py"

rule check_elastic_net_covariation_heatmap:
    input:
        confusion_matrix = "work_folder/data/dates/Covariates/check_en_model/confusion_matrix.parquet"
    output:
        heatmap = "work_folder/data/dates/Covariates/plots/check_en_model/heatmap.png"
    script:
        "../pyScripts/plotting/check_elastic_net_covariation_heatmap.py"

rule compute_mixed_feature_matrix:
    input:
        hdo_annotations = "work_folder/data/HDO/new_annotations_per_gene_with_ancestors.csv",
        go_annotations = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        bp_df = "work_folder/data/intact/bait_prey_publications_complete_dates.pq"
    output:
        feature_matrix = "work_folder/data/dates/Covariates/feature_matrix/{aspect}xHDO_feature_matrix.parquet"
    script:
        "../pyScripts/dates/Covariates/compute_mixed_feature_matrix.py"

rule obtain_mixed_feature_matrix_parameters:
    input:
        feature_matrix = "work_folder/data/dates/Covariates/feature_matrix/{aspect}xHDO_feature_matrix.parquet"
    output:
        parameters = "work_folder/data/dates/Covariates/feature_matrix/{aspect}xHDO_parameters.parquet",
        counts_plot = "work_folder/data/dates/Covariates/plots/feature_matrix/{aspect}xHDO_counts_plot.png",
        terms_drop_plot = "work_folder/data/dates/Covariates/plots/feature_matrix/{aspect}xHDO_genes_drop_plot.png",
        dropped_genes_df = "work_folder/data/dates/Covariates/feature_matrix/{aspect}xHDO_dropped_genes_percentages.parquet"
    script:
        "../pyScripts/dates/Covariates/obtain_mixed_feature_matrix_parameters.py"

rule fit_elastic_net_model_to_mixed_feature_matrix_old:
    input:
        feature_matrix = "work_folder/data/dates/Covariates/feature_matrix/{aspect}xHDO_feature_matrix.parquet",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net_old/coefficients/{aspect}xHDO_coefficients.parquet",
        en_metrics = "work_folder/data/dates/Covariates/elastic_net_old/metrics/{aspect}xHDO_metrics.tsv",
        en_model = "work_folder/data/dates/Covariates/elastic_net_old/models/{aspect}xHDO_model.rds"
    params:
        min_genes = 20
    script:
        "../pyScripts/dates/Covariates/fit_elastic_net_model_to_mixed_feature_matrix_old.R"

rule fit_elastic_net_model_to_mixed_feature_matrix_with_n_annot:
    input:
        feature_matrix = "work_folder/data/dates/Covariates/feature_matrix/{aspect}xHDO_feature_matrix.parquet",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_w_n_annot/{aspect}xHDO_coefficients.parquet",
        en_metrics = "work_folder/data/dates/Covariates/elastic_net/metrics_w_n_annot/{aspect}xHDO_metrics.tsv",
        en_model = "work_folder/data/dates/Covariates/elastic_net/models_w_n_annot/{aspect}xHDO_model.rds",
        cv_curves = "work_folder/data/dates/Covariates/plots/elastic_net/{aspect}xHDO_cv_curves_w_n_annot.pdf"
    params:
        min_genes = 20,
        use_annot_covariate = True,
        standardize = True,
        seed = 42
    script:
        "../pyScripts/dates/Covariates/fit_elastic_net_model_to_mixed_feature_matrix_with_n_annot.R"

rule visualize_mixed_elastic_net_model_old:
    input:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net_old/coefficients/{aspect}xHDO_coefficients.parquet",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        top20_plot = "work_folder/data/dates/Covariates/plots/elastic_net_old/{aspect}xHDO_top20_coefficients.png"
    script:
        "../pyScripts/plotting/visualize_mixed_elastic_net_model.py"

rule visualize_mixed_elastic_net_model:
    input:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_w_n_annot/{aspect}xHDO_coefficients.parquet",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        top20_plot = "work_folder/data/dates/Covariates/plots/elastic_net/{aspect}xHDO_top20_coefficients.png"
    script:
        "../pyScripts/plotting/visualize_mixed_elastic_net_model.py"

rule fit_elastic_net_model_to_mixed_feature_matrix_depth_specific:
    input:
        feature_matrix_hdo = "work_folder/data/HDO/feature_matrices/full_feature_matrix_with_depth_5.csv",
        feature_matrix_go = "work_folder/data/GO/feature_matrices/single_depth/{aspect}_depth_5_feature_matrix.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_single_depth/{aspect}xHDO_coefficients.parquet",
        en_metrics = "work_folder/data/dates/Covariates/elastic_net/metrics_single_depth/{aspect}xHDO_metrics.tsv",
        en_model = "work_folder/data/dates/Covariates/elastic_net/models_single_depth/{aspect}xHDO_model.rds",
        cv_curves = "work_folder/data/dates/Covariates/plots/elastic_net/{aspect}xHDO_cv_curves_single_depth.pdf"
    params:
        min_genes = 20,
        use_annot_covariate = True,
        standardize = True,
        seed = 42
    script:
        "../pyScripts/dates/Covariates/fit_elastic_net_model_to_mixed_feature_matrix_single_depth.R"

rule visualize_mixed_elastic_net_model_depth_specific:
    input:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_single_depth/{aspect}xHDO_coefficients.parquet",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        top20_plot = "work_folder/data/dates/Covariates/plots/elastic_net_single_depth/{aspect}xHDO_top20_coefficients.png"
    script:
        "../pyScripts/plotting/visualize_mixed_elastic_net_model.py"

rule fit_elastic_net_model_to_mixed_feature_matrix_with_n_annot_hdo_only:
    input:
        feature_matrix = "work_folder/data/HDO/feature_matrices/full_feature_matrix_with_depth_5.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_hdo_only/HDO_coefficients.parquet",
        en_metrics = "work_folder/data/dates/Covariates/elastic_net/metrics_hdo_only/HDO_metrics.tsv",
        en_model = "work_folder/data/dates/Covariates/elastic_net/models_hdo_only/HDO_model.rds",
        cv_curves = "work_folder/data/dates/Covariates/plots/elastic_net_hdo_only/HDO_cv_curves_w_n_annot.pdf"
    params:
        min_genes = 20,
        use_annot_covariate = True,
        standardize = True,
        seed = 42
    script:
        "../pyScripts/dates/Covariates/fit_elastic_net_model_to_mixed_feature_matrix_with_n_annot_single_ontology.R"

rule fit_elastic_net_model_to_mixed_feature_matrix_with_n_annot_go_only:
    input:
        feature_matrix = "work_folder/data/GO/feature_matrices/single_depth/{aspect}_depth_5_feature_matrix.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_go_only/{aspect}_coefficients.parquet",
        en_metrics = "work_folder/data/dates/Covariates/elastic_net/metrics_go_only/{aspect}_metrics.tsv",
        en_model = "work_folder/data/dates/Covariates/elastic_net/models_go_only/{aspect}_model.rds",
        cv_curves = "work_folder/data/dates/Covariates/plots/elastic_net_go_only/{aspect}_cv_curves_w_n_annot.pdf"
    params:
        min_genes = 20,
        use_annot_covariate = True,
        standardize = True,
        seed = 42
    script:
        "../pyScripts/dates/Covariates/fit_elastic_net_model_to_mixed_feature_matrix_with_n_annot_single_ontology.R"

rule visualize_mixed_elastic_net_model_hdo_only:
    input:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_hdo_only/HDO_coefficients.parquet",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        top20_plot = "work_folder/data/dates/Covariates/plots/elastic_net_hdo_only/HDO_top20_coefficients.png"
    script:
        "../pyScripts/plotting/visualize_mixed_elastic_net_model.py"

rule visualize_mixed_elastic_net_model_go_only:
    input:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_go_only/{aspect}_coefficients.parquet",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        top20_plot = "work_folder/data/dates/Covariates/plots/elastic_net_go_only/{aspect}_top20_coefficients.png"
    script:
        "../pyScripts/plotting/visualize_mixed_elastic_net_model.py"

rule test_en_model:
    input:
        feature_matrix = "work_folder/data/GO/feature_matrices/single_depth/{aspect}_depth_{depth}_feature_matrix.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        single_depth_elastic_net_coefficients = "work_folder/data/dates/Covariates/elastic_net_test/{aspect}_depth_{depth}coefficients.parquet"
    script:
        "../pyScripts/dates/Covariates/test_en_model.py"

rule visualize_test_en:
    input:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net_test/{aspect}_depth_{depth}coefficients.parquet",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        top20_plot = "work_folder/data/dates/Covariates/plots/elastic_net_test/{aspect}_depth_{depth}_top20_coefficients.png"
    script:
        "../pyScripts/plotting/visualize_mixed_elastic_net_model.py"

rule adapt_new_mixed_model_to_old_one:
    input:
        feature_matrix_hdo = "work_folder/data/HDO/feature_matrices/full_feature_matrix_with_depth_5.csv",
        feature_matrix_go = "work_folder/data/GO/feature_matrices/single_depth/{aspect}_depth_5_feature_matrix.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_adapted/{aspect}xHDO_coefficients.parquet",
        en_metrics = "work_folder/data/dates/Covariates/elastic_net/metrics_adapted/{aspect}xHDO_metrics.tsv",
        en_model = "work_folder/data/dates/Covariates/elastic_net/models_adapted/{aspect}xHDO_model.rds",
        cv_curves = "work_folder/data/dates/Covariates/plots/elastic_net_adapted/{aspect}xHDO_cv_curves.pdf"
    params:
        min_genes = 20,
        standardize = True,
        seed = 42
    script:
        "../pyScripts/dates/Covariates/adapt_new_mixed_model_to_old_one.R"

rule visualize_adapted_mixed_model:
    input:
        en_coefficients = "work_folder/data/dates/Covariates/elastic_net/coefficients_adapted/{aspect}xHDO_coefficients.parquet",
        go_obo = "work_folder/data/GO/go-basic.obo",
        hdo_obo = "work_folder/data/HDO/doid.obo"
    output:
        top20_plot = "work_folder/data/dates/Covariates/plots/elastic_net_adapted/{aspect}xHDO_top20_coefficients.png"
    script:
        "../pyScripts/plotting/visualize_mixed_elastic_net_model.py"