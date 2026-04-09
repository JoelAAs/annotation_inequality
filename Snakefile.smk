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

rule all:
    input:
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
        get_all_GO_single_depth_en_plots,
        get_all_GO_lost_ids_plots,

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

        # --- HDO DENDROGRAM OF COEFFICIENTS SECTION ---
         expand("work_folder/data/dendrograms/HDO/all_coefficients/all_coefficients_cutoff_{cutoff}.csv",
               cutoff = CUTOFFS),
       