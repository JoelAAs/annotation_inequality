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
include: "src/CutoffsUpdated.smk"

rule all:
    input:
        # --- GO EN SECTION ---
        [
            f"work_folder/data/ElasticNet/GO/EN_coefficients/{aspect}/elastic_net_coefficients_depth_{depth}.csv"
            for aspect in ASPECTS
            for depth in get_depths_for_aspect(aspect)
        ],
        [
            f"work_folder/data/ElasticNet/GO/plots/{aspect}/Top/top_coefficients_depth_{depth}.png"
            for aspect in ASPECTS
            for depth in get_depths_for_aspect(aspect)
        ],
        [
            f"work_folder/data/ElasticNet/GO/plots/{aspect}/Distribution/coefficients_distribution_depth_{depth}.png"
            for aspect in ASPECTS
            for depth in get_depths_for_aspect(aspect)
        ],
        [
            f"work_folder/data/ElasticNet/GO/plots/{aspect}/annotations_per_depth.png"
            for aspect in ASPECTS 
        ],

        # --- COMPLETE HDO CUTOFF SECTION ---
        expand("work_folder/data/HDO/cutoff/feature_matrices/complete/complete_feature_matrix_cutoff_{cutoff}.csv", cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete/complete_elastic_net_coefficients_cutoff_{cutoff}.csv", cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Top/complete/complete_top_coefficients_cutoff_{cutoff}.png", cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/complete/complete_coefficients_distribution_cutoff_{cutoff}.png", cutoff=CUTOFFS),
        "work_folder/data/ElasticNet/HDO_cutoff/plots/doids_lost/complete_HDO_doids_lost.png",

        # --- SINGLE DEPTH SECTION ---
        expand("work_folder/data/HDO/cutoff/feature_matrices/single_depth/depth_{depth}_feature_matrix_cutoff_{cutoff}.csv", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Top/single_depth/depth_{depth}_top_coefficients_cutoff_{cutoff}.png", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/single_depth/depth_{depth}_coefficients_distribution_cutoff_{cutoff}.png", 
               depth=HDO_DEPTHS_WITH_ANCESTORS, cutoff=CUTOFFS),
        expand("work_folder/data/ElasticNet/HDO_cutoff/plots/doids_lost/depth_{depth}_HDO_doids_lost.png", 
               depth=HDO_DEPTHS_WITH_ANCESTORS)
       