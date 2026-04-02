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
        # All GO Elastic Net necessary files
        lambda wc: [
            # All depths coefficients
            f"work_folder/data/ElasticNet/GO/EN_coefficients/{aspect}/elastic_net_coefficients_depth_{depth}.csv"
            for aspect in ASPECTS
            for depth in get_depths_for_aspect(aspect)
        ] + [
            # All top 20 plots
            f"work_folder/data/ElasticNet/GO/plots/{aspect}/Top/top_coefficients_depth_{depth}.png"
            for aspect in ASPECTS
            for depth in get_depths_for_aspect(aspect)
        ] + [
            # All distribution coefficients
            f"work_folder/data/ElasticNet/GO/plots/{aspect}/Distribution/coefficients_distribution_depth_{depth}.png"
            for aspect in ASPECTS
            for depth in get_depths_for_aspect(aspect)
        ] + [
            # All number of coefficients per depth level
            f"work_folder/data/ElasticNet/GO/plots/{aspect}/annotations_per_depth.png"
            for aspect in ASPECTS 
        ],
        # All HDO cutoff necessary files
        lambda wc: [
            # All complete feature matrices with all cutoffs
            f"work_folder/data/HDO/cutoff/feature_matrices/complete_feature_matrix_cutoff_{cutoff}.csv"
            for cutoff in CUTOFFS
        ] + [
            # All complete elastic net coefficients with all cutoffs
            f"work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete_elastic_net_coefficients_cutoff_{cutoff}.csv"
            for cutoff in CUTOFFS
        ] + [
            # All complete top 20 plots with all cutoffs
            f"work_folder/data/ElasticNet/HDO_cutoff/plots/Top/complete_top_coefficients_cutoff_{cutoff}.png"
            for cutoff in CUTOFFS
        ] + [
            # All complete distribution coefficients with all cutoffs
            f"work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/complete_coefficients_distribution_cutoff_{cutoff}.png"
            for cutoff in CUTOFFS
        ] + [
            # DOIDs lost plot
            "work_folder/data/ElasticNet/HDO_cutoff/plots/complete_HDO_doids_lost.png"
        ]
