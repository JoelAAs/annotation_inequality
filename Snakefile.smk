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
        ]
        
