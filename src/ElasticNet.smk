import pandas as pd

HDO_DEPTHS = list(range(0, 11))
HDO_DEPTHS_WITH_ANCESTORS = list(range(0, 13))
ASPECTS = ['BP', 'MF', 'CC']

rule fit_elastic_net_HDO:
    input:
        feature_matrix = expand("work_folder/data/HDO/feature_matrices/feature_matrix_with_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS),
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        elastic_net_coefficients = expand("work_folder/data/ElasticNet/HDO/EN_coefficients/elastic_net_coefficients_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS)
    script:
        "compute_HDO_elastic_net_coefficients.py"

rule visualize_HDO_coefficients:
    input:
        elastic_net_coefficients = expand("work_folder/data/ElasticNet/HDO/EN_coefficients/elastic_net_coefficients_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS),
        annotations_per_depth = "work_folder/data/HDO/annotations_per_depth.csv"
    output:
        top_coefficients = expand("work_folder/data/ElasticNet/HDO/plots/Top/top_coefficients_depth_{hdo_depth}.png", hdo_depth  =HDO_DEPTHS),
        coefficients_distribution = expand("work_folder/data/ElasticNet/HDO/plots/Distribution/coefficients_distribution_depth_{hdo_depth}.png", hdo_depth = HDO_DEPTHS),
        annotations_per_depth_plot = "work_folder/data/ElasticNet/HDO/plots/annotations_per_depth.png"
    script:
        "plot_HDO_EN_coefficients.py"

rule fit_elastic_net_HDO_with_ancestors:
    input:
        feature_matrix = expand("work_folder/data/HDO/feature_matrices/full_feature_matrix_with_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS_WITH_ANCESTORS),
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        elastic_net_coefficients = expand("work_folder/data/ElasticNet/HDO_full/EN_coefficients/elastic_net_coefficients_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS_WITH_ANCESTORS)
    script:
        "compute_HDO_elastic_net_coefficients.py"

rule visualize_HDO_coefficients_with_ancestors:
    input:
        elastic_net_coefficients = expand("work_folder/data/ElasticNet/HDO_full/EN_coefficients/elastic_net_coefficients_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS_WITH_ANCESTORS),
        annotations_per_depth = "work_folder/data/HDO/full_annotations_per_depth.csv",
        genes_per_depth = "work_folder/data/HDO/genes_per_depth.csv"
    output:
        top_coefficients = expand("work_folder/data/ElasticNet/HDO_full/plots/Top/top_coefficients_depth_{hdo_depth}.png", hdo_depth = HDO_DEPTHS_WITH_ANCESTORS),
        coefficients_distribution = expand("work_folder/data/ElasticNet/HDO_full/plots/Distribution/coefficients_distribution_depth_{hdo_depth}.png", hdo_depth = HDO_DEPTHS_WITH_ANCESTORS),
        annotations_per_depth_plot = "work_folder/data/ElasticNet/HDO_full/plots/annotations_per_depth.png",
        genes_per_depth_plot = "work_folder/data/ElasticNet/HDO_full/plots/genes_per_depth.png"
    script:
        "plot_HDO_EN_coefficients_with_ancestors.py"

rule fit_complete_elastic_net_with_ancestors:
    input:
        complete_matrix = "work_folder/data/HDO/feature_matrices/complete_matrix_with_ancestors.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/HDO_full/EN_coefficients/complete_elastic_net_coefficients.csv"
    script:
        "compute_complete_HDO_elastic_net_coefficients.py"

rule visualize_complete_HDO_coefficients_with_ancestors:
    input: 
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/HDO_full/EN_coefficients/complete_elastic_net_coefficients.csv"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/HDO_full/plots/Top/complete_top_coefficients.png",
        coefficients_distribution = "work_folder/data/ElasticNet/HDO_full/plots/Distribution/complete_coefficients_distribution.png"
    script:
        "plot_complete_HDO_EN_coefficients_with_ancestors.py"

def get_depths_for_aspect(aspect):
    depth_file = f"work_folder/data/GO/{aspect}_min_aspect_depths.csv"
    with open(depth_file) as f:
        next(f) 
        line = next(f)
        _, min_depth = line.strip().split(":")
        return [str(d) for d in range(int(min_depth) + 1)]

rule fit_elastic_net_GO:
    input:
        feature_matrix = lambda wc: f"work_folder/data/GO/feature_matrices/{wc.aspect}/feature_matrix_with_depth_{wc.depth}.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        outputfile = "work_folder/data/ElasticNet/GO/EN_coefficients/{aspect}/elastic_net_coefficients_depth_{depth}.csv"
    script:
        "compute_GO_elastic_net_coefficients.py"

rule visualize_GO_coefficients:
    input: 
        elastic_net_coefficients = lambda wc: f"work_folder/data/ElasticNet/GO/EN_coefficients/{wc.aspect}/elastic_net_coefficients_depth_{wc.depth}.csv"
    output:
        top_plot = "work_folder/data/ElasticNet/GO/plots/{aspect}/Top/top_coefficients_depth_{depth}.png", 
        distribution_plot = "work_folder/data/ElasticNet/GO/plots/{aspect}/Distribution/coefficients_distribution_depth_{depth}.png"
    script:
        "plot_GO_EN_coefficients.py"


rule visualize_GO_annotations_per_depth:
    input: 
        annotations_per_depth = lambda wc: f"work_folder/data/GO/annotations_per_depth_{wc.aspect}.csv"
    output:
        annotations_per_depth_plot = "work_folder/data/ElasticNet/GO/plots/{aspect}/annotations_per_depth.png" 
    script:
        "plot_GO_annotations_per_depth.py" 