HDO_DEPTHS = list(range(0, 11))
ASPECTS = ['BP', 'MF', 'CC']
INPUT_FOLDERS_GO = [directory(f"work_folder/data/GO/feature_matrices/{aspect}") for aspect in ASPECTS]
OUTPUT_FOLDERS_GO = [directory(f"work_folder/data/ElasticNet/GO/EN_coefficients/{aspect}") for aspect in ASPECTS]

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
        top_coefficients = expand("work_folder/data/ElasticNet/HDO/plots/top_coefficients_depth_{hdo_depth}.png", hdo_depth  =HDO_DEPTHS),
        coefficients_distribution = expand("work_folder/data/ElasticNet/HDO/plots/coefficients_distribution_depth_{hdo_depth}.png", hdo_depth = HDO_DEPTHS),
        annotations_per_depth_plot = "work_folder/data/ElasticNet/HDO/plots/annotations_per_depth.png"
    script:
        "plot_HDO_EN_coefficients.py"

rule fit_elastic_net_GO:
    input:
        INPUT_FOLDERS_GO,
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        OUTPUT_FOLDERS_GO
    script:
        "compute_GO_elastic_net_coefficients.py"
