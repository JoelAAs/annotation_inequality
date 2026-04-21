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

rule compute_GO_complete_elastic_net_coefficients:
    input:
        feature_matrix = "work_folder/data/GO/feature_matrices/complete/{aspect}_complete_feature_matrix.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/GO/EN_coefficients/complete/{aspect}_complete_elastic_net_coefficients.csv"
    script:
        "compute_GO_complete_elastic_net_coefficients.py"

rule visualize_complete_GO_coefficients:
    input: 
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/GO/EN_coefficients/complete/{aspect}_complete_elastic_net_coefficients.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/GO/plots/Top/complete/complete_{aspect}_top_coefficients.png",
        coefficients_distribution = "work_folder/data/ElasticNet/GO/plots/Distribution/complete/complete_{aspect}_coefficients_distribution.png"
    script: 
        "plot_GO_EN_coefficients.py"

rule plot_go_annotations_per_depth:
    input: 
        annotations_per_depth = "work_folder/data/GO/annotations_per_depth_{aspect}.csv"
    output: 
        annotations_per_depth_plot =  "work_folder/data/ElasticNet/GO/plots/{aspect}_annotations_per_depth.png"
    script: 
        "plot_GO_annotations_per_depth.py"

rule plot_go_genes_per_depth:
    input: 
        genes_per_depth = "work_folder/data/GO/genes_per_depth_{aspect}.csv"
    output: 
        genes_per_depth_plot =  "work_folder/data/ElasticNet/GO/plots/{aspect}_genes_per_depth.png"
    script: 
        "plot_go_genes_per_depth.py"

def get_all_en_coefficients(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep = ':')
    
    all_en_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        for d in range(max_d + 1):
            path = (
                f"work_folder/data/ElasticNet/GO/EN_coefficients/single_depth/"
                f"{asp}_depth_{d}_elastic_net_coefficients.csv"
            )
            all_en_paths.append(path)
            
    return all_en_paths

rule compute_GO_single_depth_elastic_net_coefficients:
    input:
        feature_matrix = "work_folder/data/GO/feature_matrices/single_depth/{aspect}_depth_{depth}_feature_matrix.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        single_depth_elastic_net_coefficients = "work_folder/data/ElasticNet/GO/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients.csv"
    script:
        "compute_GO_single_depth_elastic_net_coefficients.py"

def get_all_en_plots(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep = ':')
    all_plots = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        for d in range(max_d + 1):
            all_plots.append(
                f"work_folder/data/ElasticNet/GO/plots/Top/single_depth/{asp}_depth_{d}_top_coefficients.png"
            )
            all_plots.append(
                f"work_folder/data/ElasticNet/GO/plots/Distribution/single_depth/{asp}_depth_{d}_coefficients_distribution.png"
            )
            
    return all_plots

rule visualize_single_depth_GO_coefficients:
    input: 
        single_depth_elastic_net_coefficients = "work_folder/data/ElasticNet/GO/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/GO/plots/Top/single_depth/{aspect}_depth_{depth}_top_coefficients.png",
        coefficients_distribution = "work_folder/data/ElasticNet/GO/plots/Distribution/single_depth/{aspect}_depth_{depth}_coefficients_distribution.png"
    script: 
        "plot_GO_single_depth_EN_coefficients.py"