ASPECTS = ['BP', 'MF', 'CC']
CUTOFFS = [5, 10, 15, 20, 25, 30, 40, 50, 100]

## Complete cutoff section
rule compute_complete_GO_feature_matrix_with_cutoff:
    input: 
        complete_annotations = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv",
        count_df = "work_folder/data/GO/annotations_counts/{aspect}_annotations_counts.csv"
    output: 
        complete_feature_matrix_with_cutoff = "work_folder/data/GO/cutoff/feature_matrices/complete/complete_{aspect}_feature_matrix_cutoff_{cutoff}.csv",
        cutoff_file = "work_folder/data/GO/cutoff/cutoff_files/complete/complete_{aspect}_cutoff_{cutoff}_file.csv"
    script:
        "compute_complete_GO_feature_matrix_with_cutoff.py"

rule aggregate_GO_complete_cutoff_files:
    input: 
        lambda wildcards: expand(
                    "work_folder/data/GO/cutoff/cutoff_files/complete/complete_{aspect}_cutoff_{cutoff}_file.csv", 
                    aspect = wildcards.aspect, 
                    cutoff = CUTOFFS
                )
    output: 
        cutoff_file = "work_folder/data/GO/cutoff/cutoff_files/complete/complete_{aspect}_cutoff_file.csv"
    run: 
        import pandas as pd

        all_frames = []

        for f in input:
            temp_df = pd.read_csv(f, sep = ':')
            all_frames.append(temp_df)

        combined = pd.concat(all_frames, ignore_index = True)
        combined = combined.sort_values('cutoff')
        combined.to_csv(output.cutoff_file, sep = ':', index= False)

rule aggregate_GO_complete_feature_matrices:
    input:
        expand("work_folder/data/GO/cutoff/feature_matrices/complete/complete_{aspect}_feature_matrix_cutoff_{cutoff}.csv", 
               aspect = ASPECTS, 
               cutoff = CUTOFFS
            )
    output:
        touch("work_folder/data/GO/cutoff/done_files/complete_matrices_done.txt")

rule compute_complete_GO_elastic_net_coefficients_with_cutoff:
    input: 
        wait = "work_folder/data/GO/cutoff/done_files/complete_matrices_done.txt",
        complete_matrix = "work_folder/data/GO/cutoff/feature_matrices/complete/complete_{aspect}_feature_matrix_cutoff_{cutoff}.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output: 
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/complete/complete_{aspect}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    script: 
        "compute_complete_GO_elastic_net_coefficients_with_cutoff.py"

rule aggregate_complete_GO_elastic_net_coefficients:
    input:
        expand("work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/complete/complete_{aspect}_elastic_net_coefficients_cutoff_{cutoff}.csv", 
               aspect = ASPECTS, 
               cutoff = CUTOFFS
            )
    output:
        touch("work_folder/data/GO/cutoff/done_files/complete_en_coefficients_done.txt")

rule visualize_complete_GO_elastic_net_coefficients_with_cutoff:
    input: 
        wait = "work_folder/data/GO/cutoff/done_files/complete_en_coefficients_done.txt",
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/complete/complete_{aspect}_elastic_net_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/GO_cutoff/plots/Top/complete/complete_{aspect}_top_coefficients_cutoff_{cutoff}.png",
        coefficients_distribution = "work_folder/data/ElasticNet/GO_cutoff/plots/Distribution/complete/complete_{aspect}_coefficients_distribution_cutoff_{cutoff}.png"
    script: 
        "plot_complete_GO_EN_coefficients_with_cutoff.py" 

rule plot_complete_GO_ids_lost_with_each_cutoff:
    input:
        cutoff_file = "work_folder/data/GO/cutoff/cutoff_files/complete/complete_{aspect}_cutoff_file.csv"
    output:
        cutoff_losing_plot = "work_folder/data/ElasticNet/GO_cutoff/plots/go_ids_lost/complete/complete_{aspect}_GO_ids_lost.png"
    script:
        "plot_complete_GO_ids_lost_with_each_cutoff.py"


## Single depth cutoff section

rule compute_single_depth_GO_feature_matrix_with_cutoff:
    input: 
        complete_annotations = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv",
        count_df = "work_folder/data/GO/annotations_counts/{aspect}_annotations_counts.csv"
    output: 
        single_depth_feature_matrix_with_cutoff = "work_folder/data/GO/cutoff/feature_matrices/single_depth/{aspect}_feature_matrix_depth_{depth}_cutoff_{cutoff}.csv",
        cutoff_file = "work_folder/data/GO/cutoff/cutoff_files/single_depth/{aspect}_depth_{depth}_cutoff_{cutoff}_file.csv"
    script:
        "compute_single_depth_GO_feature_matrix_with_cutoff.py"

def get_all_GO_depth_cutoff_paths(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    
    # Safety check for dry-runs
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep=':')
    all_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        # We need to expand across TWO variables: depth and cutoff
        # This nested loop creates the grid for this specific aspect
        for d in range(max_d + 1):
            for c in CUTOFFS:
                path = (
                    f"work_folder/data/GO/cutoff/feature_matrices/single_depth/"
                    f"{asp}_feature_matrix_depth_{d}_cutoff_{c}.csv"
                )
                all_paths.append(path)
                
    return all_paths

def get_all_GO_depth_summary_files(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep=':')
    summary_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        # This targets the aggregated file for EACH depth
        for d in range(max_d + 1):
            path = f"work_folder/data/GO/cutoff/cutoff_files/single_depth/{asp}_depth_{d}_cutoff_file.csv"
            summary_paths.append(path)
            
    return summary_paths

rule aggregate_GO_single_depth_cutoff_files:
    input:
        lambda wildcards: expand(
            "work_folder/data/GO/cutoff/cutoff_files/single_depth/{aspect}_depth_{depth}_cutoff_{cutoff}_file.csv", 
            aspect = wildcards.aspect, depth = wildcards.depth, cutoff = CUTOFFS
        )
    output:
        cutoff_file = "work_folder/data/GO/cutoff/cutoff_files/single_depth/{aspect}_depth_{depth}_cutoff_file.csv"
    run:
        import pandas as pd
        combined = pd.concat([pd.read_csv(f, sep=':') for f in input], ignore_index=True)
        combined.sort_values('cutoff').to_csv(output.cutoff_file, sep=':', index=False)

rule aggregate_GO_single_depth_feature_matrices:
    input:
        get_all_GO_depth_cutoff_paths
    output:
        touch("work_folder/data/GO/cutoff/done_files/single_depth_matrices_done.txt")

rule compute_single_depth_GO_elastic_net_coefficients_with_cutoff:
    input: 
        wait = "work_folder/data/GO/cutoff/done_files/single_depth_matrices_done.txt",
        single_depth_feature_matrix_with_cutoff = "work_folder/data/GO/cutoff/feature_matrices/single_depth/{aspect}_feature_matrix_depth_{depth}_cutoff_{cutoff}.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output: 
        single_depth_elastic_net_coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv",
        adj_r2_file = "work_folder/data/GO/cutoff/adj_r2_files/single_depth/{aspect}_depth_{depth}_cutoff_{cutoff}_adj_r2_file.csv"
    script: 
        "compute_single_depth_GO_elastic_net_coefficients_with_cutoff.py"

## TODO function to get all adj_r2_file based on aspect, also create function to aggregate them into one for each aspect-cutoff with all depths for that cutoff then plots as for HDO

def get_all_GO_single_depth_en_coefficients(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    
    # Safety check for dry-runs
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep=':')
    all_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        for d in range(max_d + 1):
            for c in CUTOFFS:
                path = (
                    f"work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/"
                    f"{asp}_depth_{d}_elastic_net_coefficients_cutoff_{c}.csv"
                )
                all_paths.append(path)
                
    return all_paths

rule aggregate_single_depth_GO_elastic_net_coefficients:
    input:
        get_all_GO_single_depth_en_coefficients
    output:
        touch("work_folder/data/GO/cutoff/done_files/single_depth_en_coefficients_done.txt")

def get_all_GO_single_depth_adj_r2_files(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep=':')
    all_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        for c in CUTOFFS:
            # This targets the final aggregated files for each aspect/cutoff
            path = f"work_folder/data/GO/cutoff/adj_r2_files/full/{asp}_cutoff_{c}_all_depths_adj_r2.csv"
            all_paths.append(path)
                
    return all_paths

def get_depth_files_for_aggregation(wildcards):
    """Filters the adj_r2 files for a specific aspect and cutoff across all depths."""
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    if not os.path.exists(depth_file): return []
    
    df = pd.read_csv(depth_file, sep=':')
    # Get max depth for the current aspect in wildcards
    max_d = int(df.loc[df['aspect'] == wildcards.aspect, 'max_depth'].values[0])
    
    return [
        f"work_folder/data/GO/cutoff/adj_r2_files/single_depth/{wildcards.aspect}_depth_{d}_cutoff_{wildcards.cutoff}_adj_r2_file.csv"
        for d in range(max_d + 1)
    ]

rule aggregate_GO_adj_r2_by_aspect_cutoff:
    input:
        files = get_depth_files_for_aggregation
    output:
        combined = "work_folder/data/GO/cutoff/adj_r2_files/full/{aspect}_cutoff_{cutoff}_all_depths_adj_r2.csv"
    run:
        import pandas as pd        
        # Load all depth-specific files
        # Using a list comprehension is efficient here
        df_list = [pd.read_csv(f) for f in input.files]        
        # Combine and sort by 'depth' so the resulting file is organized
        combined_df = pd.concat(df_list, ignore_index=True)        
        # Assuming your CSV has a column named 'depth'
        if 'depth' in combined_df.columns:
            combined_df = combined_df.sort_values('depth')            
        combined_df.to_csv(output.combined, index=False)


rule visualize_single_depth_GO_elastic_net_coefficients_with_cutoff:
    input: 
        wait = "work_folder/data/GO/cutoff/done_files/single_depth_en_coefficients_done.txt",
        single_depth_elastic_net_coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/GO_cutoff/plots/Top/single_depth/{aspect}_depth_{depth}_top_coefficients_cutoff_{cutoff}.png",
        coefficients_distribution = "work_folder/data/ElasticNet/GO_cutoff/plots/Distribution/single_depth/{aspect}_depth_{depth}_coefficients_distribution_cutoff_{cutoff}.png"
    script: 
        "plot_single_depth_GO_EN_coefficients_with_cutoff.py" 

def get_all_GO_single_depth_en_plots(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    
    # Safety check for dry-runs
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep=':')
    all_plot_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        for d in range(max_d + 1):
            for c in CUTOFFS:
                # Add both types of plots to the list
                top_path = (
                    f"work_folder/data/ElasticNet/GO_cutoff/plots/Top/single_depth/"
                    f"{asp}_depth_{d}_top_coefficients_cutoff_{c}.png"
                )
                dist_path = (
                    f"work_folder/data/ElasticNet/GO_cutoff/plots/Distribution/single_depth/"
                    f"{asp}_depth_{d}_coefficients_distribution_cutoff_{c}.png"
                )
                all_plot_paths.extend([top_path, dist_path])
                
    return all_plot_paths

rule plot_single_depth_GO_ids_lost_with_each_cutoff:
    input:
        cutoff_file = "work_folder/data/GO/cutoff/cutoff_files/single_depth/{aspect}_depth_{depth}_cutoff_file.csv"
    output:
        cutoff_losing_plot = "work_folder/data/ElasticNet/GO_cutoff/plots/go_ids_lost/single_depth/{aspect}_depth_{depth}_GO_ids_lost.png"
    script:
        "plot_single_depth_GO_ids_lost_with_each_cutoff.py"

def get_all_GO_lost_ids_plots(wildcards):
    depth_file = checkpoints.compute_GO_max_depths.get(**wildcards).output[0]
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep=':')
    plot_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        # One plot per aspect/depth combination
        for d in range(max_d + 1):
            path = f"work_folder/data/ElasticNet/GO_cutoff/plots/go_ids_lost/single_depth/{asp}_depth_{d}_GO_ids_lost.png"
            plot_paths.append(path)
            
    return plot_paths