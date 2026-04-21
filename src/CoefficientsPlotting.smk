rule plot_HDO_highest_abs_value_coefficients_per_cutoff:
    input: 
        coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output: 
        top_abs_coeffs_plot = "work_folder/data/ElasticNet/HDO_cutoff/plots/Top_abs/top_abs_depth_{depth}_cutoff_{cutoff}.png"
    script: 
        "plot_HDO_highest_abs_value_coefficients.py"

## TODO get all highest abs value coefficients

rule plot_GO_highest_abs_value_coefficients_per_cutoff:
    input:
        coefficients = "work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/{aspect}_depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output:
        top_abs_coeffs_plot = "work_folder/data/ElasticNet/GO_cutoff/plots/Top_abs/top_abs_{aspect}_depth_{depth}_cutoff_{cutoff}.png"
    script:
        "pyScripts/plotting/plot_GO_highest_abs_value_coefficients.py"

def get_all_GO_top_abs_plots(wildcards):
    depth_file = "work_folder/data/GO/max_depths_file.csv"
    
    # Safety check for dry-runs or first-time runs
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep=':')
    all_plot_paths = []
    
    for _, row in df.iterrows():
        asp = row['aspect']
        max_d = int(row['max_depth'])
        
        # We loop from 0 to the specific max_depth for this aspect
        for d in range(max_d + 1):
            for c in CUTOFFS:
                # Mirroring your rule's output path structure
                path = (
                    f"work_folder/data/ElasticNet/GO_cutoff/plots/Top_abs/"
                    f"top_abs_{asp}_depth_{d}_cutoff_{c}.png"
                )
                all_plot_paths.append(path)
                
    return all_plot_paths