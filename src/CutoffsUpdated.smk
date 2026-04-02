CUTOFFS = [5, 10, 15, 20, 25, 30, 40, 50, 100]
HDO_DEPTHS_WITH_ANCESTORS = list(range(0, 13))

## All the cutoffs are made over the "complete" annotation dataframes, thus taking into account also the ancestor doids

rule compute_complete_HDO_feature_matrix_with_cutoff:
    input: 
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv",
        count_df = "work_folder/data/HDO/full_annotations_gene_counts.csv"
    output: 
        complete_feature_matrix_with_cutoff = "work_folder/data/HDO/cutoff/feature_matrices/complete/complete_feature_matrix_cutoff_{cutoff}.csv",
        cutoff_file = "work_folder/data/HDO/cutoff/cutoff_files/complete/complete_cutoff_{cutoff}_file.csv"
    script:
        "compute_complete_HDO_feature_matrix_with_cutoff.py"


rule aggregate_complete_cutoff_files:
    input: 
        expand("work_folder/data/HDO/cutoff/cutoff_files/complete/complete_cutoff_{cutoff}_file.csv", cutoff = CUTOFFS)
    output: 
        cutoff_file = "work_folder/data/HDO/cutoff/cutoff_files/complete/complete_cutoff_file.csv"
    run: 
        import pandas as pd

        all_frames = []

        for f in input:
            temp_df = pd.read_csv(f, sep = ':')
            all_frames.append(temp_df)

        combined = pd.concat(all_frames, ignore_index = True)
        combined = combined.sort_values('cutoff')
        combined.to_csv(output.cutoff_file, sep = ':', index= False)

rule aggregate_complete_feature_matrices:
    input:
        [f"work_folder/data/HDO/cutoff/feature_matrices/complete/complete_feature_matrix_cutoff_{cutoff}.csv" for cutoff in CUTOFFS]
    output:
        touch("work_folder/data/HDO/cutoff/done_files/matrices_done.txt")


rule compute_complete_HDO_elastic_net_coefficients_with_cutoff:
    input: 
        wait = "work_folder/data/HDO/cutoff/done_files/matrices_done.txt",
        complete_matrix = "work_folder/data/HDO/cutoff/feature_matrices/complete/complete_feature_matrix_cutoff_{cutoff}.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output: 
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete/complete_elastic_net_coefficients_cutoff_{cutoff}.csv"
    script: 
        "compute_complete_HDO_elastic_net_coefficients_with_cutoff.py"

rule aggregate_complete_elastic_net_coefficients:
    input:
        [f"work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete/complete_elastic_net_coefficients_cutoff_{cutoff}.csv" for cutoff in CUTOFFS]
    output:
        touch("work_folder/data/HDO/cutoff/done_files/en_coefficients_done.txt")

rule visualize_complete_HDO_elastic_net_coefficients_with_cutoff:
    input: 
        wait = "work_folder/data/HDO/cutoff/done_files/en_coefficients_done.txt",
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete/complete_elastic_net_coefficients_cutoff_{cutoff}.csv"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/HDO_cutoff/plots/Top/complete/complete_top_coefficients_cutoff_{cutoff}.png",
        coefficients_distribution = "work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/complete/complete_coefficients_distribution_cutoff_{cutoff}.png"
    script: 
        "plot_complete_HDO_EN_coefficients_with_cutoff.py" 

rule plot_complete_HDO_doids_lost_with_each_cutoff:
    input:
        cutoff_file = "work_folder/data/HDO/cutoff/cutoff_files/complete/complete_cutoff_file.csv"
    output:
        cutoff_losing_plot = "work_folder/data/ElasticNet/HDO_cutoff/plots/doids_lost/complete_HDO_doids_lost.png"
    script:
        "plot_complete_HDO_doids_lost_with_each_cutoff.py"