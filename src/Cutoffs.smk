CUTOFFS = [5, 10, 15, 20, 25, 30, 40, 50, 100]
HDO_DEPTHS_WITH_ANCESTORS = list(range(0, 13))

## All the cutoffs are made over the "complete" annotation dataframes, thus taking into account also the ancestor doids

rule create_complete_cutoff_file:
    output:
        cutoff_file = "work_folder/data/HDO/cutoff/cutoff_file.csv"
    shell:
        """
        echo "cutoff:removed_doids:remaining_percentage" > {output}
        """

rule compute_complete_HDO_feature_matrix_with_cutoff:
    input: 
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv",
        cutoff_file = "work_folder/data/HDO/cutoff/cutoff_file.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv",
        count_df = "work_folder/data/HDO/full_annotations_gene_counts.csv"
    output: 
        complete_feature_matrix_with_cutoff = "work_folder/data/HDO/cutoff/feature_matrices/complete_feature_matrix_cutoff_{cutoff}.csv"
    script:
        "compute_complete_HDO_feature_matrix_with_cutoff.py"

rule compute_complete_HDO_elastic_net_coefficients_with_cutoff:
    input: 
        complete_matrix = lambda wc: f"work_folder/data/HDO/cutoff/feature_matrices/complete_feature_matrix_cutoff_{wc.cutoff}.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output: 
        complete_elastic_net_coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete_elastic_net_coefficients_cutoff_{cutoff}.csv"
    script: 
        "compute_complete_HDO_elastic_net_coefficients_with_cutoff.py"

rule visualize_complete_HDO_elastic_net_coefficients_with_cutoff:
    input: 
        complete_elastic_net_coefficients = lambda wc: f"work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/complete_elastic_net_coefficients_cutoff_{wc.cutoff}.csv"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/HDO_cutoff/plots/Top/complete_top_coefficients_cutoff_{cutoff}.png",
        coefficients_distribution = "work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/complete_coefficients_distribution_cutoff_{cutoff}.png"
    script: 
        "plot_complete_HDO_EN_coefficients_with_cutoff.py" 

rule plot_complete_HDO_doids_lost_with_each_cutoff:
    input:
        cutoff_file = "work_folder/data/HDO/cutoff/cutoff_file.csv"
    output:
        cutoff_losing_plot = "work_folder/data/ElasticNet/HDO_cutoff/plots/complete_HDO_doids_lost.png"
    script:
        "plot_complete_HDO_doids_lost_with_each_cutoff.py"

rule create_single_depth_cutoff_file:
    output:
        cutoff_file = "work_folder/data/HDO/cutoff/cutoff_file_depth_{depth}.csv"
    shell:
        """
        echo "cutoff:removed_doids:remaining_percentage" > {output}
        """

rule compute_HDO_feature_matrix_with_cutoff_single_depth:
    input:
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv",
        cutoff_file = lambda wc: f"work_folder/data/HDO/cutoff/cutoff_file_depth_{wc.depth}.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv",
        count_df = "work_folder/data/HDO/full_annotations_gene_counts.csv"
    output:
        single_depth_feature_matrix_with_cutoff = "work_folder/data/HDO/cutoff/feature_matrices/feature_matrix_cutoff_{cutoff}_depth_{depth}.csv"
    script:
        "compute_single_depth_HDO_feature_matrix_with_cutoff.py"

rule compute_single_depth_HDO_elastic_net_coefficients_with_cutoff:
    input: 
        single_feature_matrix_with_cutoff = lambda wc: f"work_folder/data/HDO/cutoff/feature_matrices/feature_matrix_cutoff_{wc.cutoff}_depth_{wc.depth}.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output: 
        single_depth_elastic_net_coefficients = "work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv"
    script: 
        "compute_single_depth_HDO_elastic_net_coefficients_with_cutoff.py"

rule visualize_single_depth_HDO_elastic_net_coefficients_with_cutoff:
    input: 
        single_depth_elastic_net_coefficients = lambda wc: f"work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/depth_{wc.depth}_elastic_net_coefficients_cutoff_{wc.cutoff}.csv"
    output: 
        top_coefficients = "work_folder/data/ElasticNet/HDO_cutoff/plots/Top/depth_{depth}_top_coefficients_cutoff_{cutoff}.png",
        coefficients_distribution = "work_folder/data/ElasticNet/HDO_cutoff/plots/Distribution/depth_{depth}_coefficients_distribution_cutoff_{cutoff}.png"
    script: 
        "plot_single_depth_HDO_EN_coefficients_with_cutoff.py" 

rule plot_single_depth_HDO_doids_lost_with_each_cutoff:
    input:
        cutoff_file = lambda wc: f"work_folder/data/HDO/cutoff/cutoff_file_depth_{wc.depth}.csv"
    output:
        cutoff_losing_plot = "work_folder/data/ElasticNet/HDO_cutoff/plots/depth_{depth}_HDO_doids_lost.png"
    script:
        "plot_single_depth_HDO_doids_lost_with_each_cutoff.py"


         '''
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
               depth=HDO_DEPTHS_WITH_ANCESTORS)'''