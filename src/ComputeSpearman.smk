rule get_spearman_correlation:
    input:
        input1 = "work_folder/data/joined/interactions_annotations_baits.pq",
        input2 = "work_folder/data/joined/interactions_annotations_preys.pq"   
    output:
        plot1 = "work_folder/data/plots/interactions_annotations_baits_spearman.png",
        plot2 = "work_folder/data/plots/interactions_annotations_preys_spearman.png"
    script:
        "compute_spearman.py"