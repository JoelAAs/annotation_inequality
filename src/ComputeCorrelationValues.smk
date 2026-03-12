rule get_spearman_correlation:
    input:
        input1 = "work_folder/data/{database}/annotation_per_entrez_baits.csv",
        input2 = "work_folder/data/{database}/annotation_per_entrez_preys.csv"  
    output:
        plot1 = "work_folder/data/plots/{database}_plots/baits_spearman.png",
        plot2 = "work_folder/data/plots/{database}_plots/preys_spearman.png"
    script:
        "compute_spearman.py"

rule get_pearson_correlation:
    input:
        input1 = "work_folder/data/{database}/annotation_per_entrez_baits.csv",
        input2 = "work_folder/data/{database}/annotation_per_entrez_preys.csv" 
    output:
        plot3 = "work_folder/data/plots/{database}_plots/baits_pearson.png",
        plot4 = "work_folder/data/plots/{database}_plots/preys_pearson.png"
    script:
        "compute_pearson.py"