rule get_spearman_correlation_plots:
    input:
        input1 = "work_folder/data/{database}/annotation_per_entrez_baits.csv",
        input2 = "work_folder/data/{database}/annotation_per_entrez_preys.csv"  
    output:
        plot1 = "work_folder/data/plots/{database}_plots/baits_spearman.png",
        plot2 = "work_folder/data/plots/{database}_plots/preys_spearman.png"
    script:
        "compute_spearman.py"

rule get_pearson_correlation_plots:
    input:
        input1 = "work_folder/data/{database}/annotation_per_entrez_baits.csv",
        input2 = "work_folder/data/{database}/annotation_per_entrez_preys.csv" 
    output:
        plot1 = "work_folder/data/plots/{database}_plots/baits_pearson.png",
        plot2 = "work_folder/data/plots/{database}_plots/preys_pearson.png"
    script:
        "compute_pearson.py"

rule get_spearman_correlation_plots_GO:
    input:
        input1 = "work_folder/data/GO/annotation_per_entrez_{aspect}_baits.csv",
        input2 = "work_folder/data/GO/annotation_per_entrez_{aspect}_preys.csv"  
    output:
        plot1 = "work_folder/data/plots/GO_plots/{aspect}_baits_spearman.png",
        plot2 = "work_folder/data/plots/GO_plots/{aspect}_preys_spearman.png"
    script:
        "compute_spearman.py"

rule get_pearson_correlation_plots_GO:
    input:
        input1 = "work_folder/data/GO/annotation_per_entrez_{aspect}_baits.csv",
        input2 = "work_folder/data/GO/annotation_per_entrez_{aspect}_preys.csv" 
    output:
        plot1 = "work_folder/data/plots/GO_plots/{aspect}_baits_pearson.png",
        plot2 = "work_folder/data/plots/GO_plots/{aspect}_preys_pearson.png"
    script:
        "compute_pearson.py"