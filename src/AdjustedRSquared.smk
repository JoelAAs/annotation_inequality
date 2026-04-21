rule plot_HDO_adjusted_R_squared_among_cutoffs:
    input: 
        full_adj_r2_file = "work_folder/data/HDO/cutoff/adj_r2_files/full/cutoff_{cutoff}_adj_r2_file.csv"
    output: 
        adj_r2_plot = "work_folder/data/HDO/cutoff/adj_r2_plots/cutoff_{cutoff}_adj_r2.png"
    script: 
        "pyScripts/plotting/plot_HDO_adj_R_squared.py"

rule plot_GO_adjusted_R_squared_among_cutoffs:
    input:
        full_adj_r2_file = "work_folder/data/GO/cutoff/adj_r2_files/full/{aspect}_cutoff_{cutoff}_all_depths_adj_r2.csv"
    output:
        adj_r2_plot = "work_folder/data/GO/cutoff/adj_r2_plots/{aspect}_cutoff_{cutoff}_adj_r2.png"
    script:
        "pyScripts/plotting/plot_GO_adj_R_squared.py"