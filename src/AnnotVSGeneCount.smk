rule plot_HDO_annotations_vs_gene_counts:
    input: 
        annotation_df = "work_folder/data/HDO/annotations_per_gene_with_depth.csv"
    output: 
        output_plot = "work_folder/data/plots/HDO_plots/annotations_vs_gene_counts_distrib.png",
        count_df = "work_folder/data/HDO/annotations_gene_counts.csv"
    script:
        "plot_HDO_annot_vs_gene_count.py"

rule plot_HDO_with_ancestors_annotations_vs_gene_counts:
    input: 
        annotation_df = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv"
    output: 
        output_plot = "work_folder/data/plots/HDO_plots/full_annotations_vs_gene_counts_distrib.png",
        count_df = "work_folder/data/HDO/full_annotations_gene_counts.csv"
    script:
        "plot_HDO_annot_vs_gene_count.py"