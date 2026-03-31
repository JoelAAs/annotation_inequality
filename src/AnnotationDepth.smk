ASPECTS = ['BP', 'MF', 'CC']

rule add_annotation_depth_HDO:
    input:
        annot_df = "work_folder/data/HDO/annotations_per_gene.csv"
    output:
        annot_df_depth = "work_folder/data/HDO/annotations_per_gene_with_depth.csv" 
    script:
        "get_HDO_annotations_depth.py"

rule add_annotation_depth_GO:
    input:
        annot_df = expand("work_folder/data/GO/{aspect}_annotations_per_gene.csv", aspect = ASPECTS)
    output:
        annot_df_depth = expand("work_folder/data/GO/{aspect}_annotations_per_gene_with_depth.csv", aspect = ASPECTS)
    script:
        "get_GO_annotations_depth.py"
