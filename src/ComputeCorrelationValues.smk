import pandas as pd

ASPECTS = glob_wildcards(
    "work_folder/data/GO/annotation_per_entrez_{aspect}_baits.csv"
).aspect

rule merge_annotation_dfs:
    input:
        go_dfs = expand(
            "work_folder/data/GO/annotation_per_entrez_{aspect}_baits.csv",
            aspect = ASPECTS
        ),
        hdo_df = "work_folder/data/HDO/annotation_per_entrez_baits.csv",
        disgenet_df = "work_folder/data/disgenet/annotation_per_entrez_baits.csv"
    output:
        merged_df = "work_folder/data/correlaton/all_annotations.csv"
    script:
        "merge_dfs.py"


rule get_correlation_matrix:
    input:
        merged_df = "work_folder/data/correlaton/all_annotations.csv"
    output:
        correlation_values = "work_folder/data/correlaton/corr_matrix_{method}.csv",
    params:
        method = "{method}"
    script:
        "compute_correlation_matrix.py"

rule plot_correlation_matrix:
    input:
        correlation_values = "work_folder/data/correlaton/corr_matrix_{method}.csv"
    output:
        correlation_matrix = "work_folder/data/correlaton/corr_matrix_{method}.png"
    script:
        "plot_correlation_matrix.py"