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
        hdo_df = "work_folder/data/HDO/annotation_per_entrez_baits.csv"
    output:
        merged_df = "work_folder/data/correlaton/all_annotations.csv"
    script:
        "merge_dfs.py"