ASPECTS = ['BP', 'MF', 'CC']

rule get_HDO_annotations:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        annotation_df = "work_folder/data/HDO/annotations_per_gene.csv",
        annotation_list = "work_folder/data/HDO/all_annotations.csv"
    script:
        "get_HDO_annotations_per_gene.R"

rule compute_HDO_feature_matrix:
    input:
        annotation_df = "work_folder/data/HDO/annotations_per_gene.csv",
        annotation_list = "work_folder/data/HDO/all_annotations.csv"
    output:
        feature_matrix = "work_folder/data/HDO/feature_matrix.csv"
    script:
        "compute_HDO_feature_matrix.py"

rule get_GO_annotations:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq"
    output:
        annotation_df = expand("work_folder/data/GO/{aspect}_annotations_per_gene.csv", aspect = ASPECTS),
        annotation_list = expand("work_folder/data/GO/{aspect}_all_annotations.csv", aspect = ASPECTS)
    script:
        "get_GO_annotations_per_gene.py"

rule compute_GO_feature_matrix:
    input:
        annotation_df = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        annotation_list = "work_folder/data/GO/{aspect}_all_annotations.csv"
    output:
        feature_matrix = "work_folder/data/GO/{aspect}_feature_matrix.csv"
    script:
        "compute_Go_feature_matrix.py"
        