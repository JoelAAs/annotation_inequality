ASPECTS = ['BP', 'MF', 'CC']
INPUT_FILES_GO = [f"work_folder/data/GO/{aspect}_annotations_per_gene_with_depth.csv" for aspect in ASPECTS]
OUTPUT_FOLDERS_GO = [directory(f"work_folder/data/GO/feature_matrices/{aspect}") for aspect in ASPECTS]
HDO_DEPTHS = list(range(0, 11))
HDO_DEPTHS_WITH_ANCESTORS = list(range(0, 13))

rule get_disgenet_annotations:
    input:
        bait_ids = "work_folder/data/disgenet/annotation_per_entrez_baits.csv"
    output:
        annotation_df = "work_folder/data/disgenet/annotations_per_gene.csv",
        annotation_list = "work_folder/data/disgenet/all_annotations.csv"
    script:
        "get_disgenet_annotations_per_gene.py"

rule compute_disgenet_feature_matrix:
    input: 
        annotation_df = "work_folder/data/disgenet/annotations_per_gene.csv"
    output:
        feature_matrix = "work_folder/data/disgenet/feature_matrix.csv"
    script:
        "compute_disgenet_feature_matrix.py"
        
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
        annotation_df = "work_folder/data/HDO/annotations_per_gene_with_depth.csv"
    output:
        feature_matrix = expand("work_folder/data/HDO/feature_matrices/feature_matrix_with_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS),
        annotations_per_depth = "work_folder/data/HDO/annotations_per_depth.csv",
        counts_per_annot = expand("work_folder/data/HDO/annot_counts_per_depth/annot_counts_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS)
    script:
        "compute_HDO_feature_matrix.py"

rule compute_HDO_feature_matrix_with_ancestors:
    input:
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv"
    output:
        feature_matrix = expand("work_folder/data/HDO/feature_matrices/full_feature_matrix_with_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS_WITH_ANCESTORS),
        annotations_per_depth = "work_folder/data/HDO/full_annotations_per_depth.csv",
        counts_per_annot = expand("work_folder/data/HDO/annot_counts_per_depth/full_annot_counts_depth_{hdo_depth}.csv", hdo_depth = HDO_DEPTHS_WITH_ANCESTORS)
    script:
        "compute_HDO_feature_matrix_with_ancestors.py"

rule compute_complete_HDO_feature_matrix_with_ancestors:
    input:
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv"
    output:
        complete_matrix = "work_folder/data/HDO/feature_matrices/complete_matrix_with_ancestors.csv"
    script:
        "compute_complete_HDO_feature_matrix_with_ancestors.py"

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
        INPUT_FILES_GO
    output:
        feature_dirs = OUTPUT_FOLDERS_GO,
        annotations_per_depth = expand("work_folder/data/GO/annotations_per_depth_{aspect}.csv", aspect = ASPECTS),
        min_aspect_depths = expand("work_folder/data/GO/{aspect}_min_aspect_depths.csv", aspect = ASPECTS)
    script:
        "compute_GO_feature_matrices.py"
        