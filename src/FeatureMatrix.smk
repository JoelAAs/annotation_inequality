import pandas as pd
import os

ASPECTS = ['BP', 'MF', 'CC']
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

## HDO SECTION      

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
        

## GO SECTION
rule download_go_ontology:
    output:
        "work_folder/data/GO/go-basic.obo"
    shell:
        "wget http://purl.obolibrary.org/obo/go/go-basic.obo -O {output}"

rule get_GO_annotations:
    input:
        bp_frequencies = "work_folder/data/intact/bait_prey_frequencies.pq",
        go_obo = "work_folder/data/GO/go-basic.obo" 
    output:
        annotation_df = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        annotation_list = "work_folder/data/GO/{aspect}_all_annotations.csv"
    script:
        "get_GO_annotations_per_gene.py"

rule compute_go_annotations_counts:
    input:
        annotation_df = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        counts_per_annot = "work_folder/data/GO/annotations_counts/{aspect}_annotations_counts.csv"
    script:
        "compute_go_annotations_counts.py"

rule compute_go_annotations_per_depth:
    input: 
        annotation_df = "work_folder/data/GO/{aspect}_annotations_per_gene.csv"
    output: 
        annotations_per_depth = "work_folder/data/GO/annotations_per_depth_{aspect}.csv"
    run:
        import pandas as pd

        df = pd.read_csv(input.annotation_df, sep='\t')
        # Ensure depth is treated consistently (as float to match your example "0.0")
        counts = df.dropna().groupby('depth')['go_id'].nunique()

        with open(output.annotations_per_depth, "w") as f:
            f.write("depth:n_of_coefficients\n")
            for depth, val in counts.items():
                # Using f-string to ensure the 0.0 format if desired
                f.write(f"{float(depth)}:{val}\n")
        
        print(f"GO {wildcards.aspect} annotations per depth saved!")

rule compute_go_genes_per_depth:
    input:
        annotation_df = "work_folder/data/GO/{aspect}_annotations_per_gene.csv"
    output:
        genes_per_depth = "work_folder/data/GO/genes_per_depth_{aspect}.csv"
    run:
        import pandas as pd

        df = pd.read_csv(input.annotation_df, sep='\t')
        gene_counts = df.dropna().groupby('depth')['entrez_id'].nunique()

        with open(output.genes_per_depth, "w") as f:
            f.write("depth:n_of_genes\n")
            for depth, val in gene_counts.items():
                f.write(f"{float(depth)}:{val}\n")
        
        print(f"GO {wildcards.aspect} genes per depth saved!")


rule compute_GO_complete_feature_matrix:
    input:
        annotation_df = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        feature_matrix = "work_folder/data/GO/feature_matrices/complete/{aspect}_complete_feature_matrix.csv"
    script:
        "compute_GO_complete_feature_matrix.py"

rule compute_GO_max_depth:
    input:
        aspect_files = expand("work_folder/data/GO/{aspect}_all_annotations.csv", aspect=ASPECTS)
    output:
        max_depth_file = "work_folder/data/GO/max_depths_file.csv"
    run:
        import pandas as pd
        import os

        results = []
        
        # Iterate through the list of 3 files (BP, MF, CC)
        for filepath in input.aspect_files:
            # Extract the aspect name from the filename (e.g., 'BP')
            aspect_name = os.path.basename(filepath).split('_')[0]
            
            # Read the data
            df = pd.read_csv(filepath, sep='\t')
            
            # Calculate max depth if the file isn't empty
            if not df.empty and 'depth' in df.columns:
                max_d = int(df['depth'].max())
            else:
                max_d = 0
            
            results.append({
                'aspect': aspect_name,
                'max_depth': max_d
            })
            print(f"Computed max depth for {aspect_name}: {max_d}")

        # Create the summary dataframe
        summary_df = pd.DataFrame(results)
        
        # Save to the output path defined in the rule
        summary_df.to_csv(output.max_depth_file, sep = ':', index=False)

def get_all_single_depth_matrices(wildcards):
    depth_file = "work_folder/data/GO/max_depths_file.csv"
    
    # Safety check: if the file isn't there yet (e.g., first dry-run)
    if not os.path.exists(depth_file):
        return []

    df = pd.read_csv(depth_file, sep = ':')
    
    all_outputs = []
    for _, row in df.iterrows():
        asp = row['aspect']      
        max_d = int(row['max_depth'])
        
        # Request files from depth 0 up to max_depth
        for d in range(max_d + 1):
            all_outputs.append(
                f"work_folder/data/GO/feature_matrices/single_depth/{asp}_depth_{d}_feature_matrix.csv"
            )
            
    return all_outputs

rule compute_GO_single_depth_feature_matrix:
    input:
        annotation_df = "work_folder/data/GO/{aspect}_annotations_per_gene.csv",
        bait_usage = "work_folder/data/intact/bait_count.csv"
    output:
        feature_matrix = "work_folder/data/GO/feature_matrices/single_depth/{aspect}_depth_{depth}_feature_matrix.csv"
    script:
        "compute_GO_single_depth_feature_matrix.py"