rule add_HDO_ancestors:
    input: 
        annotation_df = "work_folder/data/HDO/annotations_per_gene.csv"
    output: 
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv"
    script:
        "add_HDO_ancestors.R"

rule count_genes_per_depth:
    input: 
        complete_annotations = "work_folder/data/HDO/annotations_per_gene_with_ancestors.csv"
    output: 
        genes_per_depth = "work_folder/data/HDO/genes_per_depth.csv"
    run: 
        import pandas as pd

        # Load the data
        df = pd.read_csv(input[0], sep='\t')
        
        # Drop rows where depth is NaN (like the 'No_doid' entries)
        df = df.dropna(subset=['depth'])
        
        # Count unique genes per depth
        # We use nunique() in case a gene appears multiple times at the same depth
        counts = df.groupby('depth')['entrez_id'].nunique()

        with open(output[0], 'w') as f:
            f.write("depth:n_of_genes\n")
        
        # Write to file in the format "depth:count"
        with open(output[0], 'a') as f:
            for depth, count in counts.items():
                # Formatting as depth:n_genes
                f.write(f"{int(depth)}:{count}\n")