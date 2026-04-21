CUTOFFS = [5, 10, 15, 20, 25, 30, 40, 50, 100]
HDO_DEPTHS_WITH_ANCESTORS = list(range(0, 13))

## HDO SECTION

rule merge_HDO_single_depth_coefficients_files:
    input: 
        coefficients = lambda wildcards: expand("work_folder/data/ElasticNet/HDO_cutoff/EN_coefficients/single_depth/depth_{depth}_elastic_net_coefficients_cutoff_{cutoff}.csv",
            depth = HDO_DEPTHS_WITH_ANCESTORS,
            cutoff = wildcards.cutoff
        )
    output: 
        all_coefficients = "work_folder/data/dendrograms/HDO/all_coefficients/all_coefficients_cutoff_{cutoff}.csv"
    run: 
        import pandas as pd

        print(f"Merging all HDO single depth coefficients for cutoff {wildcards.cutoff}...")

        df_list = [pd.read_csv(f, sep = '\t') for f in input.coefficients]

        merged_df = pd.concat(df_list, ignore_index = True)
        merged_df = merged_df.drop_duplicates(subset = ['HDO_doid'])
        merged_df = merged_df.sort_values(by = 'HDO_doid')

        merged_df.to_csv(output.all_coefficients, sep = '\t', index = False)

        print(f"All coefficients file for HDO cutoff {wildcards.cutoff} ready!\n")

rule create_HDO_dendrogram:
    input: 
        wait = "work_folder/data/HDO/cutoff/done_files/complete_en_coefficients_done.txt",
        all_coefficients = "work_folder/data/dendrograms/HDO/all_coefficients/all_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output: 
        dendrogram = "work_folder/data/dendrograms/HDO/visualization/dendrogram/dendrogram_cutoff_{cutoff}.pdf"
    script: 
        "create_HDO_dendrogram.py"

rule create_HDO_treemap:
    input:
        wait = "work_folder/data/HDO/cutoff/done_files/complete_en_coefficients_done.txt",
        all_coefficients = "work_folder/data/dendrograms/HDO/all_coefficients/all_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output: 
        treemap = "work_folder/data/dendrograms/HDO/visualization/treemap/treemap_cutoff_{cutoff}.pdf",
        treemap_html = "work_folder/data/dendrograms/HDO/visualization/treemap/treemap_cutoff_{cutoff}.html"
    script: 
        "create_HDO_treemap.py"

## GO SECTION

def get_GO_single_depth_coefficients_files_by_aspect(wildcards):
    checkpoint_output = checkpoints.compute_GO_max_depths.get(**wildcards).output.max_depth_file
    
    import pandas as pd
    df_depths = pd.read_csv(checkpoint_output, sep=':')
    
    aspect_row = df_depths[df_depths['aspect'] == wildcards.aspect]
    
    if aspect_row.empty:
        return []

    max_d = int(aspect_row['max_depth'].iloc[0])
    
    all_files = []
    for d in range(max_d + 1):
        file_path = (
            f"work_folder/data/ElasticNet/GO_cutoff/EN_coefficients/single_depth/"
            f"{wildcards.aspect}_depth_{d}_elastic_net_coefficients_cutoff_{wildcards.cutoff}.csv"
        )
        all_files.append(file_path)
            
    return all_files

rule merge_GO_single_depth_coefficients_files:
    input: 
        wait = "work_folder/data/GO/cutoff/done_files/complete_en_coefficients_done.txt",
        coefficients = get_GO_single_depth_coefficients_files_by_aspect
    output: 
        all_coefficients = "work_folder/data/dendrograms/GO/all_coefficients/{aspect}_all_coefficients_cutoff_{cutoff}.csv"
    run: 
        import pandas as pd

        print(f"Merging all GO {wildcards.aspect} single depth coefficients for cutoff {wildcards.cutoff}...")

        df_list = [pd.read_csv(f, sep = '\t') for f in input.coefficients]

        merged_df = pd.concat(df_list, ignore_index = True)
        
        id_col = 'GO_id' 
        merged_df = merged_df.drop_duplicates(subset = [id_col])
        merged_df = merged_df.sort_values(by = id_col)

        merged_df.to_csv(output.all_coefficients, sep = '\t', index = False)

        print(f"All coefficients file for GO cutoff {wildcards.cutoff} ready!\n")

rule create_GO_dendrogram:
    input: 
        wait = "work_folder/data/GO/cutoff/done_files/complete_en_coefficients_done.txt",
        all_coefficients = "work_folder/data/dendrograms/GO/all_coefficients/{aspect}_all_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        dendrogram = "work_folder/data/dendrograms/GO/visualization/dendrogram/{aspect}_dendrogram_cutoff_{cutoff}.pdf"
    script: 
        "pyScripts/dendrogram/create_GO_dendrogram.py"

rule create_GO_treemap:
    input:
        wait = "work_folder/data/GO/cutoff/done_files/complete_en_coefficients_done.txt",
        all_coefficients = "work_folder/data/dendrograms/GO/all_coefficients/{aspect}_all_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/GO/go-basic.obo"
    output: 
        treemap = "work_folder/data/dendrograms/GO/visualization/treemap/{aspect}_treemap_cutoff_{cutoff}.pdf",
        treemap_html = "work_folder/data/dendrograms/GO/visualization/treemap/{aspect}_treemap_cutoff_{cutoff}.html"
    script: 
        "pyScripts/dendrogram/create_GO_treemap.py"