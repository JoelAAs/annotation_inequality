CUTOFFS = [5, 10, 15, 20, 25, 30, 40, 50, 100]
HDO_DEPTHS_WITH_ANCESTORS = list(range(0, 13))

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
        all_coefficients = "work_folder/data/dendrograms/HDO/all_coefficients/all_coefficients_cutoff_{cutoff}.csv",
        ontology = "work_folder/data/HDO/doid.obo"
    output: 
        dendrogram = "work_folder/data/dendrograms/HDO/visualization/dendrogram_cutoff_{cutoff}.png"
    script: 
        "create_HDO_dendrogram.py"