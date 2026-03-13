import pandas as pd 
import mygene
from pathlib import Path
import os

# Args
input_pod = snakemake.input.bp_frequencies
ontologies = {
    'BP': snakemake.output.annotation_df_BP,
    'MF': snakemake.output.annotation_df_MF,
    'CC': snakemake.output.annotation_df_CC
}

# Input
pod_df = pd.read_parquet(input_pod)

# Processing
mg = mygene.MyGeneInfo()

for ontology, outputfile in ontologies.items():
    # Query mygene
    results = mg.querymany(
        pod_df['entrez_id'].astype(str),
        scopes = 'entrezgene',
        fields = f'entrezgene,go.{ontology}.id',
        species = 'human',
        as_dataframe = True
    )

    # Save results, drop na values, and clean df
    go_df = results
    go_df['go_id'] = go_df[f'go.{ontology}'].apply(
        lambda x: [d['id'] for d in x] if isinstance(x, list) else []
    )
    go_df = go_df[['entrezgene', 'go_id']].explode('go_id')
    go_df = go_df.dropna()
    go_df = go_df.rename(columns = {
        'entrezgene': 'entrez_id',
    })

    # Count GO annotations for each entrez_id
    pod_df = pod_df.drop_duplicates(subset = ['entrez_id'])
    go_df = go_df.drop_duplicates(subset = ['entrez_id', 'go_id'])
    go_annot_count = go_df.groupby('entrez_id').size().reset_index(name = 'count')

    # Merge, including genes with zero annotations if present
    final_df = (
        pod_df[['entrez_id']]
        .astype(str)
        .merge(go_annot_count, on = 'entrez_id', how = 'left')
        .fillna({'count':0}).astype(int)
    )

    print(f'GO {ontology} annotation counts done!')

    final_df.to_csv(outputfile, sep = '\t', index = False)