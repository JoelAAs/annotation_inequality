import pandas as pd 
import mygene
from pathlib import Path

# Args
input_pod = snakemake.input.bp_frequencies
output_pod = snakemake.output.annotation_df_BP

# Input
pod_df = pd.read_parquet(input_pod)

# Processing
mg = mygene.MyGeneInfo()
print('Here!!')

# Query mygene
results = mg.querymany(
    pod_df['entrez_id'].astype(str),
    scopes = 'entrezgene',
    fields = 'entrezgene,go.BP.id',
    species = 'human',
    as_dataframe = True
)

# Save results, drop na values, and clean df
go_df = results[['entrezgene', 'go.BP.id']].explode('go.BP.id')
go_df = go_df.dropna()
go_df = go_df.rename(columns = {
    'entrezgene': 'entrez_id',
    'go.BP.id': 'go_id'
})

# Count GO annotations for each entrez_id
go_annot_count = go_df.groupby('entrez_id').size().reset_index(name = 'go_count')

# Merge, including genes with zero annotations if present
final_df = (
    pod_df[['entrez_id']]
    .astype(str)
    .merge(go_annot_count, on = 'entrez_id', how = 'right')
    .fillna({'go_count':0})
)

# Ensure parent folder exists
Path(output_pod).parent.mkdir(parents = True, exist_ok = True)
print('Final dataframe shape:', final_df.shape)

final_df.to_csv(output_pod, index = False)