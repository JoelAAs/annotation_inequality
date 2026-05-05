import pandas as pd
import mygene

evidence_dir = snakemake.input.evidence_folder
targets_dir = snakemake.input.target_folder
disease_dir = snakemake.input.disease_folder
outputdf = snakemake.output.first_annotation_dates

print(f"Assigning HDO disease-gene first association dates...")

mg = mygene.MyGeneInfo()

print(f"Loading and filtering HDO evidence data...")

df_evidence = pd.read_parquet(
    evidence_dir,
    columns = ['targetId', 'diseaseId', 'publicationDate']
)

df_evidence['publicationDate'] = pd.to_datetime(df_evidence['publicationDate'], errors = 'coerce')
df_evidence = df_evidence.dropna(subset = ['publicationDate'])

first_dates = df_evidence.groupby(['targetId', 'diseaseId'])['publicationDate'].min().reset_index()
first_dates = first_dates.rename(columns = {'publicationDate': 'first_publication_date'})
first_dates['first_publication_date'] = first_dates['first_publication_date'].dt.strftime('%Y%m%d')

print(f"Loading and filtering HDO target data...")

df_targets = pd.read_parquet(targets_dir, columns = ['id', 'approvedSymbol'])
df_targets = df_targets.rename(columns = {
    'id': 'targetId',
    'approvedSymbol': 'Gene_Symbol'
})

print(f"Loading and filtering HDO disease data...")

df_diseases = pd.read_parquet(disease_dir, columns = ['id', 'name'])
df_diseases = df_diseases.rename(columns = {
    'id': 'diseaseId',
    'name': 'Disease_Name'
})

print(f"Merging targets and diseases ids with dates...")

final_df = first_dates.merge(df_targets, on = 'targetId', how = 'left')
final_df = final_df.merge(df_diseases, on = 'diseaseId', how = 'left')

print(f"Translating Ensembl to Entrez IDs...")

# Retrieve entrez ids of the genes starting from ENSG format through myGene
ensembl_list = final_df['targetId'].dropna().unique().tolist()
results = mg.querymany(
    ensembl_list, 
    scopes='ensembl.gene', 
    fields='entrezgene', 
    species='human', 
    as_dataframe=True
)
results = results.reset_index()

# Create the mapping table to later merge, dropping duplicates mistakes
mapping_table = results[['query', 'entrezgene']].rename(columns={
    'query': 'targetId', 
    'entrezgene': 'entrez_targetId'  
})
mapping_table = mapping_table.drop_duplicates(subset=['targetId'])

final_df = final_df.merge(mapping_table, on='targetId', how='left')
final_df['entrez_targetId'] = pd.to_numeric(final_df['entrez_targetId'], errors='coerce').astype('Int64')

print(f"Standardizing disease IDs to DOID...")

# Use the dbXrefs attribute in the disease df to retrieve when possible the DOID format of the disease to later merge
df_diseases_xrefs = pd.read_parquet(disease_dir, columns=['id', 'dbXRefs'])
df_exploded = df_diseases_xrefs.explode('dbXRefs')
df_xrefs_doid = df_exploded[df_exploded['dbXRefs'].str.startswith('DOID:', na=False)]

mapping_series = df_xrefs_doid.drop_duplicates(subset=['id']).set_index('id')['dbXRefs']
final_df['standardized_DOID'] = final_df['diseaseId'].map(mapping_series)

is_already_doid = final_df['diseaseId'].str.startswith('DOID', na=False)
final_df.loc[is_already_doid, 'standardized_DOID'] = final_df.loc[is_already_doid, 'diseaseId'].str.replace('_', ':')

final_df = final_df[[
    'Gene_Symbol', 'targetId', 'entrez_targetId', 
    'Disease_Name', 'diseaseId', 'standardized_DOID', 
    'first_publication_date'
]]

output_df = final_df[['entrez_targetId', 'standardized_DOID', 'first_publication_date']]
output_df = output_df.dropna(subset = ['entrez_targetId', 'standardized_DOID'])

print(f"Saving HDO disease-gene first association dates...")

output_df.to_csv(outputdf, sep = '\t', index = False)

print(f"HDO disease-gene first association dates computed!")