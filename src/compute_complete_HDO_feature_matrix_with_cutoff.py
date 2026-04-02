import pandas as pd

cutoff = int(snakemake.wildcards.cutoff)
input_df = snakemake.input.complete_annotations
cutoff_file = snakemake.output.cutoff_file
input_bait_usage = snakemake.input.bait_usage
input_count_file = snakemake.input.count_df
output_matrix = snakemake.output.complete_feature_matrix_with_cutoff

print(f"Processing complete full feature matrix with cutoff {cutoff}...\n")

print("Loading data...")

df = pd.read_csv(input_df, sep = '\t')
df_copy = df.copy()
df_copy.fillna({
    'doid': 'No_doid'
}, inplace=True)

bait_usage = pd.read_csv(input_bait_usage, sep = '\t')

count_df = pd.read_csv(input_count_file, sep = '\t')

print("Data loaded!\n")

print("Merging to obtain common genes...")

all_genes = df_copy['entrez_id'].unique()
print(f"Total gene universe: {len(all_genes)}")

bait_ids = bait_usage['entrez_id_bait'].unique()
df_copy = df_copy[df_copy['entrez_id'].isin(bait_ids)]
n_common_genes = df_copy['entrez_id'].nunique()

print(f"{n_common_genes} common genes obtained!\n")

print(f"Removing doids with gene count < {cutoff}...\n")

unique_doids_before = df_copy['doid'].nunique()

print(f"unique genes before feature matrix = {df_copy['entrez_id'].nunique()}")

print(f"valid doids = {df_copy['doid'].nunique()}")

valid_doids = count_df[count_df['gene_count'] >= cutoff]['doid'].unique()
df_copy = df_copy[df_copy['doid'].isin(valid_doids)].copy()
df_copy['value'] = 1

unique_doids_after = df_copy['doid'].nunique()

print(f"valid doids after cutoff = {unique_doids_after}\n")

print("Doids removed!\n")

# Removed doids
removed_doids = unique_doids_before - unique_doids_after

# Calculating the remaining doid percentage
if unique_doids_before > 0:
    percentage_remaining = (unique_doids_after / unique_doids_before) * 100
    percentage_remaining = round(percentage_remaining, 2)
else:
    percentage_remaining = 0

print("Creating feature matrix...")

df_onehot = df_copy.pivot_table(
    index = 'entrez_id',
    columns = 'doid',
    values = 'value',
    fill_value = 0
)

print(f"unique genes after feature matrix and cutoff doid removal before reindexing = {df_onehot.index.nunique()}")

df_onehot = df_onehot.reindex(all_genes, fill_value = 0)
df_onehot = df_onehot.reset_index()

print(f"unique genes after feature matrix and cutoff doid removal after reindexing = {df_onehot['entrez_id'].nunique()}")

print(f"Unique doids in the feature matrix = {len(df_onehot.columns) - 1}")

print("Feature_matrix ready!\n")

print("Saving feature matrix...")

df_onehot.to_csv(output_matrix, sep = '\t', index = False)

print("Feature matrix saved!\n")

print(f"Creating cutoff {cutoff} file...")

stats_data = {
    'cutoff': [cutoff],
    'removed_doids': [removed_doids],
    'remaining_percentage': [percentage_remaining]
}

df_stats = pd.DataFrame(stats_data)
df_stats.to_csv(cutoff_file, sep = ':', index = False)

print(f"Cutoff {cutoff} file created!\n")

print(f"Complete full feature matrix with cutoff {cutoff} processed!")