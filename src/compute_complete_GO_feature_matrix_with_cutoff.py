import pandas as pd

aspect = snakemake.wildcards.aspect
cutoff = int(snakemake.wildcards.cutoff)
input_df = snakemake.input.complete_annotations
cutoff_file = snakemake.output.cutoff_file
input_bait_usage = snakemake.input.bait_usage
input_count_file = snakemake.input.count_df
output_matrix = snakemake.output.complete_feature_matrix_with_cutoff

print(f"Processing complete GO {aspect} feature matrix with cutoff {cutoff}...\n")

print("Loading data...")

df = pd.read_csv(input_df, sep = '\t')
df_copy = df.copy()
df_copy.fillna({
    'go_id': 'No_go_id'
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
common_genes = df_copy['entrez_id'].unique()

print(f"{n_common_genes} common genes obtained!\n")

print(f"Removing go ids with gene count < {cutoff}...\n")

unique_go_ids_before = df_copy['go_id'].nunique()

print(f"unique genes before feature matrix = {df_copy['entrez_id'].nunique()}")

print(f"valid go ids = {df_copy['go_id'].nunique()}")

valid_go_ids = count_df[count_df['count'] >= cutoff]['go_id'].unique()
df_copy = df_copy[df_copy['go_id'].isin(valid_go_ids)].copy()
df_copy['value'] = 1

unique_go_ids_after = df_copy['go_id'].nunique()

print(f"valid go ids after cutoff = {unique_go_ids_after}\n")

print(f"Go {aspect} ids removed!\n")

# Removed doids
removed_go_ids = unique_go_ids_before - unique_go_ids_after

# Calculating the remaining doid percentage
if unique_go_ids_before > 0:
    percentage_remaining = (unique_go_ids_after / unique_go_ids_before) * 100
    percentage_remaining = round(percentage_remaining, 2)
else:
    percentage_remaining = 0

print("Creating feature matrix...")

df_onehot = df_copy.pivot_table(
    index = 'entrez_id',
    columns = 'go_id',
    values = 'value',
    fill_value = 0
)

print(f"unique genes after feature matrix and cutoff doid removal before reindexing = {df_onehot.index.nunique()}")

df_onehot = df_onehot.reindex(common_genes, fill_value = 0)
df_onehot = df_onehot.reset_index()

print(f"unique genes after feature matrix and cutoff doid removal after reindexing = {df_onehot['entrez_id'].nunique()}")

print(f"Unique go_ids in the feature matrix = {len(df_onehot.columns) - 1}")

print("Feature_matrix ready!\n")

print("Saving feature matrix...")

df_onehot.to_csv(output_matrix, sep = '\t', index = False)

print("Feature matrix saved!\n")

print(f"Creating cutoff {cutoff} file...")

stats_data = {
    'cutoff': [cutoff],
    'removed_go_ids': [removed_go_ids],
    'remaining_percentage': [percentage_remaining]
}

df_stats = pd.DataFrame(stats_data)
df_stats.to_csv(cutoff_file, sep = ':', index = False)

print(f"Cutoff {cutoff} file created!\n")

print(f"Complete GO {aspect} feature matrix with cutoff {cutoff} processed!")