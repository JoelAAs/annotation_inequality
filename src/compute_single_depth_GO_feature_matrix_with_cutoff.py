import pandas as pd
import numpy as np

aspect = snakemake.wildcards.aspect
cutoff = int(snakemake.wildcards.cutoff)
depth = int(snakemake.wildcards.depth)

input_df = snakemake.input.complete_annotations
input_bait_usage = snakemake.input.bait_usage
input_count_file = snakemake.input.count_df

output_matrix = snakemake.output.single_depth_feature_matrix_with_cutoff
output_stats = snakemake.output.cutoff_file  
print(f"Processing feature matrix for GO {aspect} depth {depth} with cutoff {cutoff}...")

print("Loading data...")

df = pd.read_csv(input_df, sep='\t')
df_copy = df.copy()
df_copy.fillna({'go_id': 'No_go_id', 'depth': -1}, inplace=True)

bait_usage = pd.read_csv(input_bait_usage, sep='\t')
count_df = pd.read_csv(input_count_file, sep='\t')

print("Data loaded!\n")

print("Finding genes in common...")

all_genes = df_copy['entrez_id'].unique()
bait_ids = bait_usage['entrez_id_bait'].unique()
df_copy = df_copy[df_copy['entrez_id'].isin(bait_ids)]
n_common_genes = df_copy['entrez_id'].nunique()
common_genes = df_copy['entrez_id'].unique()

print(f"{n_common_genes} common genes obtained!\n")

print(f"Filtering for depth {depth}...")

df_copy = df_copy[df_copy['depth'] == depth]

print(f"Keeping {df_copy['entrez_id'].nunique()} genes!\n")

print(f"Removing go ids with gene count < {cutoff}")

unique_goids_before = df_copy['go_id'].nunique()

valid_goids = count_df[count_df['count'] >= cutoff]['go_id'].unique()
df_copy = df_copy[df_copy['go_id'].isin(valid_goids)].copy()

unique_goids_after = df_copy['go_id'].nunique()

print(f"valid go ids after cutoff = {unique_goids_after}\n")

removed_goids = unique_goids_before - unique_goids_after

print(f"Removed go ids = {removed_goids}\n")

if unique_goids_before > 0:
    percentage_remaining = round((unique_goids_after / unique_goids_before) * 100, 2)
else:
    percentage_remaining = 0.0

print(f"Creating {aspect} depth {depth} cutoff {cutoff} feature matrix...")
if not df_copy.empty:
    df_copy['value'] = 1
    df_onehot = df_copy.pivot_table(
        index='entrez_id',
        columns='go_id',
        values='value',
        fill_value=0
    )
    df_onehot = df_onehot.reindex(common_genes, fill_value=0)
else:
    df_onehot = pd.DataFrame(index=common_genes)

df_onehot = df_onehot.reset_index().rename(columns={'index': 'entrez_id'})
df_onehot.to_csv(output_matrix, sep='\t', index=False)

print(f"unique genes after feature matrix and cutoff go ids removal after reindexing = {df_onehot['entrez_id'].nunique()}")

print(f"Unique go ids in the feature matrix = {len(df_onehot.columns) - 1}")

print(f"Feature matrix saved to {output_matrix}")

print(f"Saving {aspect} depth {depth} feature matrix...")

stats_df = pd.DataFrame({
    'cutoff': [cutoff],
    'removed_go_ids': [removed_goids],
    'remaining_percentage': [percentage_remaining]
})
stats_df.to_csv(output_stats, sep=':', index=False)

print(f"Individual stats saved to {output_stats}")
print(f"GO {aspect} Depth {depth} feature matrix with cutoff {cutoff} processed!")