import pandas as pd

cutoff = int(snakemake.wildcards.cutoff)
depth = int(snakemake.wildcards.depth)
input_df = snakemake.input.complete_annotations
cutoff_file = snakemake.input.cutoff_file
input_bait_usage = snakemake.input.bait_usage
input_count_file = snakemake.input.count_df
output_matrix = snakemake.output.single_depth_feature_matrix_with_cutoff

print(f"Processing full depth {depth} feature matrix with cutoff {cutoff}...\n")

print("Loading data...")

df = pd.read_csv(input_df, sep = '\t')
df_copy = df.copy()
df_copy.fillna({
    'doid': 'No_doid',
    'depth': -1
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

print(f"Keeping only doids with depth {depth}...")

df_copy = df_copy[df_copy['depth'] == depth]

print(f"Doids not with depth {depth} removed!\n")

print(f"Removing doids with gene count < {cutoff}...\n")

unique_doids_before = df_copy['doid'].nunique()

print(f"unique genes before feature matrix = {df_copy['entrez_id'].nunique()}")

print(f"valid doids = {df_copy['doid'].nunique()}")

valid_doids = count_df[count_df['gene_count'] >= cutoff]['doid'].unique()
df_copy = df_copy[df_copy['doid'].isin(valid_doids)].copy()

print(f"Doids with gene count < {cutoff} removed!\n")

df_copy['value'] = 1
unique_doids_after = df_copy['doid'].nunique()

print(f"valid doids after cutoff = {unique_doids_after}\n")

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

print("Updating cutoff file...")

'''new_row = {
    'cutoff': cutoff,
    'removed_doids': removed_doids,
    'remaining_percentage': percentage_remaining
}

log_df = pd.read_csv(cutoff_file, sep = ':')
log_df = log_df[log_df['cutoff'] != cutoff]
log_df = pd.concat([log_df, pd.DataFrame([new_row])], ignore_index = True)
target_columns = ['cutoff', 'removed_doids', 'remaining_percentage']
log_df = log_df[target_columns]
log_df = log_df.sort_values('cutoff')
log_df.to_csv(cutoff_file, sep = ':', index = False)'''

import fcntl  # Standard library for file locking on Linux

# ... (your previous code) ...

# Use a context manager to handle the file safely
with open(cutoff_file, 'a+') as f:
    # Acquire an EXCLUSIVE lock. If another job has the lock, this job will WAIT here.
    fcntl.flock(f, fcntl.LOCK_EX)
    
    try:
        # 1. Seek to the start of the file to read existing data
        f.seek(0)
        content = f.read()
        
        if not content.strip():
            # Handle the case where the file might be truly empty (first run)
            log_df = pd.DataFrame(columns=['cutoff', 'removed_doids', 'remaining_percentage'])
        else:
            from io import StringIO
            log_df = pd.read_csv(StringIO(content), sep=':')

        # 2. Update the data
        new_row = pd.DataFrame([{
            'cutoff': cutoff,
            'removed_doids': removed_doids,
            'remaining_percentage': percentage_remaining
        }])
        
        log_df = log_df[log_df['cutoff'] != cutoff]
        log_df = pd.concat([log_df, new_row], ignore_index=True)
        log_df = log_df.sort_values('cutoff')

        # 3. Wipe the file and write the updated DataFrame
        f.seek(0)
        f.truncate()
        log_df.to_csv(f, sep=':', index=False)
        
    finally:
        # 4. Release the lock so the next waiting job can proceed
        fcntl.flock(f, fcntl.LOCK_UN)

print("Cutoff file updated!\n")

print(f"Complete full feature matrix with cutoff {cutoff} processed!")