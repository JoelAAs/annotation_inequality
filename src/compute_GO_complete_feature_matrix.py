import pandas as pd

aspect = snakemake.wildcards.aspect
input_df = snakemake.input.annotation_df
bait_usage_df = snakemake.input.bait_usage
output_matrix = snakemake.output.feature_matrix

print(f"Processing GO {aspect} complete matrix...\n")

print(f"Loading {aspect} input dfs and filling missing values...")

df = pd.read_csv(input_df, sep = '\t')
df_copy = df.copy()
df_copy.fillna({
    'go_id': 'No_go_id',
    'annotation': 'No_annot',
    'depth': -1
}, inplace=True)

bait_usage = pd.read_csv(bait_usage_df, sep = '\t')

print(f"{aspect} Data loaded!\n")

print("Finding genes in common...")

all_genes = df_copy['entrez_id'].unique()
print(f"{aspect} Total gene universe: {len(all_genes)}")

bait_ids = bait_usage['entrez_id_bait'].unique()
df_copy = df_copy[df_copy['entrez_id'].isin(bait_ids)]
n_common_genes = df_copy['entrez_id'].nunique()

print(f"{n_common_genes} genes in common found for {aspect}!\n")

print(f"Unique GO {aspect} ids = {df_copy['go_id'].nunique() - 1}")

print(f"Creating {aspect} complete feature matrix...")

df_copy['value'] = 1

df_onehot = df_copy.pivot_table(
    index = 'entrez_id',
    columns = 'go_id',
    values = 'value',
    fill_value = 0
)

if 'No_go_id' in df_onehot.columns:
    df_onehot = df_onehot.drop(columns=['No_go_id'])

print(f"unique genes after feature matrix = {df_onehot.index.nunique()}")

print(f"GO {aspect} complete matrix ready!\n")

print(f"Saving {aspect} complete feature matrix...")

df_onehot.to_csv(output_matrix, sep = '\t')

print(f"{aspect} Complete feature matrix saved!\n")

print(f"GO {aspect} complete matrix processed!\n")