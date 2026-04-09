import pandas as pd

aspect = snakemake.wildcards.aspect
annot_df = snakemake.input.annotation_df
input_bait_usage = snakemake.input.bait_usage
output_counts = snakemake.output.counts_per_annot

print(f"Processing GO {aspect} annotations counts...")

print(f"Loading {aspect} data...")

df = pd.read_csv(annot_df, sep = '\t')
df_copy = df.copy()
df_copy = df_copy.dropna()

bait_usage = pd.read_csv(input_bait_usage, sep = '\t')

print(f"{aspect} Data loaded!\n")

print("Merging to obtain common genes...")

all_genes = df_copy['entrez_id'].unique()
print(f"Total gene universe: {len(all_genes)}")

bait_ids = bait_usage['entrez_id_bait'].unique()
df_copy = df_copy[df_copy['entrez_id'].isin(bait_ids)]
n_common_genes = df_copy['entrez_id'].nunique()

print(f"{n_common_genes} common genes obtained!\n")

counts_df = (
    df_copy.drop_duplicates(subset=['entrez_id', 'go_id'])
    .groupby('go_id')
    .size()
    .reset_index(name='count')
)
counts_df.columns = ['go_id', 'count']

print(f"Saving {aspect} annotations counts...")

counts_df.to_csv(output_counts, sep = '\t', index = False)

print(f"{aspect} annotations counts saved!\n")

print(f"GO {aspect} annotations counts obtained!")