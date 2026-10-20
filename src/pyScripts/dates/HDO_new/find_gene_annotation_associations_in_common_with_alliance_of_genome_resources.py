import os
import pandas as pd

file_agr = snakemake.input.agr_file
file_target = snakemake.input.annot_file
file_hgnc = snakemake.input.hgnc_file
file_output = snakemake.output.gene_annot_associations

# Load AGR dataset
print(f"--- [INFO] Reading AGR dataset from: {file_agr} ---")
skip_idx = 0
with open(file_agr, "r") as f:
    for i, line in enumerate(f):
        if line.startswith("UniqueID"):
            skip_idx = i
            break

# Read AGR, dropping commented headers correctly
df_agr = pd.read_csv(file_agr, sep="\t", skiprows=skip_idx, dtype=str)

# HGNC to ENTREZ mapping
print("--- [INFO] Fetching HGNC to Entrez ID mapping... ---")
df_mapping = pd.read_csv(file_hgnc, sep="\t", usecols=['hgnc_id', 'entrez_id'], dtype=str)
df_mapping = df_mapping.dropna(subset=['entrez_id'])

# Merge Entrez IDs into AGR
df_agr = pd.merge(df_agr, df_mapping, left_on='Gene IDs', right_on='hgnc_id', how='left')

# Load target dataset
print(f"--- [INFO] Reading Target dataset from: {file_target} ---")
df_target = pd.read_csv(file_target, sep="\t", dtype=str)

# Filter out 'No_doid'
df_target_valid = df_target[df_target['doid'] != 'No_doid'].copy()
initial_target_count = len(df_target_valid.drop_duplicates(subset=['entrez_id', 'doid']))
print(f"--- [INFO] Found {initial_target_count} unique valid target associations ---")

# Merge target with AGR
print("--- [INFO] Merging Target and AGR datasets... ---")
# Left join ensures we keep all valid target associations even if they are not in AGR
df_merge = pd.merge(
    df_target_valid,
    df_agr[['entrez_id', 'Disease ID', 'Reference']],
    left_on=['entrez_id', 'doid'],
    right_on=['entrez_id', 'Disease ID'],
    how='left'
)

# Explode references in case of multiple reference for one single gene-annotation association
print("--- [INFO] Processing and exploding references... ---")
# Replace pandas NaN or "nan" string with 'no_reference'
df_merge['Reference'] = df_merge['Reference'].fillna('no_reference').replace('nan', 'no_reference')

# Split multiple references separated by '|' into lists, then explode into separate rows
df_merge['Reference'] = df_merge['Reference'].str.split('|')
df_exploded = df_merge.explode('Reference')

# Clean up any potential whitespace
df_exploded['Reference'] = df_exploded['Reference'].str.strip()

# Drop exact duplicates created by the merge/explode process
df_final = df_exploded.drop_duplicates(subset=['entrez_id', 'doid', 'depth', 'term', 'Reference'])

# Clean up
print("--- [INFO] Cleaning up redundant 'no_reference' entries... ---")
# If a gene-disease pair has valid references AND a 'no_reference' row, drop the 'no_reference' row.
valid_refs = df_final[df_final['Reference'] != 'no_reference'][['entrez_id', 'doid']].drop_duplicates()
valid_refs['has_valid_ref'] = True

df_final = pd.merge(df_final, valid_refs, on=['entrez_id', 'doid'], how='left')

# Keep the row if it's a valid reference OR if it truly has no valid references at all
mask = (df_final['Reference'] != 'no_reference') | (df_final['has_valid_ref'].isna())
df_final = df_final[mask].drop(columns=['has_valid_ref', 'Disease ID'])

# Ensure all columns are strings (avoids Arrow/Parquet mixed-type schema errors)
df_final = df_final.astype(str)

# Save
total_rows = len(df_final)
pmid_rows = len(df_final[df_final['Reference'].str.startswith('PMID', na=False)])
no_ref_rows = len(df_final[df_final['Reference'] == 'no_reference'])

print(f"--- [INFO] Final dataset contains {total_rows} total rows ---")
print(f"--- [INFO] Rows with a PMID: {pmid_rows} ---")
print(f"--- [INFO] Rows with no_reference: {no_ref_rows} ---")

print("--- [INFO] Breakdown of all Reference types: ---")
ref_prefixes = df_final['Reference'].str.split(':').str[0]
ref_counts = ref_prefixes.value_counts()

for ref_type, count in ref_counts.items():
    print(f"      > {ref_type}: {count}")

# Save to PARQUET output file
print(f"--- [INFO] Saving exploded associations to Parquet format... ---")
df_final.to_parquet(file_output, index=False)
print(f"--- [INFO] Process completed successfully ---")