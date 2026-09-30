import pandas as pd

hdo_annotations_path = snakemake.input.hdo_annotations
go_annotations_path = snakemake.input.go_annotations
bp_df_path = snakemake.input.bp_df
feature_matrix_output_path = snakemake.output.feature_matrix

print(f"--- [START] MIXED FEATURE MATRIX COMPUTATION STARTED ---", flush=True)

# Load data
print(f"--- [RUNNING] LOADING INPUT DATA ---", flush=True)
bp_df = pd.read_parquet(bp_df_path)
hdo_df = pd.read_csv(hdo_annotations_path, sep='\t')
go_df = pd.read_csv(go_annotations_path, sep='\t')
print(f"--- [INFO] INPUT DATA LOADED SUCCESSFULLY ---", flush=True)

# Extract unique genes universe
print(f"--- [RUNNING] EXTRACTING UNIQUE GENES FROM BAIT AND PREY ---", flush=True)
all_genes = pd.concat([bp_df['entrez_id_bait'], bp_df['entrez_id_prey']], ignore_index=True)
# Remove NaNs and take unique ids
universe_genes = all_genes.dropna().astype(str).unique()
print(f"--- [INFO] TOTAL UNIQUE GENES IN UNIVERSE: {len(universe_genes)} ---", flush=True)

# Process HDO annotations
print(f"--- [RUNNING] PROCESSING HDO ANNOTATIONS ---", flush=True)
hdo_df['entrez_id'] = hdo_df['entrez_id'].astype(str)
hdo_df = hdo_df[hdo_df['doid'] != 'No_doid']
hdo_df = hdo_df[['entrez_id', 'doid']].rename(columns={'doid': 'annotation_id'})

# Process GO annotations
print(f"--- [RUNNING] PROCESSING GO ANNOTATIONS ---", flush=True)
go_df['entrez_id'] = go_df['entrez_id'].astype(str)
go_df = go_df[['entrez_id', 'go_id']].rename(columns={'go_id': 'annotation_id'})

# Merge and clean annotations
print(f"--- [RUNNING] MERGING HDO AND GO ANNOTATIONS ---", flush=True)
combined_annotations = pd.concat([hdo_df, go_df], ignore_index=True)
combined_annotations.dropna(subset=['annotation_id'], inplace=True)
combined_annotations.drop_duplicates(inplace=True)

# Filter annotations to keep only genes present in the universe
combined_annotations = combined_annotations[combined_annotations['entrez_id'].isin(universe_genes)]
print(f"--- [INFO] ANNOTATIONS FILTERED TO UNIVERSE: {combined_annotations['annotation_id'].nunique()} TERMS FOR {combined_annotations['entrez_id'].nunique()} GENES ---", flush=True)

# Build binary feature matrix
print(f"--- [RUNNING] BUILDING BINARY FEATURE MATRIX ---", flush=True)
# Assign 1 to all valid gene-annotation associations
combined_annotations['value'] = 1
feature_matrix = combined_annotations.pivot_table(
    index='entrez_id', 
    columns='annotation_id', 
    values='value', 
    fill_value=0
)

# Add missing genes (genes in universe but with absolutely no valid annotations)
print(f"--- [RUNNING] ADDING MISSING GENES (ZERO ANNOTATIONS) ---", flush=True)
missing_genes = list(set(universe_genes) - set(feature_matrix.index))

if missing_genes:
    # Zero DF for missing genes
    missing_df = pd.DataFrame(0, index=missing_genes, columns=feature_matrix.columns)
    # Add it to feature matrixx
    feature_matrix = pd.concat([feature_matrix, missing_df])

# Sort genes
feature_matrix.sort_index(inplace=True)
print(f"--- [INFO] FINAL FEATURE MATRIX SHAPE: {feature_matrix.shape[0]} GENES, {feature_matrix.shape[1]} ANNOTATIONS ---", flush=True)

# Save to Parquet
print(f"--- [RUNNING] SAVING FEATURE MATRIX TO PARQUET ---", flush=True)
feature_matrix.columns = feature_matrix.columns.astype(str)
feature_matrix.to_parquet(feature_matrix_output_path)

print(f"--- [COMPLETE] MIXED FEATURE MATRIX COMPUTATION FINISHED SUCCESSFULLY ---", flush=True)