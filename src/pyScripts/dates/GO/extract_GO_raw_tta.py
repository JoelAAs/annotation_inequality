import os
import pandas as pd
import time

input_adj_file = snakemake.input.mean_adj_file
input_dates_file = snakemake.input.annot_dates_file
output_file = snakemake.output.raw_tta_file
term = snakemake.wildcards.term
term_formatted = term.replace('_', ':')

start_time = time.time()
print(f"--- Extracting Raw TTA for {term} ---", flush=True)

# Load the data (we only need Date and Future_Gene from the adj file)
adj_df = pd.read_parquet(input_adj_file)[['Date', 'Future_Gene']]
dates_df = pd.read_csv(input_dates_file, sep='\t')
dates_df = dates_df[dates_df['GO_id'] == term_formatted]

# Format and Merge
adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
dates_df['gene_id'] = dates_df['gene_id'].astype(str)

df = pd.merge(
    adj_df, 
    dates_df[['gene_id', 'first_annotation_date']], 
    left_on='Future_Gene', 
    right_on='gene_id', 
    how='inner'
)

# Time Math
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')

# Filter for strictly future annotations
df = df[df['first_annotation_date'] > df['Date']].copy()

if df.empty:
    print(f"[{term}] No valid data. Creating empty output.")
    pd.DataFrame(columns=['Date', 'Future_Gene', 'TTA_Years']).to_parquet(output_file, index=False)
    exit(0)

# Calculate exact TTA in years
df['TTA_Years'] = (df['first_annotation_date'] - df['Date']).dt.days / 365.25

# Save
os.makedirs(os.path.dirname(output_file), exist_ok=True)
df[['Date', 'Future_Gene', 'TTA_Years']].to_parquet(output_file, engine='pyarrow', index=False)

elapsed = round(time.time() - start_time, 2)
print(f"--- Saved {len(df)} predictions for {term} in {elapsed}s ---", flush=True)