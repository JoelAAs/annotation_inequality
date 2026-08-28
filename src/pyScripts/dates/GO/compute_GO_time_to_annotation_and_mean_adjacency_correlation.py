import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from pathlib import Path
import time

adj_dir = Path(snakemake.input.mean_adj_dir)
dates_dir = Path(snakemake.input.annot_dates_dir)
output_file = snakemake.output.global_stats

aspect = snakemake.wildcards.aspect

start_time = time.time()
print(f"--- [GLOBAL {aspect.upper()}] Calculating Spearman Correlations for All Terms ---")

results = []
adj_files = list(adj_dir.glob("*_mean_adjacencies.parquet"))
print(f"Found {len(adj_files)} adjacency files. Processing one by one...")

for adj_file in adj_files:
    term_str = adj_file.name.replace("_mean_adjacencies.parquet", "")
    go_id = term_str.replace("_", ":")
    
    # Locate matching date file
    date_file = dates_dir / f"{term_str}_first_annotation_dates.csv"
    if not date_file.exists():
        continue
        
    # LOAD DATA FOR THIS SPECIFIC TERM
    adj_df = pd.read_parquet(adj_file)[['Date', 'Future_Gene', 'PID0_mean_adj']]
    dates_df = pd.read_csv(date_file, sep='\t')[['gene_id', 'first_annotation_date']]
    
    # MERGE
    adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
    dates_df['gene_id'] = dates_df['gene_id'].astype(str)
    
    df = pd.merge(
        adj_df, 
        dates_df, 
        left_on='Future_Gene', 
        right_on='gene_id', 
        how='inner'
    )
    
    # CALCULATE DELTA T
    df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
    df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')
    df = df[df['first_annotation_date'] > df['Date']].copy()
    
    n_preds = len(df)
    
    # COMPUTE CORRELATION
    if n_preds >= 12:
        df['Delta_T'] = (df['first_annotation_date'] - df['Date']).dt.days
        rho, pval = spearmanr(df['Delta_T'], df['PID0_mean_adj'])
    else:
        rho, pval = np.nan, np.nan
        
    # STORE RESULT
    results.append({
        'GO_id': go_id,
        'N_predictions': n_preds,
        'Spearman_rho': rho,
        'p_value': pval
    })

# Compile master table
print(f"Processing complete. Compiling master table...")
results_df = pd.DataFrame(results)

# Drop terms that failed the N>=12 check, and sort from most negative correlation to least
results_df = results_df.dropna(subset=['Spearman_rho']).sort_values('GO_id')

# Save as CSV
print(f"Saving master table to {output_file}...")
results_df.to_csv(output_file, sep='\t', index=False)

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [GLOBAL {aspect.upper()}] Master table successfully generated in {elapsed_time} seconds! ---")