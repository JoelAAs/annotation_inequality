import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import time

adj_dir = Path(snakemake.input.mean_adj_dir)
dates_dir = Path(snakemake.input.annot_dates_dir)
output_file = snakemake.output.global_stats

start_time = time.time()
print("--- [GLOBAL] Calculating HDO Spearman Correlations for All Terms ---")

results = []
adj_files = list(adj_dir.glob("*_mean_adjacencies.parquet"))
print(f"Found {len(adj_files)} adjacency files. Processing one by one...")

for adj_file in adj_files:
    term_str = adj_file.name.replace("_mean_adjacencies.parquet", "")
    hdo_id = term_str.replace("_", ":")
    
    # Helper to record empty/missing terms with NaN before skipping
    def record_missing(reason):
        print(f"  [{term_str}] Skipping: {reason}. Recording as N/A in table.")
        results.append({
            'HDO_id': hdo_id,
            'N_predictions': 0,
            'Spearman_rho': np.nan,
            'p_value': np.nan
        })

    # Locate matching date file
    date_file = dates_dir / f"{term_str}_first_annotation_dates.csv"
    if not date_file.exists():
        record_missing("No date file found")
        continue
        
    # LOAD DATA FOR THIS SPECIFIC TERM (with defensive checks)
    try:
        adj_df = pd.read_parquet(adj_file)
        if adj_df.empty or 'PID0_mean_adj' not in adj_df.columns:
            record_missing("Matrix is empty")
            continue
        adj_df = adj_df[['Date', 'Future_Gene', 'PID0_mean_adj']]
        
        dates_df = pd.read_csv(date_file, sep='\t')
        if dates_df.empty:
            record_missing("Dates file is empty")
            continue
        dates_df = dates_df[['gene_id', 'first_annotation_date']]
    except Exception as e:
        record_missing(f"Error reading files: {e}")
        continue
        
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
    
    if df.empty:
        record_missing("No matching genes after merge")
        continue

    # CALCULATE DELTA T SAFELY
    df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d', errors='coerce')
    df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d', errors='coerce')
    df = df.dropna(subset=['Date', 'first_annotation_date'])
    
    df = df[df['first_annotation_date'] > df['Date']].copy()
    
    n_preds = len(df)
    
    # COMPUTE CORRELATION
    if n_preds >= 12:
        df['Delta_T'] = (df['first_annotation_date'] - df['Date']).dt.days
        rho, pval = spearmanr(df['Delta_T'], df['PID0_mean_adj'])
    else:
        rho, pval = np.nan, np.nan
        
    # STORE RESULT (For valid terms)
    results.append({
        'HDO_id': hdo_id,
        'N_predictions': n_preds,
        'Spearman_rho': rho,
        'p_value': pval
    })

# Compile master table
print("\nProcessing complete. Compiling master table...")

# Ensure output directory exists
output_path = Path(output_file)
output_path.parent.mkdir(parents=True, exist_ok=True)

if not results:
    print("--- [WARNING] No valid results for any terms. Creating placeholder CSV. ---")
    empty_df = pd.DataFrame(columns=['HDO_id', 'N_predictions', 'Spearman_rho', 'p_value', 'FDR'])
    empty_df.to_csv(output_file, sep='\t', index=False)
    exit(0)

results_df = pd.DataFrame(results)

# INSTEAD OF DROPPING NAs, KEEP THEM FOR THE PLOT
valid_mask = results_df['p_value'].notna()
results_df['FDR'] = np.nan

# Calculate FDR (Benjamini-Hochberg) ONLY on valid rows
if valid_mask.any():
    _, fdr_values, _, _ = multipletests(results_df.loc[valid_mask, 'p_value'], method='fdr_bh')
    results_df.loc[valid_mask, 'FDR'] = fdr_values

# Sort: valid Spearman_rho values first, N/A terms pushed to the right (last)
results_df = results_df.sort_values(by='Spearman_rho', na_position='last')

# Save as CSV
print(f"Saving master table to {output_file}...")
results_df.to_csv(output_file, sep='\t', index=False)

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [GLOBAL] HDO Master table successfully generated in {elapsed_time} seconds! ---")