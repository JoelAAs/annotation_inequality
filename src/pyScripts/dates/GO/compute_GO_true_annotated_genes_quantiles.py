import os
import pandas as pd
import numpy as np

input_file = snakemake.input.mean_adj_file
output_file = snakemake.output.quantile_file
term = snakemake.wildcards.term

print(f"\n--- [Core] Computing Quantiles for {term} ---\n", flush=True)

# Load the master mean adjacency matrix
df = pd.read_parquet(input_file)

# Remove the "Future_Gene" column
df = df.drop('Future_Gene', axis=1)

# Group by date to get the global temporal signal
df = df.groupby('Date', as_index=False).mean()

# Isolate columns dynamically
true_col = 'PID0_mean_adj'
decoy_cols = [c for c in df.columns if c.startswith('PID') and c != true_col]
n_decoys = len(decoy_cols)

# QUANTILE COMPUTATION (Vectorized with Pseudo-Count)
# We add +1 to both numerator and denominator per to include the "true" observation (Laplace correction)
df['quantile'] = (df[decoy_cols].lt(df[true_col], axis=0).sum(axis=1) + 1) / (n_decoys + 1)

# Final dataframe
results_df = pd.DataFrame({
    'Date': df['Date'],
    'True_Mean_Adj': df[true_col].astype(np.float32), 
    'Quantile': df['quantile'].astype(np.float32)
})

# Save the resulting quantiles
os.makedirs(os.path.dirname(output_file), exist_ok=True)
results_df.to_parquet(output_file, engine='pyarrow', index=False)

print(f"--- [Core] Successfully saved Quantile Matrix for term {term} to: {output_file} ---\n", flush=True)