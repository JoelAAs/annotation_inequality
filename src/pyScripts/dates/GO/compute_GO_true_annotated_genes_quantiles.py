import os
import pandas as pd
import numpy as np

input_file = snakemake.input.mean_adj_file
output_file = snakemake.output.quantile_file
term = snakemake.wildcards.term

print(f"\n--- [Core] Computing Quantiles for {term} ---\n", flush=True)

# Load the master mean adjacency matrix
df = pd.read_parquet(input_file)

# Isolate the true values (1D array)
true_values = df['PID0_mean_adj'].values

# Isolate the 1,000 decoy permutations (2D matrix)
decoy_cols = [f"PID{i}_mean_adj" for i in range(1, 1001)]
decoy_matrix = df[decoy_cols].values

# QUANTILE COMPUTATION
# For every row, how many decoys are < the true value? (We could also do <= but there will be a lot of 1s)
# Should we filter out where there are many 0s?
quantiles = (decoy_matrix < true_values[:, None]).sum(axis=1) / 1000.0

# Final dataframe
results_df = pd.DataFrame({
    'Date': df['Date'],
    'Future_Gene': df['Future_Gene'],
    'True_Mean_Adj': true_values.astype(np.float32), 
    'Quantile': quantiles.astype(np.float32)
})

# Save the resulting quantiles
os.makedirs(os.path.dirname(output_file), exist_ok=True)
results_df.to_parquet(output_file, engine='pyarrow', index=False)

print(f"--- [Core] Successfully saved Quantile Matrix for term {term} to: {output_file} ---\n", flush=True)