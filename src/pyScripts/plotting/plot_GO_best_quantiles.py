import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

input_file = snakemake.input.mean_adj_file
output_file = snakemake.output.plot_file
term_name = snakemake.wildcards.term

print(f"\n--- [Core] Generating best quantile plot for {term_name} ---", flush=True)

# Load the master mean adjacency matrix
df = pd.read_parquet(input_file)

# Isolate the arrays
true_values = df['PID0_mean_adj'].values
decoy_cols = [f"PID{i}_mean_adj" for i in range(1, 1001)]
decoy_matrix = df[decoy_cols].values

# Recalculate quantiles instantly to find the best prediction
quantiles = (decoy_matrix <= true_values[:, None]).sum(axis=1) / 1000.0
df['Quantile'] = quantiles

# Sort to find the absolute best row
df_sorted = df.sort_values(by=['Quantile', 'PID0_mean_adj'], ascending=[False, False])
best_row = df_sorted.iloc[0]

date = best_row['Date']
gene = best_row['Future_Gene']
best_true_val = best_row['PID0_mean_adj']
best_quantile = best_row['Quantile']

# Extract the 1,000 decoy values for this specific Date/Gene combination
null_distribution = best_row[decoy_cols].values.astype(float)

# PLOTTING
plt.figure(figsize=(10, 6))

# Plot the density of the null distribution (the 1000 permutations)
sns.kdeplot(null_distribution, fill=True, color="skyblue", alpha=0.6, linewidth=2, label="Null Distribution (1,000 Permutations)")

# Draw a bold vertical line for the True Annotation
plt.axvline(x=best_true_val, color="red", linestyle="--", linewidth=2.5, label=f"True Annotation (Adj: {best_true_val:.4f})")

# Formatting the plot
plt.title(f"Best quantile for {term_name}\nGene {gene} Date {date} (Quantile: {best_quantile})", fontsize=14, pad=15)
plt.xlabel("Mean Adjacency Value", fontsize=12)
plt.ylabel("Density", fontsize=12)

# Remove the top and right spines for a cleaner look
sns.despine()

plt.legend(loc='upper right', frameon=False)
plt.tight_layout()

# Save the plot securely to the Snakemake output directory
os.makedirs(os.path.dirname(output_file), exist_ok=True)
plt.savefig(output_file, dpi=300, bbox_inches = 'tight')
plt.close() # Free up memory 

print(f"--- [Core] Successfully saved best quantile plot for term {term_name} to: {output_file} ---\n", flush=True)