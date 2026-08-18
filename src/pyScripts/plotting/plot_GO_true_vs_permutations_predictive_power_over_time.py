import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

input_file = snakemake.input.mean_adj_file
output_file = snakemake.output.plot_file
term = snakemake.wildcards.term

start_time = time.time()
print(f"--- [{term}] Starting True vs Permutation plotting ---")
print(f"[{term}] [LOAD] Reading data from {input_file}...")

df = pd.read_parquet(input_file)
df.columns = df.columns.str.strip() # Clean column names just in case!

# Date formating
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')

print(f"[{term}] [PROCESS] Grouping {len(df)} rows by Date to prevent RAM overload...")

# Drop 'Future_Gene' as we just want the averages per date
df_daily_mean = df.drop(columns=['Future_Gene']).groupby('Date').mean()

# Separate True vs Permutations
true_line = df_daily_mean['PID0_mean_adj']

# Grab all the other 1000 PID columns
perm_cols = [col for col in df_daily_mean.columns if col.startswith('PID') and col != 'PID0_mean_adj']
perm_data = df_daily_mean[perm_cols]

print(f"[{term}] [MATH] Calculating 95% Confidence Intervals across {len(perm_cols)} permutations...")

# Calculate the mean and the 95% bounds of the random permutations for the shaded area
perm_mean = perm_data.mean(axis=1)
perm_lower = perm_data.quantile(0.025, axis=1) # 2.5% percentile
perm_upper = perm_data.quantile(0.975, axis=1) # 97.5% percentile

# Plotting
print(f"[{term}] [PLOT] Drawing the visualization...")
plt.figure(figsize=(12, 6))

# Plot the 1000 Permutations (Mean + Shaded CI)
plt.plot(perm_mean.index, perm_mean, label='Random Permutations (Mean)', color='gray', linestyle='--', linewidth=2)
plt.fill_between(
    perm_mean.index, 
    perm_lower, 
    perm_upper, 
    color='gray', 
    alpha=0.3, 
    label='Permutations 95% CI'
)

# Plot the True Annotations (PID0) on top
plt.plot(true_line.index, true_line, label='True Annotations (PID0)', color='crimson', linewidth=2.5, marker='o', markersize=4)

# Formatting
plt.title(f"Predictive Power Over Time: True Annotations vs. Random Permutations\nTerm: {term}", fontsize=14, pad=15)
plt.xlabel("Date", fontsize=12)
plt.ylabel("Mean Adjacency Score", fontsize=12)
plt.legend(loc="upper left")
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# 6. Save
print(f"[{term}] [SAVE] Saving image to disk...")
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [{term}] Plot successfully generated in {elapsed_time} seconds! ---")