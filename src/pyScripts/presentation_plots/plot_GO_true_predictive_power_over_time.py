import pandas as pd
import matplotlib.pyplot as plt
import time

input_file = snakemake.input.mean_adj_file
output_file = snakemake.output.plot_file
term = snakemake.wildcards.term

start_time = time.time()
print(f"--- [{term}] Starting True Annotation plotting ---")
print(f"[{term}] [LOAD] Reading data from {input_file}...")

df = pd.read_parquet(input_file)
df.columns = df.columns.str.strip() # Clean column names just in case!

# Date formatting
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')

print(f"[{term}] [PROCESS] Grouping {len(df)} rows by Date...")

# Drop 'Future_Gene' as we just want the averages per date
df_daily_mean = df.drop(columns=['Future_Gene']).groupby('Date').mean()

# Extract only the True Annotations line
true_line = df_daily_mean['PID0_mean_adj']

# Plotting
print(f"[{term}] [PLOT] Drawing the visualization...")
plt.figure(figsize=(12, 6))

# Plot ONLY the True Annotations (PID0) - Thicker line for presentation
plt.plot(true_line.index, true_line, label='True Annotations', color='crimson', linewidth=2.5)

# Formatting - Updated title to reflect the single line
plt.title(f"Predictive Power Over Time: True Annotations\nTerm: {term}", fontsize=16, pad=15)
plt.xlabel("Date", fontsize=14)
plt.ylabel("Mean Adjacency Score", fontsize=14)
plt.legend(loc="upper left", fontsize=12)
plt.xticks(rotation=45, fontsize=11)
plt.yticks(fontsize=11)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

# Save
print(f"[{term}] [SAVE] Saving image to disk...")
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [{term}] Plot successfully generated in {elapsed_time} seconds! ---")