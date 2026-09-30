import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import time

input_stats = snakemake.input.global_stats
output_plot = snakemake.output.waterfall_plot
aspect = snakemake.wildcards.aspect

start_time = time.time()
print(f"--- [GLOBAL {aspect.upper()}] Plotting Spearman Bar Chart (p-value & labels) ---")

# Load Data
df = pd.read_csv(input_stats, sep='\t')

if df.empty or 'Spearman_rho' not in df.columns:
    print("No valid data available. Exiting.")
    open(output_plot, 'w').close()
    exit(0)

# Sort data
df = df.sort_values(by='Spearman_rho', ascending=True).reset_index(drop=True)
df['Rank'] = df.index + 1

# Define significance based on raw p-value
PVAL_THRESHOLD = 0.05
if 'p_value' in df.columns:
    df['Significant'] = df['p_value'] < PVAL_THRESHOLD
else:
    df['Significant'] = False

# Plotting Setup
plt.figure(figsize=(10, 7))
sns.set_context("talk")
sns.set_style("white")

colors = ['#2b8cbe' if sig else '#a6bddb' for sig in df['Significant']]

plt.grid(True, axis='y', linestyle='--', alpha=0.6, zorder=0)
plt.bar(df['Rank'], df['Spearman_rho'], color=colors, width=0.8, zorder=3)
plt.axhline(y=0, color='black', linewidth=1.5, linestyle='-', zorder=4)

for i, val in enumerate(df['Spearman_rho']):
    y_offset = -0.05 if val < 0 else 0.05
    va_align = 'top' if val < 0 else 'bottom'
    
    plt.text(
        df['Rank'][i], 
        val + y_offset, 
        f"{val:.5f}", 
        ha='center', 
        va=va_align, 
        fontsize=11, 
        fontweight='bold',
        color='#333333',
        zorder=5
    )

# Custom Legend
legend_elements = [
    Patch(facecolor='#2b8cbe', label=f'Significant ($p < {PVAL_THRESHOLD}$)'),
    Patch(facecolor='#a6bddb', label='Not Significant')
]
plt.legend(handles=legend_elements, loc='best', fontsize=12, framealpha=0.9)

# Formatting
plt.title(f"Correlation Between Mean Adjacency and TTA\nGO Aspect: {aspect.upper()}", fontsize=18, fontweight='bold', pad=15)
plt.xlabel("GO Terms", fontsize=16, fontweight='bold', labelpad=15)
plt.ylabel("Spearman Correlation ($\\rho$)", fontsize=16, fontweight='bold', labelpad=10)

plt.yticks(np.arange(-1.0, 1.2, 0.2))

if 'GO_id' in df.columns:
    plt.xticks(df['Rank'], df['GO_id'], ha='center', fontsize=10)
else:
    plt.xticks(df['Rank'], [f"Term {i}" for i in df['Rank']], rotation=90, ha='center')

plt.xlim(0.4, len(df) + 0.6)
plt.ylim(-1.15, 1.15) 

sns.despine(bottom=True)
plt.tight_layout()

# Save
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- Plot successfully generated in {elapsed}s! ---")