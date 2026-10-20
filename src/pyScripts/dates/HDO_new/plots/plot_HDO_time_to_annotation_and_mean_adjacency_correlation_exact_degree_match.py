import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import time

input_stats = snakemake.input.global_stats
output_plot = snakemake.output.waterfall_plot

start_time = time.time()
print("--- [GLOBAL] Plotting HDO Spearman Bar Chart (p-value & labels) ---")

os.makedirs(os.path.dirname(output_plot), exist_ok=True)

# Load Data
df = pd.read_csv(input_stats, sep='\t')

if df.empty or 'Spearman_rho' not in df.columns:
    print("--- [WARNING] No valid correlation data available. Creating placeholder plot. ---")
    plt.figure(figsize=(10, 7))
    plt.text(0.5, 0.5, "No valid correlation data available\nfor HDO at this depth/cutoff", 
             horizontalalignment='center', verticalalignment='center', 
             fontsize=16, transform=plt.gca().transAxes)
    plt.title("HDO Spearman Correlation (No Data)", fontsize=18, fontweight='bold', pad=15)
    plt.axis('off')
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    exit(0)

# Sort data (just to be safe, NaNs at the end)
df = df.sort_values(by='Spearman_rho', na_position='last').reset_index(drop=True)
df['Rank'] = df.index + 1

# Define significance based on raw p-value, treating NaNs as False
PVAL_THRESHOLD = 0.05
df['Significant'] = (df['p_value'] < PVAL_THRESHOLD) & (df['p_value'].notna())

# Identify missing data for formatting
is_missing = df['Spearman_rho'].isna()
plot_rho = df['Spearman_rho'].fillna(0) # Temporarily fill with 0 just so plt.bar doesn't crash

# Plotting Setup
plt.figure(figsize=(10, 7))
sns.set_context("talk")
sns.set_style("white")

# Colors: blue if sig, light blue if not sig, light grey if missing
colors = ['#2b8cbe' if sig else '#e0e0e0' if miss else '#a6bddb' for sig, miss in zip(df['Significant'], is_missing)]

plt.grid(True, axis='y', linestyle='--', alpha=0.6, zorder=0)
plt.bar(df['Rank'], plot_rho, color=colors, width=0.8, zorder=3)
plt.axhline(y=0, color='black', linewidth=1.5, linestyle='-', zorder=4)

for i, val in enumerate(df['Spearman_rho']):
    if pd.isna(val):
        # Missing data: write N/A above the zero line
        plt.text(df['Rank'][i], 0.05, "N/A", ha='center', va='bottom', fontsize=11, fontweight='bold', color='#999999', zorder=5)
    else:
        # Valid data: write the actual correlation value
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
    Patch(facecolor='#a6bddb', label='Not Significant'),
    Patch(facecolor='#e0e0e0', label='N/A (Insufficient Data)')
]
plt.legend(handles=legend_elements, loc='best', fontsize=12, framealpha=0.9)

# Formatting
plt.title("Correlation Between Mean Adjacency and TTA\nHuman Disease Ontology (HDO)", fontsize=18, fontweight='bold', pad=15)
plt.xlabel("HDO Terms", fontsize=16, fontweight='bold', labelpad=15)
plt.ylabel("Spearman Correlation ($\\rho$)", fontsize=16, fontweight='bold', labelpad=10)

plt.yticks(np.arange(-1.0, 1.2, 0.2))

if 'HDO_id' in df.columns:
    rotation_angle = 45 if len(df) > 5 else 0
    align = 'right' if len(df) > 5 else 'center'
    label_size = 10 if len(df) <= 15 else 8
    plt.xticks(df['Rank'], df['HDO_id'], rotation=rotation_angle, ha=align, fontsize=label_size)
else:
    plt.xticks(df['Rank'], [f"Term {i}" for i in df['Rank']], rotation=90, ha='center')

plt.xlim(0.4, len(df) + 0.6)
plt.ylim(-1.15, 1.15) 

sns.despine(bottom=True)
plt.tight_layout()

# Save
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- [GLOBAL] HDO Plot successfully generated in {elapsed}s! ---")