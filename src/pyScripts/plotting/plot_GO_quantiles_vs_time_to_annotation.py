import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from scipy.stats import binned_statistic
import itertools

# 1. Map directly to Snakemake rule inputs/outputs
adj_dir = Path(snakemake.input.mean_adj_dir)
dates_dir = Path(snakemake.input.annot_dates_dir)
output_plot = Path(snakemake.output.plot_quantiles)

# 2. Extract wildcards for logging and the plot title
aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff

start_time = time.time()
print(f"--- [GLOBAL {aspect.upper()} | Depth:{depth} | Cutoff:{cutoff}] Plotting Quantile Trajectories ---")

plt.figure(figsize=(14, 8))

global_x = []
global_y = []

adj_files = list(adj_dir.glob("*_mean_adjacencies.parquet"))
print(f"Processing {len(adj_files)} terms...")

tab20 = plt.get_cmap('tab20').colors
safe_colors = [tab20[i] for i in [0, 1, 4, 5, 8, 9, 14, 15, 18, 19]]
color_cycle = itertools.cycle(safe_colors)

# Generate the list of permutation column names just once
perm_cols = [f'PID{i}_mean_adj' for i in range(1, 1001)]
columns_to_load = ['Date', 'Future_Gene', 'PID0_mean_adj'] + perm_cols

for adj_file in adj_files:
    term_str = adj_file.name.replace("_mean_adjacencies.parquet", "")
    date_file = dates_dir / f"{term_str}_first_annotation_dates.csv"
    
    if not date_file.exists():
        continue
        
    # Load the true value AND all 1000 permutations
    adj_df = pd.read_parquet(adj_file)[columns_to_load]
    dates_df = pd.read_csv(date_file, sep='\t')[['gene_id', 'first_annotation_date']]
    
    adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
    dates_df['gene_id'] = dates_df['gene_id'].astype(str)
    
    # Merge
    df = pd.merge(adj_df, dates_df, left_on='Future_Gene', right_on='gene_id', how='inner')
    
    df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
    df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')
    df = df[df['first_annotation_date'] > df['Date']].copy()
    
    if df.empty:
        continue
        
    df['Delta_T_years'] = (df['first_annotation_date'] - df['Date']).dt.days / 365.25
    
    # --- CALCULATE THE QUANTILE ---
    true_vals = df['PID0_mean_adj'].values[:, None]
    perm_vals = df[perm_cols].values
    
    # Strictly less than avoids scoring ties at 0 as artificially high
    df['Quantile'] = (perm_vals < true_vals).sum(axis=1) + 1 / 1000.0 + 1
    
    global_x.extend(df['Delta_T_years'])
    global_y.extend(df['Quantile'])
    
    # CALCULATE THE TRENDLINE FOR THIS SPECIFIC TERM (0.5 = 6 months)
    bins = np.arange(0, np.ceil(df['Delta_T_years'].max()) + 1, 0.5)
    
    bin_means, bin_edges, _ = binned_statistic(
        df['Delta_T_years'], df['Quantile'], statistic='mean', bins=bins
    )
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    valid_bins = ~np.isnan(bin_means)
    
    if valid_bins.sum() > 1:
        term_color = next(color_cycle)
        plt.plot(
            bin_centers[valid_bins], 
            bin_means[valid_bins], 
            color=term_color, 
            linewidth=1.5, 
            alpha=0.4,  
            zorder=1,
            label=term_str.replace('_', ':')
        )

# Layer 2: Plot the Global Average on top
if len(global_x) > 0:
    global_bins = np.arange(0, np.ceil(max(global_x)) + 1, 0.5)
    
    global_means, global_edges, _ = binned_statistic(
        global_x, global_y, statistic='mean', bins=global_bins
    )
    global_centers = (global_edges[:-1] + global_edges[1:]) / 2
    valid_global = ~np.isnan(global_means)
    
    plt.plot(
        global_centers[valid_global], 
        global_means[valid_global], 
        color='#d73027', 
        linewidth=5.0, 
        label='Global Average Trend',
        zorder=2
    )

# Formatting: added depth and cutoff to title
title_str = (
    f"Mean Quantile Score vs Time to Annotation\n"
    f"GO Aspect: {aspect.upper()} | Depth: {depth} | Cutoff: {cutoff}"
)
plt.title(title_str, fontsize=24, fontweight='bold', pad=20)
plt.xlabel("Time to Annotation ($\\Delta$T in Years)", fontsize=22, labelpad=15)
plt.ylabel("Mean Quantile (vs 1000 permutations)", fontsize=22, labelpad=15)

# Lock the Y-axis between 0 and 1 since it's a quantile
plt.ylim(-0.05, 1.05) 

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.grid(True, linestyle="--", alpha=0.4, zorder=0)

num_terms = len(adj_files)

if num_terms <= 20:
    cols = 1
elif num_terms <= 40:
    cols = 2
else:
    cols = 1
    handles, labels = plt.gca().get_legend_handles_labels()
    handles = [h for h, l in zip(handles, labels) if l == 'Global Average Trend']
    labels = ['Global Average Trend']
    plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=12)

if num_terms <= 40:
    plt.legend(
        loc="upper left", 
        bbox_to_anchor=(1.02, 1), 
        fontsize=10, 
        ncol=cols, 
        framealpha=0.9, 
        edgecolor='#cccccc'
    )

print(f"Saving quantile spaghetti plot to {output_plot}...")
output_plot.parent.mkdir(parents=True, exist_ok=True) # Ensure output directory exists just in case
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [GLOBAL {aspect.upper()}] Spaghetti plot generated in {elapsed_time} seconds! ---")