import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from scipy.stats import binned_statistic

adj_dir = Path(snakemake.input.mean_adj_dir)
dates_dir = Path(snakemake.input.annot_dates_dir)
output_plot = snakemake.output.plot_file
aspect = snakemake.wildcards.aspect

start_time = time.time()
print(f"--- [GLOBAL {aspect.upper()}] Plotting Individual Term Trajectories ---")

# Setup the figure for presentation
plt.figure(figsize=(14, 8))

# Lists to hold global data so we can still draw ONE bold average line on top
global_x = []
global_y = []

adj_files = list(adj_dir.glob("*_mean_adjacencies.parquet"))
print(f"Processing {len(adj_files)} terms...")

lines_plotted = 0

for adj_file in adj_files:
    term_str = adj_file.name.replace("_mean_adjacencies.parquet", "")
    date_file = dates_dir / f"{term_str}_first_annotation_dates.csv"
    
    if not date_file.exists():
        continue
        
    # Load and clean
    adj_df = pd.read_parquet(adj_file)[['Date', 'Future_Gene', 'PID0_mean_adj']]
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
        
    # Calculate Delta T in years
    df['Delta_T_years'] = (df['first_annotation_date'] - df['Date']).dt.days / 365.25
    
    # Save data for the global average later
    global_x.extend(df['Delta_T_years'])
    global_y.extend(df['PID0_mean_adj'])
    
    # CALCULATE THE TRENDLINE FOR THIS SPECIFIC TERM
    bins = np.arange(0, np.ceil(df['Delta_T_years'].max()) + 1, 1)
    bin_means, bin_edges, _ = binned_statistic(
        df['Delta_T_years'], df['PID0_mean_adj'], statistic='mean', bins=bins
    )
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    valid_bins = ~np.isnan(bin_means)
    
    # Only plot if this term has enough data to form a line (at least 2 points)
    if valid_bins.sum() > 1:
        # Layer 1: Plot this individual term (zorder=1, highly transparent)
        plt.plot(
            bin_centers[valid_bins], 
            bin_means[valid_bins], 
            color='#0571b0', 
            linewidth=1.5, 
            alpha=0.2,
            zorder=1
        )
        lines_plotted += 1

# Add a dummy line for the legend so the audience knows what the faint lines are
plt.plot([], [], color='#0571b0', linewidth=2.0, alpha=0.5, label='Individual GO Terms')

# Layer 2: Plot the Global Average on top to anchor the visual (zorder=2)
if len(global_x) > 0:
    global_bins = np.arange(0, np.ceil(max(global_x)) + 1, 1)
    global_means, global_edges, _ = binned_statistic(
        global_x, global_y, statistic='mean', bins=global_bins
    )
    global_centers = (global_edges[:-1] + global_edges[1:]) / 2
    valid_global = ~np.isnan(global_means)
    
    plt.plot(
        global_centers[valid_global], 
        global_means[valid_global], 
        color='#d73027', 
        linewidth=4.0, 
        label='Global Average Trend',
        zorder=2
    )

# Formatting
plt.title(f"Term-by-Term Network Predictive Power\nGO Aspect: {aspect.upper()}", fontsize=28, fontweight='bold', pad=20)
plt.xlabel("Time to Annotation ($\Delta$T in Years)", fontsize=22, labelpad=15)
plt.ylabel("Mean Adjacency Score", fontsize=22, labelpad=15)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Presentation Grid (zorder=0)
plt.grid(True, linestyle="--", alpha=0.4, zorder=0)

# Statistics Box
stats_text = f"Total Pathways Plotted: {lines_plotted:,}"
plt.text(
    0.95, 0.95, stats_text, 
    transform=plt.gca().transAxes,
    fontsize=20, verticalalignment='top', horizontalalignment='right',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#cccccc'),
    zorder=3
)

# Legend
plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(0.95, 0.82), framealpha=0.9)
plt.tight_layout()

# Save
print(f"Saving spaghetti plot to {output_plot}...")
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [GLOBAL {aspect.upper()}] Spaghetti plot generated in {elapsed_time} seconds! ---")