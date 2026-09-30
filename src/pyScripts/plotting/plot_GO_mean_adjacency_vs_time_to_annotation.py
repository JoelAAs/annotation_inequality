import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import time
from scipy.stats import binned_statistic
import itertools

adj_dir = Path(snakemake.input.mean_adj_dir)
dates_dir = Path(snakemake.input.annot_dates_dir)
output_plot = snakemake.output.plot_file
aspect = snakemake.wildcards.aspect

start_time = time.time()
print(f"--- [GLOBAL {aspect.upper()}] Plotting Individual Term Trajectories (6-Month Bins, Cool Colors, Legend) ---")

# Setup the figure for presentation
plt.figure(figsize=(14, 8))

# Lists to hold global data so we can still draw ONE bold average line on top
global_x = []
global_y = []

adj_files = list(adj_dir.glob("*_mean_adjacencies.parquet"))
print(f"Processing {len(adj_files)} terms...")

# Extract Tab20 colors and strictly KEEP only Blues, Greens, Purples, Grays, and Cyans
# Indices: 0,1 (Blue), 4,5 (Green), 8,9 (Purple), 14,15 (Gray), 18,19 (Cyan)
tab20 = plt.get_cmap('tab20').colors
safe_colors = [tab20[i] for i in [0, 1, 4, 5, 8, 9, 14, 15, 18, 19]]
color_cycle = itertools.cycle(safe_colors)

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
    
    # CALCULATE THE TRENDLINE FOR THIS SPECIFIC TERM (0.5 = 6 months)
    bins = np.arange(0, np.ceil(df['Delta_T_years'].max()) + 1, 0.5)
    bin_means, bin_edges, _ = binned_statistic(
        df['Delta_T_years'], df['PID0_mean_adj'], statistic='mean', bins=bins
    )
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    valid_bins = ~np.isnan(bin_means)
    
    # Only plot if this term has enough data to form a line (at least 2 points)
    if valid_bins.sum() > 1:
        # Layer 1: Plot this individual term with a safe cool color
        term_color = next(color_cycle)
        plt.plot(
            bin_centers[valid_bins], 
            bin_means[valid_bins], 
            color=term_color, 
            linewidth=1.5, 
            alpha=0.4,  
            zorder=1,
            label=term_str.replace('_', ':') # Add the formatted term ID to the legend
        )


# Layer 2: Plot the Global Average on top to anchor the visual (zorder=2)
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
        color='#d73027',  # The distinct red line
        linewidth=5.0, 
        label='Global Average Trend',
        zorder=2
    )

# Formatting
plt.title(f"Term-by-Term Network Mean Adjacency Score vs TTA\nGO Aspect: {aspect.upper()}", fontsize=28, fontweight='bold', pad=20)
plt.xlabel("Time to Annotation ($\\Delta$T in Years)", fontsize=22, labelpad=15)
plt.ylabel("Mean Adjacency Score", fontsize=22, labelpad=15)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Presentation Grid (zorder=0)
plt.grid(True, linestyle="--", alpha=0.4, zorder=0)

# Legend (Positioned outside the plot to the right in one column)
num_terms = len(adj_files)

if num_terms <= 20:
    # Normal behavior when only few terms
    cols = 1
elif num_terms <= 40:
    # Double columns
    cols = 2
else:
    # Remove labels if there are too many terms, showing only the global mean
    cols = 1
    handles, labels = plt.gca().get_legend_handles_labels()
    # Keep only the global average
    handles = [h for h, l in zip(handles, labels) if l == 'Global Average Trend']
    labels = ['Global Average Trend']
    plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=12)

# If we didn't deactivate the legend, then build it with the correct number of columns
if num_terms <= 40:
    plt.legend(
        loc="upper left", 
        bbox_to_anchor=(1.02, 1), 
        fontsize=10, 
        ncol=cols, 
        framealpha=0.9, 
        edgecolor='#cccccc'
    )

# Save
print(f"Saving spaghetti plot to {output_plot}...")
# bbox_inches='tight' is crucial so the external legend isn't cut off
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [GLOBAL {aspect.upper()}] Spaghetti plot generated in {elapsed_time} seconds! ---")