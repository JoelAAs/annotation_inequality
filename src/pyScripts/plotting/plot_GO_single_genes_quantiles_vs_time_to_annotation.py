import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools

quantiles_dir = snakemake.input.quantiles_dir
annot_dates_dir = snakemake.input.annot_dates_dir
output_file = snakemake.output.plot_file

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff

start_time = time.time()
print(f"--- Starting Single Genes Quantiles vs Time-to-Annotation Plotting (Multi-Color) ---", flush=True)

# Gather all files in the directory
quantile_files = glob.glob(os.path.join(quantiles_dir, "*_quantiles.parquet"))
print(f"[LOAD] Found {len(quantile_files)} quantile files to process...", flush=True)

all_term_data = []

# Iterate through every GO term to process its data
for q_file in quantile_files:
    # Extract GO ID from filename
    basename = os.path.basename(q_file)
    term_id = basename.split('_single')[0] if '_single' in basename else "_".join(basename.split('_')[:2])
    term_formatted = term_id.replace('_', ':')
    
    # Locate the corresponding dates file
    date_file = os.path.join(annot_dates_dir, f"{term_id}_first_annotation_dates.csv")
    if not os.path.exists(date_file):
        continue
        
    # Load data
    q_df = pd.read_parquet(q_file)
    d_df = pd.read_csv(date_file, sep='\t')
    
    # Filter dates to this specific term just to be safe
    d_df = d_df[d_df['GO_id'] == term_formatted]
    
    # Format for merge
    q_df['Future_Gene'] = q_df['Future_Gene'].astype(str)
    d_df['gene_id'] = d_df['gene_id'].astype(str)
    
    # Merge on Gene ID
    merged = pd.merge(
        q_df[['Date', 'Future_Gene', 'Quantile']], 
        d_df[['gene_id', 'first_annotation_date']], 
        left_on='Future_Gene', 
        right_on='gene_id', 
        how='inner'
    )
    
    if merged.empty:
        continue
        
    merged['Date'] = pd.to_datetime(merged['Date'].astype(str), format='%Y%m%d')
    merged['first_annotation_date'] = pd.to_datetime(merged['first_annotation_date'].astype(str), format='%Y%m%d')
    
    # Filter strictly for future annotations
    merged = merged[merged['first_annotation_date'] > merged['Date']].copy()
    
    # Calculate Delta T in Years
    merged['Delta_T_Years'] = (merged['first_annotation_date'] - merged['Date']).dt.days / 365.25
    
    # Bin to the nearest half-year (0.5) to turn raw scatter points into drawable continuous lines
    merged['Delta_T_Binned'] = np.round(merged['Delta_T_Years'] * 2) / 2
    
    # Average the quantiles per time bin for this specific term
    term_agg = merged.groupby('Delta_T_Binned')['Quantile'].mean().reset_index()
    term_agg['Term'] = term_id
    
    all_term_data.append(term_agg)

print(f"[PROCESS] Aggregating data across all valid GO terms...", flush=True)
full_df = pd.concat(all_term_data, ignore_index=True)

# Calculate the global mean across all terms for the thick trendline
global_mean = full_df.groupby('Delta_T_Binned')['Quantile'].mean().reset_index()
# Smooth the global mean slightly to ensure it looks elegant
global_mean['Quantile_Smoothed'] = global_mean['Quantile'].rolling(window=3, min_periods=1, center=True).mean()


# PLOT
print(f"[PLOT] Drawing visualization...", flush=True)
plt.figure(figsize=(14, 8))

# Extract Tab20 colors and strictly keep cool/neutral colors so the red trendline stands out
tab20 = plt.get_cmap('tab20').colors
safe_colors = [tab20[i] for i in [0, 1, 4, 5, 8, 9, 14, 15, 18, 19]]
color_cycle = itertools.cycle(safe_colors)

# Draw the background: all individual GO terms
for term, group in full_df.groupby('Term'):
    term_color = next(color_cycle)
    plt.plot(
        group['Delta_T_Binned'], 
        group['Quantile'], 
        color=term_color, 
        alpha=0.4,
        linewidth=1.5, 
        zorder=1,
        label=term # Add the term ID to the legend
    )

# Draw the foreground: Global Mean Trendline
plt.plot(
    global_mean['Delta_T_Binned'], 
    global_mean['Quantile_Smoothed'], 
    color='#d73027', 
    linewidth=5.0, 
    label='Global Mean Trend', 
    zorder=2
)

# FORMATTING
plt.title(f"Predictive Quantiles vs. Time to Annotation (Single Genes)\nAspect: {aspect} | Depth: {depth} | Cutoff: {cutoff}", fontsize=22, fontweight='bold', pad=20)
plt.xlabel("Time to Annotation ($\\Delta$T in Years)", fontsize=18, labelpad=15)
plt.ylabel("Mean Predictive Quantile", fontsize=18, labelpad=15)

# Start X-axis cleanly at 0 years, lock Y-axis to probability bounds
plt.xlim(left=0)
plt.ylim(0, 1.05) 

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Place the legend outside the plot to the right, using one column
plt.legend(
    loc="upper left", 
    bbox_to_anchor=(1.02, 1), 
    fontsize=10, 
    ncol=1, 
    framealpha=0.9, 
    edgecolor='#cccccc'
)

plt.grid(True, linestyle='--', alpha=0.5, zorder=0)

# SAVE
print(f"[SAVE] Saving plot to disk...", flush=True)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
# bbox_inches='tight' is crucial here so the external legend isn't cut off in the saved file
plt.savefig(output_file, dpi=300, bbox_inches='tight') 
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- Plot successfully generated in {elapsed} seconds! ---", flush=True)