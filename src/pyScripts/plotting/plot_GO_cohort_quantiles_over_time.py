import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import time

quantiles_dir = snakemake.input.quantiles_dir
output_file = snakemake.output.plot_file
aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff

start_time = time.time()
print(f"--- Starting Chronological Evolution Plot for Aspect: {aspect.upper()} ---", flush=True)

# 1. Gather and Load Data
quantile_files = glob.glob(os.path.join(quantiles_dir, "*_cohort_tta_quantiles.parquet"))
print(f"[LOAD] Found {len(quantile_files)} cohort quantile files to process...", flush=True)

all_data = []

for q_file in quantile_files:
    df = pd.read_parquet(q_file)
    if not df.empty:
        all_data.append(df)

full_df = pd.concat(all_data, ignore_index=True)

# Ensure Date is a proper datetime object for chronological plotting
full_df['Date'] = pd.to_datetime(full_df['Date'])
# Ensure bins are strings to avoid categorical grouping issues
full_df['Delta_T_Binned'] = full_df['Delta_T_Binned'].astype(str)

# 2. Aggregate across all GO Terms using a Size-Weighted Average
print(f"[PROCESS] Aggregating temporal trends using weighted averages...", flush=True)

def weighted_quantile(group):
    total_size = group['Group_Size'].sum()
    if total_size == 0:
        return np.nan
    return (group['Quantile'] * group['Group_Size']).sum() / total_size

# Apply the weighted average to properly anchor the baseline
global_trend = full_df.groupby(['Date', 'Delta_T_Binned']).apply(weighted_quantile).reset_index(name='Quantile')

# 3. Plotting Setup
print(f"[PLOT] Drawing visualization...", flush=True)
plt.figure(figsize=(14, 8))
sns.set_context("talk")
sns.set_style("white")

# Define specific, presentation-safe colors for your 3 bins
bin_colors = {
    '< 1 year': '#d73027', 
    '1-5 years': '#74add1', 
    '>= 5 years': '#313695'
}

ordered_bins = ['< 1 year', '1-5 years', '>= 5 years']

# 4. Draw the Trendlines (12-Month Rolling Window)
for bin_label in ordered_bins:
    # Extract temporal data, ensuring chronological order
    subset = global_trend[global_trend['Delta_T_Binned'] == bin_label].sort_values('Date').copy()
    
    if not subset.empty:
        # Apply a 12-snapshot centered rolling average to absorb curation bursts
        subset['Quantile_Smoothed'] = subset['Quantile'].rolling(
            window=24, 
            min_periods=3, 
            center=True
        ).mean()
        
        plt.plot(
            subset['Date'], 
            subset['Quantile_Smoothed'], 
            label=bin_label, 
            color=bin_colors.get(bin_label, 'black'),
            linewidth=4.0,
            zorder=3
        )

# 5. Formatting
plt.title(f"Evolution of Network Predictive Power Over Time\nAspect: {aspect.upper()} | Depth: {depth} | Cutoff: {cutoff}", 
          fontsize=22, fontweight='bold', pad=20)
plt.xlabel("Network Snapshot Date", fontsize=18, labelpad=15, fontweight='bold')
plt.ylabel("Mean Cohort Quantile Score (Weighted)", fontsize=18, labelpad=15, fontweight='bold')

# Y-axis bounds for probability
plt.ylim(0.0, 1.05)

# Format X-axis to cleanly display years (e.g., "2010", "2015")
plt.gca().xaxis.set_major_locator(mdates.YearLocator(base=2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xticks(fontsize=14, rotation=45, ha='right')
plt.yticks(fontsize=14)

# Legend configuration
plt.legend(
    title="Time to Annotation ($\\Delta$T)",
    title_fontsize=14,
    fontsize=14,
    loc="lower right",
    framealpha=0.9,
    edgecolor='#cccccc'
)

# Clean borders and grid
sns.despine()
plt.grid(True, axis='y', linestyle='-', alpha=0.3, color='gray', zorder=0)
plt.grid(True, axis='x', linestyle='--', alpha=0.2, color='gray', zorder=0)

# 6. Save
print(f"[SAVE] Saving plot to disk...", flush=True)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- Plot successfully generated in {elapsed} seconds! ---", flush=True)