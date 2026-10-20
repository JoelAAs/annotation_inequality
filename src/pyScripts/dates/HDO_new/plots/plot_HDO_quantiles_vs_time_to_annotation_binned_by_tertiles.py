import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import time
import pickle
import warnings

warnings.filterwarnings('ignore')

adj_file = snakemake.input.mean_adj_file
dates_file = snakemake.input.annot_dates_file
network_file = snakemake.input.network_file 
output_plot_quantiles = snakemake.output.plot_quantiles
output_plot_degrees = snakemake.output.plot_degrees

term = snakemake.wildcards.term
term_formatted = term.replace("_", ":")

start_time = time.time()
print(f"--- [TERM {term_formatted} HDO] Generating Quantile and Degree Plots (Tertiles 33-33-33 Split) ---")

# Ensure output directories exist early
os.makedirs(os.path.dirname(output_plot_quantiles), exist_ok=True)
os.makedirs(os.path.dirname(output_plot_degrees), exist_ok=True)

# Helper function to generate placeholder plots if data is missing
def create_placeholder_plots(message):
    print(f"--- [WARNING] {message}. Creating placeholder plots for {term_formatted}. ---")
    
    for out_file, title in [(output_plot_quantiles, "Quantiles"), (output_plot_degrees, "Degrees")]:
        plt.figure(figsize=(14, 8))
        plt.text(0.5, 0.5, f"No valid temporal data available\nfor {term_formatted}\n({title})", 
                 horizontalalignment='center', verticalalignment='center', 
                 fontsize=18, transform=plt.gca().transAxes)
        plt.title(f"HDO Temporal Analysis - {term_formatted}", fontsize=22, fontweight='bold', pad=15)
        plt.axis('off')
        plt.savefig(out_file, dpi=300, bbox_inches='tight')
        plt.close()
    exit(0)

# Load Network
print("Loading static global network from pickle...")
try:
    with open(network_file, 'rb') as f:
        G = pickle.load(f)
except Exception as e:
    create_placeholder_plots(f"Failed to load network: {e}")

# Load Data with defensive checks
print("Loading adjacency and dates data...")
try:
    adj_df = pd.read_parquet(adj_file)
    if adj_df.empty or 'PID0_mean_adj' not in adj_df.columns:
        create_placeholder_plots("Matrix is empty")
    
    dates_df = pd.read_csv(dates_file, sep='\t')
    if dates_df.empty:
        create_placeholder_plots("Dates file is empty")
        
    dates_df = dates_df[['gene_id', 'first_annotation_date']]
except Exception as e:
    create_placeholder_plots(f"Error reading files: {e}")

adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
# FIX: Convert adj_df['Date'] to datetime BEFORE the merge to prevent ValueError
adj_df['Date'] = pd.to_datetime(adj_df['Date'].astype(str), format='%Y%m%d', errors='coerce')

dates_df['gene_id'] = dates_df['gene_id'].astype(str)

df = pd.merge(adj_df, dates_df, left_on='Future_Gene', right_on='gene_id', how='inner')

if df.empty:
    create_placeholder_plots("No matching genes after merge")

# Format dates and Apply Sanity Check safely
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d', errors='coerce')
df = df.dropna(subset=['Date', 'first_annotation_date'])

df = df[df['first_annotation_date'] > df['Date']].copy()

if df.empty:
    create_placeholder_plots("No valid future annotations available after date filtering")

# Calculate Historical Degree
EDGE_DATE_ATTR = 'discovery_date'

def get_historical_degree(row):
    node = str(row['Future_Gene'])
    if node not in G:
        return np.nan
        
    current_date_int = int(row['Date'].strftime('%Y%m%d'))
    historical_degree = 0
    
    for neighbor, edge_data in G[node].items():
        edge_date = edge_data.get(EDGE_DATE_ATTR, 99999999) 
        if edge_date <= current_date_int:
            historical_degree += 1
            
    return historical_degree

print("Calculating time-resolved historical degree for each target gene...")
df['Historical_Degree'] = df.apply(get_historical_degree, axis=1)
df = df.dropna(subset=['Historical_Degree'])

# --- DATA-DRIVEN TERTILE BINNING (33% - 33% - 33%) ---
df['Delta_T_days'] = (df['first_annotation_date'] - df['Date']).dt.days

t1 = df['Delta_T_days'].quantile(0.3333)
t2 = df['Delta_T_days'].quantile(0.6667)

print("\n--- Temporal Tertiles ---")
print(f"First Tertile: 0 -> {t1/365.25:.2f} years")
print(f"Middle Tertile: {t1/365.25:.2f} -> {t2/365.25:.2f} years")
print(f"Last Tertile: > {t2/365.25:.2f} years\n")

def assign_tertile_bin(days):
    if days <= t1:
        return 'First Tertile'
    elif days <= t2:
        return 'Middle Tertile'
    else:
        return 'Last Tertile'

df['TTA_Bin'] = df['Delta_T_days'].apply(assign_tertile_bin)
# ---------------------------------------------------

# Identify Permutation Columns
true_col = 'PID0_mean_adj'
# Ensure we only grab valid permutation columns
perm_cols = [col for col in adj_df.columns if col.startswith('PID') and col.endswith('_mean_adj') and col != true_col]

# If there are no permutations (e.g., if you generated the matrix without them), skip quantile calc safely
num_perms = len(perm_cols)
if num_perms == 0:
    create_placeholder_plots("No permutation columns found in the matrix for quantile calculation")

# Compute Quantiles and Mean Degree per Date/Bin
def compute_stats(group):    
    true_mean = group[true_col].mean()
    perm_means = group[perm_cols].mean()
    count_greater = (true_mean > perm_means).sum()
    quantile = (count_greater + 1) / (num_perms + 1)
    
    mean_degree = group['Historical_Degree'].mean()
    
    return pd.Series({
        'Quantile': quantile,
        'Mean_Degree': mean_degree,
        'Gene_Count': len(group)
    })

print("Aggregating quantiles and degrees per date and bin...")
agg_df = df.groupby(['Date', 'TTA_Bin']).apply(compute_stats).reset_index()

# Drop time-points with fewer than 5 genes
MIN_GENES = 5
agg_df = agg_df[agg_df['Gene_Count'] >= MIN_GENES].copy()

if agg_df.empty:
    create_placeholder_plots(f"No dates remained with >= {MIN_GENES} target genes")

# --- Prepare Global Stats Text Boxes ---
quantiles_stats_text = "Mean Quantile:\n"
degrees_stats_text = "Mean Degree:\n"

bin_order = ['First Tertile', 'Middle Tertile', 'Last Tertile']

for b in bin_order:
    # Safely calculate means handling cases where a bin might be completely empty
    subset = agg_df[agg_df['TTA_Bin'] == b]
    mean_q = subset['Quantile'].mean() if not subset.empty else np.nan
    mean_deg = subset['Mean_Degree'].mean() if not subset.empty else np.nan
    
    short_b = b.split(' ')[0] 
    quantiles_stats_text += f"• {short_b}: {mean_q:.3f}\n" if not pd.isna(mean_q) else f"• {short_b}: N/A\n"
    degrees_stats_text += f"• {short_b}: {mean_deg:.1f}\n" if not pd.isna(mean_deg) else f"• {short_b}: N/A\n"
    
quantiles_stats_text = quantiles_stats_text.strip()
degrees_stats_text = degrees_stats_text.strip()

# General Plotting Setup
sns.set_context("talk")
sns.set_style("whitegrid")

colors = {
    'First Tertile': '#7cb3e8',       
    'Middle Tertile': '#fce473',    
    'Last Tertile': '#f08b8b'      
}
props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='lightgray')

# PLOT 1: QUANTILES
plt.figure(figsize=(14, 8))

plotted_any = False
for b in bin_order:
    subset = agg_df[agg_df['TTA_Bin'] == b].sort_values('Date').copy()
    if subset.empty:
        continue
    subset.set_index('Date', inplace=True)
    subset['Smoothed_Quantile'] = subset['Quantile'].rolling('30D', min_periods=1).mean()
    subset.reset_index(inplace=True)
    
    plt.plot(subset['Date'], subset['Smoothed_Quantile'], color=colors[b], linestyle='-', linewidth=3.0, label=f'{b}')
    plotted_any = True

if not plotted_any:
    create_placeholder_plots("Data grouping failed to produce valid plot lines")

plt.title(f"Tertile-Binned Quantiles Over Time (HDO)\nTerm: {term_formatted}", fontsize=18, fontweight='bold', pad=15)
plt.xlabel("Year", fontsize=16, fontweight='bold', labelpad=10)
plt.ylabel("Quantile", fontsize=16, fontweight='bold', labelpad=10)
plt.ylim(0, 1.05)

ax1 = plt.gca()
ax1.xaxis.set_major_locator(mdates.YearLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xticks(rotation=45)

plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Random Expectation (0.5)')
plt.axhline(y=0.95, color='red', linestyle=':', alpha=0.5, label='p < 0.05 Threshold (0.95)')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title="Time to Annotation", fontsize=12, title_fontsize=14)
ax1.text(1.02, 0.4, quantiles_stats_text, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props, family='monospace')

sns.despine()
plt.tight_layout()
plt.savefig(output_plot_quantiles, dpi=300, bbox_inches='tight')
plt.close()

# PLOT 2: HISTORICAL DEGREES
plt.figure(figsize=(14, 8))

for b in bin_order:
    subset = agg_df[agg_df['TTA_Bin'] == b].sort_values('Date').copy()
    if subset.empty:
        continue
    subset.set_index('Date', inplace=True)
    subset['Smoothed_Degree'] = subset['Mean_Degree'].rolling('30D', min_periods=1).mean()
    subset.reset_index(inplace=True)
    
    plt.plot(subset['Date'], subset['Smoothed_Degree'], color=colors[b], linestyle='-', linewidth=3.0, label=f'{b}')

plt.title(f"Tertile-Binned Node Degree Average Over Time (HDO)\nTerm: {term_formatted}", fontsize=18, fontweight='bold', pad=15)
plt.xlabel("Year", fontsize=16, fontweight='bold', labelpad=10)
plt.ylabel("Average Node Degree", fontsize=16, fontweight='bold', labelpad=10)

ax2 = plt.gca()
ax2.xaxis.set_major_locator(mdates.YearLocator())
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xticks(rotation=45)

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title="Time to Annotation", fontsize=12, title_fontsize=14)
ax2.text(1.02, 0.4, degrees_stats_text, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props, family='monospace')

sns.despine()
plt.tight_layout()
plt.savefig(output_plot_degrees, dpi=300, bbox_inches='tight')
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- [TERM {term_formatted} HDO] Done! Both plots generated in {elapsed}s ---")