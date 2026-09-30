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
print(f"--- [TERM {term_formatted}] Generating Quantile and Degree Plots ---")

# Load Network
print("Loading static global network from pickle...")
with open(network_file, 'rb') as f:
    G = pickle.load(f)

# Load and Merge Data
print("Loading adjacency and dates data...")
adj_df = pd.read_parquet(adj_file)
dates_df = pd.read_csv(dates_file, sep='\t')[['gene_id', 'first_annotation_date']]

adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
dates_df['gene_id'] = dates_df['gene_id'].astype(str)

df = pd.merge(adj_df, dates_df, left_on='Future_Gene', right_on='gene_id', how='inner')

# Format dates and Apply Sanity Check
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')

df = df[df['first_annotation_date'] > df['Date']].copy()

# Safety measure in case of empty data
if df.empty:
    print("No valid data available after date filtering. Exiting.")
    os.makedirs(os.path.dirname(output_plot_quantiles), exist_ok=True)
    os.makedirs(os.path.dirname(output_plot_degrees), exist_ok=True)
    open(output_plot_quantiles, 'w').close()
    open(output_plot_degrees, 'w').close()
    exit(0)

# Calculate Historical Degree
EDGE_DATE_ATTR = 'discovery_date'

def get_historical_degree(row):
    node = str(row['Future_Gene'])
    if node not in G:
        return np.nan
        
    current_date_int = int(row['Date'].strftime('%Y%m%d'))
    historical_degree = 0
    
    for neighbor, edge_data in G[node].items():
        # Obtain edge date, if missing we assume it is later in the future
        edge_date = edge_data.get(EDGE_DATE_ATTR, 99999999) 
        if edge_date <= current_date_int:
            historical_degree += 1
            
    return historical_degree

print("Calculating time-resolved historical degree for each target gene...")
df['Historical_Degree'] = df.apply(get_historical_degree, axis=1)
df = df.dropna(subset=['Historical_Degree'])

# Compute Time to Annotation and Binning
df['Delta_T_days'] = (df['first_annotation_date'] - df['Date']).dt.days

def assign_bin(days):
    if days < 365:
        return '< 1 Year'
    elif days < 1825:
        return '1 - 5 Years'
    else:
        return '>= 5 Years'

df['TTA_Bin'] = df['Delta_T_days'].apply(assign_bin)

# Identify Permutation Columns
true_col = 'PID0_mean_adj'
perm_cols = [col for col in df.columns if col.startswith('PID') and col.endswith('_mean_adj') and col != true_col]
if not perm_cols:
    raise ValueError("Error: No permutation columns found!")
num_perms = len(perm_cols)

# Compute Quantiles and Mean Degree per Date/Bin
def compute_stats(group):
    # Quantile Calculation
    true_mean = group[true_col].mean()
    perm_means = group[perm_cols].mean()
    count_greater = (true_mean > perm_means).sum()
    quantile = (count_greater + 1) / (num_perms + 1)
    
    # Historical Degree Calculation
    mean_degree = group['Historical_Degree'].mean()
    
    return pd.Series({
        'Quantile': quantile,
        'Mean_Degree': mean_degree,
        'Gene_Count': len(group)
    })

print("Aggregating quantiles and degrees per date and bin...")
agg_df = df.groupby(['Date', 'TTA_Bin']).apply(compute_stats).reset_index()

# FILTER: Drop time-points with fewer than 5 genes
MIN_GENES = 5
agg_df = agg_df[agg_df['Gene_Count'] >= MIN_GENES].copy()

# --- Prepare Global Stats Text Boxes ---
quantiles_stats_text = "Mean Quantile:\n"
degrees_stats_text = "Mean Degree:\n"

for b in ['< 1 Year', '1 - 5 Years', '>= 5 Years']:
    mean_q = agg_df[agg_df['TTA_Bin'] == b]['Quantile'].mean()
    mean_deg = agg_df[agg_df['TTA_Bin'] == b]['Mean_Degree'].mean()
    
    quantiles_stats_text += f"• {b}: {mean_q:.3f}\n" if not pd.isna(mean_q) else f"• {b}: N/A\n"
    degrees_stats_text += f"• {b}: {mean_deg:.1f}\n" if not pd.isna(mean_deg) else f"• {b}: N/A\n"
    
quantiles_stats_text = quantiles_stats_text.strip()
degrees_stats_text = degrees_stats_text.strip()

# General Plotting Setup
sns.set_context("talk")
sns.set_style("whitegrid")

bin_order = ['< 1 Year', '1 - 5 Years', '>= 5 Years']
colors = {
    '< 1 Year': '#7cb3e8',       
    '1 - 5 Years': '#fce473',    
    '>= 5 Years': '#f08b8b'      
}
props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='lightgray')

os.makedirs(os.path.dirname(output_plot_quantiles), exist_ok=True)
os.makedirs(os.path.dirname(output_plot_degrees), exist_ok=True)

# PLOT 1: QUANTILES
plt.figure(figsize=(14, 8))

for b in bin_order:
    subset = agg_df[agg_df['TTA_Bin'] == b].sort_values('Date').copy()
    if subset.empty:
        continue
    subset.set_index('Date', inplace=True)
    subset['Smoothed_Quantile'] = subset['Quantile'].rolling('30D', min_periods=1).mean()
    subset.reset_index(inplace=True)
    
    plt.plot(subset['Date'], subset['Smoothed_Quantile'], color=colors[b], linestyle='-', linewidth=3.0, label=f'{b}')

plt.title(f"Binned Quantiles Over Time\nTerm: {term_formatted}", fontsize=18, fontweight='bold', pad=15)
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

plt.title(f"Binned Node Degree Average Over Time\nTerm: {term_formatted}", fontsize=18, fontweight='bold', pad=15)
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
print(f"--- [TERM {term_formatted}] Done! Both plots generated in {elapsed}s ---")