import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import time
import warnings

# Suppress seaborn layout warnings for clean logs
warnings.filterwarnings('ignore')

adj_file = snakemake.input.mean_adj_file
dates_file = snakemake.input.annot_dates_file
output_plot = snakemake.output.plot_file
term = snakemake.wildcards.term
term_formatted = term.replace("_", ":")

start_time = time.time()
print(f"--- [TERM {term_formatted}] Plotting Binned TTA Lineplot ---")

# Load and Merge Data
adj_df = pd.read_parquet(adj_file)
dates_df = pd.read_csv(dates_file, sep='\t')[['gene_id', 'first_annotation_date']]

adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
dates_df['gene_id'] = dates_df['gene_id'].astype(str)

df = pd.merge(adj_df, dates_df, left_on='Future_Gene', right_on='gene_id', how='inner')

# Format dates and Apply Sanity Check
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')

initial_rows = len(df)
df = df[df['first_annotation_date'] > df['Date']].copy()
dropped_rows = initial_rows - len(df)

if dropped_rows > 0:
    print(f"Sanity Check Warning: Dropped {dropped_rows} rows where annotation was already present.")

if df.empty:
    print("No valid data available after date filtering. Exiting.")
    os.makedirs(os.path.dirname(output_plot), exist_ok=True)
    open(output_plot, 'w').close()
    exit(0)

# Compute Time to Annotation (Delta T) and Binning
df['Delta_T_days'] = (df['first_annotation_date'] - df['Date']).dt.days

def assign_bin(days):
    if days < 365:
        return '< 1 Year'
    elif days < 1825:
        return '1 - 5 Years'
    else:
        return '>= 5 Years'

df['TTA_Bin'] = df['Delta_T_days'].apply(assign_bin)

# Handle Permutations
true_col = 'PID0_mean_adj'
perm_cols = [col for col in df.columns if col.startswith('PID') and col.endswith('_mean_adj') and col != true_col]

if not perm_cols:
    print("Warning: No permutation columns found in the dataset!")
    df['Perm_mean_adj'] = np.nan
else:
    df['Perm_mean_adj'] = df[perm_cols].mean(axis=1)

# True - Perm difference mean
df['True_vs_Perm_Diff'] = df[true_col] - df['Perm_mean_adj']
diff_stats = df.groupby('TTA_Bin')['True_vs_Perm_Diff'].mean()

d_short = diff_stats.get('< 1 Year', np.nan)
d_med = diff_stats.get('1 - 5 Years', np.nan)
d_long = diff_stats.get('>= 5 Years', np.nan)

stats_text = "Mean Δ (True - Permuted):\n"
stats_text += f"• < 1 Year:   {d_short:+.2e}\n" if not pd.isna(d_short) else "• < 1 Year:   N/A\n"
stats_text += f"• 1-5 Years:  {d_med:+.2e}\n" if not pd.isna(d_med) else "• 1-5 Years:  N/A\n"
stats_text += f"• >= 5 Years: {d_long:+.2e}" if not pd.isna(d_long) else "• >= 5 Years: N/A"

# Aggregate data per Date and Bin
agg_df = df.groupby(['Date', 'TTA_Bin']).agg(
    True_Adj_Mean=(true_col, 'mean'),
    Perm_Adj_Mean=('Perm_mean_adj', 'mean'),
    Gene_Count=(true_col, 'size')
).reset_index()

MIN_GENES = 5
agg_df = agg_df[agg_df['Gene_Count'] >= MIN_GENES].copy()

# Plotting Setup
plt.figure(figsize=(14, 8))
sns.set_context("talk")
sns.set_style("whitegrid")

bin_order = ['< 1 Year', '1 - 5 Years', '>= 5 Years']
colors = {
    '< 1 Year': '#7cb3e8',       
    '1 - 5 Years': '#fce473',    
    '>= 5 Years': '#f08b8b'      
}

# Draw the 6 Lines with Temporal Rolling Average
for b in bin_order:
    subset = agg_df[agg_df['TTA_Bin'] == b].sort_values('Date').copy()
    
    if subset.empty:
        continue
        
    subset.set_index('Date', inplace=True)
    subset['Smoothed_True'] = subset['True_Adj_Mean'].rolling('30D', min_periods=1).mean()
    subset['Smoothed_Perm'] = subset['Perm_Adj_Mean'].rolling('30D', min_periods=1).mean()
    subset.reset_index(inplace=True)
    
    plt.plot(
        subset['Date'], 
        subset['Smoothed_True'], 
        color=colors[b], 
        linestyle='-', 
        linewidth=2.5,
        label=f'True Adj ({b})'
    )
    
    if not subset['Smoothed_Perm'].isna().all():
        plt.plot(
            subset['Date'], 
            subset['Smoothed_Perm'], 
            color=colors[b], 
            linestyle='--', 
            linewidth=2,
            alpha=0.85,
            label=f'Perm Mean ({b})'
        )

# Plot Formatting
plt.title(f"Predictive Power Over Time by TTA Bins\nTerm: {term_formatted}", fontsize=18, fontweight='bold', pad=15)
plt.xlabel("Year", fontsize=16, fontweight='bold', labelpad=10)
plt.ylabel("Mean Adjacency Score", fontsize=16, fontweight='bold', labelpad=10)

ax = plt.gca()
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xticks(rotation=45)

# Main legend
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title="Condition & Bin", fontsize=12, title_fontsize=14)

# Stats differences box
props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='lightgray')
ax.text(1.02, 0.5, stats_text, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props, family='monospace')

sns.despine()
plt.tight_layout()

# Save and Close
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- [TERM {term_formatted}] Binned Lineplot generated in {elapsed}s! ---")