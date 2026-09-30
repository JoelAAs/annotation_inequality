import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import time
import warnings

warnings.filterwarnings('ignore')

adj_file = snakemake.input.mean_adj_file
dates_file = snakemake.input.annot_dates_file

output_super_predictors = snakemake.output.plot_super_predictors
output_density_2d = snakemake.output.plot_density_2d

term = snakemake.wildcards.term
term_formatted = term.replace("_", ":")

start_time = time.time()
print(f"--- [TERM {term_formatted}] Plotting Super-Predictors & 2D Continuous Density ---")

# Load and Merge Data
adj_df = pd.read_parquet(adj_file)
dates_df = pd.read_csv(dates_file, sep='\t')[['gene_id', 'first_annotation_date']]

adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
dates_df['gene_id'] = dates_df['gene_id'].astype(str)

df = pd.merge(adj_df, dates_df, left_on='Future_Gene', right_on='gene_id', how='inner')

# Format dates and Apply Sanity Check
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')

df = df[df['first_annotation_date'] > df['Date']].copy()

if df.empty:
    print("No valid data available after date filtering. Exiting.")
    os.makedirs(os.path.dirname(output_super_predictors), exist_ok=True)
    os.makedirs(os.path.dirname(output_density_2d), exist_ok=True)
    open(output_super_predictors, 'w').close()
    open(output_density_2d, 'w').close()
    exit(0)

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

# Gene-Level Quantile Calculation (Vettorizzato)
true_col = 'PID0_mean_adj'
perm_cols = [col for col in df.columns if col.startswith('PID') and col.endswith('_mean_adj') and col != true_col]

if not perm_cols:
    raise ValueError("Error: No permutation columns found in the dataset!")

num_perms = len(perm_cols)

print("Computing gene-level empirical quantiles...")
# Confronto vettorializzato riga per riga: quante colonne di permutazione sono inferiori a true_col
counts_greater = (df[perm_cols].values < df[[true_col]].values).sum(axis=1)
df['Gene_Quantile'] = (counts_greater + 1.0) / (num_perms + 1.0)

# PLOT 1: SUPER-PREDICTORS OVER TIME (>= 0.95)
SUPER_THRESHOLD = 0.95
df['Is_Super_Predictor'] = df['Gene_Quantile'] >= SUPER_THRESHOLD

print("Aggregating super-predictor fractions...")
agg_super = df.groupby(['Date', 'TTA_Bin']).agg(
    Fraction_Super=('Is_Super_Predictor', 'mean'),
    Gene_Count=('Future_Gene', 'size')
).reset_index()

# Minimum genes filter
MIN_GENES = 5
agg_super = agg_super[agg_super['Gene_Count'] >= MIN_GENES].copy()

# Box Text: Global Mean Super-Predictors Fraction
global_stats_text = f"Mean Fraction >= {SUPER_THRESHOLD}:\n"
for b in ['< 1 Year', '1 - 5 Years', '>= 5 Years']:
    mean_frac = agg_super[agg_super['TTA_Bin'] == b]['Fraction_Super'].mean()
    global_stats_text += f"• {b}: {mean_frac*100:.1f}%\n" if not pd.isna(mean_frac) else f"• {b}: N/A\n"
global_stats_text = global_stats_text.strip()

plt.figure(figsize=(14, 8))
sns.set_context("talk")
sns.set_style("whitegrid")

bin_order = ['< 1 Year', '1 - 5 Years', '>= 5 Years']
colors = {
    '< 1 Year': '#7cb3e8',       
    '1 - 5 Years': '#fce473',    
    '>= 5 Years': '#f08b8b'      
}

for b in bin_order:
    subset = agg_super[agg_super['TTA_Bin'] == b].sort_values('Date').copy()
    if subset.empty:
        continue
    subset.set_index('Date', inplace=True)
    subset['Smoothed_Frac'] = subset['Fraction_Super'].rolling('30D', min_periods=1).mean()
    subset.reset_index(inplace=True)
    
    plt.plot(
        subset['Date'], 
        subset['Smoothed_Frac'] * 100, 
        color=colors[b], 
        linestyle='-', 
        linewidth=3.0, 
        label=f'{b}'
    )

plt.title(f"Percentage of High-Confidence Target Genes (Quantile >= {SUPER_THRESHOLD})\nTerm: {term_formatted}", fontsize=18, fontweight='bold', pad=15)
plt.xlabel("Year", fontsize=16, fontweight='bold', labelpad=10)
plt.ylabel("Genes with Empirical Quantile >= 0.95 (%)", fontsize=16, fontweight='bold', labelpad=10)
plt.ylim(-2, 102)

ax1 = plt.gca()
ax1.xaxis.set_major_locator(mdates.YearLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xticks(rotation=45)

plt.axhline(y=5.0, color='gray', linestyle='--', alpha=0.5, label='Random Expectation (5%)')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title="Time to Annotation", fontsize=12, title_fontsize=14)

props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='lightgray')
ax1.text(1.02, 0.4, global_stats_text, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props, family='monospace')

sns.despine()
plt.tight_layout()
os.makedirs(os.path.dirname(output_super_predictors), exist_ok=True)
plt.savefig(output_super_predictors, dpi=300, bbox_inches='tight')
plt.close()

# PLOT 2: 2D CONTINUOUS DENSITY (HEXBIN MULTI-PANEL)
print("Generating 2D continuous density landscape...")
fig, axes = plt.subplots(1, 3, figsize=(22, 7), sharey=True, sharex=True)

# Convert dates in numeric coordinates for the hexbin
df['Date_Num'] = mdates.date2num(df['Date'])
date_min, date_max = df['Date_Num'].min(), df['Date_Num'].max()

for idx, b in enumerate(bin_order):
    ax = axes[idx]
    sub = df[df['TTA_Bin'] == b]
    
    if sub.empty:
        ax.set_title(f"{b} (No data)", fontsize=16, fontweight='bold')
        continue
        
    # Hexbin 2D with log scale
    hb = ax.hexbin(
        sub['Date_Num'], 
        sub['Gene_Quantile'], 
        gridsize=(40, 25), 
        cmap='magma', 
        bins='log',
        mincnt=1,
        extent=[date_min, date_max, 0, 1.05]
    )
    
    ax.set_title(f"Time to Annotation: {b}\n(N genes = {len(sub):,})", fontsize=16, fontweight='bold', pad=10)
    ax.set_xlabel("Year", fontsize=14, fontweight='bold', labelpad=8)
    if idx == 0:
        ax.set_ylabel("Empirical Quantile (Per Gene)", fontsize=14, fontweight='bold', labelpad=8)
        
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.tick_params(axis='x', rotation=45)
    ax.set_ylim(-0.02, 1.05)
    
    # 0.95 threshold line
    ax.axhline(0.95, color='cyan', linestyle=':', linewidth=1.5, alpha=0.8)
    
    cb = fig.colorbar(hb, ax=ax, orientation='horizontal', pad=0.2, shrink=0.7)
    cb.set_label('Log10(Gene Count)', fontsize=11)

plt.suptitle(f"Gene-Level Quantile Density Landscape Over Time\nTerm: {term_formatted}", fontsize=20, fontweight='bold', y=1.05)
sns.despine()
os.makedirs(os.path.dirname(output_density_2d), exist_ok=True)
plt.savefig(output_density_2d, dpi=300, bbox_inches='tight')
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- [TERM {term_formatted}] Done! Both exploratory plots generated in {elapsed}s ---")