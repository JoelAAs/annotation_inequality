import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
import warnings

warnings.filterwarnings('ignore')

adj_file = snakemake.input.mean_adj_file
dates_file = snakemake.input.annot_dates_file
output_plot = snakemake.output.plot_quantiles

term = snakemake.wildcards.term
term_formatted = term.replace("_", ":")

start_time = time.time()
print(f"--- [TERM {term_formatted}] Generating Continuous Sliding Window Plot ---")

# Load Data
adj_df = pd.read_parquet(adj_file)
dates_df = pd.read_csv(dates_file, sep='\t')[['gene_id', 'first_annotation_date']]

adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
dates_df['gene_id'] = dates_df['gene_id'].astype(str)

df = pd.merge(adj_df, dates_df, left_on='Future_Gene', right_on='gene_id', how='inner')

# Format dates and calculate continuous TTA
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')

df['Delta_T_days'] = (df['first_annotation_date'] - df['Date']).dt.days
df = df[df['Delta_T_days'] > 0].copy()

if df.empty:
    print("No valid data. Exiting.")
    os.makedirs(os.path.dirname(output_plot), exist_ok=True)
    open(output_plot, 'w').close()
    exit(0)

df['TTA_Years'] = df['Delta_T_days'] / 365.25

# Gene-Level Quantile Calculation
true_col = 'PID0_mean_adj'
perm_cols = [col for col in df.columns if col.startswith('PID') and col.endswith('_mean_adj') and col != true_col]
num_perms = len(perm_cols)

counts_greater = (df[perm_cols].values < df[[true_col]].values).sum(axis=1)
df['Gene_Quantile'] = (counts_greater + 1.0) / (num_perms + 1.0)

# Sort and apply Sliding Window
df = df.sort_values('TTA_Years').reset_index(drop=True)

# Window at 2%, or minimum 500 genes
window_size = max(500, int(len(df) * 0.02))

df['Rolling_Quantile'] = df['Gene_Quantile'].rolling(window=window_size, center=True).mean()
df['Rolling_TTA'] = df['TTA_Years'].rolling(window=window_size, center=True).mean()

plot_df = df.dropna(subset=['Rolling_Quantile', 'Rolling_TTA'])

# Plotting
plt.figure(figsize=(12, 7))
sns.set_context("talk")
sns.set_style("whitegrid")

color_q = '#4c72b0'

plt.plot(plot_df['Rolling_TTA'], plot_df['Rolling_Quantile'], color=color_q, linewidth=3.5, label='Mean Quantile')
plt.xlabel("Time to Annotation (Years)", fontsize=16, fontweight='bold', labelpad=10)
plt.ylabel("Mean Empirical Quantile", fontsize=16, fontweight='bold', labelpad=10)

plt.ylim(0.0, 1.05)

plt.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='Random Expectation')

plt.title(f"Predictive Power Landscape across Time to Annotation\nTerm: {term_formatted} (Window: {window_size} genes)", fontsize=18, fontweight='bold', pad=15)
plt.xlim(plot_df['Rolling_TTA'].min(), plot_df['Rolling_TTA'].max())
plt.legend(loc='lower right', fontsize=12)

sns.despine()
plt.tight_layout()
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()