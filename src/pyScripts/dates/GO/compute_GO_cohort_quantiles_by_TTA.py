import os
import pandas as pd
import numpy as np
import time

input_adj_file = snakemake.input.mean_adj_file
input_dates_file = snakemake.input.annot_dates_file
output_file = snakemake.output.quantile_file
term = snakemake.wildcards.term
term_formatted = term.replace('_', ':')

print(f"\n--- [Core] Computing Cohort Quantiles (Date + TTA Grouped) for {term} ---\n", flush=True)
start_time = time.time()

# LOAD DATA
print(f"[{term}] Loading mean adjacency and annotation date files...")
adj_df = pd.read_parquet(input_adj_file)
dates_df = pd.read_csv(input_dates_file, sep='\t')

# Filter dates specifically for this term
dates_df = dates_df[dates_df['GO_id'] == term_formatted]

# MERGE & CALCULATE DELTA T
print(f"[{term}] Merging and calculating Time-to-Annotation...")
adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
dates_df['gene_id'] = dates_df['gene_id'].astype(str)

df = pd.merge(
    adj_df, 
    dates_df[['gene_id', 'first_annotation_date']], 
    left_on='Future_Gene', 
    right_on='gene_id', 
    how='inner'
)

# Parse dates
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')

# Strictly filter for future annotations
df = df[df['first_annotation_date'] > df['Date']].copy()

if df.empty:
    print(f"[{term}] No valid future annotations found. Creating empty output.")
    empty_df = pd.DataFrame(columns=['Date', 'Delta_T_Binned', 'Group_Size', 'True_Mean_Adj', 'Quantile'])
    empty_df.to_parquet(output_file, index=False)
    exit(0)

# Calculate Delta T in years
df['Delta_T_Years'] = (df['first_annotation_date'] - df['Date']).dt.days / 365.25

# Bin to exactly 3 biological phases: <1 yr, 1-5 yrs, >=5 yrs
bins = [0, 1, 5, np.inf]
labels = ['< 1 year', '1-5 years', '>= 5 years']

df['Delta_T_Binned'] = pd.cut(
    df['Delta_T_Years'], 
    bins=bins, 
    labels=labels, 
    right=False  # Ensures the intervals are [0,1), [1,5), [5, inf)
).astype(str)

# DUAL-ISOLATION (GROUP BY DATE & TTA BIN)
print(f"[{term}] Grouping by Date and TTA Bin to compute cohort means...")

# Track how many genes fall into each specific Date/TTA bucket
group_sizes = df.groupby(['Date', 'Delta_T_Binned']).size().reset_index(name='Group_Size')

# Calculate the mean for PID0 and all 1000 decoys simultaneously for every bucket
numeric_cols = [f"PID{i}_mean_adj" for i in range(1001)]
grouped_means = df.groupby(['Date', 'Delta_T_Binned'])[numeric_cols].mean().reset_index()

# Merge the group sizes back into the grouped dataframe
grouped_means = pd.merge(grouped_means, group_sizes, on=['Date', 'Delta_T_Binned'])

# COMPUTE THE QUANTILE
print(f"[{term}] Computing statistical quantiles against 1,000 decoys...")
true_vals = grouped_means['PID0_mean_adj'].values
decoy_cols = [f"PID{i}_mean_adj" for i in range(1, 1001)]
decoy_mat = grouped_means[decoy_cols].values

# Vectorized quantile math: How many decoy means are smaller than the true mean?
quantiles = (decoy_mat < true_vals[:, None]).sum(axis=1) / 1000.0

# ASSEMBLE FINAL DATAFRAME
results_df = pd.DataFrame({
    'Date': grouped_means['Date'],
    'Delta_T_Binned': grouped_means['Delta_T_Binned'],
    'Group_Size': grouped_means['Group_Size'],
    'True_Mean_Adj': true_vals.astype(np.float32),
    'Quantile': quantiles.astype(np.float32)
})

# Save to disk
os.makedirs(os.path.dirname(output_file), exist_ok=True)
results_df.to_parquet(output_file, engine='pyarrow', index=False)

elapsed = round(time.time() - start_time, 2)
print(f"--- [Core] Successfully saved Cohort Quantiles to {output_file} in {elapsed}s ---\n", flush=True)