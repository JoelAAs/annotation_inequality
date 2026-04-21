import pandas as pd
import numpy as np
import pickle
import networkx as nx
from scipy.stats import wilcoxon

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
observed_df = snakemake.input.observed
baseline_df = snakemake.input.baseline
outputstats = snakemake.output.summary_stats

print(f"Processing top HDO coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline...")

print(f"Loading data for top HDO coefficients depth {depth} cutoff {cutoff}...")

with open(observed_df, 'rb') as f:
    obs_df = pickle.load(f)
with open(baseline_df, 'rb') as f:
    base_df = pickle.load(f)

print(f"Computing top HDO coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline...")

results = []

for i in range(len(obs_df)):
    annot = obs_df.iloc[i]['annotation']

    # Extract the list of neighbor bait count sums (matched 1 to 1)
    obs_vals = np.array([g['annotated_neighbors_count'] for g in obs_df.iloc[i]['genes_with_annotated_neighbors']])
    baseline_vals = np.array([g['annotated_neighbors_count'] for g in base_df.iloc[i]['genes_with_annotated_neighbors']])

    # Calculate statistics
    obs_mean = np.mean(obs_vals)
    baseline_mean = np.mean(baseline_vals)
    fold_change = obs_mean / baseline_mean if baseline_mean != 0 else np.inf

    # Wilcoxon signed-rank test (paired)
    try:
        stat, p_val = wilcoxon(obs_vals, baseline_vals, alternative = 'two-sided')
    except ValueError:
        stat, p_val = np.nan, 1.0

    # Saving results
    results.append({
        'HDO_annotation': annot,
        'N_genes': len(obs_vals),
        'Observed_mean': round(obs_mean, 3),
        'Baseline_mean': round(baseline_mean, 3),
        'Fold_enrichment': round(fold_change, 3),
        'Wilcoxon_p_value': p_val
    })

print(f"Saving results for top HDO coefficients annotated neighbors depth {depth} cutoff {cutoff}...")

summary_df = pd.DataFrame(results)
summary_df.to_csv(outputstats, sep = '\t', index = False)

print(f"Top HDO coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline ready!")