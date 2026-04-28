import pandas as pd
import scipy.stats as stats
import numpy as np
import pronto

# 1. Snakemake variables
aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.networks_statistics
go_path = snakemake.input.ontology
output_stats = snakemake.output.correlation_results

print(f"Calculating temporal statistics (Spearman & Wilcoxon) for GO {aspect} depth {depth} cutoff {cutoff}...")

print(f"Loading data for GO {aspect} depth {depth} cutoff {cutoff}...")

df = pd.read_csv(statistics_df, sep='\t')
go = pronto.Ontology(go_path)

results = []

print(f"Starting temporal statistics computation for GO {aspect} depth {depth} cutoff {cutoff}...")

# Analyze each GO term individually
for go_id, group in df.groupby('GO_id'):
    go_name = go[go_id].name if go_id in go else go_id
    
    # ==========================================
    # TOTAL NEIGHBORS
    # ==========================================
    past_total = group['total_neighbors_past'].values
    present_total = group['total_neighbors_present'].values
    future_total = group['total_neighbors_future'].values
    
    # -- Spearman (Trend across all 3 timepoints) --
    total_counts_seq = np.concatenate([past_total, present_total, future_total])
    time_seq = np.concatenate([np.repeat(1, len(past_total)), 
                               np.repeat(2, len(present_total)), 
                               np.repeat(3, len(future_total))])
    total_rho, total_spearman_p = stats.spearmanr(time_seq, total_counts_seq)
    
    # -- Wilcoxon (Magnitude diff Past vs Future) --
    total_diff = future_total - past_total
    total_mean_diff = np.mean(total_diff)
    
    # -- Log2 Fold Change (Pseudocounts added to avoid Division by Zero) --
    total_log2fc = np.mean(np.log2((future_total + 1) / (past_total + 1)))
    
    if np.all(total_diff == 0):
        total_wilcoxon_p = 1.0
    else:
        _, total_wilcoxon_p = stats.wilcoxon(past_total, future_total, alternative='two-sided')

    # ==========================================
    # ANNOTATED NEIGHBORS
    # ==========================================
    past_annot = group['annotated_neighbors_past'].values
    present_annot = group['annotated_neighbors_present'].values
    future_annot = group['annotated_neighbors_future'].values
    
    # -- Spearman --
    annot_counts_seq = np.concatenate([past_annot, present_annot, future_annot])
    annot_rho, annot_spearman_p = stats.spearmanr(time_seq, annot_counts_seq)
    
    # -- Wilcoxon --
    annot_diff = future_annot - past_annot
    annot_mean_diff = np.mean(annot_diff)
    
    # -- Log2 Fold Change --
    annot_log2fc = np.mean(np.log2((future_annot + 1) / (past_annot + 1)))
    
    if np.all(annot_diff == 0):
        annot_wilcoxon_p = 1.0
    else:
        _, annot_wilcoxon_p = stats.wilcoxon(past_annot, future_annot, alternative='two-sided')

    # ==========================================
    # STORE RESULTS
    # ==========================================
    results.append({
        'GO_id': go_id,
        'GO_Name': go_name,
        # Total Stats
        'Total_Spearman_Rho': round(total_rho, 4),
        'Total_Spearman_P': '{:.2e}'.format(total_spearman_p),
        'Total_Wilcoxon_MeanDiff': round(total_mean_diff, 4),
        'Total_Log2FC': round(total_log2fc, 4),
        'Total_Wilcoxon_P': '{:.2e}'.format(total_wilcoxon_p),
        # Annotated Stats
        'Annot_Spearman_Rho': round(annot_rho, 4),
        'Annot_Spearman_P': '{:.2e}'.format(annot_spearman_p),
        'Annot_Wilcoxon_MeanDiff': round(annot_mean_diff, 4),
        'Annot_Log2FC': round(annot_log2fc, 4),
        'Annot_Wilcoxon_P': '{:.2e}'.format(annot_wilcoxon_p)
    })

print(f"Saving temporal statistics for GO {aspect} depth {depth} cutoff {cutoff}...")

results_df = pd.DataFrame(results)
results_df.to_csv(output_stats, sep='\t', index=False)

print(f"Temporal statistics results successfully saved to {output_stats}!")

print(f"Temporal statistics for GO {aspect} depth {depth} cutoff {cutoff} ready!")