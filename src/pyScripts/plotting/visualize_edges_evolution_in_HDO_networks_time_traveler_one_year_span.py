import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import pronto
import textwrap

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.edges_evolution_statistics
obo_file = snakemake.input.ontology
outputplot = snakemake.output.fractions_plot

do_ontology = pronto.Ontology(obo_file)
df = pd.read_csv(statistics_df, sep='\t')

plt.figure(figsize=(16, 10))

do_ids = df['DO_id'].unique()
colors = sns.color_palette("viridis", len(do_ids))

for i, do_id in enumerate(do_ids):
    do_df = df[df['DO_id'] == do_id].sort_values('relative_month').copy()
    
    if do_df.empty:
        continue

    do_df['cum_total'] = do_df['total_edges_added'].cumsum()
    do_df['cum_annot'] = do_df['annotated_edges_added'].cumsum()
    
    do_df = do_df.set_index('relative_month')
    min_m = do_df.index.min()
    max_m = do_df.index.max()
    
    full_range = range(int(min_m), int(max_m) + 1)
    do_df = do_df.reindex(full_range).ffill()
    do_df['fraction'] = do_df['cum_annot'] / do_df['cum_total'].replace(0, np.nan)
    
    do_df = do_df.reset_index()
    do_df.rename(columns={'index': 'relative_month'}, inplace=True)

    window_mask = (do_df['relative_month'] >= -12) & (do_df['relative_month'] <= 12)
    window_df = do_df[window_mask]
    
    if len(window_df) < 3:
        continue 
        
    x_vals = window_df['relative_month'].values
    y_vals = window_df['fraction'].values
    
    try:
        do_name = do_ontology[do_id].name
    except KeyError:
        do_name = do_id
    wrapped_label = "\n".join(textwrap.wrap(f"{do_id}: {do_name}", width=30))
    
    plt.step(x_vals, y_vals, where='post', color=colors[i], alpha=0.3, linewidth=1.5, linestyle='--')
    
    plt.scatter(x_vals, y_vals, color=colors[i], s=30, alpha=0.6, zorder=5)
    
    lowess = sm.nonparametric.lowess(y_vals, x_vals, frac=0.3)
    plt.plot(lowess[:, 0], lowess[:, 1], label=wrapped_label, color=colors[i], linewidth=3.5, zorder=4)

plt.axvline(x=0, color='red', linestyle='-', linewidth=2.5, zorder=1, label="Annotation Event (Month 0)")

plt.yticks(fontsize=16)
plt.xticks(ticks=range(-12, 13), labels=[f"{m}" if m != 0 else "0\n(Annotated)" for m in range(-12, 13)], fontsize=14)

plt.xlim(-12, 12)

plt.title(f"Time-Traveler HDO Edges Evolution\n(Depth: {depth} - Cutoff: {cutoff})", fontsize=24, fontweight='bold', pad=25)
plt.ylabel("Cumulative Fraction (Annotated / Total)", fontsize=22, labelpad=20)
plt.xlabel("Months Relative to Annotation Event", fontsize=22, labelpad=20)

plt.grid(True, linestyle=':', alpha=0.6, zorder=0)

legend = plt.legend(title="DOID", bbox_to_anchor=(1.02, 1), loc='upper left', title_fontsize=16, fontsize=14)
legend.get_frame().set_alpha(0.8)

plt.tight_layout()
plt.savefig(outputplot, dpi=300, bbox_inches='tight')
plt.close()