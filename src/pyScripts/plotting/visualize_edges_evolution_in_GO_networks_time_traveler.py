import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import pronto
import textwrap

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.edges_evolution_statistics
obo_file = snakemake.input.ontology
outputplot = snakemake.output.fractions_plot

go_ontology = pronto.Ontology(obo_file)
df = pd.read_csv(statistics_df, sep='\t')

plt.figure(figsize=(16, 10))

go_ids = df['GO_id'].unique()
colors = sns.color_palette("viridis", len(go_ids))

global_min_month = df['relative_month'].min()
global_max_month = df['relative_month'].max()

for i, go_id in enumerate(go_ids):
    go_df = df[df['GO_id'] == go_id].sort_values('relative_month').copy()
    
    if go_df.empty:
        continue

    go_df['cum_total'] = go_df['total_edges_added'].cumsum()
    go_df['cum_annot'] = go_df['annotated_edges_added'].cumsum()
    
    valid_mask = go_df['cum_total'] >= 15
    go_df = go_df[valid_mask]
    
    if len(go_df) < 3:
        continue
        
    go_df = go_df.set_index('relative_month')
    min_m = go_df.index.min()
    max_m = go_df.index.max()
    
    full_range = range(int(min_m), int(max_m) + 1)
    go_df = go_df.reindex(full_range).ffill()
    go_df['fraction'] = go_df['cum_annot'] / go_df['cum_total'].replace(0, np.nan)
    
    go_df = go_df.reset_index()
    go_df.rename(columns={'index': 'relative_month'}, inplace=True)
    
    x_vals = go_df['relative_month'].values / 12.0
    y_vals = go_df['fraction'].values
    
    try:
        go_name = go_ontology[go_id].name
    except KeyError:
        go_name = go_id
    wrapped_label = "\n".join(textwrap.wrap(f"{go_id}: {go_name}", width=30))
    
    plt.step(x_vals, y_vals, where='post', color=colors[i], alpha=0.25, linewidth=1.5, linestyle='--')
    
    lowess = sm.nonparametric.lowess(y_vals, x_vals, frac=0.2)
    plt.plot(lowess[:, 0], lowess[:, 1], label=wrapped_label, color=colors[i], linewidth=3.5, zorder=4)

plt.axvline(x=0, color='red', linestyle='-', linewidth=2.5, zorder=5, label="Annotation Event (Time 0)")

ax = plt.gca()

ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))

plt.yticks(fontsize=16)
plt.xticks(fontsize=12, rotation=45)

ax.grid(which='major', axis='x', linestyle='-', alpha=0.5, color='gray', zorder=0)
ax.grid(which='minor', axis='x', linestyle=':', alpha=0.3, color='gray', zorder=0)
ax.grid(which='major', axis='y', linestyle=':', alpha=0.6, zorder=0)

plt.xlim(global_min_month / 12.0, global_max_month / 12.0)

plt.title(f"Time Traveler GO {aspect} Lifespan Edges Evolution\n(Depth: {depth} - Cutoff: {cutoff})", fontsize=26, fontweight='bold', pad=25)
plt.ylabel("Cumulative Fraction (Annotated / Total)", fontsize=22, labelpad=20)
plt.xlabel("Timeline: Years Relative to Annotation Event\n(- = Before Annotation | + = After Annotation)", fontsize=22, labelpad=20)

legend = plt.legend(title="GO ID", bbox_to_anchor=(1.02, 1), loc='upper left', title_fontsize=16, fontsize=14)
legend.get_frame().set_alpha(0.8)

plt.tight_layout()
plt.savefig(outputplot, dpi=300, bbox_inches='tight')
plt.close()