import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pronto
import textwrap

# --- Catching Wildcards ---
aspect = snakemake.wildcards.aspect if hasattr(snakemake.wildcards, 'aspect') else ""
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.edges_evolution_statistics
obo_file = snakemake.input.ontology
outputplot = snakemake.output.fractions_plot

# --- Load Data ---
ontology = pronto.Ontology(obo_file)
df = pd.read_csv(statistics_df, sep='\t')

plt.figure(figsize=(16, 10))

# Automatically detect if this is GO or HDO
id_col = 'GO_id' if 'GO_id' in df.columns else 'DO_id'
ids = df[id_col].unique()
colors = sns.color_palette("viridis", len(ids))

for i, term_id in enumerate(ids):
    term_df = df[df[id_col] == term_id].sort_values('relative_month').copy()
    
    if term_df.empty:
        continue

    # Filter strictly to the 1-year window
    window_mask = (term_df['relative_month'] >= -12) & (term_df['relative_month'] <= 12)
    window_df = term_df[window_mask]
    
    if len(window_df) < 3:
        continue 
        
    x_vals = window_df['relative_month'].values
    y_vals = window_df['mean_fraction'].values
    
    try:
        term_name = ontology[term_id].name
    except KeyError:
        term_name = term_id
    wrapped_label = "\n".join(textwrap.wrap(f"{term_id}: {term_name}", width=30))
    
    # Plotting
    plt.step(x_vals, y_vals, where='post', label=wrapped_label, color=colors[i], linewidth=3.5, zorder=4)
    plt.scatter(x_vals, y_vals, color=colors[i], s=50, alpha=1.0, zorder=5)

# --- Formatting ---
plt.axvline(x=0, color='red', linestyle='-', linewidth=2.5, zorder=1, label="Annotation Event (Month 0)")

plt.yticks(fontsize=16)
plt.xticks(ticks=range(-12, 13), labels=[f"{m}" if m != 0 else "0\n(Annotated)" for m in range(-12, 13)], fontsize=14)
plt.xlim(-12, 12)

onto_label = f"GO {aspect}" if id_col == 'GO_id' else "HDO"
plt.title(f"{onto_label} Edges Evolution (1-Year Window)\n(Depth: {depth} - Cutoff: {cutoff})", fontsize=24, fontweight='bold', pad=25)
plt.ylabel("Mean Annotated Neighborhood Fraction per Protein", fontsize=22, labelpad=20)
plt.xlabel("Months Relative to Annotation Event", fontsize=22, labelpad=20)

plt.grid(True, linestyle=':', alpha=0.6, zorder=0)

legend = plt.legend(title=id_col.replace('_id', ' ID'), bbox_to_anchor=(1.02, 1), loc='upper left', title_fontsize=16, fontsize=14)
legend.get_frame().set_alpha(0.8)

plt.tight_layout()
plt.savefig(outputplot, dpi=300, bbox_inches='tight')
plt.close()