import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pronto
import textwrap
from matplotlib.ticker import MaxNLocator

# Catching wildcards
aspect = snakemake.wildcards.aspect if hasattr(snakemake.wildcards, 'aspect') else ""
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.edges_evolution_statistics
obo_file = snakemake.input.ontology
outputplot = snakemake.output.fractions_plot

ontology = pronto.Ontology(obo_file)
df = pd.read_csv(statistics_df, sep='\t')

plt.figure(figsize=(16, 10))

# Automatically detect if this is GO or HDO based on the dataframe columns
id_col = 'GO_id' if 'GO_id' in df.columns else 'DO_id'

ids = df[id_col].unique()
colors = sns.color_palette("viridis", len(ids))

# Find global temporal limits to set the X-axis boundaries properly
global_min_month = df['relative_month'].min()
global_max_month = df['relative_month'].max()

for i, term_id in enumerate(ids):
    # Data is already computed as mean_fraction per month!
    term_df = df[df[id_col] == term_id].sort_values('relative_month').copy()
    
    if len(term_df) < 3:
        continue 
        
    x_vals = term_df['relative_month'].values
    y_vals = term_df['mean_fraction'].values
    
    try:
        term_name = ontology[term_id].name
    except KeyError:
        term_name = term_id
        
    wrapped_label = "\n".join(textwrap.wrap(f"{term_id}: {term_name}", width=30))
    
    # 1. Plot the main step line (slightly thinner for the zoomed-out view)
    plt.step(x_vals, y_vals, where='post', label=wrapped_label, color=colors[i], linewidth=2.5, zorder=4)
    
    # 2. Plot the scatter dots (made smaller since the timeline is now much longer)
    plt.scatter(x_vals, y_vals, color=colors[i], s=20, alpha=0.8, zorder=5)

# Add the red reference line for when the annotation actually happened
plt.axvline(x=0, color='red', linestyle='-', linewidth=2.5, zorder=1, label="Annotation Event (Month 0)")

# Format the axes
plt.yticks(fontsize=16)
plt.xticks(fontsize=14)

# Force the x-axis to show sensible, clean intervals automatically
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True, nbins=15))

# Set the limits to the absolute min and max of your entire dataset
plt.xlim(global_min_month, global_max_month)

# Dynamic Title depending on the ontology
onto_label = f"GO {aspect}" if id_col == 'GO_id' else "HDO"
plt.title(f"Time Traveler {onto_label} Edges Evolution\n(Depth: {depth} - Cutoff: {cutoff})", fontsize=24, fontweight='bold', pad=25)
plt.ylabel("Mean Annotated Neighborhood Fraction per Protein", fontsize=22, labelpad=20)
plt.xlabel("Months Relative to Annotation Event", fontsize=22, labelpad=20)

plt.grid(True, linestyle=':', alpha=0.6, zorder=0)

legend = plt.legend(title=id_col.replace('_id', ' ID'), bbox_to_anchor=(1.02, 1), loc='upper left', title_fontsize=16, fontsize=14)
legend.get_frame().set_alpha(0.8)

plt.tight_layout()
plt.savefig(outputplot, dpi=300, bbox_inches='tight')
plt.close()