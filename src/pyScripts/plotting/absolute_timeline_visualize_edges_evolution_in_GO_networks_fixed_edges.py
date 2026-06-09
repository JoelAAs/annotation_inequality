import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pronto
import textwrap
import matplotlib.dates as mdates

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.edges_evolution_statistics
obo_file = snakemake.input.ontology
outputplot = snakemake.output.fractions_plot

print(f"Plotting absolute timeline edges evolution with fixed edges for GO {aspect} network...")

ontology = pronto.Ontology(obo_file)
df = pd.read_csv(statistics_df, sep='\t')

df['exact_date'] = pd.to_datetime(df['exact_date'])

plt.figure(figsize=(16, 10))

# Automatically detect if this is GO or HDO
id_col = 'GO_id' if 'GO_id' in df.columns else 'DO_id'
ids = df[id_col].unique()
colors = sns.color_palette("viridis", len(ids))

for i, term_id in enumerate(ids):
    term_df = df[df[id_col] == term_id].sort_values('exact_date').copy()
    
    if len(term_df) < 3:
        continue 
        
    x_vals = term_df['exact_date']
    y_vals = term_df['mean_fraction'].values
    
    try:
        term_name = ontology[term_id].name
    except KeyError:
        term_name = term_id
        
    wrapped_label = "\n".join(textwrap.wrap(f"{term_id}: {term_name}", width=30))
    
    plt.step(x_vals, y_vals, where='post', label=wrapped_label, color=colors[i], linewidth=3.5, zorder=4)
    plt.scatter(x_vals, y_vals, color=colors[i], s=20, alpha=0.8, zorder=5)

plt.yticks(fontsize=16)
plt.xticks(fontsize=14, rotation=90, ha='center')

# Smart X-Axis Formatting: Matplotlib will automatically figure out how to label 
# the axis with Years (e.g., 2000, 2005, 2010) cleanly without overlapping text.
locator = mdates.AutoDateLocator(minticks=15, maxticks=30)
ax = plt.gca()
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))

onto_label = f"GO {aspect}" if id_col == 'GO_id' else "HDO"

plt.title(f"Historical {onto_label} Network Evolution\n(Fixed Edges | Daily Precision | Depth: {depth} - Cutoff: {cutoff})", fontsize=24, fontweight='bold', pad=25)
plt.ylabel("Mean Annotated Neighborhood Fraction", fontsize=22, labelpad=20)
plt.xlabel("Timeline", fontsize=22, labelpad=20)

plt.grid(True, linestyle=':', alpha=0.6, zorder=0)

legend = plt.legend(title=id_col.replace('_id', ' ID'), bbox_to_anchor=(1.02, 1), loc='upper left', title_fontsize=16, fontsize=14)
legend.get_frame().set_alpha(0.8)

plt.tight_layout()
plt.savefig(outputplot, dpi=300, bbox_inches='tight')
plt.close()

print(f"Absolute timeline edges evolution plot with fixed edges for GO {aspect} network ready!")