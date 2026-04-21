import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pronto
import textwrap

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
stats_df = snakemake.input.summary_stats
go = snakemake.input.ontology
outputtable = snakemake.output.summary_stats_table
outputplot = snakemake.output.fold_enrichment_plot

print(f"Plotting stats for GO {aspect} top coefficients depth {depth} cutoff {cutoff}...")

print(f"Loading data for GO {aspect} top coefficients depth {depth} cutoff {cutoff}...")

stats = pd.read_csv(stats_df, sep = '\t')
ontology = pronto.Ontology(go)

print(f"Plotting GO {aspect} top coefficients depth {depth} cutoff {cutoff} stats table...")

def get_goid_name(goid):
    try:
        return ontology[goid].name
    except KeyError:
        return goid 
    
stats['GO_annotation'] = stats['GO_annotation'].apply(get_goid_name)
stats['GO_annotation'] = stats['GO_annotation'].apply(lambda x: "\n".join(textwrap.wrap(x, 20)))

fig, ax = plt.subplots(figsize = (12, 6))
ax.axis('off')

table = ax.table(cellText = stats.values, 
                 colLabels = stats.columns, 
                 cellLoc = 'center', 
                 loc = 'center')

table.auto_set_font_size(False)
table.set_fontsize(9)
table.auto_set_column_width(col=list(range(len(stats.columns))))
table.scale(1.2, 3.5)

for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_text_props(weight='bold', color='white')
        cell.set_facecolor('#4c72b0')

plt.title(f"Top Coefficients GO {aspect} Neighborhood Comparison Statistics", fontweight = 'bold', pad = 30)

print(f"Saving stats table for GO {aspect} top coefficients depth {depth} cutoff {cutoff}...")

plt.savefig(outputtable, bbox_inches = 'tight', dpi = 300)
plt.close()

print(f"Plotting GO {aspect} top coefficients depth {depth} cutoff {cutoff} fold change plot...")

plt.figure(figsize=(12, 6))
sns.set_style("whitegrid")

stats = stats.sort_values('Fold_enrichment', ascending = False)

ax = sns.barplot(data = stats, x = 'GO_annotation', y = 'Fold_enrichment', palette='viridis')

# Add a dashed line at 1.0 (null hypothesis)
plt.axhline(1.0, color = 'red', linestyle = '--', alpha = 0.6)

# Annotate with p-values
for i, p in enumerate(stats['Wilcoxon_p_value']):
    ax.text(i, stats['Fold_enrichment'].iloc[i] + 0.01, f"p = {p:.1e}", 
            ha = 'center', fontweight = 'bold', color = 'black')

plt.title(f"Fold Enrichment of GO {aspect} Top Coefficients Neighbor Bait Counts (Observed vs Baseline)")
plt.ylabel("Fold Enrichment ($Obs/Base$)")
plt.xlabel("GO Annotation")
plt.tight_layout()

print(f"Saving fold change plot for GO {aspect} top coefficients depth {depth} cutoff {cutoff}...")

plt.savefig(outputplot)

print(f"Stats plots for GO {aspect} top coefficients depth {depth} cutoff {cutoff} ready!")