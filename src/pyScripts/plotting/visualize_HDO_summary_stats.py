import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pronto
import textwrap

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
stats_df = snakemake.input.summary_stats
do = snakemake.input.ontology
outputtable = snakemake.output.summary_stats_table
outputplot = snakemake.output.fold_enrichment_plot

print(f"Plotting stats for HDO top coefficients depth {depth} cutoff {cutoff}...")

print(f"Loading data for HDO top coefficients depth {depth} cutoff {cutoff}...")

stats = pd.read_csv(stats_df, sep = '\t')
ontology = pronto.Ontology(do)

print(f"Plotting HDO top coefficients depth {depth} cutoff {cutoff} stats table...")

def get_doid_name(doid):
    try:
        return ontology[doid].name
    except KeyError:
        return doid 
    
stats['HDO_annotation'] = stats['HDO_annotation'].apply(get_doid_name)
stats['HDO_annotation'] = stats['HDO_annotation'].apply(lambda x: "\n".join(textwrap.wrap(x, 20)))

fig, ax = plt.subplots(figsize = (12, 6))
ax.axis('off')

table = ax.table(cellText = stats.values, 
                 colLabels = stats.columns, 
                 cellLoc = 'center', 
                 loc = 'center')

table.auto_set_font_size(False)
table.set_fontsize(9)
table.auto_set_column_width(col=list(range(len(stats.columns))))
table.scale(1.2, 3.0)

for (row, col), cell in table.get_celld().items():
    if row == 0:
        cell.set_text_props(weight='bold', color='white')
        cell.set_facecolor('#4c72b0')

plt.title("Top Coefficients HDO Neighborhood Comparison Statistics", fontweight = 'bold', pad = 30)

print(f"Saving stats table for HDO top coefficients depth {depth} cutoff {cutoff}...")

plt.savefig(outputtable, bbox_inches = 'tight', dpi = 300)
plt.close()

print(f"Plotting HDO top coefficients depth {depth} cutoff {cutoff} fold change plot...")

plt.figure(figsize=(16, 8))
sns.set_style("whitegrid")

stats = stats.sort_values('Fold_enrichment', ascending = False)

ax = sns.barplot(data = stats, x = 'HDO_annotation', y = 'Fold_enrichment', palette='viridis', edgecolor='0.2', linewidth=2)

max_y = stats['Fold_enrichment'].max()
ax.set_ylim(0, max_y * 1.3)

# Add a dashed line at 1.0 (null hypothesis)
plt.axhline(1.0, color = 'red', linestyle = '--', alpha = 0.8)

# Annotate with p-values
for i, p in enumerate(stats['Wilcoxon_p_value']):
    ax.text(i, stats['Fold_enrichment'].iloc[i] + (max_y * 0.03), f"p = {p:.1e}", 
            ha = 'center', va='bottom', fontweight = 'bold', color = 'black', fontsize = 18)

plt.title("Fold Enrichment of HDO Top Coefficients Neighbor Bait Counts (Observed vs Baseline)", fontsize=26, fontweight='bold', pad=25)
plt.ylabel("Fold Enrichment ($Obs/Base$)", fontsize=22, labelpad=20)
plt.xlabel("HDO Annotation", fontsize=22, labelpad=20)
plt.xticks(fontsize = 20, fontweight = 'medium')
plt.yticks(fontsize = 20, fontweight = 'medium')
plt.tight_layout()

print(f"Saving fold change plot for HDO top coefficients depth {depth} cutoff {cutoff}...")

plt.savefig(outputplot, dpi=300, bbox_inches='tight')

print(f"Stats plots for HDO top coefficients depth {depth} cutoff {cutoff} ready!")