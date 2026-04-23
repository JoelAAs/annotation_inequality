import pandas as pd
import pickle
import textwrap
import seaborn as sns
import matplotlib.pyplot as plt
import pronto
import networkx as nx

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
stats_df = snakemake.input.summary_stats
observed_df = snakemake.input.observed
baseline_df = snakemake.input.baseline
go = snakemake.input.ontology
outputplots = snakemake.output.summary_stats_plots

print(f"Processing top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline plotting...")

print(f"Loading top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline data...")

stats = pd.read_csv(stats_df, sep  = '\t')
with open(observed_df, 'rb') as f:
    obs_df = pickle.load(f)
with open(baseline_df, 'rb') as f:
    base_df = pickle.load(f)
ontology = pronto.Ontology(go)

print(f"Setting up top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline data...")

def get_go_id_name(goid):
    try:
        return ontology[goid].name
    except KeyError:
        return goid 
    
stats['GO_term'] = stats['GO_annotation'].apply(get_go_id_name)

dist_data = []

for i in range(len(obs_df)):
    annot = obs_df.iloc[i]['annotation']
    name = get_go_id_name(annot)
    obs_vals = [g['annotated_neighbors_count'] for g in obs_df.iloc[i]['genes_with_annotated_neighbors']]
    base_vals = [g['annotated_neighbors_count'] for g in base_df.iloc[i]['genes_with_annotated_neighbors']]

    for v in obs_vals: dist_data.append({
        'GO_term': name,
        'Group': 'Observed',
        'Value': v
    })
    for v in base_vals: dist_data.append({
        'GO_term': name,
        'Group': 'Baseline',
        'Value': v
    })

dist_df = pd.DataFrame(dist_data)

def format_name(name):
    return "\n".join(textwrap.wrap(name, 20))

stats['GO_term_plot'] = stats['GO_term'].apply(format_name)
dist_df['Annotation_plot'] = dist_df['GO_term'].apply(format_name)

print(f"Plotting top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline...")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (12, 14), sharex = True)
sns.set_style('whitegrid')

plot_bar_df = stats.melt(
    id_vars = ['GO_term', 'Wilcoxon_p_value'],
    value_vars = ['Observed_mean', 'Baseline_mean'],
    var_name = 'Group',
    value_name = 'Mean'
)
plot_bar_df['Group'] = plot_bar_df['Group'].str.replace('_mean', '')

sns.barplot(data = plot_bar_df,
            x = 'GO_term',
            y = 'Mean',
            hue = 'Group',
            palette = 'Set2',
            ax = ax1,
            edgecolor='0.3',
            width = 0.7,
            zorder=3
)

max_val = stats['Observed_mean'].max()
ax1.set_ylim(0, max_val * 1.4)
ax1.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.3)

for i, row in stats.iterrows():
    p = row['Wilcoxon_p_value']
    stars = "****" if p < 0.0001 else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    p_text = f"{stars}\n(p = {p:.1e})"
    h = max(row['Observed_mean'], row['Baseline_mean'])
    ax1.text(i, h + (max_val * 0.02), p_text, ha = 'center', va = 'bottom', fontsize = 16, fontweight = 'bold')

ax1.set_title(f"A. Mean GO {aspect} Annotated Neighbors (Paired Wilcoxon)", loc = 'left', fontweight = 'bold', fontsize = 26, pad = 25)
ax1.set_ylabel("Mean count", fontsize = 22, labelpad= 20)
ax1.tick_params(axis='y', labelsize=18)
ax1.legend(title = "Group", title_fontsize=20, fontsize=18, loc='upper right')

sns.violinplot(data = dist_df,
               x = 'GO_term',
               y = 'Value',
               hue = 'Group',
               split = True,
               inner = 'quart',
               palette = 'Set2',
               ax = ax2,
               linewidth=2
)

ax2.set_yscale('log')
ax2.set_title(f"B. Distribution of GO {aspect} Annotated Neighbors (Log Scale)", loc = 'left', fontweight='bold', fontsize = 26, pad = 25)
ax2.set_ylabel("Annotated neighbors (log scale)", fontsize = 22, labelpad = 20)
ax2.set_xlabel("GO Annotation", fontsize = 22, labelpad = 20)
ax2.set_xticklabels(
    [format_name(t.get_text()) for t in ax2.get_xticklabels()],
    rotation=45, 
    ha='right',          
    rotation_mode='anchor', 
    fontsize=18, 
    fontweight='medium'
)
ax2.tick_params(axis='y', labelsize=18)
ax2.legend(title = "Group", title_fontsize=20, fontsize=18, loc='upper right')

plt.tight_layout()

print(f"Saving Top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline plots...")

plt.subplots_adjust(bottom=0.2, hspace=0.3)
plt.savefig(outputplots, dpi = 300, bbox_inches='tight')

print(f"Top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline plotting done!")