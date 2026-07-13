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

print(f"Plotting top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline...")

# Switched to a single plot layout
fig, ax = plt.subplots(figsize = (14, 10))
sns.set_style('whitegrid')

# Plot the split violin plot
sns.violinplot(
    data = dist_df,
    x = 'GO_term',
    y = 'Value',
    hue = 'Group',
    split = True,
    inner = 'quart',
    palette = 'Set2',
    ax = ax,
    linewidth=2,
    order=stats['GO_term'] # Ensures the x-axis order exactly matches the iteration order for p-values
)

ax.set_yscale('log')

# Expand the upper y-limit to make room for the p-value text on a log scale
max_dist_val = dist_df['Value'].max()
ax.set_ylim(bottom=None, top=max_dist_val * 20)

# Iterate through the stats to add p-value annotations
for i, row in stats.iterrows():
    p = row['Wilcoxon_p_value']
    stars = "****" if p < 0.0001 else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    p_text = f"{stars}\n(p = {p:.1e})"
    
    # Find the maximum value for this specific GO term to place text right above the violin
    term_max = dist_df[dist_df['GO_term'] == row['GO_term']]['Value'].max()
    
    # Multiply for offset due to the logarithmic scale
    h = term_max * 2.5 if term_max > 0 else 1.5 
    
    ax.text(i, h, p_text, ha='center', va='bottom', fontsize=14, fontweight='bold')

# Removed "B." from title since it's the only plot now
ax.set_title(f"Distribution of GO {aspect} Annotated Neighbors (Log Scale)", loc = 'left', fontweight='bold', fontsize = 26, pad = 25)
ax.set_ylabel("Annotated neighbors (log scale)", fontsize = 22, labelpad = 20)
ax.set_xlabel("GO Annotation", fontsize = 22, labelpad = 20)

# Apply formatting to xticklabels
ax.set_xticklabels(
    [format_name(t.get_text()) for t in ax.get_xticklabels()],
    rotation=45, 
    ha='right',          
    rotation_mode='anchor', 
    fontsize=18, 
    fontweight='medium'
)
ax.tick_params(axis='y', labelsize=18)
ax.legend(title = "Group", title_fontsize=20, fontsize=18, loc='upper right')

plt.tight_layout()

print(f"Saving Top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline plots...")

plt.savefig(outputplots, dpi = 300, bbox_inches='tight')

print(f"Top GO {aspect} coefficients depth {depth} cutoff {cutoff} annotated neighbors correlation with baseline plotting done!")