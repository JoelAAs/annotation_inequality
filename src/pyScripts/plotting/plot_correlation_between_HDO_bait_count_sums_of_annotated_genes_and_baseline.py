import pandas as pd
import pickle
import textwrap
import seaborn as sns
import matplotlib.pyplot as plt
import pronto

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
stats_df = snakemake.input.summary_stats
observed_df = snakemake.input.observed
baseline_df = snakemake.input.baseline
do = snakemake.input.ontology
outputplots = snakemake.output.summary_stats_plots

print(f"Processing top HDO coefficients depth {depth} cutoff {cutoff} correlation with baseline plotting...")

print(f"Loading top HDO coefficients depth {depth} cutoff {cutoff} correlation with baseline data...")

stats = pd.read_csv(stats_df, sep  = '\t')
with open(observed_df, 'rb') as f:
    obs_df = pickle.load(f)
with open(baseline_df, 'rb') as f:
    base_df = pickle.load(f)
ontology = pronto.Ontology(do)

print(f"Setting up top HDO coefficients depth {depth} cutoff {cutoff} correlation with baseline data...")

def get_doid_name(doid):
    try:
        return ontology[doid].name
    except KeyError:
        return doid 
    
stats['HDO_term'] = stats['HDO_annotation'].apply(get_doid_name)

dist_data = []

for i in range(len(obs_df)):
    annot = obs_df.iloc[i]['annotation']
    name = get_doid_name(annot)
    obs_vals = [g['neighbor_sum'] for g in obs_df.iloc[i]['genes_with_sums']]
    base_vals = [g['neighbor_sum'] for g in base_df.iloc[i]['genes_with_sums']]

    for v in obs_vals: dist_data.append({
        'HDO_term': name,
        'Group': 'Observed',
        'Value': v
    })
    for v in base_vals: dist_data.append({
        'HDO_term': name,
        'Group': 'Baseline',
        'Value': v
    })

dist_df = pd.DataFrame(dist_data)

def format_name(name):
    return "\n".join(textwrap.wrap(name, 20))

stats['HDO_term_plot'] = stats['HDO_term'].apply(format_name)
dist_df['Annotation_plot'] = dist_df['HDO_term'].apply(format_name)

print(f"Plotting top HDO coefficients depth {depth} cutoff {cutoff} correlation with baseline...")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (12, 14), sharex = True)
sns.set_style('whitegrid')

plot_bar_df = stats.melt(
    id_vars = ['HDO_term', 'Wilcoxon_p_value'],
    value_vars = ['Observed_mean', 'Baseline_mean'],
    var_name = 'Group',
    value_name = 'Mean'
)
plot_bar_df['Group'] = plot_bar_df['Group'].str.replace('_mean', '')

sns.barplot(data = plot_bar_df,
            x = 'HDO_term',
            y = 'Mean',
            hue = 'Group',
            palette = 'Set2',
            ax = ax1)

max_val = stats['Observed_mean'].max()
ax1.set_ylim(0, max_val * 1.25)

for i, row in stats.iterrows():
    p = row['Wilcoxon_p_value']
    stars = "****" if p < 0.0001 else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    p_text = f"{stars}\n(p = {p:.1e})"
    h = max(row['Observed_mean'], row['Baseline_mean'])
    ax1.text(i, h + (max_val * 0.02), p_text, ha = 'center', va = 'bottom', fontsize = 10, fontweight = 'bold')

ax1.set_title(f"A. Mean HDO Neighbor Bait Counts (Paired Wilcoxon)", loc = 'left', fontweight = 'bold')
ax1.set_ylabel("Mean Sum")
ax1.legend(title = "Group")

sns.violinplot(data = dist_df,
               x = 'HDO_term',
               y = 'Value',
               hue = 'Group',
               split = True,
               inner = 'quart',
               palette = 'Set2',
               ax = ax2)

ax2.set_yscale('log')
ax2.set_title(f"B. Distribution of HDO Neighbor Sums (Log Scale)", loc = 'left', fontweight='bold')
ax2.set_ylabel("Neighbors Bait Count Sum (log - scale)")
ax2.set_xlabel("HDO Annotation")
ax2.set_xticklabels([format_name(t.get_text()) for t in ax2.get_xticklabels()])
ax2.legend(title = "Group")

plt.tight_layout()

print(f"Saving Top HDO coefficients depth {depth} cutoff {cutoff} correlation with baseline plots...")

plt.savefig(outputplots, dpi = 300)

print(f"Top HDO coefficients depth {depth} cutoff {cutoff} correlation with baseline plotting done!")