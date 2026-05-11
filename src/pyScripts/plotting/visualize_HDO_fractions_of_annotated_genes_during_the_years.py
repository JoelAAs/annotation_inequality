import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm
import pronto
import textwrap

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.networks_statistics
obo_file = snakemake.input.ontology
outputplot = snakemake.output.fractions_plot
outputplot_general = snakemake.output.fractions_plot_general

do_ontology = pronto.Ontology(obo_file)

df = pd.read_csv(statistics_df, sep='\t')
df['date'] = pd.to_datetime(df['node_annotation_date'], format='%Y%m%d')

# --- PLOT 1: Past / Present / Future Step Plot ---
agg_df = df.groupby('DO_id')[[
    'annotated_neighbors_past', 'total_neighbors_past',
    'annotated_neighbors_present', 'total_neighbors_present',
    'annotated_neighbors_future', 'total_neighbors_future'
]].sum().reset_index()

agg_df['Past'] = agg_df['annotated_neighbors_past'] / agg_df['total_neighbors_past'].replace(0, np.nan)
agg_df['Present'] = agg_df['annotated_neighbors_present'] / agg_df['total_neighbors_present'].replace(0, np.nan)
agg_df['Future'] = agg_df['annotated_neighbors_future'] / agg_df['total_neighbors_future'].replace(0, np.nan)

plt.figure(figsize=(14, 10))

x_discrete = np.array([0, 1, 2])
x_dense = np.linspace(0, 2, 500)

do_ids = agg_df['DO_id'].unique()
colors = sns.color_palette("viridis", len(do_ids))
color_map = {do_id: colors[i] for i, do_id in enumerate(do_ids)}

for i, do_id in enumerate(do_ids):
    row = agg_df[agg_df['DO_id'] == do_id].iloc[0]
    y_discrete = np.array([row['Past'], row['Present'], row['Future']])

    try:
        do_name = do_ontology[do_id].name
    except KeyError:
        do_name = do_id

    wrapped_label = "\n".join(textwrap.wrap(f"{do_id}: {do_name}", width=30))
    valid_mask = ~np.isnan(y_discrete)

    if valid_mask.sum() == 3:
        plt.step(x_discrete, y_discrete, where='mid', color=color_map[do_id], alpha=0.3, linestyle='--')
        
        y_step_dense = np.piecewise(
            x_dense, 
            [x_dense < 0.5, (x_dense >= 0.5) & (x_dense < 1.5), x_dense >= 1.5], 
            [y_discrete[0], y_discrete[1], y_discrete[2]]
        )
        
        lowess = sm.nonparametric.lowess(y_step_dense, x_dense, frac=0.2, it=0)
        
        plt.plot(lowess[:, 0], lowess[:, 1], label=wrapped_label, color=color_map[do_id], linewidth=2.5, zorder=4)
        plt.scatter(x_discrete, y_discrete, color=color_map[do_id], s=50, zorder=5)

    elif valid_mask.sum() > 1:
        plt.step(x_discrete[valid_mask], y_discrete[valid_mask], where='mid', color=color_map[do_id], alpha=0.3, linestyle='--')
        plt.plot(x_discrete[valid_mask], y_discrete[valid_mask], label=wrapped_label, color=color_map[do_id], linewidth=2.5, zorder=4)
        plt.scatter(x_discrete[valid_mask], y_discrete[valid_mask], color=color_map[do_id], s=50, zorder=5)

plt.xticks(ticks=[0, 1, 2], labels=['Past (1 Year Prior)', 'Present (Annotation Date)', 'Future (1 Year After)'], fontsize=16)
plt.yticks(fontsize=16)
plt.title(f"HDO Trend of Annotated Neighbors Fractions\n(Depth: {depth} - Cutoff: {cutoff})", fontsize=26, fontweight='bold', pad=25)
plt.ylabel("Fraction (Sum Annotated / Sum Total)", fontsize=22, labelpad=20)
plt.xlabel("Time Period", fontsize=22, labelpad=20)
plt.grid(True, linestyle=':', alpha=0.6, zorder=0)
plt.legend(title="DOID", bbox_to_anchor=(1.05, 1), loc='upper left', title_fontsize=16, fontsize=14)
plt.tight_layout()

plt.savefig(outputplot, dpi=300, bbox_inches='tight')
plt.close()

# --- PLOT 2: Cumulative Progression Over Years ---
min_year = df['date'].dt.year.min()
max_year = df['date'].dt.year.max()
years = np.arange(min_year, max_year + 1)

cumulative_records = []

for year in years:
    cutoff_date = pd.Timestamp(year, 12, 31)
    subset_df = df[df['date'] <= cutoff_date]
    
    yearly_agg = subset_df.groupby('DO_id')[['annotated_neighbors_present', 'total_neighbors_present']].sum()
    
    for do_id in do_ids:
        if do_id in yearly_agg.index:
            ann_sum = yearly_agg.loc[do_id, 'annotated_neighbors_present']
            tot_sum = yearly_agg.loc[do_id, 'total_neighbors_present']
            frac = ann_sum / tot_sum if tot_sum > 0 else np.nan
        else:
            frac = np.nan
            
        cumulative_records.append({'Year': year, 'DO_id': do_id, 'Fraction': frac})

cum_df = pd.DataFrame(cumulative_records)

plt.figure(figsize=(14, 10))

for do_id in do_ids:
    do_data = cum_df[cum_df['DO_id'] == do_id].dropna(subset=['Fraction'])
    
    if not do_data.empty:
        try:
            do_name = do_ontology[do_id].name
        except KeyError:
            do_name = do_id
            
        wrapped_label = "\n".join(textwrap.wrap(f"{do_id}: {do_name}", width=30))
        
        plt.plot(do_data['Year'], do_data['Fraction'], marker='o', markersize=8, 
                 label=wrapped_label, color=color_map[do_id], linewidth=2.5, zorder=4)

plt.xticks(ticks=years, labels=[str(y) for y in years], fontsize=14, rotation=45)
plt.yticks(fontsize=16)
plt.title(f"HDO Cumulative Trend of Annotated Neighbors Fractions Over Time\n(Depth: {depth} - Cutoff: {cutoff})", fontsize=26, fontweight='bold', pad=25)
plt.ylabel("Cumulative Fraction (Annotated / Total)", fontsize=22, labelpad=20)
plt.xlabel("Year", fontsize=22, labelpad=20)
plt.grid(True, linestyle=':', alpha=0.6, zorder=0)
plt.legend(title="DOID", bbox_to_anchor=(1.05, 1), loc='upper left', title_fontsize=16, fontsize=14)
plt.tight_layout()

plt.savefig(outputplot_general, dpi=300, bbox_inches='tight')
plt.close()