import pandas as pd
import networkx as nx
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pronto
import textwrap

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
statistics_df = snakemake.input.networks_statistics
go_path = snakemake.input.ontology
outputplot_total = snakemake.output.network_statistics_total_plots
outputplot_annot = snakemake.output.network_statistics_annotated_plots

print(f"Plotting GO {aspect} depth {depth} cutoff {cutoff} networks statistics...")
print(f"Loading GO {aspect} depth {depth} cutoff {cutoff} networks statistics data...")

df = pd.read_csv(statistics_df, sep = '\t')

print(f"Loading GO ontology from local file...")
go = pronto.Ontology(go_path)

print(f"Translating GO {aspect} depth {depth} cutoff {cutoff} IDs to Names...")

# Wrap at 15 characters to stack them horizontally without overlapping
df['GO_Name'] = df['GO_id'].apply(lambda x: go[x].name if x in go else x)
df['GO_Name'] = df['GO_Name'].apply(lambda x: textwrap.fill(x, width=15))

print(f"Processing GO {aspect} depth {depth} cutoff {cutoff} networks statistics plots...")

## --- Melt Dataframes ---
# Total neighbors
melted_total = df.melt(
    id_vars = ['GO_Name', 'node_id'],
    value_vars = ['total_neighbors_past', 'total_neighbors_present', 'total_neighbors_future'],
    var_name = 'Time_Point',
    value_name = 'Count'
)

time_map_total = {
    'total_neighbors_past': 'Year Prior',
    'total_neighbors_present': 'Discovery Year',
    'total_neighbors_future': 'Year After'
}

melted_total['Time_Point'] = melted_total['Time_Point'].map(time_map_total)
melted_total['Count_Plot'] = melted_total['Count'] + 1 ## Pseudocounts to avoid log(0) values

# Annotated_neighbors
melted_annot = df.melt(
    id_vars = ['GO_Name', 'node_id'],
    value_vars = ['annotated_neighbors_past', 'annotated_neighbors_present', 'annotated_neighbors_future'],
    var_name = 'Time_Point',
    value_name = 'Count'
)

time_map_annot = {
    'annotated_neighbors_past': 'Year Prior',
    'annotated_neighbors_present': 'Discovery Year',
    'annotated_neighbors_future': 'Year After'
}

melted_annot['Time_Point'] = melted_annot['Time_Point'].map(time_map_annot)
melted_annot['Count_Plot'] = melted_annot['Count'] + 1

## --- Plot settings ---
sns.set_theme(style='whitegrid')
palette = sns.color_palette('viridis', 3)

## --- Figure 1 (Total Neighbors) ---
fig1, axes1 = plt.subplots(2, 1, figsize=(12, 14))

# Top: Bar Plot
sns.barplot(data=melted_total, x='GO_Name', y='Count', hue='Time_Point', ax=axes1[0], palette=palette, errorbar='ci', capsize=0.1)
axes1[0].set_title('A. Mean Total Neighbors Over Time', fontsize=26, fontweight='bold', pad=25)
axes1[0].set_ylabel('Mean Count', fontsize=22, labelpad=20)
axes1[0].set_xlabel('', fontsize=22, labelpad=20)
axes1[0].tick_params(axis='y', labelsize=16)
axes1[0].tick_params(axis='x', labelsize=14) # Keeps X-axis horizontal and readable
axes1[0].legend(title='Timeline', loc='upper left', fontsize=14, title_fontsize=16)

# Bottom: Violin Plot
sns.violinplot(data=melted_total, x='GO_Name', y='Count_Plot', hue='Time_Point', ax=axes1[1], palette=palette, inner='quartile', cut=0, density_norm='width')
axes1[1].set_title('B. Distribution of Total Neighbors (Log Scale)', fontsize=26, fontweight='bold', pad=25)
axes1[1].set_ylabel('Count (Log Scale)', fontsize=22, labelpad=20)
axes1[1].set_xlabel('GO Annotation', fontsize=22, labelpad=20)
axes1[1].set_yscale('log')
axes1[1].tick_params(axis='y', labelsize=16)
axes1[1].tick_params(axis='x', labelsize=14) # Keeps X-axis horizontal and readable
axes1[1].legend(title='Timeline', loc='upper left', fontsize=14, title_fontsize=16)

plt.tight_layout(pad=3.0)
fig1.savefig(outputplot_total, dpi=300, bbox_inches='tight')
plt.close(fig1)

# --- Figure 2 (Annotated Neighbors) ---
fig2, axes2 = plt.subplots(2, 1, figsize=(12, 14))

# Top: Bar Plot
sns.barplot(data=melted_annot, x='GO_Name', y='Count', hue='Time_Point', ax=axes2[0], palette=palette, errorbar='ci', capsize=0.1)
axes2[0].set_title('A. Mean Annotated Neighbors Over Time', fontsize=26, fontweight='bold', pad=25)
axes2[0].set_ylabel('Mean Count', fontsize=22, labelpad=20)
axes2[0].set_xlabel('', fontsize=22, labelpad=20)
axes2[0].tick_params(axis='y', labelsize=16)
axes2[0].tick_params(axis='x', labelsize=14) # Keeps X-axis horizontal and readable
axes2[0].legend(title='Timeline', loc='upper left', fontsize=14, title_fontsize=16)

# Bottom: Violin Plot
sns.violinplot(data=melted_annot, x='GO_Name', y='Count_Plot', hue='Time_Point', ax=axes2[1], palette=palette, inner='quartile', cut=0, density_norm='width')
axes2[1].set_title('B. Distribution of Annotated Neighbors (Log Scale)', fontsize=26, fontweight='bold', pad=25)
axes2[1].set_ylabel('Count (Log Scale)', fontsize=22, labelpad=20)
axes2[1].set_xlabel('GO Annotation', fontsize=22, labelpad=20)
axes2[1].set_yscale('log')
axes2[1].tick_params(axis='y', labelsize=16)
axes2[1].tick_params(axis='x', labelsize=14) # Keeps X-axis horizontal and readable
axes2[1].legend(title='Timeline', loc='upper left', fontsize=14, title_fontsize=16)

plt.tight_layout(pad=3.0)
fig2.savefig(outputplot_annot, dpi=300, bbox_inches='tight')
plt.close(fig2)

print(f"Saving GO {aspect} depth {depth} cutoff {cutoff} networks statistics plots...")
print(f"GO {aspect} depth {depth} cutoff {cutoff} networks statistics plotted!")