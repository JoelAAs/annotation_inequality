import pandas as pd
import matplotlib.pyplot as plt

degree_frequencies_df = snakemake.input.degree_frequencies
outputplot = snakemake.output.degree_frequencies_plot

print(f"Processing raw network degree frequencies plot...\n")

print(f"Loading data...")

df = pd.read_csv(degree_frequencies_df, sep = '\t')

print(f"Data loaded!\n")

print(f"Plotting raw network degree frequencies...")

total_nodes = df['entrez_id'].nunique()
df_distrib = df['degree'].value_counts().reset_index()
df_distrib.columns = ['degree', 'frequency']
df_distrib = df_distrib.sort_values('degree')
max_frequency = df_distrib['frequency'].max()
max_degree = df_distrib['degree'].max()

x = df_distrib['degree']
y = df_distrib['frequency']

plt.figure(figsize = (10,6))
plt.bar(x, y, color = 'skyblue', edgecolor = 'black', zorder = 0, label = total_nodes)
plt.xscale('log')
plt.yscale('log')
plt.xlim(0, max_degree * 1.1)
plt.ylim(0, max_frequency * 1.1)
plt.xlabel('Degree (log-scale)')
plt.ylabel('Frequency (log-scale)')
plt.title("Bait-Prey Publications Network's Node Degree Distribution")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc = 'upper right', frameon = False, fontsize = 'small', title = 'Total nodes')
plt.tight_layout()

print(f"Raw network degree frequencies plotting done!\n")

print(f"Saving raw network degree frequencies plot...")

plt.savefig(outputplot)
plt.close()

print(f"Raw network degree frequencies plot saved!\n")

print(f"Raw network degree frequencies plot ready!\n")