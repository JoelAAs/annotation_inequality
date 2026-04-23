import pandas as pd
import matplotlib.pyplot as plt

cutoff = snakemake.wildcards.cutoff
input_df = snakemake.input.full_adj_r2_file
outputplot = snakemake.output.adj_r2_plot

print(f"Plotting HDO cutoff {cutoff} R^2 per depth...\n")

print(f"Loading HDO cutoff {cutoff} data...")

df = pd.read_csv(input_df, sep = ':')

x = df['depth']
y = df['adjusted_r2']

print(f"HDO cutoff {cutoff} data loaded!\n")

plt.figure(figsize = (12, 8))
plt.bar(x, y, color = 'skyblue', edgecolor = 'black')
plt.grid(axis='y', linestyle='--', alpha=0.8)
plt.xticks(df['depth'])
plt.xlabel('Depth', fontsize = 22, labelpad = 20)
plt.ylabel('Adjusted R^2', fontsize = 22, labelpad = 20)
plt.title(f'Adjusted HDO R^2 value per depth, cutoff = {cutoff}', fontsize = 26, fontweight = 'bold', pad = 25)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.margins(y=0)
plt.ylim(0, y.max() * 1.05)
plt.tight_layout()
plt.savefig(outputplot, dpi = 300, bbox_inches='tight')
plt.close()

print(f"HDO cutoff {cutoff} R^2 per depth plot ready!\n")