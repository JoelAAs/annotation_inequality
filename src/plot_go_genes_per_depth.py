import pandas as pd
import matplotlib.pyplot as plt

aspect = snakemake.wildcards.aspect
input_textfile = snakemake.input.genes_per_depth
outputplot = snakemake.output.genes_per_depth_plot

print(f'Plotting GO {aspect} number of genes per depth...')

df_textfile = pd.read_csv(input_textfile, sep = ':')

x = df_textfile['depth']
y = df_textfile['n_of_genes']

plt.figure(figsize = (10, 6))
plt.bar(x, y, color = 'skyblue', edgecolor = 'black')
plt.yscale('log')
plt.xlabel('Depth')
plt.ylabel('Number of Genes (log-scale)')
plt.title(f'Number of GO {aspect} Genes per Depth level')
plt.savefig(outputplot, dpi = 300)
plt.close()

print(f'GO {aspect} number of genes per depth plotted!')