import pandas as pd
import matplotlib.pyplot as plt

aspect = snakemake.wildcards.aspect
input_textfile = snakemake.input.annotations_per_depth
outputplot = snakemake.output.annotations_per_depth_plot

print(f'Plotting GO {aspect} number of coefficients per depth...')

df_textfile = pd.read_csv(input_textfile, sep = ':')

x = df_textfile['depth']
y = df_textfile['n_of_coefficients']

plt.figure(figsize = (10, 6))
plt.bar(x, y, color = 'skyblue', edgecolor = 'black')
plt.yscale('log')
plt.xlabel('Depth')
plt.ylabel('Number of Coefficients (log-scale)')
plt.title(f'Number of GO {aspect} Coefficients per Depth level')
plt.savefig(outputplot, dpi = 300)
plt.close()

print(f'GO {aspect} number of coefficients per depth plotted!')