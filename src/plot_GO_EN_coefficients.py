import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
coefficients = snakemake.input.elastic_net_coefficients
top_plot = snakemake.output.top_plot
distrib_plot = snakemake.output.distribution_plot
    
print(f'Plotting plots for GO {aspect} depth {depth} coefficients...')

df = pd.read_csv(coefficients, sep = '\t')

# Plot the top 20 positive and negative coefficients
top_n = 20
df_pos = df[df['Coefficient'] > 0].sort_values(by = 'Coefficient', ascending = False).head(top_n)
df_neg = df[df['Coefficient'] < 0].sort_values(by = 'Coefficient').head(top_n)
df_plot = pd.concat([df_neg, df_pos])

df_plot['color'] = df_plot['Coefficient'].apply(lambda x: 'green' if x > 0 else 'red')

go_aspect_term = f'GO_{aspect}_term'
plt.figure(figsize = (16 ,8))
plt.barh(df_plot['GO_term'], df_plot['Coefficient'], color = df_plot['color'])
plt.xlabel('Coefficient')
plt.title(f'Top depth {depth} GO {aspect} Elastic Net Coefficients (Green = Positive, Red = Negative)')
plt.gca().invert_yaxis()
plt.xticks(rotation = 0)
plt.yticks(rotation = 0)
plt.tight_layout()
plt.savefig(top_plot, dpi = 300)
plt.close()

# Plot the coefficients distribution with a histogram
coeff_values = df['Coefficient'].round(3)
counts = Counter(coeff_values)

values = list(counts.keys())
freqs = list(counts.values())
total_coefs = len(df)

plt.figure(figsize = (10, 6))
plt.hist(df['Coefficient'], bins = 30, color = 'skyblue', edgecolor = 'black')
plt.yscale('log')
plt.xlabel('Coefficient value')
plt.ylabel('Frequency (log-scale)')
plt.title(f'Histogram of depth {depth} GO {aspect} Elastic Net Coefficients')
plt.legend([f'Total coefficients = {total_coefs}'], loc='upper right')
plt.savefig(distrib_plot, dpi = 300)
plt.close()

print(f'GO {aspect} depth {depth} plots done!')