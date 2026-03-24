import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

input_coeff = snakemake.input.elastic_net_coefficients
input_textfile = snakemake.input.annotations_per_depth
output_top = snakemake.output.top_coefficients
output_distrib = snakemake.output.coefficients_distribution
output_plot = snakemake.output.annotations_per_depth_plot

for coefficients, top_coeff, distrib in zip(input_coeff, output_top, output_distrib):
    df = pd.read_csv(coefficients, sep = '\t')
    depth = coefficients.split("_")[-1].replace(".csv", "")
    
    print(f'Plotting plots for depth {depth} HDO coefficients...')

    # Plot the top 20 positive and negative coefficients
    top_n = 20
    df_pos = df[df['Coefficient'] > 0].sort_values(by = 'Coefficient', ascending = False).head(top_n)
    df_neg = df[df['Coefficient'] < 0].sort_values(by = 'Coefficient').head(top_n)
    df_plot = pd.concat([df_neg, df_pos])

    df_plot['color'] = df_plot['Coefficient'].apply(lambda x: 'green' if x > 0 else 'red')

    plt.figure(figsize = (8 ,6))
    plt.barh(df_plot['HDO_Term'], df_plot['Coefficient'], color = df_plot['color'])
    plt.xlabel('Coefficient')
    plt.title(f'Top depth {depth} HDO Elastic Net Coefficients (Green = Positive, Red = Negative)')
    plt.gca().invert_yaxis()
    plt.xticks(rotation = 0)
    plt.yticks(rotation = 0)
    plt.savefig(top_coeff, dpi = 300)
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
    plt.title(f'Histogram of depth {depth} HDO Elastic Net Coefficients')
    plt.legend([f'Total coefficients = {total_coefs}'], loc='upper right')
    plt.savefig(distrib, dpi = 300)
    plt.close()
    
    print(f'Plots for depth {depth} done!')
    
print(f'Plotting HDO coefficient counts per depth...')

df_textfile = pd.read_csv(input_textfile, sep = ':')

x = df_textfile['depth']
y = df_textfile['n_of_coefficients']

plt.figure(figsize = (8, 6))
plt.bar(x, y, color = 'skyblue', edgecolor = 'black')
plt.yscale('log')
plt.xlabel('Depth')
plt.ylabel('Number of Coefficients (log-scale)')
plt.title('Number of Coefficients per Depth Level')
plt.xticks(x)
plt.savefig(output_plot, dpi = 300)
plt.close()

print(f'Coefficient counts per depth plotted!')