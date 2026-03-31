import pandas as pd
import matplotlib.pyplot as plt

input_df = snakemake.input.annotation_df
outputdf = snakemake.output.count_df
outputfile = snakemake.output.output_plot

df = pd.read_csv(input_df, sep = '\t')
df = df.dropna()

print('Obtaining HDO annotations gene counts...')

counts = df['doid'].value_counts().reset_index()
counts.columns = ['doid', 'gene_count']

counts.to_csv(outputdf, sep = '\t', index = False)

print('HDO annotations gene counts obtained\n')

print('Plotting HDO distribution of doid sizes...')

plt.figure(figsize = (10, 6))
plt.bar(counts['gene_count'], counts['doid'])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Gene Count (log-scale)')
plt.ylabel('DOIDs (log-scale)')
plt.title('Distribution of Annotation Sizes')
plt.tight_layout()
plt.savefig(outputfile, dpi = 300)
plt.close()

print('HDO distribution of annotation sizes plotted!')