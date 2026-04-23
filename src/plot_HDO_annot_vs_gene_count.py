import pandas as pd
import matplotlib.pyplot as plt

input_df = snakemake.input.annotation_df
bait_usage_df = snakemake.input.bait_usage
outputdf = snakemake.output.count_df
outputfile = snakemake.output.output_plot

df = pd.read_csv(input_df, sep = '\t')
df = df.dropna()

bait_usage = pd.read_csv(bait_usage_df, sep = '\t')

bait_ids = bait_usage['entrez_id_bait'].unique()
df = df[df['entrez_id'].isin(bait_ids)]

print('Obtaining HDO annotations gene counts...')

counts = df['doid'].value_counts().reset_index()
counts.columns = ['doid', 'gene_count']

counts.to_csv(outputdf, sep = '\t', index = False)

print('HDO annotations gene counts obtained\n')

print('Plotting HDO distribution of doid sizes...')

plt.figure(figsize = (12, 8))
plt.bar(counts['gene_count'], counts['doid'])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Gene Count (log-scale)', fontsize = 22, labelpad = 15)
plt.ylabel('DOIDs (log-scale)', fontsize = 22, labelpad = 15)
plt.title('Distribution of Annotation Sizes', fontsize = 28, fontweight = 'bold', pad = 25)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(outputfile, dpi = 300, bbox_inches='tight')
plt.close()

print('HDO distribution of annotation sizes plotted!')