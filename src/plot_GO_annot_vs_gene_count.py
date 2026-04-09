import pandas as pd
import matplotlib.pyplot as plt

aspect = snakemake.wildcards.aspect
input_df = snakemake.input.annotation_df
bait_usage_df = snakemake.input.bait_usage
outputdf = snakemake.output.count_df
outputfile = snakemake.output.output_plot

print(f'Loading GO {aspect} data...')

df = pd.read_csv(input_df, sep = '\t')
df = df.dropna()

bait_usage = pd.read_csv(bait_usage_df, sep = '\t')

bait_ids = bait_usage['entrez_id_bait'].unique()
df = df[df['entrez_id'].isin(bait_ids)]

print(f'GO {aspect} data loaded!')

print(f'Obtaining GO {aspect} annotations gene counts...')

counts = df['go_id'].value_counts().reset_index()
counts.columns = ['go_id', 'gene_count']

counts.to_csv(outputdf, sep = '\t', index = False)

print(f'GO {aspect} annotations gene counts obtained\n')

print(f'Plotting GO {aspect} distribution of GO ids sizes...')

plt.figure(figsize = (10, 6))
plt.bar(counts['gene_count'], counts['go_id'])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Gene Count (log-scale)')
plt.ylabel('GO ids (log-scale)')
plt.title(f'Distribution of GO {aspect} Annotation Sizes')
plt.tight_layout()
plt.savefig(outputfile, dpi = 300)
plt.close()

print(f'GO {aspect} distribution of annotation sizes plotted!')