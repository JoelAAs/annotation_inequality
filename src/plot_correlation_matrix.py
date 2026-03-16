import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

inputfile = snakemake.input.correlation_values
outputfile = snakemake.output.correlation_matrix
method = snakemake.wildcards.method

corr_df = pd.read_csv(inputfile, sep = '\t', index_col = 0)
corr_df = corr_df.rename(columns = {'count_studies': 'Bait count',
                                    'count_annot_hdo': 'HDO',
                                    'count_annot_disgenet': 'DISGENET',
                                    'count_annot_bp': 'GO:BP',
                                    'count_annot_mf': 'GO:MF',
                                    'count_annot_cc': 'GO:CC'})
corr_df.index = corr_df.columns

colors = ['#4D7902', '#F9A602']
custom_cmap = LinearSegmentedColormap.from_list('green_orange', colors)

plt.figure(figsize = (10, 6))
sns.heatmap(corr_df.iloc[::-1],
            annot = True,
            cmap = custom_cmap,
            square = True,
            fmt = '.2g')

plt.xticks(rotation = 0)
plt.yticks(rotation = 0)
plt.title(f'{method} correlation matrix')
plt.tight_layout()
plt.savefig(outputfile, dpi = 300)