import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

inputfile = snakemake.input.correlation_values
n_of_baits_df = snakemake.input.n_of_baits
outputfile = snakemake.output.correlation_matrix
method = snakemake.wildcards.method

method = method.capitalize()
n_of_baits = pd.read_csv(n_of_baits_df, sep = '\t')
corr_df = pd.read_csv(inputfile, sep = '\t', index_col = 0)
corr_df = corr_df.rename(columns = {'count_studies': 'Bait count',
                                    'count_annot_hdo': 'HDO',
                                    'count_annot_disgenet': 'DISGENET',
                                    'count_annot_bp': 'GO:BP',
                                    'count_annot_mf': 'GO:MF',
                                    'count_annot_cc': 'GO:CC'})
corr_df.index = corr_df.columns

baits = n_of_baits.iloc[0]['n_of_baits']

colors = ['#4D7902', '#F9A602']
custom_cmap = LinearSegmentedColormap.from_list('green_orange', colors)

plt.figure(figsize = (12, 10))
plot = sns.heatmap(
    corr_df.iloc[::-1],
    annot = True,
    cmap = custom_cmap,
    square = True,
    fmt = '.2f',
    annot_kws={"size": 16, "weight": "bold"}, 
    cbar_kws={'label': f'{method} Correlation Coefficient', 'shrink': 0.8, 'aspect': 20}
)

plt.xticks(rotation = 0, fontsize = 16)
plt.yticks(rotation = 0, fontsize = 16)
plt.title(f'{method} Correlation Matrix | N = {baits} Baits', fontsize = 26, fontweight='bold', pad = 25)
cbar = plot.collections[0].colorbar
cbar.ax.tick_params(labelsize=16)
cbar.set_label(f'{method} Correlation Value', fontsize=18, labelpad=15)
plt.tight_layout()
plt.savefig(outputfile, dpi = 300, bbox_inches = 'tight')