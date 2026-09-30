import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

matrix = snakemake.input.confusion_matrix
heatmap = snakemake.output.heatmap

X_matrix = pd.read_parquet(matrix)

# Calculate the correlation matrix
# For binary data (0s and 1s), Pearson correlation is mathematically equivalent to the Phi coefficient
covariation_matrix = X_matrix.corr()

# Create a clustered heatmap
plt.figure(figsize=(10, 10))
clustermap = sns.clustermap(
    covariation_matrix, 
    annot=True,
    cmap="vlag",
    center=0,
    vmin=-1,
    vmax=1,
    fmt=".2f",
    figsize=(10, 8),
    cbar_kws={'label': 'Correlation (Phi coefficient)'}
)

clustermap.fig.suptitle("Covariation of GO and HDO Terms\n(Based on shared gene annotations)", y=1.02, fontsize=14)

clustermap.savefig(heatmap, dpi=300, bbox_inches='tight')

print(f"Heatmap salvata con successo in {heatmap}", flush=True)