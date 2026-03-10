import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

def get_spearman_correlation(df):
    x = df['count']
    y = df['annotation_count']

    rho, pval = spearmanr(x, y)
    print(f'rho = {rho}\tp-value = {pval}')

    return rho, pval

def make_scatterplot():
    inputs = snakemake.input
    outputs = snakemake.output

    for inputfile, outputfile in zip(inputs, outputs):
        df = pd.read_parquet(inputfile)
        print(f'Opening file {inputfile}')
        print(f'Svaing image on {outputfile}')

        rho, pval = get_spearman_correlation(df)

        plt.figure(figsize = (10, 8))
        plt.scatter(df['count'], df['annotation_count'], alpha = 0.6, s = 50)
        plt.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\np-value = {pval:.3f}',
                 transform = plt.gca().transAxes, fontsize = 12,
                 bbox = dict(boxstyle = "round,pad=0.3", facecolor = "white"))
        plt.xlabel('Interactions_count')
        plt.ylabel('Annotations_count')
        plt.grid(True, alpha = 0.3)
        plt.title(f'Spearman Correlation\nρ = {rho:.3f}')

        plt.savefig(outputfile, dpi = 300, bbox_inches = 'tight')
        plt.close()

make_scatterplot()