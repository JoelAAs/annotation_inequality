import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pathlib import Path

def get_spearman_correlation(df):
    '''
    col_1 = df['annotation_count']
    col_2 = df['count']
    rho_test = col_1.corr(col_2, method = 'spearman')
    print(f'Spearman correlation with pandas: {rho_test}')
    '''

    df_clean = df.dropna(subset = ['count', 'annotation_count'])

    x = df_clean['annotation_count']
    y = df_clean['count']

    rho, pval = spearmanr(x, y)
    print(f'rho = {rho}\tp-value = {pval}')

    return rho, pval

def get_name(input_path):
    return Path(input_path).stem.split('_')[-1]

def make_scatterplot():
    inputs = snakemake.input
    outputs = snakemake.output

    for inputfile, outputfile in zip(inputs, outputs):
        df = pd.read_parquet(inputfile)
        print(f'Opening file {inputfile}')
        print(f'Saving image on {outputfile}')

        rho, pval = get_spearman_correlation(df)
        bait_or_prey = get_name(inputfile)

        plt.figure(figsize = (10, 8))
        plt.scatter(df['annotation_count'], df['count'], alpha = 0.6, s = 50)
        plt.text(0.05, 0.95, f'Spearman ρ = {rho:.5f}\np-value = {pval:.3e}',
                 transform = plt.gca().transAxes, fontsize = 12,
                 bbox = dict(boxstyle = "round,pad=0.3", facecolor = "white"))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Annotations Count (log scale)')
        plt.ylabel('Interactions Count (log scale)')
        plt.grid(True, alpha = 0.3)
        plt.title(f'Spearman Correlation Annotations vs Interactions counts ({bait_or_prey})\nρ = {rho:.3f}')

        plt.savefig(outputfile, dpi = 300, bbox_inches = 'tight')
        plt.close()

make_scatterplot()