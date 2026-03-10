import pandas as pd
import matplotlib.pyplot as plt

def make_scatterplot():
    inputfile = snakemake.input
    outputfile = snakemake.output

    df = pd.read_parquet(inputfile)

    plt.figure(figsize = (10, 6))
    plt.scatter(df[x_col], df[y_col], alpha = 0.6, s = 50)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.grid(True, alpha = 0.3)
    plt.title(f'{x_col} vs {y_col}')

    plt.savefig(outputfile, dpi = 300, bbox_inches = 'tight')
    plt.close()

make_scatterplot()