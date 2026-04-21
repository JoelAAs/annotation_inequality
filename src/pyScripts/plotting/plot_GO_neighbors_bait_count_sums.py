import pickle
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
neighbors_sums_df = snakemake.input.neighbors_bait_count_sums_pickle
outputpdf = snakemake.output.plots_pdf

print(f"Generating multi-page PDF for GO {aspect} depth {depth} cutoff {cutoff} neighbors bait count sums...\n")

print(f"Loading data for GO {aspect} depth {depth} cutoff {cutoff} neighbors bait count sums...")

with open(neighbors_sums_df, 'rb') as f:
    df = pickle.load(f)

print(f"GO {aspect} depth {depth} cutoff {cutoff} neighbors bait count sums data loaded!\n")

print(f"Plotting GO {aspect} depth {depth} cutoff {cutoff} neighbors bait count sums...")

sns.set_style('whitegrid')

with PdfPages(outputpdf) as pdf:
    d = pdf.infodict()
    d['Title'] = 'Neighbor Bait Count Distributions'
    
    for index, row in df.iterrows():
        # Each time select the annotation and the list of genes' neighbors bait count sums
        goid = row['annotation']
        gene_data = row['genes_with_sums']

        # Extract the sum values
        sums = [item['neighbor_sum'] for item in gene_data]

        # Skip annotations if no genes are found in the network for them
        if not sums:
            continue

        # Create a new figure for each page:
        plt.figure(figsize = (10, 7))

        # Plot and format
        sns.histplot(sums, kde = (len(sums) > 1), color = 'teal', bins = 20)
        plt.title(f'Frequency of GO {aspect} Neighbor Bait Count Sums\nAnnotation: {goid} | Genes with annotation: {len(sums)}', fontsize = 14)
        plt.xlabel('Bait Count Sum Of Neighbors (log-scale)', fontsize = 12)
        plt.ylabel('Frequency', fontsize = 12)
        plt.xscale('log')
        plt.tight_layout()
        pdf.savefig()
        plt.close()

print(f"GO {aspect} depth {depth} cutoff {cutoff} neighbors bait count sums plotted!\n")

print(f"Multi-page PDF for GO {aspect} depth {depth} cutoff {cutoff} neighbors bait count sums ready!\n")