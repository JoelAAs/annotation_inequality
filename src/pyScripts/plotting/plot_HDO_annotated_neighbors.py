import pickle
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
annotated_neighbors_df = snakemake.input.neighbors_with_annotation_pickle
outputpdf = snakemake.output.plots_pdf

print(f"Generating multi-page PDF for HDO depth {depth} cutoff {cutoff} annotated neighbors...\n")

print(f"Loading data for HDO depth {depth} cutoff {cutoff} annotated neighbors...")

with open(annotated_neighbors_df, 'rb') as f:
    df = pickle.load(f)

print(f"HDO depth {depth} cutoff {cutoff} annotated neighbors data loaded!\n")

print(f"Plotting HDO depth {depth} cutoff {cutoff} annotated neighbors...")

sns.set_style('whitegrid')

with PdfPages(outputpdf) as pdf:
    d = pdf.infodict()
    d['Title'] = 'Annotated Neighbors Distributions'
    
    for index, row in df.iterrows():
        # Each time select the annotation and the list of genes' annotated neighbors
        doid = row['annotation']
        gene_data = row['genes_with_annotated_neighbors']

        # Extract the sum values
        sums = [item['annotated_neighbors_count'] for item in gene_data]

        # Skip annotations if no genes are found in the network for them
        if not sums:
            continue

        # Create a new figure for each page:
        plt.figure(figsize = (10, 7))

        # Plot and format
        sns.histplot(sums, kde = (len(sums) > 1), color = 'teal', bins = 20)
        plt.title(f'Frequency of HDO Annotated Neighbors\nAnnotation: {doid} | Genes with annotation: {len(sums)}', fontsize = 14)
        plt.xlabel('Annotated Neighbors (log-scale)', fontsize = 12)
        plt.ylabel('Frequency', fontsize = 12)
        plt.xscale('log')
        plt.tight_layout()
        pdf.savefig()
        plt.close()

print(f"HDO depth {depth} cutoff {cutoff} annotated neighbors plotted!\n")

print(f"Multi-page PDF for HDO depth {depth} cutoff {cutoff} annotated neighbors ready!\n")