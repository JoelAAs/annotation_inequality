import pandas as pd
import pickle
import networkx as nx

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.network_with_all_dates
coeff_df = snakemake.input.coefficients
outputdf = snakemake.output.nodes_with_top_5_annotations_df
outputpickle = snakemake.output.nodes_with_top_5_annotations_pickle

print(f"Computing lists of annotated HDO genes with the top 5 annotations for depth {depth} and cutoff {cutoff}...")

print(f"Loading data for HDO depth {depth} cutoff {cutoff}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
coefficients = pd.read_csv(coeff_df, sep = '\t')

print(f"Processing HDO depth {depth} cutoff {cutoff} lists...")

coefficients['abs_coeff'] = coefficients['Coefficient'].abs()
top_5_do = coefficients.nlargest(5, 'abs_coeff')['HDO_doid'].tolist()

grouped_genes = {
    do_id: [] for do_id in top_5_do
}

for n, attr in G.nodes(data = True):
    annotations = attr.get('do_annotations', [])

    for d in annotations:
        do_id = d.get('doid')
        # If it's one of our top 5, append the gene to that specific list!
        if do_id in grouped_genes:
            grouped_genes[do_id].append(n)

final_df = pd.DataFrame(list(grouped_genes.items()), columns = ['DO_id', 'annotated_genes'])

print(f"Saving lists of annotated genes with top 5 HDO depth {depth} cutoff {cutoff} annotations...")

final_df.to_csv(outputdf, sep  = '\t', index = False)

with open(outputpickle, 'wb') as f:
    pickle.dump(final_df, f)

print(f"Lists of annotated HDO genes with the top 5 annotations for depth {depth} and cutoff {cutoff} ready!")