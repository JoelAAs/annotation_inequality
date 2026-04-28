import pandas as pd
import pickle
import networkx as nx

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.network_with_all_dates
coeff_df = snakemake.input.coefficients
outputdf = snakemake.output.nodes_with_top_5_annotations_df
outputpickle = snakemake.output.nodes_with_top_5_annotations_pickle

print(f"Computing lists of annotated GO {aspect} genes with the top 5 annotations for depth {depth} and cutoff {cutoff}...")

print(f"Loading data for GO {aspect} depth {depth} cutoff {cutoff}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
coefficients = pd.read_csv(coeff_df, sep = '\t')

print(f"Processing GO {aspect} depth {depth} cutoff {cutoff} lists...")

coefficients['abs_coeff'] = coefficients['Coefficient'].abs()
top_5_go = coefficients.nlargest(5, 'abs_coeff')['GO_id'].tolist()

grouped_genes = {
    go_id: [] for go_id in top_5_go
}

for n, attr in G.nodes(data = True):
    annotations = attr.get('go_annotations', [])

    for d in annotations:
        go_id = d.get('go_id')
        # If it's one of our top 5, append the gene to that specific list!
        if go_id in grouped_genes:
            grouped_genes[go_id].append(n)

final_df = pd.DataFrame(list(grouped_genes.items()), columns = ['GO_id', 'annotated_genes'])

print(f"Saving lists of annotated genes with top 5 GO {aspect} depth {depth} cutoff {cutoff} annotations...")

final_df.to_csv(outputdf, sep  = '\t', index = False)

with open(outputpickle, 'wb') as f:
    pickle.dump(final_df, f)

print(f"Lists of annotated GO {aspect} genes with the top 5 annotations for depth {depth} and cutoff {cutoff} ready!")