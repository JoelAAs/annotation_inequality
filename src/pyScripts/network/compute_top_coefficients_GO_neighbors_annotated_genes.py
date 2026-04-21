import pandas as pd
import pickle
import networkx as nx

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
inputnetwork = snakemake.input.bp_network_with_attributes
coefficients_df = snakemake.input.coefficients
outputdf = snakemake.output.neighbors_with_annotation_df
outputpickle = snakemake.output.neighbors_with_annotation_pickle

print(f"Computing annotated neighbors for GO {aspect} top coefficients depth {depth} cutoff {cutoff}...")

print(f"Loading data for GO {aspect} top coefficients depth {depth} cutoff {cutoff}...")

with open(inputnetwork, 'rb') as f:
    G = pickle.load(f)
df_coeffs = pd.read_csv(coefficients_df, sep = '\t')

print(f"GO {aspect} top coefficients depth {depth} cutoff {cutoff} data loaded!")

print(f"Processing GO {aspect} top coefficients depth {depth} cutoff {cutoff} annotated neighbors...")

df_coeffs['abs_coeff'] = df_coeffs['Coefficient'].abs()
top_5_go = df_coeffs.nlargest(5, 'abs_coeff')['GO_id'].tolist()
    
results = []

for go_id in top_5_go:
    # Find all genes with the specific GO id
    annotated_genes = [
        n for n, attr in G.nodes(data = True)
        if 'go_ids_with_depth' in attr and any(d['go_id'] == go_id for d in attr['go_ids_with_depth'])
    ]

    genes_with_annotated_neighbors = []

    for gene in annotated_genes:
        # Look at all neighbors and count how many have the same GO id we are looking at
        annotated_neighbors_count = 0

        for neighbor in G.neighbors(gene):
            neighbor_attr = G.nodes[neighbor]
            if 'go_ids_with_depth' in neighbor_attr:
                # Check if the neighbor has the speific GO id
                if any(d['go_id'] == go_id for d in neighbor_attr['go_ids_with_depth']):
                    annotated_neighbors_count += 1

        genes_with_annotated_neighbors.append({
            'gene': str(gene),
            'annotated_neighbors_count': int(annotated_neighbors_count)
        })

    results.append({
        'annotation': go_id,
        'genes_with_annotated_neighbors': genes_with_annotated_neighbors
    })

print(f"GO {aspect} top coefficients depth {depth} cutoff {cutoff} annotated neighbors processed!")

print(f"Saving GO {aspect} top coefficients depth {depth} cutoff {cutoff} annotated neighbors...")

final_df = pd.DataFrame(results)
final_df.to_csv(outputdf, sep = '\t', index = False)

with open(outputpickle, 'wb') as f:
    pickle.dump(final_df, f)

print(f"GO {aspect} top coefficients depth {depth} cutoff {cutoff} annotated neighbors saved!")

print(f"Annotated neighbors for GO {aspect} top coefficients depth {depth} cutoff {cutoff} ready!")