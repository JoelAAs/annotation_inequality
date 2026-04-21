import pandas as pd
import pickle
import networkx as nx

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
inputnetwork = snakemake.input.bp_network_with_attributes
coefficients_df = snakemake.input.coefficients
outputdf = snakemake.output.neighbors_with_annotation_df
outputpickle = snakemake.output.neighbors_with_annotation_pickle

print(f"Computing annotated neighbors for HDO top coefficients depth {depth} cutoff {cutoff}...")

print(f"Loading data for HDO top coefficients depth {depth} cutoff {cutoff}...")

with open(inputnetwork, 'rb') as f:
    G = pickle.load(f)
df_coeffs = pd.read_csv(coefficients_df, sep = '\t')

print(f"HDO top coefficients depth {depth} cutoff {cutoff} data loaded!")

print(f"Processing HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors...")

df_coeffs['abs_coeff'] = df_coeffs['Coefficient'].abs()
top_5_hdo = df_coeffs.nlargest(5, 'abs_coeff')['HDO_doid'].tolist()
    
results = []

for hdo_id in top_5_hdo:
    # Find all genes with the specific HDO id
    annotated_genes = [
        n for n, attr in G.nodes(data = True)
        if 'doids_with_depth' in attr and any(d['doid'] == hdo_id for d in attr['doids_with_depth'])
    ]

    genes_with_annotated_neighbors = []

    for gene in annotated_genes:
        # Look at all neighbors and count how many have the same HDO id we are looking at
        annotated_neighbors_count = 0

        for neighbor in G.neighbors(gene):
            neighbor_attr = G.nodes[neighbor]
            if 'doids_with_depth' in neighbor_attr:
                # Check if the neighbor has the speific HDO id
                if any(d['doid'] == hdo_id for d in neighbor_attr['doids_with_depth']):
                    annotated_neighbors_count += 1

        genes_with_annotated_neighbors.append({
            'gene': str(gene),
            'annotated_neighbors_count': int(annotated_neighbors_count)
        })

    results.append({
        'annotation': hdo_id,
        'genes_with_annotated_neighbors': genes_with_annotated_neighbors
    })

print(f"HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors processed!")

print(f"Saving HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors...")

final_df = pd.DataFrame(results)
final_df.to_csv(outputdf, sep = '\t', index = False)

with open(outputpickle, 'wb') as f:
    pickle.dump(final_df, f)

print(f"HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors saved!")

print(f"Annotated neighbors for HDO top coefficients depth {depth} cutoff {cutoff} ready!")