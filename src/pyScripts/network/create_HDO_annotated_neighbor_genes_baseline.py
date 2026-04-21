import pandas as pd
import pickle
import networkx as nx
import random

depth = snakemake.wildcards. depth
cutoff = snakemake.wildcards.cutoff
inputnetwork = snakemake.input.bp_network_with_attributes
coefficients_df = snakemake.input.coefficients
baseline_df = snakemake.output.baseline_neighbors_with_annotation_df
baseline_pickle = snakemake.output.baseline_neighbors_with_annotation_pickle

print(f"Building annotated neighbors baseline for HDO top coefficients depth {depth} cutoff {cutoff}...")

print(f"Loading data for HDO top coefficients depth {depth} cutoff {cutoff}...")

with open(inputnetwork, 'rb') as f:
    G = pickle.load(f)
coefficients = pd.read_csv(coefficients_df, sep='\t')

print(f"HDO top coefficients depth {depth} cutoff {cutoff} data loaded!")

print(f"Processing HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors baseline...")

# Identify top 5 HDO IDs
coefficients['abs_coeff'] = coefficients['Coefficient'].abs()
top_5_hdo = coefficients.nlargest(5, 'abs_coeff')['HDO_doid'].tolist()

# Pre-compute degree-to-nodes 
degree_map = {}
for node, degree in G.degree():
    if degree not in degree_map:
        degree_map[degree] = []
    degree_map[degree].append(node)

results = []

# Generate baseline
for doid in top_5_hdo:
    # Find the real genes annotated with this specific DOID
    annotated_genes = [
        n for n, attr in G.nodes(data=True)
        if 'doids_with_depth' in attr and any(d['doid'] == doid for d in attr['doids_with_depth'])
    ]

    baseline_genes_with_counts = []

    for real_gene in annotated_genes:
        # Get degree of the real gene
        real_degree = G.degree(real_gene)

        # Select a random node with the same degree
        random_node = random.choice(degree_map[real_degree])

        # Count neighbors of the random node that have the specific annotation
        annotated_neighbors_count = 0
        for neighbor in G.neighbors(random_node):
            neighbor_attr = G.nodes[neighbor]
            if 'doids_with_depth' in neighbor_attr:
                if any(d['doid'] == doid for d in neighbor_attr['doids_with_depth']):
                    annotated_neighbors_count += 1
        
        baseline_genes_with_counts.append({
            'gene': str(random_node),
            'annotated_neighbors_count': int(annotated_neighbors_count)
        })

    results.append({
        'annotation': doid,
        'genes_with_annotated_neighbors': baseline_genes_with_counts
    })

print(f"HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors baseline processed!")

print(f"Saving HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors baseline...")

final_df = pd.DataFrame(results)
final_df.to_csv(baseline_df, sep='\t', index=False)

with open(baseline_pickle, 'wb') as f:
    pickle.dump(final_df, f)

print(f"HDO top coefficients depth {depth} cutoff {cutoff} annotated neighbors baseline saved!")

print(f"Annotated neighbors baseline for HDO top coefficients depth {depth} cutoff {cutoff} ready!")
