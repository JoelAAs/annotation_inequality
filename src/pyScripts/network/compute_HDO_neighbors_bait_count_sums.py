import pickle
import networkx as nx
import pandas as pd

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.bp_network_with_attributes
annotations_df = snakemake.input.annotations
outputdf = snakemake.output.neighbors_bait_count_sums_df
outputpickle = snakemake.output.neighbors_bait_count_sums_pickle

print(f"Processing HDO neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}...\n")

print(f"Loading data for HDO depth {depth} cutoff {cutoff}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
annotations = pd.read_csv(annotations_df, sep = '\t', dtype={'entrez_id': str})
doid_cols = [col for col in annotations.columns if col != 'entrez_id']

print(f"Data for HDO depth {depth} cutoff {cutoff} loaded!\n")

print(f"Computing HDO neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}...")

# Pre-calculating the bait count sums for all nodes, in order not to recompute them each time!
node_neighbor_sums = {
    str(node): sum(G.nodes[neighbor].get('bait_count', 0) for neighbor in G.neighbors(node))
    for node in G.nodes()
}

annotation_to_gene_data = {}

# For each node 
for node, attrs in G.nodes(data = True):
    node_str = str(node)
    doids = attrs.get('doids_with_depth', [])
    neighbor_sum = node_neighbor_sums.get(node_str, 0)

    # For each doid in the node
    for item in doids:
        doid = item['doid']
        if doid not in annotation_to_gene_data:
            annotation_to_gene_data[doid] = []

        annotation_to_gene_data[doid].append({
            'gene': node_str,
            'neighbor_sum': neighbor_sum
        })

results = []

for doid in doid_cols:
    gene_list = annotation_to_gene_data.get(doid, [])
    results.append({
        'annotation': doid,
        'genes_with_sums': gene_list
    })

print(f"Counts for HDO neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff} ready!\n")

print(f"Saving HDO neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}...")

compressed_df = pd.DataFrame(results)
compressed_df = compressed_df.sort_values(by = 'annotation')
compressed_df.to_csv(outputdf, sep = '\t', index = False)

with open(outputpickle, 'wb') as f:
    pickle.dump(compressed_df, f)

print(f"HDO neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}saved!\n")

print(f"HDO neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff} processed!\n")