import pickle
import networkx as nx
import pandas as pd

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.bp_network_with_attributes
annotations_df = snakemake.input.annotations
outputdf = snakemake.output.neighbors_bait_count_sums_df
outputpickle = snakemake.output.neighbors_bait_count_sums_pickle

print(f"Processing GO {aspect} neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}...\n")

print(f"Loading data for GO {aspect} depth {depth} cutoff {cutoff}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
annotations = pd.read_csv(annotations_df, sep = '\t', dtype={'entrez_id': str})
go_id_cols = [col for col in annotations.columns if col != 'entrez_id']

print(f"Data for GO {aspect} depth {depth} cutoff {cutoff} loaded!\n")

print(f"Computing GO {aspect} neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}...")

# Pre-calculating the bait count sums for all nodes, in order not to recompute them each time!
node_neighbor_sums = {
    str(node): sum(G.nodes[neighbor].get('bait_count', 0) for neighbor in G.neighbors(node))
    for node in G.nodes()
}

annotation_to_gene_data = {}

# For each node 
for node, attrs in G.nodes(data = True):
    node_str = str(node)
    goids = attrs.get('go_ids_with_depth', [])
    neighbor_sum = node_neighbor_sums.get(node_str, 0)

    # For each doid in the node
    for item in goids:
        goid = item['go_id']
        if goid not in annotation_to_gene_data:
            annotation_to_gene_data[goid] = []

        annotation_to_gene_data[goid].append({
            'gene': node_str,
            'neighbor_sum': neighbor_sum
        })

results = []

for goid in go_id_cols:
    gene_list = annotation_to_gene_data.get(goid, [])
    results.append({
        'annotation': goid,
        'genes_with_sums': gene_list
    })

print(f"Counts for GO {aspect} neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff} ready!\n")

print(f"Saving GO {aspect} neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}...")

compressed_df = pd.DataFrame(results)
compressed_df = compressed_df.sort_values(by = 'annotation')
compressed_df.to_csv(outputdf, sep = '\t', index = False)

with open(outputpickle, 'wb') as f:
    pickle.dump(compressed_df, f)

print(f"GO {aspect} neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff}saved!\n")

print(f"GO {aspect} neighbors bait count sums for each annotation at depth {depth} cutoff {cutoff} processed!\n")