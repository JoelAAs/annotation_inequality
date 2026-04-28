import pandas as pd
import networkx as nx
import pickle

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.final_network
top_annot_df = snakemake.input.nodes_with_top_5_annotations_pickle
outputstatistics = snakemake.output.networks_statistics

print(f"Computing across-time networks statistics for GO {aspect} depth {depth} cutoff {cutoff}...")

print(f"Loading data for GO {aspect} depth {depth} cutoff {cutoff}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
with open(top_annot_df, 'rb') as f:
    top_annot = pickle.load(f)

print(f"Processing GO {aspect} depth {depth} cutoff {cutoff} temporal neighbor statistics...")

# HELPER FUNCTION
def get_annotation_date(node_id, target_go):
    """Returns the full YYYYMMDD integer for a given node and GO term."""
    annotations = G.nodes[node_id].get('go_annotations', [])
    for ann in annotations:
        if ann.get('go_id') == target_go:
            date = ann.get('first_annotation_date')
            if date is not None:
                return int(date) 
    return None

statistics_data = []

# For every Top 5 GO term
for index, row in top_annot.iterrows():
    go_id = row['GO_id']
    genes = row['annotated_genes']

    # For every gene annotated with the specific annotation
    for target_node in genes:
        target_date = get_annotation_date(target_node, go_id)
        if target_date is None:
            continue

        # Annotated neighbors counters
        annot_count_past = 0
        annot_count_present = 0
        annot_count_future = 0

        # Number of edges counters
        edges_count_past = 0
        edges_count_present = 0
        edges_count_future = 0

        # For every neighbor of the specific node
        for neighbor in G.neighbors(target_node):
            
            edge_date = G.edges[target_node, neighbor].get('discovery_date')
            if edge_date is None:
                continue

            # ----------------------------------------------------
            # TOTAL EDGES (Before checking neighbor GO)
            # ----------------------------------------------------
            # Past (Target Date - 10,000 = One Year Prior)
            if edge_date <= (target_date - 10000):
                edges_count_past += 1

            # Present (Exact Date)
            if edge_date <= target_date:
                edges_count_present += 1

            # Future (Target Date + 10,000 = One Year Later)
            if edge_date <= (target_date + 10000):
                edges_count_future += 1

            # ----------------------------------------------------
            # NEIGHBOR GO ANNOTATIONS
            # ----------------------------------------------------
            neighbor_date = get_annotation_date(neighbor, go_id)
            if neighbor_date is None:
                continue # Now it only skips the annotated counters, which is correct!

            # Check one year before
            if edge_date <= (target_date - 10000) and neighbor_date <= (target_date - 10000):
                annot_count_past += 1

            # Check the current date
            if edge_date <= target_date and neighbor_date <= target_date:
                annot_count_present += 1

            # Check one year after
            if edge_date <= (target_date + 10000) and neighbor_date <= (target_date + 10000):
                annot_count_future += 1            

        # Save the results for the specific gene
        statistics_data.append({
            'GO_id': go_id,
            'node_id': target_node,
            'node_annotation_date': target_date,
            'total_neighbors_past': edges_count_past,
            'total_neighbors_present': edges_count_present,
            'total_neighbors_future': edges_count_future,
            'annotated_neighbors_past': annot_count_past,
            'annotated_neighbors_present': annot_count_present,
            'annotated_neighbors_future': annot_count_future
        })

stats_df = pd.DataFrame(statistics_data)

print(f"Saving across-time networks statistics for GO {aspect} depth {depth} cutoff {cutoff}...")
stats_df.to_csv(outputstatistics, sep='\t', index=False)
print(f"Across-time networks statistics for GO {aspect} depth {depth} cutoff {cutoff} ready!")