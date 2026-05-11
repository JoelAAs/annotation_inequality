import pandas as pd
import networkx as nx
import pickle
from datetime import datetime
from collections import defaultdict

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.final_network
top_annot_df = snakemake.input.nodes_with_top_5_annotations_pickle
outputstatistics = snakemake.output.edges_evolution_statistics

print(f"Computing Full Time-Zero (Time-Traveler) statistics for HDO depth {depth} cutoff {cutoff}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
with open(top_annot_df, 'rb') as f:
    top_annot = pickle.load(f)

def get_annotation_date(node_id, target_do):
    annotations = G.nodes[node_id].get('do_annotations', [])
    for ann in annotations:
        if ann.get('doid') == target_do:
            date = ann.get('first_annotation_date')
            if date is not None:
                return int(date) 
    return None

statistics_data = []

for index, row in top_annot.iterrows():
    do_id = row['DO_id']
    genes = row['annotated_genes']
    
    # Track the ENTIRE relative timeline
    relative_timeline = defaultdict(lambda: {'total_added': 0, 'annot_added': 0})

    for target_node in genes:
        target_date_int = get_annotation_date(target_node, do_id)
        if target_date_int is None:
            continue
            
        target_dt = datetime.strptime(str(target_date_int), '%Y%m%d')

        for neighbor in G.neighbors(target_node):
            edge_date_int = G.edges[target_node, neighbor].get('discovery_date')
            if edge_date_int is None:
                continue
                
            edge_dt = datetime.strptime(str(edge_date_int), '%Y%m%d')
            
            # Calculate precise month difference
            delta_months = (edge_dt.year - target_dt.year) * 12 + (edge_dt.month - target_dt.month)
            
            # TIME TRAVELER LOGIC: Does this neighbor EVER get annotated?
            is_annotated = get_annotation_date(neighbor, do_id) is not None
            
            relative_timeline[delta_months]['total_added'] += 1
            if is_annotated:
                relative_timeline[delta_months]['annot_added'] += 1

    # Save the aggregated timeline
    for m, counts in sorted(relative_timeline.items()):
        statistics_data.append({
            'DO_id': do_id,
            'relative_month': m,
            'total_edges_added': counts['total_added'],
            'annotated_edges_added': counts['annot_added']
        })

stats_df = pd.DataFrame(statistics_data)
stats_df.to_csv(outputstatistics, sep='\t', index=False)
print(f"HDO Full Time-Zero statistics ready!")