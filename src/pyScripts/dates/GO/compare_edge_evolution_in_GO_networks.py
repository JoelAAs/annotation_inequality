import pandas as pd
import networkx as nx
import pickle
from datetime import datetime
from collections import defaultdict

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.final_network
top_annot_df = snakemake.input.nodes_with_top_5_annotations_pickle
outputstatistics = snakemake.output.edges_evolution_statistics

print(f"Computing Full Time-Zero (True Historian Logic) statistics for GO {aspect} depth {depth} cutoff {cutoff}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
with open(top_annot_df, 'rb') as f:
    top_annot = pickle.load(f)

def get_annotation_date(node_id, target_go):
    annotations = G.nodes[node_id].get('go_annotations', [])
    for ann in annotations:
        if ann.get('go_id') == target_go:
            date = ann.get('first_annotation_date')
            if date is not None:
                return int(date) 
    return None

def get_month_diff(date1_int, date2_int):
    d1 = datetime.strptime(str(date1_int), '%Y%m%d')
    d2 = datetime.strptime(str(date2_int), '%Y%m%d')
    return (d1.year - d2.year) * 12 + (d1.month - d2.month)

statistics_data = []

for index, row in top_annot.iterrows():
    go_id = row['GO_id']
    genes = row['annotated_genes']
    
    # Track the ENTIRE relative timeline
    relative_timeline = defaultdict(lambda: {'total_added': 0, 'annot_added': 0})

    for target_node in genes:
        target_date_int = get_annotation_date(target_node, go_id)
        if target_date_int is None:
            continue

        for neighbor in G.neighbors(target_node):
            edge_date_int = G.edges[target_node, neighbor].get('discovery_date')
            if edge_date_int is None:
                continue
                
            # 1. Get the relative month the edge was discovered
            m_edge = get_month_diff(edge_date_int, target_date_int)
            
            # Event A: The edge is added to the network. Denominator increases.
            relative_timeline[m_edge]['total_added'] += 1
            
            # 2. TRUE HISTORIAN LOGIC: When does the neighbor actually get annotated?
            neighbor_annot_int = get_annotation_date(neighbor, go_id)
            if neighbor_annot_int is not None:
                m_annot = get_month_diff(neighbor_annot_int, target_date_int)
                
                # Event B: The neighbor counts as annotated ONLY when BOTH conditions are met.
                # It joins the numerator at whichever event happens LAST.
                m_activation = max(m_edge, m_annot)
                relative_timeline[m_activation]['annot_added'] += 1

    # Save the aggregated timeline
    for m, counts in sorted(relative_timeline.items()):
        statistics_data.append({
            'GO_id': go_id,
            'relative_month': m,
            'total_edges_added': counts['total_added'],
            'annotated_edges_added': counts['annot_added']
        })

stats_df = pd.DataFrame(statistics_data)
stats_df.to_csv(outputstatistics, sep='\t', index=False)
print(f"Full Time-Zero (True Historian) statistics ready!")