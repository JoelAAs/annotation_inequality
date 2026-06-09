import pandas as pd
import networkx as nx
import pickle
import numpy as np
from datetime import datetime
from collections import defaultdict

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.final_network
top_annot_df = snakemake.input.nodes_with_top_5_annotations_pickle
outputstatistics = snakemake.output.edges_evolution_statistics

print(f"Computing Node-Centric Mean Fractions (True Historian Logic) for GO {aspect} depth {depth} cutoff {cutoff}...")

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
    
    # Store events for each individual gene: {gene_id: {month: {'total': x, 'annot': y}}}
    all_gene_events = {}
    
    # Keep track of the absolute min and max months for this specific GO_id
    global_min_m = float('inf')
    global_max_m = float('-inf')
    
    for target_node in genes:
        target_date_int = get_annotation_date(target_node, go_id)
        if target_date_int is None:
            continue
            
        gene_events = defaultdict(lambda: {'total_added': 0, 'annot_added': 0})
        has_valid_edges = False
        
        for neighbor in G.neighbors(target_node):
            edge_date_int = G.edges[target_node, neighbor].get('discovery_date')
            if edge_date_int is None:
                continue
                
            has_valid_edges = True
            m_edge = get_month_diff(edge_date_int, target_date_int)
            gene_events[m_edge]['total_added'] += 1
            
            # Update global timeline bounds
            global_min_m = min(global_min_m, m_edge)
            global_max_m = max(global_max_m, m_edge)
            
            neighbor_annot_int = get_annotation_date(neighbor, go_id)
            if neighbor_annot_int is not None:
                m_annot = get_month_diff(neighbor_annot_int, target_date_int)
                m_activation = max(m_edge, m_annot)
                gene_events[m_activation]['annot_added'] += 1
                
                global_max_m = max(global_max_m, m_activation)

        if has_valid_edges:
            all_gene_events[target_node] = gene_events
            
    if not all_gene_events or global_min_m == float('inf'):
        continue
        
    # Dictionary to collect all gene fractions at a specific month
    fractions_per_month = defaultdict(list)
    
    # Process each gene's timeline individually
    for gene, events in all_gene_events.items():
        cum_total = 0
        cum_annot = 0
        
        # Walk through the timeline month by month for this gene
        for m in range(global_min_m, global_max_m + 1):
            if m in events:
                cum_total += events[m]['total_added']
                cum_annot += events[m]['annot_added']
            
            # Only calculate a fraction if the gene has at least one edge by this month
            if cum_total > 0:
                fractions_per_month[m].append(cum_annot / cum_total)
                
    # Calculate the mean of those fractions for each month
    for m in range(global_min_m, global_max_m + 1):
        if len(fractions_per_month[m]) > 0:
            mean_frac = np.mean(fractions_per_month[m])
            statistics_data.append({
                'GO_id': go_id,
                'relative_month': m,
                'mean_fraction': mean_frac
            })

stats_df = pd.DataFrame(statistics_data)
stats_df.to_csv(outputstatistics, sep='\t', index=False)
print("Node-Centric Mean Time-Zero statistics ready!")