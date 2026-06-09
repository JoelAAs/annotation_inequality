import pandas as pd
import networkx as nx
import pickle
from collections import defaultdict

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.final_network
top_annot_df = snakemake.input.nodes_with_top_5_annotations_pickle
outputstatistics = snakemake.output.edges_evolution_statistics

print(f"Computing Exact Daily Fractions (Evolving Edges / Fixed Annotations) for GO {aspect}...")

with open(network, 'rb') as f:
    G = pickle.load(f)
with open(top_annot_df, 'rb') as f:
    top_annot = pickle.load(f)

def is_ever_annotated(node_id, target_go):
    annotations = G.nodes[node_id].get('go_annotations', [])
    for ann in annotations:
        if ann.get('go_id') == target_go:
            return True
    return False

statistics_data = []

for index, row in top_annot.iterrows():
    go_id = row['GO_id']
    genes = row['annotated_genes']
    
    all_gene_dfs = []
    global_dates = set()
    
    for target_node in genes:
        # Check if the target node itself belongs to this cohort
        if not is_ever_annotated(target_node, go_id):
            continue
            
        gene_events = defaultdict(lambda: {'total_added': 0, 'annot_added': 0})
        has_valid_edges = False
        
        for neighbor in G.neighbors(target_node):
            edge_date_str = G.edges[target_node, neighbor].get('discovery_date')
            if edge_date_str is None:
                continue
                
            has_valid_edges = True
            dt_edge = pd.to_datetime(str(edge_date_str), format='%Y%m%d')
            
            # The denominator increases on the exact date the edge is discovered
            gene_events[dt_edge]['total_added'] += 1
            global_dates.add(dt_edge)
            
            # If the neighbor is EVER annotated, it counts immediately on the edge date
            if is_ever_annotated(neighbor, go_id):
                gene_events[dt_edge]['annot_added'] += 1

        if has_valid_edges:
            df_gene = pd.DataFrame.from_dict(gene_events, orient='index').sort_index()
            df_gene['cum_total'] = df_gene['total_added'].cumsum()
            df_gene['cum_annot'] = df_gene['annot_added'].cumsum()
            
            # Drop rows before the gene had any edges to prevent dividing by zero
            df_gene = df_gene[df_gene['cum_total'] > 0].copy()
            df_gene['fraction'] = df_gene['cum_annot'] / df_gene['cum_total']
            
            # Save just the fraction column for this gene
            all_gene_dfs.append(df_gene[['fraction']])
            
    if not all_gene_dfs or not global_dates:
        continue
        
    # Create a sorted global timeline of every day an edge was discovered for this GO cohort
    sorted_global_dates = sorted(list(global_dates))
    
    reindexed_fractions = []
    for df_g in all_gene_dfs:
        # Align this gene's timeline to the global dates. 
        # .ffill() carries the fraction forward on days where this specific gene had no new edges.
        r_df = df_g.reindex(sorted_global_dates, method='ffill')
        reindexed_fractions.append(r_df['fraction'])
        
    # Stack all genes and calculate the daily mean across the cohort (Node-Centric)
    combined_df = pd.concat(reindexed_fractions, axis=1)
    
    daily_means = combined_df.mean(axis=1)
    
    for date, mean_val in daily_means.items():
        if pd.notna(mean_val):
            statistics_data.append({
                'GO_id': go_id,
                'exact_date': date.strftime('%Y-%m-%d'),
                'mean_fraction': mean_val
            })

stats_df = pd.DataFrame(statistics_data)
stats_df.to_csv(outputstatistics, sep='\t', index=False)
print(f"Absolute Daily Timeline (Evolving Edges / Fixed Annotations) for GO {aspect} statistics ready!")