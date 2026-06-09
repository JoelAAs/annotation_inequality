import pandas as pd
import networkx as nx
import pickle

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
network = snakemake.input.final_network
top_annot_df = snakemake.input.nodes_with_top_5_annotations_pickle
outputstatistics = snakemake.output.edges_evolution_statistics

print(f"Computing Exact Daily Fractions (Static Network / Evolving Annotations) for GO {aspect}...")

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
                return str(date) 
    return None

statistics_data = []

for index, row in top_annot.iterrows():
    go_id = row['GO_id']
    genes = row['annotated_genes']
    
    all_gene_dfs = []
    global_dates = set()
    
    for target_node in genes:
        # Check if the target node itself is properly annotated
        if get_annotation_date(target_node, go_id) is None:
            continue
            
        # Get the total number of neighbors in the modern network
        neighbors = list(G.neighbors(target_node))
        total_neighbors = len(neighbors)
        
        # If the gene has no neighbors in the modern network, skip it
        if total_neighbors == 0:
            continue
            
        # Collect the exact dates neighbors got annotated
        annot_dates = []
        for neighbor in neighbors:
            annot_date_str = get_annotation_date(neighbor, go_id)
            if annot_date_str is not None:
                dt = pd.to_datetime(annot_date_str, format='%Y%m%d')
                annot_dates.append(dt)
                global_dates.add(dt)

        # CALCULATE THIS GENE'S TIMELINE
        if not annot_dates:
            # If the gene NEVER gets an annotated neighbor, its fraction is always 0.0
            df_gene = pd.DataFrame(columns=['fraction']) 
        else:
            # Count how many neighbors get annotated on each specific day
            date_counts = pd.Series(annot_dates).value_counts().sort_index()
            cum_annot = date_counts.cumsum()
            
            df_gene = pd.DataFrame({'fraction': cum_annot / total_neighbors})
            
        all_gene_dfs.append(df_gene)
            
    if not all_gene_dfs or not global_dates:
        continue
        
    # Create a sorted global timeline of every day an annotation occurred for this GO term
    sorted_global_dates = sorted(list(global_dates))
    
    reindexed_fractions = []
    for df_g in all_gene_dfs:
        if df_g.empty:
            # Gene with 0 annotated neighbors gets a flatline of 0.0
            r_df = pd.Series(0.0, index=sorted_global_dates, name='fraction')
        else:
            # Align the gene's timeline to the global dates.
            # .ffill() carries the fraction forward on days with no new annotations.
            # .fillna(0.0) ensures that the time *before* the first neighbor was annotated starts at 0.
            r_df = df_g['fraction'].reindex(sorted_global_dates).ffill().fillna(0.0)
            
        reindexed_fractions.append(r_df)
        
    # Combine all genes and calculate the daily mean across the cohort
    combined_df = pd.concat(reindexed_fractions, axis=1)
    daily_means = combined_df.mean(axis=1)
    
    for date, mean_val in daily_means.items():
        statistics_data.append({
            'GO_id': go_id,
            'exact_date': date.strftime('%Y-%m-%d'),
            'mean_fraction': mean_val
        })

stats_df = pd.DataFrame(statistics_data)
stats_df.to_csv(outputstatistics, sep='\t', index=False)
print(f"Absolute Daily Timeline (Static Network) for GO {aspect} statistics ready!")