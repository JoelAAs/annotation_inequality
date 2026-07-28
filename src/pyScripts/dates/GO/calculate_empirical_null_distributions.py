import os
import gc
import pickle
import pandas as pd
import numpy as np
import graph_tool.all as gt
import concurrent.futures
import traceback
import re

aspect = snakemake.wildcards.aspect
network_file = snakemake.input.final_network
top_annot_file = snakemake.input.top_annot_df
distances_dir = snakemake.input.distances_dir
output_file = snakemake.output.null_master_file 
num_threads = snakemake.threads 

# Ensure the parent directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)

NUM_PERMUTATIONS = 100  # 100 static random sets, carried across all dates
NUM_WALKS = 100
STEPS = 2

print(f"--- Loading Global GO {aspect} Network for Null Permutations ---")
with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

with open(top_annot_file, 'rb') as f:
    top_annot_df = pickle.load(f)
    
print("Preparing graph_tool canonical mappings and degree bins...")
unique_gids = list(G_nx.nodes())
n_total_nodes = len(unique_gids)
name_to_gt_id = {name: i for i, name in enumerate(unique_gids)}

# Inspect node annotations attribute name safely
sample_node_data = list(G_nx.nodes(data=True))[0][1]
annot_key = next((k for k in ['go_annotations', 'annotations', 'GO_annotations'] if k in sample_node_data), None)

node_annotations = {
    name_to_gt_id[node]: data.get(annot_key, []) 
    for node, data in G_nx.nodes(data=True)
}

# --- GLOBAL DEGREE BINNING ---
global_degrees = np.array([G_nx.degree(n) for n in unique_gids])
bins = np.unique(np.percentile(global_degrees, np.linspace(0, 100, 100)).astype(int))
node_bins = np.digitize(global_degrees, bins)

bin_to_nodes = {}
for i, b in enumerate(node_bins):
    if b not in bin_to_nodes:
        bin_to_nodes[b] = []
    bin_to_nodes[b].append(i)
for b in bin_to_nodes:
    bin_to_nodes[b] = np.array(bin_to_nodes[b])

# Extract edges with dates for the temporal graph rebuild
edges_by_date = []
for u, v, data in G_nx.edges(data=True):
    date_str = data.get('discovery_date') or data.get('date')
    if date_str is not None:
        try:
            edges_by_date.append((int(date_str), name_to_gt_id[u], name_to_gt_id[v]))
        except ValueError:
            pass

edges_by_date.sort(key=lambda x: x[0]) 

# Memory wipe
print("Wiping NetworkX graph from memory...")
del G_nx
gc.collect()

# --- WORKER FUNCTION ---
def process_go_term_null(row_data):
    go_id, raw_genes = row_data
    numeric_id = int(go_id.split(':')[1])
    rng = np.random.default_rng(42 + numeric_id)
    
    go_id_safe = go_id.replace(':', '_')
    true_dist_file = os.path.join(distances_dir, f'{go_id_safe}_mean_distances_from_annotated.csv')
    
    if not os.path.exists(true_dist_file):
        return f"Skipped id: {go_id} (No true distance file found)"
        
    df_true = pd.read_csv(true_dist_file, sep='\t')
    
    # INDESTRUCTIBLE PARSING: Force extract numbers in case Pandas loaded the array as a string
    genes = re.findall(r"['\"]?(\d+)['\"]?", str(raw_genes))
    
    gene_to_date = {}
    for gene in genes:
        matched_gene = None
        if gene in name_to_gt_id:
            matched_gene = gene
        elif str(gene) in name_to_gt_id:
            matched_gene = str(gene)
        elif str(gene).isdigit() and int(gene) in name_to_gt_id:
            matched_gene = int(gene)
            
        if matched_gene is not None:
            gt_id = name_to_gt_id[matched_gene]
            annotations = node_annotations.get(gt_id, [])
            if annotations: 
                for ann in annotations:
                    if isinstance(ann, dict):
                        ann_go = ann.get('go_id') or ann.get('GO_id') or ann.get('GO')
                        if ann_go == go_id:
                            date = ann.get('first_annotation_date') or ann.get('date') or ann.get('discovery_date')
                            if date is not None:
                                gene_to_date[matched_gene] = int(date) 
                            break
                        
    if not gene_to_date:
        print(f"[{go_id}] Warning: No matching dates found for extracted genes {genes[:5]}...")
        return f"Skipped id: {go_id} (No dates)"

    unique_dates = sorted(df_true['Date'].unique().tolist())
    
    # --- STATIC TEMPORAL MAPPINGS ---
    all_go_genes = list(gene_to_date.keys())
    null_mappings = []
    
    for p in range(NUM_PERMUTATIONS):
        mapping_for_p = {}
        for a_gene in all_go_genes:
            global_idx = name_to_gt_id[a_gene]
            b = node_bins[global_idx]
            pool = bin_to_nodes[b]
            mapping_for_p[a_gene] = rng.choice(pool)
        null_mappings.append(mapping_for_p)

    # Initialize dynamic graph
    g_curr = gt.Graph(directed=False)
    g_curr.add_vertex(n_total_nodes)
    
    def _step(j, from_node, remaining_steps, non_value):
        if remaining_steps < 1:
            return [j]
        else:
            nbrs_j = g_curr.get_out_neighbors(j)
            nbrs_j = nbrs_j[nbrs_j != from_node]
            if not nbrs_j.any(): 
                return [j] + [non_value] * remaining_steps
            else:
                next_j = nbrs_j[rng.integers(len(nbrs_j))]
                return [j] + _step(next_j, j, remaining_steps - 1, non_value)

    def walks(j, N, max_steps):
        non_value = n_total_nodes 
        random_paths = np.full((N, max_steps + 1), fill_value=non_value, dtype=int)
        for i in range(N):
            random_paths[i, :] = _step(j, -1, max_steps, non_value)
        return random_paths

    edge_index = 0
    total_edges = len(edges_by_date)
    final_results = []

    print(f"\n--- [{go_id}] Processing GO {aspect} id: {go_id} ---")

    for current_date in unique_dates:
        edges_to_add = []
        while edge_index < total_edges and edges_by_date[edge_index][0] <= current_date:
            _, u, v = edges_by_date[edge_index]
            edges_to_add.append([u, v])
            edge_index += 1
            
        if edges_to_add:
            g_curr.add_edge_list(np.array(edges_to_add))
            
        df_true_date = df_true[df_true['Date'] == current_date]
        current_genes = [g for g, d in gene_to_date.items() if d <= current_date and g in name_to_gt_id]
        
        future_genes_raw = df_true_date['Future_Gene'].tolist()
        future_genes = []
        future_gt_ids = []
        valid_indices = []
        
        for i, f_gene in enumerate(future_genes_raw):
            matched_f = None
            if f_gene in name_to_gt_id:
                matched_f = f_gene
            elif str(f_gene) in name_to_gt_id:
                matched_f = str(f_gene)
            elif str(f_gene).isdigit() and int(f_gene) in name_to_gt_id:
                matched_f = int(f_gene)
                
            if matched_f is not None:
                future_genes.append(f_gene)
                future_gt_ids.append(name_to_gt_id[matched_f])
                valid_indices.append(i)
                
        true_adj_values = df_true_date['Mean_Probability_From_Annotated'].values[valid_indices]
        
        if not current_genes or not future_genes:
            print(f"[{go_id}] Skipped Date {current_date}: Current Genes={len(current_genes)}, Future Genes={len(future_genes)}")
            continue
            
        null_matrices = np.zeros((NUM_PERMUTATIONS, len(future_genes)))
        
        for p in range(NUM_PERMUTATIONS):
            sum_probs = np.zeros(len(future_genes))
            num_random_annotated = len(current_genes)
            
            random_annotated_gt_ids = [null_mappings[p][a_gene] for a_gene in current_genes]
            
            for r_gt_id in random_annotated_gt_ids:
                if g_curr.get_out_degrees([r_gt_id])[0] > 0:
                    random_paths = walks(r_gt_id, NUM_WALKS, STEPS)
                    non_value = n_total_nodes
                    valid_mask = random_paths != non_value
                    path_lengths = valid_mask.sum(axis=1)
                    path_weights = np.repeat((1.0 / path_lengths) / float(NUM_WALKS), path_lengths)
                    valid_nodes = random_paths[valid_mask]
                    node_probs = np.bincount(valid_nodes, weights=path_weights, minlength=n_total_nodes + 1)
                else:
                    node_probs = np.zeros(n_total_nodes + 1)
                    node_probs[r_gt_id] = 1.0
                    
                sum_probs += node_probs[future_gt_ids]
                
            null_matrices[p, :] = sum_probs / num_random_annotated
            
        null_means = null_matrices.mean(axis=0)
        null_stds = null_matrices.std(axis=0) + 1e-12 
        z_scores = (true_adj_values - null_means) / null_stds
        p_values = (null_matrices >= true_adj_values[None, :]).mean(axis=0)
        
        for idx, f_gene in enumerate(future_genes):
            final_results.append({
                'GO_ID': go_id, 
                'Date': current_date,
                'Future_Gene': f_gene,
                'True_Adjacency': true_adj_values[idx],
                'Null_Mean': null_means[idx],
                'Null_Std': null_stds[idx],
                'Z_Score': z_scores[idx],
                'Empirical_P_Value': p_values[idx]
            })

    if final_results:
        print(f"[{go_id}] Null permutations completed! Generated {len(final_results)} records.")
        return pd.DataFrame(final_results)
    else:
        print(f"[{go_id}] Error: No valid data generated for any dates.")
        return f"[{go_id}] No valid data."


# --- MULTIPROCESSING POOL CONTROL ---
if __name__ == '__main__':
    tasks = [(row['GO_id'], row['annotated_genes']) for index, row in top_annot_df.iterrows()]
    
    print(f"\nLaunching parallel processing with {num_threads} cores...\n")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        all_dfs = list(executor.map(process_go_term_null, tasks))
        
    valid_dfs = [df for df in all_dfs if isinstance(df, pd.DataFrame)]
    
    if valid_dfs:
        master_df = pd.concat(valid_dfs, ignore_index=True)
        print(f"--- [Core] SUCCESS: Master Null Distribution compiled successfully! ---")
    else:
        master_df = pd.DataFrame(columns=[
            'GO_ID', 'Date', 'Future_Gene', 'True_Adjacency', 
            'Null_Mean', 'Null_Std', 'Z_Score', 'Empirical_P_Value'
        ])
        print(f"--- [Core] WARNING: No valid target/future pairs were processed. Writing empty template file. ---")
        
    try:
        master_df.to_csv(output_file, sep='\t', index=False)
        print(f"--- [Core] File successfully written to {output_file} ---")
    except Exception as e:
        print(f"--- [Core] ERROR writing file to {output_file}: {e} ---")
        traceback.print_exc()
        raise e