import os
import gc
import pickle
import pandas as pd
import numpy as np
import graph_tool.all as gt
import concurrent.futures

aspect = snakemake.wildcards.aspect
network_file = snakemake.input.final_network
top_annot_file = snakemake.input.top_annot_df
output_dir = snakemake.output.matrix_dir
num_threads = snakemake.threads 

print(f"Loading GO {aspect} NetworkX graph...")
with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

with open(top_annot_file, 'rb') as f:
    top_annot_df = pickle.load(f)
    
os.makedirs(output_dir, exist_ok=True)
NUM_WALKS = 1000
STEPS = 2

print("Preparing graph_tool canonical mappings...")
unique_gids = list(G_nx.nodes())
n_total_nodes = len(unique_gids)
name_to_gt_id = {name: i for i, name in enumerate(unique_gids)}

node_annotations = {
    name_to_gt_id[node]: data.get('go_annotations', []) 
    for node, data in G_nx.nodes(data=True)
}

edges_by_date = []
for u, v, data in G_nx.edges(data=True):
    date_str = data.get('discovery_date')
    if date_str is not None:
        try:
            edges_by_date.append((int(date_str), name_to_gt_id[u], name_to_gt_id[v]))
        except ValueError:
            pass

edges_by_date.sort(key=lambda x: x[0]) 

# --- MEMORY WIPE ---
print("Wiping NetworkX graph from memory to free up RAM...")
del G_nx
gc.collect()

# --- WORKER FUNCTION ---
def process_go_term(row_data):
    go_id, genes = row_data
    
    numeric_id = int(go_id.split(':')[1])
    rng = np.random.default_rng(42 + numeric_id)
    
    print(f"\n--- [Core] Calculating temporal matrices for GO {aspect} id: {go_id} ---")
    
    gene_to_date = {}
    for gene in genes:
        if gene in name_to_gt_id:
            gt_id = name_to_gt_id[gene]
            annotations = node_annotations.get(gt_id, [])
            
            if annotations: 
                for ann in annotations:
                    if ann.get('go_id') == go_id:
                        date = ann.get('first_annotation_date')
                        if date is not None:
                            gene_to_date[gene] = int(date) 
                        break
                    
    if not gene_to_date:
        return f"Skipped id: {go_id} (No dates)"

    safe_go_id = go_id.replace(':', '_')
    go_dir = os.path.join(output_dir, safe_go_id)
    os.makedirs(go_dir, exist_ok=True)

    unique_dates = sorted(list(set(gene_to_date.values())))

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

    for current_date in unique_dates:
        current_genes = [g for g, d in gene_to_date.items() if d <= current_date]
        future_genes = [g for g, d in gene_to_date.items() if d > current_date]
        
        if not current_genes:
            continue 
            
        all_target_genes = current_genes + future_genes
        print(f"  [{go_id}] Processing date {current_date} | Walking {len(all_target_genes)} total genes...")
        
        edges_to_add = []
        while edge_index < total_edges and edges_by_date[edge_index][0] <= current_date:
            _, u, v = edges_by_date[edge_index]
            edges_to_add.append([u, v])
            edge_index += 1
            
        if edges_to_add:
            g_curr.add_edge_list(np.array(edges_to_add))
            
        col_tuples = [('Annotated', g) for g in current_genes] + [('Future', g) for g in future_genes]
        multi_index = pd.MultiIndex.from_tuples(col_tuples, names=['Status', 'Gene_ID'])
        
        matrix_data = []
        
        for row_gene in all_target_genes:
            row_gt_id = name_to_gt_id[row_gene]
            
            if g_curr.get_out_degrees([row_gt_id])[0] > 0:
                random_paths = walks(row_gt_id, NUM_WALKS, STEPS)
                
                non_value = n_total_nodes
                valid_mask = random_paths != non_value
                
                # Get the exact path length (|pi|) for every single walk
                path_lengths = valid_mask.sum(axis=1)
                
                # Calculate the theoretical weight: (1 / |pi|) * (1 / NUM_WALKS)
                # We repeat this weight for every valid node in that specific path
                path_weights = np.repeat((1.0 / path_lengths) / float(NUM_WALKS), path_lengths)
                
                # Sum the weighted probabilities globally
                valid_nodes = random_paths[valid_mask]
                node_probs = np.bincount(valid_nodes, weights=path_weights, minlength=n_total_nodes + 1)
                
            else:
                # Isolated Node
                node_probs = np.zeros(n_total_nodes + 1)
                node_probs[row_gt_id] = 1.0
                
            row_data = []
            for status, col_gene in col_tuples:
                col_gt_id = name_to_gt_id[col_gene]
                row_data.append(node_probs[col_gt_id])
                
            matrix_data.append(row_data)
            
        df_matrix = pd.DataFrame(matrix_data, index=multi_index, columns=multi_index)
        file_path = os.path.join(go_dir, f"{current_date}.csv")
        df_matrix.to_csv(file_path)

    return f"[{go_id}] temporal matrices completed!"

# --- MULTIPROCESSING POOL CONTROL ---
if __name__ == '__main__':
    tasks = [(row['GO_id'], row['annotated_genes']) for index, row in top_annot_df.iterrows()]
    
    print(f"\nLaunching parallel processing with {num_threads} cores...\n")
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = executor.map(process_go_term, tasks)
        
    print(f"--- [Core] SUCCESS: All GO {aspect} random walk temporal matrices saved to {output_dir}! ---")