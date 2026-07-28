import os
import gc
import json
import pickle
import ast
import pandas as pd
import numpy as np
import graph_tool.all as gt
import concurrent.futures

aspect = snakemake.wildcards.aspect
target_term_wildcard = snakemake.wildcards.term  
target_term = target_term_wildcard.replace("_", ":") 

network_file = snakemake.input.final_network
top_annot_file = snakemake.input.top_annot_df
bins_json_file = snakemake.input.bins_json
output_matrix_dir = snakemake.output.matrix_dir

num_threads = snakemake.threads 

NUM_WALKS = 1000
STEPS = 2
NUM_PERMUTATIONS = 1000

os.makedirs(output_matrix_dir, exist_ok=True)
print(f"\n--- [Core] Initializing Matrix Generation for {target_term} ({aspect}) ---\n", flush=True)

with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

top_annot_df = pd.read_pickle(top_annot_file)
    
with open(bins_json_file, 'r') as f:
    degree_bins = json.load(f)

# Strict String Casting: Ensures dataset comparison IDs perfectly match graph nodes
unique_gids_raw = list(G_nx.nodes())
n_total_nodes = len(unique_gids_raw)
name_to_gt_id = {str(name): i for i, name in enumerate(unique_gids_raw)}

gene_to_bin = {}
for bin_key, nodes in degree_bins.items():
    for node in nodes:
        gene_to_bin[str(node)] = bin_key

# --- EXTRACT GLOBAL EDGES ONCE ---
edges_by_date = []
for u, v, data in G_nx.edges(data=True):
    date_str = data.get('discovery_date')
    if date_str is not None:
        try:
            edges_by_date.append((int(date_str), name_to_gt_id[str(u)], name_to_gt_id[str(v)]))
        except KeyError:
            pass

edges_by_date.sort(key=lambda x: x[0]) 
del G_nx
gc.collect()

# --- WORKER FUNCTION ---
def worker_compute_walks(nodes_chunk, edges_up_to_date, n_nodes, num_walks, steps):
    g_curr = gt.Graph(directed=False)
    g_curr.add_vertex(n_nodes)
    if len(edges_up_to_date) > 0:
        g_curr.add_edge_list(np.array(edges_up_to_date))

    rng = np.random.default_rng()

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
        non_value = n_nodes 
        random_paths = np.full((N, max_steps + 1), fill_value=non_value, dtype=int)
        for i in range(N):
            random_paths[i, :] = _step(j, -1, max_steps, non_value)
        return random_paths

    chunk_results = {}
    for row_gt_id in nodes_chunk:
        if g_curr.get_out_degrees([row_gt_id])[0] > 0:
            random_paths = walks(row_gt_id, num_walks, steps)
            non_value = n_nodes
            valid_mask = random_paths != non_value
            
            path_lengths = valid_mask.sum(axis=1)
            path_weights = np.repeat((1.0 / path_lengths) / float(num_walks), path_lengths)
            valid_nodes = random_paths[valid_mask]
            
            node_probs = np.bincount(valid_nodes, weights=path_weights, minlength=n_nodes + 1).astype(np.float16)
        else:
            node_probs = np.zeros(n_nodes + 1, dtype=np.float16)
            node_probs[row_gt_id] = 1.0
            
        chunk_results[row_gt_id] = node_probs
    return chunk_results


# --- PROCESS THE SINGLE TARGET TERM ---
term_data = top_annot_df[top_annot_df['GO_id'] == target_term]

if term_data.empty:
    raise ValueError(f"Term {target_term} not found in the annotation dataframe!\n")

term_row = term_data.iloc[0]
annotated_genes_raw = term_row['annotated_genes']

# Safely evaluate if the dataframe stored the list as a raw string
if isinstance(annotated_genes_raw, str):
    try:
        annotated_genes_raw = ast.literal_eval(annotated_genes_raw)
    except (ValueError, SyntaxError):
        pass

# Dynamic Timeline Construction
gene_to_date = {}
for item in annotated_genes_raw:
    if isinstance(item, (tuple, list)) and len(item) == 2:
        g_name, date = str(item[0]), int(item[1])
        gene_to_date[g_name] = min(gene_to_date.get(g_name, date), date)
    else:
        g_name = str(item)
        gt_id = name_to_gt_id.get(g_name)
        if gt_id is not None:
            for e_date, u, v in edges_by_date:
                if u == gt_id or v == gt_id:
                    gene_to_date[g_name] = min(gene_to_date.get(g_name, e_date), e_date)
                    break

unique_dates = sorted(list(set(gene_to_date.values())))
all_term_genes = set(gene_to_date.keys())

# Reset network chronological state
edge_index = 0
total_edges = len(edges_by_date)
current_edges = []

for current_date in unique_dates:
    print(f"--- [{target_term}] Starting date date: {current_date} ---\n", flush=True)
    
    true_annotated_genes = [g for g, d in gene_to_date.items() if d <= current_date]
    future_genes = [g for g, d in gene_to_date.items() if d > current_date]
    
    if not true_annotated_genes or not future_genes:
        print(f"--- [{target_term}] Skipping {current_date}: Lacking source or future target genes.\n", flush=True)
        continue
        
    rng_decoy = np.random.default_rng(seed=int(current_date)) 
    permutation_rows = [] 
    all_required_nodes = set(true_annotated_genes)
    
    for g in true_annotated_genes:
        permutation_rows.append((g, 0))
        
    # 1. OPTIMIZATION: Precompute valid pools ONCE per gene
    precomputed_pools = {}
    for true_gene in true_annotated_genes:
        bin_key = gene_to_bin.get(true_gene)
        if bin_key and bin_key in degree_bins:
            bin_pool = degree_bins[bin_key]
            valid_pool = [n for n in bin_pool if n not in all_term_genes]
            precomputed_pools[true_gene] = valid_pool if valid_pool else bin_pool
        else:
            precomputed_pools[true_gene] = []

    # 2. OPTIMIZATION: Vectorize 1,000 random choices instantly
    for true_gene in true_annotated_genes:
        pool = precomputed_pools[true_gene]
        if pool:
            # Draw all 1,000 permutations in a single C-optimized NumPy call
            decoys = rng_decoy.choice(pool, size=NUM_PERMUTATIONS, replace=True)
            for i, decoy in enumerate(decoys):
                perm_id = i + 1
                permutation_rows.append((decoy, perm_id))
                all_required_nodes.add(decoy)

    unique_gt_ids = list(set([name_to_gt_id[g] for g in all_required_nodes if g in name_to_gt_id]))

    # Advance network state to the current date
    while edge_index < total_edges and edges_by_date[edge_index][0] <= current_date:
        _, u, v = edges_by_date[edge_index]
        current_edges.append([u, v])
        edge_index += 1

    chunk_size = max(1, len(unique_gt_ids) // num_threads)
    node_chunks = [unique_gt_ids[i:i + chunk_size] for i in range(0, len(unique_gt_ids), chunk_size)]
    
    master_row_data = {}
    total_matrix_rows = len(permutation_rows)
    print(f"--- [{target_term}] Computing {len(unique_gt_ids)} unique nodes to populate {total_matrix_rows:,} matrix rows across {num_threads} cores... ---\n", flush=True)
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(worker_compute_walks, chunk, current_edges, n_total_nodes, NUM_WALKS, STEPS) for chunk in node_chunks]
        for future in concurrent.futures.as_completed(futures):
            master_row_data.update(future.result())
            
    # Assemble 3-Tier MultiIndex DataFrame
    row_multi_index = pd.MultiIndex.from_tuples(permutation_rows, names=['Gene_ID', 'Permutation_ID'])
    col_tuples = [(g, 'Future_Target') for g in future_genes]
    col_multi_index = pd.MultiIndex.from_tuples(col_tuples, names=['Target_ID', 'Status'])
    col_gt_ids = [name_to_gt_id[g] for g in future_genes if g in name_to_gt_id]
    
    matrix_data = []
    for row_gene, _ in permutation_rows:
        if row_gene in name_to_gt_id:
            row_gt_id = name_to_gt_id[row_gene]
            probs = master_row_data.get(row_gt_id, np.zeros(n_total_nodes + 1, dtype=np.float16))
            matrix_data.append([probs[target_id] for target_id in col_gt_ids])
        else:
            matrix_data.append([0.0] * len(col_gt_ids))

    df_matrix = pd.DataFrame(matrix_data, index=row_multi_index, columns=col_multi_index, dtype=np.float16)
    
    # Save directly into Snakemake's term directory
    file_path = os.path.join(output_matrix_dir, f"{current_date}.parquet")
    df_matrix.to_parquet(file_path, engine='pyarrow')
    print(f"--- [{target_term}] Saved Matrix: {file_path}! ---\n", flush=True)

print(f"--- [Core] TERM {target_term} COMPLETED! ---\n", flush=True)