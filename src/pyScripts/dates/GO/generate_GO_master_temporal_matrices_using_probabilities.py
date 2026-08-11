import os
import gc
import json
import pickle
import ast
import time
import shutil
import pandas as pd
import numpy as np
import scipy.sparse as sp
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
NUM_PERMUTATIONS = 1000

os.makedirs(output_matrix_dir, exist_ok=True)
print(f"\n=======================================================", flush=True)
print(f"--- [START] Generating Matrices for {target_term} ({aspect})", flush=True)
print(f"=======================================================\n", flush=True)

with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

top_annot_df = pd.read_pickle(top_annot_file)
    
with open(bins_json_file, 'r') as f:
    degree_bins = json.load(f)

# Strict String Casting
unique_gids_raw = list(G_nx.nodes())
n_total_nodes = len(unique_gids_raw)
name_to_gt_id = {str(name): i for i, name in enumerate(unique_gids_raw)}

# Timeline Fix
node_annotations = {
    name_to_gt_id[str(node)]: data.get('go_annotations', []) 
    for node, data in G_nx.nodes(data=True)
}

gene_to_bin = {}
for bin_key, nodes in degree_bins.items():
    for node in nodes:
        gene_to_bin[str(node)] = bin_key

# Extract Edges
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

# --- WORKER FUNCTION: SCIPY SPARSE MATRICES ---
def worker_compute_exact_walks(nodes_chunk, edges_up_to_date, n_nodes):
    if len(edges_up_to_date) > 0:
        edges = np.array(edges_up_to_date)
        rows = np.concatenate([edges[:, 0], edges[:, 1]])
        cols = np.concatenate([edges[:, 1], edges[:, 0]])
        data = np.ones(len(rows), dtype=np.float64)
        A = sp.csr_matrix((data, (rows, cols)), shape=(n_nodes + 1, n_nodes + 1))
        A.data = np.ones_like(A.data)
    else:
        A = sp.csr_matrix((n_nodes + 1, n_nodes + 1), dtype=np.float64)

    D = np.array(A.sum(axis=1)).flatten()
    chunk_results = {}
    
    for u in nodes_chunk:
        probs = np.zeros(n_nodes + 1, dtype=np.float64) 
        d_u = D[u]
        
        if d_u == 0:
            probs[u] = 1.0
            chunk_results[int(u)] = probs.astype(np.float16)
            continue
            
        nbrs_u = A.indices[A.indptr[u]:A.indptr[u+1]]
        d_v_all = D[nbrs_u]
        d_v_sub = d_v_all - 1
        
        prob_u_trans = 1.0 / d_u
        dead_mask = (d_v_sub == 0)
        cont_mask = (d_v_sub > 0)
        
        val_dead = prob_u_trans / 2.0
        val_cont = prob_u_trans / 3.0
        
        probs[u] += np.sum(dead_mask) * val_dead + np.sum(cont_mask) * val_cont
        
        if np.any(dead_mask):
            probs[nbrs_u[dead_mask]] += val_dead
        if np.any(cont_mask):
            probs[nbrs_u[cont_mask]] += val_cont
            x = np.zeros(n_nodes + 1, dtype=np.float64)
            x[nbrs_u[cont_mask]] = prob_u_trans / (d_v_sub[cont_mask] * 3.0)
            w_probs = A.dot(x)
            w_probs[u] = 0.0
            probs += w_probs
            
        chunk_results[int(u)] = probs.astype(np.float16)
        
    return chunk_results
# --------------------------------------------------

term_data = top_annot_df[top_annot_df['GO_id'] == target_term]
if term_data.empty:
    raise ValueError(f"Term {target_term} not found in the annotation dataframe!\n")

term_row = term_data.iloc[0]
annotated_genes_raw = term_row['annotated_genes']

if isinstance(annotated_genes_raw, str):
    try:
        annotated_genes_raw = ast.literal_eval(annotated_genes_raw)
    except (ValueError, SyntaxError):
        pass

gene_to_date = {}
for g_name in annotated_genes_raw:
    g_name = str(g_name)
    if g_name in name_to_gt_id:
        gt_id = name_to_gt_id[g_name]
        annotations = node_annotations.get(gt_id, [])
        if annotations: 
            for ann in annotations:
                if ann.get('go_id') == target_term:
                    date = ann.get('first_annotation_date')
                    if date is not None:
                        gene_to_date[g_name] = int(date)
                    break

unique_dates = sorted(list(set(gene_to_date.values())))
all_term_genes = set(gene_to_date.keys())
edge_index = 0
total_edges = len(edges_by_date)
current_edges = []

for current_date in unique_dates:
    # Adding the prefix to easily identify every line in a crowded log file
    log_prefix = f"[{target_term} | {current_date}]"
    
    print(f"\n--- {log_prefix} Initializing... ---", flush=True)
    
    true_annotated_genes = [g for g, d in gene_to_date.items() if d <= current_date]
    future_genes = [g for g, d in gene_to_date.items() if d > current_date]
    
    if not true_annotated_genes or not future_genes:
        print(f"    {log_prefix} [SKIP] Lacking source ({len(true_annotated_genes)}) or future target ({len(future_genes)}) genes.", flush=True)
        continue
        
    print(f"    {log_prefix} [INFO] Source Genes: {len(true_annotated_genes)} | Future Targets: {len(future_genes)}", flush=True)
        
    rng_decoy = np.random.default_rng(seed=int(current_date)) 
    permutation_rows = [] 
    all_required_nodes = set(true_annotated_genes)
    
    for g in true_annotated_genes:
        permutation_rows.append((g, 0))
        
    precomputed_pools = {}
    for true_gene in true_annotated_genes:
        bin_key = gene_to_bin.get(true_gene)
        if bin_key and bin_key in degree_bins:
            bin_pool = degree_bins[bin_key]
            valid_pool = [str(n) for n in bin_pool if str(n) not in all_term_genes]
            precomputed_pools[true_gene] = valid_pool if valid_pool else [str(n) for n in bin_pool]
        else:
            precomputed_pools[true_gene] = []

    for true_gene in true_annotated_genes:
        pool = precomputed_pools[true_gene]
        if pool:
            decoys = rng_decoy.choice(pool, size=NUM_PERMUTATIONS, replace=True)
            for i, decoy in enumerate(decoys):
                perm_id = i + 1
                decoy_str = str(decoy)
                permutation_rows.append((decoy_str, perm_id))
                all_required_nodes.add(decoy_str)

    unique_gt_ids = list(set([name_to_gt_id[g] for g in all_required_nodes if g in name_to_gt_id]))

    while edge_index < total_edges and edges_by_date[edge_index][0] <= current_date:
        _, u, v = edges_by_date[edge_index]
        current_edges.append([u, v])
        edge_index += 1

    g_balance = gt.Graph(directed=False)
    g_balance.add_vertex(n_total_nodes)
    if len(current_edges) > 0:
        g_balance.add_edge_list(np.array(current_edges))
        
    degrees = g_balance.get_out_degrees(unique_gt_ids)
    sorted_indices = np.argsort(degrees)[::-1]
    sorted_gt_ids = np.array(unique_gt_ids)[sorted_indices]
    
    node_chunks = [[] for _ in range(num_threads)]
    for i, node in enumerate(sorted_gt_ids):
        node_chunks[i % num_threads].append(node)
        
    node_chunks = [c for c in node_chunks if len(c) > 0]
    
    master_row_data = {}
    total_matrix_rows = len(permutation_rows)
    
    print(f"    {log_prefix} [COMPUTE] Launching SciPy math on {len(unique_gt_ids)} nodes / {total_matrix_rows:,} rows across {len(node_chunks)} cores...", flush=True)
    start_time = time.time()
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(worker_compute_exact_walks, chunk, current_edges, n_total_nodes) for chunk in node_chunks]
        for future in concurrent.futures.as_completed(futures):
            master_row_data.update(future.result())
            
    elapsed = time.time() - start_time
    print(f"    {log_prefix} [COMPUTE] Finished in {elapsed:.2f} seconds.", flush=True)
            
    row_multi_index = pd.MultiIndex.from_tuples(permutation_rows, names=['Gene_ID', 'Permutation_ID'])
    col_tuples = [(g, 'Future_Target') for g in future_genes]
    col_multi_index = pd.MultiIndex.from_tuples(col_tuples, names=['Target_ID', 'Status'])
    col_gt_ids = [name_to_gt_id[g] for g in future_genes if g in name_to_gt_id]
    
    print(f"    {log_prefix} [BUILD] Constructing float16 DataFrame...", flush=True)
    matrix_data = []
    for row_gene, _ in permutation_rows:
        if row_gene in name_to_gt_id:
            row_gt_id = name_to_gt_id[row_gene]
            probs = master_row_data.get(row_gt_id, np.zeros(n_total_nodes + 1, dtype=np.float16))
            matrix_data.append([probs[target_id] for target_id in col_gt_ids])
        else:
            matrix_data.append([0.0] * len(col_gt_ids))

    df_matrix = pd.DataFrame(matrix_data, index=row_multi_index, columns=col_multi_index, dtype=np.float16)
    file_path = os.path.join(output_matrix_dir, f"{current_date}.parquet")
    
    # Check Disk Space Before Save
    total_bytes, used_bytes, free_bytes = shutil.disk_usage(output_matrix_dir)
    free_gb = free_bytes / (1024**3)
    print(f"    {log_prefix} [DISK] Space available before save: {free_gb:.2f} GB", flush=True)
    
    if free_gb < 0.1:
        print(f"    {log_prefix} [WARNING] DANGER: Less than 100MB remaining on drive! Save may fail.", flush=True)

    print(f"    {log_prefix} [SAVE] Compressing with zstd level 19...", flush=True)
    df_matrix.to_parquet(file_path, engine='pyarrow', compression='zstd', compression_level=19)
    
    # Report File Size
    file_size_mb = os.path.getsize(file_path) / (1024**2)
    print(f"    {log_prefix} [SUCCESS] Saved Matrix: {file_path} | Size: {file_size_mb:.2f} MB", flush=True)

print(f"\n=======================================================", flush=True)
print(f"--- [COMPLETE] ALL DATES FINISHED FOR {target_term} ---", flush=True)
print(f"=======================================================\n", flush=True)