import os
# --- IMPLICIT MULTITHREADING BLOCK ---
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import gc
import pickle
import ast
import time
import shutil
import pandas as pd
import numpy as np
import scipy.sparse as sp
import graph_tool.all as gt
import concurrent.futures
import pyarrow as pa
import pyarrow.parquet as pq

# --- SNAKEMAKE I/O ---
aspect = snakemake.wildcards.aspect
target_term_wildcard = snakemake.wildcards.term  
target_term = target_term_wildcard.replace("_", ":") 

network_file = snakemake.input.final_network
top_annot_file = snakemake.input.top_annot_df
output_matrix_dir = snakemake.output.matrix_dir

num_threads = snakemake.threads 
NUM_PERMUTATIONS = 1000

os.makedirs(output_matrix_dir, exist_ok=True)
print(f"\n=======================================================", flush=True)
print(f"--- [START] Generating Matrices for {target_term} ({aspect})", flush=True)
print(f"=======================================================\n", flush=True)

# --- LOAD DATA (Shared read-only for all forks) ---
with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

top_annot_df = pd.read_pickle(top_annot_file)
    
# Strict String Casting
unique_gids_raw = list(G_nx.nodes())
n_total_nodes = len(unique_gids_raw)
name_to_gt_id = {str(name): i for i, name in enumerate(unique_gids_raw)}

# Timeline Extraction & SELF-LOOP REMOVAL
node_annotations = {
    name_to_gt_id[str(node)]: data.get('go_annotations', []) 
    for node, data in G_nx.nodes(data=True)
}

edges_by_date = []
for u, v, data in G_nx.edges(data=True):
    if str(u) == str(v): 
        continue

    date_str = data.get('discovery_date')
    if date_str is not None:
        try:
            edges_by_date.append((int(date_str), name_to_gt_id[str(u)], name_to_gt_id[str(v)]))
        except KeyError:
            pass

edges_by_date.sort(key=lambda x: x[0]) 
del G_nx
gc.collect()

# --- TARGET TERM PREPARATION ---
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

# --- WORKER FUNCTION: PROCESS A SINGLE FULL DATE ---
def process_single_date(current_date):
    log_prefix = f"[{target_term} | {current_date}]"
    
    true_annotated_genes = [g for g, d in gene_to_date.items() if d <= current_date]
    future_genes = [g for g, d in gene_to_date.items() if d > current_date]
    
    if not true_annotated_genes or not future_genes:
        return f"    {log_prefix} [SKIP] Lacking source or future target genes."
        
    print(f"    {log_prefix} [INFO] Processing {len(true_annotated_genes)} sources -> {len(future_genes)} targets...", flush=True)

    # 1. TIMELINE EXTRACTION (Fast slice from sorted list)
    current_edges = []
    for d, u, v in edges_by_date:
        if d <= current_date:
            current_edges.append([u, v])
        else:
            break

    # 2. BUILD HISTORICAL GRAPH & GET EXACT DEGREES
    g_balance = gt.Graph(directed=False)
    g_balance.add_vertex(n_total_nodes)
    if len(current_edges) > 0:
        g_balance.add_edge_list(np.array(current_edges))
        
    current_degrees = g_balance.get_out_degrees(np.arange(n_total_nodes))

    # 3. EXACT DEGREE MATCHING
    degree_dict = {}
    for i, deg in enumerate(current_degrees):
        d = int(deg)
        if d not in degree_dict:
            degree_dict[d] = []
        degree_dict[d].append(str(unique_gids_raw[i]))

    gt_id_to_degree = {i: int(deg) for i, deg in enumerate(current_degrees)}

    # 4. DRAW HISTORICAL PERMUTATIONS
    term_numeric_id = int(target_term.split(":")[-1])
    unique_seed = int(current_date) + term_numeric_id
    rng_decoy = np.random.default_rng(seed=unique_seed)
    
    permutation_rows = [] 
    all_required_nodes = set(true_annotated_genes)
    future_genes_set = set(future_genes)
    
    for g in true_annotated_genes:
        permutation_rows.append((g, 0))
        
    for true_gene in true_annotated_genes:
        true_gt_id = name_to_gt_id[true_gene]
        exact_degree = gt_id_to_degree[true_gt_id]
        
        exact_pool = degree_dict.get(exact_degree, [])
        pool = [n for n in exact_pool if n not in future_genes_set]
        if not pool:
             pool = [str(n) for n in unique_gids_raw if str(n) not in future_genes_set]

        decoys = rng_decoy.choice(pool, size=NUM_PERMUTATIONS, replace=True)
        for i, decoy in enumerate(decoys):
            perm_id = i + 1
            decoy_str = str(decoy)
            permutation_rows.append((decoy_str, perm_id))
            all_required_nodes.add(decoy_str)

    unique_gt_ids = list(set([name_to_gt_id[g] for g in all_required_nodes if g in name_to_gt_id]))

    # 5. BUILD LOCAL SPARSE MATRIX
    if len(current_edges) > 0:
        edges_arr = np.array(current_edges)
        rows = np.concatenate([edges_arr[:, 0], edges_arr[:, 1]])
        cols = np.concatenate([edges_arr[:, 1], edges_arr[:, 0]])
        data = np.ones(len(rows), dtype=np.float64)
        A_main = sp.csr_matrix((data, (rows, cols)), shape=(n_total_nodes + 1, n_total_nodes + 1))
        A_main.data = np.ones_like(A_main.data)
    else:
        A_main = sp.csr_matrix((n_total_nodes + 1, n_total_nodes + 1), dtype=np.float64)

    D_main = np.array(A_main.sum(axis=1)).flatten()

    # 6. SEQUENTIAL MATH EXECUTION (Single Core per Date)
    start_time = time.time()
    master_row_data = {}
    
    for u in unique_gt_ids:
        probs = np.zeros(n_total_nodes + 1, dtype=np.float64) 
        d_u = D_main[u]
        
        if d_u == 0:
            probs[u] = 1.0
            master_row_data[int(u)] = probs.astype(np.float32)
            continue
            
        nbrs_u = A_main.indices[A_main.indptr[u]:A_main.indptr[u+1]]
        d_v_all = D_main[nbrs_u]
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
            x = np.zeros(n_total_nodes + 1, dtype=np.float64)
            x[nbrs_u[cont_mask]] = prob_u_trans / (d_v_sub[cont_mask] * 3.0)
            
            w_probs = A_main.dot(x)
            w_probs[u] = 0.0
            probs += w_probs
            
        master_row_data[int(u)] = probs.astype(np.float32)
        
    math_elapsed = time.time() - start_time
            
    # 7. BUILD DATAFRAME & EXPORT
    total_matrix_rows = len(permutation_rows)
    col_gt_ids = [name_to_gt_id[g] for g in future_genes if g in name_to_gt_id]
    num_targets = len(col_gt_ids)
    
    matrix_array = np.zeros((total_matrix_rows, num_targets), dtype=np.float32, order='F')
    
    for row_idx, (row_gene, _) in enumerate(permutation_rows):
        if row_gene in name_to_gt_id:
            row_gt_id = name_to_gt_id[row_gene]
            probs = master_row_data.get(row_gt_id)
            if probs is not None:
                matrix_array[row_idx, :] = probs[col_gt_ids]
                
    min_prob = np.min(matrix_array)
    max_prob = np.max(matrix_array)

    gene_ids = [str(r[0]) for r in permutation_rows]
    perm_ids = np.array([r[1] for r in permutation_rows], dtype=np.uint16) 

    pa_arrays = [pa.array(gene_ids), pa.array(perm_ids)]
    pa_names = ['Gene_ID', 'Permutation_ID']

    for i, target_gene in enumerate(future_genes):
        pa_arrays.append(pa.array(matrix_array[:, i]))
        pa_names.append(str(target_gene))

    table = pa.Table.from_arrays(pa_arrays, names=pa_names)
    file_path = os.path.join(output_matrix_dir, f"{current_date}.parquet")
    
    save_start_time = time.time()
    pq.write_table(table, file_path, compression='snappy')
    save_elapsed = time.time() - save_start_time
    
    file_size_mb = os.path.getsize(file_path) / (1024**2)
    
    # RAM Cleanup per worker
    del table
    del matrix_array
    del master_row_data
    del A_main
    del g_balance
    gc.collect()
    
    return f"    {log_prefix} [SUCCESS] Math: {math_elapsed:.2f}s | Save: {save_elapsed:.2f}s | File: {file_size_mb:.2f} MB | Max Prob: {max_prob:.4f}"

# --- EXECUTE IN PARALLEL OVER DATES ---
print(f"--- Launching Date-Level Parallelization across {num_threads} cores ---", flush=True)

with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
    # Map dates to workers
    future_to_date = {executor.submit(process_single_date, date): date for date in unique_dates}
    
    for future in concurrent.futures.as_completed(future_to_date):
        date = future_to_date[future]
        try:
            result_msg = future.result()
            print(result_msg, flush=True)
        except Exception as exc:
            print(f"    [{target_term} | {date}] [ERROR] Generated an exception: {exc}", flush=True)

print(f"\n=======================================================", flush=True)
print(f"--- [COMPLETE] ALL DATES FINISHED FOR {target_term} ---", flush=True)
print(f"=======================================================\n", flush=True)