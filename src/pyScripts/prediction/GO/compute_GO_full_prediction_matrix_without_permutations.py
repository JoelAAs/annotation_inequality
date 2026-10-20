import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import gc
import time
import pickle
import numpy as np
import scipy.sparse as sp
import pyarrow as pa
import pyarrow.parquet as pq
import concurrent.futures

network_file = snakemake.input.final_network
output_matrix_file = snakemake.output.full_matrix
aspect = snakemake.wildcards.aspect
num_threads = snakemake.threads

os.makedirs(os.path.dirname(output_matrix_file), exist_ok=True)
print(f"\n=======================================================", flush=True)
print(f"--- [START] Generating Full {aspect} Prediction Matrix", flush=True)
print(f"=======================================================\n", flush=True)

# --- LOAD DATA ---
load_start = time.time()
with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

unique_gids_raw = list(G_nx.nodes())
n_total_nodes = len(unique_gids_raw)
name_to_id = {str(name): i for i, name in enumerate(unique_gids_raw)}

print(f"    [INFO] Extracted {n_total_nodes} nodes from network in {time.time() - load_start:.2f}s.", flush=True)

edges = []
for u, v in G_nx.edges():
    if str(u) != str(v):
        edges.append([name_to_id[str(u)], name_to_id[str(v)]])

del G_nx
gc.collect()

# --- BUILD GLOBAL SPARSE MATRIX ---
print("    [INFO] Building Global Sparse Adjacency Matrix...", flush=True)
sparse_start = time.time()

if len(edges) > 0:
    edges_arr = np.array(edges)
    rows = np.concatenate([edges_arr[:, 0], edges_arr[:, 1]])
    cols = np.concatenate([edges_arr[:, 1], edges_arr[:, 0]])
    data = np.ones(len(rows), dtype=np.float64)
    A_main = sp.csr_matrix((data, (rows, cols)), shape=(n_total_nodes, n_total_nodes))
    A_main.data = np.ones_like(A_main.data)
else:
    A_main = sp.csr_matrix((n_total_nodes, n_total_nodes), dtype=np.float64)

D_main = np.array(A_main.sum(axis=1)).flatten()
print(f"    [INFO] Sparse Matrix built in {time.time() - sparse_start:.2f}s.", flush=True)

# --- WORKER FUNCTION ---
def process_chunk(args):
    """Processes a chunk of source nodes with progress logging."""
    chunk_id, chunk_nodes = args
    total_in_chunk = len(chunk_nodes)
    chunk_results = np.zeros((total_in_chunk, n_total_nodes), dtype=np.float16)
    
    print(f"    [Worker {chunk_id}] Started processing {total_in_chunk} nodes.", flush=True)
    
    for idx, u in enumerate(chunk_nodes):
        if (idx + 1) % 500 == 0 or (idx + 1) == total_in_chunk:
            print(f"    [Worker {chunk_id}] Processed {idx + 1}/{total_in_chunk} nodes...", flush=True)

        d_u = D_main[u]
        if d_u == 0:
            # Isolated node: retains the entire probability budget
            probs = np.zeros(n_total_nodes, dtype=np.float64)
            probs[u] = 1.0
            chunk_results[idx, :] = probs.astype(np.float16)
            continue 

        probs = np.zeros(n_total_nodes, dtype=np.float64) 
        nbrs_u = A_main.indices[A_main.indptr[u]:A_main.indptr[u+1]]
        d_v_all = D_main[nbrs_u]
        d_v_sub = d_v_all - 1
        
        prob_u_trans = 1.0 / d_u
        dead_mask = (d_v_sub == 0)
        cont_mask = (d_v_sub > 0)
        
        val_dead = prob_u_trans / 2.0
        val_cont = prob_u_trans / 3.0
        
        # 1-Hop Math
        if np.any(dead_mask):
            probs[nbrs_u[dead_mask]] += val_dead
            
        if np.any(cont_mask):
            probs[nbrs_u[cont_mask]] += val_cont
            
            # 2-Hop Math Push
            x = np.zeros(n_total_nodes, dtype=np.float64)
            x[nbrs_u[cont_mask]] = prob_u_trans / (d_v_sub[cont_mask] * 3.0)
            
            w_probs = A_main.dot(x)
            w_probs[u] = 0.0 # Prevent backtracking
            probs += w_probs
        
        # --- DIAGONAL CALCULATION ---
        # Assign all undistributed probability back to the source node itself
        probs[u] = max(0.0, 1.0 - np.sum(probs))
        
        chunk_results[idx, :] = probs.astype(np.float16)
        
    print(f"    [Worker {chunk_id}] Completed all {total_in_chunk} nodes.", flush=True)
    return chunk_id, chunk_nodes, chunk_results

# --- PARALLEL EXECUTION ---
print(f"\n--- Launching Parallelization across {num_threads} cores ---", flush=True)

all_nodes = np.arange(n_total_nodes)
chunks = np.array_split(all_nodes, num_threads)

final_matrix = np.zeros((n_total_nodes, n_total_nodes), dtype=np.float16)

math_start = time.time()
with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
    futures = {executor.submit(process_chunk, (i, chunk)): i for i, chunk in enumerate(chunks)}
    
    for future in concurrent.futures.as_completed(futures):
        try:
            chunk_id, chunk_nodes, chunk_res = future.result()
            print(f"    [MAIN] Merging Worker {chunk_id} results into final matrix...", flush=True)
            final_matrix[chunk_nodes, :] = chunk_res
            print(f"    [MAIN] Successfully merged Worker {chunk_id}.", flush=True)
        except Exception as exc:
            chunk_id = futures[future]
            print(f"    [ERROR] Worker {chunk_id} generated an exception: {exc}", flush=True)

math_elapsed = time.time() - math_start
print(f"--- Math completed in {math_elapsed:.2f}s ---\n", flush=True)

# --- MATRIX VALIDATION & SPARSITY ---
print("    [INFO] Validating Matrix and Calculating Sparsity...", flush=True)
max_prob = np.max(final_matrix)
print(f"           -> Global Maximum Probability: {max_prob:.4f}", flush=True)

if max_prob == 0.0:
    print("           -> [WARNING] The entire matrix is composed of zeros!", flush=True)
elif max_prob > 1.01: # Limit set to 1.0 because isolated nodes retain exactly 1.0
    print("           -> [WARNING] Maximum probability exceeds the theoretical limit of 1.0!", flush=True)
else:
    print("           -> [OK] Matrix values are within expected mathematical bounds.", flush=True)

total_cells = n_total_nodes * n_total_nodes
non_zero_cells = np.count_nonzero(final_matrix)
zero_cells = total_cells - non_zero_cells
sparsity_pct = (zero_cells / total_cells) * 100

print(f"           -> Total Cells: {total_cells:,}", flush=True)
print(f"           -> Zero Cells:  {zero_cells:,} ({sparsity_pct:.2f}%)", flush=True)
print("")

# --- BUILD DATAFRAME & EXPORT ---
print("    [INFO] Building PyArrow Table for Parquet Export...", flush=True)
table_start = time.time()

str_gene_ids = [str(g) for g in unique_gids_raw]
pa_arrays = [pa.array(str_gene_ids)]
pa_names = ['Gene_ID']

for i in range(n_total_nodes):
    pa_arrays.append(pa.array(final_matrix[:, i]))
    pa_names.append(str_gene_ids[i])

table = pa.Table.from_arrays(pa_arrays, names=pa_names)
print(f"    [INFO] PyArrow Table built in {time.time() - table_start:.2f}s.", flush=True)

print("    [INFO] Saving to Parquet (Snappy compression)...", flush=True)
save_start = time.time()
pq.write_table(table, output_matrix_file, compression='snappy')
save_elapsed = time.time() - save_start

file_size_mb = os.path.getsize(output_matrix_file) / (1024**2)

print(f"\n=======================================================", flush=True)
print(f"--- [COMPLETE] FULL MATRIX SAVED ---", flush=True)
print(f"--- Save Time: {save_elapsed:.2f}s | File Size: {file_size_mb:.2f} MB", flush=True)
print(f"=======================================================\n", flush=True)