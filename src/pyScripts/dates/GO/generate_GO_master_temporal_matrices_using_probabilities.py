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
NUM_PERMUTATIONS = 1000

os.makedirs(output_matrix_dir, exist_ok=True)
print(f"\n--- [Core] Initializing Exact Math Matrix Generation for {target_term} ({aspect}) ---\n", flush=True)

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

# --- WORKER FUNCTION: EXACT ANALYTICAL MATH ---
def worker_compute_exact_walks(nodes_chunk, edges_up_to_date, n_nodes):
    g_curr = gt.Graph(directed=False)
    g_curr.add_vertex(n_nodes)
    if len(edges_up_to_date) > 0:
        g_curr.add_edge_list(np.array(edges_up_to_date))

    chunk_results = {}
    for u in nodes_chunk:
        # Compute in 64-bit to perfectly capture tiny graph fractions
        probs = np.zeros(n_nodes + 1, dtype=np.float64) 
        nbrs_u = g_curr.get_out_neighbors(u)
        d_u = len(nbrs_u)

        if d_u == 0:
            probs[u] = 1.0
        else:
            prob_u_transition = 1.0 / d_u
            for v in nbrs_u:
                nbrs_v = g_curr.get_out_neighbors(v)
                nbrs_v = nbrs_v[nbrs_v != u] # Exclude immediate backtracking (Step 2 logic)
                d_v_sub = len(nbrs_v)

                if d_v_sub == 0:
                    # Dead end at step 1: Path length is 2 (nodes u, v)
                    val = prob_u_transition / 2.0
                    probs[u] += val
                    probs[v] += val
                else:
                    # Full 2 steps: Path length is 3 (nodes u, v, w)
                    P = prob_u_transition / d_v_sub
                    val = P / 3.0
                    
                    probs[u] += val * d_v_sub
                    probs[v] += val * d_v_sub
                    probs[nbrs_v] += val

        # Downcast to float16 exactly like original code to save memory
        chunk_results[u] = probs.astype(np.float16)
        
    return chunk_results


# --- PROCESS THE SINGLE TARGET TERM ---
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
    print(f"--- [{target_term}] Starting date: {current_date} ---\n", flush=True)
    
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
        
    precomputed_pools = {}
    for true_gene in true_annotated_genes:
        bin_key = gene_to_bin.get(true_gene)
        if bin_key and bin_key in degree_bins:
            bin_pool = degree_bins[bin_key]
            valid_pool = [n for n in bin_pool if n not in all_term_genes]
            precomputed_pools[true_gene] = valid_pool if valid_pool else bin_pool
        else:
            precomputed_pools[true_gene] = []

    for true_gene in true_annotated_genes:
        pool = precomputed_pools[true_gene]
        if pool:
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

    # --- ADVANCED LOAD BALANCING (ANTI-STRAGGLER) ---
    # Build a quick temporary graph to count node degrees
    g_balance = gt.Graph(directed=False)
    g_balance.add_vertex(n_total_nodes)
    if len(current_edges) > 0:
        g_balance.add_edge_list(np.array(current_edges))
        
    # Sort all required nodes from highest degree to lowest degree
    degrees = g_balance.get_out_degrees(unique_gt_ids)
    sorted_indices = np.argsort(degrees)[::-1]
    sorted_gt_ids = np.array(unique_gt_ids)[sorted_indices]
    
    # Deal nodes into chunks round-robin style to perfectly balance massive hubs across cores
    node_chunks = [[] for _ in range(num_threads)]
    for i, node in enumerate(sorted_gt_ids):
        node_chunks[i % num_threads].append(node)
        
    node_chunks = [c for c in node_chunks if len(c) > 0]
    
    master_row_data = {}
    total_matrix_rows = len(permutation_rows)
    print(f"--- [{target_term}] Computing {len(unique_gt_ids)} nodes (Exact Math) for {total_matrix_rows:,} rows across {len(node_chunks)} balanced cores... ---\n", flush=True)
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(worker_compute_exact_walks, chunk, current_edges, n_total_nodes) for chunk in node_chunks]
        for future in concurrent.futures.as_completed(futures):
            master_row_data.update(future.result())
            
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
    
    file_path = os.path.join(output_matrix_dir, f"{current_date}.parquet")
    df_matrix.to_parquet(file_path, engine='pyarrow')
    print(f"--- [{target_term}] Saved Matrix: {file_path}! ---\n", flush=True)

print(f"--- [Core] TERM {target_term} COMPLETED! ---\n", flush=True)