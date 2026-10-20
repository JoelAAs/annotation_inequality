import os
import pickle
import pandas as pd

network_file = snakemake.input.final_network
output_y_file = snakemake.output.y_file

# Reconstruct the GO ID
target_term_wildcard = snakemake.wildcards.term
target_term = target_term_wildcard.replace("_", ":")

os.makedirs(os.path.dirname(output_y_file), exist_ok=True)
print(f"\n=======================================================", flush=True)
print(f"--- [START] Extracting Labels (Y) for {target_term}", flush=True)
print(f"=======================================================\n", flush=True)

# --- LOAD NETWORK ---
print(f"    [INFO] Loading final network topology...", flush=True)
with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

# --- EXTRACT LABELS ---
print(f"    [INFO] Parsing node attributes for {target_term}...", flush=True)
unique_gids = list(G_nx.nodes())
labels = []

for node in unique_gids:
    node_data = G_nx.nodes[node]
    annotations = node_data.get('go_annotations', [])
    
    # Check if the exact target term is present in the node's annotations
    has_term = any(ann.get('go_id') == target_term for ann in annotations)
    labels.append(1 if has_term else 0)

# Free up RAM immediately
del G_nx

# --- BUILD DATAFRAME & EXPORT ---
df_y = pd.DataFrame({
    "Gene_ID": [str(g) for g in unique_gids],
    "Label": labels
})

# Cast to memory-efficient types
df_y['Gene_ID'] = df_y['Gene_ID'].astype(str)
df_y['Label'] = df_y['Label'].astype('int8')

num_positives = df_y['Label'].sum()
num_negatives = len(df_y) - num_positives

print(f"           -> Total Nodes: {len(df_y)}", flush=True)
print(f"           -> Positive Labels: {num_positives} | Negative Labels: {num_negatives}", flush=True)

if num_positives == 0:
    print(f"           -> [WARNING] No positive annotations found for {target_term}!", flush=True)

print(f"    [INFO] Saving to Parquet (Snappy compression)...", flush=True)
df_y.to_parquet(output_y_file, engine='pyarrow', compression='snappy', index=False)

file_size_kb = os.path.getsize(output_y_file) / 1024

print(f"\n=======================================================", flush=True)
print(f"--- [COMPLETE] LABELS SAVED ---", flush=True)
print(f"--- File Size: {file_size_kb:.2f} KB", flush=True)
print(f"=======================================================\n", flush=True)