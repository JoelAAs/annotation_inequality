import os
import json
import pickle
import numpy as np

aspect = snakemake.wildcards.aspect
network_file = snakemake.input.final_network
output_bins_json = snakemake.output.bins_json

print(f"\n--- [Core] Initializing percentile-based degree binning for GO {aspect}... ---\n")

with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

unique_gids = list(G_nx.nodes())
print(f"--- [Core] Loaded network with {len(unique_gids)} unique nodes.\n")

# --- PERCENTILE BINS ---
global_degrees = np.array([G_nx.degree(n) for n in unique_gids])
bins = np.unique(np.percentile(global_degrees, np.linspace(0, 100, 100)).astype(int))
node_bins = np.digitize(global_degrees, bins)

degree_bins = {}
for i, b in enumerate(node_bins):
    bin_key = f"bin_{int(b)}"
    if bin_key not in degree_bins:
        degree_bins[bin_key] = []
    degree_bins[bin_key].append(unique_gids[i])

sorted_degree_bins = {k: degree_bins[k] for k in sorted(degree_bins.keys(), key=lambda x: int(x.split('_')[1]))}

os.makedirs(os.path.dirname(output_bins_json), exist_ok=True)
with open(output_bins_json, 'w') as f:
    json.dump(sorted_degree_bins, f, indent=4)

print(f"--- [Core] SUCCESS: Exported {len(degree_bins)} unique degree bins to {output_bins_json} ---\n")