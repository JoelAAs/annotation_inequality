import pickle
import networkx as nx
import pandas as pd

graph_path = snakemake.input.network
ancestors_path = snakemake.input.ancestors
output_path = snakemake.output.pairs

with open(graph_path, "rb") as f:
    G = pickle.load(f)

with open(ancestors_path, "rb") as f:
    ancestors_map = pickle.load(f)
    
print(f"--- Extracting leaf pairs... ---\n")

def get_leaf_doids(doid_list, ancestors_map):
    doid_set = set(doid_list)
    leaves = set(doid_list)
    for doid in doid_list:
        ancestors = ancestors_map.get(doid, set())
        # remove from `leaves` any DOID that is an ancestor of `doid`
        leaves -= (ancestors & doid_set)
    return leaves

# --- Extract leaf DOIDs per node and flatten into pairs ---
rows = []
skipped_no_doids = 0

for entrez_id, data in G.nodes(data=True):
    doid_entries = data.get("doids_with_depth", [])
    if not doid_entries:
        skipped_no_doids += 1
        continue

    doids = [entry["doid"] for entry in doid_entries]
    leaves = get_leaf_doids(doids, ancestors_map)

    for doid in leaves:
        rows.append({"entrez_id": str(entrez_id), "doid": doid})

pairs_df = pd.DataFrame(rows).drop_duplicates()

# --- Basic sanity logging (goes to Snakemake's log if configured) ---
print(f"Nodes with no DOID annotations skipped: {skipped_no_doids}")
print(f"Total (entrez_id, doid) leaf pairs: {len(pairs_df)}")
print(f"Unique genes: {pairs_df['entrez_id'].nunique()}")
print(f"Unique DOIDs: {pairs_df['doid'].nunique()}")

pairs_df.to_csv(output_path, index=False)