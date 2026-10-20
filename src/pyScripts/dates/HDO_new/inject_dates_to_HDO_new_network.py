import pandas as pd
import networkx as nx
import pickle

file_network = snakemake.input.network
file_dates = snakemake.input.only_associated
file_output = snakemake.output.network_with_dates

print("--- [INFO] Loading network and temporal data... ---")
# Load NetworkX graph
with open(file_network, 'rb') as f:
    G = pickle.load(f)

# Load Parquet file with dates
df_dates = pd.read_parquet(file_dates)

# Prepare lookup dictionary-
print("--- [INFO] Processing temporal associations (finding the oldest date per gene-DOID)... ---")

# Ensure 'entrez_id' is treated as a string to match NetworkX node IDs
df_dates['entrez_id'] = df_dates['entrez_id'].astype(str)

# Convert Publication_Date to integer for correct min() evaluation
df_dates['Publication_Date'] = df_dates['Publication_Date'].astype(int)

# Group by gene and DOID, taking the minimum (oldest) date
oldest_dates = df_dates.groupby(['entrez_id', 'doid'])['Publication_Date'].min().reset_index()

# Create a nested dictionary for fast O(1) lookups: lookup[entrez_id][doid] = oldest_date
date_lookup = {}
for _, row in oldest_dates.iterrows():
    gene = row['entrez_id']
    doid = row['doid']
    date = int(row['Publication_Date'])
    
    if gene not in date_lookup:
        date_lookup[gene] = {}
    date_lookup[gene][doid] = date

# Inject dates into network
print("--- [INFO] Injecting dates into the network nodes... ---")
updated_nodes = 0
total_annotations = 0
dated_annotations = 0

for node_id, node_data in G.nodes(data=True):
    # Only process nodes that have the original annotation list
    if 'doids_with_depth' in node_data:
        old_list = node_data['doids_with_depth']
        new_list = []
        
        for annot in old_list:
            doid = annot['doid']
            depth = annot['depth']
            
            # Check if we have a date for this specific gene-DOID pair
            found_date = None
            if node_id in date_lookup and doid in date_lookup[node_id]:
                found_date = date_lookup[node_id][doid]
                dated_annotations += 1
                
            # Append to the new structure requested by the user
            new_list.append({
                'doid': doid,
                'depth': depth,
                'first_annotation_date': found_date
            })
            total_annotations += 1
        
        # Replace the old key with the new key 'do_annotations'
        node_data['do_annotations'] = new_list
        
        # Remove the old key
        del node_data['doids_with_depth']
        
        updated_nodes += 1

# Final states and Save
print(f"    > Nodes updated: {updated_nodes}")
print(f"    > Total DOID annotations processed: {total_annotations}")
print(f"    > Annotations successfully dated: {dated_annotations}")
print(f"    > Annotations left with 'None': {total_annotations - dated_annotations}")

print("\n--- [INFO] Saving the updated network... ---")
with open(file_output, 'wb') as f:
    pickle.dump(G, f)

print("--- [INFO] Process completed successfully. ---")