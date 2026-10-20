import networkx as nx
import pickle
import pronto
import pandas as pd

file_network = snakemake.input.network_with_dates
file_ontology = snakemake.input.ontology
file_bp_pubs = snakemake.input.bp_publications
file_output = snakemake.output.final_network

# Load data
print("--- [INFO] Loading network... ---")
with open(file_network, 'rb') as f:
    G = pickle.load(f)

print("--- [INFO] Loading Human Disease Ontology (HDO)... ---")
ontology = pronto.Ontology(file_ontology)

# Dates propagation
print("--- [INFO] Mapping ancestor relationships... ---")
ancestor_map = {}
for term in ontology.terms():
    ancestor_map[term.id] = {t.id for t in term.superclasses()}

print("--- [INFO] Propagating dates from descendants to ancestors... ---")
none_resolved = 0
dates_improved = 0
remaining_nones = 0

for node_id, node_data in G.nodes(data=True):
    if 'do_annotations' in node_data:
        annots = node_data['do_annotations']
        original_dates = {a['doid']: a['first_annotation_date'] for a in annots}
        new_dates = dict(original_dates)
        
        for doid, date in original_dates.items():
            if date is not None and doid in ancestor_map:
                for ancestor_doid in ancestor_map[doid]:
                    if ancestor_doid in new_dates:
                        curr_anc_date = new_dates[ancestor_doid]
                        if curr_anc_date is None or date < curr_anc_date:
                            new_dates[ancestor_doid] = date
                            
        for a in annots:
            doid = a['doid']
            old_val = a['first_annotation_date']
            new_val = new_dates[doid]
            
            if old_val is None and new_val is not None:
                none_resolved += 1
            elif old_val is not None and new_val is not None and new_val < old_val:
                dates_improved += 1
                
            a['first_annotation_date'] = new_val
            
            if new_val is None:
                remaining_nones += 1

print("\n--- [INFO] Temporal Propagation Results ---")
print(f"    > Missing dates resolved ('None' -> Date): {none_resolved}")
print(f"    > Ancestor dates pushed back (Improved): {dates_improved}")
print(f"    > Annotations still left with 'None': {remaining_nones}")

# Edge dates injection
print("\n--- [INFO] Injecting PPI interaction discovery dates into edges... ---")
df_edges = pd.read_parquet(file_bp_pubs)

COL_NODE1 = 'entrez_id_bait'
COL_NODE2 = 'entrez_id_prey'
COL_DATE  = 'date_created'

edge_date_lookup = {}
print(f"    > Scanning {len(df_edges)} publication records to find first discovery dates...")

# Extract the earliest date for each interacting pair
for _, row in df_edges.iterrows():
    u = str(row[COL_NODE1])
    v = str(row[COL_NODE2])
    
    # Skip invalid rows
    if pd.isna(row[COL_DATE]):
        continue
        
    # The 'date_created' format in your DataFrame is already YYYYMMDD as an integer
    date = int(row[COL_DATE]) 
    
    # Sort to ensure (A, B) and (B, A) are treated as the same undirected edge
    edge_key = tuple(sorted([u, v]))
    
    # Keep strictly the OLDEST (minimum) date for the interaction
    if edge_key not in edge_date_lookup or date < edge_date_lookup[edge_key]:
        edge_date_lookup[edge_key] = date

edges_updated = 0
edges_missing_dates = 0

# Apply the earliest dates to the NetworkX graph edges
for u, v in G.edges():
    edge_key = tuple(sorted([str(u), str(v)]))
    if edge_key in edge_date_lookup:
        G[u][v]['discovery_date'] = edge_date_lookup[edge_key]
        edges_updated += 1
    else:
        edges_missing_dates += 1

print(f"    > Edges successfully updated with discovery dates: {edges_updated}")
if edges_missing_dates > 0:
    print(f"    > WARNING: {edges_missing_dates} edges in the network have no matching date in the intact file.")

# Save final network
print("\n--- [INFO] Saving the final time-resolved network... ---")
with open(file_output, 'wb') as f:
    pickle.dump(G, f)

print("--- [INFO] Process completed successfully. ---")