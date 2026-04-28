import pandas as pd
import pickle
import networkx as nx

aspect = snakemake.wildcards.aspect
network_file = snakemake.input.network
dates_file = snakemake.input.dates
mapping_file = snakemake.input.mapping 
output_file = snakemake.output.network_with_dates

# Build the HGNC Reverse-Lookup Dictionary 
print(f"[STATUS] Loading HGNC mapping for GO {aspect} ID resolution...")
df_hgnc = pd.read_csv(mapping_file, sep="\t", dtype=str, low_memory=False)
entrez_col = "entrez_id" if "entrez_id" in df_hgnc.columns else "ncbi_id"

name_to_entrez = {}
for _, row in df_hgnc.iterrows():
    eid = str(row[entrez_col]).split('.')[0]
    if pd.isna(eid) or eid == "nan": continue
    
    if pd.notna(row["symbol"]):
        name_to_entrez[row["symbol"].upper()] = eid
    if pd.notna(row["uniprot_ids"]):
        for u in str(row["uniprot_ids"]).replace('"', '').split('|'):
            name_to_entrez[u.strip().upper()] = eid

# Build the Date Lookup Dictionary
print(f"[STATUS] Building GO {aspect} date lookup table...")
df_dates = pd.read_csv(dates_file, sep="\t")
date_lookup = {}
for _, row in df_dates.iterrows():
    gene_alias = str(row["Gene"]).upper()
    go_id = str(row["GO_ID"])
    date = int(row["First_Date_Annotated"])
    
    # Resolve the alias back to its Entrez ID
    entrez_id = name_to_entrez.get(gene_alias)
    if entrez_id:
        # Keep the oldest date if multiple aliases map to the same Entrez ID
        key = (entrez_id, go_id)
        if key not in date_lookup or date < date_lookup[key]:
            date_lookup[key] = date

# Load Network and Build Comprehensive Data Structure
with open(network_file, "rb") as f:
    G = pickle.load(f)

print(f"[STATUS] Integrating comprehensive dates into {G.number_of_nodes()} GO {aspect} nodes...")
for node_id, attr in G.nodes(data=True):
    eid = str(node_id)
    
    # Grab the original list of dictionaries
    original_annotations = attr.get("go_ids_with_depth", []) 
    
    # Build a new, clean list to hold our unified data
    comprehensive_annotations = []
    
    for item in original_annotations:
        go_id = item.get('go_id')
        if go_id:
            # Create a unified dictionary for this specific GO term
            comprehensive_dict = {
                'go_id': go_id,
                'depth': item.get('depth'),
                'first_annotation_date': date_lookup.get((eid, str(go_id)))
            }
            comprehensive_annotations.append(comprehensive_dict)

    # Assign the unified list to a brand new, clearly named attribute
    attr["go_annotations"] = comprehensive_annotations
    
    # Cleanup --> Delete the old keys so your network nodes stay perfectly clean
    if "go_ids_with_depth" in attr:
        del attr["go_ids_with_depth"]
    if "annotation_dates" in attr:
        del attr["annotation_dates"]

# Save 
with open(output_file, "wb") as f:
    pickle.dump(G, f)

print(f"[STATUS] Successfully saved fully integrated GO {aspect} network!")