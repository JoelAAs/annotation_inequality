import pandas as pd
import mygene
import pronto
import sys
import time

# 1. Setup Wildcards and Paths
# Accessing wildcards and files from the Snakemake object
aspect = snakemake.wildcards.aspect  # Expected values: 'BP', 'MF', or 'CC'
input_df_path = snakemake.input.bp_frequencies
obo_path = snakemake.input.go_obo
output_df = snakemake.output.annotation_df
output_list = snakemake.output.annotation_list

# --- ROOT MAPPING ---
# Maps your wildcards directly to the official GO Root IDs
ROOT_MAP = {
    'BP': 'GO:0008150',
    'MF': 'GO:0003674',
    'CC': 'GO:0005575'
}

target_root = ROOT_MAP.get(aspect)

if not target_root:
    print(f"ERROR: Wildcard '{aspect}' not recognized. Use BP, MF, or CC.")
    sys.exit(1)

print(f"--- Processing GO {aspect} (Target Root: {target_root}) ---")

# 2. Load Gene Universe
# Using parquet for speed as discussed previously
df = pd.read_parquet(input_df_path)
all_genes_list = df['entrez_id'].astype(str).unique().tolist()
all_genes_df = pd.DataFrame({'entrez_id': all_genes_list})

# 3. Load Ontology & Calculate Depths
print(f"Loading Ontology and calculating depths...")
ontology = pronto.Ontology(obo_path)
depths = {}

# Explicitly initialize the root to 0
if target_root in ontology:
    depths[target_root] = 0
else:
    print(f"ERROR: Root {target_root} not found in the OBO file!")
    sys.exit(1)

def get_depth_recursive(term):
    # Base case: if we already calculated it
    if term.id in depths:
        return depths[term.id]
    
    # Get direct parents
    parents = list(term.superclasses(distance=1, with_self=False))
    
    if not parents:
        # This belongs to a different branch of the GO tree
        depths[term.id] = -1
        return -1
    
    p_depths = [get_depth_recursive(p) for p in parents]
    valid_p_depths = [d for d in p_depths if d != -1]
    
    if not valid_p_depths:
        depths[term.id] = -1
        return -1
    
    # Depth = 1 + maximum distance to our target root (Longest Path)
    d = 1 + max(valid_p_depths)
    depths[term.id] = d
    return d

print("Starting depth calculation...")
for term in ontology.terms():
    get_depth_recursive(term)

# Cleanup: Filter out terms that don't lead to our root
final_depth_map = {k: v for k, v in depths.items() if v != -1}
# Ensure the root is definitely in there as 0
final_depth_map[target_root] = 0

print(f"Depths calculated for {len(final_depth_map)} terms. Root {target_root} is Depth {final_depth_map.get(target_root)}")

# 4. Query MyGene with Retry Logic
mg = mygene.MyGeneInfo()
results = None
for i in range(3):
    try:
        print(f"Querying MyGene (Attempt {i+1})...")
        results = mg.querymany(
            all_genes_list,
            scopes='entrezgene',
            fields=f'go.{aspect}',
            species='human',
            as_dataframe=True
        )
        break
    except Exception as e:
        print(f"Server error: {e}. Retrying in 10 seconds...")
        time.sleep(10)

if results is None:
    sys.exit(1)

# 5. Process Results
rows = []
go_col = f'go.{aspect}'
if go_col not in results.columns:
    results[go_col] = None

root_term_obj = ontology[target_root]

print("Climbing the hierarchy for gene annotations...")
for entrez_id, data in results.iterrows():
    go_entries = data.get(go_col)
    
    # Skip genes with no annotations
    if go_entries is None or (isinstance(go_entries, float) and pd.isna(go_entries)):
        continue

    entries_to_process = go_entries if isinstance(go_entries, list) else [go_entries]
    gene_annotations = set()
    
    # --- FORCE ROOT (DEPTH 0) ---
    # If the gene has any valid GO data, we force the Root into the set
    gene_annotations.add((root_term_obj.id, root_term_obj.name, 0))

    for entry in entries_to_process:
        direct_id = entry.get('id')
        if direct_id in final_depth_map:
            target_term = ontology[direct_id]
            # Add the term and all its ancestors that are in our target branch
            for ancestor in target_term.superclasses():
                if ancestor.id in final_depth_map:
                    gene_annotations.add((ancestor.id, ancestor.name, final_depth_map[ancestor.id]))
    
    for anc_id, anc_name, anc_depth in gene_annotations:
        rows.append({
            'entrez_id': entrez_id, 
            'go_id': anc_id, 
            'go_term_name': anc_name, 
            'depth': anc_depth
        })

# 6. Save Outputs
if not rows:
    final_df = all_genes_df.assign(go_id=None, go_term_name=None, depth=None)
else:
    raw_annot_df = pd.DataFrame(rows).drop_duplicates()
    final_df = pd.merge(all_genes_df, raw_annot_df, on='entrez_id', how='left')

# Output 1: Full mapping for matrix generation
final_df.to_csv(output_df, sep='\t', index=False)

# Output 2: Unique list for depth summary calculation
final_list = final_df[['go_id', 'go_term_name', 'depth']].dropna().drop_duplicates()
final_list.to_csv(output_list, sep='\t', index=False)

print(f"Successfully processed {aspect}. Check {output_list} for Depth 0.")