import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import pronto

cutoff = snakemake.wildcards.cutoff
coefficients_df = snakemake.input.all_coefficients
ontology = snakemake.input.ontology
outputimage = snakemake.output.treemap
outputimage_html = snakemake.output.treemap_html

print(f"Creating treemap for HDO cutoff {cutoff}...\n")

print(f"Loading data for HDO cutoff {cutoff}...")

coefficients = pd.read_csv(coefficients_df, sep = '\t')
coeff_map = dict(zip(coefficients['HDO_doid'], coefficients['Coefficient']))

print(f"Data for HDO cutoff {cutoff} loaded!\n")

print("Loading HDO ontology from local file...")

hdo = pronto.Ontology(ontology)

print("HDO ontology loaded!\n")

print(f"Tracing HDO cutoff {cutoff} lineage and coefficients...")

paths = {}

# Root
root_id = 'DOID:4'
root_val = coeff_map.get(root_id, 0.0)
paths[root_id] = {
    'label': 'Disease',
    'parent': '',
    'coefficient': root_val,
    'value': 15
}

for doid in coeff_map.keys():
    if doid in hdo:
        # Get lineage: pronto superclasses usually go [Leaf -> ... -> Root]
        lineage = [ancestor.id for ancestor in hdo[doid].superclasses()]
        lineage = [n for n in lineage if str(n).startswith("DOID:")]
        
        # ENSURE DOID:4 is at the START
        if root_id in lineage:
            lineage.remove(root_id)
        lineage = [root_id] + list(reversed(lineage)) 
        # Now lineage is guaranteed to be [DOID:4, Category, ..., doid]

        for i in range(len(lineage)):
            node_id = lineage[i]
            u_path = "/".join(lineage[:i+1])
            u_parent = "/".join(lineage[:i]) if i > 0 else ''
            
            # Root is the only one with an empty parent
            if node_id == root_id:
                continue

            if u_path not in paths:
                name = hdo[node_id].name if node_id in hdo else node_id
                val = coeff_map.get(node_id, 0.0)
                # i is now the correct depth from the root
                paths[u_path] = {
                    'label': name,
                    'parent': u_parent,
                    'coefficient': val,
                    'value': max(15 - i, 1)
                }

# Final Conversion
tree_df = pd.DataFrame.from_dict(paths, orient='index').reset_index()
tree_df.rename(columns={'index': 'id'}, inplace=True)

# --- SAFETY CHECKS ---
# A. Remove any exact duplicates
tree_df = tree_df.drop_duplicates(subset=['id'])

# B. Ensure 'value' is numeric and > 0 (Plotly crashes on 0 or NaN values)
tree_df['value'] = pd.to_numeric(tree_df['value'], errors='coerce').fillna(1)
tree_df.loc[tree_df['value'] <= 0, 'value'] = 1

# C. Prevent self-parenting (id must not equal parent)
tree_df = tree_df[tree_df['id'] != tree_df['parent']]

# D. Final check: if 'id' is DOID:4, parent MUST be ''
tree_df.loc[tree_df['id'] == 'DOID:4', 'parent'] = ''
# -------------------------

all_ids = set(tree_df['id'])
missing_parents = [p for p in tree_df['parent'] if p != '' and p not in all_ids]
if missing_parents:
    print(f"CRITICAL ERROR: These parents are missing from the ID list: {set(missing_parents)}")

print(f"HDO cutoff {cutoff} lineage and coefficients traced!\n")

print(f"Creating HDO cutoff {cutoff} treemap visualization...")

dynamic_max = max(tree_df['coefficient'].abs().max(), 0.001)

fig = go.Figure(go.Treemap(
    ids = tree_df['id'],
    labels = tree_df['label'],
    parents = tree_df['parent'],
    values = tree_df['value'],
    branchvalues = "remainder",
    marker = dict(
        colors = tree_df['coefficient'],
        colorscale = ["#d63031", "#ffffff", "#27ae60"],
        cmid = 0.0,
        cmin = -dynamic_max,
        cmax = dynamic_max,
        showscale = True
    ),
    hovertemplate = "<b>%{label}</b><br>ID: %{id}<br>Coeff: %{color:.4f}<extra></extra>"
))

fig.update_layout(margin=dict(t=50, l=10, r=10, b=10), 
                  uniformtext=dict(minsize=10, mode='hide')
)

print(f"HDO cutoff {cutoff} treemap visualization done!\n")

print(f"Saving HDO cutoff {cutoff} treemap into pdf and HTML...")

fig.write_image(outputimage, engine = 'kaleido')
fig.write_html(outputimage_html)

print(f"HDO cutoff {cutoff} treemap in pdf and HTML saved!\n")

print(f"HDO cutoff {cutoff} treemap ready!\n")