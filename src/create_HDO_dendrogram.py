import pandas as pd
import matplotlib.colors as mcolors
import pronto
import pydot

cutoff = snakemake.wildcards.cutoff
coefficients_df = snakemake.input.all_coefficients
ontology = snakemake.input.ontology
outputimage = snakemake.output.dendrogram

print(f"Creating HDO dendrogram for cutoff {cutoff}...\n")

print(f"Loading data for HDO cutoff {cutoff}...")

coefficients = pd.read_csv(coefficients_df, sep = '\t')
coeff_map = dict(zip(coefficients['HDO_doid'], coefficients['Coefficient']))

print(f"Data for HDO cutoff {cutoff} loaded!\n")

print("Loading HDO ontology from local file...")

hdo = pronto.Ontology(ontology)

print("HDO ontology loaded!\n")

# Setup color mapping (
#   Red -> low
#   White -> 0
#   Green -> high 
# )

colors = ["#d63031", "#ffffff", "#27ae60"]
custom_cmap = mcolors.LinearSegmentedColormap.from_list("VibrantRdGn", colors)
max_val = max(abs(coefficients['Coefficient'].min()), abs(coefficients['Coefficient'].max()), 0.01)
norm = mcolors.Normalize(vmin=-max_val, vmax=max_val)

dendrogram_style = {
    'rankdir': 'TB',
    'splines': 'ortho',
    'nodesep': '1.0',       
    'ranksep': '1.5',    
    'fontname': 'Arial',
    'fontsize': '16',    
    'dpi': '300'         
}

graph = pydot.Dot(graph_type='digraph', **dendrogram_style)

print(f"Building HDO cutoff {cutoff} Hierarchy...")

nodes_to_plot = set()
edges_to_plot = set()

for doid in coeff_map.keys():
    if doid in hdo:
        nodes_to_plot.add(doid)
        # Trace the lineage to ensure everything is connected to the root
        for ancestor in hdo[doid].superclasses():
            nodes_to_plot.add(ancestor.id)
            # Find the immediate parent to draw the tree branch
            for parent in hdo[ancestor.id].superclasses(distance = 1):
                if parent.id != ancestor.id:
                    edges_to_plot.add((parent.id, ancestor.id))

print(f"HDO cutoff {cutoff} Hierarchy computed!\n")

def clean_id(doid):
    return doid.replace(':', '_')

print(f"Adding HDO cutoff {cutoff} nodes to the graph...")

for node_id in nodes_to_plot:
    safe_id = clean_id(node_id)
    val = coeff_map.get(node_id, 0)
    color = mcolors.to_hex(custom_cmap(norm(val)))
    name = hdo[node_id].name if node_id in hdo else node_id
    label = f'"{node_id}\n{name}\n({val:.3f})"'

    node = pydot.Node(
        safe_id, 
        label = label, 
        fillcolor = color, 
        style = 'filled', 
        shape = 'box', 
        fontsize = '14', 
        margin = '0.2', 
        width = '1.5'
    )
    graph.add_node(node)

print(f"HDO cutoff {cutoff} nodes added to the graph!\n")

print(f"Adding HDO cutoff {cutoff} edges to the graph...")

for parent_id, child_id in edges_to_plot:
    graph.add_edge(pydot.Edge(clean_id(parent_id), clean_id(child_id)))

print(f"HDO cutoff {cutoff} edges added to the graph!\n")

print(f"Saving HDO cutoff {cutoff} to pdf...")

graph.write_pdf(outputimage)

print(f"HDO cutoff {cutoff} pdf saved!\n")

print(f"HDO cutoff {cutoff} dendrogram ready!\n")