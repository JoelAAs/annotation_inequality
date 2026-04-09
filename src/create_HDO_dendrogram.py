import pandas as pd
import matplotlib.pyplot as plt
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

colors = ["#ffb3b3", "#ffffff", "#b3ffb3"]
custom_cmap = mcolors.LinearSegmentedColormap.from_list("LightRdGn", colors)
norm = mcolors.Normalize(vmin=-1, vmax=1)

dendrogram_style = {
    'rankdir': 'TB',    
    'splines': 'ortho', 
    'nodesep': '0.8',   
    'ranksep': '1.0',   
    'fontname': 'Arial'
}

graph = pydot.Dot(graph_type='digraph', **dendrogram_style)

def get_node_attributes(doid):
    if doid not in coeff_map:
        return {'fillcolor': '#f0f0f0', 'label': doid, 'style': 'filled'}
    
    val = coeff_map[doid]
    color_hex = mcolors.to_hex(custom_cmap(norm(val)))
    
    # Get term name
    name = hdo[doid].name if doid in hdo else "Unknown"
    label = f"{doid}\n{name}\n({val:.4f})"
    
    return {
        'fillcolor': color_hex, 
        'label': label, 
        'style': 'filled',
        'shape': 'box',
        'fontname': 'Helvetica'
    }

nodes_to_plot = set()
edges_to_plot = set()

## TODO finish this

print(f"Building HDO cutoff {cutoff} DAG lineage...")

for doid in coeff_map.keys():
    if doid in hdo:
        nodes_to_plot.add(doid)
        # Add lineage
        for ancestor in hdo[doid].superclasses():
            nodes_to_plot.add(ancestor.id)
            # Add relationships
            for parent in hdo[ancestor.id].superclasses(distance = 1):
                edges_to_plot.add((parent.id, ancestor.id))

print(f"HDO cutoff {cutoff} DAG lineage computed!\n")

print(f"Creating HDO cutoff {cutoff} visualization...")

for node_id in nodes_to_plot:
    # Only plot if node is part of the DOID:4 (Disease) hierarchy to avoid floating nodes
    attrs = get_node_attributes(node_id)
    graph.add_node(pydot.Node(node_id, **attrs))

for edge in edges_to_plot:
    if edge[0] in nodes_to_plot and edge[1] in nodes_to_plot:
        graph.add_edge(pydot.Edge(edge[0], edge[1]))

print(f"HDO cutoff {cutoff} visualization done!\n")

print(f"Saving HDO cutoff {cutoff} to png...")

graph.write_png(outputimage)

print(f"HDO cutoff {cutoff} png saved!\n")

print(f"HDO cutoff {cutoff} dendrogram ready!\n")