import pandas as pd
import matplotlib.colors as mcolors
import pronto
import pydot

aspect = snakemake.wildcards.aspect
cutoff = snakemake.wildcards.cutoff
coefficients_df = snakemake.input.all_coefficients
ontology = snakemake.input.ontology
outputimage = snakemake.output.dendrogram

print(f"Creating GO {aspect} dendrogram for cutoff {cutoff}...\n")

print(f"Loading data for GO {aspect} cutoff {cutoff}...")

coefficients = pd.read_csv(coefficients_df, sep = '\t')
coeff_map = dict(zip(coefficients['GO_id'], coefficients['Coefficient']))

print(f"Data for GO {aspect} cutoff {cutoff} loaded!\n")

print("Loading GO ontology from local file...")

go = pronto.Ontology(ontology)

print("GO ontology loaded!\n")

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

print(f"Building GO {aspect} cutoff {cutoff} Hierarchy...")

nodes_to_plot = set()
edges_to_plot = set()

for goid in coeff_map.keys():
    if goid in go:
        nodes_to_plot.add(goid)
        # Trace the lineage to ensure everything is connected to the root
        for ancestor in go[goid].superclasses():
            nodes_to_plot.add(ancestor.id)
            # Find the immediate parent to draw the tree branch
            for parent in go[ancestor.id].superclasses(distance = 1):
                if parent.id != ancestor.id:
                    edges_to_plot.add((parent.id, ancestor.id))

print(f"GO {aspect} cutoff {cutoff} Hierarchy computed!\n")

def clean_id(goid):
    return goid.replace(':', '_')

print(f"Adding GO {aspect} cutoff {cutoff} nodes to the graph...")

for node_id in nodes_to_plot:
    safe_id = clean_id(node_id)
    val = coeff_map.get(node_id, 0)
    color = mcolors.to_hex(custom_cmap(norm(val)))
    name = go[node_id].name if node_id in go else node_id
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

print(f"GO {aspect} cutoff {cutoff} nodes added to the graph!\n")

print(f"Adding GO {aspect} cutoff {cutoff} edges to the graph...")

for parent_id, child_id in edges_to_plot:
    graph.add_edge(pydot.Edge(clean_id(parent_id), clean_id(child_id)))

print(f"GO {aspect} cutoff {cutoff} edges added to the graph!\n")

print(f"Saving GO {aspect} cutoff {cutoff} to pdf...")

graph.write_pdf(outputimage)

print(f"GO {aspect} cutoff {cutoff} pdf saved!\n")

print(f"GO {aspect} cutoff {cutoff} dendrogram ready!\n")