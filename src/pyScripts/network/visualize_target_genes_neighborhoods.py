import pandas as pd
import networkx as nx
import numpy as np
import pickle
from pyvis.network import Network

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
observed_df = snakemake.input.observed
network = snakemake.input.bp_network_with_attributes
outputnetwork = snakemake.output.interactive_network

print(f"Creating interactive network for HDO top coefficients depth {depth} cutoff {cutoff}...")

print(f"Loading data for HDO top coefficients depth {depth} cutoff {cutoff}...")

with open(observed_df, 'rb') as f:
    obs_df = pickle.load(f)
with open(network, 'rb') as f:
    G = pickle.load(f)

print(f"Extracting subgraph for HDO top coefficients depth {depth} cutoff {cutoff}...")

# Extract a manageable Subgraph (focusing on top predictors)
target_genes = []

for i in range(len(obs_df)):
    # Get the top 15 genes for each HDO annotation from your pickle
    genes = [str(g['gene']) for g in obs_df.iloc[i]['genes_with_sums'][:15]]
    target_genes.extend(genes)

target_genes = list(set(target_genes)) # Remove duplicates

# Get neighbors to see the "neighborhood" context
nodes_to_plot = set(target_genes)
for node in target_genes:
    if node in G:
        nodes_to_plot.update(G.neighbors(node))

sub = G.subgraph(nodes_to_plot)

print(f"Initializing Pyvis and mapping NetworkX to it for HDO top coefficients depth {depth} cutoff {cutoff}...")

# Initialize Pyvis
net = Network(height="850px", width="100%", bgcolor="#ffffff", font_color="black")

# Map NetworkX to Pyvis with Custom Styling
for node_id in sub.nodes():
    # Get metadata from original Graph
    b_count = G.nodes[node_id].get('bait_count', 0)
    deg = G.degree(node_id)
    
    # Visual Properties
    # Size scales with the square root of bait count to keep it visible
    size = np.sqrt(b_count + 1) * 3 
    
    # Color: Green for disease predictors, Gray for neighbors
    color = '#66c2a5' if node_id in target_genes else '#ced4da'
    
    # Label and Hover-over title
    label = str(node_id)
    hover_title = f"Gene ID: {node_id}\nBait Count: {b_count}\nDegree: {deg}"
    
    net.add_node(node_id, 
                 label=label, 
                 title=hover_title, 
                 size=size, 
                 color=color)

# Add edges from the subgraph
for source, target in sub.edges():
    net.add_edge(source, target, color='#e9ecef', width=0.5)

# Physics and Interaction Settings
# This makes the "hubs" naturally float to the center
net.barnes_hut(gravity=-8000, central_gravity=0.3, spring_length=100)
net.show_buttons(filter_=['physics']) # Allows you to tweak live in browser

print(f"Saving HTML for HDO top coefficients depth {depth} cutoff {cutoff}...")

# Save to HTML
net.save_graph(outputnetwork)

print(f"Interactive network for HDO top coefficients depth {depth} cutoff {cutoff} ready!")