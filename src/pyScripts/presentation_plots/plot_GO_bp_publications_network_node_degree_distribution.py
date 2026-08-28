import networkx as nx
import matplotlib.pyplot as plt
import collections
import pickle

# --- Snakemake I/O ---
input_network = snakemake.input.bp_network
output_plot = snakemake.output.plot_file

print("Loading bait-prey network...")
with open(input_network, "rb") as f:
    G = pickle.load(f)

# Extract degrees (exclude 0-degree isolated nodes)
degrees = [d for n, d in G.degree() if d > 0]
total_nodes = G.number_of_nodes()
total_edges = G.number_of_edges()

# ==========================================
# 1. CALCULATE EXACT FREQUENCIES
# ==========================================
degree_counts = collections.Counter(degrees)
# Sort by degree (X axis)
x_raw, y_raw = zip(*sorted(degree_counts.items()))

# ==========================================
# 2. DRAW THE PLOT (LINEAR SCALE)
# ==========================================
print("Drawing linear presentation plot...")
plt.figure(figsize=(14, 8))

# Layer 1 (Bottom): The slightly thicker line running UNDER the dots
plt.plot(
    x_raw, y_raw, 
    linestyle='-', color='#0571b0', linewidth=2.5, 
    zorder=1
)

# Layer 2 (Top): The dots with a very thin black border OVER the line
plt.plot(
    x_raw, y_raw, 
    linestyle='None', 
    marker='o', markersize=5, 
    color='#5ab4e5',  # Slightly lighter blue for the dots
    markeredgecolor='black', markeredgewidth=0.5,
    label='Node Count',
    zorder=2
)

# Title and labels scaled massively up for presentation slides
plt.title("Bait-Prey Publications Network Node Degree Distribution", fontsize=28, fontweight='bold', pad=20)
plt.xlabel("Degree", fontsize=22, labelpad=15)
plt.ylabel("Number of Nodes", fontsize=22, labelpad=15)

# Enlarge tick marks for readability from the back of the room
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Standard linear grid (pushed behind the data using zorder=0)
plt.grid(True, linestyle="--", alpha=0.5, zorder=0)

# Add statistics box with larger font
stats_text = f"Total Nodes: {total_nodes}\nTotal Edges: {total_edges}"
plt.text(
    0.95, 0.95, stats_text, 
    transform=plt.gca().transAxes,
    fontsize=20, verticalalignment='top', horizontalalignment='right',
    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#cccccc')
)

# Enlarge the legend and place it safely under the stats box
plt.legend(fontsize=18, loc='upper right', bbox_to_anchor=(0.95, 0.78))
plt.tight_layout()

# Save the plot
print(f"Saving presentation plot to {output_plot}...")
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

print("Plot successfully generated!")