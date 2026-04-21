import networkx as nx
import pandas as pd
import pickle

bp_publications_df = snakemake.input.bp_publications
bait_usage_df = snakemake.input.bait_usage
outputdf = snakemake.output.degree_frequencies
outputgraph = snakemake.output.bp_network

print(f"Processng bait prey publications raw network...\n")

print(f"Loading raw data...")

df = pd.read_parquet(bp_publications_df)
bait_usage = pd.read_csv(bait_usage_df, sep = '\t')

print(f"Raw data ready!\n")

print(f"Creating network...")

G = nx.from_pandas_edgelist(df, "entrez_id_bait",  "entrez_id_prey")

print(f"Network ready!\n")

print("-----")
print(f"Unique nodes = {G.number_of_nodes()}")
print(f"Unique edges = {G.number_of_edges()}")
print("-----\n")

print("Adding bait count attribute to the nodes...")

# Transform the bait usage df into a dict for faster computation
bait_map = {str(k): v for k, v in bait_usage.set_index('entrez_id_bait')['count'].to_dict().items()}
# Look inside the dictionary for each node in the graph, if the entrez id is found in the df
# return its count value, otherwise assign 0
bait_attr = {node: bait_map.get(str(node), 0) for node in G.nodes}

# Assign the bait count as attribute to each node
nx.set_node_attributes(G, bait_attr, name = 'bait_count')

print("Bait counts added to the nodes!\n")

print(f"Saving degree frequencies...")

degree_to_entrez = dict()
for node, degree in dict(G.degree()).items():
    if degree in degree_to_entrez:
        degree_to_entrez[degree].append(node)
    else:
        degree_to_entrez[degree] = [node,]

rows = []

for degree, genes in degree_to_entrez.items():
    for gene in genes:
        rows.append({'degree': degree, 'entrez_id': gene})

degree_df = pd.DataFrame(rows)
degree_df = degree_df.sort_values(by = 'degree', ascending = True)
degree_df.to_csv(outputdf, sep = '\t', index = False)

print(f"Degree frequencies saved!\n")

print(f"Saving network...")

with open(outputgraph, "wb") as f:
    pickle.dump(G, f)

print(f"Network saved!\n")

print(f"Bait prey publications raw network ready!\n")