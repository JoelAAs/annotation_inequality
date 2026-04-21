import pandas as pd
import networkx as nx
import pickle

raw_network = snakemake.input.bp_network
annot_df = snakemake.input.complete_annotations
HDO_network = snakemake.output.bp_network_with_attributes

print(f"Adding HDO annotations to raw network...\n")

print(f"Loading raw data...")

with open(raw_network, 'rb') as f:
    G = pickle.load(f)
annotations = pd.read_csv(annot_df, sep = '\t')
annotations = annotations.dropna(subset = ['depth'])

print(f"Raw data loaded!\n")

print(f"Processing HDO network annotations...")

annotations['entrez_id'] = annotations['entrez_id'].astype(str)
annotation_map = annotations.groupby('entrez_id').apply(
    lambda x: x[['doid', 'depth']].to_dict('records'),
    include_groups = False
).to_dict()

nx.set_node_attributes(G, annotation_map, name = 'doids_with_depth')

print(f"Saving HDO network...")

with open(HDO_network, 'wb') as f:
    pickle.dump(G, f)

print(f"HDO network saved!\n")

print(f"Network with HDO annotations ready!\n")