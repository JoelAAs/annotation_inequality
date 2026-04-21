import pandas as pd
import networkx as nx
import pickle

aspect = snakemake.wildcards.aspect
raw_network = snakemake.input.bp_network
annot_df = snakemake.input.complete_annotations
GO_network = snakemake.output.bp_network_with_attributes

print(f"Adding GO {aspect} annotations to raw network...\n")

print(f"Loading raw data...")

with open(raw_network, 'rb') as f:
    G = pickle.load(f)
annotations = pd.read_csv(annot_df, sep = '\t')
annotations = annotations.dropna(subset = ['depth'])

print(f"Raw data loaded!\n")

print(f"Processing GO {aspect} network annotations...")

annotations['entrez_id'] = annotations['entrez_id'].astype(str)
annotation_map = annotations.groupby('entrez_id').apply(
    lambda x: x[['go_id', 'depth']].to_dict('records'),
    include_groups = False
).to_dict()

nx.set_node_attributes(G, annotation_map, name = 'go_ids_with_depth')

print(f"Saving GO {aspect} network...")

with open(GO_network, 'wb') as f:
    pickle.dump(G, f)

print(f"GO {aspect} network saved!\n")

print(f"Network with GO {aspect} annotations ready!\n")