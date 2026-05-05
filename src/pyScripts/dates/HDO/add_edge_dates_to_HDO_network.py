import pandas as pd
import networkx as nx
import pickle

network = snakemake.input.network_with_all_dates
bp_publications_df = snakemake.input.bp_publications
outputnetwork = snakemake.output.final_network

print(f"Adding edge dates to HDO network...")

print(f"Loading HDO network data...")

with open(network, 'rb') as f:
    G = pickle.load(f)
df = pd.read_parquet(bp_publications_df)

print(f"Processing HDO edge dates...")

edge_dates = {}

for _, row in df.iterrows():
    bait = str(row['entrez_id_bait'])
    prey = str(row['entrez_id_prey'])
    date = int(row['date_created'])

    # Check if the edge already exist in the dictionary
    if (bait, prey) in edge_dates:
        oldest_date = min(edge_dates[(bait, prey)], date)
        edge_dates[(bait, prey)] = oldest_date
        edge_dates[(prey, bait)] = oldest_date
    else:
        edge_dates[(bait, prey)] = date
        edge_dates[(prey, bait)] = date    

print(f"Injecting HDO edge dates into {G.number_of_edges()} edges...")
updated_count = 0

for u, v, attr in G.edges(data = True):
    discovery_date = edge_dates.get((str(u), str(v)))

    if discovery_date is not None:
        attr['discovery_date'] = discovery_date
        updated_count += 1
    else:
        attr['discovery_date'] = None

print(f"Successfully injected {updated_count} into HDO network!")

print(f"Saving HDO network with edge dates...")

with open(outputnetwork, 'wb') as f:
    pickle.dump(G, f)

print(f"HDO edge dates added!")