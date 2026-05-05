import pandas as pd
import networkx as nx
import pickle

input_network = snakemake.input.network
dates_df = snakemake.input.dates
outputnetwork = snakemake.output.network_with_dates

print(f"Processing data addition to HDO network nodes...")

print(f"Loading HDO network and dates...")

with open(input_network, 'rb') as f:
    G = pickle.load(f)

df = pd.read_csv(dates_df, sep = '\t')

print(f"Building HDO fast date lookup dictionary...")

df = df.dropna(subset = ['entrez_targetId', 'standardized_DOID', 'first_publication_date'])

date_lookup = {}

for row in df.itertuples():
    entrezid = str(int(row.entrez_targetId))
    doid = str(row.standardized_DOID)
    date = int(row.first_publication_date)

    # In case there are multiple publications for the same gene-disease pair, keep the oldest
    key = (entrezid, doid)
    if key not in date_lookup or date < date_lookup[key]:
        date_lookup[key] = date

print(f"Adding dates to HDO network nodes...")

for node_id, attr in G.nodes(data = True):
    entrezid = str(node_id)

    # Take the original list of dictionaries and create a new one
    original_annotations = attr.get('doids_with_depth', [])
    new_annotations = []

    for item in original_annotations:
        doid = item.get('doid')
        if doid:
            new_dict = {
                'doid': doid,
                'depth': item.get('depth'),
                'first_annotation_date': date_lookup.get((entrezid, str(doid)))
            }
            new_annotations.append(new_dict)

    attr['do_annotations'] = new_annotations

    # Delete the old attribute
    if 'doids_with_depth' in attr:
        del attr['doids_with_depth']

print(f"Saving HDO network with dates...")

with open(outputnetwork, 'wb') as f:
    pickle.dump(G, f)

print(f"Generating HDO Network Dates Injection Quality Report...")

# Initialize our counters
total_annotations = 0
dates_found = 0
dates_missing = 0

# Loop through the updated network
for node_id, attr in G.nodes(data=True):
    # Safely grab the new list we just created
    annotations = attr.get('do_annotations', [])
    
    for item in annotations:
        total_annotations += 1
        
        # Check if the date is actually there
        if item.get('first_annotation_date') is not None:
            dates_found += 1
        else:
            dates_missing += 1

# Print the final statistics
print(f"--- HDO INJECTION REPORT ---")
print(f"Total DO Annotations in Network: {total_annotations}")
print(f"Successfully DO injected dates:     {dates_found}")
print(f"Missing DO dates (assigned None):   {dates_missing}")
print(f"----------------------------")

print(f"HDO network nodes dates added!")