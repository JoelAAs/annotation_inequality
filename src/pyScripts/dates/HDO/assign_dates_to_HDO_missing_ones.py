import networkx as nx
import pickle
import pronto

network_input = snakemake.input.network_with_dates
do = snakemake.input.ontology
network_output = snakemake.output.network_with_all_dates

print(f"Loading HDO network with missing dates...")
with open(network_input, 'rb') as f:
    G = pickle.load(f)

print("Loading Gene Ontology database from local file (this may take a few seconds)...")
do_onto = pronto.Ontology(do)

print(f"Propagating dates from HDO child terms up to parent terms...")
updated_nodes = 0
total_dates_injected = 0

# Loop through every gene in the network
for n, attr in G.nodes(data=True):
    annotations = attr.get('do_annotations', [])
    if not annotations:
        continue
        
    # Collect all the dates we DO have for this specific gene
    node_dates = {}
    for ann in annotations:
        date = ann.get('first_annotation_date')
        if date is not None:
            node_dates[ann['doid']] = int(date)
            
    if not node_dates:
        continue # If this gene has absolutely zero dates, we can't infer anything.
        
    # Time-travel up the tree!
    # For every term with a date, tell all of its ancestors about that date.
    inferred_dates = {}
    for doid, date in node_dates.items():
        if doid in do_onto:
            # .superclasses() gets the term itself PLUS all parents/grandparents up to the root
            for ancestor in do_onto[doid].superclasses():
                anc_id = ancestor.id
                
                if anc_id not in inferred_dates:
                    inferred_dates[anc_id] = date
                else:
                    # TRUE PATH RULE: The parent gets the OLDEST (minimum) date among its children
                    inferred_dates[anc_id] = min(inferred_dates[anc_id], date)
                    
    # Fill in the blanks in the graph's memory
    node_updated = False
    for ann in annotations:
        # If this term is missing a date...
        if ann.get('first_annotation_date') is None:
            # ...check if we managed to infer one from its children!
            inferred = inferred_dates.get(ann['doid'])
            if inferred is not None:
                ann['first_annotation_date'] = inferred
                node_updated = True
                total_dates_injected += 1
                
    if node_updated:
        updated_nodes += 1

print(f"Successfully injected {total_dates_injected} missing dates across {updated_nodes} nodes in the HDO network!")

print(f"Saving completed HDO network...")

with open(network_output, 'wb') as f:
    pickle.dump(G, f)

print(f"Missing dates added to HDO network!")