import pandas as pd
import networkx as nx
import pickle

aspect = snakemake.wildcards.aspect
network_file = snakemake.input.final_network
top_annot_file = snakemake.input.nodes_with_top_5_annotations_pickle
output_csv = snakemake.output.top_5_annotations_dates_df
output_pkl = snakemake.output.top_5_annotations_dates_pkl

print(f"Loading GO {aspect} network and top annotations...")

with open(network_file, 'rb') as f:
    G = pickle.load(f)

with open(top_annot_file, 'rb') as f:
    top_annot_df = pickle.load(f)

results = []

for index, row in top_annot_df.iterrows():
    go_id = row['GO_id']
    genes = row['annotated_genes']
    
    dates_list = []
    
    for gene in genes:
        if gene in G.nodes:
            annotations = G.nodes[gene].get('go_annotations', [])
            
            for ann in annotations:
                if ann.get('go_id') == go_id:
                    date = ann.get('first_annotation_date')
                    if date is not None:
                        dates_list.append(date)
                    break 

    results.append({
        'GO_id': go_id,
        'annotation_dates': dates_list
    })

new_df = pd.DataFrame(results)
new_df.to_csv(output_csv, sep = '\t', index=False)

with open(output_pkl, 'wb') as f:
    pickle.dump(new_df, f)

print(f"Done computing top 5 GO {aspect} annotation dates!")