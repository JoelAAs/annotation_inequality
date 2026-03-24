import pandas as pd
import obonet
import networkx as nx

input_df = snakemake.input.annot_df
output_df = snakemake.output.annot_df_depth

df = pd.read_csv(input_df, sep  ='\t')

url = "http://purl.obolibrary.org/obo/doid.obo"
graph = obonet.read_obo(url)
root = 'DOID:4'
graph_rev = graph.reverse(copy = True)
depths = nx.shortest_path_length(graph_rev, source = root)
print(list(depths.keys())[:10]) 

df['depth'] = df['doid'].map(lambda x: depths.get(x, None))

df.to_csv(output_df, sep = '\t', index = False)