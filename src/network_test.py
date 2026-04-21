import networkx as nx
import pandas as pd

bait_prey = pd.read_parquet("work_folder/data/intact/bait_prey_publications.pq")

G = nx.from_pandas_edgelist(bait_prey, "entrez_id_bait"  "entrez_id_prey")

degree_to_entrez = dict()
for node, degree in dict(G.degree()).items():
    if degree in degree_to_entrez:
        degree_to_entrez[degree].append(node)
    else:
        degree_to_entrez[degree] =[node,]
