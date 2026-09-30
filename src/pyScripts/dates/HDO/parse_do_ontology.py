import obonet
import networkx as nx
import pickle
import json

obo = snakemake.input.obo
ancestors_file = snakemake.output.ancestors
xrefs_file = snakemake.output.xrefs

print(f"--- Parsing DO ontology... ---\n")
print(f"--- Loading obo... ---\n")

do_graph = obonet.read_obo(obo)

ancestors_map = {
    doid: set(nx.descendants(do_graph, doid))
    for doid in do_graph.nodes
}
with open(ancestors_file, "wb") as f:
    pickle.dump(ancestors_map, f)

doid_to_mesh_omim = {}
for doid, data in do_graph.nodes(data=True):
    xrefs = data.get("xref", [])
    ids = [x for x in xrefs if x.startswith("MESH:") or x.startswith("OMIM:")]
    if ids:
        doid_to_mesh_omim[doid] = ids
with open(xrefs_file, "w") as f:
    json.dump(doid_to_mesh_omim, f)
    
print(f"--- Ancestors and xrefs files ready! ---")