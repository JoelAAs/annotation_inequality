import pickle
import pandas as pd 
import random
import networkx as nx

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
inputnetwork = snakemake.input.bp_network_with_attributes
coefficients_df = snakemake.input.coefficients
baseline_df = snakemake.output.baseline_bait_count_sums_df
baseline_pickle = snakemake.output.baseline_bait_count_sums_pickle

print(f"Processing baseline for top GO {aspect} depth {depth} cutoff {cutoff} coefficients...")

print(f"Loading GO {aspect} depth {depth} cutoff {cutoff} data...")

with open(inputnetwork, 'rb') as f:
    G = pickle.load(f)
coefficients = pd.read_csv(coefficients_df, sep = '\t')

print(f"Computing baseline neighbors bait sum counts for GO {aspect} depth {depth} cutoff {cutoff}...")

# Fid top coefficients
coefficients['abs_coeff'] = coefficients['Coefficient'].abs()
top_5_go = coefficients.nlargest(5, 'abs_coeff')['GO_id'].tolist()

# Pre-compute degree-to-nodes mapping for all the nodes in the network
degree_map = {}

for node, degree in G.degree():
    if degree not in degree_map:
        degree_map[degree] = []
    degree_map[degree].append(node)

results = []

# Generate baseline
for goid in top_5_go:
    # Find the annotated genes for the top coefficients
    annotated_genes = [
        n for n, attr in G.nodes(data = True)
        if 'go_ids_with_depth' in attr and any(d['go_id'] == goid for d in attr['go_ids_with_depth'])
    ]

    baseline_genes = []

    for real_gene in annotated_genes:
        # Get degree of the real gene
        real_degree = G.degree(real_gene)

        # Select a random node (the node itself can be picked)
        random_node = random.choice(degree_map[real_degree])

        # Compute the neighbor bait sum of the randomly picked gene
        neighbor_sum = sum(G.nodes[neighbor].get('bait_count', 0)
                           for neighbor in G.neighbors(random_node))
        
        baseline_genes.append({
            'gene': str(random_node),
            'neighbor_sum': int(neighbor_sum)
        })

    results.append({
        'annotation': goid,
        'genes_with_sums': baseline_genes
    })

print(f"Saving outputs for GO {aspect} depth {depth} cutoff {cutoff} top coefficients baseline...")

final_df = pd.DataFrame(results)

final_df.to_csv(baseline_df, sep = '\t', index = False)
with open(baseline_pickle, 'wb') as f:
    pickle.dump(final_df, f)

print(f"Baseline for GO {aspect} depth {depth} cutoff {cutoff} top coefficients ready!")