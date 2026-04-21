import pandas as pd
import pickle
import networkx as nx

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
inputnetwork = snakemake.input.bp_network_with_attributes
coefficients_df = snakemake.input.coefficients
outputdf = snakemake.output.neighbors_bait_count_sums_df
outputpickle = snakemake.output.neighbors_bait_count_sums_pickle

# Load the inputs from Snakemake
with open(inputnetwork, 'rb') as f:
    G = pickle.load(f)

# Assume the network nodes have an attribute 'bait_count' 
# and an attribute 'annotations' which is a list of HDO IDs
df_coeffs = pd.read_csv(coefficients_df, sep = '\t')

# Identify Top 5 Coefficients (by absolute value to include strong negative ones)
df_coeffs['abs_coeff'] = df_coeffs['Coefficient'].abs()
top_5_go = df_coeffs.nlargest(5, 'abs_coeff')['GO_id'].tolist()

results = []

# Process each top annotation
for go_id in top_5_go:
    annotated_genes = [
        n for n, attr in G.nodes(data=True) 
        if 'go_ids_with_depth' in attr and any(g['go_id'] == go_id for g in attr['go_ids_with_depth'])
    ]
    
    genes_with_sums = []
    
    for gene in annotated_genes:
        # Sum the 'bait_count' of all immediate neighbors in the graph
        neighbor_sum = sum(G.nodes[neighbor].get('bait_count', 0) 
                          for neighbor in G.neighbors(gene))
        
        genes_with_sums.append({
            'gene': str(gene),
            'neighbor_sum': int(neighbor_sum)
        })
    
    results.append({
        'annotation': go_id,
        'genes_with_sums': genes_with_sums
    })

# Save Outputs
final_df = pd.DataFrame(results)

final_df.to_csv(outputdf, sep = '\t', index=False)

with open(outputpickle, 'wb') as f:
    pickle.dump(final_df, f)

print(f"Computed neighbor sums for top 5 HDO annotations at depth {depth}")