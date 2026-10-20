import pandas as pd
import pickle
import networkx as nx

file_network = snakemake.input.final_network
file_coeff = snakemake.input.coefficients
file_output = snakemake.output.nodes_with_top_5_annotations

# Find top 5 DOIDs
print(f"--- [INFO] Loading coefficients... ---")
df_coeff = pd.read_csv(file_coeff, sep='\t')

# Compute absolute value and take top 5
df_coeff['abs_coeff'] = df_coeff['Coefficient'].abs()
top_5_df = df_coeff.sort_values(by='abs_coeff', ascending=False).head(5)
top_5_doids = top_5_df['HDO_doid'].tolist()

print(f"--- [INFO] Top 5 DOIDs by absolute coefficient: ---")
for idx, row in top_5_df.iterrows():
    print(f"    > {row['HDO_doid']}: {row['Coefficient']} (Abs: {row['abs_coeff']})")

# Initialize a dictionary, removing eventual suplicates (shouldn't be present though)
doid_to_genes = {doid: set() for doid in top_5_doids}

# Scan the network to find the annotated nodes
print("\n--- [INFO] Loading network... ---")
with open(file_network, 'rb') as f:
    G = pickle.load(f)

print("--- [INFO] Extracting annotated genes from network... ---")
for node_id, node_data in G.nodes(data=True):
    if 'do_annotations' in node_data:
        for annot in node_data['do_annotations']:
            doid = annot['doid']
            if doid in doid_to_genes:
                doid_to_genes[doid].add(str(node_id))

# Format and Save
print("--- [INFO] Formatting output DataFrame... ---")
output_data = []
for doid in top_5_doids:
    genes_list = sorted(list(doid_to_genes[doid]))
    output_data.append({
        'HDO_doid': doid,
        'annotated_genes': genes_list
    })

df_output = pd.DataFrame(output_data)

print(f"--- [INFO] Saving to {file_output}... ---")
df_output.to_pickle(file_output)

print("--- [INFO] Process completed successfully. ---")