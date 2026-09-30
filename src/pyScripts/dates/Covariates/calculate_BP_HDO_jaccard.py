import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pronto

input_go = snakemake.input.go_file
input_hdo = snakemake.input.hdo_file
obo_go = snakemake.input.go_obo
obo_hdo = snakemake.input.hdo_obo

output_tsv = snakemake.output.jaccard_tsv
output_plot = snakemake.output.jaccard_plot

# Load ontologies using pronto
go_ontology = pronto.Ontology(obo_go)
hdo_ontology = pronto.Ontology(obo_hdo)

def get_term_name(term_id, ontology):
    try:
        name = ontology[term_id].name
        return f"{term_id} ({name})"
    except KeyError:
        return term_id

def jaccard_index(set1, set2):
    if not set1 and not set2:
        return 0.0
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union

# Load data
df_go = pd.read_pickle(input_go).head(5)
df_hdo = pd.read_pickle(input_hdo).head(5)

# Extract gene sets, translate IDs, and append gene counts
go_dict = {}
for _, row in df_go.iterrows():
    genes = set(row['annotated_genes'])
    base_name = get_term_name(row['GO_id'], go_ontology)
    label = f"{base_name}\n(n={len(genes)})"
    go_dict[label] = genes

hdo_dict = {}
for _, row in df_hdo.iterrows():
    genes = set(row['annotated_genes'])
    base_name = get_term_name(row['DO_id'], hdo_ontology)
    label = f"{base_name}\n(n={len(genes)})"
    hdo_dict[label] = genes

# Calculate global Jaccard index
all_go_genes = set.union(*go_dict.values()) if go_dict else set()
all_hdo_genes = set.union(*hdo_dict.values()) if hdo_dict else set()
global_jaccard = jaccard_index(all_go_genes, all_hdo_genes)

# Calculate Pairwise Jaccard Matrix
matrix_data = []
for go_name, go_genes in go_dict.items():
    row_data = {'GO_Term': go_name}
    for hdo_name, hdo_genes in hdo_dict.items():
        row_data[hdo_name] = jaccard_index(go_genes, hdo_genes)
    matrix_data.append(row_data)

df_matrix = pd.DataFrame(matrix_data)

# Save tabular results
with open(output_tsv, 'w') as f:
    f.write(f"# Global Jaccard Index (Union of Top 5 GO vs Union of Top 5 HDO): {global_jaccard:.4f}\n")
    f.write("# Pairwise Jaccard Matrix:\n")
    df_matrix.to_csv(f, sep='\t', index=False)

# Generate and save the heatmap plot
df_plot = df_matrix.set_index('GO_Term')

plt.figure(figsize=(11, 9))
sns.heatmap(df_plot, annot=True, cmap="YlGnBu", vmin=0, vmax=1, fmt=".3f", 
            cbar_kws={'label': 'Jaccard Index'})

plt.title(f"Jaccard Overlap: GO vs HDO\nGlobal Overlap = {global_jaccard:.3f}", pad=20)
plt.xlabel("HDO Terms", labelpad=15)
plt.ylabel("GO Terms", labelpad=15)
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()