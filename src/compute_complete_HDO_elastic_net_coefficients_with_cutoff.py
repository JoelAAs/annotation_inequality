import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler

cutoff = snakemake.wildcards.cutoff
input_matrix = snakemake.input.complete_matrix
input_df_bait_usage = snakemake.input.bait_usage
output_coefficients = snakemake.output.complete_elastic_net_coefficients

print(f"Processing complete HDO Elastic Net coefficients with cutoff {cutoff}...\n")

print("Loading input dfs...")

df = pd.read_csv(input_matrix, sep = '\t')
df = df.set_index('entrez_id')

bait_usage = pd.read_csv(input_df_bait_usage, sep = '\t')
bait_usage = bait_usage.set_index('entrez_id_bait')

print("Input dfs loaded!\n")

print("Finding genes in common...")

print(f"Unique IDs in df = {len(df.index)}")
print(f"Unique IDs in bait_usage = {len(bait_usage.index)}")
common_genes = df.index.intersection(bait_usage.index)
n_common_genes = len(common_genes)

print(f"{n_common_genes} genes in common found!\n")

print("Scaling data...")

X = df.loc[common_genes].values
y = np.log1p(bait_usage.loc[common_genes, 'count'].values)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("Data scaled!\n")

print("Starting coefficients calculation...")

model = ElasticNetCV(
    l1_ratio=[.1, .5, .7, .9, .95, .99, 1],
    cv=5, 
    n_jobs=-1, 
    random_state=42,
    max_iter=5000,
    tol=1e-3
)

model.fit(X_scaled,y)

coef_df = pd.DataFrame({
    'HDO_doid': df.columns,
    'Coefficient': model.coef_
})

print("Coefficients obtained\n")

print("Saving coefficients to csv...")

coef_df.to_csv(output_coefficients, sep = '\t', index = False)

print("Coefficients saved!\n")

print(f"Complete HDO Elastic Net coefficients with cutoff {cutoff} ready!")