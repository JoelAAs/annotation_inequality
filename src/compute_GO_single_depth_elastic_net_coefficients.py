import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
input_matrix = snakemake.input.feature_matrix
input_df_bait_usage = snakemake.input.bait_usage
output_coefficients = snakemake.output.single_depth_elastic_net_coefficients

print(f"Processing depth {depth} GO {aspect} Elastic Net coefficients...\n")

print(f"Loading {aspect} input dfs...")

df = pd.read_csv(input_matrix, sep = '\t')
df = df.set_index('entrez_id')

bait_usage = pd.read_csv(input_df_bait_usage, sep = '\t')
bait_usage = bait_usage.set_index('entrez_id_bait')

print("Input dfs loaded!\n")

print("Finding genes in common...")

common_genes = df.index.intersection(bait_usage.index)
n_common_genes = len(common_genes)

X_raw = df.loc[common_genes]
y_raw = bait_usage.loc[common_genes, 'count']

print(f"{n_common_genes} genes in common found!\n")

print(f"Scaling {aspect} data...")

X = X_raw.values
y = np.log1p(y_raw.values)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"{aspect} Data scaled!\n")

print(f"Starting GO {aspect} depth {depth} coefficients calculation...")

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
    'GO_id': X_raw.columns,
    'Coefficient': model.coef_
})

print(f"{aspect} Coefficients obtained\n")

print("Saving coefficients to csv...")

coef_df.to_csv(output_coefficients, sep = '\t', index = False)

print("Coefficients saved!\n")

print(f"Depth {depth} GO {aspect} Elastic Net coefficients ready!")