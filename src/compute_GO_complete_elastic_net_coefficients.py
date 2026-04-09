import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler

aspect = snakemake.wildcards.aspect
feature_matrix = snakemake.input.feature_matrix
input_df_bait_usage = snakemake.input.bait_usage
coefficients = snakemake.output.complete_elastic_net_coefficients

print(f"Processing complete GO {aspect} Elastic Net coefficients...\n")

print(f"Loading {aspect} input dfs...")

df = pd.read_csv(feature_matrix, sep = '\t')
df = df.set_index('entrez_id')

bait_usage = pd.read_csv(input_df_bait_usage, sep = '\t')
bait_usage = bait_usage.set_index('entrez_id_bait')

print(f"{aspect} Input dfs loaded!\n")

print("Finding genes in common...")

common_genes = df.index.intersection(bait_usage.index)
n_common_genes = len(common_genes)

print(f"{n_common_genes} genes in common found!\n")

print(f"Scaling {aspect} data...")

X_df = df.loc[common_genes]
print(f"{aspect} features before removing constant columns = {X_df.shape[1]}")
X_df = X_df.loc[:, X_df.std() > 0]
print(f"{aspect} features remaining after removing constant columns = {X_df.shape[1]}")
X = X_df.values
y = np.log1p(bait_usage.loc[common_genes, 'count'].values)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"{aspect} Data scaled!\n")

print(f"Starting {aspect} coefficients calculation...")

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
    'GO_id': X_df.columns,
    'Coefficient': model.coef_
})

print(f"{aspect} Coefficients obtained\n")

print(f"Saving {aspect} coefficients to csv...")

coef_df.to_csv(coefficients, sep = '\t', index = False)

print(f"{aspect} Coefficients saved!\n")

print(f"Complete GO {aspect} Elastic Net coefficients ready!")