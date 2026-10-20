import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV

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

X_raw = df.loc[common_genes]
y_raw = bait_usage.loc[common_genes, 'count']

# log requires count > 0
tested = y_raw.values > 0
X_raw = X_raw.loc[tested]
y_raw = y_raw.loc[tested]

n_common_genes = len(X_raw)
print(f"{n_common_genes} genes in common found!\n")

print(f"Filtering {aspect} terms (count >= 20)...")

X_raw = X_raw.loc[:, X_raw.sum(axis=0) >= 20]

print(f"{X_raw.shape[1]} terms retained!\n")

X = X_raw.values
y = np.log(y_raw.values)

print(f"Starting GO {aspect} depth {depth} coefficients calculation...")

model = ElasticNetCV(
    l1_ratio=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
    cv=10,
    n_jobs=-1,
    random_state=42,
    max_iter=5000
)

model.fit(X, y)

print(f"Selected alpha (l1_ratio) = {model.l1_ratio_} | lambda (sklearn alpha) = {model.alpha_}")

coef_df = pd.DataFrame({
    'annotation_id': X_raw.columns,
    'coefficient': model.coef_,
    'alpha': model.l1_ratio_,
    'lambda': model.alpha_
})

print(f"{aspect} Coefficients obtained\n")

print("Saving coefficients to parquet")

coef_df.to_parquet(output_coefficients)

print("Coefficients saved!\n")

print(f"Depth {depth} GO {aspect} Elastic Net coefficients ready!")