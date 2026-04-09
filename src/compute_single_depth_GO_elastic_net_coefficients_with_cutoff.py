import pandas as pd
import numpy as np
import sys
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler

aspect = snakemake.wildcards.aspect
cutoff = snakemake.wildcards.cutoff
depth = snakemake.wildcards.depth
input_matrix = snakemake.input.single_depth_feature_matrix_with_cutoff
input_df_bait_usage = snakemake.input.bait_usage
output_coefficients = snakemake.output.single_depth_elastic_net_coefficients

print(f"Processing depth {depth} GO {aspect} Elastic Net coefficients with cutoff {cutoff}...\n")

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

X_raw = df.loc[common_genes]
y_raw = bait_usage.loc[common_genes, 'count']

print(f"{n_common_genes} genes in common found!\n")

print(f"Scaling {aspect} depth {depth} cutoff {cutoff} data...")

if X_raw.shape[1] == 0:
    print(f"No features available for {aspect} Depth {depth}, Cutoff {cutoff}. Skipping...")
    # Create an empty DataFrame with the expected columns so Snakemake doesn't fail
    pd.DataFrame(columns=['GO_id', 'Coefficient']).to_csv(output_coefficients, sep = '\t', index=False)
    sys.exit(0) 

X = X_raw.values
y = np.log1p(y_raw.values)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"{aspect} depth {depth} cutoff {cutoff} Data scaled!\n")

print(f"Starting {aspect} depth {depth} cutoff {cutoff} coefficients calculation...")

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
    'GO_id': df.columns,
    'Coefficient': model.coef_
})

print(f"{aspect} depth {depth} cutoff {cutoff} Coefficients obtained\n")

print("Saving coefficients to csv...")

coef_df.to_csv(output_coefficients, sep = '\t', index = False)

print("Coefficients saved!\n")

print(f"Depth {depth} GO {aspect} Elastic Net coefficients with cutoff {cutoff} ready!")