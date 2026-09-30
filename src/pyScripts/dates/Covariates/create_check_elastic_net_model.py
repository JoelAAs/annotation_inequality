import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
import ast

go_file_path = snakemake.input.go_file
hdo_file_path = snakemake.input.hdo_file
bait_usage_path = snakemake.input.bait_usage

confusion_matrix_out = snakemake.output.confusion_matrix
en_coefficients_out = snakemake.output.en_coefficients

# Load Data
print("Loading GO, HDO, and Bait usage files for check EN model", flush=True)
go_df = pd.read_pickle(go_file_path)
hdo_df = pd.read_pickle(hdo_file_path)
bait_df = pd.read_csv(bait_usage_path, sep='\t')

# Format Bait Counts
print("Formatting target variable (log bait counts)", flush=True)
bait_df['entrez_id_bait'] = bait_df['entrez_id_bait'].astype(str)
bait_df = bait_df[bait_df['count'] > 0]
bait_df['log_bait_count'] = np.log(bait_df['count'])

# Construct the Incidence Matrix
print("Constructing the incidence matrix", flush=True)
go_df = go_df.rename(columns={'GO_id': 'term_id'})
hdo_df = hdo_df.rename(columns={'DO_id': 'term_id'})

terms_df = pd.concat([
    go_df[['term_id', 'annotated_genes']], 
    hdo_df[['term_id', 'annotated_genes']]
])

print(f"Exploding {len(terms_df)} terms into individual gene rows", flush=True)
exploded = terms_df.explode('annotated_genes').rename(columns={'annotated_genes': 'gene_id'})
exploded['gene_id'] = exploded['gene_id'].astype(str)
exploded['presence'] = 1

# CALCOLO E PRINT DEI GENI UNICI PRIMA DELLA MATRICE
unique_genes_count = exploded['gene_id'].nunique()
print(f"Trovati {unique_genes_count} geni unici totali all'interno dei due file", flush=True)

print("Pivoting to create final binary matrix", flush=True)
X_matrix = exploded.pivot_table(
    index='gene_id', 
    columns='term_id', 
    values='presence', 
    fill_value=0
)

print(f"Saving incidence matrix to {confusion_matrix_out}", flush=True)
X_matrix.to_parquet(confusion_matrix_out)

# Merge Features and Target
print("Merging feature matrix with target variable", flush=True)
merged = X_matrix.merge(
    bait_df[['entrez_id_bait', 'log_bait_count']], 
    left_index=True, 
    right_on='entrez_id_bait', 
    how='inner'
)

term_columns = X_matrix.columns
X = merged[term_columns]
y = merged['log_bait_count']

# Fit the Elastic Net Model
print(f"Fitting Elastic Net CV on matrix of shape {X.shape}", flush=True)
en_cv = ElasticNetCV(
    l1_ratio=[.1, .5, .7, .9, .95, .99, 1],
    cv=5, 
    random_state=42, 
    n_jobs=-1
)
en_cv.fit(X, y)

# Extract and Save Coefficients
print("Extracting model coefficients", flush=True)
coef_df = pd.DataFrame({
    'term_id': term_columns,
    'coefficient': en_cv.coef_
})

intercept_df = pd.DataFrame({
    'term_id': ['Intercept'], 
    'coefficient': [en_cv.intercept_]
})
final_out = pd.concat([intercept_df, coef_df], ignore_index=True)

print(f"Saving coefficients to {en_coefficients_out}", flush=True)
final_out.to_csv(en_coefficients_out, index=False)

print("Check EN model ready!", flush=True)