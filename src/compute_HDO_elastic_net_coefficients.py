import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from scipy.sparse import csr_matrix

input_df_annotations = snakemake.input.annotation_df
input_df_bait_usage = snakemake.input.bait_usage
output_df = snakemake.output.elastic_net_coefficients

df_annot = pd.read_csv(input_df_annotations, sep = '\t')
bait_usage = pd.read_csv(input_df_bait_usage, sep = '\t')

df_annot = df_annot.set_index('entrez_id')
bait_usage = bait_usage.set_index('entrez_id_bait')

common_genes = df_annot.index.intersection(bait_usage.index)
X = df_annot.loc[common_genes].values
y = np.log1p(bait_usage.loc[common_genes, 'count'].values)

model = ElasticNetCV(
    l1_ratio=[.1, .5, .7, .9, .95, .99, 1],
    cv=5, 
    n_jobs=-1, 
    random_state=42
)

model.fit(X,y)

coef_df = pd.DataFrame({
    'HDO_Term': df_annot.columns,
    'Coefficient': model.coef_
})

coef_df.to_csv(output_df, sep = '\t')