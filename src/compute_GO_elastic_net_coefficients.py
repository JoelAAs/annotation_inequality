import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
import os

aspect = snakemake.wildcards.aspect
depth = int(snakemake.wildcards.depth)
feature_matrix = snakemake.input.feature_matrix
bait_file = snakemake.input.bait_usage
outputfile = snakemake.output.outputfile

bait_usage = pd.read_csv(bait_file, sep = '\t')
bait_usage = bait_usage.set_index('entrez_id_bait')
df = pd.read_csv(feature_matrix, sep = '\t')
df = df.set_index('entrez_id')

print(f'Processing GO {aspect} depth {depth} Elastic Net Coefficients...')

common_genes = df.index.intersection(bait_usage.index)
n_common_genes = len(common_genes)

X = df.loc[common_genes].values
y = np.log1p(bait_usage.loc[common_genes, 'count'].values)

model = ElasticNetCV(
    l1_ratio=[.1, .5, .7, .9, .95, .99, 1],
    cv=5, 
    n_jobs=-1, 
    random_state=42
)

model.fit(X,y)

coef_df = pd.DataFrame({
    f'GO_term': df.columns,
    'Coefficient': model.coef_
})

coef_df.to_csv(outputfile, sep = '\t', index = False)

print(f'GO {aspect} depth {depth} coefficients obtained!')
