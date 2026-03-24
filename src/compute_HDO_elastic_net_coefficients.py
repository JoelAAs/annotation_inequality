import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV

input_feat_matrix = snakemake.input.feature_matrix
input_df_bait_usage = snakemake.input.bait_usage
output_coefficients = snakemake.output.elastic_net_coefficients

bait_usage = pd.read_csv(input_df_bait_usage, sep = '\t')
bait_usage = bait_usage.set_index('entrez_id_bait')

for feat_matrix, outputfile in zip(input_feat_matrix, output_coefficients):
    df_annot = pd.read_csv(feat_matrix, sep = '\t')
    df_annot = df_annot.set_index('entrez_id')
    depth = outputfile.split("_")[-1].replace(".csv", "")

    print(f'Starting depth {depth} HDO elastic net coefficients calculation...')

    common_genes = df_annot.index.intersection(bait_usage.index)
    n_common_genes = len(common_genes)

    print(f'Number of common genes = {n_common_genes}')

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

    print(f'Depth {depth} HDO elastic net coefficients obtained!')

    coef_df.to_csv(outputfile, sep = '\t', index = False)