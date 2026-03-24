import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
import os

input_df_bait_usage = snakemake.input.bait_usage

bait_usage = pd.read_csv(input_df_bait_usage, sep = '\t')
bait_usage = bait_usage.set_index('entrez_id_bait')

for inputdir, outputdir in zip(snakemake.input[0], snakemake.output[0]):
    os.makedirs(outputdir, exist_ok = True)
    matrices = os.listdir(inputdir)
    
    for feature_matrix in matrices:
        df_annot = pd.read_csv(feature_matrix, sep = '\t')
        df_annot = df_annot.set_index('entrez_id_bait')
        aspect = os.path.basename(os.path.dirname(feature_matrix))
        fname = os.path.basename(feature_matrix)
        depth = fname.split('_')[-1].replace('.csv', '')
        outputfile = os.path.join(outputdir, f'elastic_net_coefficients_depth_{depth}.csv')

        print(f'Starting depth {depth} GO {aspect} elastic net coefficients calculation...')

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
            f'GO_{aspect}_Term': df_annot.columns,
            'Coefficient': model.coef_
        })

        coef_df.to_csv(outputfile, sep = '\t', index = False)

        print(f'Depth {depth} GO {aspect} elastic net coefficients obtained!')
