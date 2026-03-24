import pandas as pd
import os

for input_df, output_dir in zip(snakemake.input, snakemake.output):
    aspect = os.path.basename(input_df).split('_')[0]
    df = pd.read_csv(input_df, sep = '\t')
    df = df.drop(columns = 'go_id')
    df = df.rename(columns = {'go_term_name': 'annotation'})
    df.fillna({
        'annotation': 'No_annot',
        'depth': -1    
    }, inplace = True)
    
    max_depth = int(df['depth'].max())
    os.makedirs(output_dir, exist_ok = True)

    for depth in range(max_depth + 1):
        print(f'Processing GO {aspect} feature matrix with depth {depth}...')

        outputfile = os.path.join(output_dir, f'feature_matrix_with_depth_{depth}.csv')

        df_mod = df.copy()
        df_mod.loc[df_mod['depth'] != depth, 'annotation'] = 'No_annot'
        df_mod['value'] = 1

        df_onehot = df_mod.pivot_table(
            index = 'entrez_id',
            columns = 'annotation',
            values = 'value',
            fill_value = 0
        ).reset_index()

        df_onehot = df_onehot.drop(columns = {'No_annot'})
        df_onehot.columns.name = None

        print(f'Unique annotations = {len(df_onehot.columns)}')


        df_onehot.to_csv(outputfile, sep = '\t', index = False)
        print(f'GO {aspect} feature matrix with depth {depth} ready!')