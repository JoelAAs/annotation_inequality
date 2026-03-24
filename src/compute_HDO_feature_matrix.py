import pandas as pd

input_df = snakemake.input.annotation_df
outputs = snakemake.output.feature_matrix
textfile = snakemake.output.annotations_per_depth

df = pd.read_csv(input_df, sep = '\t')
df = df.drop(columns = {'doid'})
df.fillna({
    'annotation': 'No_annot',
    'depth': -1
}, inplace=True)

with open(textfile, 'w') as f:
    f.write(f'depth:n_of_coefficients\n')

for outputfile in outputs:
    df_mod = df.copy()
    depth = int(outputfile.split("_")[-1].replace(".csv", ""))
    print(f"Processing depth {depth}")

    df_mod.loc[df_mod['depth'] != depth, 'annotation'] = 'No_annot'
    df_mod['value'] = 1

    print(f'Creating HDO feature matrix with depth {depth}...')

    df_onehot = df_mod.pivot_table(
        index = 'entrez_id',
        columns = 'annotation',
        values = 'value',
        fill_value = 0
    ).reset_index()

    df_onehot = df_onehot.drop(columns = {'No_annot'})
    df_onehot.columns.name = None

    unique_values = df_onehot['entrez_id'].nunique
    print(f'Unique entrez_id values = {unique_values}\n')
    print(f'Unique annotations = {len(df_onehot.columns)}')
    
    with open(textfile, 'a') as f:
        f.write(f'{depth}:{len(df_onehot.columns)}\n')
    
    print(f'HDO coefficient count for depth {depth} inserted in the text file!')

    df_onehot.to_csv(outputfile, sep = '\t', index = False)
    print(f'HDO feature matrix with depth {depth} ready!')
