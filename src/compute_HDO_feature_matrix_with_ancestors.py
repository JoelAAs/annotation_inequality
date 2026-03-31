import pandas as pd

input_df = snakemake.input.complete_annotations
outputs = snakemake.output.feature_matrix
textfile = snakemake.output.annotations_per_depth
counts = snakemake.output.counts_per_annot

df = pd.read_csv(input_df, sep = '\t')
df.fillna({
    'doid': 'No_doid',
    'depth': -1
}, inplace=True)

with open(textfile, 'w') as f:
    f.write(f'depth:n_of_coefficients\n')

for outputfile, count_file in zip(outputs, counts):
    depth = int(outputfile.split("_")[-1].replace(".csv", ""))
    print(f'Creating HDO feature matrix with depth {depth}...')

    df_mod = df.copy()
    print(f"Processing depth {depth}")

    df_mod.loc[df_mod['depth'] != depth, 'doid'] = 'No_doid'
    df_mod['value'] = 1

    df_onehot = df_mod.pivot_table(
        index = 'entrez_id',
        columns = 'doid',
        values = 'value',
        fill_value = 0
    ).reset_index()

    df_onehot = df_onehot.drop(columns = {'No_doid'})
    df_onehot.columns.name = None

    unique_values = df_onehot.shape[0]
    print(f'Unique entrez_id values = {unique_values}')
    print(f'Unique doids = {len(df_onehot.columns) - 1}')
    
    with open(textfile, 'a') as f:
        f.write(f'{depth}:{len(df_onehot.columns) - 1}\n')
    
    print(f'HDO coefficient count for depth {depth} inserted in the text file!')

    df_onehot.to_csv(outputfile, sep = '\t', index = False)
    print(f'HDO feature matrix with depth {depth} ready!')

    print(f'Counting HDO depth {depth} annotations counts...')

    df_count = df.copy()

    df_count = df_count[df_count['depth'] == depth]

    counts = df_count['doid'].value_counts().reset_index()
    counts.columns = ['doid', 'count']

    counts.to_csv(count_file, sep = '\t', index = False)

    print(f'HDO depth {depth} annotations counts obtained!\n')