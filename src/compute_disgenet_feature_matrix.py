import pandas as pd

input_df = snakemake.input.annotation_df
outputfile = snakemake.output.feature_matrix

df = pd.read_csv(input_df, sep = '\t')

df['value'] = 1
print(f'Creating DISGENET feature matrix...')

df_onehot = df.pivot_table(
    index = 'entrez_id',
    columns = 'disease_name',
    values = 'value',
    fill_value = 0
).reset_index()

df_onehot.columns.name = None

df_onehot.to_csv(outputfile, sep = '\t', index = False)

print(f'DISGENET feature matrix ready!')