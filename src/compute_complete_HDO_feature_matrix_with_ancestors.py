import pandas as pd

input_df = snakemake.input.complete_annotations
output_matrix = snakemake.output.complete_matrix

print("Processing HDO complete matrix with ancestors...\n")

print("Loading input df and filling missing values...")

df = pd.read_csv(input_df, sep = '\t')
df.fillna({
    'doid': 'No_doid'
})

print("Done!\n")

print("Creating feature matrix...\n")

df['value'] = 1

df_onehot = df.pivot_table(
    index = 'entrez_id',
    columns = 'doid',
    values = 'value',
    fill_value = 0
).reset_index()

df_onehot = df_onehot.drop(columns = {'No_doid'})
df_onehot.columns.name = None

unique_values = df_onehot.shape[0]
print(f'Unique entrez_id values = {unique_values}')
print(f'Unique doids = {len(df_onehot.columns) - 1}\n')

print("Feature matrix ready!\n")

df_onehot.to_csv(output_matrix, sep = '\t', index = False)

print("HDO complete matrix with ancestors ready!")