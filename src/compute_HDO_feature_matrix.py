import pandas as pd

input_df = snakemake.input.annotation_df
input_list = snakemake.input.annotation_list
output_matrix = snakemake.output.feature_matrix

df = pd.read_csv(input_df, sep = '\t')
annot_list = pd.read_csv(input_list, sep = '\t')

df['annotations_list'] = df['annotations'].str.split(';')
df = df.explode('annotations_list')
df['value'] = 1

df_onehot = df.pivot_table(
    index = 'entrez_id',
    columns = 'annotations_list',
    values = 'value',
    fill_value = 0
).reset_index()

# TODO review this!!!

df_onehot.columns.name = None

df_onehot.to_csv(output_matrix, sep = '\t', index = False)