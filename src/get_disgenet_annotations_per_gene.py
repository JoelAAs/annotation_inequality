import pandas as pd

input_df_annotations = snakemake.input.annotation_df
output_annotation_df = snakemake.output.annotation_df
output_annotations_list = snakemake.output.annotation_list

df_annotations = pd.read_csv(input_df_annotations, sep = '\t')

df_annotations = df_annotations.dropna() 
df_annotations = df_annotations.drop_duplicates(subset = ['geneid', 'diseaseUMLSCUI'])
df_annotations = df_annotations[['geneid', 'disease_name']]
df_annotations = df_annotations.rename(columns = {
    'geneid': 'entrez_id'
})

print("DISGENET annotations obtained!")

annot_list = df_annotations[['disease_name']]
annot_list = annot_list.drop_duplicates()

df_annotations.to_csv(output_annotation_df, sep = '\t')
annot_list.to_csv(output_annotations_list, sep = '\t')
