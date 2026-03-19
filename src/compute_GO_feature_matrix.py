import pandas as pd

ASPECTS = ['BP', 'MF', 'CC']

for aspect in ASPECTS:
    input_df = [p for p in snakemake.input.annotation_df if f"{aspect}_" in p][0]    
    output_df = [p for p in snakemake.output.feature_matrix if f"{aspect}_" in p][0]

    df = pd.read_csv(input_df, sep  ='\t')

    df['value'] = 1

    df_onehot = df.pivot_table(
        index = 'entrez_id',
        columns = 'go_term_name',
        values = 'value',
        fill_value = 0
    ).reset_index()
    
    df_onehot.columns.name = None

    df_onehot.to_csv(output_df, sep = '\t', index = False)

    print(f'{aspect} feature matrix ready!')