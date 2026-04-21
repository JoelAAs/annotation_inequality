import pandas as pd

input_df = snakemake.input.merged_df
output_corr = snakemake.output.correlation_values
file_path = snakemake.output.n_of_baits
method = snakemake.wildcards.method

merged_df = pd.read_csv(input_df, sep = '\t', index_col = 0)
merged_df = merged_df.drop(columns = {'entrez_id'})
n_of_baits = len(merged_df)

if method == 'spearman':    
    corr = merged_df.select_dtypes('number').corr(method = 'spearman')

elif method == 'pearson':
    corr = merged_df.select_dtypes('number').corr(method = 'pearson')

new_df = pd.DataFrame({
        'n_of_baits': [n_of_baits]
    })

new_df.to_csv(file_path, sep = '\t', index = False)

corr_df = pd.DataFrame(corr, index = merged_df.columns, columns = merged_df.columns)
corr_df = corr_df.round(3)
corr_df.to_csv(output_corr, sep = '\t')

