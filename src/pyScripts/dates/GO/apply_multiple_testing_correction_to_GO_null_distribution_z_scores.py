import pandas as pd
from statsmodels.stats.multitest import multipletests

input_file = snakemake.input.null_master_file
output_file = snakemake.output.corrected_file

print(f"\n--- [Core] Computing multiple-testing corrected p-values... ---\n")

print(f"--- [Core] Loading the master null adjacency distribution file... ---\n")

df = pd.read_csv(input_file, sep = '\t')

print(f"--- [Core] Applying Benjamini-Hochberg FDR correction to the empirical p-values obtained... ---\n")

reject, pvals_corrected, _, _ = multipletests(
    df['Empirical_P_Value'],
    alpha = 0.05,
    method = 'fdr_bh'
)

df['FDR_Adjusted_P_Value'] = pvals_corrected
df['Is_Significant_FDR'] = reject

print(f"--- [Core] Saving to: {output_file} ---\n")

df.to_csv(output_file, sep = '\t', index = False)

print(f"--- [Core] Succesfully saved to: {output_file}! ---\n")