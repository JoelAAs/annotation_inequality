import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

matrix_path = snakemake.input.confusion_matrix
coefs_path = snakemake.input.en_coefficients

matrix_plot_out = snakemake.output.matrix_plot
coef_plot_out = snakemake.output.coef_plot

# Load Data
print("--- Loading data for visualization ---", flush=True)
matrix = pd.read_parquet(matrix_path)
coefs = pd.read_csv(coefs_path)

# Plot Incidence Matrix (Sparsity Pattern)
print("--- Plotting incidence matrix sparsity ---", flush=True)
plt.figure(figsize=(10, 8))

# plt.spy efficiently plots non-zero elements in a matrix
plt.spy(matrix, aspect='auto', markersize=2, color='navy')

plt.title("Gene-Term Incidence Matrix", fontsize=16, pad=15)
plt.xlabel("Terms (GO & HDO)", fontsize=14)
plt.ylabel(f"Genes (n={matrix.shape[0]})", fontsize=14)
plt.xticks(ticks=range(len(matrix.columns)), labels=matrix.columns, rotation=45, ha='right')

plt.tight_layout()
plt.savefig(matrix_plot_out, dpi=300)
plt.close()

# Plot Elastic Net Coefficients
print("--- Plotting Elastic Net coefficients ---", flush=True)

# Separate the intercept so it doesn't distort the bar chart scale
intercept_row = coefs[coefs['term_id'] == 'Intercept']
intercept_val = intercept_row['coefficient'].values[0] if not intercept_row.empty else 0

# Sort the remaining terms by their coefficient values
terms_only = coefs[coefs['term_id'] != 'Intercept'].sort_values(by='coefficient', ascending=False)

plt.figure(figsize=(10, 6))
# Using a diverging palette to clearly show positive vs negative weights
sns.barplot(data=terms_only, x='coefficient', y='term_id', palette='vlag')

plt.title(f"Elastic Net Coefficients\n(Intercept: {intercept_val:.3f})", fontsize=16, pad=15)
plt.xlabel("Coefficient Value (Impact on Log Bait Count)", fontsize=14)
plt.ylabel("Terms", fontsize=14)

# Add a vertical line at 0 for reference
plt.axvline(x=0, color='black', linestyle='-', linewidth=1)

plt.tight_layout()
plt.savefig(coef_plot_out, dpi=300)
plt.close()

print(f"--- Visualizations saved to {matrix_plot_out} and {coef_plot_out} ---", flush=True)