import pandas as pd
import matplotlib.pyplot as plt

feature_matrix_path = snakemake.input.feature_matrix
parameters_output_path = snakemake.output.parameters
counts_plot_path = snakemake.output.counts_plot
terms_drop_plot_path = snakemake.output.terms_drop_plot 
excluded_terms_df_path = snakemake.output.dropped_genes_df

print(f"--- [START] MIXED FEATURE MATRIX PARAMETERS COMPUTATION STARTED ---", flush=True)

# Load feature matrix
print(f"--- [RUNNING] LOADING FEATURE MATRIX ---", flush=True)
feature_matrix = pd.read_parquet(feature_matrix_path)
total_genes = feature_matrix.shape[0]
total_terms = feature_matrix.shape[1]
print(f"--- [INFO] LOADED MATRIX: {total_genes} GENES AND {total_terms} ANNOTATIONS ---", flush=True)

# Calculate counts and fractions
print(f"--- [RUNNING] CALCULATING COUNTS AND FRACTIONS ---", flush=True)
annotation_counts = feature_matrix.sum(axis=0)
annotation_fractions = annotation_counts / total_genes

parameters_df = pd.DataFrame({
    'count': annotation_counts,
    'fraction': annotation_fractions
}).reset_index()

parameters_df.rename(columns={'index': 'annotation_id', 'Features': 'annotation_id', 'annotation_id': 'annotation_id'}, inplace=True)

parameters_df.sort_values(by='count', ascending=False, inplace=True)
parameters_df.reset_index(drop=True, inplace=True)
print(f"--- [INFO] PARAMETERS CALCULATED AND SORTED ---", flush=True)

# Generate counts histogram
print(f"--- [RUNNING] GENERATING COUNTS HISTOGRAM ---", flush=True)
plt.figure(figsize=(10, 6))
plt.hist(parameters_df['count'], bins=100, color='blue', edgecolor='black', log=True)
plt.title('Distribution of Genes per Annotation (Histogram)')
plt.xlabel('Number of Genes (Count)')
plt.ylabel('Number of Annotations (Log Scale)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(counts_plot_path, dpi=300)
plt.close()
print(f"--- [INFO] COUNTS HISTOGRAM SAVED ---", flush=True)

# Calculate excluded terms across thresholds
print(f"--- [RUNNING] CALCULATING EXCLUDED TERMS FOR THRESHOLDS 1-10 ---", flush=True)
thresholds = list(range(1, 51))
excluded_terms_counts = []
excluded_terms_percentages = []

for t in thresholds:
    # A term is excluded if it appears in fewer genes than the threshold 't'
    excluded_count = (parameters_df['count'] < t).sum()
    excluded_percentage = (excluded_count / total_terms) * 100
    
    excluded_terms_counts.append(excluded_count)
    excluded_terms_percentages.append(excluded_percentage)
    print(f"--- [INFO] THRESHOLD {t}: {excluded_count} TERMS EXCLUDED ({excluded_percentage:.2f}%) ---", flush=True)

# Generate excluded terms plot
print(f"--- [RUNNING] GENERATING EXCLUDED TERMS PLOT ---", flush=True)
plt.figure(figsize=(10, 6))
plt.plot(thresholds, excluded_terms_percentages, marker='o', color='red', linewidth=2, markersize=8)

# Add text annotations above each dot with a white background box
for x, y in zip(thresholds, excluded_terms_percentages):
    plt.annotate(
        f"{y:.3f}%", 
        (x, y), 
        textcoords="offset points", 
        xytext=(0, 12),
        ha='center', 
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8)
    )

plt.title('Impact of Minimum Gene Threshold on Annotation Retention')
plt.xlabel('Genes Count Threshold')
plt.ylabel('Percentage of Annotations Excluded (%)')
plt.xticks(thresholds)

plt.ylim(bottom=-3, top=max(excluded_terms_percentages) + 10)

plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(terms_drop_plot_path, dpi=300)
plt.close()
print(f"--- [INFO] EXCLUDED TERMS PLOT SAVED ---", flush=True)

# Save excluded terms percentages to Parquet
print(f"--- [RUNNING] SAVING EXCLUDED TERMS PERCENTAGES TO PARQUET ---", flush=True)
excluded_terms_df = pd.DataFrame({
    'threshold': thresholds,
    'excluded_terms_count': excluded_terms_counts,
    'excluded_terms_percentage': excluded_terms_percentages
})
excluded_terms_df.to_parquet(excluded_terms_df_path)
print(f"--- [INFO] EXCLUDED TERMS DATAFRAME SAVED ---", flush=True)

# Save parameters to Parquet
print(f"--- [RUNNING] SAVING PARAMETERS TO PARQUET ---", flush=True)
parameters_df.to_parquet(parameters_output_path)

print(f"--- [COMPLETE] MIXED FEATURE MATRIX PARAMETERS COMPUTATION FINISHED SUCCESSFULLY ---", flush=True)