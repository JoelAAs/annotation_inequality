import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pronto
from matplotlib.patches import Patch

coef_path = snakemake.input.en_coefficients
metrics_path = snakemake.input.en_metrics
go_obo_path = snakemake.input.go_obo
hdo_obo_path = snakemake.input.hdo_obo

top20_plot_dir = snakemake.output.top20_plot
metrics_plot_out = snakemake.output.metrics_plot

print("--- [START] MIXED ELASTIC NET MODEL VISUALIZATIONS STARTED ---", flush=True)

# PLOTTING ALPHA METRICS
print("--- [RUNNING] GENERATING METRICS PLOT ---", flush=True)
metrics_df = pd.read_csv(metrics_path, sep="\t")

# Extract all tested alpha values for the x-axis ticks
alpha_ticks = metrics_df['alpha'].tolist()
# Extract baseline deviance
baseline_val = metrics_df['baseline_test_dev_explained'].iloc[0]

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Plot 1: Deviances
axes[0].plot(metrics_df['alpha'], metrics_df['cv_dev_explained_min'], marker='o', color='#99ccff', linewidth=2, markersize=8, label='CV Dev Explained (min)', zorder=3)
axes[0].plot(metrics_df['alpha'], metrics_df['cv_dev_explained_1se'], marker='s', color='#ffcc99', linewidth=2, markersize=8, label='CV Dev Explained (1se)', zorder=3)
axes[0].plot(metrics_df['alpha'], metrics_df['test_dev_explained'], marker='^', color='#ff9999', linewidth=2, markersize=8, label='Test Dev Explained (20% Hold-out set)', zorder=3)

# Add horizontal line for the baseline test deviance
axes[0].axhline(y=baseline_val, color='#666666', linestyle='--', linewidth=2, alpha=0.8, zorder=2, label=f'Baseline Test Dev Expl. ({baseline_val:.4f})')

axes[0].set_title('Deviance Explained', fontsize=14)
axes[0].set_xlabel('Alpha (0 = Ridge, 1 = Lasso)', fontsize=12)
axes[0].set_ylabel('Deviance Explained', fontsize=12)
axes[0].set_xticks(alpha_ticks)
axes[0].grid(True, linestyle='--', alpha=0.7, zorder=1)
axes[0].legend(loc='best')

# Plot 2: Non-zero features
axes[1].plot(metrics_df['alpha'], metrics_df['nonzero_min'], marker='o', color='#c2c2f0', linewidth=2, markersize=8, label='Non-zero features (min)', zorder=3)
axes[1].plot(metrics_df['alpha'], metrics_df['nonzero_1se'], marker='s', color='#ffb3e6', linewidth=2, markersize=8, label='Non-zero features (1se)', zorder=3)

axes[1].set_title('Feature Selection', fontsize=14)
axes[1].set_xlabel('Alpha (0 = Ridge, 1 = Lasso)', fontsize=12)
axes[1].set_ylabel('Non-Zero Annotations', fontsize=12)
axes[1].set_xticks(alpha_ticks)
axes[1].grid(True, linestyle='--', alpha=0.7, zorder=1)
axes[1].legend(loc='best')

plt.tight_layout()
plt.savefig(metrics_plot_out, dpi=300)
plt.close()
print(f"--- [INFO] METRICS PLOT SAVED TO {metrics_plot_out} ---", flush=True)

# PROCESSING TOP 20 COEFFICIENTS
print("--- [RUNNING] LOADING ONTOLOGIES WITH PRONTO ---", flush=True)
go_ontology = pronto.Ontology(go_obo_path)
hdo_ontology = pronto.Ontology(hdo_obo_path)

def get_term_name(term_id):
    if term_id in go_ontology:
        return go_ontology[term_id].name
    elif term_id in hdo_ontology:
        return hdo_ontology[term_id].name
    return term_id

print("--- [RUNNING] PROCESSING COEFFICIENTS ALONG ALPHAS AND THRESHOLDS ---", flush=True)
coef_df_full = pd.read_parquet(coef_path)

# Remove intercept and unpenalized covariate globally
coef_df_full = coef_df_full[~coef_df_full['annotation_id'].isin(['(Intercept)', 'n_annot'])]

# Create the output directory if it doesn't exist
os.makedirs(top20_plot_dir, exist_ok=True)

# Define the thresholds you want to plot.
freq_thresholds = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
alphas = coef_df_full['alpha'].unique()

for current_alpha in alphas:
    # Filter dataset for the current alpha
    df_alpha = coef_df_full[coef_df_full['alpha'] == current_alpha].copy()
    
    for thresh in freq_thresholds:
        # Filter for the current threshold
        df_thresh = df_alpha[df_alpha['selection_freq'].fillna(0.0) > thresh].copy()
        robust_count = len(df_thresh)
        
        # Skip plotting if no terms survived the threshold
        if robust_count == 0:
            print(f"--- [SKIP] Alpha {current_alpha} | Threshold {thresh}: No features survived. ---", flush=True)
            continue
            
        # Compute absolute coefficient for sorting
        df_thresh['abs_coefficient'] = df_thresh['coefficient'].abs()

        # Take top 20 by absolute coefficient
        top20_df = df_thresh.sort_values(by='abs_coefficient', ascending=False).head(20).copy()

        # Map names and create labels
        top20_df['term_name'] = top20_df['annotation_id'].apply(get_term_name)
        top20_df['label'] = top20_df['term_name'] + ' (' + top20_df['annotation_id'] + ')'

        # Assign colors
        top20_df['Color'] = top20_df['coefficient'].apply(lambda x: '#ff9999' if x > 0 else '#99ccff')

        # Generate Barplot
        plt.figure(figsize=(14, 10))
        ax = sns.barplot(
            x='abs_coefficient', 
            y='label', 
            data=top20_df, 
            palette=top20_df['Color'].tolist()
        )

        ax.grid(axis='x', color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
        ax.set_axisbelow(True)

        max_width = top20_df['abs_coefficient'].max()
        text_offset = max_width * 0.015 if max_width > 0 else 0.01

        for i, patch in enumerate(ax.patches):
            annot_id = top20_df.iloc[i]['annotation_id']
            freq = top20_df.iloc[i]['selection_freq']
            
            patch.set_edgecolor('white')
            patch.set_linewidth(1.2)
            
            if 'DOID' in annot_id or 'HDO' in annot_id:
                patch.set_hatch('///')
                
            width = patch.get_width()
            y_pos = patch.get_y() + patch.get_height() / 2
            
            ax.text(width + text_offset, y_pos, f"Freq: {freq:.2f}", 
                    va='center', ha='left', fontsize=10, color='#333333', fontweight='bold')

        plt.title(f'Top 20 Most Predictive Annotations (Freq > {thresh:.2f})\nα = {current_alpha:.2f} | Non-Zero Terms Surviving = {robust_count}', fontsize=16, pad=20)
        plt.xlabel('Absolute Elastic Net Coefficient', fontsize=14)
        plt.ylabel('Annotation', fontsize=14)

        ax.set_xlim(0, max_width * 1.15)

        legend_elements = [
            Patch(facecolor='#ff9999', edgecolor='white', label='Positive Coeff. (GO)'),
            Patch(facecolor='#ff9999', edgecolor='white', hatch='///', label='Positive Coeff. (HDO)'),
            Patch(facecolor='#99ccff', edgecolor='white', label='Negative Coeff. (GO)'),
            Patch(facecolor='#99ccff', edgecolor='white', hatch='///', label='Negative Coeff. (HDO)')
        ]
        plt.legend(handles=legend_elements, loc='lower right', fontsize=12)

        plt.tight_layout()
        
        # Save dynamically inside the directory as PNG
        plot_filename = f"alpha_{current_alpha:.2f}_freq_{thresh:.2f}.png"
        plot_filepath = os.path.join(top20_plot_dir, plot_filename)
        plt.savefig(plot_filepath, dpi=300)
        plt.close()

print("--- [COMPLETE] MIXED ELASTIC NET MODEL VISUALIZATIONS FINISHED SUCCESSFULLY ---", flush=True)