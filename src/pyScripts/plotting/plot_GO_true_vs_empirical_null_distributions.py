import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

input_file = snakemake.input.corrected_file
output_folder = snakemake.output.output_folder
aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff

print(f"\n--- [Core] Starting plotting for GO {aspect} depth {depth} cutoff {cutoff} true vs empirical null distributions... ---\n")

print(f"--- [Core] Loading the corrected master adjacency file... ---\n")

df = pd.read_csv(input_file, sep = '\t')
unique_gos = df['GO_ID'].unique()

os.makedirs(output_folder, exist_ok = True)

print(f"--- [Core] Starting plots generation... ---\n")

for target_go in unique_gos:
    go_data = df[df['GO_ID'] == target_go].copy()
    
    # Optional safety: skip GO terms if they have less than 2 nodes
    if len(go_data) <= 2:
        print(f"--- [Core] Skipping {target_go}: Not enough data points! ---")
        continue
    
    # Canvas setup
    fig, axes = plt.subplots(1, 2, figsize = (14, 6))
    
    # --- PARITY PLOT (Scatter - Log Scale) ---
    # We define a tiny epsilon so 0.0 values don't disappear on a log scale
    epsilon = 1e-6 
    
    sns.scatterplot(
        x = go_data['Null_Mean'] + epsilon,
        y = go_data['True_Adjacency'] + epsilon,
        alpha = 0.6,
        color = 'royalblue',
        ax = axes[0]
    )
    
    # Convert axes to logarithmic scale
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    
    # Calculate min and max for the diagonal line based on the log scale
    max_val = max(go_data['True_Adjacency'].max(), go_data['Null_Mean'].max()) + epsilon
    min_val = epsilon 
    
    axes[0].plot([min_val, max_val], [min_val, max_val], color = 'red', linestyle = '--', linewidth = 2, label = 'Equal Performance')
    
    axes[0].set_title(f'Parity Plot: True vs. Random (Log Scale)', fontsize = 13, pad = 10)
    axes[0].set_xlabel('Random Degree-Matched Annotations (Mean) + 1e-6')
    axes[0].set_ylabel('True Annotations + 1e-6')
    axes[0].legend()
    axes[0].grid(True, linestyle = ':', alpha = 0.6)

    # --- DISTRIBUTION SHIFT (Boxplot) ---
    melted_df = go_data.melt(
        id_vars = ['Future_Gene'], 
        value_vars = ['True_Adjacency', 'Null_Mean'],
        var_name = 'Annotation_Type', 
        value_name = 'Adjacency_Score'
    )

    melted_df['Annotation_Type'] = melted_df['Annotation_Type'].replace({
        'True_Adjacency': 'True GO Annotations',
        'Null_Mean': 'Random Degree-Matched'
    })

    sns.boxplot(
        data = melted_df, 
        x = 'Annotation_Type', 
        y = 'Adjacency_Score',
        hue = 'Annotation_Type',
        palette = ['#87CEEB', '#D3D3D3'],
        ax = axes[1],
        showfliers = False,
        legend = False
    )

    sns.stripplot(
        data = melted_df,
        x = 'Annotation_Type',
        y = 'Adjacency_Score',
        color = 'black',
        alpha = 0.4,
        size = 3,
        jitter = True,
        ax = axes[1]
    )

    axes[1].set_title(f'Distribution Shift', fontsize = 13, pad = 10)
    axes[1].set_ylabel('Transition Probability (Adjacency Score)')
    axes[1].set_xlabel('')
    axes[1].grid(True, linestyle = ':', alpha = 0.6)

    plt.suptitle(f'Topological Predictive Power: {target_go}', fontsize = 16, y = 1.05)
    plt.tight_layout()
    
    safe_filename = target_go.replace(':', '_')
    save_path = os.path.join(output_folder, f"{safe_filename}_comparison.png")
    
    plt.savefig(save_path, bbox_inches = 'tight', dpi = 150)
    plt.close(fig) 
    
    print(f"--- [Core] Saved plot for {target_go} ---\n")

print(f"\n --- [Core] All plots generated successfully! You can find them in: {output_folder} ---\n")