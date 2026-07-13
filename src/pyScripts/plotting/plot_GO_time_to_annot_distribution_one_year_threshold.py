import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import sys

aspect = snakemake.wildcards.aspect
distances_dir = snakemake.input.distances_dir
first_annotation_dates_dir = snakemake.input.first_annotation_dates_dir
correlation_file = snakemake.input.correlation_file
outputplot = snakemake.output.distribution_plot

print(f"\nGenerating Close vs Far distribution plots for GO {aspect} annotations...")

os.makedirs(os.path.dirname(outputplot), exist_ok = True)

distance_files = glob.glob(os.path.join(distances_dir, '*_mean_distances_from_annotated.csv'))

all_data = []
THRESHOLD_DAYS = 365 # 1 year threshold

# DATA AGGREGATION
for dist_file in distance_files:
    go_id_safe = os.path.basename(dist_file).replace('_mean_distances_from_annotated.csv', '')
    clean_go_id = go_id_safe.replace('_', ':')
    annot_file = os.path.join(first_annotation_dates_dir, f'{go_id_safe}_first_annotation_dates.csv')

    if not os.path.exists(annot_file):
        print(f"\nWARNING: File path {annot_file} doesn't exist! Skipping it...")
        continue

    df_dist = pd.read_csv(dist_file, sep = '\t')
    df_annot = pd.read_csv(annot_file, sep = '\t')

    # Safety string casting
    df_annot['gene_id'] = df_annot['gene_id'].astype(str)
    df_dist['Future_Gene'] = df_dist['Future_Gene'].astype(str)

    # Mapping the first annotation dates
    gene_to_date = df_annot.set_index('gene_id')['first_annotation_date'].to_dict()
    df_dist['first_annotation_date'] = df_dist['Future_Gene'].map(gene_to_date)
    df_dist = df_dist.dropna(subset = ['first_annotation_date'])

    # Keep only the future annotated genes
    df_future = df_dist[df_dist['first_annotation_date'] > df_dist['Date']].copy()
    df_future = df_future.groupby('Date').filter(lambda x: len(x) >= 4).copy()
    
    if df_future.empty:
        continue

    # Time to annotation in days
    dt_current = pd.to_datetime(df_future['Date'].astype(str), format = '%Y%m%d')
    dt_annot = pd.to_datetime(df_future['first_annotation_date'].astype(str).str.split('.').str[0], format = '%Y%m%d')
    df_future['time_to_annot_days'] = (dt_annot - dt_current).dt.days

    # Variance check
    if df_future['Mean_Probability_From_Annotated'].nunique() <= 1:
        continue

    # --- TAG THE DATA FOR THE DISTRIBUTIONS ---
    # New column categorizing the data into our two groups
    df_future['Time_Window'] = np.where(df_future['time_to_annot_days'] <= THRESHOLD_DAYS, '≤ 1 Year', '> 1 Year')
    df_future['GO_id'] = clean_go_id
    
    # We only need these three columns for the violin plot
    all_data.append(df_future[['GO_id', 'Mean_Probability_From_Annotated', 'Time_Window']])

# If no valid data is found create an empty plot and exit the script
if not all_data:
    print(f"\nWARNING: No valid data found! Exiting...")
    plt.figure()
    plt.text(0.5, 0.5, 'No valid data plot', ha = 'center', va = 'center')
    plt.savefig(outputplot)
    plt.close()
    sys.exit(0)

# Combine everything into one massive dataframe for Seaborn
df_master = pd.concat(all_data, ignore_index = True)

# Load pre-calculated stats
print(f"\nLoading statistical data from {correlation_file}...")
df_stats = pd.read_csv(correlation_file, sep='\t')

print(f"\nDrawing distribution plots...")
go_order = sorted(df_master['GO_id'].unique())

plt.figure(figsize = (14, 8))

# Draw the split violin plot
sns.violinplot(
    data = df_master,
    x = 'GO_id',
    y = 'Mean_Probability_From_Annotated',
    hue = 'Time_Window',
    hue_order = ['≤ 1 Year', '> 1 Year'],
    order = go_order,
    split = True,
    inner = 'quartile',
    palette = {'≤ 1 Year': '#ff6666', '> 1 Year': '#66b3ff'},
    cut = 0
)

# Y axis zoom
y_upper_limit = df_master['Mean_Probability_From_Annotated'].quantile(0.96)
plt.ylim(bottom=-y_upper_limit * 0.05, top=y_upper_limit * 1.5)

# Add fold change text annotations from stats file
custom_labels = []
for i, go_id in enumerate(go_order):
    # Lookup this specific GO term in pre-calculated stats file
    row = df_stats[df_stats['GO_id'] == go_id]
    
    if not row.empty:
        fc = row['Fold_Change'].values[0]
        # Create a string with a newline so the FC sits below the ID
        custom_labels.append(f"{go_id}\n(FC: {fc:.2f}x)")
    else:
        custom_labels.append(go_id)

plt.xlabel('GO Term', fontsize = 14)
plt.ylabel('Adjacency (Mean Probability From Annotated)', fontsize = 14)
plt.title(f'GO {aspect} Adjacency Distribution: ≤ 1 Year vs > 1 Year Prior to Annotation', fontsize = 18)
plt.xticks(ticks=range(len(go_order)), labels=custom_labels, fontsize=12)
plt.yticks(fontsize = 12)

plt.legend(title='Time to Annotation', fontsize = 12, title_fontsize = 13, loc = 'upper right')
plt.grid(True, axis = 'y', linestyle = '--', alpha = 0.5)

plt.tight_layout()
plt.savefig(outputplot, dpi = 400, bbox_inches = 'tight')
plt.close()

print(f"\nSUCCESS: Distribution plot ready and saved to {outputplot}!")