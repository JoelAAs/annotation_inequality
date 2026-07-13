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
outputplot = snakemake.output.plot_file

print(f"\nGenerating adjacency vs time to annotation plots for GO {aspect} annotations...")

os.makedirs(os.path.dirname(outputplot), exist_ok = True)

distance_files = glob.glob(os.path.join(distances_dir, '*_mean_distances_from_annotated.csv'))

all_data = []
lines_data = []

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

    # Keep only the future annotated genes and filter out dates with number of unannotated genes < 4
    df_future = df_dist[df_dist['first_annotation_date'] > df_dist['Date']].copy()
    df_future = df_future.groupby('Date').filter(lambda x: len(x) >= 4).copy()
    if df_future.empty:
        print(f"\nWARNING: {go_id_safe} df after threshold application is empty! Skipping it..")
        continue

    # Time to annotation in days
    dt_current = pd.to_datetime(df_future['Date'].astype(str), format = '%Y%m%d')
    dt_annot = pd.to_datetime(df_future['first_annotation_date'].astype(str).str.split('.').str[0], format = '%Y%m%d')
    df_future['time_to_annot_days'] = (dt_annot - dt_current).dt.days

    # Variance check
    if df_future['Mean_Probability_From_Annotated'].nunique() <= 1 or df_future['time_to_annot_days'].nunique() <= 1:
        print(f"\nNo varying data for id {go_id_safe}! Spearman correlation cannot be computed!")
        continue

    # Cronological sorting mean probability, standard deviation, and sample size calculation for each day to normalize for days with low sample size
    df_curve = df_future.groupby('time_to_annot_days')['Mean_Probability_From_Annotated'].agg(['mean', 'std', 'count']).reset_index()
    df_curve = df_curve.sort_values(by = 'time_to_annot_days')
    df_curve['std'] = df_curve['std'].fillna(0)
    df_curve['se'] = df_curve['std'] / np.sqrt(df_curve['count'])

    SMOOTHING_WINDOW = 90
    df_curve['smoothed_y'] = df_curve['mean'].rolling(window = SMOOTHING_WINDOW, min_periods = 1).mean()
    df_curve['smoothed_se'] = df_curve['se'].rolling(window = SMOOTHING_WINDOW, min_periods = 1).mean()

    # 95% confidence intervals
    df_curve['ci_upper'] = df_curve['smoothed_y'] + (1.96 * df_curve['smoothed_se'])
    df_curve['ci_lower'] = df_curve['smoothed_y'] - (1.96 * df_curve['smoothed_se'])

    lines_data.append({
        'go_id': clean_go_id,
        'obs_count': df_curve['count'].sum(),
        'x': df_curve['time_to_annot_days'].values,
        'y': df_curve['smoothed_y'].values,
        'y_upper': df_curve['ci_upper'].values,
        'y_lower': df_curve['ci_lower'].values
    })

    # Save raw dataframe for the later global trendline
    all_data.append(df_future[['time_to_annot_days', 'Mean_Probability_From_Annotated']])

# If no valid data is found create an empty plot and exit the script
if not all_data:
    print(f"\nWARNING: No valid data found! Exiting...")

    plt.figure()
    plt.text(0.5, 0.5, 'No valid data plot', ha = 'center', va = 'center')
    plt.savefig(outputplot)
    plt.close()
    sys.exit(0)

# PLOTTING
# Instead if valid data is found:
num_lines = len(lines_data)
print(f"\nDrawing {num_lines} true Go {aspect} annotation trajectories...")

plt.figure(figsize = (14, 10))

colors = sns.color_palette('tab10', num_lines)

# Draw the true curved lines
for i, line in enumerate(lines_data):
    plt.plot(
        line['x'],
        line['y'],
        color = colors[i],
        alpha = 0.9,
        linewidth = 2.5,
        label = f"{line['go_id']} (N = {int(line['obs_count']):,})"
    )

    plt.fill_between(
        line['x'],
        line['y_lower'],
        line['y_upper'],
        color = colors[i],
        alpha = 0.15,
        edgecolor = None
    )

# Draw the global trend line over everything
df_master = pd.concat(all_data, ignore_index = True)
sns.regplot(
    data = df_master,
    x = 'time_to_annot_days',
    y = 'Mean_Probability_From_Annotated',
    scatter = False,
    line_kws = {'color': 'red', 'linewidth': 6, 'linestyle': '--', 'label': f'Global {aspect} Trend (Avg)'}
)

plt.xlabel('Time To Annotation (Days)', fontsize = 14)
plt.ylabel('Adjacency (Mean Probability From Annotated)', fontsize = 14)
plt.title(f'GO {aspect} Annotations Adjacency vs Time To Annotation', fontsize = 18)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)

# Legend adjustment kept as requested to prevent overlapping out of the box
plt.legend(fontsize=12, loc='best', framealpha=0.9)
plt.grid(True, linestyle = '--', alpha = 0.4)

plt.tight_layout()
plt.savefig(outputplot, dpi = 400, bbox_inches = 'tight')
plt.close()

print(f"\nSUCCESS: Trajectory plot of GO {aspect} annotations adjacency vs time to annotation ready and saved to {outputplot}!")