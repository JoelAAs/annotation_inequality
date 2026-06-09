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
    go_id_safe = os.path.basename(dist_file).replace('_mean_distance_from_annotated.csv')
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

    # Cronological sorting
    df_future = df_future.sort_values(by = 'time_to_annot_days')

    # Save raw sorted arrays to draw the true curve for the GO annot
    lines_data.append({
        'x': df_future['time_to_annot_days'].values,
        'y': df_future['Mean_Probability_From_Annotated'].values
    })

    # Save raw dataframe for the later global trendline
    all_data.append(df_future[['time_to_annot_days', 'Mean_Probability_From_Annotated']])

# If no valid data is found create an empty plot and exit the script
if not all_data:
    print(f"\nWARNING: No valid data found! Exiting...")

    plt.figure()
    plt.text(0.5, 0.5, 'No valid data plot', ha = 'center', va = 'center')
    plt.savefig(outputplot)
    sys.exit(0)

# PLOTTING
# Instead if valid data is found:
