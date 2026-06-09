import pandas as pd
import numpy as np
import os
import glob
import concurrent.futures
from scipy.stats import spearmanr

aspect = snakemake.wildcards.aspect
distances_dir = snakemake.input.distances_dir
first_annotation_dates_dir = snakemake.input.first_annotation_dates_dir
outputfile = snakemake.output.correlation_file
num_threads = snakemake.threads

print(f"\nStarting Spearman correlation calculation between adjacency and time to annotation for GO {aspect}...\n")

os.makedirs(os.path.dirname(outputfile), exist_ok = True)

# WORKER FUNCTION
def process_go_correlation(go_id_safe):
    dist_file = os.path.join(distances_dir, f'{go_id_safe}_mean_distances_from_annotated.csv')
    annot_file = os.path.join(first_annotation_dates_dir, f'{go_id_safe}_first_annotation_dates.csv')

    if not os.path.exists(dist_file) or not os.path.exists(annot_file):
        print(f"\n[{go_id_safe}] Mean distances or first annotation file missing! Returning empty result!")

        return []
    
    df_dist = pd.read_csv(dist_file, sep = '\t')
    df_annot = pd.read_csv(annot_file, sep = '\t')

    print(f"\n--- [Core] Starting GO {aspect} id: {go_id_safe} Spearman correlation between adjacency and time to annotation calculation...")

    df_annot['gene_id'] = df_annot['gene_id'].astype(str)
    df_dist['Future_Gene'] = df_dist['Future_Gene'].astype(str)

    # Map the first annotation dates
    gene_to_date = df_annot.set_index('gene_id')['first_annotation_date'].to_dict()
    df_dist['first_annotation_date'] = df_dist['Future_Gene'].map(gene_to_date)
    df_dist = df_dist.dropna(subset = ['first_annotation_date'])

    # Exclude genes annotated till current date (included) and apply threshold for each date
    df_future = df_dist[df_dist['first_annotation_date'] > df_dist['Date']].copy()

    df_future = df_future.groupby('Date').filter(lambda x: len(x) >= 4).copy()

    if df_future.empty:
        print(f"[{go_id_safe}] Empty df after applying threshold!")

        return[]

    clean_go_id = go_id_safe.replace('_', ':')
    
    # Calculate time to annotation for the entire pool in days
    dt_current = pd.to_datetime(df_future['Date'].astype(str), format = '%Y%m%d')
    dt_annot = pd.to_datetime(df_future['first_annotation_date'].astype(str).str.split('.').str[0], format = '%Y%m%d')

    df_future['time_to_annot_days'] = (dt_annot - dt_current).dt.days

    # Variance check (we have actual varying data to correlate?)
    if df_future['Mean_Probability_From_Annotated'].nunique() <= 1 or df_future['time_to_annot_days'].nunique() <= 1:
        print(f"\n[{go_id_safe}] No varying data... Spearman correlation cannot be computed!")
        
        return []

    try:
        # Compute overall spearman correlation for the entire GO annotation's history
        rho, pval = spearmanr(df_future['Mean_Probability_From_Annotated'], df_future['time_to_annot_days'])

        if not np.isnan(rho):
            print(f"\n[{clean_go_id}] Spearman correlation ready!")

            return[{
                'GO_id': clean_go_id,
                'Spearman_rho': rho,
                'p_value': pval,
                'Total_Observations': len(df_future)
            }]

    except:
        print(f"\n[{go_id_safe}] Spearman rho is null!")

        pass

    return []

# MAIN
if __name__ == '__main__':
    distance_files = glob.glob(os.path.join(distances_dir, '*_mean_distances_from_annotated.csv'))
    go_terms_safe = [os.path.basename(f).replace('_mean_distances_from_annotated.csv', '') for f in distance_files]

    all_results = []

    print(f"\nLaunching Spearman correlation calculation between adjacency and time to annotation for GO {aspect} with {num_threads} cores...")

    with concurrent.futures.ProcessPoolExecutor(max_workers = num_threads) as executor:
        worker_outputs = list(executor.map(process_go_correlation, go_terms_safe))
        for output in worker_outputs:
            all_results.extend(output)

    if all_results:
        df_final = pd.DataFrame(all_results)
        df_final = df_final.sort_values(by = ['GO_id'])
        df_final.to_csv(outputfile, sep = '\t', index = False)

        print(f"\n--- [Core] SUCCESS: GO {aspect} Spearman correlation between adjacency and time to annotation ready and saved to {outputfile}! ---")

    else:
        # Fallback if all dates/terms fail the thresholds
        df_empty = pd.DataFrame(columns = [
            'GO_id',
            'Spearman_rho',
            'p_value',
            'Total_Observations'
        ])

        df_empty.to_csv(outputfile, sep = '\t', index = False)

        print(f"\n--- [Core] WARNING: No valid data met the thresholds. Empty file created to satisfy Snakemake! ---")