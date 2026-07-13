import pandas as pd
import os
import glob
import concurrent.futures

aspect = snakemake.wildcards.aspect
input_dir = snakemake.input.input_dir
output_dir = snakemake.output.distances_dir
num_threads = snakemake.threads

print(f"Computing GO {aspect} mean distance of future annotated genes from already annotated genes...")

os.makedirs(output_dir, exist_ok = True)

def process_go_folder(go_folder_path):
    go_id_safe = os.path.basename(go_folder_path)
    csv_files = glob.glob(os.path.join(go_folder_path, '*.csv'))
    
    print(f"\n--- [Core] Starting GO {aspect} id: {go_id_safe} mean distance of future annotated genes from already annotated genes calculation...")
    
    results = []
    
    for csv_file in csv_files:
        date_str = os.path.basename(csv_file).replace('.csv', '')
        
        # Load the temporal matrix keeping into account the multi index
        df = pd.read_csv(csv_file, header = [0, 1], index_col = [0, 1])
        
        num_annotated = (df.index.get_level_values(0) == 'Annotated').sum()
        num_future = (df.index.get_level_values(0) == 'Future').sum()
        
        print(f"  [{go_id_safe}] Processing date {date_str} | Annotated: {num_annotated} | Future: {num_future}")
        
        # Taking the Annotated -> Future quadrant as we want the probabilities of the paths starting from Annotated genes to Future ones
        # [Rows = Annotated genes, Columns = Future genes]
        try:
            annot_to_future_quadrant = df.loc['Annotated', 'Future']
            
            if not annot_to_future_quadrant.empty and num_annotated > 0:
                # Sum all the probabilities for each future gene (column wise) and divide by the total number of Annotated genes
                future_gene_sums = annot_to_future_quadrant.values.sum(axis = 0)
                future_gene_means = future_gene_sums / num_annotated
                
                # Get the Future gene ids
                future_gene_ids = annot_to_future_quadrant.columns.tolist()
                
                # Store the mean probability for each individual Future gene
                for gene_id, mean_prob in zip(future_gene_ids, future_gene_means):
                    results.append({
                        'Date': date_str,
                        'Future_Gene': gene_id,
                        'Mean_Probability_From_Annotated': mean_prob
                    })
        
        except KeyError:
            # In case there are no 'Annotated' or 'Future' genes in this date the date is skipped
            pass
        
    if results:
        df_results = pd.DataFrame(results).sort_values(by = ['Date', 'Future_Gene'])
        
        output_file = os.path.join(output_dir, f'{go_id_safe}_mean_distances_from_annotated.csv')
        df_results.to_csv(output_file, sep = '\t', index = False)
        
    return f'{go_id_safe} mean distances of future annotated genes from annotated genes completed!'

# MULTIPROCESSING CONTROL
if __name__ == '__main__':
    tasks = [f.path for f in os.scandir(input_dir) if f.is_dir()]
    
    print(f"\nStarting GO {aspect} mean distances of future annotated genes from annotated genes parallel processing with {num_threads} cores...\n")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers = num_threads) as executor:
        list(executor.map(process_go_folder, tasks))
        
    print(f"--- [Core] SUCCESS: All GO {aspect} mean distances of future annotated genes from annotated genes saved to {output_dir}! ---")