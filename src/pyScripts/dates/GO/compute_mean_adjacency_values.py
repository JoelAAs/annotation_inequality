import pandas as pd
import os
import glob
import concurrent.futures

aspect = snakemake.wildcards.aspect
matrix_dir = snakemake.input.matrix_dir
mean_adjacency_dir = snakemake.output.mean_adjacency_dir
num_threads = snakemake.threads

print(f"Computing GO {aspect} mean adjacency values")

os.makedirs(mean_adjacency_dir, exist_ok = True)

def process_go_folder(go_folder_path):
    go_id_safe = os.path.basename(go_folder_path)
    csv_files = glob.glob(os.path.join(go_folder_path, '*.csv'))
    
    print(f"\n--- [Core] Starting GO {aspect} id: {go_id_safe} mean adjacency calculations ---")
    
    results = []
    
    for csv_file in csv_files:
        date_str = os.path.basename(csv_file).replace('.csv', '')
        
        # Loading the temporal matrix, keeping into account the multi index
        df = pd.read_csv(csv_file, header = [0, 1], index_col = [0, 1])
        
        num_annotated = (df.index.get_level_values(0) == 'Annotated').sum()
        num_future = (df.index.get_level_values(0) == 'Future').sum()
        
        print(f"  [{go_id_safe}] Processing date {date_str} | Annotated: {num_annotated} | Future {num_future}...")        
        
        # Extracting the various quadrants to calculate the adjacency means
        # Annotated -> Annotated quadrant
        try:
            annotated_quadrant = df.loc['Annotated', 'Annotated']
            mean_annotated = annotated_quadrant.values.mean() if not annotated_quadrant.empty else 0.0
        except KeyError:
            # If there are NO 'Annotated' genes on this date:
            mean_annotated = 0.0
            
        # Future -> Future quadrant
        try:
            future_quadrant = df.loc['Future', 'Future']
            mean_future = future_quadrant.values.mean() if not future_quadrant.empty else 0.0
        except KeyError:
            # If there are NO 'Future' genes on this date:
            mean_future = 0.0
            
        results.append({
            'Date': date_str,
            'Mean_Adjacency_Annotated': mean_annotated,
            'Mean_Adjacency_Future': mean_future
        })
        
    if results:
        # Convert results to pandas DataFrame, sorting chronologically by date
        df_results = pd.DataFrame(results).sort_values(by = 'Date')
        
        # Save the time series for this GO term
        output_file = os.path.join(mean_adjacency_dir, f"{go_id_safe}_mean_adjacency.csv")
        df_results.to_csv(output_file, sep = '\t', index = False)
        
    return f"{go_id_safe} mean adjacency values completed!"

# MULTIPROCESSING CONTROL
if __name__ == '__main__':
    # Find all GO terms subdirectories inside the matrix directory
    tasks = [f.path for f in os.scandir(matrix_dir) if f.is_dir()]
    
    print(f"\nLaunching GO {aspect} mean adjacency parallel processing with {num_threads} cores...\n")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers = num_threads) as executor:
        list(executor.map(process_go_folder, tasks))
        
    print(f"--- [Core] SUCCESS: All GO {aspect} mean adjacencies saved to {mean_adjacency_dir}! ---")