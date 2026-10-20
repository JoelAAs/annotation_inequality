import pandas as pd
import pickle
import networkx as nx 
import os
import concurrent.futures
import gc

final_network = snakemake.input.final_network
top_5_pickle = snakemake.input.nodes_with_top_5_annotations
output_dir = snakemake.output.first_annotation_dates_dir
num_threads = snakemake.threads

print("Retrieving first annotation dates for top 5 HDO annotations...")

os.makedirs(output_dir, exist_ok=True)

# Empty global dictionary so all the workers recognize it for later
node_annotations = {}

# WORKER FUNCTION
def process_and_save_hdo(task_data):
    # Save the task's data into the proper variables
    hdo_id, gene_list, out_dir = task_data
    
    results = []

    print(f"\n--- [Core] Starting HDO id: {hdo_id} first annotation dates extraction...")
    
    for gene_id in gene_list:
        annotations_list = node_annotations.get(gene_id, [])
        
        for item in annotations_list:
            if item.get('doid') == hdo_id:
                first_date = item.get('first_annotation_date')

                if first_date is not None:
                    results.append({
                        'gene_id': gene_id,
                        'HDO_id': hdo_id,
                        'first_annotation_date': first_date
                    })
                break

    # If no valid dates are found, create an empty placeholder with correct headers
    if not results:
        df_hdo = pd.DataFrame(columns=['gene_id', 'HDO_id', 'first_annotation_date'])
        status_message = f"  [{hdo_id}] No valid dates found. Empty placeholder created."
    else:
        df_hdo = pd.DataFrame(results)
        df_hdo = df_hdo.sort_values(by='first_annotation_date', ascending=True)
        status_message = f"\n[{hdo_id}] First annotation dates extracted and saved!"

    safe_hdo_id = str(hdo_id).replace(':', '_')
    output_file = os.path.join(out_dir, f'{safe_hdo_id}_first_annotation_dates.csv')
    df_hdo.to_csv(output_file, sep='\t', index=False)

    print(status_message)
    
    return

# MAIN
if __name__ == '__main__':
    print("\nLoading HDO pickle files (final network and genes annotated with top 5 annotations)...\n")

    with open(final_network, 'rb') as f:
        G = pickle.load(f)
    with open(top_5_pickle, 'rb') as f:
        top_5_data = pickle.load(f)

    if isinstance(top_5_data, pd.DataFrame):
        # Update column names to match HDO
        top_5_dict = top_5_data.set_index('HDO_doid')['annotated_genes'].to_dict()
    else:
        top_5_dict = top_5_data

    print("\nExtracting HDO lightweight annotations and clearing graph from RAM...\n")

    # Change 'go_annotations' to 'do_annotations' for HDO
    node_annotations = {
        node: data.get('do_annotations', [])
        for node, data in G.nodes(data=True)
        if 'do_annotations' in data
    }

    # Memory wipe
    del G
    gc.collect()

    # Prepare parallel tasks
    tasks = []
    for hdo_id, gene_list in top_5_dict.items():
        tasks.append((hdo_id, gene_list, output_dir))

    print(f"\nLaunching HDO first annotation dates extraction and save on {num_threads} cores...\n")

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        list(executor.map(process_and_save_hdo, tasks))

    print(f"\n--- [Core] SUCCESS: All HDO first annotation dates CSVs saved to {output_dir}! ---")