import pandas as pd
import pickle
import networkx as nx 
import os
import concurrent.futures
import gc

aspect = snakemake.wildcards.aspect
final_network = snakemake.input.final_network
top_5_pickle = snakemake.input.nodes_with_top_5_annotations_pickle
output_dir = snakemake.output.first_annotation_dates_dir
num_threads = snakemake.threads

print(f"Retrieving first annotation dates for top 5 GO {aspect} annotations...")

os.makedirs(output_dir, exist_ok = True)

# Empty global dictionary so all the workers recognize it for later
node_annotations = {}

# WORKER FUNCTION
def process_and_save_go(task_data):
    # Save the task's data into the proper variables
    go_id, gene_list, out_dir = task_data
    
    results = []

    print(f"\n--- [Core] Starting GO {aspect} id: {go_id} first annotation dates extraction...")
    
    for gene_id in gene_list:
        annotations_list = node_annotations.get(gene_id, [])
        
        for item in annotations_list:
            if item.get('go_id')== go_id:
                first_date = item.get('first_annotation_date')

                if first_date is not None:
                    results.append({
                        'gene_id': gene_id,
                        'GO_id': go_id,
                        'first_annotation_date': first_date
                    })

                # If we found the annotation in the node there is no point in analyzing it further so we break the loop
                break

    # If no valid dates are found then we skip the saving of the annotation
    if not results:
        return f"  [{go_id}] Skipped (no valid dates found)"
    
    df_go = pd.DataFrame(results)
    df_go = df_go.sort_values(by = 'first_annotation_date', ascending = True)

    safe_go_id = str(go_id).replace(':', '_')
    output_file = os.path.join(out_dir, f'{safe_go_id}_first_annotation_dates.csv')
    df_go.to_csv(output_file, sep = '\t', index = False)

    print(f"\n[{go_id}] First annotation dates extracted and saved!")
    
    return

# MAIN
if __name__ == '__main__':
    print(f"\nLoading GO {aspect} pickle files (final network and genes annotated with top 5 annotations)...\n")

    with open(final_network, 'rb') as f:
        G = pickle.load(f)
    with open(top_5_pickle, 'rb') as f:
        top_5_data = pickle.load(f)

    if isinstance(top_5_data, pd.DataFrame):
        top_5_dict = top_5_data.set_index('GO_id')['annotated_genes'].to_dict()
    else:
        top_5_dict = top_5_data

    print(f"\nExtracting GO {aspect} lightweight annotations and clearing graph from RAM...\n")

    node_annotations = {
        node: data.get('go_annotations', [])
        for node, data in G.nodes(data = True)
        if 'go_annotations' in data
    }

    # Memory wipe
    del G
    gc.collect()

    # Prepare parallel tasks
    tasks = []
    for go_id, gene_list in top_5_dict.items():
        tasks.append((go_id, gene_list, output_dir))

    print(f"\nLaunching GO {aspect} first annotation dates extraction and save on {num_threads} cores...\n")

    with concurrent.futures.ProcessPoolExecutor(max_workers = num_threads) as executor:
        list(executor.map(process_and_save_go, tasks))

    print(f"\n--- [Core] SUCCESS: All GO {aspect} first annotation dates CSVs saved to {output_dir}! ---")