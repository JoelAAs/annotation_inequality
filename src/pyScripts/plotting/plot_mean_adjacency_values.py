import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import glob
import concurrent.futures

aspect = snakemake.wildcards.aspect
input_dir = snakemake.input.mean_adjacency_dir
output_dir = snakemake.output.plots_dir
num_threads = snakemake.threads

print(f"Generating GO {aspect} mean adjacency value plots...")

os.makedirs(output_dir, exist_ok = True)

def plot_go_timeline(csv_file):
    filename = os.path.basename(csv_file)
    go_id_safe = filename.replace('_mean_adjacency.csv', '')
    
    print(f"\n--- [Core] Starting GO {aspect} id: {go_id_safe} mean adjacency values plot ---")
    
    df = pd.read_csv(csv_file, sep = '\t')
    
    # Converting the date strings into Datetime object for corret X axis scaling
    df['Date'] = pd.to_datetime(df['Date'].astype(str), format = '%Y%m%d')
    
    # PLOTTING SECTION
    fig, ax = plt.subplots(figsize = (12, 7))
    
    ax.plot(df['Date'], df['Mean_Adjacency_Annotated'],
            label = 'Annotated Genes Mean Adjacency Probability', color = '#1f77b4', linewidth = 3.0)
    
    ax.plot(df['Date'], df['Mean_Adjacency_Future'],
            label = 'Future Genes Mean Adjacency Probability', color = '#ff7f0e', linewidth = 3.0)
    
    clean_go_id = go_id_safe.replace('_', ':')

    ax.set_title(f'GO {aspect} Mean Adjacency Values Over Time | ID: {clean_go_id}', fontsize = 20, fontweight = 'bold', pad = 15)
    ax.set_xlabel('Timeline', fontsize = 14, fontweight = 'bold', labelpad = 10)
    ax.set_ylabel('Mean Adjacency Value (Log Scale)', fontsize = 14, fontweight = 'bold', labelpad = 10)
    ax.set_yscale('log')
    
    ax.xaxis.set_major_locator(mdates.YearLocator(base = 1))    
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    plt.setp(ax.get_xticklabels(), rotation = 90, ha = 'right')
    ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
    ax.grid(True, linestyle = ':', alpha = 0.6)
    ax.legend(fontsize = 14, loc = 'upper right')
    
    fig.tight_layout()
    output_file = os.path.join(output_dir, f'{go_id_safe}_mean_adjacency_plot.png')
    fig.savefig(output_file, dpi = 300, bbox_inches='tight')
    
    plt.close(fig)
    
    return f"  [{clean_go_id}] Mean adjacency values plot ready!"
    
# MULTIPROCESSING CONTROL
if __name__ == '__main__':
    csv_files = glob.glob(os.path.join(input_dir, '*.csv'))
    
    print(f"\nLaunching GO {aspect} mean adjacency value plotting with {num_threads} cores...")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers = num_threads) as executor:
        list(executor.map(plot_go_timeline, csv_files))
        
    print(f"--- [Core] SUCCESS: All GO {aspect} mean adjacency value plots saved to {output_dir}! ---")