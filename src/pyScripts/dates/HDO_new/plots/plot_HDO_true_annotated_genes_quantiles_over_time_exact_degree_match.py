import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os

input_file = snakemake.input.quantile_file
output_file = snakemake.output.plot_file
term = snakemake.wildcards.term

start_time = time.time()
print(f"--- [{term}] Starting HDO temporal quantile plotting ---")
print(f"[{term}] [LOAD] Reading data from {input_file}...")

df = pd.read_parquet(input_file)

print(f"[{term}] [INFO] Loaded {len(df)} rows of data.")

os.makedirs(os.path.dirname(output_file), exist_ok=True)

# --- EMPTY INPUT PLACEHOLDER ---
if df.empty:
    print(f"--- [WARNING] [{term}] Input quantile file is empty. "
          f"Writing placeholder plot. ---", flush=True)

    plt.figure(figsize=(12, 6))
    plt.text(
        0.5, 0.5,
        f"No temporal matrices available\nTerm: {term}",
        ha='center', va='center', fontsize=16, color='gray',
        transform=plt.gca().transAxes
    )
    plt.title(f"Predictive Quantile of True Annotated Genes Over Time\nTerm: {term}",
              fontsize=14, pad=15)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    raise SystemExit(0)

print(f"[{term}] [PROCESS] Formatting dates...")

# Date formating
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')

plt.figure(figsize=(12, 6))

# sns.lineplot automatically aggregates the multiple Future_Gene quantiles per Date
sns.lineplot(
    data=df, 
    x='Date', 
    y='Quantile', 
    linewidth=1.8,
    color='royalblue'
)

plt.axhline(y=0.95, color='red', linestyle='--', linewidth=1.5, label='95th Percentile Threshold')
plt.legend(loc='lower right')

plt.title(f"Predictive Quantile of True Annotated Genes Over Time\nTerm: {term}", fontsize=14, pad=15)
plt.xlabel("Date", fontsize=12)
plt.ylabel("Quantile Score", fontsize=12)

# Fix Y-axis to 0-1 scale so all GO term plots are visually comparable
plt.ylim(0, 1.05) 

plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()

plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [{term}] Plot successfully generated in {elapsed_time} seconds! ---")