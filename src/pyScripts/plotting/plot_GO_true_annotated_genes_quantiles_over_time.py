import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import time

input_file = snakemake.input.quantile_file
output_file = snakemake.output.plot_file
term = snakemake.wildcards.term

start_time = time.time()
print(f"--- [{term}] Starting temporal quantile plotting ---")
print(f"[{term}] [LOAD] Reading data from {input_file}...")

df = pd.read_parquet(input_file)

print(f"[{term}] [INFO] Loaded {len(df)} rows of data.")
print(f"[{term}] [PROCESS] Formatting dates...")

# Date formating
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df = df.sort_values('Date')
unique_dates = df['Date'].nunique()

print(f"[{term}] [PLOT] Aggregating {unique_dates} unique dates for lineplot...")

plt.figure(figsize=(12, 6))

# sns.lineplot automatically aggregates the multiple Future_Gene quantiles per Date
# It plots the mean as a line and shades the 95% Confidence Interval
sns.lineplot(
    data=df, 
    x='Date', 
    y='Quantile', 
    marker='o', 
    linewidth=2,
    color='royalblue',
    errorbar=('ci', 95)
)

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