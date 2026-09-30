import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import time

network_file = snakemake.input.network
output_plot = snakemake.output.plot
aspect = snakemake.wildcards.aspect

start_time = time.time()
print("--- Starting Network Growth Plotting ---", flush=True)

# 1. Load Network
print(f"[LOAD] Reading network from {network_file}...", flush=True)
with open(network_file, 'rb') as f:
    G_nx = pickle.load(f)

# 2. Extract Dates (Excluding self-loops to match previous logic)
dates = []
for u, v, data in G_nx.edges(data=True):
    if str(u) == str(v): 
        continue
        
    date_str = data.get('discovery_date')
    if date_str is not None:
        dates.append(int(date_str))

print(f"[PROCESS] Extracted {len(dates)} valid edges with discovery dates.", flush=True)

# 3. Aggregate and Calculate Cumulative Sum
date_counts = Counter(dates)
df = pd.DataFrame(list(date_counts.items()), columns=['Date', 'New_Edges'])
df = df.sort_values('Date').reset_index(drop=True)

# The total number of edges present in the network up to each date
df['Cumulative_Edges'] = df['New_Edges'].cumsum()

# Format dates for plotting
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')

# 4. Plotting Setup
print("[PLOT] Drawing the visualization...", flush=True)
plt.figure(figsize=(12, 6))
sns.set_context("talk")
sns.set_style("whitegrid")

# Cumulative Lineplot
sns.lineplot(
    data=df, 
    x='Date', 
    y='Cumulative_Edges', 
    linewidth=3, 
    color='royalblue'
)

# Formatting
plt.title(f"GO {aspect} Network Growth Over Time\n(Total Edges Present)", fontsize=20, fontweight='bold', pad=20)
plt.xlabel("Date", fontsize=16, fontweight='bold', labelpad=15)
plt.ylabel("Total Number of Edges", fontsize=16, fontweight='bold', labelpad=15)

# Axis limits and ticks
plt.xlim(df['Date'].min(), df['Date'].max())
plt.ylim(bottom=0)
plt.xticks(rotation=45, fontsize=12)
plt.yticks(fontsize=12)

sns.despine()
plt.tight_layout()

# 5. Save
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- Plot successfully generated in {elapsed} seconds! ---", flush=True)