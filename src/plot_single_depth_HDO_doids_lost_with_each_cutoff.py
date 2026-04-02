import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

depth = snakemake.wildcards.depth
inputfile = snakemake.input.cutoff_file
outputplot = snakemake.output.cutoff_losing_plot

print(f"Plotting depth {depth} doids lost...")

df = pd.read_csv(inputfile, sep=':')

# Check if the df is empty
if df.empty:
    print(f"Depth {depth} cutoff file is empty. Creating an empty output file and exiting...")
    plt.figure()
    plt.text(0.5, 0.5, 'No data available', ha = 'center')
    plt.savefig(outputplot)
    print("Done. Exiting...")
    sys.exit()

# Calculate the doids still left for each cutoff
df['doids_left'] = 0.0
mask = df['remaining_percentage']

print("Calculating remaining DOIDs...")

df['removed_percentage'] = 100 - df['remaining_percentage']

df['doids_left'] = (df['removed_doids'] * df['remaining_percentage']) / df['removed_percentage'].replace(0, np.nan)
df['doids_left'] = df['doids_left'].fillna(0).round().astype(int)

# Final safety check, skipping plot if there are no lost DOIDs, in that case skipping the plot
if df['doids_left'].sum() == 0 and df['removed_doids'].sum() == 0:
    print("Data exists but contains no lost DOIDs. Skipping plot...")

print("Remaining DOIDs computed!\n")

print("Computing the total DOIDs...")

total_calc_series = (df['removed_doids'] / df['removed_percentage'].replace(0, np.nan)) * 100
total_doids = int(round(total_calc_series.dropna().median()))

print("Total DOIDs computed!\n")

plt.figure(figsize=(12, 7))

plt.plot(df['cutoff'], df['doids_left'], 
         marker='o', linestyle='-', color='#2c3e50', 
         linewidth=2, markersize=8,
         label=f'Initial DOIDs: {total_doids}'
        )

plt.xticks(df['cutoff'])

for i, row in df.iterrows():
    plt.annotate(
        f"{row['remaining_percentage']}%", 
        (row['cutoff'], row['doids_left']),
        textcoords="offset points", 
        xytext=(0, 12), 
        ha='center', 
        fontsize=10,
        fontweight='bold',
        color='#e67e22'
    )

plt.legend(loc='upper right', fontsize=12, frameon=True, shadow=True)

plt.title(f'Number of DOIDs Remaining by Cutoff Frequency for depth {depth}', fontsize=15, pad=20)
plt.xlabel('Cutoff (Minimum Gene Count)', fontsize=12, labelpad=10)
plt.ylabel(f'Count of DOIDs in Feature Matrix depth {depth}', fontsize=12, labelpad=10)
plt.grid(True, linestyle='--', alpha=0.6)

plt.ylim(0, max(df['doids_left']) + 300)

plt.tight_layout()
plt.savefig(outputplot)

print(f"Depth {depth} Doids lost plot saved!")