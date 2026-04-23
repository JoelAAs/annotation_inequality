import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

depth = snakemake.wildcards.depth
inputfile = snakemake.input.cutoff_file
outputplot = snakemake.output.cutoff_losing_plot

print(f"Plotting depth {depth} DOIDs remaining...")

print("Loading data...")
df = pd.read_csv(inputfile, sep=':')

if df.empty:
    plt.figure()
    plt.text(0.5, 0.5, f'No data available for Depth {depth}', ha='center')
    plt.savefig(outputplot)
    sys.exit(0)

print(f"Data loaded for depth {depth}!\n")

print(f"Computing depth {depth} metrics...")

# 1. Calculate the 'lost' percentage
df['removed_percentage'] = 100 - df['remaining_percentage']

# 2. Reconstruct the total and remaining counts
# Use a helper to avoid division by zero
def calc_remaining(row):
    if row['removed_percentage'] == 0:
        # If 0% removed, we can't mathematically derive the total from 'removed_doids'
        # But we know 100% remained. We'll handle this in the next step.
        return np.nan
    total = (row['removed_doids'] / row['removed_percentage']) * 100
    return total - row['removed_doids']

df['doids_left'] = df.apply(calc_remaining, axis=1)

# 3. Handle cases where 0 DOIDs were removed (100% remained)
# We find the total by looking at other rows where some WERE removed
valid_totals = (df['removed_doids'] / df['removed_percentage'].replace(0, np.nan)) * 100
total_doids_at_depth = int(round(valid_totals.median())) if not valid_totals.dropna().empty else 0

# Fill in the 100% remaining cases
df.loc[df['remaining_percentage'] == 100, 'doids_left'] = total_doids_at_depth
df['doids_left'] = df['doids_left'].fillna(0).round().astype(int)

print(f"Depth {depth} metrics ready!\n")

print(f"Plotting depth {depth}...")

plt.figure(figsize=(14, 8))

plt.plot(df['cutoff'], df['doids_left'], 
         marker='o', linestyle='-', color='#2c3e50', 
         linewidth=4, markersize=13,
         label=f'Derived Total DOIDs: {total_doids_at_depth}')

for i, row in df.iterrows():
    plt.annotate(
        f"{row['remaining_percentage']}%", 
        (row['cutoff'], row['doids_left']),
        textcoords="offset points", 
        xytext=(8, 8), 
        rotation=45, ha='left', va='bottom', fontsize=18, fontweight='bold', color='#e67e22')

plt.xticks(df['cutoff'], fontsize = 20)
plt.yticks(fontsize=20)
plt.ylim(0, (df['doids_left'].max() * 1.2) + 10) 
plt.grid(True, linestyle='--', alpha=0.8)
plt.margins(x=0.1)
plt.title(f'DOIDs Remaining by Cutoff Frequency (Depth {depth})', fontsize=26, fontweight='bold', pad=25)
plt.xlabel('Cutoff (Minimum Gene Count per DOID)', fontsize=22, labelpad=20)
plt.ylabel('Number of DOIDs in Matrix', fontsize=22, labelpad=20)
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(outputplot, dpi=300, bbox_inches='tight')

print(f"Depth {depth} doids lost per cutoff ready!\n")