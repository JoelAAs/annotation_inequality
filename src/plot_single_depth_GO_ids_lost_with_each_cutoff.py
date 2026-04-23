import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
inputfile = snakemake.input.cutoff_file
outputplot = snakemake.output.cutoff_losing_plot

print(f"Plotting GO {aspect} depth {depth} GO ids remaining...")

print(f"Loading {aspect} depth {depth} data...")
df = pd.read_csv(inputfile, sep=':')

if df.empty:
    plt.figure()
    plt.text(0.5, 0.5, f'No data available for Depth {depth}', ha='center')
    plt.savefig(outputplot)
    sys.exit(0)

print(f"Data loaded for {aspect} depth {depth}!\n")

print(f"Computing GO {aspect} depth {depth} metrics...")

df['removed_percentage'] = 100 - df['remaining_percentage']

def calc_remaining(row):
    if row['removed_percentage'] == 0:
        return np.nan
    total = (row['removed_go_ids'] / row['removed_percentage']) * 100
    return total - row['removed_go_ids']

df['go_ids_left'] = df.apply(calc_remaining, axis=1)

valid_totals = (df['removed_go_ids'] / df['removed_percentage'].replace(0, np.nan)) * 100
total_goids_at_depth = int(round(valid_totals.median())) if not valid_totals.dropna().empty else 0

df.loc[df['remaining_percentage'] == 100, 'go_ids_left'] = total_goids_at_depth
df['go_ids_left'] = df['go_ids_left'].fillna(0).round().astype(int)

print(f"GO {aspect} Depth {depth} metrics ready!\n")

print(f"Plotting GO {aspect} depth {depth}...")

plt.figure(figsize=(14, 8))

plt.plot(df['cutoff'], df['go_ids_left'], 
         marker='o', linestyle='-', color='#2c3e50', 
         linewidth=4, markersize=13,
         label=f'Derived Total GO ids: {total_goids_at_depth}')

for i, row in df.iterrows():
    plt.annotate(
        f"{row['remaining_percentage']}%", 
        (row['cutoff'], row['go_ids_left']),
        textcoords="offset points", 
        xytext=(8, 8), 
        rotation=45, ha='left', va='bottom', fontsize=18, fontweight='bold', color='#e67e22')

plt.xticks(df['cutoff'], fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylim(0, (df['go_ids_left'].max() * 1.2) + 10) 
plt.grid(True, linestyle='--', alpha=0.8)
plt.margins(x=0.1)
plt.title(f'GO {aspect} ids Remaining by Cutoff Frequency (Depth {depth})', fontsize=26, fontweight='bold', pad=25)
plt.xlabel('Cutoff (Minimum Gene Count per GO id)', fontsize=22, labelpad=20)
plt.ylabel('Number of GO ids in Matrix', fontsize=22, labelpad=20)
plt.legend(fontsize=20)

plt.tight_layout()
plt.savefig(outputplot, dpi=300, bbox_inches='tight')

print(f"GO {aspect} Depth {depth} GO ids lost per cutoff ready!\n")