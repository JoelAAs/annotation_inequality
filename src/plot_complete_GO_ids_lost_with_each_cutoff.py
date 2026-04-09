import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

aspect = snakemake.wildcards.aspect
inputfile = snakemake.input.cutoff_file
outputplot = snakemake.output.cutoff_losing_plot

print(f"Plotting GO {aspect} ids lost...")

print(f"Loading data for complete GO {aspect}...")

df = pd.read_csv(inputfile, sep=':')

print(f"Data loaded for complete GO {aspect}!\n")

print(f"Computing complete GO {aspect} metrics...")

df['removed_percentage'] = 100 - df['remaining_percentage']

def calc_remaining(row):
    if row['removed_percentage'] == 0:
        return np.nan
    total = (row['removed_go_ids'] / row['removed_percentage']) * 100
    return total - row['removed_go_ids']

df['go_ids_left'] = df.apply(calc_remaining, axis = 1)

valid_totals = (df['removed_go_ids'] / df['removed_percentage'].replace(0, np.nan)) * 100
total_goids = int(round(valid_totals.median())) if not valid_totals.dropna().empty else 0

df.loc[df['remaining_percentage'] == 100, 'go_ids_left'] = total_goids
df['go_ids_left'] = df['go_ids_left'].fillna(0).round().astype(int)

print(f"Complete GO {aspect} metrics ready!\n")

plt.figure(figsize=(12, 7))

plt.plot(df['cutoff'], df['go_ids_left'], 
         marker='o', linestyle='-', color='#2c3e50', 
         linewidth=2, markersize=8, 
         label=f'Initial GO ids: {total_goids}')

plt.xticks(df['cutoff'])

for i, row in df.iterrows():
    plt.annotate(
        f"{row['remaining_percentage']}%", 
        (row['cutoff'], row['go_ids_left']),
        textcoords="offset points", 
        xytext=(0, 12), 
        ha='center', 
        fontsize=10,
        fontweight='bold',
        color='#e67e22'
    )

plt.legend(loc='upper right', fontsize=12, frameon=True, shadow=True)

plt.title(f'Number of GO {aspect} ids Remaining by Cutoff Frequency', fontsize=15, pad=20)
plt.xlabel('Cutoff (Minimum Gene Count)', fontsize=12, labelpad=10)
plt.ylabel(f'Count of GO {aspect} ids in Feature Matrix (any depth)', fontsize=12, labelpad=10)
plt.grid(True, linestyle='--', alpha=0.6)

plt.ylim(0, max(df['go_ids_left']) + 300)

plt.tight_layout()
plt.savefig(outputplot)

print(f"GO {aspect} ids lost plot saved!")