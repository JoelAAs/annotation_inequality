import pandas as pd
import matplotlib.pyplot as plt

inputfile = snakemake.input.cutoff_file
outputplot = snakemake.output.cutoff_losing_plot

print("Plotting doids lost...")

df = pd.read_csv(inputfile, sep=':')

total_doids = 6335
df['doids_left'] = total_doids - df['removed_doids']

plt.figure(figsize=(12, 7))

plt.plot(df['cutoff'], df['doids_left'], 
         marker='o', linestyle='-', color='#2c3e50', 
         linewidth=2, markersize=8, 
         label=f'Initial DOIDs: {total_doids}')

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

plt.title('Number of DOIDs Remaining by Cutoff Frequency', fontsize=15, pad=20)
plt.xlabel('Cutoff (Minimum Gene Count)', fontsize=12, labelpad=10)
plt.ylabel('Count of DOIDs in Feature Matrix (any depth)', fontsize=12, labelpad=10)
plt.grid(True, linestyle='--', alpha=0.6)

plt.ylim(0, max(df['doids_left']) + 300)

plt.tight_layout()
plt.savefig(outputplot)

print("Doids lost plot saved!")