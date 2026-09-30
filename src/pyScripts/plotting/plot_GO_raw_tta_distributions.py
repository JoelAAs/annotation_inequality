import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import time

input_file = snakemake.input.raw_tta_file
output_plot = snakemake.output.plot_file
term = snakemake.wildcards.term
term_formatted = term.replace('_', ':')
aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff

start_time = time.time()
print(f"--- Starting Presentation-Ready Density Plot for {term_formatted} ---", flush=True)

# Load Data
df = pd.read_parquet(input_file)

if df.empty:
    print(f"[{term_formatted}] No data available. Creating empty plot to satisfy Snakemake.")
    plt.figure(figsize=(8, 6))
    plt.text(0.5, 0.5, 'No valid future annotations', ha='center', va='center', fontsize=20)
    plt.axis('off')
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    exit(0)

# Convert Years to Months
df['TTA_Months'] = df['TTA_Years'] * 12
median_tta = df['TTA_Months'].median()

# Plotting Setup (High Contrast, Large Canvas)
plt.figure(figsize=(12, 7))
sns.set_context("talk") 
sns.set_style("white")  

# Generate the Density Curve (KDE)
ax = sns.kdeplot(
    data=df,
    x='TTA_Months',
    fill=True,
    color="#2b8cbe", 
    alpha=0.5,
    linewidth=3.5    
)

# Add Median Anchor Line
plt.axvline(x=median_tta, color='#de2d26', linestyle='--', linewidth=2.5, zorder=3)
plt.text(median_tta + 0.5, ax.get_ylim()[1] * 0.85, f'Median: {median_tta:.1f} mo', 
         color='#de2d26', fontsize=16, fontweight='bold')

# Formatting
plt.title(f"Time-to-Annotation Density: {term_formatted}", 
          fontsize=24, fontweight='bold', pad=25)

plt.suptitle(f"Aspect: {aspect.upper()} | Depth: {depth} | Cutoff: {cutoff}", 
             fontsize=14, color='gray', y=0.92)

plt.xlabel("Time to Annotation ($\\Delta$T in Months)", fontsize=18, labelpad=15, fontweight='bold')
plt.ylabel("Density", fontsize=18, labelpad=15, fontweight='bold')

# Axis limits
plt.xlim(left=0)
plt.ylim(bottom=0)

# Ticks and Grid
plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(6))
plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(1))

# Rotate labels 90 degrees vertically to fit the 6-month scale perfectly
plt.xticks(fontsize=12, rotation=90)
plt.yticks(fontsize=14)

# Clean up spines 
sns.despine()

# Vertical gridlines only, aligned with the 6-month major ticks
plt.grid(True, axis='x', linestyle='--', alpha=0.4, color='gray')
plt.grid(False, axis='y')

# Save
os.makedirs(os.path.dirname(output_plot), exist_ok=True)
plt.savefig(output_plot, dpi=300, bbox_inches='tight', transparent=False)
plt.close()

elapsed = round(time.time() - start_time, 2)
print(f"--- Plot successfully generated in {elapsed} seconds! ---", flush=True)