import pandas as pd
import pronto
import sys
import matplotlib.pyplot as plt
from collections import Counter

cutoff = snakemake.wildcards.cutoff
depth = snakemake.wildcards.depth
input_coefficients = snakemake.input.single_depth_elastic_net_coefficients
obo_path = snakemake.input.ontology 
output_top = snakemake.output.top_coefficients
output_distribution = snakemake.output.coefficients_distribution

print(f"Processing Depth {depth}, Cutoff {cutoff} Plotting...")

print(f"Loading data for depth {depth}, cutoff {cutoff}...")

df = pd.read_csv(input_coefficients, sep='\t')

# Ensure Coefficient is numeric and drop any NaNs
df['Coefficient'] = pd.to_numeric(df['Coefficient'], errors='coerce')
df = df.dropna(subset=['Coefficient'])

print(f"Data for depth {depth}, cutoff {cutoff} loaded!\n")

# Handle Empty Case (Map-Reduce Safety)
# If the Elastic Net produced no coefficients, create placeholder plots
if df.empty or (df['Coefficient'] == 0).all():
    print(f"No valid coefficients found. Creating empty placeholder plots.")
    for out in [output_top, output_distribution]:
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, f"No significant DOIDs found\n(Depth {depth}, Cutoff {cutoff})", 
                 ha='center', va='center', fontsize=12, color='gray', fontweight='bold')
        plt.axis('off')
        plt.savefig(out, dpi=300)
        plt.close()
    sys.exit(0)

print(f"Loading HDO Ontology from local file: {obo_path}...")
try:
    ontology = pronto.Ontology(obo_path)

    print("HDO Ontology loaded!\n")

    print(f"Assigning HDO terms for depth {depth}, cutoff {cutoff} to their respective doids...")
    doid_to_name = {term.id: term.name for term in ontology.terms()} 
    df['HDO_term'] = df['HDO_doid'].map(doid_to_name).fillna(df['HDO_doid'])
    print(f"Terms for depth {depth}, cutoff {cutoff} assigned!\n")
except Exception as e:
    print(f"Error parsing ontology file: {e}. Falling back to raw IDs.")
    df['HDO_term'] = df['HDO_doid']

# Handle specific labels
df['HDO_term'] = df['HDO_term'].replace('No_doid', 'No_doid_assigned')

# Plot Top Coefficients
print(f"Generating depth {depth}, cutoff {cutoff} Top Coefficients plot...")

top_n = 20
# Get top positive and top negative
df_pos = df[df['Coefficient'] > 0].sort_values('Coefficient', ascending=False).head(top_n)
df_neg = df[df['Coefficient'] < 0].sort_values('Coefficient', ascending=True).head(top_n)

# Combine and sort for a nice visual flow
df_plot = pd.concat([df_pos, df_neg]).sort_values('Coefficient', ascending=True)

plt.figure(figsize=(16, 10))
colors = ['green' if x > 0 else 'red' for x in df_plot['Coefficient']]
plt.barh(df_plot['HDO_term'], df_plot['Coefficient'], color=colors, alpha=0.8)

# Add a vertical line at zero for clarity
plt.axvline(0, color='black', linewidth=1, alpha=0.7)

plt.xlabel('Elastic Net Coefficient Value', fontsize=12)
plt.title(f'Top Influential DOIDs: Depth {depth}, Cutoff {cutoff}\n(Green = Positive Association, Red = Negative)', 
          fontsize=16, pad=20)

plt.grid(axis='x', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(output_top, dpi=300)
plt.close()

print(f"Depth {depth}, cutoff {cutoff} Top Coefficients plot ready!\n")

# Plot Coefficient Distribution
print(f"Generating depth {depth}, cutoff {cutoff} Distribution plot...")

total_coefs = len(df)
nonzero_coefs = len(df[df['Coefficient'] != 0])

plt.figure(figsize=(10, 6))
plt.hist(df['Coefficient'], bins=40, color='skyblue', edgecolor='black', alpha=0.7)
plt.yscale('log') 

plt.xlabel('Coefficient Value', fontsize=12)
plt.ylabel('Frequency (Log Scale)', fontsize=12)
plt.title(f'Distribution of Coefficients: Depth {depth}, Cutoff {cutoff}', fontsize=15)

plt.legend([f'Total Terms: {total_coefs}\nNon-zero: {nonzero_coefs}'], loc='upper right')

plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig(output_distribution, dpi=300)
plt.close()

print(f"Depth {depth}, cutoff {cutoff} Distribution plot ready!\n")

print(f"Depth {depth}, Cutoff {cutoff} plots saved successfully.")