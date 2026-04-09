import pandas as pd
import pronto
import matplotlib.pyplot as plt
from collections import Counter

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
obo_path = snakemake.input.ontology
input_coefficients = snakemake.input.single_depth_elastic_net_coefficients
output_top = snakemake.output.top_coefficients
output_distribution = snakemake.output.coefficients_distribution

print(f"Processing complete GO {aspect} depth {depth} coefficients plotting...\n")

print("Loading data...")

df = pd.read_csv(input_coefficients, sep = '\t')

print("Data loaded!\n")

print(f"Loading GO Ontology from local file: {obo_path}")

ontology = pronto.Ontology(obo_path)

print("GO ontology loaded!\n")

print(f"Retrieving associations between GO {aspect} ids and GO {aspect} Terms...")

goid_to_name = {}
for term in ontology.terms():
    goid_to_name[str(term.id)] = str(term.name)
    # Also map alternate IDs so they don't show up as 'assigned'
    for alt in term.alternate_ids:
        goid_to_name[str(alt)] = str(term.name)

# 2. Clean the dataframe IDs
df['GO_id'] = df['GO_id'].astype(str).str.strip()

# 3. Map and check how many failed
df['GO_term'] = df['GO_id'].map(goid_to_name)

n_missing = df['GO_term'].isna().sum()
if n_missing > 0:
    sample_missing = df[df['GO_term'].isna()]['GO_id'].head(5).tolist()
    print(f"WARNING: {n_missing} IDs could not be found in the OBO file.")
    print(f"Sample missing IDs: {sample_missing}")

# 4. Final fill
df['GO_term'] = df['GO_term'].fillna('Unknown_GO_ID')

print("Associations retrieved!\n")

print(f"Plotting {aspect} depth {depth} Top coefficients...")

top_n = 20
df_pos = df[df['Coefficient'] > 0].sort_values(by = 'Coefficient', ascending = False).head(top_n)
df_neg = df[df['Coefficient'] < 0].sort_values(by = 'Coefficient').head(top_n)
df_plot = pd.concat([df_neg, df_pos])

df_plot['color'] = df_plot['Coefficient'].apply(lambda x: '#2ca02c' if x > 0 else '#d62728')

# Truncate very long names for the plot
df_plot['GO_term_short'] = df_plot['GO_term'].apply(lambda x: (x[:50] + '..') if len(x) > 50 else x)

plt.figure(figsize = (16 ,8))
plt.barh(df_plot['GO_term_short'], df_plot['Coefficient'], color = df_plot['color'], zorder = 0)
plt.axvline(0, color='black', linewidth=0.8)
plt.xlabel('Coefficient')
plt.title(f'Top GO {aspect} depth {depth} Elastic Net Coefficients (Green = Positive, Red = Negative)')
plt.gca().invert_yaxis()
plt.xticks(rotation = 0)
plt.yticks(rotation = 0)
plt.grid()
plt.tight_layout()
plt.savefig(output_top, dpi = 300)
plt.close()

print(f"Top {aspect} depth {depth} coefficients plotted!\n")

print(f"Plotting {aspect} depth {depth} coefficients distribution...\n")

coeff_values = df['Coefficient'].round(3)
counts = Counter(coeff_values)

values = list(counts.keys())
freqs = list(counts.values())
total_coefs = len(df)

plt.figure(figsize = (10, 6))
plt.hist(df['Coefficient'], bins = 30, color = 'skyblue', edgecolor = 'black')
plt.yscale('log')
plt.xlabel('Coefficient value')
plt.ylabel('Frequency (log-scale)')
plt.title(f'Histogram of depth {depth} GO {aspect} Elastic Net Coefficients Distribution')
plt.legend([f'Total coefficients = {total_coefs}'], loc='upper right')
plt.savefig(output_distribution, dpi = 300)
plt.close()

print(f"GO {aspect} depth {depth} Coefficient distribution plotted!\n")

print(f"Depth {depth} GO {aspect} coefficients plotting done!")