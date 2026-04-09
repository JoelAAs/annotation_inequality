import pandas as pd
import pronto
import matplotlib.pyplot as plt
from collections import Counter

aspect = snakemake.wildcards.aspect
cutoff = snakemake.wildcards.cutoff
input_coefficients = snakemake.input.complete_elastic_net_coefficients
obo_path = snakemake.input.ontology 
output_top = snakemake.output.top_coefficients
output_distribution = snakemake.output.coefficients_distribution

print(f"Processing complete GO {aspect} coefficients with cutoff {cutoff} plotting...\n")

print("Loading data...")

df = pd.read_csv(input_coefficients, sep = '\t')

print("Data loaded!\n")

print(f"Loading GO Ontology from local file: {obo_path}")

ontology = pronto.Ontology(obo_path)

print("GO ontology loaded!\n")

print("Retrieving associations between GO ids and GO Terms...")

goid_to_name = {term.id: term.name for term in ontology.terms()}
df['GO_term'] = df['GO_id'].map(goid_to_name)

# Handle 'No_doid' and missing values
df['GO_term'] = df['GO_term'].fillna('No_go_id_assigned')

print("Associations retrieved!\n")

print(f"Plotting Top coefficients for {aspect} cutoff {cutoff}...")

top_n = 20
df_pos = df[df['Coefficient'] > 0].sort_values(by = 'Coefficient', ascending = False).head(top_n)
df_neg = df[df['Coefficient'] < 0].sort_values(by = 'Coefficient').head(top_n)
df_plot = pd.concat([df_neg, df_pos])

df_plot['color'] = df_plot['Coefficient'].apply(lambda x: 'green' if x > 0 else 'red')

plt.figure(figsize = (16 ,8))
plt.barh(df_plot['GO_term'], df_plot['Coefficient'], color = df_plot['color'], zorder = 3)
plt.xlabel('Coefficient')
plt.title(f'Top complete GO {aspect} Elastic Net Coefficients with cutoff {cutoff} (Green = Positive, Red = Negative)')
plt.gca().invert_yaxis()
plt.xticks(rotation = 0)
plt.yticks(rotation = 0)
plt.grid()
plt.tight_layout()
plt.savefig(output_top, dpi = 300)
plt.close()

print(f"{aspect} Cutoff {cutoff} Top coefficients plotted!\n")

print(f"Plotting {aspect} coefficients distribution for cutoff {cutoff}...\n")

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
plt.title(f'Histogram of complete GO {aspect} Elastic Net Coefficients')
plt.legend([f'Total coefficients = {total_coefs}'], loc='upper right')
plt.savefig(output_distribution, dpi = 300)
plt.close()

print(f"{aspect} Cutoff {cutoff} Coefficient distribution plotted!\n")

print(f"Complete GO {aspect} coefficients with cutoff {cutoff} plotting done!")