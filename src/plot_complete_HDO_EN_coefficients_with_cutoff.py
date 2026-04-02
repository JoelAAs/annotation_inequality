import pandas as pd
import pronto
import matplotlib.pyplot as plt
from collections import Counter

cutoff = snakemake.wildcards.cutoff
input_coefficients = snakemake.input.complete_elastic_net_coefficients
output_top = snakemake.output.top_coefficients
output_distribution = snakemake.output.coefficients_distribution

print(f"Processing complete HDO coefficients with cutoff {cutoff} plotting...\n")

print("Loading data...")

df = pd.read_csv(input_coefficients, sep = '\t')

print("Data loaded!\n")

print("Retrieving associations between DOIDs and HDO Terms...")

do_url = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"
ontology = pronto.Ontology(do_url)
doid_to_name = {term.id: term.name for term in ontology.terms()}

df['HDO_term'] = df['HDO_doid'].map(doid_to_name)

# Handle 'No_doid' and missing values
df['HDO_term'] = df['HDO_term'].fillna('No_doid_assigned')

print("Associations retrieved!\n")

print(f"Plotting Top coefficients for cutoff {cutoff}...")

top_n = 20
df_pos = df[df['Coefficient'] > 0].sort_values(by = 'Coefficient', ascending = False).head(top_n)
df_neg = df[df['Coefficient'] < 0].sort_values(by = 'Coefficient').head(top_n)
df_plot = pd.concat([df_neg, df_pos])

df_plot['color'] = df_plot['Coefficient'].apply(lambda x: 'green' if x > 0 else 'red')

plt.figure(figsize = (16 ,8))
plt.barh(df_plot['HDO_term'], df_plot['Coefficient'], color = df_plot['color'])
plt.xlabel('Coefficient')
plt.title(f'Top complete HDO full Elastic Net Coefficients with cutoff {cutoff} (Green = Positive, Red = Negative)')
plt.gca().invert_yaxis()
plt.xticks(rotation = 0)
plt.yticks(rotation = 0)
plt.tight_layout()
plt.savefig(output_top, dpi = 300)
plt.close()

print(f"Cutoff {cutoff} Top coefficients plotted!\n")

print(f"Plotting coefficients distribution for cutoff {cutoff}...\n")

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
plt.title(f'Histogram of complete HDO full Elastic Net Coefficients')
plt.legend([f'Total coefficients = {total_coefs}'], loc='upper right')
plt.savefig(output_distribution, dpi = 300)
plt.close()

print(f"Cutoff {cutoff}Coefficient distribution plotted!\n")

print(f"Complete HDO coefficients with cutoff {cutoff} plotting done!")