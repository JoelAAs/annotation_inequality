import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pronto import Ontology

depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
input_coefficients = snakemake.input.coefficients
ontology = snakemake.input.ontology
outputplot = snakemake.output.top_abs_coeffs_plot

print(f"Plotting top abs value coefficients for HDO depth {depth} cutoff {cutoff}...\n")

print(f"Loading data for HDO depth {depth} cutoff {cutoff}...")

coefficients = pd.read_csv(input_coefficients, sep = '\t')

print(f"HDO depth {depth} cutoff {cutoff} data loaded!\n")

print(f"Loading HDO ontology from local file...")

do = Ontology(ontology)

print(f"HDO ontology loaded!\n")

coefficients['abs_coefficient'] = coefficients['Coefficient'].abs()
id_to_name = {term.id: term.name for term in do.terms()}
coefficients['HDO_term'] = coefficients['HDO_doid'].map(id_to_name)
top_coefs = coefficients.sort_values('abs_coefficient', ascending = False).head(50)
top_coefs['HDO_term_short'] = top_coefs['HDO_term'].apply(
    lambda x: x[:30] + '...' if len(str(x)) > 30 else x
)

plt.figure(figsize = (12, 6))
plot = sns.barplot(
    data = top_coefs,
    x = 'HDO_term_short',
    y = 'Coefficient',
    palette = 'vlag',
    hue = 'Coefficient',
    edgecolor = '.2',
    linewidth = 1.5
)
plt.grid(axis = 'y', linestyle = '--', alpha = 0.7)
plt.axhline(0, color = 'black', linewidth = 1)
plt.xticks(
    rotation = 45,
    ha = 'right',
    fontsize = 9
)
plt.title(f'Top 50 HDO Depth {depth} Coefficients by Abs Value (Cutoff = {cutoff})')
plt.xlabel('HDO Term')
plt.ylabel('Absolute Elastic Net Coefficient Value')
plt.tight_layout()

print(f"Saving HDO depth {depth} cutoff {cutoff} top coefficients by abs value plot...\n")

plt.savefig(outputplot)
plt.close()

print(f"Top abs value coefficients plot for HDO depth {depth} cutoff {cutoff} ready!\n")