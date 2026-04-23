import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pronto import Ontology

aspect = snakemake.wildcards.aspect
depth = snakemake.wildcards.depth
cutoff = snakemake.wildcards.cutoff
input_coefficients = snakemake.input.coefficients
ontology = snakemake.input.ontology
outputplot = snakemake.output.top_abs_coeffs_plot

print(f"Plotting top abs value coefficients for GO {aspect} depth {depth} cutoff {cutoff}...\n")

print(f"Loading data for GO {aspect} depth {depth} cutoff {cutoff}...")

coefficients = pd.read_csv(input_coefficients, sep = '\t')

print(f"GO {aspect} depth {depth} cutoff {cutoff} data loaded!\n")

print(f"Loading GO ontology from local file...")

go = Ontology(ontology)

print(f"GO ontology loaded!\n")

coefficients['abs_coefficient'] = coefficients['Coefficient'].abs()
id_to_name = {term.id: term.name for term in go.terms()}
coefficients['GO_term'] = coefficients['GO_id'].map(id_to_name)
top_coefs = coefficients.sort_values('abs_coefficient', ascending = False).head(20)
top_coefs['GO_term_short'] = top_coefs['GO_term'].apply(
    lambda x: x[:30] + '...' if len(str(x)) > 30 else x
)

plt.figure(figsize = (10, 14))
plot = sns.barplot(
    data = top_coefs,
    x = 'GO_term_short',
    y = 'abs_coefficient',
    hue = 'Coefficient',
    palette = 'RdBu_r',
    edgecolor = '.2',
    width = 0.8,
    linewidth = 0.8
)
plt.grid(axis = 'y', linestyle = '--', alpha = 0.6)
plt.axhline(0, color = 'black', linewidth = 1)
plt.xticks(
    rotation=45, 
    ha='right', 
    fontsize=14, 
    rotation_mode='anchor'
)
plt.yticks(fontsize=14)
plt.title(f'Top 20 GO {aspect} Depth {depth} Coefficients by Abs Value (Cutoff = {cutoff})', fontsize = 22, pad = 25)
plt.xlabel('GO Term', fontsize = 18, labelpad = 15)
plt.ylabel('Absolute Elastic Net Coefficient Value', fontsize = 18, labelpad = 15)
plt.tight_layout()

if plot.get_legend() is not None:
    plot.get_legend().remove()

print(f"Saving GO {aspect} depth {depth} cutoff {cutoff} top coefficients by abs value plot...\n")

plt.xlim(-0.5, len(top_coefs) - 0.5)
plt.savefig(outputplot, dpi = 300, bbox_inches = 'tight')
plt.close()

print(f"Top abs value coefficients plot for GO {aspect} depth {depth} cutoff {cutoff} ready!\n")