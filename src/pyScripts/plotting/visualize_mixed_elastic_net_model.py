import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pronto

coef_path = snakemake.input.en_coefficients
go_obo = snakemake.input.go_obo
hdo_obo = snakemake.input.hdo_obo
out_file = snakemake.output.top20_plot

print("--- [START] TOP 20 COEFFICIENTS VISUALIZATION STARTED ---")

print("--- [RUNNING] LOADING ONTOLOGIES AND DATA ---")
go_ont = pronto.Ontology(go_obo)
hdo_ont = pronto.Ontology(hdo_obo)

df = pd.read_parquet(coef_path)
df = df[['annotation_id', 'coefficient']]

print(f"--- [INFO] LOADED {len(df)} NON-ZERO COEFFICIENTS ---")

print("--- [RUNNING] EXTRACTING TOP 20 FEATURES ---")
df['abs_coef'] = df['coefficient'].abs()
df = df.sort_values(by='abs_coef', ascending=False).head(20).copy()
df = df.sort_values(by='abs_coef', ascending=True)

print("--- [RUNNING] MAPPING IDS TO ONTOLOGY NAMES ---")
def get_term_name(term_id):
    if term_id in go_ont:
        return f"{go_ont[term_id].name} ({term_id})"
    elif term_id in hdo_ont:
        return f"{hdo_ont[term_id].name} ({term_id})"
    return term_id

df['term_name'] = df['annotation_id'].apply(get_term_name)
df['direction'] = df['coefficient'] > 0
df['ontology'] = df['annotation_id'].apply(lambda x: 'GO' if 'GO:' in x else 'HDO')

plot_title = os.path.basename(coef_path).replace("_coefficients.parquet", "")
print(f"--- [INFO] PLOT TITLE SET TO: {plot_title} ---")

print("--- [RUNNING] GENERATING AND SAVING BARPLOT ---")
os.makedirs(os.path.dirname(out_file), exist_ok=True)

plt.figure(figsize=(10, 8))
plt.gca().set_axisbelow(True)
plt.grid(axis='x', color='lightgray', linestyle='--', alpha=0.7)

color_map = {True: '#FEAFAB', False: "#D9F1F3"}
hatch_map = {'GO': '', 'HDO': '///'}

colors = df['direction'].map(color_map).tolist()
hatches = df['ontology'].map(hatch_map).tolist()

bars = plt.barh(df['term_name'], df['abs_coef'], color=colors, edgecolor='white', linewidth=1)

for bar, hatch in zip(bars, hatches):
    bar.set_hatch(hatch)
    if hatch:
        bar.set_edgecolor('dimgray')
        bar.set_linewidth(1)

plt.title(f"Top 20 Coefficients: {plot_title}")
plt.xlabel("Absolute Coefficient Magnitude")
plt.ylabel("")

legend_elements = [
    Patch(facecolor="#D9F1F3", edgecolor='none', label='Negative Coefficient'),
    Patch(facecolor='#FEAFAB', edgecolor='none', label='Positive Coefficient'),
    Patch(facecolor='lightgray', edgecolor='none', label='GO Term'),
    Patch(facecolor='lightgray', edgecolor='dimgray', hatch='///', label='HDO Term')
]
plt.legend(handles=legend_elements, loc='lower right', framealpha=0.9)

plt.tight_layout()
plt.savefig(out_file, bbox_inches='tight', dpi=300)

print(f"--- [INFO] PLOT SAVED AT: {out_file} ---")
print("--- [COMPLETE] TOP 20 COEFFICIENTS VISUALIZATION FINISHED SUCCESSFULLY ---")