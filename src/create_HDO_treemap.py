import pandas as pd
import plotly.express as px
from goatools.obo_parser import GODag
import plotly.io as pio
import os

os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Snakemake variables
cutoff = snakemake.wildcards.cutoff
input_csv = snakemake.input.all_coefficients
obo_file = snakemake.input.ontology
output_pdf = snakemake.output.treemap
output_html = snakemake.output.treemap_html

print(f"Processing HDO cutoff {cutoff} visualization...")

print(f"Loading HDO ontology...")
# GODag works for doid.obo just as it does for go-basic.obo
hdo_dag = GODag(obo_file)
max_hdo_depth = max(term.level for term in hdo_dag.values())

print(f"Loading HDO cutoff {cutoff} data...")
# Note: Using 'HDO_doid' as the ID column per your previous HDO rules
df = pd.read_csv(input_csv, sep='\t').drop_duplicates(subset=['HDO_doid'])
target_col = 'Coefficient' 

def get_longest_path_to_root(term):
    """Trace the longest hierarchical path to ensure deep terms show their full depth."""
    path = []
    curr = term
    while curr:
        path.append(curr.name.replace("/", "-"))
        if not curr.parents: break
        # Using MAX level to force the deepest possible hierarchy
        curr = max(curr.parents, key=lambda x: x.level)
    return "/".join(reversed(path))

print(f"Processing all HDO coefficients...")
plot_data = []
for _, row in df.iterrows():
    disease_id = row['HDO_doid']
    if disease_id in hdo_dag:
        term = hdo_dag[disease_id]
        hier_path = get_longest_path_to_root(term)
        
        # Area Scaling: Root = Big Area, Leaf = Small Area
        area_weight = (max_hdo_depth - term.level + 1) ** 2
        
        # Label format: "Disease Name (DOID:ID) [D:X]"
        leaf_label = f"{term.name.replace('/', '-')} ({disease_id}) [D:{term.level}]"
        
        if "/" in hier_path:
            full_path_with_depth = f"{hier_path.rsplit('/', 1)[0]}/{leaf_label}"
        else:
            full_path_with_depth = leaf_label

        plot_data.append({
            "full_path": full_path_with_depth,
            "depth": term.level,
            "area_size": area_weight, 
            "coefficient": row[target_col],
            "abs_coeff": abs(row[target_col])
        })

plot_df = pd.DataFrame(plot_data)

# Color Balancing: 98th percentile
v_max = plot_df['coefficient'].abs().quantile(0.98)

print(f"Creating HDO cutoff {cutoff} treemap for {len(plot_df)} diseases...")

if not plot_df.empty:
    fig = px.treemap(
        plot_df,
        path=['full_path'], 
        values='area_size',
        color='coefficient',
        hover_data={'depth': True, 'coefficient': ':.6f'},
        color_continuous_scale='RdBu_r',
        color_continuous_midpoint=0,
        range_color=[-v_max, v_max],
        title=f"Full HDO Landscape | Area = Depth (Root=Large) | Cutoff = {cutoff} | N = {len(plot_df)}"
    )
    
    fig.update_layout(
        width=3000, 
        height=2000, 
        margin=dict(t=250, l=40, r=150, b=40), 
        
        
        title=dict(
            text=f"Full HDO Landscape | Area = Depth (Root=Large) | Cutoff = {cutoff} | N = {len(plot_df)}",
            x=0.5,
            y=0.95,
            xanchor='center',
            yanchor='top',
            font=dict(size=80, color='black', family="Arial, sans-serif")
        ),

        
        coloraxis_colorbar=dict(
            title="<b>Coefficient Value</b>",
            title_font=dict(size=50), 
            tickfont=dict(size=40),  
            thicknessmode="pixels", 
            thickness=120,            
            lenmode="fraction", 
            len=0.7,
            yanchor="middle",
            y=0.5,
            x=1.02,                  
            ticks="outside",
            ticklen=15                  
        ),

        
        font=dict(size=30) 
    )
    
    fig.update_layout(uniformtext=dict(minsize=6, mode='hide'))

    print(f"Saving HTML to {output_html}...")
    fig.write_html(output_html)
    
    print(f"Saving PDF to {output_pdf}...")
    try:
        fig.write_image(output_pdf, engine="kaleido", scale=3)
    except Exception as e:
        print(f"PDF export failed: {e}")

print(f"HDO cutoff {cutoff} processing done!")