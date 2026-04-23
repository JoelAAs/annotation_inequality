import pandas as pd
import plotly.express as px
from goatools.obo_parser import GODag
import plotly.io as pio
import os

os.environ["QT_QPA_PLATFORM"] = "offscreen"

aspect = snakemake.wildcards.aspect
cutoff = snakemake.wildcards.cutoff
input_csv = snakemake.input.all_coefficients
obo_file = snakemake.input.ontology
output_pdf = snakemake.output.treemap
output_html = snakemake.output.treemap_html
aspect_name = snakemake.wildcards.aspect

print(f"Processing GO {aspect} cutoff {cutoff} visualization...")

print(f"Loading GO ontology...")
godag = GODag(obo_file)
max_ontology_depth = max(term.level for term in godag.values())

print(f"Loading GO {aspect} cutoff {cutoff} data...")
df = pd.read_csv(input_csv, sep='\t').drop_duplicates(subset=['GO_id'])
target_col = 'Coefficient' 

def get_shortest_path_to_root(term):
    path = []
    curr = term
    while curr:
        # Sanitize name to prevent extra slashes breaking the path logic
        path.append(curr.name.replace("/", "-"))
        if not curr.parents: break
        # Trace shortest path via levels
        curr = min(curr.parents, key=lambda x: x.level)
    return "/".join(reversed(path))

print(f"Processing all coefficients for GO {aspect}...")
plot_data = []
for _, row in df.iterrows():
    go_id = row['GO_id']
    if go_id in godag:
        term = godag[go_id]
        hier_path = get_shortest_path_to_root(term)
        
        # Area Scaling: Root (D:0) gets max area, Leaves (High D) get min area
        area_weight = (max_ontology_depth - term.level + 1) ** 2
        
        # Leaf label including ID and Depth
        leaf_label = f"{term.name.replace('/', '-')} ({go_id}) [D:{term.level}]"
        
        # Reconstruct path with the detailed leaf label
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

# Color Balancing: Use 98th percentile to ensure outliers don't wash out the map
v_max = plot_df['coefficient'].abs().quantile(0.98)

print(f"Creating GO {aspect} cutoff {cutoff} treemap for {len(plot_df)} terms...")

# Create Treemap
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
        title=f"Full {aspect_name} Landscape | Area = Depth (Root=Large) | Cutoff = {cutoff} | N = {len(plot_df)}"
    )
    
    fig.update_layout(
        width=3000, 
        height=2000, 
        margin=dict(t=250, l=40, r=150, b=40), 
        
        
        title=dict(
            text=f"Full {aspect_name} Landscape | Area = Depth (Root=Large) | Cutoff = {cutoff} | N = {len(plot_df)}",
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
    
    # Hide text in very small boxes to keep PDF from becoming a 'black blob'
    fig.update_layout(uniformtext=dict(minsize=6, mode='hide'))

    print(f"Saving HTML to {output_html}...")
    fig.write_html(output_html)
    
    print(f"Saving PDF to {output_pdf} (this may take a moment for large datasets)...")
    try:
        # scale=3 provides 9000px resolution equivalent for deep zooming
        fig.write_image(output_pdf, engine="kaleido", scale=3)
    except Exception as e:
        print(f"PDF export failed: {e}. Check if 'kaleido' is updated.")

print(f"GO {aspect} cutoff {cutoff} processing done!")