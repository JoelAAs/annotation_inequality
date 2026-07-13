import pandas as pd
import plotly.express as px
from goatools.obo_parser import GODag
import os

os.environ["QT_QPA_PLATFORM"] = "offscreen"

# Snakemake variables
cutoff = snakemake.wildcards.cutoff
input_csv = snakemake.input.all_coefficients
obo_file = snakemake.input.ontology

# Only keeping the PNG output
output_png = snakemake.output.scatterplot_png

print(f"Processing HDO cutoff {cutoff} scatter plot visualization...")

print(f"Loading HDO ontology...")
hdo_dag = GODag(obo_file)

print(f"Loading HDO cutoff {cutoff} data...")
df = pd.read_csv(input_csv, sep='\t').drop_duplicates(subset=['HDO_doid'])
target_col = 'Coefficient' 

print(f"Processing all HDO coefficients...")
plot_data = []

for _, row in df.iterrows():
    disease_id = row['HDO_doid']
    if disease_id in hdo_dag:
        term = hdo_dag[disease_id]
        
        plot_data.append({
            "HDO_doid": disease_id,
            "name": term.name,
            "depth": term.level,
            "coefficient": row[target_col]
        })

plot_df = pd.DataFrame(plot_data)

print(f"Creating HDO cutoff {cutoff} scatter plot for {len(plot_df)} diseases...")

if not plot_df.empty:
    fig = px.scatter(
        plot_df,
        x='depth',
        y='coefficient',
        hover_name='name',
        hover_data={'HDO_doid': True, 'depth': True, 'coefficient': ':.6f'},
        opacity=1.0, # Kept solid so no dots look faded
        labels={
            'depth': 'Ontology Depth',
            'coefficient': 'Coefficient Value'
        },
        title=f"HDO Coefficients by Depth | Cutoff = {cutoff} | N = {len(plot_df)}"
    )
    
    # Set a uniform solid color (Plotly blue) with a black outline
    fig.update_traces(
        marker=dict(size=12, color='#1f77b4', line=dict(width=1, color='black'))
    )
    
    fig.update_layout(
        width=1600,  
        height=1000, 
        margin=dict(t=150, l=80, r=150, b=80), 
        
        title=dict(
            x=0.5,
            y=0.95,
            xanchor='center',
            yanchor='top',
            font=dict(size=40, color='black', family="Arial, sans-serif")
        ),
        
        xaxis=dict(
            title_font=dict(size=24),
            tickfont=dict(size=18),
            dtick=1 
        ),
        yaxis=dict(
            title_font=dict(size=24),
            tickfont=dict(size=18),
            zeroline=True, 
            zerolinewidth=2,
            zerolinecolor='black'
        )
        # Note: coloraxis_colorbar has been completely removed
    )
    
    print(f"Saving PNG to {output_png}...")
    try:
        # scale=2 gives a crisp high-res PNG
        fig.write_image(output_png, engine="kaleido", scale=2)
    except Exception as e:
        print(f"PNG export failed: {e}. Check if 'kaleido' is updated.")

print(f"HDO cutoff {cutoff} processing done!")