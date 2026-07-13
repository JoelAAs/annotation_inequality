import pandas as pd
import plotly.express as px
from goatools.obo_parser import GODag
import os

os.environ["QT_QPA_PLATFORM"] = "offscreen"

aspect = snakemake.wildcards.aspect
cutoff = snakemake.wildcards.cutoff
input_csv = snakemake.input.all_coefficients
obo_file = snakemake.input.ontology

# Only keeping the PNG output
output_png = snakemake.output.scatterplot_png
aspect_name = snakemake.wildcards.aspect

print(f"Processing GO {aspect} cutoff {cutoff} scatter plot visualization...")

print(f"Loading GO ontology...")
godag = GODag(obo_file)

print(f"Loading GO {aspect} cutoff {cutoff} data...")
df = pd.read_csv(input_csv, sep='\t').drop_duplicates(subset=['GO_id'])
target_col = 'Coefficient' 

print(f"Processing all coefficients for GO {aspect}...")
plot_data = []

for _, row in df.iterrows():
    go_id = row['GO_id']
    if go_id in godag:
        term = godag[go_id]
        
        plot_data.append({
            "GO_id": go_id,
            "name": term.name,
            "depth": term.level,
            "coefficient": row[target_col]
        })

plot_df = pd.DataFrame(plot_data)

# Color Balancing: Use 98th percentile to ensure outliers don't wash out the color scale
v_max = plot_df['coefficient'].abs().quantile(0.98)

print(f"Creating GO {aspect} cutoff {cutoff} scatter plot for {len(plot_df)} terms...")

if not plot_df.empty:
    fig = px.scatter(
        plot_df,
        x='depth',
        y='coefficient',
        hover_name='name',
        hover_data={'GO_id': True, 'depth': True, 'coefficient': ':.6f'},
        color='coefficient',
        color_continuous_scale='RdBu_r',
        color_continuous_midpoint=0,
        range_color=[-v_max, v_max],
        opacity=0.6,
        labels={
            'depth': 'Ontology Depth',
            'coefficient': 'Coefficient Value'
        },
        title=f"GO {aspect_name} Coefficients by Depth | Cutoff = {cutoff} | N = {len(plot_df)}"
    )
    
    # Adjust point size and layout for better readability
    fig.update_traces(marker=dict(size=12, line=dict(width=1, color='DarkSlateGrey')))
    
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
        ),

        coloraxis_colorbar=dict(
            title="<b>Coefficient</b>",
            title_font=dict(size=24), 
            tickfont=dict(size=18),  
            thicknessmode="pixels", 
            thickness=40,            
            yanchor="middle",
            y=0.5,
            ticks="outside"
        )
    )

    print(f"Saving PNG to {output_png}...")
    try:
        fig.write_image(output_png, engine="kaleido", scale=2)
    except Exception as e:
        print(f"PNG export failed: {e}. Check if 'kaleido' is updated.")

print(f"GO {aspect} cutoff {cutoff} processing done!")