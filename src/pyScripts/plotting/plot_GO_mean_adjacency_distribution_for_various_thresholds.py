import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mannwhitneyu
from pathlib import Path
import time
import warnings

# Suppress seaborn layout warnings for clean logs
warnings.filterwarnings('ignore')

adj_dir = Path(snakemake.input.mean_adj_dir)
dates_dir = Path(snakemake.input.annot_dates_dir)
output_pdf = snakemake.output.plot_file
aspect = snakemake.wildcards.aspect

start_time = time.time()
print(f"--- [GLOBAL {aspect.upper()}] Generating Master Distribution PDF ---")

# LOAD AND AGGREGATE DATA
all_data = []
adj_files = list(adj_dir.glob("*_mean_adjacencies.parquet"))
print(f"Aggregating {len(adj_files)} files...")

for adj_file in adj_files:
    term_str = adj_file.name.replace("_mean_adjacencies.parquet", "")
    term_formatted = term_str.replace("_", ":")
    date_file = dates_dir / f"{term_str}_first_annotation_dates.csv"
    
    if not date_file.exists():
        continue
        
    adj_df = pd.read_parquet(adj_file)[['Date', 'Future_Gene', 'PID0_mean_adj']]
    dates_df = pd.read_csv(date_file, sep='\t')[['gene_id', 'first_annotation_date']]
    
    adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
    dates_df['gene_id'] = dates_df['gene_id'].astype(str)
    
    df = pd.merge(adj_df, dates_df, left_on='Future_Gene', right_on='gene_id', how='inner')
    
    df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
    df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')
    df = df[df['first_annotation_date'] > df['Date']].copy()
    
    if not df.empty:
        df['Delta_T_days'] = (df['first_annotation_date'] - df['Date']).dt.days
        df['GO_Term'] = term_formatted
        
        all_data.append(df[['GO_Term', 'PID0_mean_adj', 'Delta_T_days']])

master_df = pd.concat(all_data, ignore_index=True)

# DEFINE THRESHOLDS
thresholds = [
    ('1 Week', 7), ('1 Month', 30), ('3 Months', 90), ('6 Months', 180),
    ('1 Year', 365), ('3 Years', 1095), ('5 Years', 1825), ('10 Years', 3650)
]

# GENERATE MULTI-PAGE PDF
print(f"Drawing Master PDF to {output_pdf}...")

with PdfPages(output_pdf) as pdf:
    for label, days in thresholds:
        print(f" -> Rendering page for {label} threshold...")
        
        # Create a copy for this specific page's logic
        page_df = master_df.copy()
        
        # Categorize predictions based on the current threshold
        page_df['Condition'] = np.where(
            page_df['Delta_T_days'] <= days, 
            f'≤ {label}', 
            f'> {label}'
        )
        
        # Calculate Fold Changes and Mann-Whitney U Test
        fc_labels = []
        go_terms = sorted(page_df['GO_Term'].unique())
        
        for term in go_terms:
            subset = page_df[page_df['GO_Term'] == term]
            short_data = subset[subset['Condition'] == f'≤ {label}']['PID0_mean_adj']
            long_data = subset[subset['Condition'] == f'> {label}']['PID0_mean_adj']
            
            mean_short = short_data.mean()
            mean_long = long_data.mean()
            
            # Failsafe for missing data or zero division
            if pd.isna(mean_short) or pd.isna(mean_long) or mean_long == 0:
                fc_labels.append(f"{term}\nFC: N/A\np = N/A")
                continue
                
            fc = mean_short / mean_long
            
            # Wilcoxon (Mann-Whitney U) Test
            if len(short_data) >= 3 and len(long_data) >= 3:
                stat, p_val = mannwhitneyu(short_data, long_data, alternative='two-sided')
                
                # Format the p-value
                if p_val < 0.001:
                    p_str = "p < 0.001 ***"
                elif p_val < 0.01:
                    p_str = f"p = {p_val:.3f} **"
                elif p_val < 0.05:
                    p_str = f"p = {p_val:.3f} *"
                else:
                    p_str = f"p = {p_val:.2f} (ns)"
            else:
                p_str = "p = N/A"
                
            fc_labels.append(f"{term}\nFC: {fc:.2f}x\n{p_str}")
        
        # -- DRAW PLOT --
        plt.figure(figsize=(16, 9))
        sns.set_theme(style="whitegrid")
        
        palette = {f'≤ {label}': '#ea7373', f'> {label}': '#7cb3e8'}
        
        ax = sns.violinplot(
            data=page_df, 
            x='GO_Term', 
            y='PID0_mean_adj', 
            hue='Condition', 
            split=True, 
            inner="quartile",
            palette=palette,
            linewidth=1.2,
            cut=0, 
            bw_adjust=0.5, 
            order=go_terms
        )
        
        # Formatting
        plt.title(f"GO {aspect.upper()} Adjacency Distribution: ≤ {label} vs > {label} Prior to Annotation", fontsize=22, fontweight='bold', pad=20)
        plt.ylabel("Adjacency (Mean Probability From Annotated)", fontsize=18)
        plt.xlabel("GO Term", fontsize=18)
        
        ax.set_xticklabels(fc_labels, fontsize=14)
        plt.yticks(fontsize=14)
        
        # Dynamic Y-Limit
        y_max = page_df['PID0_mean_adj'].quantile(0.98)
        
        if not pd.isna(y_max) and y_max > 0:
            plt.ylim(-(y_max * 0.05), y_max + (y_max * 0.05))
            
        plt.legend(title='Time to Annotation', fontsize=14, title_fontsize=16, loc='upper right')
        
        plt.grid(True, axis='y', linestyle='--', alpha=0.6)
        plt.tight_layout()
        
        # Save page to PDF and clear canvas
        pdf.savefig()
        plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [GLOBAL {aspect.upper()}] Master PDF successfully generated in {elapsed_time} seconds! ---")