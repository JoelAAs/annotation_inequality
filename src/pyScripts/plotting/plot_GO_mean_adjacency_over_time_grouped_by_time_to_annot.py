import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import time

input_adj_file = snakemake.input.mean_adj_file
input_dates_file = snakemake.input.annot_dates_file
output_file = snakemake.output.plot_file
term = snakemake.wildcards.term

term_formatted = term.replace('_', ':')

start_time = time.time()
print(f"--- [{term}] Starting Time-to-Annotation Relative Quartiles (Horizon Cutoff) ---")

# LOAD ADJACENCY DATA
print(f"[{term}] [LOAD] Reading adjacency data...")
adj_df = pd.read_parquet(input_adj_file)
adj_df.columns = adj_df.columns.str.strip()

# Filter for only the true annotations immediately
adj_df = adj_df[['Date', 'Future_Gene', 'PID0_mean_adj']].copy()

# LOAD ANNOTATION DATES DATA
print(f"[{term}] [LOAD] Reading first annotation dates...")
dates_df = pd.read_csv(input_dates_file, sep='\t')
dates_df = dates_df[dates_df['GO_id'] == term_formatted].copy()

# PREPARE FOR MERGE
print(f"[{term}] [PROCESS] Merging datasets and calculating Delta T...")
adj_df['Future_Gene'] = adj_df['Future_Gene'].astype(str)
dates_df['gene_id'] = dates_df['gene_id'].astype(str)

# Inner Join
df = pd.merge(
    adj_df, 
    dates_df[['gene_id', 'first_annotation_date']], 
    left_on='Future_Gene', 
    right_on='gene_id', 
    how='inner'
)

# Parse the dates
df['Date'] = pd.to_datetime(df['Date'].astype(str), format='%Y%m%d')
df['first_annotation_date'] = pd.to_datetime(df['first_annotation_date'].astype(str), format='%Y%m%d')

# --- VECTORIZED MATH & HORIZON CUTOFF ---
# 1. Strictly filter for genes annotated AFTER the network snapshot date
df = df[df['first_annotation_date'] > df['Date']].copy()

# 2. Calculate time-to-annotation in days
df['Delta_T'] = (df['first_annotation_date'] - df['Date']).dt.days

# 3. THE HORIZON CUTOFF: Prevent the end-of-timeline squeeze
max_annot_date = df['first_annotation_date'].max()
horizon_cutoff = max_annot_date - pd.DateOffset(years=3)

print(f"[{term}] [INFO] Global max annotation date: {max_annot_date.date()}")
print(f"[{term}] [INFO] Imposing Horizon Cutoff at: {horizon_cutoff.date()} to preserve stratification.")

df = df[df['Date'] <= horizon_cutoff].copy()


# QUARTILE MATH (Relative Percentages)
def compute_quartile_stats(group):
    # Failsafe: need a minimal number of genes to make percentiles meaningful
    if len(group) < 12: 
        return pd.Series({
            'Close_mean': np.nan, 'Close_sem': np.nan,
            'Middle_mean': np.nan, 'Middle_sem': np.nan,
            'Far_mean': np.nan, 'Far_sem': np.nan
        })
    
    group = group.copy()
    
    # Rank-based splitting to guarantee perfectly even bin sizes and prevent tie-crashing
    group['Rank'] = group['Delta_T'].rank(method='first')
    
    # Explicit 25% / 50% / 25% splits to group Q2 and Q3 together cleanly
    labels = ['Close', 'Middle', 'Far']
    group['Quartile'] = pd.qcut(group['Rank'], q=[0, 0.25, 0.75, 1.0], labels=labels)
    
    # Compute mean and Standard Error of the Mean (SEM)
    stats = {}
    for q_name in ['Close', 'Middle', 'Far']:
        subset = group[group['Quartile'] == q_name]['PID0_mean_adj']
        stats[f'{q_name}_mean'] = subset.mean()
        stats[f'{q_name}_sem'] = subset.sem() if len(subset) > 1 else 0.0
        
    return pd.Series(stats)

# Apply the math to every snapshot date
print(f"[{term}] [MATH] Slicing timelines into dynamically generated percentiles...")
time_to_annot_df = df.groupby('Date').apply(compute_quartile_stats, include_groups=False).reset_index()
time_to_annot_df = time_to_annot_df.dropna(subset=['Close_mean']).sort_values('Date')


# PLOT
print(f"[{term}] [PLOT] Drawing Time-to-Annotation trajectories...")
plt.figure(figsize=(14, 8))

# Define colors for the trajectories
categories = [
    ('Close', 'Fastest 25% of Annotations (Close)', '#d73027', '#f46d43'),
    ('Middle', 'Middle 50% of Annotations', '#fdae61', '#fee08b'),
    ('Far', 'Slowest 25% of Annotations (Far)', '#4575b4', '#91bfdb')
]

for cat_key, label, line_color, band_color in categories:
    mean_col = f'{cat_key}_mean'
    sem_col = f'{cat_key}_sem'
    
    # 95% Confidence Bounds (Mean ± 1.96 * SEM)
    upper_bound = time_to_annot_df[mean_col] + (1.96 * time_to_annot_df[sem_col])
    lower_bound = time_to_annot_df[mean_col] - (1.96 * time_to_annot_df[sem_col])
    
    # Draw the Error Band (zorder=1 to push behind lines)
    plt.fill_between(
        time_to_annot_df['Date'],
        lower_bound,
        upper_bound,
        color=band_color,
        alpha=0.35,
        zorder=1
    )
    
    # Draw the Mean Line (Markers removed for smooth, continuous lines)
    plt.plot(
        time_to_annot_df['Date'],
        time_to_annot_df[mean_col],
        label=label,
        color=line_color,
        linewidth=3.0,
        zorder=2
    )

# FORMATTING
plt.title(f"Mean Adjacency Scores Over Time Grouped by Relative $\\Delta$T\n{term_formatted}", fontsize=26, fontweight='bold', pad=20)
plt.xlabel("Network Snapshot Date", fontsize=20, labelpad=15)
plt.ylabel("Mean Adjacency Score", fontsize=20, labelpad=15)

ax = plt.gca()
ax.xaxis.set_major_locator(mdates.YearLocator(2)) 
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.xticks(fontsize=16, rotation=45)
plt.yticks(fontsize=16)

plt.legend(
    loc="upper right", # Moved to top right to avoid colliding with early-year data points
    fontsize=16, 
    framealpha=0.9, 
    edgecolor='#cccccc',
    title="Relative Time-to-Annotation",
    title_fontsize=18
)

plt.grid(True, linestyle='--', alpha=0.5, zorder=0)
plt.tight_layout()

# SAVE
print(f"[{term}] [SAVE] Saving image to disk...")
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

elapsed_time = round(time.time() - start_time, 2)
print(f"--- [{term}] Horizon-Cutoff plot successfully generated in {elapsed_time} seconds! ---")