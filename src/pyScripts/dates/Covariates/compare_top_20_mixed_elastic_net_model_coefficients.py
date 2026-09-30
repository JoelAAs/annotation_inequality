import os
import pandas as pd
import pronto

coef_path = snakemake.input.en_coefficients
go_obo_path = snakemake.input.go_obo
hdo_obo_path = snakemake.input.hdo_obo

out_dir = snakemake.output.comparisons

print("--- [START] TOP 20 COEFFICIENTS COMPARISON STARTED ---", flush=True)

# Create the output directory
os.makedirs(out_dir, exist_ok=True)

# LOAD ONTOLOGIES
print("--- [RUNNING] LOADING ONTOLOGIES ---", flush=True)
go_ontology = pronto.Ontology(go_obo_path)
hdo_ontology = pronto.Ontology(hdo_obo_path)

def get_term_name(term_id):
    if term_id in go_ontology:
        return go_ontology[term_id].name
    elif term_id in hdo_ontology:
        return hdo_ontology[term_id].name
    return term_id

# PREPARE DATASET
print("--- [RUNNING] PROCESSING COEFFICIENTS ---", flush=True)
df = pd.read_parquet(coef_path)

# Remove intercept and unpenalized covariate
df = df[~df['annotation_id'].isin(['(Intercept)', 'n_annot'])].copy()

# Ensure selection_freq has no NAs and compute absolute coefficient
df['selection_freq'] = df['selection_freq'].fillna(0.0)
df['abs_coefficient'] = df['coefficient'].abs()

# Define thresholds and get sorted alphas
freq_thresholds = [0.50, 0.70, 0.80, 0.90, 0.95]
alphas = sorted(df['alpha'].unique())

# BUILD COMPARISON DATAFRAMES
for thresh in freq_thresholds:
    print(f"--- [RUNNING] Building table for threshold: {thresh:.2f} ---", flush=True)
    
    # Dictionary to hold a column for each alpha
    alpha_columns = {}
    
    for alpha in alphas:
        # Filter for current alpha and threshold
        sub_df = df[(df['alpha'] == alpha) & (df['selection_freq'] > thresh)].copy()
        
        # Get top 20 absolute coefficients
        top20 = sub_df.sort_values(by='abs_coefficient', ascending=False).head(20)
        
        # Generate the formatted labels
        labels = []
        for _, row in top20.iterrows():
            annot_id = row['annotation_id']
            name = get_term_name(annot_id)
            coeff = row['coefficient']
            # Format: "Term Name (ID) [c: 0.123]"
            label = f"{name} ({annot_id}) [c: {coeff:.3f}]"
            labels.append(label)
            
        alpha_columns[f"Alpha_{alpha:.2f}"] = pd.Series(labels)
        
    # Combine into a single DataFrame
    comparison_df = pd.DataFrame(alpha_columns)
    
    # Clean up the index to show ranks (Rank_1, Rank_2, ..., Rank_20)
    comparison_df.index = [f"Rank_{i+1}" for i in range(len(comparison_df))]
    
    # Fill empty slots (NaNs) with a dash for cleaner TSV reading
    comparison_df.fillna("-", inplace=True)
    
    # Save the dataframe as a TSV file inside the directory
    out_filename = f"freq_{thresh:.2f}.tsv"
    out_filepath = os.path.join(out_dir, out_filename)
    
    comparison_df.to_csv(out_filepath, sep='\t', index=True)
    
print(f"--- [COMPLETE] ALL COMPARISON TABLES SAVED IN {out_dir} ---", flush=True)