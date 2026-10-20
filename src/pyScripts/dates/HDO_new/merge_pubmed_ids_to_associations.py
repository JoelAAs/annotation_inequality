import pandas as pd

file_associations = snakemake.input.associations
file_dates = snakemake.input.pmid_dates
file_output_full = snakemake.output.full_dataset
file_output_temporal = snakemake.output.only_associated

# Load datasets
print("--- [INFO] Loading datasets... ---")
df_assoc = pd.read_parquet(file_associations)
df_dates = pd.read_parquet(file_dates)

print(f"    > Original associations: {len(df_assoc)} rows")
print(f"    > Publication dates: {len(df_dates)} unique PMIDs")

# Merge and format
print("--- [INFO] Merging temporal data with the associations... ---")
# Left merge: keeps all rows from df_assoc, adds Publication_Date where Reference matches
df_full = df_assoc.merge(df_dates, on='Reference', how='left')

# Convert date format from YYYY-MM-DD to YYYYMMDD by removing hyphens
df_full['Publication_Date'] = df_full['Publication_Date'].str.replace('-', '', regex=False)

# -Extract 
print("--- [INFO] Extracting the strict time-resolved associations... ---")
# Drop rows where Publication_Date is missing (NaN)
df_temporal = df_full.dropna(subset=['Publication_Date']).copy()

# Sort chronologically to see the biological timeline unfold (oldest to newest)
df_temporal = df_temporal.sort_values(by='Publication_Date')

# Final stats and Save
print("\n--- [INFO] Final Associations Statistics ---")
print(f"    > Full Associations (Global): {len(df_full)}")
print(f"    > Time-Resolved Associations (Dated associations only): {len(df_temporal)}")

df_full.to_parquet(file_output_full, index=False)
df_temporal.to_parquet(file_output_temporal, index=False)

print("\n--- [INFO] Process completed successfully. Data saved to Parquet. ---")