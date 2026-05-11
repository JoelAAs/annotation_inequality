import pandas as pd
import networkx as nx

ot_evidences_dir = snakemake.input.evidence_folder
ot_diseases_dir = snakemake.input.disease_folder
outputdf = snakemake.output.first_annotation_dates

evidence_df = pd.read_parquet(
    ot_evidences_dir,
    columns = ['targetId', 'diseaseId', 'evidenceDate']
)
disease_df = pd.read_parquet(
    ot_diseases_dir,
    columns = ['id', 'dbXRefs']
)

print("Building the Master Translation Dictionary...")

# Drop diseases that don't have cross-references
df_refs = disease_df.dropna(subset=['dbXRefs'])

# Explode the lists so every cross-reference gets its own row
exploded_df = df_refs.explode('dbXRefs')

# Create a dictionary that maps ANY alternate ID to its Open Targets ID
# Format: { "MONDO:0016241": "EFO_0000305", "DOID:1612": "EFO_0000305", ... }
xref_to_efo = dict(zip(exploded_df['dbXRefs'], exploded_df['id']))
print(f"Built a dictionary with {len(xref_to_efo)} translations!")


print("Translating the Evidence DataFrame...")

# Create a new column by mapping the old IDs through our dictionary
# If it finds 'MONDO:123', it translates to 'EFO_456'. 
# If it doesn't find it (like if it's ALREADY an EFO ID), it returns NaN.
evidence_df['diseaseId_EFO'] = evidence_df['diseaseId'].map(xref_to_efo)

# Fill the NaNs with the original ID
evidence_df['diseaseId_EFO'] = evidence_df['diseaseId_EFO'].fillna(evidence_df['diseaseId'])

efo_count = evidence_df['diseaseId_EFO'].astype(str).str.startswith('EFO').sum()
print(f"Translation complete! {efo_count} out of {len(evidence_df)} rows are now EFO formatted.")

print("Extracting the oldest dates for each Gene-Disease pair...")

# Drop rows that don't have a date 
clean_df = evidence_df.dropna(subset=['evidenceDate'])

# Group by the Gene and EFO column, and grab the minimum date
oldest_dates_df = clean_df.groupby(['targetId', 'diseaseId_EFO'])['evidenceDate'].min().reset_index()

oldest_dates_df = oldest_dates_df.rename(columns={'evidenceDate': 'first_publication_date'})

print(f"Success! Compressed down to {len(oldest_dates_df)} unique, oldest associations.")

oldest_dates_df.to_csv(outputdf, sep = '\t', index = False)