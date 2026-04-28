import pandas as pd

entrez_file = snakemake.input.network_genes
mapping_file = snakemake.input.mapping
output_file = snakemake.output[0]

# 1. Load your network's Entrez IDs
with open(entrez_file, "r") as f:
    entrez_ids = set(line.strip() for line in f if line.strip())

# 2. Load the official HGNC table (forcing strings so Entrez IDs don't become floats)
df_hgnc = pd.read_csv(mapping_file, sep="\t", dtype=str, low_memory=False)

mapped_names = set()

# 3. Search the table for your IDs
entrez_col = "entrez_id" if "entrez_id" in df_hgnc.columns else "ncbi_id"

# Clean the column (remove decimals like 5050.0)
df_hgnc[entrez_col] = df_hgnc[entrez_col].str.split('.').str[0]

# Filter the table for rows that match your network IDs
mask = df_hgnc[entrez_col].isin(entrez_ids)
matched_df = df_hgnc[mask]

# Extract Symbols (TP53, etc.)
if "symbol" in matched_df.columns:
    mapped_names.update(matched_df["symbol"].dropna().str.upper())

# Extract UniProt IDs (P04637, etc.) - handling pipe-separated strings
if "uniprot_ids" in matched_df.columns:
    u_series = matched_df["uniprot_ids"].dropna().str.replace('"', '').str.split('|').explode()
    mapped_names.update(u_series.str.strip().str.upper())

# 4. Write the translated names to a new text file for the GAF parser
with open(output_file, "w") as f:
    for name in sorted(mapped_names):
        f.write(f"{name}\n")

print(f"[STATUS] Translated {len(entrez_ids)} Entrez IDs into {len(mapped_names)} valid Symbols/UniProt targets.")