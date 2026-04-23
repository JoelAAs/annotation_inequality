import gzip

# Access Snakemake properties
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Parse targets, filtering out empty strings if the user provided empty lists
target_genes = set(filter(None, snakemake.params.genes.split(",")))
target_go = set(filter(None, snakemake.params.go_terms.split(",")))

snapshot_name = input_file.split("/")[-1]

with gzip.open(input_file, "rt", encoding="utf-8") as f_in, open(output_file, "w", encoding="utf-8") as f_out:
    for line in f_in:
        # Skip GAF headers/metadata
        if line.startswith("!"):
            continue
            
        parts = line.strip("\n").split("\t")
        
        # A valid GAF file must have at least 15 columns (or 14 in very early formats)
        if len(parts) < 14:
            continue
            
        # GAF Format Column Indices (0-based)
        gene_symbol = parts[2]
        go_id = parts[4]
        evidence = parts[6]
        date = parts[13]  # Format: YYYYMMDD
        
        # Check criteria (if target sets are empty, it acts as a passthrough for all)
        match_gene = (not target_genes) or (gene_symbol in target_genes)
        match_go = (not target_go) or (go_id in target_go)
        
        if match_gene and match_go:
            f_out.write(f"{gene_symbol}\t{go_id}\t{date}\t{evidence}\t{snapshot_name}\n")