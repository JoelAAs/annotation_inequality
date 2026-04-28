import gzip

snapshot_file = snakemake.input.snapshot
genes_file = snakemake.input.genes
output_file = snakemake.output[0]

target_gos = set(filter(None, snakemake.params.target_gos.split(",")))

with open(genes_file, "r") as f:
    target_genes = set(line.strip() for line in f if line.strip())

with gzip.open(snapshot_file, "rt", encoding="utf-8") as f_in, open(output_file, "w", encoding="utf-8") as f_out:
    f_out.write("Gene\tGO_ID\tDate\n")
    
    for line in f_in:
        if line.startswith("!"): continue
            
        parts = line.strip("\n").split("\t")
        if len(parts) < 14: continue
            
        gene = parts[2]
        go_id = parts[4]
        date = parts[13] 
        
        if gene in target_genes:
            if not target_gos or go_id in target_gos:
                f_out.write(f"{gene}\t{go_id}\t{date}\n")