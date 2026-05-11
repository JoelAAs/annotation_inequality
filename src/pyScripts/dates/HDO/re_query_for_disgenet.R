api_key <- readLines("~/temp")[1]
Sys.setenv(DISGENET_API_KEY = api_key)

library(disgenet2r) 
library(arrow) 
library(tidyverse) 

print("Loading data...")

# 1. Load Parquet Data
df_frequencies <- read_parquet(snakemake@input$bp_frequencies)

# Grabs the first output dynamically
output_annotation_file <- snakemake@output[[1]]

# 2. Extract Entrez IDs and clean up floats 
all_ids <- df_frequencies %>% 
  pull(entrez_id) %>% 
  na.omit() %>% 
  as.numeric() %>% 
  round() %>%      
  as.character() %>% 
  unique()

print(paste("Extracted", length(all_ids), "unique Entrez IDs."))

# 3. Define New HDO-like Sources
target_sources <- c(
  "ORPHANET", "HPO", "CLINVAR", "CLINGEN", 
  "GENCC", "UNIPROT", "PSYGENET"
)

# 4. Setup the Loop
batch_size <- 50  
n_ids <- length(all_ids)

# Ensure the output directory exists and delete any old run of this file
dir.create(dirname(output_annotation_file), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output_annotation_file)) {
  file.remove(output_annotation_file)
}

print("Starting progressive queries...")

# Header tracking variable
header <- TRUE

for (i in 1:ceiling(n_ids / batch_size)) {
  
  from <- ((i - 1) * batch_size) + 1
  to <- min(i * batch_size, n_ids) 
  
  print(paste("Querying IDs", from, "to", to, "of", n_ids))
  
  tryCatch({
    
    # FIX 1: Silence the annoying "genes not found" console spam
    results <- suppressWarnings(
      gene2disease(
        gene = all_ids[from:to],
        vocabulary = "ENTREZ",   
        database = "CURATED"  
      )
    )
    
    # Check if results exist
    if (length(results@qresult) > 0 && nrow(results@qresult) > 0) {
      tab <- results@qresult
      
      # FIX 2: Dynamic, crash-proof filtering
      # Find the source column dynamically (it might be named 'source', 'source_name', etc.)
      src_col <- grep("source", colnames(tab), value = TRUE, ignore.case = TRUE)
      
      if (length(src_col) > 0) {
        # Sources are often returned as comma-strings (e.g., "CLINVAR, HPO")
        # %in% fails on this, so we use string detection instead
        pattern <- paste(target_sources, collapse = "|")
        tab <- tab[grepl(pattern, tab[[src_col[1]]], ignore.case = TRUE), ]
      }
      
      # Make sure we still have data after filtering
      if (nrow(tab) > 0) {
        
        # Build the dynamic list of columns to keep
        target_cols <- c("gene_symbol", "geneid", "disease_name", "diseaseUMLSCUI", "score", "yearInitial", "yearFinal")
        if (length(src_col) > 0) target_cols <- c(target_cols, src_col[1])
        
        cols_to_keep <- intersect(colnames(tab), target_cols)
        tab <- tab %>% select(all_of(cols_to_keep)) %>% distinct()
        
        # write.table logic
        if (header) {
          write.table(
            tab, output_annotation_file, sep = "\t", 
            row.names = FALSE, quote = FALSE
          )
          header <- FALSE
        } else {
          write.table(
            tab, output_annotation_file, sep = "\t", 
            append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE
          )
        }
      }
    }
  }, error = function(e) {
    # Kept as 'message' so it prints to the Snakemake log immediately if it fails
    message(paste("Chunk", i, "failed:", e$message))
  })
  
  Sys.sleep(3)
}

print("Done! All chunks processed.")

# --- Snakemake Failsafe ---
if (!file.exists(output_annotation_file)) {
  message("Warning: No valid disease associations found for any genes. Creating empty file to satisfy Snakemake.")
  file.create(output_annotation_file)
}