library(HDO.db)
library(DO.db)
library(AnnotationDbi)
library(tidyverse)
library(igraph)
library(arrow)

input_pod <- snakemake@input$bp_frequencies
output_pod_df <- snakemake@output$annot_df

cat("Loading data from Parquet...\n")
df <- read_parquet(input_pod)

# Extract unique genes as characters, removing any NAs
all_genes <- df %>% 
    filter(!is.na(entrez_id)) %>%
    distinct(entrez_id) %>% 
    mutate(entrez_id_char = as.character(entrez_id))

cat(sprintf("Querying offline HDO.db for %d genes...\n", nrow(all_genes)))

doids <- tryCatch({
    # EXPLICITLY call AnnotationDbi::select to avoid tidyverse masking
    AnnotationDbi::select(HDO.db, keys = all_genes$entrez_id_char, keytype = "gene", columns = c("doid"))
}, error = function(e) {
    cat("\n--- BIOCONDUCTOR ERROR ---\n")
    cat("Actual Error: ", e$message, "\n")
    stop("Query failed! See the actual error message above.")
})

# Safely rename the 'gene' column to match our data and filter out empty query results
doids <- doids %>% 
    rename(entrez_id_char = gene) %>% 
    filter(!is.na(doid))

cat("Building DO graph from DO.db to find ancestors...\n")
parent_map <- as.list(DOPARENTS)
edges <- do.call(rbind, lapply(names(parent_map), function(child) {
  parents <- parent_map[[child]]
  if (is.null(parents) || all(is.na(parents))) return(NULL)
  data.frame(from = parents, to = child)
}))
g <- graph_from_data_frame(edges, directed = TRUE)

# Calculate global depths from root
root <- "DOID:4"
dists <- distances(g, v = root, mode = "out")
depth_df <- data.frame(
  doid = colnames(dists),
  depth = as.numeric(dists[1, ]),
  stringsAsFactors = FALSE
)

cat("Expanding hierarchy to include all ancestors...\n")
unique_real_doids <- unique(doids$doid)

# subcomponent(mode='in') gets all ancestors rapidly
ancestor_list <- lapply(unique_real_doids, function(d) {
    if(!d %in% V(g)$name) return(d)
    names(subcomponent(g, d, mode = "in"))
})
names(ancestor_list) <- unique_real_doids

# Expand the dataframe to include all ancestor DOIDs for each gene
expanded_df <- doids %>%
    group_by(entrez_id_char) %>%
    do({
        all_anc <- unique(unlist(ancestor_list[.$doid]))
        data.frame(doid = all_anc, stringsAsFactors = FALSE)
    }) %>%
    ungroup()

cat("Adding depths and DO terms...\n")
terms <- as.list(DOTERM)
term_df <- data.frame(
    doid = names(terms),
    term = sapply(terms, function(x) x@Term),
    stringsAsFactors = FALSE
)

# Process the successfully mapped genes that connect back to DOID:4
mapped_df <- expanded_df %>%
    left_join(depth_df, by = "doid") %>%
    left_join(term_df, by = "doid") %>%
    filter(depth >= 0) %>% 
    dplyr::select(entrez_id_char, doid, depth, term) %>%
    distinct()

# Find any gene from our starting list that isn't in mapped_df
# This catches both genes that had NO mapping, and genes whose only mappings were disconnected/obsolete.
missing_ids <- setdiff(all_genes$entrez_id_char, mapped_df$entrez_id_char)

unmapped_df <- data.frame(
    entrez_id_char = missing_ids,
    doid = "No_doid",
    depth = -1,
    term = "No_term",
    stringsAsFactors = FALSE
)

cat(sprintf("Found %d fully mapped genes and %d unmapped/disconnected genes.\n", 
            length(unique(mapped_df$entrez_id_char)), 
            nrow(unmapped_df)))

# Combine and format
final_df <- bind_rows(mapped_df, unmapped_df) %>%
    rename(entrez_id = entrez_id_char) %>%
    mutate(entrez_id = as.integer(entrez_id)) %>%
    arrange(entrez_id, depth)

no_doid_count <- sum(final_df$doid == "No_doid")
cat(sprintf("Validation: Output contains %d rows with 'No_doid'.\n", no_doid_count))

cat(sprintf("FINAL RESULT: Saving %d total annotations to output.\n", nrow(final_df)))
write.table(final_df, output_pod_df, sep = "\t", row.names = FALSE, quote = FALSE)
cat("HDO to gene annotations ready!\n")