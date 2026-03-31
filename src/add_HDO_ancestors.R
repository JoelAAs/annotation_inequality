library(DO.db)
library(AnnotationDbi)
library(dplyr)
library(igraph)

annot_df = snakemake@input$annotation_df
output_df = snakemake@output$complete_annotations

cat("Processing ancestors from DOIDS...\n")

## Load the dataframe
cat("Loading the dataframe...\n")
df <- read.csv(
    annot_df,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
)
cat("Dataframe loaded!\n")

## Fill the Na values
cat("Filling the Na values...\n")
df <- df %>% select(-annotation)
df$doid[is.na(df$doid)] <- "No_doid"
cat("Na values filled!\n")

# Access DO.db mappings directly
cat("Loading DO.db mappings...\n")
parent_map <- as.list(get("DOPARENTS", envir = asNamespace("DO.db")))
cat("Mappings loaded!\n")

## Compute GLOBAL depth from root
cat("Computing global depth from root...\n")
edges <- do.call(rbind, lapply(names(parent_map), function(child) {
  parents <- parent_map[[child]]
  if (is.null(parents) || all(is.na(parents))) return(NULL)
  data.frame(from = parents, to = child)
}))
g <- graph_from_data_frame(edges, directed = TRUE)
root <- "DOID:4"
distances <- distances(g, v = root, mode = "out")
depth_df <- data.frame(
  doid = colnames(distances),
  depth = as.numeric(distances[1, ])
)
cat("Depth computed!\n")

## Get ancestors recursively
cat("Loading hierarchy logic...\n")
getAllDOAncestors <- function(doid, parentMap = parent_map) {
  if(doid == "No_doid") {
    return(data.frame(ancestor = doid, stringsAsFactors = FALSE))
  }
  parents <- parentMap[[doid]]
  if(is.null(parents) || all(is.na(parents))) {
    return(data.frame(ancestor = doid, stringsAsFactors = FALSE))
  }
  paths <- lapply(parents, function(p) getAllDOAncestors(p, parentMap))
  paths_df <- do.call(rbind, paths)
  rbind(paths_df, data.frame(ancestor = doid, stringsAsFactors = FALSE))
}
cat("Logic ready!\n")

## Cache unique DOIDs
cat("Caching unique DOIDs...\n")
unique_doids <- unique(df$doid)
ancestor_cache <- lapply(unique_doids, function(d) getAllDOAncestors(d))
names(ancestor_cache) <- unique_doids
cat("Ancestor cache built!\n")

# Expand dataframe with hierarchy
cat("Expanding dataframe with hierarchy...\n")
expanded_list <- lapply(seq_len(nrow(df)), function(i) {
  row <- df[i, ]
  ancestors_df <- ancestor_cache[[row$doid]]
  ancestors_df$entrez_id <- row$entrez_id
  ancestors_df
})
full_df <- do.call(rbind, expanded_list) %>%
  distinct(ancestor, entrez_id, .keep_all = TRUE) %>%
  rename(doid = ancestor)

# Add the root DOID:4 ONLY for entrez_ids that have at least one real DOID
cat("Adding root DOID:4 for valid entrez_ids...\n")

valid_entrez_ids <- df %>%
  filter(doid != "No_doid") %>%
  pull(entrez_id) %>%
  unique()

root_df <- data.frame(
  entrez_id = valid_entrez_ids,
  doid = "DOID:4",
  stringsAsFactors = FALSE
)

full_df <- bind_rows(full_df, root_df) %>%
  distinct(doid, entrez_id, .keep_all = TRUE)

# Assign depth
cat("Assigning depth...\n")
full_df <- full_df %>%
  left_join(depth_df, by = "doid") %>%
  arrange(entrez_id, depth)
cat("Depth assigned!\n")

# Write output
cat("Writing csv...\n")
write.table(full_df, output_df, sep = '\t', row.names = FALSE)
cat("Csv obtained!\n")
cat("Ancestors processed!\n")