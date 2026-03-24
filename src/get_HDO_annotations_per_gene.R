library(HDO.db)
library(tidyverse)
library(org.Hs.eg.db)
library(arrow)
library(dplyr)

### args

input_pod <- snakemake@input$bp_frequencies
output_pod_df <- snakemake@output$annotation_df
output_pod_list <- snakemake@output$annotation_list

### Input
pod_df <- read_parquet(
    input_pod
)

pod_df <- pod_df %>%
    distinct(entrez_id)
cat("Unique entrez_id values in bp_frequencies = ", length(unique(pod_df$entrez_id)))
all_genes <- data.frame(entrez_id = pod_df$entrez_id)
all_genes <- all_genes %>% mutate(entrez_id = as.character(entrez_id))

### Processing
doids <- select(
    x = HDO.db,
    keys = pod_df$entrez_id, keytype = "gene",
    columns = c("doid", "term")
)

cat("Unique entrez_id values after HDO query = ", length(unique(pod_df$entrez_id)))

colnames(doids) <- c("doid", "entrez_id", "annotation")

doids <- doids[!duplicated(doids), ]

annotations_list <- doids %>%
    distinct(annotation, .keep_all = TRUE)
annotations_list <- annotations_list[, "annotation", drop = FALSE]

full <- doids %>%
    full_join(all_genes, by = "entrez_id")

cat("Unique entrez_id values after adding ids with no annotations = ", length(unique(pod_df$entrez_id)))

print("HDO annotations obtained!")

write.table(full, output_pod_df, sep = "\t", row.names = FALSE)
write.table(annotations_list, output_pod_list, sep = '\t', row.names = FALSE)