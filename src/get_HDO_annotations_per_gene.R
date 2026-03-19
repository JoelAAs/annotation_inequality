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

### Processing
doids <- select(
    x = HDO.db,
    keys = pod_df$entrez_id, keytype = "gene",
    columns = c("doid", "term")
)

colnames(doids) <- c("doid", "entrez_id", "annotation")
pod_df <- pod_df %>%
    distinct(entrez_id, .keep_all = TRUE)

full <- merge(doids, pod_df, by = "entrez_id", all = TRUE)
full <- full[, c("entrez_id", "annotation")]
full <- na.omit(full)

annotations_list <- full %>%
    distinct(annotation, .keep_all = TRUE)
annotations_list <- annotations_list[, "annotation", drop = FALSE]

full <- full %>%
    distinct(entrez_id, annotation) %>%
    group_by(entrez_id) %>%
    summarise(annotations = list(annotation))
full$annotations <- sapply(full$annotations, function(x) paste(x, collapse = ";"))

print("HDO annotations obtained!")

write.table(full, output_pod_df, sep = "\t", row.names = FALSE)
write.table(annotations_list, output_pod_list, sep = '\t', row.names = FALSE)