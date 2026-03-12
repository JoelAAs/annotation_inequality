library(HDO.db)
library(tidyverse)
library(org.Hs.eg.db)
library(arrow)
library(dplyr)

### args

input_pod <- snakemake@input$bp_frequencies
output_pod <- snakemake@output$annotation_df

### Input
pod_df <- read_parquet(
    input_pod
)


doids <- select(
    x = HDO.db,
    keys = pod_df$entrez_id, keytype = "gene",
    columns = c("doid")
)

colnames(doids) <- c("doid", "entrez_id")

full <- merge(doids, pod_df, by = "entrez_id", all = TRUE)
full %>%
    # Counting the number of different DO annotations for each entrez_id
    group_by(entrez_id) %>%
    summarise(count = n_distinct(doid)) %>%
    arrange(desc(count)) -> doid_counts

doid_counts <- doid_counts %>%
    filter(!is.na(entrez_id)) %>%
    mutate(count = replace_na(count, 0))
 
write.table(doid_counts, output_pod, sep = "\t", row.names = FALSE)