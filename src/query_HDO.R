library(HDO.db)
library(tidyverse)
library(org.Hs.eg.db)
library(arrow)

### args

input_pod <- snakemake@input$bp_frequencies
output_pod <- snakemake@output$annotation_df

### Input
pod_df <- read_parquet(
    input_pod
)

### Processing
all_gene <- unique(pod_df$uniprot_id)
entrez_id <- mapIds(
    org.Hs.eg.db,
    keys = all_gene,
    column = "ENTREZID", keytype = "UNIPROT", multiVals = "first"
)

gene_df <- data.frame(
    uniprot_id = names(entrez_id),
    entrez_id = unlist(entrez_id)
)

doids <- select(
    x = HDO.db,
    keys = gene_df$entrez_id, keytype = "gene",
    columns = c("doid")
)

colnames(doids) <- c("doid", "entrez_id")

full <- merge(doids, gene_df, by = "entrez_id", all = TRUE)
full <- full[complete.cases(full), ]

full %>%
    group_by(uniprot_id) %>%
    summarise(count = n_distinct(doid)) %>%
    arrange(desc(count)) -> doid_counts

write.table(doid_counts, output_pod, sep = "\t", row.names = FALSE)

df <- read.table("work_folder/data/HDO/anotation_per_uniprot.csv", sep="\t", header = TRUE)

ggplot(df, aes(x = count_degree, count_annot)) +
    geom_point() +
    theme_bw()
