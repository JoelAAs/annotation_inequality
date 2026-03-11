library(org.Hs.eg.db)
library(tidyverse)
library(arrow)

bp_df <- read_parquet(snakemake@input$formated_in)

all_uniprot_id <- c(unique(bp_df$uniprot_id_bait), unique(bp_df$uniprot_id_prey))
all_uniprot_id <- unique(unlist(sapply(all_uniprot_id, function(x) str_split(x, "-")[[1]][1])))
entrez_id <- mapIds(
    org.Hs.eg.db,
    keys = all_uniprot_id,
    column = "ENTREZID", keytype = "UNIPROT", multiVals = "first"
)

gene_df <- data.frame(
    uniprot_id = names(entrez_id),
    entrez_id = unlist(entrez_id)
)

bp_df <- merge(bp_df, gene_df, by.x = "uniprot_id_bait", by.y="uniprot_id", all = TRUE)
bp_df <- merge(bp_df, gene_df, by.x = "uniprot_id_prey", by.y="uniprot_id", all = TRUE, suffixes=c("_bait", "_prey"))

bp_df <- bp_df[complete.cases(bp_df), ]

write_parquet(bp_df, snakemake@output$formated_out)
