api_key <- readLines("~/temp")[1]
Sys.setenv(DISGENET_API_KEY = api_key)

library(disgenet2r) # for disgenet
library(arrow) # for reading parquet
library(tidyverse) # general QoL stuff


df_frequencies <- read_parquet(
  snakemake@input$bp_frequencies
)

df_frequencies <- read_parquet(
  "work_folder/data/intact/bait_prey_frequecies.pq"
)

output_annotation_file <- snakemake@output$annotation_df

all_ids <- unique(df_frequencies$uniprot_id)
header <- TRUE
i <- 0
batch_size <- 99  # 1 index so 100
n_ids <- length(all_ids)
for (i in 1:ceiling(n_ids / batch_size)) {
  from <- (i - 1) * batch_size
  to <- i * batch_size
  if (to > n_ids) {
    to <- n_ids
  }

  print(paste(to, "of", n_ids))
  results <- gene2disease(
    gene = all_ids[from:to],
    vocabulary = "UNIPROT",
    database = "CURATED"
  )

  Sys.sleep(5)

  if (class(results) != "character") {
    tab <- unique(
      results@qresult[
        ,
        c(
          "gene_symbol", "uniprotids",
          "disease_name", "diseaseUMLSCUI", "score",
          "yearInitial", "yearFinal"
        )
      ]
    )
    tab <- tab[order(tab$uniprotids), ]
    if (header) {
      write.table(
        tab,
        output_annotation_file,
        sep = "\t",
        row.names = FALSE
      )
      header <- FALSE
    } else {
      write.table(
        tab,
        output_annotation_file,
        sep = "\t",
        append = TRUE,
        col.names = FALSE,
        row.names = FALSE
      )
    }
  }
}
