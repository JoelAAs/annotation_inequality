# Load required libraries
suppressMessages(library(arrow))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages(library(Matrix))

# Snakemake variables initialization
feature_matrix_path <- snakemake@input[["feature_matrix"]]
bait_usage_path <- snakemake@input[["bait_usage"]]

coef_out <- snakemake@output[["en_coefficients"]]
metrics_out <- snakemake@output[["en_metrics"]]
model_out <- snakemake@output[["en_model"]]

# Threshold parameters from Snakemake
min_genes <- snakemake@params[["min_genes"]]

cat("--- [START] ELASTIC NET MODEL FITTING WITH GRID SEARCH (OLD) STARTED ---\n")

# 1. Load Data
cat("--- [RUNNING] LOADING DATA ---\n")
feature_matrix_df <- read_parquet(feature_matrix_path)
# FIX: Lettura del file separato da tabulazione
bait_count_df <- read.table(bait_usage_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract index (Entrez ID) safely
if ("entrez_id" %in% colnames(feature_matrix_df)) {
  universe_genes <- as.character(feature_matrix_df$entrez_id)
  feature_matrix_df$entrez_id <- NULL
} else if ("__index_level_0__" %in% colnames(feature_matrix_df)) {
  universe_genes <- as.character(feature_matrix_df[["__index_level_0__"]])
  feature_matrix_df[["__index_level_0__"]] <- NULL
} else {
  universe_genes <- as.character(feature_matrix_df[[1]])
  feature_matrix_df[[1]] <- NULL
}

# 2. Filtering Features (Memory efficient)
cat("--- [RUNNING] FILTERING FEATURE MATRIX ---\n")
col_sums <- colSums(feature_matrix_df)
valid_cols <- col_sums >= min_genes

feature_matrix_df <- feature_matrix_df[, valid_cols]
cat(sprintf("--- [INFO] FILTERING RETAINED %d OUT OF %d ANNOTATIONS (MIN >= %d) ---\n", sum(valid_cols), length(valid_cols), min_genes))

# 3. Convert to Sparse Matrix
cat("--- [RUNNING] CONVERTING TO SPARSE MATRIX ---\n")
X_sparse <- Matrix(as.matrix(feature_matrix_df), sparse = TRUE)
rownames(X_sparse) <- universe_genes

# Deep memory cleanup
rm(feature_matrix_df)
gc()

# 4. Aligning Target Variable (y)
cat("--- [RUNNING] ALIGNING TARGET VARIABLE (Y) ---\n")
y <- rep(0, length(universe_genes))
names(y) <- universe_genes

bait_id_col <- ifelse("entrez_id_bait" %in% colnames(bait_count_df), "entrez_id_bait", colnames(bait_count_df)[1])
count_col <- ifelse("count" %in% colnames(bait_count_df), "count", colnames(bait_count_df)[2])

bait_count_df[[bait_id_col]] <- as.character(bait_count_df[[bait_id_col]])
valid_baits <- bait_count_df[bait_count_df[[bait_id_col]] %in% universe_genes, ]

y[valid_baits[[bait_id_col]]] <- valid_baits[[count_col]]

# Keep only genes with bait count >= 1
# cat("--- [RUNNING] FILTERING OUT GENES WITH 0 COUNTS ---\n")
# keep_genes <- names(y)[y >= 1]

# y <- y[keep_genes]
# X_sparse <- X_sparse[keep_genes, , drop = FALSE]

# cat(sprintf("--- [INFO] RETAINED %d GENES CON COUNT >= 1 (REMOVED ZEROS) ---\n", length(y)))

# Log-transform the target variable
# cat("--- [RUNNING] LOG-TRANSFORMING TARGET VARIABLE (Y = log(y)) ---\n")
# y <- log(y)

# 5. Grid Search for Alpha
cat("--- [RUNNING] FITTING ELASTIC NET CV MODELS (GRID SEARCH ALPHA) ---\n")
alphas_to_test <- seq(0.1, 0.9, by = 0.1)
all_models <- list()
metrics_list <- list()

best_r_squared <- -Inf
best_alpha <- NA

for (alpha in alphas_to_test) {
  cat(sprintf("--- [INFO] TESTING ALPHA = %.1f ---\n", alpha))
  set.seed(42) # Stesso seed in ogni iterazione per garantire comparabilità dei 10 fold
  
  # nfolds = 10 aumenta la robustezza della cross validation
  cv_fit <- cv.glmnet(x = X_sparse, y = y, alpha = alpha, family = "gaussian", nfolds = 10)
  
  # Estrazione metriche per l'alpha corrente
  lambda_min_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
  r_sq <- cv_fit$glmnet.fit$dev.ratio[lambda_min_idx]
  nz_feat <- cv_fit$nzero[lambda_min_idx]
  
  metrics_list[[length(metrics_list) + 1]] <- data.frame(
    alpha = alpha,
    lambda_min = cv_fit$lambda.min,
    r_squared = r_sq,
    non_zero_features = nz_feat
  )
  
  # Salvataggio in lista
  all_models[[as.character(alpha)]] <- cv_fit
  
  # Check for best model
  if (r_sq > best_r_squared) {
    best_r_squared <- r_sq
    best_alpha <- alpha
  }
}

# Consolidate metrics
metrics_df <- do.call(rbind, metrics_list)
cat(sprintf("--- [INFO] BEST MODEL FOUND AT ALPHA = %.1f WITH R-SQUARED = %.4f ---\n", best_alpha, best_r_squared))

# 6. Extracting the Best Model Data
cat("--- [RUNNING] EXTRACTING METRICS & COEFFICIENTS FOR BEST MODEL ---\n")
best_model <- all_models[[as.character(best_alpha)]]

# Save comprehensive metrics table
write.table(metrics_df, file = metrics_out, row.names = FALSE, sep = "\t", quote = FALSE)

# 7. Extracting Coefficients for Best Model
beta_coefs <- coef(best_model, s = "lambda.min")
coef_df <- data.frame(
  annotation_id = rownames(beta_coefs),
  coefficient = as.numeric(beta_coefs)
)

coef_df <- coef_df[coef_df$annotation_id != "(Intercept)", ]

coef_df <- coef_df[coef_df$coefficient != 0, ]
coef_df <- coef_df[order(abs(coef_df$coefficient), decreasing = TRUE), ]

write_parquet(coef_df, coef_out)

# 8. Saving Full Best Model
cat("--- [RUNNING] SAVING FULL BEST MODEL (.RDS) ---\n")
saveRDS(best_model, file = model_out)

cat("--- [COMPLETE] ELASTIC NET MODEL FITTING WITH GRID SEARCH (OLD) FINISHED SUCCESSFULLY ---\n")