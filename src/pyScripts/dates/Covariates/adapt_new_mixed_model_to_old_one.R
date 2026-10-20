suppressMessages(library(arrow))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))

# ---- Input / output ---------------------------------------------------------
feature_matrix_hdo_path <- snakemake@input[["feature_matrix_hdo"]]
feature_matrix_go_path  <- snakemake@input[["feature_matrix_go"]]
bait_usage_path         <- snakemake@input[["bait_usage"]]

coef_out      <- snakemake@output[["en_coefficients"]]
metrics_out   <- snakemake@output[["en_metrics"]]
model_out     <- snakemake@output[["en_model"]]
cv_curves_out <- snakemake@output[["cv_curves"]]

get_param <- function(name, default) {
  v <- tryCatch(snakemake@params[[name]], error = function(e) NULL)
  if (is.null(v)) default else v
}

min_genes    <- get_param("min_genes", 20)
do_standard  <- get_param("standardize", TRUE)   # StandardScaler-equivalent applied to X before the fit
seed         <- get_param("seed", 42)

alphas_to_test <- c(0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1)
n_folds <- 5

# Gaussian deviance (residual sum of squares)
gauss_dev <- function(y, mu) sum((y - mu)^2)

# Ontology classification from the term's name/ID.
# ADAPT the regexes to your actual IDs (e.g. "GO:0008150", "DOID:1612").
get_ontology <- function(ids) {
  ifelse(grepl("^GO[:_]", ids), "GO",
         ifelse(grepl("^(DOID|HDO)", ids), "HDO", "OTHER"))
}

cat("--- [START] ELASTIC NET (GAUSSIAN ON LOG COUNT, GENES IN COMMON) WITH GRID SEARCH ---\n")

# ---- Load data --------------------------------------------------------------
cat("--- [RUNNING] LOADING DATA ---\n")
bait_count_df <- read.table(bait_usage_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)

# Reads a tab-separated feature matrix and returns a numeric matrix with gene IDs as rownames
read_feature_matrix <- function(path) {
  df <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE)
  if ("entrez_id" %in% colnames(df)) {
    genes <- as.character(df$entrez_id)
    df$entrez_id <- NULL
  } else if ("__index_level_0__" %in% colnames(df)) {
    genes <- as.character(df[["__index_level_0__"]])
    df[["__index_level_0__"]] <- NULL
  } else if (.row_names_info(df) > 0) {
    # header has one field fewer than the rows: R used the first column as rownames
    genes <- rownames(df)
  } else {
    genes <- as.character(df[[1]])
    df[[1]] <- NULL
  }
  if (anyDuplicated(genes)) stop(sprintf("Duplicated gene IDs in feature matrix: %s", path))
  m <- as.matrix(df) * 1
  rownames(m) <- genes
  m
}

# Pads a matrix with zero rows so it covers all genes in `all_genes`
expand_rows <- function(m, all_genes) {
  out <- matrix(0, nrow = length(all_genes), ncol = ncol(m),
                dimnames = list(all_genes, colnames(m)))
  out[rownames(m), ] <- m
  out
}

fm1 <- read_feature_matrix(feature_matrix_hdo_path)
fm2 <- read_feature_matrix(feature_matrix_go_path)
if (length(intersect(colnames(fm1), colnames(fm2))) > 0) {
  stop("Feature matrices share term IDs")
}

# Full outer join on genes (unique IDs), all terms from both matrices; missing = 0
all_genes <- union(rownames(fm1), rownames(fm2))
feature_matrix <- cbind(expand_rows(fm1, all_genes),
                        expand_rows(fm2, all_genes))
rm(fm1, fm2); invisible(gc())
cat(sprintf("--- [INFO] JOINED FEATURE MATRICES: %d genes x %d terms ---\n",
            nrow(feature_matrix), ncol(feature_matrix)))

ont_all <- get_ontology(colnames(feature_matrix))
for (ont in unique(ont_all)) {
  cat(sprintf("--- [INFO] %-5s: %d terms ---\n", ont, sum(ont_all == ont)))
}

# ---- Genes in common with the bait table, target y = log(count) ------------
cat("--- [RUNNING] FINDING GENES IN COMMON AND ALIGNING TARGET VARIABLE (Y) ---\n")
bait_id_col <- ifelse("entrez_id_bait" %in% colnames(bait_count_df),
                      "entrez_id_bait", colnames(bait_count_df)[1])
count_col   <- ifelse("count" %in% colnames(bait_count_df),
                      "count", colnames(bait_count_df)[2])

bait_agg <- bait_count_df %>%
  transmute(gene = as.character(.data[[bait_id_col]]),
            count = .data[[count_col]]) %>%
  group_by(gene) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

n_dup <- nrow(bait_count_df) - nrow(bait_agg)
in_matrix <- bait_agg$gene %in% rownames(feature_matrix)
cat(sprintf("--- [INFO] BAIT TABLE: %d rows, %d duplicates aggregated, %d genes not in feature matrix (dropped) ---\n",
            nrow(bait_count_df), n_dup, sum(!in_matrix)))

universe_genes <- rownames(feature_matrix)[rownames(feature_matrix) %in% bait_agg$gene]
n_genes <- length(universe_genes)
if (n_genes == 0) stop("No genes in common between feature matrix and bait table")
cat(sprintf("--- [INFO] %d GENES IN COMMON FOUND ---\n", n_genes))

X_raw <- feature_matrix[universe_genes, , drop = FALSE]
rm(feature_matrix); invisible(gc())

# ---- Feature filtering (terms with >= min_genes genes among genes in common) --
cat("--- [RUNNING] FILTERING FEATURE MATRIX ---\n")
col_sums <- colSums(X_raw)
valid_cols <- col_sums >= min_genes
ont_all <- get_ontology(colnames(X_raw))
for (ont in unique(ont_all)) {
  cat(sprintf("--- [INFO] %-5s: retained %d / %d terms (min >= %d genes) ---\n",
              ont, sum(valid_cols & ont_all == ont), sum(ont_all == ont),
              min_genes))
}
if (sum(valid_cols) == 0) stop("No annotation left after filtering")

X_raw <- X_raw[, valid_cols, drop = FALSE]
term_n_genes <- col_sums[valid_cols]   # how many genes support each term, kept for diagnostics
                                       # in the coefficients output (named by term id)

y <- log(setNames(bait_agg$count, bait_agg$gene)[universe_genes])
cat(sprintf("--- [INFO] Y: %d genes; max count = %g; median count = %.1f ---\n",
            length(y), max(exp(y)), median(exp(y))))

# ---- Scaling (equivalent to sklearn StandardScaler: population sd, zero-variance scale = 1) ----
if (do_standard) {
  cat("--- [RUNNING] SCALING DATA ---\n")
  mu  <- colMeans(X_raw)
  sdv <- sqrt(colMeans(sweep(X_raw, 2, mu, "-")^2))
  sdv[sdv < 10 * .Machine$double.eps] <- 1
  X_full <- sweep(sweep(X_raw, 2, mu, "-"), 2, sdv, "/")
} else {
  sdv <- setNames(rep(1, ncol(X_raw)), colnames(X_raw))
  X_full <- X_raw
}
rm(X_raw); invisible(gc())

# ---- Folds: contiguous, not shuffled (same as sklearn KFold(n_splits = 5)) ----
fold_sizes <- rep(n_genes %/% n_folds, n_folds)
if (n_genes %% n_folds > 0) {
  idx_extra <- seq_len(n_genes %% n_folds)
  fold_sizes[idx_extra] <- fold_sizes[idx_extra] + 1
}
foldid <- rep(seq_len(n_folds), times = fold_sizes)

# ---- Grid search on alpha (CV-based selection, fixed foldid) ---------------
cat("--- [RUNNING] FITTING ELASTIC NET CV MODELS (GRID SEARCH ALPHA) ---\n")
all_models <- list()
metrics_list <- list()

# Null-model CV deviance, on the same folds. Reference for deviance explained.
cv_null <- mean(sapply(seq_len(n_folds), function(k) {
  tr <- foldid != k
  gauss_dev(y[!tr], rep(mean(y[tr]), sum(!tr))) / sum(!tr)
}))

for (alpha in alphas_to_test) {
  cat(sprintf("--- [INFO] TESTING ALPHA = %.2f ---\n", alpha))
  cv_fit <- cv.glmnet(x = X_full, y = y, family = "gaussian",
                      alpha = alpha, foldid = foldid,
                      lambda.min.ratio = 1e-3, standardize = FALSE)

  i_min <- match(cv_fit$lambda.min, cv_fit$lambda)
  i_1se <- match(cv_fit$lambda.1se, cv_fit$lambda)

  # Out-of-fold deviance explained, relative to the null model (cv_null)
  cv_dev_expl_min <- 1 - cv_fit$cvm[i_min] / cv_null
  cv_dev_expl_1se <- 1 - cv_fit$cvm[i_1se] / cv_null

  metrics_list[[length(metrics_list) + 1]] <- data.frame(
    alpha = alpha,
    lambda_min = cv_fit$lambda.min,
    lambda_1se = cv_fit$lambda.1se,
    cv_dev_explained_min = cv_dev_expl_min,
    cv_dev_explained_1se = cv_dev_expl_1se,
    cvm_min = cv_fit$cvm[i_min],
    cvsd_min = cv_fit$cvsd[i_min],
    nonzero_min = cv_fit$nzero[i_min],
    nonzero_1se = cv_fit$nzero[i_1se]
  )
  all_models[[as.character(alpha)]] <- cv_fit
}

metrics_df <- do.call(rbind, metrics_list)
best_row   <- which.max(metrics_df$cv_dev_explained_min)   # selection criterion: CV
best_alpha <- metrics_df$alpha[best_row]
cat(sprintf("--- [INFO] BEST ALPHA (by CV) = %.2f | CV dev. explained = %.5f ---\n",
            best_alpha, metrics_df$cv_dev_explained_min[best_row]))

write.table(metrics_df, file = metrics_out, row.names = FALSE, sep = "\t", quote = FALSE)

# ---- Coefficients: best model only -------------------------------------------
final_cv <- all_models[[as.character(best_alpha)]]

cat(sprintf("--- [RUNNING] EXTRACTING COEFFICIENTS (alpha = %.2f) ---\n", best_alpha))
cf_min <- as.numeric(coef(final_cv, s = "lambda.min"))
cf_1se <- as.numeric(coef(final_cv, s = "lambda.1se"))
term_ids <- rownames(coef(final_cv, s = "lambda.min"))

coef_df <- data.frame(
  alpha = best_alpha,
  annotation_id = term_ids,
  ontology = get_ontology(term_ids),
  n_genes_term = ifelse(term_ids %in% names(term_n_genes),
                        unname(term_n_genes[term_ids]), NA_integer_),
  coefficient = cf_min,
  coefficient_original_scale = unname(cf_min / sdv[term_ids]),
  rate_ratio = exp(cf_min),
  coefficient_1se = cf_1se,
  stringsAsFactors = FALSE
)
coef_df <- coef_df[coef_df$annotation_id != "(Intercept)" &
                     (coef_df$coefficient != 0 | coef_df$coefficient_1se != 0), ]
coef_df <- coef_df[order(abs(coef_df$coefficient), decreasing = TRUE), ]
rownames(coef_df) <- NULL

cat(sprintf("--- [INFO] NON-ZERO TERMS BY ONTOLOGY (lambda.min, alpha = %.2f) ---\n", best_alpha))
print(table(coef_df$ontology[coef_df$coefficient != 0]))

write_parquet(coef_df, coef_out)

# ---- Model saving -------------------------------------------------------------
cat("--- [RUNNING] SAVING MODELS (.RDS) ---\n")
saveRDS(list(
  final_cv = final_cv,             # best-alpha model (reported coefficients)
  best_alpha = best_alpha,
  cv_models = all_models,          # grid models
  seed = seed
), file = model_out)

# ---- CV curves ------------------------------------------------------------------
cat("--- [RUNNING] SAVING CV CURVES (.PDF) ---\n")
pdf(cv_curves_out, width = 8, height = 6)

plot(metrics_df$alpha, metrics_df$cv_dev_explained_min, type = "b", pch = 16,
     ylim = range(c(metrics_df$cv_dev_explained_min, 0), na.rm = TRUE),
     xlab = "alpha", ylab = "Deviance explained",
     main = "CV deviance explained")
abline(h = 0, lty = 3)

for (alpha_val in names(all_models)) {
  plot(all_models[[alpha_val]])
  title(main = sprintf("Cross-Validation Curve for Alpha = %s", alpha_val), line = 3)
}
dev.off()

cat("--- [COMPLETE] ELASTIC NET MODEL FITTING FINISHED SUCCESSFULLY ---\n")