suppressMessages(library(arrow))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages(library(Matrix))

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

min_genes    <- snakemake@params[["min_genes"]]
use_covar    <- get_param("use_annot_covariate", TRUE)   # compute n_annot, baseline GLM, AND add it to the elastic net fits
do_standard  <- get_param("standardize", TRUE)
seed         <- get_param("seed", 42)

alphas_to_test <- c(seq(0.1, 0.9, by = 0.1), 1)
n_folds <- 10

# Gaussian deviance (residual sum of squares) and deviance explained
gauss_dev <- function(y, mu) sum((y - mu)^2)
dev_explained <- function(y_obs, mu_pred, mu_null) {
  1 - gauss_dev(y_obs, mu_pred) / gauss_dev(y_obs, rep(mu_null, length(y_obs)))
}

# Ontology classification from the term's name/ID.
# ADAPT the regexes to your actual IDs (e.g. "GO:0008150", "DOID:1612").
get_ontology <- function(ids) {
  ifelse(grepl("^GO[:_]", ids), "GO",
         ifelse(grepl("^(DOID|HDO)", ids), "HDO", "OTHER"))
}

cat("--- [START] ELASTIC NET (GAUSSIAN ON LOG COUNT, TESTED GENES) WITH GRID SEARCH ---\n")

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
universe_genes <- union(rownames(fm1), rownames(fm2))
feature_matrix_df <- cbind(expand_rows(fm1, universe_genes),
                           expand_rows(fm2, universe_genes))
rm(fm1, fm2); invisible(gc())
n_genes <- length(universe_genes)
cat(sprintf("--- [INFO] JOINED FEATURE MATRICES: %d genes x %d terms ---\n",
            n_genes, ncol(feature_matrix_df)))

# ---- Feature filtering -------------------------------------------------------
cat("--- [RUNNING] FILTERING FEATURE MATRIX ---\n")
col_sums <- colSums(feature_matrix_df)
row_annot_all <- rowSums(feature_matrix_df)           # n. annotations per gene (pre-filter)

valid_cols <- col_sums >= min_genes
ont_all <- get_ontology(colnames(feature_matrix_df))
term_n_genes <- col_sums  # how many genes support each term, kept for diagnostics
                          # in the coefficients output (named by term id)

for (ont in unique(ont_all)) {
  cat(sprintf("--- [INFO] %-5s: retained %d / %d terms (min >= %d genes) ---\n",
              ont, sum(valid_cols & ont_all == ont), sum(ont_all == ont),
              min_genes))
}
if (sum(valid_cols) == 0) stop("No annotation left after filtering")

feature_matrix_df <- feature_matrix_df[, valid_cols, drop = FALSE]

# ---- Sparse matrix (numeric, dgCMatrix) --------------------------------------
cat("--- [RUNNING] CONVERTING TO SPARSE MATRIX ---\n")
X_sparse <- Matrix(as.matrix(feature_matrix_df) * 1, sparse = TRUE)
X_sparse <- as(as(X_sparse, "generalMatrix"), "CsparseMatrix")
rownames(X_sparse) <- universe_genes
rm(feature_matrix_df); invisible(gc())

# n_annot: log(1 + gene's number of annotations). Used for the standalone
# baseline GLM below (glm(y ~ n_annot, ...)) AND, in THIS script, also added
# unpenalized to X_full - i.e. it enters the elastic net fits themselves.
if (use_covar) {
  n_annot <- log1p(row_annot_all)
  X_full <- cbind(X_sparse, Matrix(n_annot, ncol = 1, sparse = TRUE))
  colnames(X_full) <- c(colnames(X_sparse), "n_annot")
  pen_factor <- c(rep(1, ncol(X_sparse)), 0)
  cat("--- [INFO] UNPENALIZED COVARIATE 'n_annot' ADDED TO THE ELASTIC NET FITS ---\n")
} else {
  n_annot <- NULL
  X_full <- X_sparse
  pen_factor <- rep(1, ncol(X_sparse))
}
rownames(X_full) <- universe_genes

# ---- Target y ---------------------------------------------------------------
cat("--- [RUNNING] ALIGNING TARGET VARIABLE (Y) ---\n")
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
in_universe <- bait_agg$gene %in% universe_genes
cat(sprintf("--- [INFO] BAIT TABLE: %d rows, %d duplicates aggregated, %d genes not in universe (dropped) ---\n",
            nrow(bait_count_df), n_dup, sum(!in_universe)))

bait_agg <- bait_agg[in_universe, ]
y <- setNames(rep(0, n_genes), universe_genes)
y[bait_agg$gene] <- bait_agg$count
if (any(y != round(y))) cat("--- [WARN] NON-INTEGER COUNTS: ROUNDING ---\n")
y <- round(y)

cat(sprintf("--- [INFO] Y: %d genes with count > 0, %d zeros (%.1f%%); max = %d; median(>0) = %.1f ---\n",
            sum(y > 0), sum(y == 0), 100 * mean(y == 0), max(y),
            ifelse(any(y > 0), median(y[y > 0]), NA)))

# ---- Restrict to tested genes (y > 0) and model log(y) -----------------------
keep <- y > 0
y <- log(y[keep])
X_full <- X_full[keep, , drop = FALSE]
if (use_covar) n_annot <- n_annot[keep]
universe_genes <- universe_genes[keep]
n_genes <- length(universe_genes)

# Re-apply min_genes among tested genes only (n_annot is always kept)
is_term <- colnames(X_full) != "n_annot"
keep_cols <- !is_term | colSums(X_full) >= min_genes
X_full <- X_full[, keep_cols, drop = FALSE]
pen_factor <- pen_factor[keep_cols]
cat(sprintf("--- [INFO] TESTED GENES ONLY: %d genes, %d terms with >= %d tested genes ---\n",
            n_genes, sum(keep_cols & is_term), min_genes))

# Fixed foldid: same for every alpha
set.seed(seed)
foldid <- sample(rep(seq_len(n_folds), length.out = n_genes))

# ---- Baseline: n_annot only --------------------------------------------------
baseline <- NULL
baseline_dev_expl <- NA
if (use_covar) {
  base_df <- data.frame(y = y, n_annot = n_annot)
  baseline <- glm(y ~ n_annot, data = base_df, family = gaussian())
  baseline_dev_expl <- dev_explained(y, fitted(baseline), mean(y))
  cat(sprintf("--- [INFO] BASELINE (n_annot only) DEVIANCE EXPLAINED = %.4f ---\n",
              baseline_dev_expl))
}

# ---- Grid search on alpha (CV-based selection, fixed foldid) ---------------
cat("--- [RUNNING] FITTING ELASTIC NET CV MODELS (GRID SEARCH ALPHA) ---\n")
all_models <- list()
metrics_list <- list()

# Null-model CV deviance, on the same folds. Used as the reference for
# deviance explained. Computed independently (not from cv_fit$cvm[1]) because
# n_annot is unpenalized and stays in the model even at the first lambda, so
# cv_fit$cvm[1] would NOT be the true null model here.
cv_null <- mean(sapply(seq_len(n_folds), function(k) {
  tr <- foldid != k
  gauss_dev(y[!tr], rep(mean(y[tr]), sum(!tr))) / sum(!tr)
}))

for (alpha in alphas_to_test) {
  cat(sprintf("--- [INFO] TESTING ALPHA = %.2f ---\n", alpha))
  cv_fit <- cv.glmnet(x = X_full, y = y, family = "gaussian",
                      alpha = alpha, foldid = foldid,
                      penalty.factor = pen_factor, standardize = do_standard)

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
    nonzero_1se = cv_fit$nzero[i_1se],
    baseline_dev_explained = baseline_dev_expl
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
# The best-alpha model is the CV fit already computed on all genes in the grid.
final_cv <- all_models[[as.character(best_alpha)]]

cat(sprintf("--- [RUNNING] EXTRACTING COEFFICIENTS (alpha = %.2f) ---\n", best_alpha))
cf_min <- as.numeric(coef(final_cv, s = "lambda.min"))
cf_1se <- as.numeric(coef(final_cv, s = "lambda.1se"))
term_ids <- rownames(coef(final_cv, s = "lambda.min"))

coef_df <- data.frame(
  alpha = best_alpha,
  annotation_id = term_ids,
  ontology = ifelse(term_ids == "n_annot", "COVARIATE", get_ontology(term_ids)),
  n_genes_term = ifelse(term_ids %in% names(term_n_genes),
                        unname(term_n_genes[term_ids]), NA_integer_),
  coefficient = cf_min,
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
  baseline = baseline,
  penalty_factor = pen_factor,
  seed = seed
), file = model_out)

# ---- CV curves ------------------------------------------------------------------
cat("--- [RUNNING] SAVING CV CURVES (.PDF) ---\n")
pdf(cv_curves_out, width = 8, height = 6)

plot(metrics_df$alpha, metrics_df$cv_dev_explained_min, type = "b", pch = 16,
     ylim = range(c(metrics_df$cv_dev_explained_min, baseline_dev_expl, 0), na.rm = TRUE),
     xlab = "alpha", ylab = "Deviance explained",
     main = "CV deviance explained")
if (!is.null(baseline)) abline(h = baseline_dev_expl, lty = 2, col = "grey40")
abline(h = 0, lty = 3)
if (!is.null(baseline)) {
  legend("bottomright", legend = c("CV", "Baseline n_annot"),
         pch = c(16, NA), lty = c(1, 2), col = c("black", "grey40"), bty = "n")
}

for (alpha_val in names(all_models)) {
  plot(all_models[[alpha_val]])
  title(main = sprintf("Cross-Validation Curve for Alpha = %s", alpha_val), line = 3)
}
dev.off()

cat("--- [COMPLETE] ELASTIC NET MODEL FITTING FINISHED SUCCESSFULLY ---\n")