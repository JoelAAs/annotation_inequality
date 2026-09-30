suppressMessages(library(arrow))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages(library(Matrix))

# ---- Input / output ---------------------------------------------------------
feature_matrix_path <- snakemake@input[["feature_matrix"]]
bait_usage_path     <- snakemake@input[["bait_usage"]]

coef_out      <- snakemake@output[["en_coefficients"]]
metrics_out   <- snakemake@output[["en_metrics"]]
model_out     <- snakemake@output[["en_model"]]
cv_curves_out <- snakemake@output[["cv_curves"]]

get_param <- function(name, default) {
  v <- tryCatch(snakemake@params[[name]], error = function(e) NULL)
  if (is.null(v)) default else v
}

min_genes    <- snakemake@params[["min_genes"]]
holdout_frac <- get_param("holdout_frac", 0.20)
n_boot       <- get_param("n_boot", 30)
use_covar    <- get_param("use_annot_covariate", TRUE)   # compute n_annot, baseline GLM, AND add it to the elastic net fits
do_standard  <- get_param("standardize", TRUE)
seed         <- get_param("seed", 42)
compare_alphas_raw <- get_param("compare_alphas", numeric(0))  # extra alphas to fully refit besides best_alpha

alphas_to_test <- c(0.05, seq(0.1, 0.9, by = 0.1), 1)
n_folds <- 10

# compare_alphas = "all" expands to the full alpha grid (alphas_to_test),
# so every tested alpha gets its coefficients fully refit and exported.
if (identical(compare_alphas_raw, "all")) {
  compare_alphas <- alphas_to_test
  cat(sprintf("--- [INFO] compare_alphas = 'all': will fully refit all %d alphas in the grid ---\n",
              length(alphas_to_test)))
} else {
  compare_alphas <- as.numeric(compare_alphas_raw)
}

# Gaussian deviance (residual sum of squares) and deviance explained on new data
gauss_dev <- function(y, mu) sum((y - mu)^2)
dev_explained <- function(y_new, mu_new, mu_null) {
  1 - gauss_dev(y_new, mu_new) / gauss_dev(y_new, rep(mu_null, length(y_new)))
}

# Ontology classification from the term's name/ID.
# ADAPT the regexes to your actual IDs (e.g. "GO:0008150", "DOID:1612").
get_ontology <- function(ids) {
  ifelse(grepl("^GO[:_]", ids), "GO",
         ifelse(grepl("^(DOID|HDO)", ids), "HDO", "OTHER"))
}

cat("--- [START] ELASTIC NET (GAUSSIAN ON LOG COUNT, TESTED GENES) WITH GRID SEARCH, HOLD-OUT AND STABILITY ---\n")

# ---- Load data --------------------------------------------------------------
cat("--- [RUNNING] LOADING DATA ---\n")
feature_matrix_df <- read_parquet(feature_matrix_path)
bait_count_df <- read.table(bait_usage_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)

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
n_genes <- length(universe_genes)
if (anyDuplicated(universe_genes)) stop("Duplicated gene IDs in feature matrix")

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

feature_matrix_df <- feature_matrix_df[, valid_cols]

# ---- Sparse matrix (numeric, dgCMatrix) --------------------------------------
cat("--- [RUNNING] CONVERTING TO SPARSE MATRIX ---\n")
X_sparse <- Matrix(as.matrix(feature_matrix_df) * 1, sparse = TRUE)
X_sparse <- as(as(X_sparse, "generalMatrix"), "CsparseMatrix")
rownames(X_sparse) <- universe_genes
rm(feature_matrix_df); invisible(gc())

# n_annot: log(1 + gene's number of annotations). Used for the standalone
# baseline GLM below (glm(y ~ n_annot, ...)) AND, in THIS script, also added
# unpenalized to X_full - i.e. it enters the elastic net fits themselves
# (grid search, final refit, stability selection). This is the "with n_annot"
# variant, meant to be compared against the version that fits on the
# annotation terms alone (n_annot used only for the baseline).
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

# ---- Hold-out split ----------------------------------------------------------
cat("--- [RUNNING] HOLD-OUT SPLIT ---\n")
set.seed(seed)
test_idx <- sample(seq_len(n_genes), size = round(holdout_frac * n_genes))
train_idx <- setdiff(seq_len(n_genes), test_idx)

X_train <- X_full[train_idx, , drop = FALSE]; y_train <- y[train_idx]
X_test  <- X_full[test_idx,  , drop = FALSE]; y_test  <- y[test_idx]
mu_null <- mean(y_train)
cat(sprintf("--- [INFO] TRAIN: %d genes | TEST: %d genes ---\n",
            length(train_idx), length(test_idx)))

# Fixed foldid: same for every alpha
set.seed(seed)
foldid <- sample(rep(seq_len(n_folds), length.out = length(train_idx)))

# ---- Baseline: n_annot only --------------------------------------------------
baseline <- NULL
if (use_covar) {
  base_df <- data.frame(y = y_train, n_annot = n_annot[train_idx])
  baseline <- glm(y ~ n_annot, data = base_df, family = gaussian())
  base_pred <- predict(baseline,
                       newdata = data.frame(n_annot = n_annot[test_idx]),
                       type = "response")
  baseline_dev_expl <- dev_explained(y_test, base_pred, mu_null)
  cat(sprintf("--- [INFO] BASELINE (n_annot only) TEST DEVIANCE EXPLAINED = %.4f ---\n",
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
  gauss_dev(y_train[!tr], rep(mean(y_train[tr]), sum(!tr))) / sum(!tr)
}))

for (alpha in alphas_to_test) {
  cat(sprintf("--- [INFO] TESTING ALPHA = %.2f ---\n", alpha))
  cv_fit <- cv.glmnet(x = X_train, y = y_train, family = "gaussian",
                      alpha = alpha, foldid = foldid,
                      penalty.factor = pen_factor, standardize = do_standard)

  i_min <- match(cv_fit$lambda.min, cv_fit$lambda)
  i_1se <- match(cv_fit$lambda.1se, cv_fit$lambda)

  # Out-of-fold deviance explained, relative to the null model (cv_null)
  cv_dev_expl_min <- 1 - cv_fit$cvm[i_min] / cv_null
  cv_dev_expl_1se <- 1 - cv_fit$cvm[i_1se] / cv_null

  # Hold-out performance (informative only: NOT used to choose alpha)
  pred_test <- as.numeric(predict(cv_fit, newx = X_test, s = "lambda.min",
                                  type = "response"))
  test_dev_expl <- dev_explained(y_test, pred_test, mu_null)
  test_spearman <- suppressWarnings(cor(pred_test, y_test, method = "spearman"))

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
    test_dev_explained = test_dev_expl,
    test_spearman = test_spearman,
    baseline_test_dev_explained = ifelse(is.null(baseline), NA, baseline_dev_expl)
  )
  all_models[[as.character(alpha)]] <- cv_fit
}

metrics_df <- do.call(rbind, metrics_list)
best_row   <- which.max(metrics_df$cv_dev_explained_min)   # selection criterion: CV
best_alpha <- metrics_df$alpha[best_row]
cat(sprintf("--- [INFO] BEST ALPHA (by CV) = %.2f | CV dev. explained = %.5f | hold-out dev. explained = %.5f ---\n",
            best_alpha, metrics_df$cv_dev_explained_min[best_row],
            metrics_df$test_dev_explained[best_row]))
if (metrics_df$test_dev_explained[best_row] <= 0) {
  cat("--- [WARN] HOLD-OUT DEVIANCE EXPLAINED <= 0: THE MODEL DOES NOT GENERALIZE, INTERPRET COEFFICIENTS WITH CAUTION ---\n")
}

write.table(metrics_df, file = metrics_out, row.names = FALSE, sep = "\t", quote = FALSE)

# ---- Final model + stability + coefficients, as a reusable function per alpha ----
# Refits cv.glmnet on ALL genes for a given alpha and runs stability selection
# via subsampling. Returns a list(final_cv, coef_df); nothing is written here -
# the coefficients of every alpha are combined and written once to coef_out,
# since the Snakemake rule only declares a single coefficients output.
fit_alpha <- function(alpha_value) {
  cat(sprintf("--- [RUNNING] REFITTING FINAL MODEL ON ALL GENES (alpha = %.2f) ---\n",
              alpha_value))
  set.seed(seed)
  foldid_all <- sample(rep(seq_len(n_folds), length.out = n_genes))
  final_cv <- cv.glmnet(x = X_full, y = y, family = "gaussian", alpha = alpha_value,
                        foldid = foldid_all, penalty.factor = pen_factor,
                        standardize = do_standard)

  # Stability: subsampling (80% of genes), coefficient interpolated at lambda.min
  stab_freq <- setNames(rep(NA_real_, ncol(X_full)), colnames(X_full))
  if (n_boot > 0) {
    cat(sprintf("--- [RUNNING] STABILITY SELECTION (%d SUBSAMPLES, alpha = %.2f) ---\n",
                n_boot, alpha_value))
    sel_counts <- setNames(numeric(ncol(X_full)), colnames(X_full))
    for (b in seq_len(n_boot)) {
      set.seed(seed + b)
      idx <- sample(seq_len(n_genes), size = round(0.8 * n_genes))
      # NB: do not pass a single lambda to glmnet() - without the full warm-start
      # path, convergence for family="poisson" is unreliable and can return
      # degenerate solutions (all coefficients zero). Let it compute the whole
      # sequence and interpolate the value at lambda.min via s=.
      fit_b <- glmnet(X_full[idx, , drop = FALSE], y[idx], family = "gaussian",
                      alpha = alpha_value,
                      penalty.factor = pen_factor, standardize = do_standard)
      cb <- as.numeric(coef(fit_b, s = final_cv$lambda.min))[-1]
      cat(sprintf("--- [INFO] STABILITY BOOT %d/%d (alpha = %.2f): %d / %d nonzero coefficients ---\n",
                  b, n_boot, alpha_value, sum(cb != 0), length(cb)))
      sel_counts <- sel_counts + (cb != 0)
    }
    stab_freq <- sel_counts / n_boot
  }

  # Coefficients
  cat(sprintf("--- [RUNNING] EXTRACTING COEFFICIENTS (alpha = %.2f) ---\n", alpha_value))
  cf_min <- as.numeric(coef(final_cv, s = "lambda.min"))
  cf_1se <- as.numeric(coef(final_cv, s = "lambda.1se"))
  term_ids <- rownames(coef(final_cv, s = "lambda.min"))

  coef_df <- data.frame(
    alpha = alpha_value,
    annotation_id = term_ids,
    ontology = ifelse(term_ids == "n_annot", "COVARIATE", get_ontology(term_ids)),
    n_genes_term = ifelse(term_ids %in% names(term_n_genes),
                          unname(term_n_genes[term_ids]), NA_integer_),
    coefficient = cf_min,
    rate_ratio = exp(cf_min),
    coefficient_1se = cf_1se,
    selection_freq = unname(stab_freq[term_ids]),
    stringsAsFactors = FALSE
  )
  coef_df <- coef_df[coef_df$annotation_id != "(Intercept)" &
                       (coef_df$coefficient != 0 | coef_df$coefficient_1se != 0), ]
  coef_df <- coef_df[order(abs(coef_df$coefficient), decreasing = TRUE), ]
  rownames(coef_df) <- NULL

  cat(sprintf("--- [INFO] NON-ZERO TERMS BY ONTOLOGY (lambda.min, alpha = %.2f) ---\n", alpha_value))
  print(table(coef_df$ontology[coef_df$coefficient != 0]))

  list(final_cv = final_cv, coef_df = coef_df)
}

# best_alpha is always included; compare_alphas adds more rows (tagged by the
# "alpha" column) to the SAME coefficients file, since the Snakemake rule
# declares a single en_coefficients output.
fit_best <- fit_alpha(best_alpha)
final_cv <- fit_best$final_cv          # kept for saveRDS() below
all_coef_list <- list(fit_best$coef_df)

extra_alphas <- setdiff(compare_alphas, best_alpha)
for (a in extra_alphas) {
  all_coef_list[[length(all_coef_list) + 1]] <- fit_alpha(a)$coef_df
}

coef_df <- do.call(rbind, all_coef_list)
write_parquet(coef_df, coef_out)

# ---- Model saving -------------------------------------------------------------
cat("--- [RUNNING] SAVING MODELS (.RDS) ---\n")
saveRDS(list(
  final_cv = final_cv,             # refit on all genes (reported coefficients)
  best_alpha = best_alpha,
  cv_train_models = all_models,    # grid models (train only)
  baseline = baseline,
  train_idx = train_idx,
  test_idx = test_idx,
  penalty_factor = pen_factor,
  seed = seed
), file = model_out)

# ---- CV curves ------------------------------------------------------------------
cat("--- [RUNNING] SAVING CV CURVES (.PDF) ---\n")
pdf(cv_curves_out, width = 8, height = 6)

plot(metrics_df$alpha, metrics_df$cv_dev_explained_min, type = "b", pch = 16,
     ylim = range(c(metrics_df$cv_dev_explained_min, metrics_df$test_dev_explained,
                    baseline_dev_expl <- if (is.null(baseline)) NULL else metrics_df$baseline_test_dev_explained[1],
                    0), na.rm = TRUE),
     xlab = "alpha", ylab = "Deviance explained",
     main = "CV (train) vs hold-out deviance explained")
lines(metrics_df$alpha, metrics_df$test_dev_explained, type = "b", pch = 17, col = "firebrick")
if (!is.null(baseline)) abline(h = metrics_df$baseline_test_dev_explained[1], lty = 2, col = "grey40")
abline(h = 0, lty = 3)
legend("bottomright", legend = c("CV (train)", "Hold-out", "Baseline n_annot"),
       pch = c(16, 17, NA), lty = c(1, 1, 2), col = c("black", "firebrick", "grey40"), bty = "n")

for (alpha_val in names(all_models)) {
  plot(all_models[[alpha_val]])
  title(main = sprintf("Cross-Validation Curve for Alpha = %s", alpha_val), line = 3)
}
dev.off()

cat("--- [COMPLETE] ELASTIC NET MODEL FITTING FINISHED SUCCESSFULLY ---\n")