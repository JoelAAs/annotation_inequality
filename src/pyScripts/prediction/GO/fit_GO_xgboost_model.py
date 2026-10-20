import os
import gc
import json
import joblib
import numpy as np
import pandas as pd
import xgboost as xgb
import shap
from sklearn.metrics import roc_auc_score, average_precision_score, log_loss
from sklearn.model_selection import StratifiedKFold

# --- SETTINGS ---
RANDOM_STATE = 1234
N_CV_FOLDS = 5

def balance_indices_1to1(y, seed):
    """
    Returns an array of indices that perfectly balances the 0s to match the 1s.
    """
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    n_pos = len(pos_idx)
    
    if n_pos == 0:
        return np.arange(len(y))
        
    rng = np.random.default_rng(seed)
    neg_sampled_idx = rng.choice(neg_idx, size=n_pos, replace=False)
    
    balanced_idx = np.concatenate([pos_idx, neg_sampled_idx])
    rng.shuffle(balanced_idx)
    return balanced_idx


if __name__ == "__main__":
    # --- SNAKEMAKE I/O ---
    matrix_file = snakemake.input.full_matrix
    y_file      = snakemake.input.y_file
    params_file = snakemake.input.best_params
    rest_file   = snakemake.input.rest_ids
    test_file   = snakemake.input.test_ids

    out_model           = snakemake.output.saved_model
    out_metrics         = snakemake.output.test_metrics
    out_predictions     = snakemake.output.test_predictions
    out_shap_values     = snakemake.output.shap_values
    out_x_test_features = snakemake.output.x_test_features

    threads = snakemake.threads
    target_term = snakemake.wildcards.term.replace("_", ":")

    for path in (out_model, out_metrics, out_predictions, out_shap_values, out_x_test_features):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    print(f"\n=======================================================", flush=True)
    print(f"--- [START] Final Fit & Evaluation for {target_term}", flush=True)
    print(f"=======================================================\n", flush=True)

    # --- LOAD PARAMETERS ---
    with open(params_file, "r") as f:
        best_params_payload = json.load(f)
        best_params = best_params_payload["params"]

    # --- LOAD IDS & LABELS ---
    rest_ids = pd.read_csv(rest_file)["Gene_ID"].astype(str).tolist()
    test_ids = pd.read_csv(test_file)["Gene_ID"].astype(str).tolist()

    y_df = pd.read_parquet(y_file)
    y_df["Gene_ID"] = y_df["Gene_ID"].astype(str)
    label_map = dict(zip(y_df["Gene_ID"], y_df["Label"]))

    y_rest = np.array([label_map[g] for g in rest_ids], dtype=np.int8)
    y_test = np.array([label_map[g] for g in test_ids], dtype=np.int8)

    # --- LOAD & SLICE MATRIX ---
    print(f"--- [{target_term}] Loading Matrix and slicing Rest/Test sets...", flush=True)
    X_df = pd.read_parquet(matrix_file)
    X_df["Gene_ID"] = X_df["Gene_ID"].astype(str)
    X_df.set_index("Gene_ID", inplace=True)

    # Extract X_rest as numpy array for training efficiency
    X_rest = X_df.loc[rest_ids].to_numpy(dtype=np.float32, copy=True)
    
    # Extract X_test as DataFrame to retain Entrez ID column names for SHAP and plotting
    X_test_df = X_df.loc[test_ids].copy()
    
    del X_df
    gc.collect()

    print(f"--- [{target_term}] -> X_rest shape: {X_rest.shape}", flush=True)
    print(f"--- [{target_term}] -> X_test shape: {X_test_df.shape}", flush=True)


    # --- STEP 1: TRAIN MASTER MODEL ON 1:1 REST ---
    print(f"\n--- [{target_term}] Step 1: Training Master Model on 1:1 balanced Rest dataset...", flush=True)
    bal_rest_idx = balance_indices_1to1(y_rest, seed=RANDOM_STATE)
    X_rest_bal = X_rest[bal_rest_idx]
    y_rest_bal = y_rest[bal_rest_idx]

    master_model = xgb.XGBClassifier(
        **best_params,
        gamma=0,
        eval_metric="logloss",
        random_state=RANDOM_STATE,
        n_jobs=threads,
        verbosity=0,
    )
    master_model.fit(X_rest_bal, y_rest_bal)


    # --- STEP 2: CROSS-VALIDATE ON 1:1 REST (AUROC) ---
    print(f"--- [{target_term}] Step 2: 5-Fold Cross-Validation on 1:1 balanced Rest...", flush=True)
    skf = StratifiedKFold(n_splits=N_CV_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    cv_auroc_scores = []

    for fold, (tr_idx, va_idx) in enumerate(skf.split(X_rest_bal, y_rest_bal), 1):
        X_tr, y_tr = X_rest_bal[tr_idx], y_rest_bal[tr_idx]
        X_va, y_va = X_rest_bal[va_idx], y_rest_bal[va_idx]
        
        tmp_model = xgb.XGBClassifier(
            **best_params,
            gamma=0,
            eval_metric="logloss",
            random_state=RANDOM_STATE,
            n_jobs=threads,
            verbosity=0,
        )
        tmp_model.fit(X_tr, y_tr)
        preds = tmp_model.predict_proba(X_va)[:, 1]
        cv_auroc_scores.append(roc_auc_score(y_va, preds))

    mean_cv_auroc = float(np.mean(cv_auroc_scores))
    print(f"--- [{target_term}] -> Mean CV AUROC (Balanced Rest): {mean_cv_auroc:.4f}", flush=True)


    # --- STEP 3: EVALUATE ON 1:1 TEST (AUROC) ---
    print(f"--- [{target_term}] Step 3: Evaluating Master Model on 1:1 balanced Test set...", flush=True)
    bal_test_idx = balance_indices_1to1(y_test, seed=RANDOM_STATE)
    X_test_bal = X_test_df.iloc[bal_test_idx]
    y_test_bal = y_test[bal_test_idx]

    bal_test_preds = master_model.predict_proba(X_test_bal)[:, 1]
    test_bal_auroc = float(roc_auc_score(y_test_bal, bal_test_preds))
    print(f"--- [{target_term}] -> Balanced Test AUROC: {test_bal_auroc:.4f}", flush=True)
    
    bal_test_aucpr = float(average_precision_score(y_test_bal, bal_test_preds))
    print(f"--- [{target_term}] -> Balanced Test AUCPR: {bal_test_aucpr:.4f}", flush=True)

    # --- STEP 4: EVALUATE ON FULL TEST (AUCPR & LOGLOSS) ---
    print(f"--- [{target_term}] Step 4: Evaluating Master Model on full imbalanced Test set...", flush=True)
    full_test_preds = master_model.predict_proba(X_test_df)[:, 1]
    
    full_test_aucpr = float(average_precision_score(y_test, full_test_preds))
    full_test_logloss = float(log_loss(y_test, np.clip(full_test_preds, 1e-15, 1 - 1e-15)))
    
    # Calculate baseline counts and enrichment
    n_total_test = len(y_test)
    n_pos_test = int(np.sum(y_test))
    fraction_positives = n_pos_test / n_total_test if n_total_test > 0 else 0.0
    percentage_positives = fraction_positives * 100.0
    fold_enrichment = float(full_test_aucpr / fraction_positives) if fraction_positives > 0 else 0.0
    
    print(f"--- [{target_term}] -> Full Test AUCPR: {full_test_aucpr:.4f} (Fold Enrichment: {fold_enrichment:.2f}x)", flush=True)
    print(f"--- [{target_term}] -> Full Test LogLoss: {full_test_logloss:.4f}", flush=True)
    
    full_test_auroc = float(roc_auc_score(y_test, full_test_preds))
    print(f"--- [{target_term}] -> Full Test AUROC: {full_test_auroc:.4f}", flush=True)

    # --- BUILD PREDICTIONS DATAFRAME ---
    is_balanced = np.zeros(len(y_test), dtype=bool)
    is_balanced[bal_test_idx] = True
    
    preds_df = pd.DataFrame({
        "Gene_ID": test_ids,
        "True_Label": y_test,
        "Predicted_Probability": full_test_preds,
        "Is_Balanced_Subset": is_balanced
    })


    # --- SHAP VALUES COMPUTATION ---
    print(f"\n--- [{target_term}] Computing SHAP values on full Test set...", flush=True)
    explainer = shap.TreeExplainer(master_model)
    shap_values = explainer.shap_values(X_test_df)
    expected_value = explainer.expected_value
    print(f"--- [{target_term}] SHAP values computed successfully.", flush=True)


    # --- EXPORTS ---
    print(f"--- [{target_term}] Saving outputs to disk...", flush=True)
    
    # 1. Save Master Model
    joblib.dump(master_model, out_model)
    
    # 2. Save Metrics
    metrics_payload = {
        "term": target_term,
        "balanced_rest_mean_cv_auroc": mean_cv_auroc,
        "balanced_test_auroc": test_bal_auroc,
        "balanced_test_aucpr": bal_test_aucpr,
        "test_total_genes": n_total_test,
        "test_positive_genes": n_pos_test,
        "test_positive_percentage": percentage_positives,
        "full_test_auroc": full_test_auroc,
        "full_test_aucpr": full_test_aucpr,
        "aucpr_fold_enrichment": fold_enrichment,
        "full_test_logloss": full_test_logloss
    }
    with open(out_metrics, "w") as f:
        json.dump(metrics_payload, f, indent=2)
        
    # 3. Save Raw Predictions
    preds_df.to_parquet(out_predictions, compression='snappy')
    
    # 4. Save Test Features DataFrame
    X_test_df.to_parquet(out_x_test_features, compression='snappy')
    
    # 5. Save SHAP array
    np.save(out_shap_values, shap_values)
    np.save(snakemake.output.shap_expected_value, expected_value)

    print(f"=======================================================", flush=True)
    print(f"--- [COMPLETE] Evaluation and SHAP generation finished.", flush=True)
    print(f"=======================================================\n", flush=True)