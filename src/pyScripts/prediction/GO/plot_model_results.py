import os
import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt
import pronto
import mygene
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

if __name__ == "__main__":
    # --- SNAKEMAKE I/O ---
    shap_file        = snakemake.input.shap_values
    features_file    = snakemake.input.x_test_features
    preds_file       = snakemake.input.test_predictions
    ontology_file    = snakemake.input.ontology

    out_shap_plot    = snakemake.output.shap_summary_plot
    out_curves_plot  = snakemake.output.performance_curves
    out_csv          = snakemake.output.shap_top_genes

    # Convert wildcard back to standard GO ID
    target_term = snakemake.wildcards.term.replace("_", ":")

    for path in [out_shap_plot, out_curves_plot, out_csv]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    print(f"\n=======================================================", flush=True)
    print(f"--- [START] Global Results & Visualization for {target_term}", flush=True)
    print(f"=======================================================\n", flush=True)

    # --- LOAD DATA ---
    print(f"--- [{target_term}] Loading predictions, SHAP values, and test features...", flush=True)
    shap_values = np.load(shap_file)
    X_test_features = pd.read_parquet(features_file)
    entrez_ids = X_test_features.columns.tolist()
    
    preds_df = pd.read_parquet(preds_file)
    
    # Full dataset arrays
    y_true = preds_df["True_Label"].values
    y_probs = preds_df["Predicted_Probability"].values
    
    # 1:1 Balanced subset arrays
    bal_df = preds_df[preds_df["Is_Balanced_Subset"] == True]
    y_true_bal = bal_df["True_Label"].values
    y_probs_bal = bal_df["Predicted_Probability"].values

    # --- GET GO TERM NAME (PRONTO) ---
    print(f"--- [{target_term}] Fetching term name from ontology...", flush=True)
    try:
        onto = pronto.Ontology(ontology_file)
        term_name = onto[target_term].name
        base_title = f"{term_name} ({target_term})"
    except Exception as e:
        print(f"--- [{target_term}] Warning: Could not fetch term name ({e}). Defaulting to ID.", flush=True)
        base_title = f"{target_term}"

    # --- GENERATE 2x3 PERFORMANCE DASHBOARD ---
    print(f"--- [{target_term}] Generating ROC, PRC, and Probability Histograms...", flush=True)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(f"Model Performance: {base_title}", fontsize=18, y=1.02)
    
    # Row 0: Full Test Set
    ax_roc_full, ax_prc_full, ax_hist_full = axes[0]
    # Row 1: Balanced Test Set
    ax_roc_bal, ax_prc_bal, ax_hist_bal = axes[1]

    # ==========================================
    # ROW 1: FULL IMBALANCED TEST SET
    # ==========================================
    # A1. ROC Curve (Full)
    fpr, tpr, _ = roc_curve(y_true, y_probs)
    roc_auc = auc(fpr, tpr)
    ax_roc_full.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC (AUC = {roc_auc:.3f})')
    ax_roc_full.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax_roc_full.set(xlim=[0.0, 1.0], ylim=[0.0, 1.05], xlabel='False Positive Rate', ylabel='True Positive Rate', title='ROC (Full Test)')
    ax_roc_full.legend(loc="lower right")
    ax_roc_full.grid(True, alpha=0.3)

    # B1. Precision-Recall Curve (Full)
    precision, recall, _ = precision_recall_curve(y_true, y_probs)
    pr_auc = average_precision_score(y_true, y_probs)
    baseline = np.sum(y_true) / len(y_true)
    ax_prc_full.plot(recall, precision, color='forestgreen', lw=2, label=f'PRC (AUC = {pr_auc:.3f})')
    ax_prc_full.plot([0, 1], [baseline, baseline], color='firebrick', lw=2, linestyle='--', label=f'Baseline ({baseline:.3f})')
    ax_prc_full.set(xlim=[0.0, 1.0], ylim=[0.0, 1.05], xlabel='Recall', ylabel='Precision', title='PRC (Full Test)')
    ax_prc_full.legend(loc="upper right")
    ax_prc_full.grid(True, alpha=0.3)

    # C1. Histogram (Full)
    ax_hist_full.hist(y_probs[y_true == 0], bins=50, alpha=0.6, color='firebrick', label='Unannotated (0)', log=True)
    ax_hist_full.hist(y_probs[y_true == 1], bins=50, alpha=0.6, color='forestgreen', label='Annotated (1)', log=True)
    ax_hist_full.set(xlim=[0.0, 1.0], xlabel='Predicted Probability', ylabel='Count (Log)', title='Probabilities (Full Test)')
    ax_hist_full.legend(loc="upper center")
    ax_hist_full.grid(True, alpha=0.3)

    # ==========================================
    # ROW 2: 1:1 BALANCED TEST SET
    # ==========================================
    # A2. ROC Curve (Balanced)
    fpr_bal, tpr_bal, _ = roc_curve(y_true_bal, y_probs_bal)
    roc_auc_bal = auc(fpr_bal, tpr_bal)
    ax_roc_bal.plot(fpr_bal, tpr_bal, color='darkorange', lw=2, label=f'ROC (AUC = {roc_auc_bal:.3f})')
    ax_roc_bal.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax_roc_bal.set(xlim=[0.0, 1.0], ylim=[0.0, 1.05], xlabel='False Positive Rate', ylabel='True Positive Rate', title='ROC (Balanced Test)')
    ax_roc_bal.legend(loc="lower right")
    ax_roc_bal.grid(True, alpha=0.3)

    # B2. Precision-Recall Curve (Balanced)
    precision_bal, recall_bal, _ = precision_recall_curve(y_true_bal, y_probs_bal)
    pr_auc_bal = average_precision_score(y_true_bal, y_probs_bal)
    baseline_bal = 0.50 # By definition, the balanced set is 50% positive
    ax_prc_bal.plot(recall_bal, precision_bal, color='forestgreen', lw=2, label=f'PRC (AUC = {pr_auc_bal:.3f})')
    ax_prc_bal.plot([0, 1], [baseline_bal, baseline_bal], color='firebrick', lw=2, linestyle='--', label=f'Baseline ({baseline_bal:.2f})')
    ax_prc_bal.set(xlim=[0.0, 1.0], ylim=[0.0, 1.05], xlabel='Recall', ylabel='Precision', title='PRC (Balanced Test)')
    ax_prc_bal.legend(loc="upper right")
    ax_prc_bal.grid(True, alpha=0.3)

    # C2. Histogram (Balanced)
    # No log scale needed here because classes are perfectly equal
    ax_hist_bal.hist(y_probs_bal[y_true_bal == 0], bins=25, alpha=0.6, color='firebrick', label='Unannotated (0)')
    ax_hist_bal.hist(y_probs_bal[y_true_bal == 1], bins=25, alpha=0.6, color='forestgreen', label='Annotated (1)')
    ax_hist_bal.set(xlim=[0.0, 1.0], xlabel='Predicted Probability', ylabel='Count', title='Probabilities (Balanced Test)')
    ax_hist_bal.legend(loc="upper center")
    ax_hist_bal.grid(True, alpha=0.3)

    # Save Dashboard
    plt.tight_layout()
    plt.savefig(out_curves_plot, format='pdf', bbox_inches='tight')
    plt.close()

    # --- TRANSLATE ENTREZ IDS TO GENE SYMBOLS (MYGENE) ---
    print(f"--- [{target_term}] Translating {len(entrez_ids)} Entrez IDs to Gene Symbols...", flush=True)
    mg = mygene.MyGeneInfo()
    results = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species='human', verbose=False)
    
    rename_dict = {}
    for res in results:
        query_id = str(res['query'])
        # If mygene doesn't find a symbol, default back to the Entrez ID
        symbol = res.get('symbol', query_id)
        rename_dict[query_id] = symbol
        
    X_test_features.rename(columns=rename_dict, inplace=True)
    
    # Safely close the mygene connection to prevent Python teardown errors in logs
    del mg

    # --- GENERATE SHAP SUMMARY PLOT ---
    print(f"--- [{target_term}] Generating SHAP summary plot...", flush=True)
    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_test_features, show=False)
    plt.title(f"SHAP Values: {base_title}", fontsize=14, pad=20)
    plt.tight_layout()
    
    plt.savefig(out_shap_plot, format='pdf', bbox_inches='tight')
    plt.close()

    # --- EXPORT RANKED TOP GENES ---
    print(f"--- [{target_term}] Exporting ranked driving genes to CSV...", flush=True)
    # SHAP importance is the mean absolute impact across all test predictions
    mean_abs_shap = np.abs(shap_values).mean(axis=0)
    
    top_genes_df = pd.DataFrame({
        'Entrez_ID': entrez_ids,
        'Gene_Symbol': [rename_dict.get(str(eid), str(eid)) for eid in entrez_ids],
        'Mean_Abs_SHAP': mean_abs_shap
    })
    
    # Sort from highest impact to lowest
    top_genes_df = top_genes_df.sort_values(by='Mean_Abs_SHAP', ascending=False).reset_index(drop=True)
    top_genes_df.to_csv(out_csv, index=False)

    print(f"--- [{target_term}] Plots and tables saved successfully.", flush=True)
    print(f"=======================================================", flush=True)
    print(f"--- [COMPLETE] Visualization finished.", flush=True)
    print(f"=======================================================\n", flush=True)