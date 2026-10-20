import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
import mygene

if __name__ == "__main__":
    shap_file        = snakemake.input.shap_values
    expected_file    = snakemake.input.shap_expected_value
    features_file    = snakemake.input.x_test_features
    preds_file       = snakemake.input.test_predictions

    out_csv          = snakemake.output.novel_candidates_csv
    out_plot_dir     = snakemake.output.waterfall_dir

    target_term = snakemake.wildcards.term.replace("_", ":")

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    os.makedirs(out_plot_dir, exist_ok=True)

    print(f"\n=======================================================", flush=True)
    print(f"--- [START] Novel Candidate Extraction for {target_term}", flush=True)
    print(f"=======================================================\n", flush=True)

    # --- LOAD DATA ---
    print(f"--- [{target_term}] Loading predictions, SHAP values, expected value, and features...", flush=True)
    shap_values = np.load(shap_file)
    expected_value = np.load(expected_file)
    
    # Ensure expected_value is a clean float (depending on numpy save format)
    if expected_value.ndim == 0:
        expected_value = float(expected_value)
    elif expected_value.size == 1:
        expected_value = float(expected_value[0])
        
    X_test_features = pd.read_parquet(features_file)
    preds_df = pd.read_parquet(preds_file)
    entrez_ids = X_test_features.columns.tolist()

    # --- TRANSLATE NETWORK NODES (MYGENE) ---
    print(f"--- [{target_term}] Translating {len(entrez_ids)} Entrez IDs to Gene Symbols...", flush=True)
    mg = mygene.MyGeneInfo()
    results = mg.querymany(entrez_ids, scopes='entrezgene', fields='symbol', species='human', verbose=False)
    
    rename_dict = {}
    for res in results:
        query_id = str(res['query'])
        symbol = res.get('symbol', query_id)
        rename_dict[query_id] = symbol
        
    X_test_features.rename(columns=rename_dict, inplace=True)
    feature_names = X_test_features.columns.tolist()

    # Translate the test Gene_IDs themselves for the CSV
    test_gene_ids = preds_df["Gene_ID"].astype(str).tolist()
    results_test = mg.querymany(test_gene_ids, scopes='entrezgene', fields='symbol', species='human', verbose=False)
    test_rename_dict = {str(r['query']): r.get('symbol', str(r['query'])) for r in results_test}
    preds_df["Gene_Symbol"] = preds_df["Gene_ID"].astype(str).map(test_rename_dict)
    
    del mg

    # --- FILTER FOR NOVEL CANDIDATES (FALSE POSITIVES) ---
    print(f"--- [{target_term}] Isolating unannotated genes with high probabilities...", flush=True)
    
    # Keep only unannotated genes (Label == 0)
    candidates_df = preds_df[preds_df["True_Label"] == 0].copy()
    
    # Preserve original row index so we can match it back to the shap_values matrix
    candidates_df['original_index'] = candidates_df.index
    
    # Sort by probability descending
    candidates_df = candidates_df.sort_values(by="Predicted_Probability", ascending=False).reset_index(drop=True)
    original_indices = candidates_df["original_index"].values

    # --- EXTRACT TOP 3 DRIVERS FOR CSV ---
    top_driver_1, top_driver_2, top_driver_3 = [], [], []

    for row_idx in original_indices:
        local_shap = shap_values[row_idx]
        # Get indices of the top 3 highest POSITIVE SHAP values (features pushing the probability up)
        top_3_idx = np.argsort(local_shap)[-3:][::-1]
        
        top_driver_1.append(f"{feature_names[top_3_idx[0]]} (+{local_shap[top_3_idx[0]]:.3f})")
        top_driver_2.append(f"{feature_names[top_3_idx[1]]} (+{local_shap[top_3_idx[1]]:.3f})")
        top_driver_3.append(f"{feature_names[top_3_idx[2]]} (+{local_shap[top_3_idx[2]]:.3f})")

    candidates_df["Top_Driver_1"] = top_driver_1
    candidates_df["Top_Driver_2"] = top_driver_2
    candidates_df["Top_Driver_3"] = top_driver_3

    # Clean up and export the top 50 candidates
    export_df = candidates_df[["Gene_ID", "Gene_Symbol", "Predicted_Probability", 
                               "Top_Driver_1", "Top_Driver_2", "Top_Driver_3"]].head(50)
    export_df.to_csv(out_csv, index=False)
    print(f"--- [{target_term}] Saved top 50 novel candidates to {out_csv}")

    # --- PLOT LOCAL IMPACT FOR TOP 3 CANDIDATES (WATERFALL) ---
    print(f"--- [{target_term}] Generating local SHAP waterfall plots for Top 3 candidates...", flush=True)
    
    for rank in range(min(3, len(candidates_df))):
        row_idx = candidates_df.loc[rank, "original_index"]
        target_symbol = candidates_df.loc[rank, "Gene_Symbol"]
        prob = candidates_df.loc[rank, "Predicted_Probability"]

        # Create the Explanation object required by shap.plots.waterfall
        explanation = shap.Explanation(
            values=shap_values[row_idx], 
            base_values=expected_value, 
            data=X_test_features.iloc[row_idx].values,
            feature_names=feature_names
        )

        plt.figure(figsize=(10, 6))
        # Generate the native waterfall plot (max_display=10 keeps it readable)
        shap.plots.waterfall(explanation, max_display=10, show=False)
        
        plt.title(f"Why was {target_symbol} flagged?\nPredicted Prob: {prob:.2f} | True Label: 0", pad=20)
        plt.tight_layout()

        clean_symbol = str(target_symbol).replace("/", "_")
        plot_path = os.path.join(out_plot_dir, f"Rank_{rank+1}_{clean_symbol}_waterfall.pdf")
        plt.savefig(plot_path, format='pdf', bbox_inches='tight')
        plt.close()

    print(f"--- [{target_term}] Local waterfall plots saved to {out_plot_dir}")
    print(f"=======================================================", flush=True)
    print(f"--- [COMPLETE] Candidate Extraction finished.", flush=True)
    print(f"=======================================================\n", flush=True)