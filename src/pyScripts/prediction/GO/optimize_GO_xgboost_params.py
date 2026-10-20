import datetime
import gc
import json
import os
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, log_loss
from sklearn.model_selection import train_test_split
from skopt import Optimizer
from skopt.space import Integer, Real
from xgboost import XGBClassifier

# --- SETTINGS ---
RANDOM_STATE = 1234
TEST_SIZE = 0.20          # 20% of full data held out as test set (80% goes to rest)
VAL_SIZE = 0.20           # Inside tuning: 20% of rest used for validation (80% for train)
N_ITERS = 50              # Number of skopt Bayesian optimization iterations


def balance_dataset_1to1(X, y, seed):
    """
    Subsamples the non-annotated class (0) to match the number of annotated class (1).
    Returns a shuffled 1:1 balanced dataset.
    """
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    
    n_pos = len(pos_idx)
    
    # If there are no positive genes in this split, return as is to avoid crashing
    if n_pos == 0:
        return X, y
        
    rng = np.random.default_rng(seed)
    
    # Randomly sample negative indices to match the number of positives
    neg_sampled_idx = rng.choice(neg_idx, size=n_pos, replace=False)
    
    # Combine and shuffle the indices
    balanced_idx = np.concatenate([pos_idx, neg_sampled_idx])
    rng.shuffle(balanced_idx)
    
    return X[balanced_idx], y[balanced_idx]


def hyperparameter_tuned_model(X_rest, y_rest, n_threads, fileout, target_term, n_iters=N_ITERS):
    print(f"--- [{target_term}] Hyperparameter tuning started", flush=True)

    # Define the hyperparameter space
    param_dist = [
        Integer(3, 10, name="max_depth"),
        Real(0.01, 0.3, prior="log-uniform", name="learning_rate"),
        Integer(1, 50, name="min_child_weight"),
        Real(0.3, 1.0, name="colsample_bytree"),
        Real(0.5, 1.0, name="subsample"),
        Integer(50, 400, name="n_estimators"),
    ]

    best_score = np.inf
    best_model = None
    best_params = None

    hyper_optimizer = Optimizer(
        dimensions=param_dist,
        base_estimator="RF",
        acq_func="EI",
        random_state=RANDOM_STATE,
    )

    for i in range(n_iters):
        s = datetime.datetime.now()

        # Fresh stratified train/validation split of the rest set
        X_train_full, X_validation_full, y_train_full, y_validation_full = train_test_split(
            X_rest,
            y_rest,
            test_size=VAL_SIZE,
            stratify=y_rest,
            random_state=RANDOM_STATE + i,
        )
        
        # Subsample the non-annotated to equal proportions (1:1)
        X_train, y_train = balance_dataset_1to1(X_train_full, y_train_full, seed=RANDOM_STATE + i)
        X_validation, y_validation = balance_dataset_1to1(X_validation_full, y_validation_full, seed=RANDOM_STATE + i + 1000)

        # "Ask" the hyperparameter optimizer for parameters
        params = hyper_optimizer.ask()
        params_dict = dict(zip([d.name for d in param_dist], params))

        # Train model on train and evaluate logloss performance on validation
        model = XGBClassifier(
            **params_dict,
            gamma=0,
            eval_metric="logloss",
            random_state=RANDOM_STATE,
            n_jobs=n_threads,
            verbosity=0,
        )
        model.fit(X_train, y_train)

        probs = model.predict_proba(X_validation)[:, 1]
        probs = np.clip(probs, 1e-15, 1 - 1e-15)
        c_log_loss = log_loss(y_validation, probs)
        val_acc = balanced_accuracy_score(y_validation, (probs > 0.5).astype(int))

        # Provide parameters:logloss to the hyperparameter optimizer
        hyper_optimizer.tell(params, c_log_loss)

        e = datetime.datetime.now()
        print("------------------------------------------------")
        print(f"--- [{target_term}] {i + 1} iteration of {n_iters} in {e - s}")
        print(f"--- [{target_term}] Current params: " + str(params_dict))
        print(f"--- [{target_term}] log loss: {c_log_loss:.5f}\tAcc: {val_acc:.5f} for Validation (1:1 Balanced)", flush=True)

        fileout.write("---------------------\n")
        fileout.write(f"Iteration {i + 1}/{n_iters} | Training took {e - s} using {n_threads} threads\n")
        fileout.write(f"Current params: " + str(params_dict) + "\n")
        fileout.write(f"Validation log_loss: {c_log_loss}\tBalanced_Acc: {val_acc}\n")
        fileout.flush()

        # If current parameters perform better, save parameters as the best ones
        if c_log_loss < best_score:
            best_score = c_log_loss
            best_model = model
            best_params = params

    # Convert best_params values to native Python types for clean serialization
    best_params_dict = {
        d.name: (int(val) if isinstance(d, Integer) else float(val))
        for d, val in zip(param_dist, best_params)
    }

    fileout.write("=====================\n")
    fileout.write("Best validation score: " + str(best_score) + "\n")
    fileout.write("Best params: " + str(best_params_dict) + "\n")
    fileout.flush()

    return best_model, float(best_score), best_params_dict


if __name__ == "__main__":
    matrix_file = snakemake.input.full_matrix
    y_file = snakemake.input.y_file

    out_model = snakemake.output.saved_model
    out_params_json = snakemake.output.best_params
    out_log = snakemake.output.tuning_log
    out_rest_ids = snakemake.output.rest_ids
    out_test_ids = snakemake.output.test_ids

    threads = snakemake.threads
    target_term = snakemake.wildcards.term.replace("_", ":")

    for path in (out_model, out_params_json, out_log, out_rest_ids, out_test_ids):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    print("\n=======================================================", flush=True)
    print(f"--- [START] Tuning & Model Fitting for {target_term}", flush=True)
    print("=======================================================\n", flush=True)

    # --- Take the full gene x feature matrix -> split 80-20% keeping proportions ---
    print(f"--- [{target_term}] Reading labels and performing 80/20 stratified rest/test split...", flush=True)
    y_df = pd.read_parquet(y_file)
    y_df["Gene_ID"] = y_df["Gene_ID"].astype(str)
    y_df = y_df.sort_values("Gene_ID").reset_index(drop=True)

    rest_df, test_df = train_test_split(
        y_df,
        test_size=TEST_SIZE,
        stratify=y_df["Label"],
        random_state=RANDOM_STATE,
    )

    # Save gene IDs for rest and test sets
    rest_df[["Gene_ID"]].to_csv(out_rest_ids, index=False)
    test_df[["Gene_ID"]].to_csv(out_test_ids, index=False)

    print(f"--- [{target_term}] -> Rest set: {len(rest_df)} genes | Positives: {rest_df['Label'].sum()} ({rest_df['Label'].mean():.4f})", flush=True)
    print(f"--- [{target_term}] -> Test set: {len(test_df)} genes | Positives: {test_df['Label'].sum()} ({test_df['Label'].mean():.4f})", flush=True)

    rest_ids = rest_df["Gene_ID"].tolist()
    y_rest = rest_df["Label"].to_numpy(dtype=np.int8)
    del y_df, rest_df, test_df
    gc.collect()

    # --- LOAD FULL CLOSENESS MATRIX (REMOVE TEST ROWS, KEEP ALL COLUMNS & DIAGONAL) ---
    print(f"--- [{target_term}] Loading closeness matrix (keeping all columns and diagonal intact)...", flush=True)
    X_df = pd.read_parquet(matrix_file)
    X_df["Gene_ID"] = X_df["Gene_ID"].astype(str)
    X_df = X_df.set_index("Gene_ID")

    X_rest = X_df.loc[rest_ids].to_numpy(dtype=np.float32, copy=True)
    del X_df
    gc.collect()

    print(f"--- [{target_term}] -> X_rest shape: {X_rest.shape}", flush=True)

    # --- HYPERPARAMETER TUNING ---
    with open(out_log, "w") as log_file:
        _, best_val_score, best_parameters = hyperparameter_tuned_model(
            X_rest,
            y_rest,
            n_threads=threads,
            fileout=log_file,
            target_term=target_term,
            n_iters=N_ITERS,
        )

    # Save best parameters as JSON for easy inspection or downstream use
    best_payload = {
        "term": target_term,
        "best_validation_logloss": best_val_score,
        "n_iters": N_ITERS,
        "random_state": RANDOM_STATE,
        "params": best_parameters,
    }
    with open(out_params_json, "w") as f:
        json.dump(best_payload, f, indent=2)

    # --- FINAL MODEL FIT (Just in case) ---
    print(f"\n--- [{target_term}] Balancing the final Rest set 1:1 and fitting final model...", flush=True)
    X_rest_balanced, y_rest_balanced = balance_dataset_1to1(X_rest, y_rest, seed=RANDOM_STATE)
    
    final_xgb = XGBClassifier(
        **best_parameters,
        gamma=0,
        eval_metric="logloss",
        random_state=RANDOM_STATE,
        n_jobs=threads,
        verbosity=0,
    )
    final_xgb.fit(X_rest_balanced, y_rest_balanced)

    joblib.dump(final_xgb, out_model)
    print(f"--- [{target_term}] Final model saved to: {out_model}", flush=True)
    print("=======================================================\n", flush=True)