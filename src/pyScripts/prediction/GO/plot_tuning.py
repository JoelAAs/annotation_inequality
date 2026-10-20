import ast
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

log_path = snakemake.input.tuning_log
out_png = snakemake.output.plot
out_csv = snakemake.output.history
n_initial = snakemake.params.n_initial

os.makedirs(os.path.dirname(out_png), exist_ok=True)
os.makedirs(os.path.dirname(out_csv), exist_ok=True)

# --- Parse the tuning log ---
rows, current = [], {}
with open(log_path) as f:
    for line in f:
        line = line.strip()
        m = re.match(r"Iteration (\d+)/(\d+)", line)
        if m:
            current = {"iter": int(m.group(1))}
        elif line.startswith("Current params:"):
            txt = line.split(":", 1)[1].strip()
            txt = re.sub(r"np\.\w+\(([^)]*)\)", r"\1", txt)  # np.int64(5) -> 5
            current.update(ast.literal_eval(txt))
        elif line.startswith("Validation log_loss:"):
            m = re.match(r"Validation log_loss: ([\d.eE+-]+)\s+Balanced_Acc: ([\d.eE+-]+)", line)
            current["log_loss"] = float(m.group(1))
            current["bal_acc"] = float(m.group(2))
            rows.append(current)

df = pd.DataFrame(rows)
df["best_so_far"] = df["log_loss"].cummin()
df.to_csv(out_csv, sep="\t", index=False)
best = df.loc[df["log_loss"].idxmin()]

# Search-space bounds, same as in the tuning script
bounds = {
    "learning_rate": (0.01, 0.3),
    "max_depth": (3, 10),
}

star = dict(marker="*", s=250, color="gold", edgecolor="black", zorder=5)

# --- Figure: 2x2 grid, constrained layout prevents overlaps ---
fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")

# Convergence (top left)
ax = axes[0, 0]
ax.plot(df["iter"], df["log_loss"], "o", alpha=0.6, label="each try")
ax.plot(df["iter"], df["best_so_far"], "-", color="red", label="best so far")
ax.scatter(best["iter"], best["log_loss"], label="best", **star)
ax.axvline(n_initial + 0.5, ls="--", color="gray")
ax.set(xlabel="iteration", ylabel="validation log loss",
       title="Convergence (dashed line = start of guided search)")
ax.legend()

# Balanced accuracy (top right)
ax = axes[0, 1]
ax.plot(df["iter"], df["bal_acc"], "o-", alpha=0.7)
ax.scatter(best["iter"], best["bal_acc"], **star)
ax.axvline(n_initial + 0.5, ls="--", color="gray")
ax.set(xlabel="iteration", ylabel="balanced accuracy",
       title="Balanced accuracy (threshold 0.5)")

# Most influential hyperparameters vs log loss (bottom row)
for ax, p in zip(axes[1, :], bounds):
    sc = ax.scatter(df[p], df["log_loss"], c=df["iter"], cmap="viridis")
    ax.scatter(best[p], best["log_loss"], **star)
    for b in bounds[p]:
        ax.axvline(b, ls=":", color="gray")
    ax.set(xlabel=p, ylabel="validation log loss", title=f"{p} vs log loss")
    if p == "learning_rate":
        ax.set_xscale("log")

# One shared colorbar for the bottom row; constrained layout reserves space for it
fig.colorbar(sc, ax=axes[1, :], label="iteration", shrink=0.9)

fig.suptitle(
    f"{snakemake.wildcards.term.replace('_', ':')}  |  best log loss {best['log_loss']:.4f} "
    f"at iter {int(best['iter'])}"
)
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close(fig)