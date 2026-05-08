"""Post-hoc eval of VQC predictions saved by train_vqc.py (without re-training).
Handles NaN merge correctly.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
WAVE1 = HERE.parent / "wave1"

vqc = pd.read_csv(HERE / "vqc_predictions.tsv", sep="\t")
print(f"[load] vqc preds: {len(vqc)}, cols={list(vqc.columns)}")

itsn_pred = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
m = vqc.merge(itsn_pred[["peptide", "HLA_norm", "split", "pred_mean", "pred_std"]],
              on=["peptide", "HLA_norm", "split"], how="left")
print(f"[merge] {len(m)} rows; n with wave1 pred: {(~m['pred_mean'].isna()).sum()}")

def auc_or_nan(y, p):
    p = np.asarray(p, dtype=float)
    y = np.asarray(y, dtype=int)
    valid = ~np.isnan(p)
    if valid.sum() == 0 or len(np.unique(y[valid])) < 2:
        return float("nan"), int(valid.sum())
    return float(roc_auc_score(y[valid], p[valid])), int(valid.sum())

rows = []
in_master = m["in_master"].astype(bool).to_numpy()
for subset, mask in [
    ("ITSNdb_combined", np.ones(len(m), dtype=bool)),
    ("ITSNdb_no_overlap", ~in_master),
    ("ITSNdb_in_master", in_master),
]:
    yy = m.loc[mask, "label"].to_numpy()
    for col, name in [("p_vqc", "VQC_8q3l"),
                      ("p_lr_8d_classical", "LR_8d_PCA_classical"),
                      ("pred_mean", "Wave1_Bayesian")]:
        a, nv = auc_or_nan(yy, m.loc[mask, col].to_numpy())
        rows.append(dict(subset=subset, model=name, n_total=int(mask.sum()), n_valid=nv, AUROC=a))

res = pd.DataFrame(rows)
res.to_csv(HERE / "vqc_results.tsv", sep="\t", index=False)
print(res.to_string(index=False))

# Figure
fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
order = ["ITSNdb_combined", "ITSNdb_no_overlap", "ITSNdb_in_master"]
models = ["VQC_8q3l", "LR_8d_PCA_classical", "Wave1_Bayesian"]
colors = ["#9467bd", "#1f77b4", "#2ca02c"]
W = 0.27
for i, mod in enumerate(models):
    sub = res[res["model"] == mod].set_index("subset").reindex(order)
    ax.bar(np.arange(len(order)) + (i - 1) * W, sub["AUROC"], width=W, label=mod, color=colors[i])
ax.axhline(0.5, ls="--", color="grey", lw=1)
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels(order, rotation=15, ha="right")
ax.set_ylabel("AUROC")
ax.set_ylim(0.0, 1.0)
ax.set_title("Wave 2 Track C — VQC (8q,3L) vs classical LR(8d-PCA) vs Wave 1 Bayesian")
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3, axis="y")
fig.savefig(HERE / "fig_quantum_vs_classical.png", dpi=160)
fig.savefig(HERE / "fig_quantum_vs_classical.pdf")
plt.close(fig)
print("[done] VQC eval")
