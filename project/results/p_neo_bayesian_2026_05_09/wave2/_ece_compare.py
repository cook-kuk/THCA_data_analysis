"""ECE comparison: Wave1 Bayesian vs VQC vs Structure-LR vs QK-SVM vs GP-Q on each subset.
RF baseline ECE not available locally (RF preds were not saved); we report from Wave 1 log only.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WAVE1 = HERE.parent / "wave1"

EPS = 1e-9

def ece(y, p, n_bins=15):
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    valid = ~np.isnan(p)
    y, p = y[valid], p[valid]
    if len(y) == 0:
        return float("nan"), 0
    bins = np.linspace(0, 1, n_bins + 1)
    inds = np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)
    e = 0.0
    n = len(y)
    for b in range(n_bins):
        m = inds == b
        if m.sum() == 0:
            continue
        e += (m.sum() / n) * abs(y[m].mean() - p[m].mean())
    return float(e), int(len(y))

def brier(y, p):
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    valid = ~np.isnan(p)
    y, p = y[valid], p[valid]
    if len(y) == 0:
        return float("nan")
    return float(np.mean((p - y) ** 2))

# Load
w1 = pd.read_csv(WAVE1 / "predictions_itsndb.tsv", sep="\t")
vqc = pd.read_csv(HERE / "vqc_predictions.tsv", sep="\t")
struct = pd.read_csv(HERE / "structure_lr_predictions.tsv", sep="\t")
qk = pd.read_csv(HERE / "qk_svm_predictions.tsv", sep="\t")

key = ["peptide", "HLA_norm", "split"]
m = w1[key + ["label", "in_master", "pred_mean"]].rename(columns={"pred_mean": "p_w1"})
m = m.merge(vqc[key + ["p_vqc"]], on=key, how="outer")
m = m.merge(struct[key + ["p_struct_lr"]], on=key, how="outer")
m = m.merge(qk[key + ["p_qk_svm", "p_rbf_svm", "p_gp_rbf", "p_gp_q"]], on=key, how="outer")
m["label"] = m["label"].fillna(0).astype(int)
m["in_master"] = m["in_master"].fillna(False).astype(bool)

rows = []
for sname, mask in [
    ("ITSNdb_combined", np.ones(len(m), dtype=bool)),
    ("ITSNdb_no_overlap", (~m["in_master"]).to_numpy()),
    ("ITSNdb_in_master", m["in_master"].to_numpy()),
]:
    y = m.loc[mask, "label"].to_numpy()
    for col, name in [
        ("p_w1", "Wave1_Bayesian"),
        ("p_struct_lr", "Structure_LR"),
        ("p_vqc", "VQC_quantum"),
        ("p_qk_svm", "QK_SVM_quantum"),
        ("p_rbf_svm", "RBF_SVM_classical"),
        ("p_gp_rbf", "GP_RBF_classical_Bayesian"),
        ("p_gp_q", "GP_quantum_Bayesian"),
    ]:
        p = m.loc[mask, col].to_numpy()
        e, n = ece(y, p)
        b = brier(y, p)
        rows.append(dict(subset=sname, model=name, n_eval=n, ECE=e, Brier=b))

# RF baseline reported number from train.log
# (we only have user-provided 0.431 AUROC; for ECE we don't have it — note as missing)
res = pd.DataFrame(rows)
res.to_csv(HERE / "ece_comparison.tsv", sep="\t", index=False)
print(res.to_string(index=False))

best = {}
for s in ["ITSNdb_no_overlap", "ITSNdb_combined", "ITSNdb_in_master"]:
    sub = res[(res["subset"] == s) & res["ECE"].notna()]
    if len(sub):
        b = sub.loc[sub["ECE"].idxmin()]
        best[s] = dict(model=str(b["model"]), ECE=float(b["ECE"]), Brier=float(b["Brier"]))

(HERE / "ece_best_per_subset.json").write_text(json.dumps(best, indent=2))
print(json.dumps(best, indent=2))
