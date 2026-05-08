#!/usr/bin/env python
"""
Compare ITSNdb AUROC: (1) Bayesian-only baseline, (2) 3D-only LR, (3) Bayesian + 3D combo LR.
Combo LR trained on train-pool 3D features + a re-derived bayesian-like score?
We don't have predictions on the train pool here, so we approximate by:
- score_bayes_test = pred_mean from predictions_itsndb.tsv on the matching peptide/HLA
- combo on test side via simple stacking with prior LR coefs (use logreg on test split CV)
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

FEATS = ["mean_pLDDT_peptide", "min_pLDDT_peptide", "anchor_pLDDT", "mean_pLDDT_HLA",
         "interface_contacts_8A", "n_buried_residues",
         "radius_of_gyration_peptide", "peptide_helicity"]


def boot_auc(y, p, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    if len(y) < 5 or len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, len(y), size=len(y))
        ys = y[idx]; ps = p[idx]
        if len(np.unique(ys)) < 2: continue
        aucs.append(roc_auc_score(ys, ps))
    if not aucs: return float("nan"), float("nan"), float("nan")
    return float(np.mean(aucs)), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


feat = pd.read_csv("pmhc_3d_features.tsv", sep="\t")
preds = pd.read_csv("predictions_itsndb.tsv", sep="\t")
m = feat[feat["set"] == "itsndb"].merge(preds[["peptide", "HLA_norm", "pred_mean", "pred_std", "in_master"]],
                                        on=["peptide", "HLA_norm"], suffixes=("", "_pred"))
m = m.dropna(subset=FEATS + ["pred_mean", "label"]).copy()
m["label"] = m["label"].astype(int)
print(f"merged itsndb rows: {len(m)}")

rows = []
for tag, mask in [("itsndb_all", np.ones(len(m), dtype=bool)),
                  ("in_master_TRUE", (m["in_master"].astype(str) == "True").values),
                  ("in_master_FALSE_no_overlap", (m["in_master"].astype(str) == "False").values)]:
    sub = m.loc[mask].copy()
    if len(sub) < 10 or sub["label"].nunique() < 2:
        continue
    y = sub["label"].values
    p_bayes = sub["pred_mean"].values
    auc_b = roc_auc_score(y, p_bayes)
    a_b, lo_b, hi_b = boot_auc(y, p_bayes)
    rows.append({"stratum": tag, "model": "bayesian_only", "n": len(sub), "auroc": auc_b,
                 "boot_mean": a_b, "ci_lo": lo_b, "ci_hi": hi_b})

    # 3D-only via 5-fold CV inside this stratum
    Xf = sub[FEATS].values
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    p_3d = np.zeros(len(sub))
    for tr, te in skf.split(Xf, y):
        sc = StandardScaler().fit(Xf[tr])
        clf = LogisticRegression(max_iter=2000, C=1.0).fit(sc.transform(Xf[tr]), y[tr])
        p_3d[te] = clf.predict_proba(sc.transform(Xf[te]))[:, 1]
    auc_3 = roc_auc_score(y, p_3d)
    a_3, lo_3, hi_3 = boot_auc(y, p_3d)
    rows.append({"stratum": tag, "model": "3d_only_cv", "n": len(sub), "auroc": auc_3,
                 "boot_mean": a_3, "ci_lo": lo_3, "ci_hi": hi_3})

    # combo: bayesian + 3D, stacked via 5-fold
    Xc = np.concatenate([p_bayes.reshape(-1, 1), Xf], axis=1)
    p_co = np.zeros(len(sub))
    for tr, te in skf.split(Xc, y):
        sc = StandardScaler().fit(Xc[tr])
        clf = LogisticRegression(max_iter=2000, C=1.0).fit(sc.transform(Xc[tr]), y[tr])
        p_co[te] = clf.predict_proba(sc.transform(Xc[te]))[:, 1]
    auc_c = roc_auc_score(y, p_co)
    a_c, lo_c, hi_c = boot_auc(y, p_co)
    rows.append({"stratum": tag, "model": "bayes_plus_3d_cv", "n": len(sub), "auroc": auc_c,
                 "boot_mean": a_c, "ci_lo": lo_c, "ci_hi": hi_c})

res = pd.DataFrame(rows)
print(res.to_string(index=False))
res.to_csv("auroc_3d_combined.tsv", sep="\t", index=False)
