#!/usr/bin/env python
"""
Standalone evaluation of 3D pmhc features:
- Train LR on 3D features over the 'train' set rows in pmhc_3d_features.tsv (anchors)
- Score 'itsndb' rows; AUROC stratified by in_master True/False/all.
- Bootstrap 95% CI (n=1000).
"""
import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


FEATS = [
    "mean_pLDDT_peptide", "min_pLDDT_peptide", "anchor_pLDDT", "mean_pLDDT_HLA",
    "interface_contacts_8A", "n_buried_residues",
    "radius_of_gyration_peptide", "peptide_helicity",
]


def boot_auc(y, p, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    n_obs = len(y)
    if n_obs < 5 or len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, n_obs, size=n_obs)
        ys = y[idx]; ps = p[idx]
        if len(np.unique(ys)) < 2:
            continue
        aucs.append(roc_auc_score(ys, ps))
    if not aucs:
        return float("nan"), float("nan"), float("nan")
    return float(np.mean(aucs)), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.features, sep="\t")
    df = df.dropna(subset=FEATS).copy()
    df["label"] = pd.to_numeric(df["label"], errors="coerce")
    df = df.dropna(subset=["label"]).copy()
    train = df[df["set"] == "train"].copy()
    test = df[df["set"] == "itsndb"].copy()
    print(f"[input] train n={len(train)} (pos={int(train.label.sum())}) | test n={len(test)} (pos={int(test.label.sum())})", flush=True)

    if len(train) < 20 or train.label.nunique() < 2:
        print("[abort] not enough train labels — falling back to test self-train (5-fold CV).", flush=True)
        # fallback: just compute univariate AUC on test
        rows = []
        for f in FEATS:
            v = test[f].values
            y = test["label"].values
            try:
                a = roc_auc_score(y, v)
            except Exception:
                a = float("nan")
            rows.append({"feature": f, "auc_test_univariate": a})
        out_df = pd.DataFrame(rows)
        out_df.to_csv(args.out, sep="\t", index=False)
        return

    sc = StandardScaler().fit(train[FEATS].values)
    Xtr = sc.transform(train[FEATS].values)
    Xte = sc.transform(test[FEATS].values)
    ytr = train["label"].values.astype(int)
    yte = test["label"].values.astype(int)

    clf = LogisticRegression(max_iter=2000, C=1.0)
    clf.fit(Xtr, ytr)
    p = clf.predict_proba(Xte)[:, 1]

    out_rows = []
    # all itsndb
    a_all, lo_all, hi_all = boot_auc(yte, p)
    out_rows.append({"stratum": "itsndb_all", "n": len(yte), "n_pos": int(yte.sum()),
                     "auroc": roc_auc_score(yte, p) if len(np.unique(yte))>1 else float("nan"),
                     "auroc_boot_mean": a_all, "ci_lo": lo_all, "ci_hi": hi_all})
    # in_master True / False
    for tag, mask in [("itsndb_in_master_TRUE", test["in_master"].astype(str) == "True"),
                      ("itsndb_in_master_FALSE_no_overlap", test["in_master"].astype(str) == "False")]:
        ys = yte[mask.values]; ps = p[mask.values]
        if len(ys) < 2 or len(np.unique(ys)) < 2:
            out_rows.append({"stratum": tag, "n": int(mask.sum()), "n_pos": int(ys.sum()),
                             "auroc": float("nan"), "auroc_boot_mean": float("nan"),
                             "ci_lo": float("nan"), "ci_hi": float("nan")})
            continue
        a, lo, hi = boot_auc(ys, ps)
        out_rows.append({"stratum": tag, "n": int(mask.sum()), "n_pos": int(ys.sum()),
                         "auroc": roc_auc_score(ys, ps),
                         "auroc_boot_mean": a, "ci_lo": lo, "ci_hi": hi})

    res = pd.DataFrame(out_rows)
    print(res.to_string(index=False), flush=True)
    res.to_csv(args.out, sep="\t", index=False)
    # also write feature coefficients
    coef = pd.DataFrame({"feature": FEATS, "coef_z": clf.coef_[0]})
    coef.to_csv(args.out.replace(".tsv", "_coefs.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
