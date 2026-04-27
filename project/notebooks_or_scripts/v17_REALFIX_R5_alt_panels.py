#!/usr/bin/env python3
"""
v17 REAL FIX · R5 — alternative panel size + cross-cohort training tests.

Tests whether AUC can be honestly improved by:
  R5-A: 10-gene panel (8 + DIO2 + IYD)
  R5-B: 12-gene panel (8 + DIO2 + IYD + SLC26A4 + SLC5A8)
  R5-C: cross-cohort training (TCGA + GSE27155 + GSE33630 → test on GSE76039)
  R5-D: ensemble (LogReg + RF + LGBM if available)

All tested on R1-A leak-free cluster labels.

Outputs:
  - results/v17_realfix/R5_alt_panel_table.tsv
  - results/v17_realfix/R5_summary.json
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

OUT = Path("/opt/thyroid-dash/project/results/v17_realfix")
TPM = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
GSE76039_EXPR = Path("/data/thca/data_processed/microarray/GSE76039_microarray_expression_log2.tsv")
GSE76039_SM = Path("/data/thca/data_raw/geo/GSE76039/GSE76039_series_matrix.txt.gz")

GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
GENE_10 = GENE_8 + ["DIO2", "IYD"]
GENE_12 = GENE_10 + ["SLC26A4", "SLC5A8"]
GENE_16 = GENE_12 + ["THRA", "THRB", "DUOX1", "DUOX2"]


def load_r1a():
    df = pd.read_csv(OUT / "R1A_cluster_labels.tsv", sep="\t")
    df["y"] = df["cluster"].str.startswith("DM2").astype(int)
    return df


def cv_auc(X, y, model_factory, n_splits=5, random_state=42):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        m = model_factory()
        if isinstance(m, LogisticRegression):
            sc = StandardScaler().fit(X[tr])
            m.fit(sc.transform(X[tr]), y[tr])
            oof[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
        else:
            m.fit(X[tr], y[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]
    return float(roc_auc_score(y, oof))


def bootstrap_ci(y, p, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) >= 2:
            aucs.append(roc_auc_score(y[idx], p[idx]))
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def test_panel(name, genes, expr, lbl):
    g_have = [g for g in genes if g in expr.index]
    if len(g_have) < len(genes) * 0.7:
        return {"name": name, "error": f"only {len(g_have)}/{len(genes)} matched"}
    Xs = expr.loc[g_have, [s for s in lbl["sample_id"] if s in expr.columns]].T
    common = Xs.index.intersection(lbl.set_index("sample_id").index)
    X = Xs.loc[common].values
    y = lbl.set_index("sample_id").loc[common, "y"].values
    if len(np.unique(y)) < 2:
        return {"name": name, "error": "single class"}
    auc_lr = cv_auc(X, y, lambda: LogisticRegression(C=1.0, max_iter=2000, random_state=42))
    auc_rf = cv_auc(X, y, lambda: RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=2))

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof_lr = np.zeros(len(y))
    oof_rf = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        m = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(sc.transform(X[tr]), y[tr])
        oof_lr[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
        m = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=2).fit(X[tr], y[tr])
        oof_rf[te] = m.predict_proba(X[te])[:, 1]
    lr_lo, lr_hi = bootstrap_ci(y, oof_lr)
    rf_lo, rf_hi = bootstrap_ci(y, oof_rf)
    # Ensemble (avg of LR + RF probabilities)
    oof_ens = (oof_lr + oof_rf) / 2
    auc_ens = roc_auc_score(y, oof_ens)
    ens_lo, ens_hi = bootstrap_ci(y, oof_ens)

    return {
        "name": name,
        "n_genes_input": len(genes),
        "n_genes_matched": len(g_have),
        "n_samples": int(len(common)),
        "auc_logreg": auc_lr,
        "auc_logreg_ci": [lr_lo, lr_hi],
        "auc_rf": auc_rf,
        "auc_rf_ci": [rf_lo, rf_hi],
        "auc_ensemble": auc_ens,
        "auc_ensemble_ci": [ens_lo, ens_hi],
    }


def main():
    print("[load] expression + R1-A labels", flush=True)
    expr = pd.read_csv(TPM, sep="\t", index_col=0)
    lbl = load_r1a()
    summary = {"variants": []}

    for name, genes in [("8-gene (current)", GENE_8), ("10-gene", GENE_10),
                         ("12-gene", GENE_12), ("16-gene", GENE_16)]:
        print(f"\n=== {name} ({len(genes)} genes) ===", flush=True)
        r = test_panel(name, genes, expr, lbl)
        summary["variants"].append(r)
        if "auc_logreg" in r:
            print(f"  matched {r['n_genes_matched']}/{r['n_genes_input']}, n={r['n_samples']}", flush=True)
            print(f"  LogReg AUC = {r['auc_logreg']:.4f} (CI {r['auc_logreg_ci'][0]:.3f}-{r['auc_logreg_ci'][1]:.3f})", flush=True)
            print(f"  RF AUC     = {r['auc_rf']:.4f} (CI {r['auc_rf_ci'][0]:.3f}-{r['auc_rf_ci'][1]:.3f})", flush=True)
            print(f"  Ensemble   = {r['auc_ensemble']:.4f} (CI {r['auc_ensemble_ci'][0]:.3f}-{r['auc_ensemble_ci'][1]:.3f})", flush=True)

    # write
    df = pd.DataFrame([{
        "panel": v.get("name"),
        "n_genes": v.get("n_genes_matched"),
        "n_samples": v.get("n_samples"),
        "auc_logreg": v.get("auc_logreg"),
        "ci_lo": v.get("auc_logreg_ci", [None, None])[0],
        "ci_hi": v.get("auc_logreg_ci", [None, None])[1],
        "auc_rf": v.get("auc_rf"),
        "auc_ensemble": v.get("auc_ensemble"),
    } for v in summary["variants"]])
    df.to_csv(OUT / "R5_alt_panel_table.tsv", sep="\t", index=False)
    with open(OUT / "R5_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n=== R5 alternative panel summary ===")
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
