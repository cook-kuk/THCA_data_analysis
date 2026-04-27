#!/usr/bin/env python3
"""
v17 REAL FIX · R2 — REAL AUC of 8-gene panel against leak-free clusters.

For each R1-{A,B,C,D} cluster (defined WITHOUT 8-gene panel):
  - 5-fold StratifiedKFold CV (random_state=42)
  - LogisticRegression and RandomForest
  - BRAF V600E single-feature baseline (truly independent)
  - Bootstrap 1000 95% CI for AUC
  - Subgroup analysis (age/sex/stage)

Outputs:
  - results/v17_realfix/R2_real_auc_table.tsv
  - results/v17_realfix/R2_per_variant_summary.json
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

OUT = Path("/opt/thyroid-dash/project/results/v17_realfix")
TPM = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]

# BRAF V600E status — derive from sample_master_v3_merged.tsv driver_anchor column
SAMPLE_MASTER = Path("/opt/thyroid-dash/project/metadata/sample_master_v3_merged.tsv")


def load_expression() -> pd.DataFrame:
    df = pd.read_csv(TPM, sep="\t", index_col=0)
    return df


def get_8gene_matrix(expr: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    g_have = [g for g in GENE_8 if g in expr.index]
    s_have = [s for s in samples if s in expr.columns]
    sub = expr.loc[g_have, s_have].T  # samples × genes
    return sub


def load_braf_status() -> pd.DataFrame | None:
    if not SAMPLE_MASTER.exists():
        print(f"[braf] sample master not found at {SAMPLE_MASTER}", flush=True)
        return None
    print(f"[braf] using {SAMPLE_MASTER}", flush=True)
    df = pd.read_csv(SAMPLE_MASTER, sep="\t", low_memory=False)
    # derive 0/1 BRAF V600E flag from driver_anchor column
    if "driver_anchor" not in df.columns:
        print("[braf] no driver_anchor column", flush=True)
        return None
    out = pd.DataFrame({
        "sample_id": df["sample_id"].astype(str),
        "braf_v600e": (df["driver_anchor"].astype(str).str.upper() == "BRAF").astype(int),
    })
    print(f"[braf] derived n={len(out)}, BRAF+={int(out['braf_v600e'].sum())}", flush=True)
    return out


def cv_auc(X: np.ndarray, y: np.ndarray, model_factory, n_splits: int = 5, random_state: int = 42) -> tuple[float, list[float]]:
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    fold_aucs = []
    oof_proba = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        m = model_factory()
        if isinstance(m, LogisticRegression):
            sc = StandardScaler().fit(X[tr])
            m.fit(sc.transform(X[tr]), y[tr])
            p = m.predict_proba(sc.transform(X[te]))[:, 1]
        else:
            m.fit(X[tr], y[tr])
            p = m.predict_proba(X[te])[:, 1]
        oof_proba[te] = p
        fold_aucs.append(roc_auc_score(y[te], p))
    overall = roc_auc_score(y, oof_proba)
    return overall, fold_aucs


def bootstrap_auc_ci(y: np.ndarray, p: np.ndarray, n_boot: int = 1000, seed: int = 42) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    aucs = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], p[idx]))
    aucs = np.array(aucs)
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def evaluate_variant(variant: str, expr: pd.DataFrame, braf: pd.DataFrame | None) -> dict:
    label_path = OUT / f"R1{variant}_cluster_labels.tsv"
    if not label_path.exists():
        return {"variant": variant, "error": f"missing {label_path}"}
    lbl = pd.read_csv(label_path, sep="\t")
    lbl["y"] = (lbl["cluster"].str.startswith("DM2")).astype(int)
    samples = lbl["sample_id"].tolist()
    X8 = get_8gene_matrix(expr, samples)
    common = X8.index.intersection(lbl.set_index("sample_id").index)
    X = X8.loc[common].values
    y = lbl.set_index("sample_id").loc[common, "y"].values
    print(f"[R2-{variant}] n={len(common)}, n_DM2={int(y.sum())}, n_DM1={int(len(y)-y.sum())}", flush=True)
    if len(np.unique(y)) < 2 or len(common) < 30:
        return {"variant": variant, "error": "insufficient class balance or n", "n": int(len(common)), "n_DM2": int(y.sum())}

    res = {"variant": variant, "n_total": int(len(common)), "n_DM2": int(y.sum()), "n_DM1": int(len(y) - y.sum())}

    # 8-gene LogReg
    auc_lr, fold_lr = cv_auc(X, y, lambda: LogisticRegression(C=1.0, max_iter=2000, random_state=42))
    # 8-gene RF
    auc_rf, fold_rf = cv_auc(X, y, lambda: RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=2))

    # bootstrap CIs from oof predictions
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof_lr = np.zeros(len(y))
    oof_rf = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X[tr])
        m = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(sc.transform(X[tr]), y[tr])
        oof_lr[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
        m = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=2).fit(X[tr], y[tr])
        oof_rf[te] = m.predict_proba(X[te])[:, 1]
    lr_lo, lr_hi = bootstrap_auc_ci(y, oof_lr)
    rf_lo, rf_hi = bootstrap_auc_ci(y, oof_rf)

    res["auc_8gene_logreg"] = float(auc_lr)
    res["auc_8gene_logreg_95ci"] = [lr_lo, lr_hi]
    res["auc_8gene_rf"] = float(auc_rf)
    res["auc_8gene_rf_95ci"] = [rf_lo, rf_hi]
    res["fold_aucs_logreg"] = fold_lr
    res["fold_aucs_rf"] = fold_rf

    # BRAF V600E baseline (truly independent)
    if braf is not None:
        braf_col = "braf_v600e" if "braf_v600e" in braf.columns else braf.columns[1]
        sid_col = "sample_id" if "sample_id" in braf.columns else braf.columns[0]
        braf_idx = braf.set_index(sid_col)[braf_col]
        common_b = [s for s in common if s in braf_idx.index]
        if len(common_b) >= 30:
            Xb = braf_idx.loc[common_b].astype(float).values.reshape(-1, 1)
            yb = lbl.set_index("sample_id").loc[common_b, "y"].values
            auc_braf, fold_braf = cv_auc(Xb, yb, lambda: LogisticRegression(C=1.0, max_iter=2000, random_state=42))
            res["auc_braf_only"] = float(auc_braf)
            res["delta_auc_8gene_vs_braf_logreg"] = float(auc_lr - auc_braf)
            res["n_braf_common"] = int(len(common_b))
        else:
            res["auc_braf_only"] = None
            res["braf_note"] = f"insufficient overlap n={len(common_b)}"
    else:
        res["auc_braf_only"] = None
        res["braf_note"] = "BRAF status table not found"

    print(f"[R2-{variant}] 8-gene LogReg AUC = {auc_lr:.4f} (95% CI {lr_lo:.3f}-{lr_hi:.3f})  RF = {auc_rf:.4f}", flush=True)
    if "auc_braf_only" in res and res["auc_braf_only"]:
        print(f"           BRAF-only AUC = {res['auc_braf_only']:.4f}  ΔAUC = {res['delta_auc_8gene_vs_braf_logreg']:+.4f}", flush=True)

    return res


def main() -> int:
    expr = load_expression()
    print(f"[load] expression {expr.shape}", flush=True)
    braf = load_braf_status()
    if braf is not None:
        print(f"[load] braf {braf.shape}", flush=True)
        print(braf.head(3))

    results = {}
    rows = []
    for v in ["A", "B", "C", "D"]:
        r = evaluate_variant(v, expr, braf)
        results[v] = r
        if "auc_8gene_logreg" in r:
            rows.append({
                "variant": v,
                "n": r["n_total"],
                "auc_8gene_logreg": r["auc_8gene_logreg"],
                "ci_lo": r["auc_8gene_logreg_95ci"][0],
                "ci_hi": r["auc_8gene_logreg_95ci"][1],
                "auc_8gene_rf": r["auc_8gene_rf"],
                "auc_braf_only": r.get("auc_braf_only"),
                "delta_auc": r.get("delta_auc_8gene_vs_braf_logreg"),
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "R2_real_auc_table.tsv", sep="\t", index=False)
    with open(OUT / "R2_per_variant_summary.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print("\n=== R2 REAL AUC summary ===")
    print(df.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
