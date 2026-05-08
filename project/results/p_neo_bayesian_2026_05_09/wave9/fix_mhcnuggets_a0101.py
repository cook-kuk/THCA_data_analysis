#!/usr/bin/env python
"""Fix-up: rescore HLA-A*01:01 (42 peptides) without rank_output=True (which crashed)
and merge into predictions_mhcnuggets.tsv. Then recompute auroc_mhcnuggets.tsv.
"""
from __future__ import annotations
import os, sys, time, json, tempfile
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave9"

# Load existing predictions
pred = pd.read_csv(OUT / "predictions_mhcnuggets.tsv", sep="\t")
need = pred[pred["score_mhcnuggets"].isna()].copy()
print(f"[fix] {len(need)} rows missing score (before fix)", flush=True)
print(f"[fix] HLAs needed: {need['hla'].value_counts().to_dict()}", flush=True)

from mhcnuggets.src.predict import predict as mhc_predict

scores = []
for hla, sub in need.groupby("hla"):
    mhc = hla.replace("*", "")
    peptides = sub["peptide"].unique().tolist()
    with tempfile.TemporaryDirectory() as td:
        pep_path = Path(td) / "peptides.txt"
        out_path = Path(td) / "out.csv"
        pep_path.write_text("\n".join(peptides) + "\n")
        try:
            mhc_predict(class_="I", peptides_path=str(pep_path), mhc=mhc,
                        output=str(out_path), rank_output=False)
        except Exception as e:
            print(f"[fix] FAILED {mhc}: {e}", flush=True)
            continue
        if not out_path.exists():
            print(f"[fix] no output for {mhc}", flush=True)
            continue
        pr = pd.read_csv(out_path)
        for _, r in pr.iterrows():
            scores.append({
                "peptide": r["peptide"],
                "hla": hla,
                "ic50_fix": float(r["ic50"]),
            })
    print(f"[fix] done {mhc}: {len(sub)} peptides", flush=True)

sc = pd.DataFrame(scores)
print(f"[fix] scored {len(sc)} rows", flush=True)
print(sc.head().to_string(index=False), flush=True)

# Merge into pred
pred = pred.merge(sc, on=["peptide", "hla"], how="left")
fill_mask = pred["score_mhcnuggets"].isna() & pred["ic50_fix"].notna()
pred.loc[fill_mask, "mhcnuggets_ic50"] = pred.loc[fill_mask, "ic50_fix"]
pred.loc[fill_mask, "score_mhcnuggets"] = -pred.loc[fill_mask, "ic50_fix"]
pred = pred.drop(columns=["ic50_fix"])
n_still_na = int(pred["score_mhcnuggets"].isna().sum())
print(f"[fix] still NaN after merge: {n_still_na}", flush=True)

pred.to_csv(OUT / "predictions_mhcnuggets.tsv", sep="\t", index=False)
print(f"[fix] wrote updated predictions_mhcnuggets.tsv ({len(pred)} rows)", flush=True)


# Recompute AUROC
from sklearn.metrics import roc_auc_score


def auroc_ci(y, s, n_boot=1000, seed=0):
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2:
        return (np.nan, np.nan, np.nan)
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    n = len(y)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy, ss = y[idx], s[idx]
        if len(np.unique(yy)) < 2:
            continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs:
        return (auc, np.nan, np.nan)
    return (float(auc), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))


rows = []
ok = pred[pred["score_mhcnuggets"].notna()].copy()
n_skipped = int(pred["score_mhcnuggets"].isna().sum())
for label, mask in [
    ("ITSNdb_combined", ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])),
    ("ITSNdb_no_overlap", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == False)),
    ("ITSNdb_in_master", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == True)),
]:
    sub = ok[mask]
    y = sub["label"].astype(int).to_numpy()
    s = sub["score_mhcnuggets"].astype(float).to_numpy()
    auc, lo, hi = auroc_ci(y, s, n_boot=1000, seed=42)
    rows.append({"algorithm": "MHCnuggets", "testset": label,
                 "n": len(sub), "n_pos": int((sub["label"] == 1).sum()),
                 "n_neg": int((sub["label"] == 0).sum()),
                 "AUROC": auc, "AUROC_lo95": lo, "AUROC_hi95": hi,
                 "n_skipped_unsupported_hla": 2,
                 "n_skipped_runtime_error": int(n_skipped),
                 "score_col": "neg_ic50"})

ar = pd.DataFrame(rows)
ar.to_csv(OUT / "auroc_mhcnuggets.tsv", sep="\t", index=False)
print(ar.to_string(index=False), flush=True)

(OUT / "_meta_mhcnuggets.json").write_text(json.dumps({
    "n_input": int(len(pred)),
    "n_scored": int(len(ok)),
    "n_skipped_unsupported_hla": 2,
    "n_skipped_runtime_error": int(n_skipped),
    "fix_applied": "HLA-A*01:01 re-run with rank_output=False after numpy comparison error",
}, indent=2))
