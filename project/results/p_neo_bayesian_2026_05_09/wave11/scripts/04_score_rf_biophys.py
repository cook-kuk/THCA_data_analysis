#!/usr/bin/env python3
"""Wave 11 — train RF biophys on master pool n=2396 and score all bundles.

This is the in-house RF baseline (`fit_rf` in cancer_vaccine_robustness/_common.py).
Features: biophys (24-d) + HLA one-hot.
Trained on bundle.tsv split=='train'; scored on every wave11 test bundle.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import time

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
BUN = WAVE11 / "test_bundles"
OUT = WAVE11 / "wave11_predictions"

sys.path.insert(0, "/home/seungho/personal/THCA_data_analysis/project/results/cancer_vaccine_robustness_2026_05_09")
from _common import biophys, build_hla_onehot_factory, fit_rf, fit_logreg, normalize_hla

t0 = time.time()
print("[load] master pool n=2396 from bundle.tsv ...", flush=True)
bun = pd.read_csv(ROOT / "bundle.tsv", sep="\t")
train = bun[bun["split"] == "train"].dropna(subset=["HLA_norm"]).copy()
print(f"  train rows: {len(train)}; pos: {int(train['label'].sum())}", flush=True)

# Build HLA encoder
encode, n_dim, hla_idx = build_hla_onehot_factory(train["HLA_norm"].tolist())
def feats(df):
    bp = biophys(df["peptide"].tolist())
    hla = encode(df["hla"].tolist() if "hla" in df.columns else df["HLA_norm"].tolist())
    return np.hstack([bp, hla]).astype(np.float32)

Xtr = feats(train.rename(columns={"HLA_norm": "hla"}))
ytr = train["label"].astype(int).to_numpy()

print(f"[fit] RF on Xtr={Xtr.shape}", flush=True)
rf = fit_rf(Xtr, ytr)
print(f"[fit] LR on Xtr={Xtr.shape}", flush=True)
lr, sc = fit_logreg(Xtr, ytr)

# Score every bundle
bundles = sorted(BUN.glob("*.tsv"))
for bp_path in bundles:
    bname = bp_path.stem
    df = pd.read_csv(bp_path, sep="\t")
    if len(df) == 0:
        print(f"[skip] {bname} empty"); continue
    Xt = feats(df)
    score_rf = rf.predict_proba(Xt)[:, 1]
    score_lr = lr.predict_proba(sc.transform(Xt))[:, 1]
    out = df[["peptide", "hla", "label", "source", "length",
              "n_overlap_with_other_sources"]].copy()
    out["score_rf_biophys"] = score_rf
    out["score_lr_biophys"] = score_lr
    out.to_csv(OUT / f"RF_biophys__{bname}.tsv", sep="\t", index=False)
    out_lr = out.drop(columns=["score_rf_biophys"]).rename(columns={"score_lr_biophys": "score_lr_biophys"})
    # also dump LR separately for completeness, but we mostly use RF
    print(f"[score] {bname}: n={len(df)} -> RF + LR scores written", flush=True)

print(f"[done] total {time.time()-t0:.1f}s")
