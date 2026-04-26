#!/usr/bin/env python3
"""v8 deep pathway DIAL — per-pathway LODO DIAL distribution across the 50
MSigDB Hallmark pathways for THCA.

Question: is the pooled DIAL=0.000 at 50-pathway level because every pathway
individually has DIAL≈0 (coherent biology) or because aggregation across
pathways dilutes individual flips (some pathways flip, others don't, pool
cancels out)?

For each of the 50 Hallmark pathways j:
  - Feature vector: pathway_j activity score (392 samples)
  - LogReg_l2 inside a LODO loop on B (TCGA-THCA vs GSE27155)
  - AUC_pre uses the raw pathway score
  - AUC_post uses the pathway score after _combat_preserve with Y as covariate
  - DIAL_j = AUC_pre - max(AUC_post, 1 - AUC_post)  (same as v5.1 compute_dial)

Outputs:
  results/v8_statgen/v8_pathway_perpath_dial.tsv   — one row per pathway
  results/v8_statgen/v8_pathway_deep_analysis.md   — narrative interpretation
"""
from __future__ import annotations
import sys, os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score

ROOT = Path("/opt/thyroid-dash/project")
sys.path.insert(0, str(ROOT / "notebooks_or_scripts"))
from v5p1_common import _combat_preserve  # noqa: E402

RES  = ROOT / "results" / "v8_statgen"
DATA = ROOT / "data_processed" / "v5_cross_cancer" / "THCA"

def load_pathway_activity():
    df = pd.read_csv(RES / "v8_pathway_activity_THCA.tsv", sep="\t", index_col=0)
    # rows = samples, columns = 50 pathways
    sample_names = df.index.values
    pathways = df.columns.values
    X = df.values.astype(np.float32)  # (n_samples, 50)
    return sample_names, pathways, X


def load_labels_and_batch(sample_names: np.ndarray):
    Y_all = np.loadtxt(DATA / "Y.tsv", dtype=str)
    B_all = np.loadtxt(DATA / "B.tsv", dtype=str)
    sn_all = np.loadtxt(DATA / "sample_names.txt", dtype=str)
    idx = {s: i for i, s in enumerate(sn_all)}
    order = np.array([idx[s] for s in sample_names])
    return Y_all[order], B_all[order]


def per_pathway_dial(X: np.ndarray, Y: np.ndarray, B: np.ndarray) -> pd.DataFrame:
    """Per-pathway DIAL using the *jointly* ComBat-corrected 50-D matrix.

    We first run ComBat on the full 50-pathway matrix (matching v8 Task 3 so
    the pooled DIAL=0.000 result is directly decomposable), then compute
    per-column AUC_pre (from raw activity) and AUC_post (from ComBat column)
    under LODO. DIAL = AUC_pre - AUC_post_flipped per the v5.1 rule. This
    shows which pathways individually carried the flip and which did not."""
    classes = sorted(np.unique(Y).tolist())
    Ybin = (Y == classes[0]).astype(int)
    cv = list(LeaveOneGroupOut().split(X, Ybin, groups=B))

    # ComBat once on the full 50-D pathway matrix (same as v8 Task 3).
    X_post = _combat_preserve(X, Y, B)
    assert np.isfinite(X_post).all(), "ComBat on 50-D pathway matrix produced NaN"

    n_path = X.shape[1]
    out = []
    for j in range(n_path):
        xj_raw  = X[:, j:j+1]
        xj_post = X_post[:, j:j+1]

        pre, post = [], []
        for tr, te in cv:
            try:
                c_pre = LogisticRegression(max_iter=2000, C=1.0)
                c_pre.fit(xj_raw[tr], Ybin[tr])
                pre.append(roc_auc_score(Ybin[te], c_pre.predict_proba(xj_raw[te])[:, 1]))
                c_post = LogisticRegression(max_iter=2000, C=1.0)
                c_post.fit(xj_post[tr], Ybin[tr])
                post.append(roc_auc_score(Ybin[te], c_post.predict_proba(xj_post[te])[:, 1]))
            except Exception as e:
                pre.append(np.nan); post.append(np.nan)

        auc_pre = float(np.mean(pre))
        auc_post = float(np.mean(post))
        auc_post_flip = max(auc_post, 1 - auc_post)
        dial = float(max(0.0, auc_pre - auc_post_flip))
        flipped = auc_post < 0.5

        if dial >= 0.3:
            interp = "batch_entangled"
        elif auc_pre >= 0.7 and auc_post >= 0.7:
            interp = "true_biology"
        elif auc_pre < 0.6:
            interp = "no_signal"
        else:
            interp = "ambiguous"

        out.append(dict(pathway="", auc_pre=auc_pre, auc_post=auc_post,
                        auc_post_flip=auc_post_flip, dial=dial,
                        flipped_direction=flipped, interpretation=interp))
    return out


def main():
    sample_names, pathways, X = load_pathway_activity()
    Y, B = load_labels_and_batch(sample_names)
    print(f"[deep] X={X.shape} Y_uniq={dict(zip(*np.unique(Y, return_counts=True)))} "
          f"B_uniq={dict(zip(*np.unique(B, return_counts=True)))}", flush=True)

    rows = per_pathway_dial(X, Y, B)
    df = pd.DataFrame(rows)
    df["pathway"] = pathways
    df = df[["pathway","auc_pre","auc_post","auc_post_flip","dial",
             "flipped_direction","interpretation"]]
    df = df.sort_values("dial", ascending=False).reset_index(drop=True)

    out_tsv = RES / "v8_pathway_perpath_dial.tsv"
    df.to_csv(out_tsv, sep="\t", index=False,
              float_format="%.4f")
    print(f"[deep] wrote {out_tsv}", flush=True)

    # Summary stats
    n = len(df)
    n_flip = int((df["dial"] >= 0.3).sum())
    n_bio  = int((df["interpretation"] == "true_biology").sum())
    n_null = int((df["interpretation"] == "no_signal").sum())
    n_amb  = int((df["interpretation"] == "ambiguous").sum())
    print(f"[deep] n_pathways={n}  batch_entangled={n_flip}  true_biology={n_bio}  "
          f"no_signal={n_null}  ambiguous={n_amb}", flush=True)
    print(f"[deep] DIAL distribution: min={df.dial.min():.4f} "
          f"median={df.dial.median():.4f} mean={df.dial.mean():.4f} "
          f"max={df.dial.max():.4f}", flush=True)
    print(f"[deep] top-5 DIAL pathways:", flush=True)
    print(df[["pathway","auc_pre","auc_post","dial","interpretation"]].head().to_string(index=False),
          flush=True)

if __name__ == "__main__":
    main()
