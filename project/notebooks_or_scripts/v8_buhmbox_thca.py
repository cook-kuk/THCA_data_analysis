#!/usr/bin/env python3
"""v8 Task 8B — BUHMBOX concept check on THCA.

Analogical to Han et al. (2016) Nat Genet BUHMBOX: if two cohorts of the
same clinical label share the same biological sub-group, their principal
component distributions should coincide. We compute PC1 on pooled class-
specific samples (TCGA + GSE for each class), then run a 2-sample
Kolmogorov-Smirnov test comparing PC1 of TCGA vs GSE27155 within each
class. p>0.05 = same sub-group; p<=0.05 = hidden heterogeneity.

Note: this is conceptual, not the original BUHMBOX statistic.

Outputs: /opt/thyroid-dash/project/results/v8_statgen/v8_buhmbox_concept_check.tsv
Columns: cohort_A, cohort_B, class, n_A, n_B, KS_stat, KS_pvalue, interpretation.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

PROJECT = Path("/opt/thyroid-dash/project")
THCA_DIR = PROJECT / "data_processed" / "v5_cross_cancer" / "THCA"
RESULTS = PROJECT / "results" / "v8_statgen"
LOGFILE = PROJECT / "logs" / "v8_hanlab.log"

N_TOP = 3000
MIN_N = 10


def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] [T8B] {msg}"
    print(line, flush=True)
    with open(LOGFILE, "a") as fh:
        fh.write(line + "\n")


def main() -> None:
    log("loading X/Y/B for THCA")
    X = np.load(THCA_DIR / "X_combined.npz")["X"]
    Y = np.loadtxt(THCA_DIR / "Y.tsv", dtype=str)
    B = np.loadtxt(THCA_DIR / "B.tsv", dtype=str)
    log(f"X {X.shape}")

    # top-variance filter
    var_per_gene = X.var(axis=0)
    top_idx = np.argsort(var_per_gene)[-N_TOP:]
    X = X[:, top_idx]
    log(f"top-variance filter: {N_TOP} genes")

    cohorts = sorted(np.unique(B).tolist())
    classes = sorted(np.unique(Y).tolist())
    log(f"cohorts={cohorts}  classes={classes}")

    rows = []
    for cls in classes:
        # PC1 on pooled class-specific samples (both cohorts together)
        sel = Y == cls
        Xc = X[sel]
        Bc = B[sel]
        log(f"class={cls}: pooled n={len(Xc)}  per-cohort={dict(zip(*np.unique(Bc, return_counts=True)))}")
        if len(Xc) < 4:
            log(f"  skip: too few pooled samples")
            continue
        # center + scale then PC1
        Xs = StandardScaler(with_mean=True, with_std=True).fit_transform(Xc)
        pca = PCA(n_components=1, random_state=42).fit(Xs)
        pc1 = pca.transform(Xs).ravel()
        # within-class PC1 by cohort
        for i, a in enumerate(cohorts):
            for b in cohorts[i + 1:]:
                pc_a = pc1[Bc == a]
                pc_b = pc1[Bc == b]
                n_a, n_b = len(pc_a), len(pc_b)
                if n_a < MIN_N or n_b < MIN_N:
                    ks_stat = float("nan")
                    pv = float("nan")
                    interp = "insufficient_n"
                    log(f"  {a} vs {b} ({cls}): n_A={n_a}, n_B={n_b} -> NaN")
                else:
                    r = stats.ks_2samp(pc_a, pc_b)
                    ks_stat = float(r.statistic)
                    pv = float(r.pvalue)
                    interp = "same_subgroup" if pv > 0.05 else "hidden_heterogeneity"
                    log(f"  {a} vs {b} ({cls}): n_A={n_a}, n_B={n_b}  KS={ks_stat:.4f}  p={pv:.4g}  -> {interp}")
                rows.append(dict(
                    cohort_A=a,
                    cohort_B=b,
                    **{"class": cls},
                    n_A=n_a,
                    n_B=n_b,
                    KS_stat=ks_stat,
                    KS_pvalue=pv,
                    interpretation=interp,
                ))

    out = pd.DataFrame(rows)
    out_path = RESULTS / "v8_buhmbox_concept_check.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    log(f"wrote {out_path}  rows={len(out)}")
    log("T8B DONE")


if __name__ == "__main__":
    main()
