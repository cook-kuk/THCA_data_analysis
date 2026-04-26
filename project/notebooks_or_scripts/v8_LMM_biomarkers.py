#!/usr/bin/env python3
"""v8 Task 8C — LMM-based biomarker correction on THCA.

Per-gene statsmodels MixedLM: expression ~ subtype + (1 | cohort)
Compare naive OLS (no random effect) vs. LMM; apply BH-FDR to each.
Status: same / lost_with_LMM / gained_with_LMM.

Reference: Sul et al. (2018) Genome Biology — LMMs for multi-tissue /
multi-cohort DE analysis.

Outputs: /opt/thyroid-dash/project/results/v8_statgen/v8_LMM_corrected_biomarkers.tsv
Columns: gene, beta_subtype, pvalue_naive, fdr_naive, pvalue_LMM, fdr_LMM, status.
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")

PROJECT = Path("/opt/thyroid-dash/project")
THCA_DIR = PROJECT / "data_processed" / "v5_cross_cancer" / "THCA"
RESULTS = PROJECT / "results" / "v8_statgen"
LOGFILE = PROJECT / "logs" / "v8_hanlab.log"

N_TOP = 3000
FDR_THRESH = 0.05


def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] [T8C] {msg}"
    print(line, flush=True)
    with open(LOGFILE, "a") as fh:
        fh.write(line + "\n")


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    fdr = np.empty(n, dtype=float)
    fdr[:] = np.nan
    valid = ~np.isnan(ranked)
    ranks = np.arange(1, n + 1)
    adj = ranked * n / ranks
    # enforce monotonicity from the right
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    fdr_ordered = adj
    fdr[order] = fdr_ordered
    fdr[np.isnan(pvals)] = np.nan
    return fdr


def fit_one_gene(args):
    """Return (beta, p_naive, p_lmm) for a single gene."""
    g_idx, expr, subtype_bin, cohort_codes = args
    beta_lmm = np.nan
    p_naive = np.nan
    p_lmm = np.nan
    try:
        df = pd.DataFrame(dict(expr=expr, subtype=subtype_bin, cohort=cohort_codes))
        # Naive OLS: expr ~ subtype
        X_ols = sm.add_constant(df[["subtype"]].astype(float))
        ols = sm.OLS(df["expr"].astype(float), X_ols).fit()
        p_naive = float(ols.pvalues.get("subtype", np.nan))
    except Exception:
        pass
    try:
        df = pd.DataFrame(dict(expr=expr, subtype=subtype_bin, cohort=cohort_codes))
        # LMM: random intercept for cohort
        md = smf.mixedlm("expr ~ subtype", df, groups=df["cohort"])
        mdf = md.fit(method="lbfgs", reml=True, disp=False)
        beta_lmm = float(mdf.params.get("subtype", np.nan))
        p_lmm = float(mdf.pvalues.get("subtype", np.nan))
    except Exception:
        pass
    return g_idx, beta_lmm, p_naive, p_lmm


def main() -> None:
    log("loading X/Y/B for THCA")
    X = np.load(THCA_DIR / "X_combined.npz")["X"]
    Y = np.loadtxt(THCA_DIR / "Y.tsv", dtype=str)
    B = np.loadtxt(THCA_DIR / "B.tsv", dtype=str)
    genes = np.loadtxt(THCA_DIR / "shared_genes.txt", dtype=str)
    log(f"X {X.shape}, genes {genes.shape}")

    # top-variance filter
    var_per_gene = X.var(axis=0)
    top_idx = np.argsort(var_per_gene)[-N_TOP:]
    X = X[:, top_idx]
    genes = genes[top_idx]
    log(f"top-variance filter: {N_TOP} genes")

    # BRAF=1, RAS=0 (binary)
    subtype_bin = (Y == "BRAF").astype(int)
    cohort_codes = pd.Categorical(B).codes.astype(int)
    log(f"subtype counts: BRAF={int((subtype_bin==1).sum())}, RAS={int((subtype_bin==0).sum())}")
    log(f"cohort codes: {np.unique(cohort_codes, return_counts=True)}")

    tasks = []
    for j in range(X.shape[1]):
        tasks.append((j, X[:, j].astype(float), subtype_bin, cohort_codes))

    log(f"fitting {len(tasks)} gene models in parallel (joblib n_jobs=-1)")
    t0 = time.time()
    results = Parallel(n_jobs=-1, backend="loky", verbose=5)(
        delayed(fit_one_gene)(t) for t in tasks
    )
    dt = time.time() - t0
    log(f"LMM+OLS fit done in {dt:.1f}s")

    betas = np.full(X.shape[1], np.nan)
    p_naive = np.full(X.shape[1], np.nan)
    p_lmm = np.full(X.shape[1], np.nan)
    for g_idx, b, pn, pl in results:
        betas[g_idx] = b
        p_naive[g_idx] = pn
        p_lmm[g_idx] = pl

    fdr_naive = bh_fdr(p_naive)
    fdr_lmm = bh_fdr(p_lmm)

    status = []
    for fn, fl in zip(fdr_naive, fdr_lmm):
        sig_n = (not np.isnan(fn)) and (fn < FDR_THRESH)
        sig_l = (not np.isnan(fl)) and (fl < FDR_THRESH)
        if sig_n and sig_l:
            status.append("same")
        elif sig_n and not sig_l:
            status.append("lost_with_LMM")
        elif not sig_n and sig_l:
            status.append("gained_with_LMM")
        else:
            status.append("same")  # both non-sig
    status = np.asarray(status)

    df = pd.DataFrame(dict(
        gene=genes,
        beta_subtype=betas,
        pvalue_naive=p_naive,
        fdr_naive=fdr_naive,
        pvalue_LMM=p_lmm,
        fdr_LMM=fdr_lmm,
        status=status,
    ))
    out_path = RESULTS / "v8_LMM_corrected_biomarkers.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    n_naive_sig = int(((fdr_naive < FDR_THRESH) & ~np.isnan(fdr_naive)).sum())
    n_lmm_sig = int(((fdr_lmm < FDR_THRESH) & ~np.isnan(fdr_lmm)).sum())
    n_same_sig = int(((fdr_naive < FDR_THRESH) & (fdr_lmm < FDR_THRESH)).sum())
    n_lost = int(((fdr_naive < FDR_THRESH) & ~(fdr_lmm < FDR_THRESH)).sum())
    n_gained = int((~(fdr_naive < FDR_THRESH) & (fdr_lmm < FDR_THRESH)).sum())

    log(f"n_naive_sig(FDR<{FDR_THRESH}) = {n_naive_sig}")
    log(f"n_LMM_sig(FDR<{FDR_THRESH})  = {n_lmm_sig}")
    log(f"same_sig={n_same_sig}  lost_with_LMM={n_lost}  gained_with_LMM={n_gained}")
    log(f"wrote {out_path}")
    log("T8C DONE")


if __name__ == "__main__":
    main()
