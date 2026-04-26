#!/usr/bin/env python3
"""v8.1 Task B — full-gene LMM (no 3 000-gene pre-filter).

v8 S6C fit `expression ~ subtype + (1|cohort)` on 3 000 top-variance genes
from 11 710 shared. Report claimed 99.4 % biomarker survival — but the
denominator was misleading: the 99.4 % was of the SUBSET testable. The
honest number requires testing every one of the 11 710 shared genes.

This script fits MixedLM per gene across all 11 710 shared genes (no
variance pre-filter) and re-computes biomarker survival against the
published 2 773 list.

Outputs (results/v8p1_rigor/b_full_lmm/):
  lmm_all_11710_genes.tsv       — per-gene beta, p_OLS, p_LMM, FDR_OLS, FDR_LMM
  biomarker_survival_honest.tsv — crossover vs 2 773 biomarkers
  summary.md                    — narrative

Log: logs/v8p1_full_lmm.log
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.regression.linear_model import OLS
import statsmodels.api as sm

ROOT = Path("/opt/thyroid-dash/project")
DATA = ROOT / "data_processed" / "v5_cross_cancer" / "THCA"
OUT  = ROOT / "results" / "v8p1_rigor" / "b_full_lmm"
LOG  = ROOT / "logs" / "v8p1_full_lmm.log"
OUT.mkdir(parents=True, exist_ok=True)

_log_f = open(LOG, "w")
def log(msg: str) -> None:
    stamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
    line = f"{stamp} {msg}"
    print(line, flush=True)
    _log_f.write(line + "\n")
    _log_f.flush()

log("v8.1 Task B — full-gene LMM on 11 710 shared THCA genes")
log("=" * 60)

# ------------------------------------------------------------------
# Load
# ------------------------------------------------------------------
X = np.load(DATA / "X_combined.npz")["X"].astype(np.float32)  # (392, 11710)
Y = np.loadtxt(DATA / "Y.tsv", dtype=str)      # subtype
B = np.loadtxt(DATA / "B.tsv", dtype=str)      # cohort
genes = np.loadtxt(DATA / "shared_genes.txt", dtype=str)
sample_names = np.loadtxt(DATA / "sample_names.txt", dtype=str)

log(f"[load] X={X.shape} genes={len(genes)} Y_uniq={dict(zip(*np.unique(Y, return_counts=True)))} "
    f"B_uniq={dict(zip(*np.unique(B, return_counts=True)))}")

# Encode subtype as 0/1 (BRAF=1, RAS=0); cohort integer codes
subtype_bin = (Y == "BRAF").astype(int)
cohort_codes, cohort_labels = pd.factorize(B)
log(f"[encode] subtype mean={subtype_bin.mean():.3f}  "
    f"cohort codes: {list(zip(cohort_labels, np.bincount(cohort_codes)))}")

# ------------------------------------------------------------------
# Per-gene fit
# ------------------------------------------------------------------
def fit_one(idx: int) -> tuple:
    """Return (beta_subtype, p_OLS, p_LMM) for gene idx."""
    y = X[:, idx].astype(np.float64)
    # Drop NaNs defensively
    if not np.isfinite(y).all():
        return (np.nan, np.nan, np.nan)
    # OLS: y ~ subtype
    x_ols = np.column_stack([np.ones_like(subtype_bin), subtype_bin]).astype(float)
    try:
        ols = OLS(y, x_ols).fit()
        beta = float(ols.params[1])
        p_ols = float(ols.pvalues[1])
    except Exception:
        beta, p_ols = np.nan, np.nan
    # LMM: y ~ subtype + (1|cohort)
    try:
        md = MixedLM(
            endog=y,
            exog=x_ols,
            groups=cohort_codes,
        )
        mdf = md.fit(reml=True, method="lbfgs", disp=False)
        # subtype coef is index 1 (intercept is 0)
        p_lmm = float(mdf.pvalues[1])
        beta_lmm = float(mdf.params[1])
        # Use LMM beta as authoritative beta
        return (beta_lmm, p_ols, p_lmm)
    except Exception:
        return (beta, p_ols, np.nan)

log("[fit] starting per-gene LMM/OLS, n_jobs=-1, batch=11710")
t0 = time.time()
results = Parallel(n_jobs=-1, verbose=5, batch_size=200)(
    delayed(fit_one)(i) for i in range(len(genes))
)
dt = time.time() - t0
log(f"[fit] done in {dt:.1f}s "
    f"({dt/len(genes)*1000:.1f}ms/gene)")

betas, pvals_ols, pvals_lmm = zip(*results)
df = pd.DataFrame({
    "gene": genes,
    "beta_subtype": betas,
    "pvalue_OLS": pvals_ols,
    "pvalue_LMM": pvals_lmm,
})
# BH-FDR per column
from statsmodels.stats.multitest import multipletests
for col_p, col_fdr in [("pvalue_OLS", "fdr_OLS"),
                       ("pvalue_LMM", "fdr_LMM")]:
    mask = df[col_p].notna()
    fdr = np.full(len(df), np.nan)
    if mask.any():
        _, p_adj, _, _ = multipletests(df.loc[mask, col_p].values, method="fdr_bh")
        fdr[mask.values] = p_adj
    df[col_fdr] = fdr

# Status classification
def classify(row):
    s_ols = (row["fdr_OLS"] < 0.05)
    s_lmm = (row["fdr_LMM"] < 0.05)
    if s_ols and s_lmm: return "same"
    if (not s_ols) and s_lmm: return "gained_with_LMM"
    if s_ols and (not s_lmm): return "lost_with_LMM"
    return "nonsig_both"

df["status"] = df.apply(classify, axis=1)

n_ols = int((df["fdr_OLS"] < 0.05).sum())
n_lmm = int((df["fdr_LMM"] < 0.05).sum())
status_counts = df["status"].value_counts().to_dict()
log(f"[stat] OLS sig = {n_ols}   LMM sig = {n_lmm}")
log(f"[stat] status breakdown: {status_counts}")

out_tsv = OUT / "lmm_all_11710_genes.tsv"
df.to_csv(out_tsv, sep="\t", index=False, float_format="%.6g")
log(f"[write] {out_tsv}")

# ------------------------------------------------------------------
# Biomarker survival vs 2 773 slate
# ------------------------------------------------------------------
V8DE = ROOT / "results" / "v8_statgen" / "v8_biomarker_DE_recomputation.tsv"
v8_de = pd.read_csv(V8DE, sep="\t")
orig_list = set(v8_de.loc[v8_de["old_significant"], "gene"])
log(f"[biomarker] original published list: {len(orig_list)}")

lmm_universe = set(df["gene"])
lmm_sig = set(df.loc[df["fdr_LMM"] < 0.05, "gene"])

orig_in_universe = orig_list & lmm_universe
orig_surviving = orig_list & lmm_sig

log(f"[biomarker] original ∩ 11 710 LMM universe: {len(orig_in_universe)}"
    f"  ({100*len(orig_in_universe)/len(orig_list):.1f}% of 2 773)")
log(f"[biomarker] of those, LMM-significant: {len(orig_surviving)}"
    f"  ({100*len(orig_surviving)/len(orig_in_universe):.1f}% of in-universe)")
log(f"[biomarker] HONEST survival (surviving / 2 773): "
    f"{100*len(orig_surviving)/len(orig_list):.1f}%")

# Compare to v8 3 000-gene result
V8LMM = ROOT / "results" / "v8_statgen" / "v8_LMM_corrected_biomarkers.tsv"
v8_lmm = pd.read_csv(V8LMM, sep="\t")
v8_lmm_universe = set(v8_lmm["gene"])
v8_lmm_sig = set(v8_lmm.loc[v8_lmm["fdr_LMM"] < 0.05, "gene"])
log(f"[v8-compare] v8 3 000-gene LMM universe: {len(v8_lmm_universe)}")
log(f"[v8-compare] v8 3 000-gene LMM sig: {len(v8_lmm_sig)}")
log(f"[v8-compare] v8 claimed in-universe survival: "
    f"{100*len(orig_list & v8_lmm_sig)/max(1,len(orig_list & v8_lmm_universe)):.1f}%")

# Druggable retention
druggable = {"TACSTD2","TMPRSS4","PLEKHA6","CYP1B1","LDLR","GABRB2","B3GNT3","PTPRE"}
drug_in_universe = druggable & lmm_universe
drug_surviving = druggable & lmm_sig
log(f"[biomarker] druggable in LMM universe: "
    f"{sorted(drug_in_universe)} ({len(drug_in_universe)}/8)")
log(f"[biomarker] druggable LMM-significant:  "
    f"{sorted(drug_surviving)} ({len(drug_surviving)}/8)")

# Build summary
surv_df = pd.DataFrame({
    "metric": [
        "original_biomarkers_total",
        "full_LMM_universe",
        "original_in_LMM_universe",
        "original_surviving_LMM_FDR_0p05",
        "in_universe_survival_rate",
        "HONEST_survival_rate_of_2773",
        "v8_3000_gene_universe_size",
        "v8_3000_gene_in_universe_survival",
        "druggable_targets_in_universe",
        "druggable_targets_surviving_LMM",
        "LMM_sig_total",
        "OLS_sig_total",
        "gained_with_LMM",
        "lost_with_LMM",
        "same_both",
    ],
    "value": [
        len(orig_list),
        len(lmm_universe),
        len(orig_in_universe),
        len(orig_surviving),
        f"{100*len(orig_surviving)/max(1,len(orig_in_universe)):.2f}%",
        f"{100*len(orig_surviving)/len(orig_list):.2f}%",
        len(v8_lmm_universe),
        f"{100*len(orig_list & v8_lmm_sig)/max(1,len(orig_list & v8_lmm_universe)):.2f}%",
        f"{len(drug_in_universe)}/8",
        f"{len(drug_surviving)}/8",
        n_lmm,
        n_ols,
        status_counts.get("gained_with_LMM", 0),
        status_counts.get("lost_with_LMM", 0),
        status_counts.get("same", 0),
    ],
})
surv_df.to_csv(OUT / "biomarker_survival_honest.tsv", sep="\t", index=False)
log(f"[write] {OUT / 'biomarker_survival_honest.tsv'}")

log("v8.1 Task B complete")
_log_f.close()
