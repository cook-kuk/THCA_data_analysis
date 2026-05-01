#!/usr/bin/env python3
"""P5 — Finding 2 autocorrelation: is the 8-gene panel a 'disguised HLA-II / immune' score?

If 8-gene RAI score and HLA-II module are highly Spearman-correlated:
  - either they're independent biological axes that happen to correlate
  - or 8-gene panel is partially a marker of immune infiltration

Method:
  1. TCGA-THCA: Spearman ρ between 8-gene RAI score and HLA-II module
  2. Residualize HLA-II out of 8-gene; re-test DM1 vs DM2 with residual 8-gene
  3. Same for purity (CPE/ABSOLUTE)
  4. GSE286332: same panel ρ; verify direction
  5. 3-cohort meta forest plot table (TCGA + K2 + GSE286332 8-gene effect on Hashimoto-like)
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/p5_8gene_vs_hla_autocorr"
RES.mkdir(parents=True, exist_ok=True)

GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
HLA_I = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"]
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
          "HLA-DMA", "HLA-DMB", "CIITA", "HLA-DOB"]
IMMUNE_PROXY = ["CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD20", "MS4A1", "GZMB", "PRF1",
                "IFNG", "TNF", "IL10", "FOXP3", "PTPRC", "CXCL13"]

def zmean(table, genes):
    keep = [g for g in genes if g in table.index]
    if len(keep) < 3:
        return None, keep
    z = table.loc[keep].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
    return z.mean(axis=0), keep

# ============================================================
# 1. TCGA — load + compute scores
# ============================================================
print("=== 1. TCGA: 8-gene RAI vs HLA-II module ===")
expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                   sep="\t", index_col=0)
labels = pd.read_csv(PROJ / "results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
labels["DM"] = labels["cluster"].str[:3]
samps = [s for s in labels["sample_id"] if s in expr.columns]
labels = labels[labels["sample_id"].isin(samps)].set_index("sample_id").loc[samps]
E = expr[samps]
y = (labels["DM"] == "DM1").astype(int).values  # DM1 = 1

g8_score, g8_used = zmean(E, GENE_8)
hla1_score, _ = zmean(E, HLA_I)
hla2_score, _ = zmean(E, HLA_II)
imm_score, _ = zmean(E, IMMUNE_PROXY)

scores = pd.DataFrame({
    "g8_RAI": g8_score, "HLA_I": hla1_score, "HLA_II": hla2_score, "immune": imm_score,
}, index=samps)
scores["DM"] = labels["DM"].values

print(f"  TCGA samples: {len(samps)}")
print(f"  g8 RAI used: {g8_used}")

corr_table = scores[["g8_RAI", "HLA_I", "HLA_II", "immune"]].corr(method="spearman")
print("\n  Spearman ρ matrix:")
print(corr_table.round(3).to_string())

# ============================================================
# 2. DM1 vs DM2 — raw vs residualized 8-gene
# ============================================================
print("\n=== 2. DM1 vs DM2 — raw vs HLA-II/immune-residualized 8-gene ===")
def cohen_compare(score, y):
    a = score[y == 1]; b = score[y == 0]
    sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    d = (a.mean() - b.mean()) / max(sp, 1e-9)
    _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(d), float(p)

# Raw 8-gene
d_raw, p_raw = cohen_compare(scores["g8_RAI"].values, y)
# Residualize HLA-II out
X = sm.add_constant(scores[["HLA_II"]])
res = sm.OLS(scores["g8_RAI"], X).fit()
g8_resid_hla2 = scores["g8_RAI"].values - res.fittedvalues.values
d_r2, p_r2 = cohen_compare(g8_resid_hla2, y)

# Residualize immune-proxy
X = sm.add_constant(scores[["immune"]])
res2 = sm.OLS(scores["g8_RAI"], X).fit()
g8_resid_imm = scores["g8_RAI"].values - res2.fittedvalues.values
d_ri, p_ri = cohen_compare(g8_resid_imm, y)

# Residualize HLA-II + immune
X = sm.add_constant(scores[["HLA_II", "immune"]])
res3 = sm.OLS(scores["g8_RAI"], X).fit()
g8_resid_both = scores["g8_RAI"].values - res3.fittedvalues.values
d_rb, p_rb = cohen_compare(g8_resid_both, y)

# Reverse direction — HLA-II raw vs residualized for 8-gene
d_h_raw, p_h_raw = cohen_compare(scores["HLA_II"].values, y)
X = sm.add_constant(scores[["g8_RAI"]])
res4 = sm.OLS(scores["HLA_II"], X).fit()
hla2_resid_g8 = scores["HLA_II"].values - res4.fittedvalues.values
d_h_r, p_h_r = cohen_compare(hla2_resid_g8, y)

residual_rows = [
    dict(model="raw 8-gene RAI", d=round(d_raw,3), mw_p=p_raw),
    dict(model="8-gene | HLA-II residualized", d=round(d_r2,3), mw_p=p_r2),
    dict(model="8-gene | immune-proxy residualized", d=round(d_ri,3), mw_p=p_ri),
    dict(model="8-gene | HLA-II + immune residualized", d=round(d_rb,3), mw_p=p_rb),
    dict(model="raw HLA-II", d=round(d_h_raw,3), mw_p=p_h_raw),
    dict(model="HLA-II | 8-gene residualized", d=round(d_h_r,3), mw_p=p_h_r),
]
res_df = pd.DataFrame(residual_rows)
print(res_df.to_string(index=False))
res_df.to_csv(RES / "residualization_TCGA.tsv", sep="\t", index=False)
scores.to_csv(RES / "scores_TCGA.tsv", sep="\t")

# ============================================================
# 3. GSE286332 — same correlation analysis
# ============================================================
print("\n=== 3. GSE286332: same 8-gene vs HLA-II correlation ===")
gse_path = "/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz"
df = pd.read_csv(gse_path, sep="\t", low_memory=False)
fpkm_cols = [c for c in df.columns if c.endswith("_FPKM")]
fpkm = df.groupby("Gene_Symbol", as_index=True)[fpkm_cols].sum()
fpkm.columns = [c.replace("_FPKM", "") for c in fpkm.columns]
log_fpkm = np.log2(fpkm + 1.0)

g8_g, _ = zmean(log_fpkm, GENE_8)
h1_g, _ = zmean(log_fpkm, HLA_I)
h2_g, _ = zmean(log_fpkm, HLA_II)
imm_g, _ = zmean(log_fpkm, IMMUNE_PROXY)

gse_scores = pd.DataFrame({"g8_RAI": g8_g, "HLA_I": h1_g, "HLA_II": h2_g, "immune": imm_g})
gse_scores["group"] = ["PTC" if s.startswith("NG_") else "PTC_HT" for s in gse_scores.index]

print(f"  GSE286332 samples: {len(gse_scores)}")
print("\n  Spearman ρ:")
print(gse_scores[["g8_RAI", "HLA_I", "HLA_II", "immune"]].corr(method="spearman").round(3))

# Per-group residualization for PTC+HT effect
yg = (gse_scores["group"] == "PTC_HT").astype(int).values
d_raw_g, p_raw_g = cohen_compare(gse_scores["g8_RAI"].values, yg)
# Note: in GSE286332 PTC_HT = 1 (Hashimoto-overlap)
# So d should be NEGATIVE since PTC+HT shows lower RAI
X = sm.add_constant(gse_scores[["HLA_II"]])
fit = sm.OLS(gse_scores["g8_RAI"], X).fit()
g8_resid = gse_scores["g8_RAI"].values - fit.fittedvalues.values
d_r_g, p_r_g = cohen_compare(g8_resid, yg)

# HLA-II raw vs residualized for 8-gene
d_h_raw_g, p_h_raw_g = cohen_compare(gse_scores["HLA_II"].values, yg)
X = sm.add_constant(gse_scores[["g8_RAI"]])
fit2 = sm.OLS(gse_scores["HLA_II"], X).fit()
hla_resid = gse_scores["HLA_II"].values - fit2.fittedvalues.values
d_h_r_g, p_h_r_g = cohen_compare(hla_resid, yg)

gse_resid_rows = [
    dict(model="raw 8-gene RAI (PTC_HT vs PTC)", d=round(d_raw_g,3), mw_p=p_raw_g),
    dict(model="8-gene | HLA-II residualized", d=round(d_r_g,3), mw_p=p_r_g),
    dict(model="raw HLA-II", d=round(d_h_raw_g,3), mw_p=p_h_raw_g),
    dict(model="HLA-II | 8-gene residualized", d=round(d_h_r_g,3), mw_p=p_h_r_g),
]
gse_res_df = pd.DataFrame(gse_resid_rows)
print(gse_res_df.to_string(index=False))
gse_res_df.to_csv(RES / "residualization_GSE286332.tsv", sep="\t", index=False)
gse_scores.to_csv(RES / "scores_GSE286332.tsv", sep="\t")

# ============================================================
# 4. Cross-cohort meta — 3-cohort effect on Hashimoto-like vs PTC
# ============================================================
print("\n=== 4. 3-cohort meta — 8-gene effect on Hashimoto-like vs PTC ===")
# TCGA: use Q12 Hashimoto-like sub-cluster vs rest
# Try to load Q12 labels
q12_files = [
    PROJ / "results/q12_hashimoto/hashimoto_like_calls.tsv",
    PROJ / "results/v17_realfix/Q12_hashimoto_calls.tsv",
    PROJ / "results/v17_realfix/Hashimoto_classified.tsv",
]
tcga_meta = None
for f in q12_files:
    if f.exists():
        tcga_meta = pd.read_csv(f, sep="\t")
        print(f"  Q12 file found: {f}")
        break
if tcga_meta is None:
    # Use DM1 as proxy for "differentiated" and DM2 for "Hashimoto-like-included"
    print("  Q12 file not found — using DM1 vs DM2 as TCGA proxy effect")
    a = scores.loc[scores["DM"]=="DM2", "g8_RAI"].values  # Hashimoto-included
    b = scores.loc[scores["DM"]=="DM1", "g8_RAI"].values  # differentiated
    sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    tcga_d = (a.mean() - b.mean()) / max(sp, 1e-9)
    tcga_n_case = len(a); tcga_n_ctrl = len(b)
    tcga_label = "TCGA-THCA (DM2 vs DM1, proxy)"
else:
    # Use real Q12 calls
    tcga_meta = tcga_meta.set_index("sample_id") if "sample_id" in tcga_meta.columns else tcga_meta
    common = [s for s in scores.index if s in tcga_meta.index]
    if "hashimoto_like" in tcga_meta.columns:
        is_hashi = tcga_meta.loc[common, "hashimoto_like"].astype(bool)
        a = scores.loc[common, "g8_RAI"].loc[is_hashi].values
        b = scores.loc[common, "g8_RAI"].loc[~is_hashi].values
        sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
        tcga_d = (a.mean() - b.mean()) / max(sp, 1e-9)
        tcga_n_case = len(a); tcga_n_ctrl = len(b)
        tcga_label = "TCGA-THCA (Q12 Hashimoto-like vs rest)"
    else:
        tcga_d = None

# K2 cohort: load K2 8-gene scores + Hashimoto calls
k2_paths = [PROJ / "results/v17_korean/K2_korean_predictions_v4.tsv",
            PROJ / "results/v17_korean/K2_korean_predictions_v3.tsv"]
k2_d = None; k2_n_case = k2_n_ctrl = 0; k2_label = "K2 (PRJEB11591)"
for f in k2_paths:
    if f.exists():
        k2 = pd.read_csv(f, sep="\t")
        print(f"  K2 prediction file: {f}, cols={k2.columns.tolist()[:10]}")
        # Try to find a hashimoto-like column or DM call
        if "DM_call" in k2.columns and "g8_score" in k2.columns:
            is_dm2 = k2["DM_call"].str.startswith("DM2")
            a = k2.loc[is_dm2, "g8_score"].values
            b = k2.loc[~is_dm2, "g8_score"].values
            if len(a) > 5 and len(b) > 5:
                sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
                k2_d = (a.mean() - b.mean()) / max(sp, 1e-9)
                k2_n_case = len(a); k2_n_ctrl = len(b)
                k2_label = "K2 (DM2 vs DM1)"
        break

# GSE286332
gse_d = d_raw_g
gse_n_case = int(yg.sum()); gse_n_ctrl = int((1 - yg).sum())
gse_label = "GSE286332 (PTC+HT vs PTC)"

meta_rows = []
if tcga_d is not None:
    meta_rows.append(dict(cohort=tcga_label, n_case=tcga_n_case, n_ctrl=tcga_n_ctrl, cohen_d=round(tcga_d, 3)))
if k2_d is not None:
    meta_rows.append(dict(cohort=k2_label, n_case=k2_n_case, n_ctrl=k2_n_ctrl, cohen_d=round(k2_d, 3)))
meta_rows.append(dict(cohort=gse_label, n_case=gse_n_case, n_ctrl=gse_n_ctrl, cohen_d=round(gse_d, 3)))
meta_df = pd.DataFrame(meta_rows)
print("\n  3-cohort meta (8-gene RAI score effect):")
print(meta_df.to_string(index=False))
meta_df.to_csv(RES / "meta_3cohort.tsv", sep="\t", index=False)

# Random-effects pooled (DerSimonian-Laird approximation)
def pool_d(rows):
    """Random-effects pooled Cohen's d (approx)."""
    if len(rows) < 2:
        return None
    ds = np.array([r["cohen_d"] for r in rows])
    ns = np.array([r["n_case"] + r["n_ctrl"] for r in rows])
    # var(d) ≈ (n1+n2)/(n1*n2) + d²/(2*(n1+n2))
    n1 = np.array([r["n_case"] for r in rows])
    n2 = np.array([r["n_ctrl"] for r in rows])
    v_d = (n1+n2)/(n1*n2) + (ds**2)/(2*(n1+n2))
    w = 1.0 / v_d
    d_fix = (w * ds).sum() / w.sum()
    Q = (w * (ds - d_fix) ** 2).sum()
    df_h = len(ds) - 1
    tau2 = max(0.0, (Q - df_h) / (w.sum() - (w**2).sum() / w.sum()))
    w_re = 1.0 / (v_d + tau2)
    d_re = (w_re * ds).sum() / w_re.sum()
    var_re = 1.0 / w_re.sum()
    return dict(d_pooled=round(d_re, 3),
                se=round(np.sqrt(var_re), 3),
                ci_lo=round(d_re - 1.96 * np.sqrt(var_re), 3),
                ci_hi=round(d_re + 1.96 * np.sqrt(var_re), 3),
                tau2=round(tau2, 4), Q=round(Q, 3), df=df_h)

pooled = pool_d(meta_rows)
if pooled is not None:
    print(f"\n  Random-effects pooled d: {pooled['d_pooled']:+.3f} [95% CI {pooled['ci_lo']:+.3f}, {pooled['ci_hi']:+.3f}]")
    print(f"  tau²={pooled['tau2']:.4f}, Cochran's Q={pooled['Q']:.3f} (df={pooled['df']})")

# ============================================================
# 5. Summary
# ============================================================
summary = {
    "TCGA": {
        "n_samples": len(samps),
        "spearman_rho_g8_vs_HLA2": round(stats.spearmanr(scores["g8_RAI"], scores["HLA_II"]).correlation, 3),
        "spearman_rho_g8_vs_immune": round(stats.spearmanr(scores["g8_RAI"], scores["immune"]).correlation, 3),
        "residualization": residual_rows,
    },
    "GSE286332": {
        "n_samples": len(gse_scores),
        "spearman_rho_g8_vs_HLA2": round(stats.spearmanr(gse_scores["g8_RAI"], gse_scores["HLA_II"]).correlation, 3),
        "spearman_rho_g8_vs_immune": round(stats.spearmanr(gse_scores["g8_RAI"], gse_scores["immune"]).correlation, 3),
        "residualization": gse_resid_rows,
    },
    "meta_3cohort": meta_rows,
    "pooled": pooled,
}
(RES / "P5_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ All outputs saved to {RES}")
