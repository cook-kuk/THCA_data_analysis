#!/usr/bin/env python3
"""
R9 — Pan-immune methylation deep-dive within BRAF-cPTC stratum.

Goal: close reviewer Q "is DM1 > DM2 just 8-gene methylation, or is there
broader immune-gene promoter de-methylation?"

Disk audit (2026-05-09):
  - 8-gene HM450 (TCGA-THCA, n=503): on disk
        /data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv
  - Full TCGA-THCA HM450 (~485k probe x ~510 sample matrix): NOT on disk
  - GSE97466 HM450: only top-5000 variable probes on disk (no full beta matrix)
  - HM450 manifest (GPL13534): NOT on disk — cannot map probe -> gene without network

Pivot: per-task spec section (6), fall back to gene-expression-based
"indirect promoter activity" proxy (RNA z-score; promoter de-methylation
typically correlates with RNA up). We:

  (A) Re-run 8-gene HM450 DM1 vs DM2 within BRAF-cPTC (positive control;
      existing finding).
  (B) Build RNA-proxy promoter-activity scores for an immune-gene panel
      (checkpoint, Treg, effector, antigen presentation, B cell,
      cytokine/chemokine).
  (C) Composite immune-promoter activity score per sample, with sign
      inverted to mimic "de-methylation" semantics.
  (D) Test composite separation DM1 vs DM2 INDEPENDENT of 8-gene mean beta
      (residualize 8-gene panel out, then re-test).
  (E) Direction/concordance check: for each immune gene, RNA z direction
      consistent with DM1 hypothesis (effectors/antigen UP, FOXP3/IDO1 may
      go either way).
  (F) Causal layering: regress HT-axis surrogate ~ 8-gene + immune-proxy +
      MAPK output to test mediation channels.

Output: project/results/p2_braf_nature_sprint_2026_05_09/r9_hm450_immune/
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/p2_braf_nature_sprint_2026_05_09/r9_hm450_immune"
OUT.mkdir(parents=True, exist_ok=True)

CLIN = Path(
    "/data/thca/repo_results/p2_braf_nature_sprint_2026_05_09/h6_survival/h6_merged_clinical.tsv"
)
MASTER = Path("/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv")
METH8 = Path("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv")
TCGA_Z = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv")
GSE97466_TOP5K = Path("/data/thca/data_processed/methylation/GSE97466_beta_top5000.tsv")


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    n1, n2 = len(a), len(b)
    s1, s2 = a.var(ddof=1), b.var(ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    if sp == 0 or np.isnan(sp):
        return float("nan")
    return (a.mean() - b.mean()) / sp


def wilc_p(a, b) -> float:
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def bh_fdr(pvec: np.ndarray) -> np.ndarray:
    p = np.asarray(pvec, float)
    n = len(p)
    out = np.full(n, np.nan)
    mask = ~np.isnan(p)
    pm = p[mask]; m = len(pm)
    if m == 0:
        return out
    order = np.argsort(pm)
    ranked = pm[order]
    fdr = ranked * m / (np.arange(1, m + 1))
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    res = np.empty(m)
    res[order] = fdr
    out[mask] = res
    return out


# ---------------------------------------------------------------------------
# load: clinical + 8-gene HM450 + RNA-z, build BRAF-cPTC stratum
# ---------------------------------------------------------------------------
print("=" * 78)
print("R9 — pan-immune methylation deep-dive within BRAF-cPTC")
print("=" * 78)

clin = pd.read_csv(CLIN, sep="\t")
print(f"  clinical n={len(clin)}; cols={len(clin.columns)}")

meth8 = pd.read_csv(METH8, sep="\t")
print(f"  8-gene HM450: n={len(meth8)}, cols={list(meth8.columns)}")

# clin sample_short = "TCGA-XX-XXXX-01"; meth8 sample_short = "TCGA-XX-XXXX" (12-char)
# Use 12-char patient barcode as join key.
clin["sample_short_match"] = clin["sample_short"].astype(str).str.slice(0, 12)
meth8["sample_short_match"] = meth8["sample_short"].astype(str).str.slice(0, 12)

mer = clin.merge(meth8, on="sample_short_match", how="left", suffixes=("", "_meth"))
print(f"  merged clin x 8-gene HM450: {len(mer)} (with HM450: {mer['mean_8g_beta'].notna().sum()})")

# stratum
braf_cptc_full = mer[(mer["molecular_subtype"] == "BRAF_like") & (mer["histology_subtype"] == "cPTC")].copy()
braf_cptc = braf_cptc_full[braf_cptc_full["dm"].isin(["DM1", "DM2"])].copy()
print(f"  BRAF-cPTC stratum: total={len(braf_cptc_full)} | DM1+DM2={len(braf_cptc)} "
      f"(DM1={(braf_cptc.dm == 'DM1').sum()}, DM2={(braf_cptc.dm == 'DM2').sum()})")

# only those with HM450 8-gene
braf_meth = braf_cptc.dropna(subset=["mean_8g_beta"]).copy()
print(f"  BRAF-cPTC with HM450 8-gene panel: n={len(braf_meth)} "
      f"(DM1={(braf_meth.dm == 'DM1').sum()}, DM2={(braf_meth.dm == 'DM2').sum()})")

# ---------------------------------------------------------------------------
# (A) 8-gene HM450 — DM1 vs DM2 within BRAF-cPTC (positive control)
# ---------------------------------------------------------------------------
print("\n[A] 8-gene HM450 DM1 vs DM2, BRAF-cPTC (positive control)")
genes_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR", "mean_8g_beta"]
rows_8 = []
for g in genes_8:
    if g not in braf_meth.columns:
        continue
    a = braf_meth.loc[braf_meth.dm == "DM1", g].values
    b = braf_meth.loc[braf_meth.dm == "DM2", g].values
    d = cohen_d(a, b)
    p = wilc_p(a, b)
    rows_8.append({
        "layer": "HM450_beta_8gene",
        "gene": g,
        "n_DM1": len(a) - np.isnan(a).sum(),
        "n_DM2": len(b) - np.isnan(b).sum(),
        "mean_DM1": float(np.nanmean(a)) if len(a) else np.nan,
        "mean_DM2": float(np.nanmean(b)) if len(b) else np.nan,
        "cohen_d_DM1_minus_DM2": d,
        "wilcoxon_p": p,
    })
df_8 = pd.DataFrame(rows_8)
df_8["fdr_bh"] = bh_fdr(df_8["wilcoxon_p"].values)
df_8.to_csv(OUT / "r9_8gene_HM450_DM1_vs_DM2.tsv", sep="\t", index=False)
print(df_8.to_string(index=False))

# ---------------------------------------------------------------------------
# (B) RNA-proxy promoter activity — immune gene panel
#     RNA up <=> promoter de-methylation (proxy; cannot distinguish promoter
#     vs enhancer activity from RNA alone — honest caveat)
# ---------------------------------------------------------------------------
print("\n[B] RNA-proxy promoter activity — immune gene panel")
panel: dict[str, list[str]] = {
    "checkpoint_inh": ["CTLA4", "PDCD1", "CD274", "HAVCR2", "LAG3", "TIGIT", "IDO1"],
    "treg_supp":      ["FOXP3", "IL10"],
    "effector":       ["GZMA", "GZMB", "PRF1", "IFNG", "NKG7"],
    "ag_pres_I":      ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2"],
    "ag_pres_II":     ["HLA-DRA", "HLA-DRB1", "HLA-DPB1", "CIITA"],
    "b_cell":         ["CD79A", "CD79B", "MS4A1", "AICDA"],
    "chemokine":      ["CXCL9", "CXCL10", "CXCL13", "CCL5"],
}
all_genes = sorted({g for v in panel.values() for g in v})
gene_to_module = {g: m for m, lst in panel.items() for g in lst}

# load RNA-z, subset to panel rows
print(f"  reading TCGA RNA-z subset for {len(all_genes)} immune genes ...")
rna = pd.read_csv(TCGA_Z, sep="\t")
rna = rna.rename(columns={rna.columns[0]: "gene_symbol"})
rna_p = rna[rna["gene_symbol"].isin(all_genes)].copy()
missing = sorted(set(all_genes) - set(rna_p["gene_symbol"]))
print(f"  panel genes found in RNA-z: {len(rna_p)}/{len(all_genes)}; missing={missing}")

# transpose to sample x gene; index = TCGA sample_id (e.g. TCGA-DJ-A2Q6-01A)
rna_p = rna_p.drop_duplicates("gene_symbol").set_index("gene_symbol").T.reset_index().rename(columns={"index": "sample_id"})
print(f"  RNA-z samples: {len(rna_p)}")

# braf_cptc has dm + sample_id (TCGA-XX-XXXX-01A); merge
braf_rna = braf_cptc.merge(rna_p, on="sample_id", how="left")
present_genes = [g for g in all_genes if g in braf_rna.columns]
print(f"  BRAF-cPTC merged with RNA: {len(braf_rna)}; "
      f"with z available: {braf_rna[present_genes].notna().any(axis=1).sum() if present_genes else 0}")

# also merge 8-gene HM450 onto this (for residualization)
braf_rna = braf_rna.merge(meth8[["sample_short_match", "mean_8g_beta"]].rename(
    columns={"mean_8g_beta": "mean_8g_beta_meth"}),
    on="sample_short_match", how="left", suffixes=("", "_meth8"))

# per-gene RNA DM1 vs DM2
rows_rna = []
for g in present_genes:
    a = braf_rna.loc[braf_rna.dm == "DM1", g].values
    b = braf_rna.loc[braf_rna.dm == "DM2", g].values
    d = cohen_d(a, b)
    p = wilc_p(a, b)
    rows_rna.append({
        "layer": "RNA_z_proxy",
        "module": gene_to_module.get(g, ""),
        "gene": g,
        "n_DM1": int(np.sum(~np.isnan(a))),
        "n_DM2": int(np.sum(~np.isnan(b))),
        "mean_z_DM1": float(np.nanmean(a)) if np.sum(~np.isnan(a)) else np.nan,
        "mean_z_DM2": float(np.nanmean(b)) if np.sum(~np.isnan(b)) else np.nan,
        "cohen_d_DM1_minus_DM2_RNAz": d,
        "wilcoxon_p": p,
        # convention: promoter activity proxy = +RNA z; "de-methylation proxy" = +RNA z
        "interp_promoter_demeth_in_DM1": "yes" if (not np.isnan(d) and d > 0) else ("no" if not np.isnan(d) else "na"),
    })
df_rna = pd.DataFrame(rows_rna)
df_rna["fdr_bh"] = bh_fdr(df_rna["wilcoxon_p"].values)
df_rna = df_rna.sort_values("cohen_d_DM1_minus_DM2_RNAz", ascending=False)
df_rna.to_csv(OUT / "r9_per_gene_methylation_d.tsv", sep="\t", index=False)
print(f"  per-gene RNA proxy table: {len(df_rna)} genes; "
      f"sig FDR<0.05: {(df_rna['fdr_bh'] < 0.05).sum()}; "
      f"DM1-up (d>0): {(df_rna['cohen_d_DM1_minus_DM2_RNAz'] > 0).sum()}")

# ---------------------------------------------------------------------------
# (C) Composite immune-promoter activity score
#     up-direction expected ("hot in DM1"): effectors, antigen, chemokines
#     bidirectional / not signed: IDO1, FOXP3, checkpoints — keep raw mean
# ---------------------------------------------------------------------------
print("\n[C] Composite immune-promoter activity scores")
modules_up = ["effector", "ag_pres_I", "ag_pres_II", "b_cell", "chemokine"]
modules_amb = ["checkpoint_inh", "treg_supp"]

def module_mean(df: pd.DataFrame, mod: str) -> pd.Series:
    gs = [g for g in panel[mod] if g in df.columns]
    return df[gs].mean(axis=1) if gs else pd.Series(np.nan, index=df.index)

for m in panel:
    braf_rna[f"mod_{m}_z"] = module_mean(braf_rna, m)

braf_rna["immune_demeth_proxy_up"] = braf_rna[[f"mod_{m}_z" for m in modules_up]].mean(axis=1)
braf_rna["immune_proxy_amb"] = braf_rna[[f"mod_{m}_z" for m in modules_amb]].mean(axis=1)
braf_rna["immune_demeth_proxy_full"] = braf_rna[[f"mod_{m}_z" for m in panel]].mean(axis=1)

# per-sample export
keep_cols = ["sample_id", "sample_short", "dm", "molecular_subtype", "histology_subtype",
             "mean_8g_beta_meth", "rai_score_v17", "tds16_score_v17"] + \
            [f"mod_{m}_z" for m in panel] + \
            ["immune_demeth_proxy_up", "immune_proxy_amb", "immune_demeth_proxy_full"]
keep_cols = [c for c in keep_cols if c in braf_rna.columns]
braf_rna[keep_cols].to_csv(OUT / "r9_immune_promoter_composite.tsv", sep="\t", index=False)

# DM1 vs DM2 on composites
rows_comp = []
for col in ["immune_demeth_proxy_up", "immune_proxy_amb", "immune_demeth_proxy_full",
            "mean_8g_beta_meth"] + [f"mod_{m}_z" for m in panel]:
    if col not in braf_rna.columns:
        continue
    a = braf_rna.loc[braf_rna.dm == "DM1", col].values
    b = braf_rna.loc[braf_rna.dm == "DM2", col].values
    rows_comp.append({
        "feature": col,
        "n_DM1": int(np.sum(~np.isnan(a))),
        "n_DM2": int(np.sum(~np.isnan(b))),
        "mean_DM1": float(np.nanmean(a)) if np.sum(~np.isnan(a)) else np.nan,
        "mean_DM2": float(np.nanmean(b)) if np.sum(~np.isnan(b)) else np.nan,
        "cohen_d_DM1_minus_DM2": cohen_d(a, b),
        "wilcoxon_p": wilc_p(a, b),
    })
df_comp = pd.DataFrame(rows_comp)
df_comp["fdr_bh"] = bh_fdr(df_comp["wilcoxon_p"].values)
df_comp.to_csv(OUT / "r9_composite_DM1_vs_DM2.tsv", sep="\t", index=False)
print(df_comp.to_string(index=False))

# ---------------------------------------------------------------------------
# (D) Independence from 8-gene mean beta — residualize, re-test
#     If immune-demeth proxy survives after regressing out mean_8g_beta_meth,
#     it adds a layer beyond thyroid-differentiation methylation.
# ---------------------------------------------------------------------------
print("\n[D] Independence from 8-gene mean beta")
rows_indep = []
have_8g = braf_rna.dropna(subset=["mean_8g_beta_meth"]).copy()
print(f"  BRAF-cPTC with both RNA + 8-gene HM450: n={len(have_8g)} "
      f"(DM1={(have_8g.dm == 'DM1').sum()}, DM2={(have_8g.dm == 'DM2').sum()})")
from numpy.linalg import lstsq

def residualize(y: pd.Series, x: pd.Series) -> pd.Series:
    m = (~y.isna()) & (~x.isna())
    if m.sum() < 3:
        return pd.Series(np.nan, index=y.index)
    X = np.column_stack([np.ones(m.sum()), x[m].values])
    beta, *_ = lstsq(X, y[m].values, rcond=None)
    pred = X @ beta
    res = pd.Series(np.nan, index=y.index)
    res[m] = y[m].values - pred
    return res

for col in ["immune_demeth_proxy_up", "immune_proxy_amb", "immune_demeth_proxy_full"] + \
           [f"mod_{m}_z" for m in panel]:
    if col not in have_8g.columns:
        continue
    raw = have_8g[col]
    res = residualize(raw, have_8g["mean_8g_beta_meth"])
    have_8g[f"{col}__resid8g"] = res
    a = res[have_8g.dm == "DM1"].values
    b = res[have_8g.dm == "DM2"].values
    a_raw = raw[have_8g.dm == "DM1"].values
    b_raw = raw[have_8g.dm == "DM2"].values
    rows_indep.append({
        "feature": col,
        "raw_d": cohen_d(a_raw, b_raw),
        "raw_p": wilc_p(a_raw, b_raw),
        "resid_d_8gMethBeta_out": cohen_d(a, b),
        "resid_p_8gMethBeta_out": wilc_p(a, b),
        "n_DM1": int(np.sum(~np.isnan(a))),
        "n_DM2": int(np.sum(~np.isnan(b))),
    })
df_indep = pd.DataFrame(rows_indep)
df_indep["raw_fdr"] = bh_fdr(df_indep["raw_p"].values)
df_indep["resid_fdr"] = bh_fdr(df_indep["resid_p_8gMethBeta_out"].values)
df_indep.to_csv(OUT / "r9_independence_from_8gene_meth.tsv", sep="\t", index=False)
print(df_indep.to_string(index=False))

# ---------------------------------------------------------------------------
# (E) RNA <-> methylation concordance check (8-gene only — only layer with
#     real HM450 beta on disk)
# ---------------------------------------------------------------------------
print("\n[E] RNA z vs HM450 beta concordance (8-gene panel, BRAF-cPTC)")
rows_conc = []
rna_full = rna.set_index("gene_symbol")
gene_8_real = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
rna_full_T = rna_full.loc[[g for g in gene_8_real if g in rna_full.index]].T.reset_index().rename(
    columns={"index": "sample_id"})
braf_join = braf_meth.merge(rna_full_T, on="sample_id", how="left", suffixes=("_meth", "_rna"))
for g in gene_8_real:
    rna_col = f"{g}_rna" if f"{g}_rna" in braf_join.columns else g
    meth_col = f"{g}_meth" if f"{g}_meth" in braf_join.columns else g
    if rna_col not in braf_join.columns or meth_col not in braf_join.columns:
        continue
    x = pd.to_numeric(braf_join[meth_col], errors="coerce")
    y = pd.to_numeric(braf_join[rna_col], errors="coerce")
    m = (~x.isna()) & (~y.isna())
    if m.sum() < 5:
        continue
    r, p = stats.spearmanr(x[m], y[m])
    rows_conc.append({
        "gene": g,
        "n": int(m.sum()),
        "spearman_meth_vs_rna": float(r),
        "p": float(p),
        "interp_expected_negative": "yes" if r < 0 else "no",
    })
df_conc = pd.DataFrame(rows_conc)
df_conc.to_csv(OUT / "r9_methylation_rna_concordance.tsv", sep="\t", index=False)
print(df_conc.to_string(index=False))

# ---------------------------------------------------------------------------
# (F) Causal layering — does immune-promoter proxy mediate MAPK -> HT axis?
#     Use:
#       MAPK output: tds16_score_v17 NEGATIVELY tracks MAPK (TDS = thyroid
#         differentiation; lower = more dediff/MAPK-driven downstream).
#         We use -tds16 as a MAPK output proxy; sanity check below.
#       HT axis: rai_score_v17 (RAI-likeness of immune/HT signature in master).
#       Mediator candidates: mean_8g_beta_meth (thyroid-diff methylation),
#                            immune_demeth_proxy_up (immune-promoter activity).
#     We compare three nested OLS:
#       M0:  rai ~ MAPK_proxy
#       M1:  rai ~ MAPK_proxy + mean_8g_beta
#       M2:  rai ~ MAPK_proxy + mean_8g_beta + immune_demeth_proxy_up
#     and report partial R^2 for each mediator.
# ---------------------------------------------------------------------------
print("\n[F] Causal layering — partial R^2 of methylation layers on HT axis (rai_score)")
mediation = have_8g.dropna(subset=["mean_8g_beta_meth", "immune_demeth_proxy_up",
                                    "tds16_score_v17", "rai_score_v17"]).copy()
print(f"  mediation n={len(mediation)} BRAF-cPTC with 8g HM450 + RNA + axes")

def ols_r2(y: np.ndarray, X: np.ndarray) -> float:
    if len(y) < X.shape[1] + 2:
        return float("nan")
    X1 = np.column_stack([np.ones(len(y)), X])
    beta, *_ = lstsq(X1, y, rcond=None)
    pred = X1 @ beta
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

rai = mediation["rai_score_v17"].values
mapk = -mediation["tds16_score_v17"].values
m8g = mediation["mean_8g_beta_meth"].values
imm = mediation["immune_demeth_proxy_up"].values

r2_M0 = ols_r2(rai, mapk.reshape(-1, 1))
r2_M1 = ols_r2(rai, np.column_stack([mapk, m8g]))
r2_M2 = ols_r2(rai, np.column_stack([mapk, m8g, imm]))
r2_M_imm_only = ols_r2(rai, np.column_stack([mapk, imm]))

mediation_rows = [
    {"model": "M0_rai~MAPK", "R2": r2_M0, "delta": np.nan, "mediator_added": ""},
    {"model": "M1_rai~MAPK+8gMeth", "R2": r2_M1, "delta": r2_M1 - r2_M0, "mediator_added": "mean_8g_beta_meth"},
    {"model": "M_imm_only_rai~MAPK+immProxy", "R2": r2_M_imm_only,
     "delta": r2_M_imm_only - r2_M0, "mediator_added": "immune_demeth_proxy_up"},
    {"model": "M2_rai~MAPK+8gMeth+immProxy", "R2": r2_M2,
     "delta": r2_M2 - r2_M1, "mediator_added": "immune_demeth_proxy_up_after_8gMeth"},
]
df_med = pd.DataFrame(mediation_rows)
df_med.to_csv(OUT / "r9_causal_layering.tsv", sep="\t", index=False)
print(df_med.to_string(index=False))

# ---------------------------------------------------------------------------
# (G) GSE97466 top-5000 sanity: any panel-immune probe gene names parseable
#     from probe_id alone? No (HM450 probe IDs are cgXXXXXXXX not gene names),
#     and no GPL13534 manifest on disk -> we cannot map to immune genes here.
#     Documented honestly in report.
# ---------------------------------------------------------------------------
gse_n_probes = sum(1 for _ in open(GSE97466_TOP5K)) - 1
print(f"\n[G] GSE97466 top-5000 sanity: {gse_n_probes} probes; no GPL13534 manifest "
      "on disk -> probe->gene mapping for immune panel skipped (network not used).")

# ---------------------------------------------------------------------------
# Summary JSON
# ---------------------------------------------------------------------------
sig_rna = df_rna.query("fdr_bh < 0.05")
demeth_DM1_count = int((df_rna["cohen_d_DM1_minus_DM2_RNAz"] > 0).sum())
demeth_DM1_sig = int(((df_rna["cohen_d_DM1_minus_DM2_RNAz"] > 0) & (df_rna["fdr_bh"] < 0.05)).sum())

summary = {
    "stratum": "BRAF-cPTC (TCGA-THCA, molecular_subtype=BRAF_like, histology=cPTC)",
    "n_total": int(len(braf_cptc_full)),
    "n_DM1_DM2_contrast": int(len(braf_cptc)),
    "n_DM1_DM2_with_HM450_8gene": int(len(braf_meth)),
    "n_DM1_DM2_with_RNA": int(len(braf_rna)),
    "n_DM1_DM2_with_both": int(len(have_8g)),
    "data_audit": {
        "full_HM450_TCGA": "NOT_ON_DISK",
        "8gene_HM450_TCGA": str(METH8),
        "GSE97466_full_HM450": "TOP_5000_ONLY (no GPL13534 manifest on disk)",
        "RNA_proxy_used": True,
    },
    "panel_genes_total": len(all_genes),
    "panel_genes_RNA_found": len(present_genes),
    "panel_RNA_DM1_up_count": int((df_rna["cohen_d_DM1_minus_DM2_RNAz"] > 0).sum()),
    "panel_RNA_DM1_up_FDR05_count": demeth_DM1_sig,
    "composite_immune_demeth_proxy_up_DM1_d": float(df_comp.set_index("feature").loc["immune_demeth_proxy_up", "cohen_d_DM1_minus_DM2"]) if "immune_demeth_proxy_up" in df_comp["feature"].values else None,
    "composite_resid8g_d": float(df_indep.set_index("feature").loc["immune_demeth_proxy_up", "resid_d_8gMethBeta_out"]) if "immune_demeth_proxy_up" in df_indep["feature"].values else None,
    "mediation_R2_M0_MAPK_only": r2_M0,
    "mediation_R2_M1_plus_8gMeth": r2_M1,
    "mediation_R2_M2_plus_immProxy": r2_M2,
    "mediation_delta_immProxy_after_8gMeth": r2_M2 - r2_M1,
    "concordance_8gene_meth_vs_RNA_median_r": float(np.nanmedian(df_conc["spearman_meth_vs_rna"].values)) if len(df_conc) else None,
    "concordance_8gene_negative_count": int((df_conc["spearman_meth_vs_rna"] < 0).sum()) if len(df_conc) else 0,
}
with open(OUT / "r9_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print("\nWrote summary:", OUT / "r9_summary.json")

# Top 10 RNA-proxy hits up in DM1
print("\nTop 10 immune RNA-proxy hits UP in DM1 (BRAF-cPTC):")
print(df_rna.head(10)[["module", "gene", "cohen_d_DM1_minus_DM2_RNAz", "wilcoxon_p", "fdr_bh"]].to_string(index=False))

print("\nDONE.")
