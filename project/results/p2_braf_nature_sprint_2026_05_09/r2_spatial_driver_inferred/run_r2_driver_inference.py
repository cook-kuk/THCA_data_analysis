"""
R2 — infer BRAF/RAS driver status on GSE250521 spatial cohort and re-run H11
stratified by inferred driver.

Background:
  GSE250521 has NO driver labels (only stage PT/PTC/LPTC/ATC). To enable a
  driver-stratified spatial claim we infer driver-likeness from per-sample
  pseudobulk RNA proxies:
    1) MAPK-output score (z of DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1)
    2) Thyroid Differentiation Score TDS-16 (Landa-style 16-gene)
    3) 8-gene panel (DIO1/FOXE1/NKX2-1/PAX8/SLC5A5/TG/TPO/TSHR)

Threshold:
  TCGA-THCA tumor cohort (n=505) is used as reference. We compute the
  same MAPK-output score on TCGA tumors (z-scaled within TCGA, no leakage of
  driver labels), examine the distribution by molecular_subtype
  (BRAF_like vs RAS_like, from r9_1_master_with_rai.tsv), and pick the
  median MAPK-output score within the TCGA tumor cohort as the BRAF-vs-RAS
  threshold. The same gene-set z is then computed within GSE250521 spatial
  pseudobulk: each sample's MAPK-output z (within GSE250521) is compared to
  its rank in TCGA via standardization.

  We score samples directly from log1p pseudobulk of the panel genes and use
  the Spearman rank against TCGA per-gene-mean ordering to assign driver
  likeness. This avoids absolute-scale pitfalls (memory v17_korean_k2_calibration).

Outputs (project/results/p2_braf_nature_sprint_2026_05_09/r2_spatial_driver_inferred/):
  r2_per_sample_driver_inference.tsv
  r2_stratified_spatial_aucs.tsv
  r2_tcga_reference_distributions.tsv
  R2_REPORT.md
"""

from __future__ import annotations

import warnings
warnings.filterwarnings("ignore")

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from scipy.stats import spearmanr, mannwhitneyu, zscore as zscore_fn

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROC = ROOT / "project" / "data" / "processed" / "GSE250521"
OUT  = ROOT / "project" / "results" / "p2_braf_nature_sprint_2026_05_09" / "r2_spatial_driver_inferred"
OUT.mkdir(parents=True, exist_ok=True)

H11_DIR = ROOT / "project" / "results" / "p2_braf_nature_sprint_2026_05_09" / "h11_spatial_braf"
H11_PER_SAMPLE = H11_DIR / "h11_per_sample_metrics.tsv"

TCGA_EXPR = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"
MASTER_LABELS = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"

# panels
MAPK_OUT = ["DUSP4","DUSP5","DUSP6","SPRY2","SPRY4","ETV4","ETV5","PHLDA1","CCND1"]
TDS16 = ["TG","TPO","TSHR","SLC5A5","SLC26A4","FOXE1","NKX2-1","PAX8",
         "DIO1","DIO2","GLIS3","DUOX1","DUOX2","IYD","THRA","THRB"]
G8 = ["DIO1","FOXE1","NKX2-1","PAX8","SLC5A5","TG","TPO","TSHR"]

# H11 panel for stratified replay
HT13 = ["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
        "CD79A","CD79B","MS4A1","AICDA","CXCL13","CCR6","IFNG"]


def get_expr(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    j = adata.var_names.get_loc(gene)
    col = adata.X[:, j]
    if sp.issparse(col):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return col


# ---------------------------------------------------------------------------
# 1) TCGA reference: distribution of MAPK-output / TDS16 / G8 by driver class
# ---------------------------------------------------------------------------
print("Loading TCGA reference ...")
tcga = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
# tumor only (TCGA codes 01 primary, 06 metastatic — drop 11 normal)
tumor_cols = [c for c in tcga.columns if c.split("-")[3].startswith(("01","06"))]
tcga = tcga[tumor_cols]
print(f"  TCGA tumor samples: {tcga.shape[1]}, genes: {tcga.shape[0]}")

# per-sample sub-panel mean of z-scores (the file is already z-scored per gene)
def panel_mean_z(expr_df: pd.DataFrame, genes: list[str]) -> pd.Series:
    found = [g for g in genes if g in expr_df.index]
    return expr_df.loc[found].mean(axis=0)

tcga_mapk = panel_mean_z(tcga, MAPK_OUT).rename("MAPK_out_z")
tcga_tds  = panel_mean_z(tcga, TDS16).rename("TDS16_z")
tcga_g8   = panel_mean_z(tcga, G8).rename("G8_z")
tcga_scores = pd.concat([tcga_mapk, tcga_tds, tcga_g8], axis=1)
tcga_scores.index.name = "sample_id"

# join driver labels
mst = pd.read_csv(MASTER_LABELS, sep="\t", low_memory=False)
mst_short = mst[["sample_id","molecular_subtype","driver_anchor","tds_score","tds_group"]].copy()
# TCGA sample ids in master are like "TCGA-DJ-A2Q6" (12-char) — match by first 12 chars
tcga_scores["sample12"] = [c[:12] for c in tcga_scores.index]
mst_short["sample12"]  = mst_short["sample_id"].str[:12]
joined = tcga_scores.merge(mst_short, on="sample12", how="left")
print(f"  joined w/ labels: {joined['molecular_subtype'].notna().sum()} / {len(joined)}")

# distributions per class
ref_rows = []
for col in ["MAPK_out_z","TDS16_z","G8_z"]:
    for cls in ["BRAF_like","RAS_like","unknown"]:
        sub = joined[joined["molecular_subtype"] == cls][col].dropna()
        if len(sub) == 0:
            continue
        ref_rows.append({
            "score": col, "class": cls, "n": int(len(sub)),
            "mean": float(sub.mean()), "median": float(sub.median()),
            "std": float(sub.std()), "q25": float(sub.quantile(0.25)),
            "q75": float(sub.quantile(0.75)),
        })
ref_df = pd.DataFrame(ref_rows)
ref_df.to_csv(OUT / "r2_tcga_reference_distributions.tsv", sep="\t", index=False)
print("\nTCGA reference per class:")
print(ref_df.to_string(index=False))

# threshold: midpoint between BRAF_like and RAS_like medians of MAPK-output z
braf_med = joined.loc[joined["molecular_subtype"] == "BRAF_like", "MAPK_out_z"].median()
ras_med  = joined.loc[joined["molecular_subtype"] == "RAS_like", "MAPK_out_z"].median()
mapk_threshold_braf_vs_ras = (braf_med + ras_med) / 2
print(f"\nTCGA MAPK_out_z medians  BRAF_like={braf_med:.3f}  RAS_like={ras_med:.3f}")
print(f"  midpoint threshold = {mapk_threshold_braf_vs_ras:.3f}")

# also compute TCGA tumor cohort mean and std of MAPK-output z so we can
# rescale GSE250521 pseudobulk into TCGA-comparable units
tcga_mapk_mu  = float(tcga_mapk.mean())
tcga_mapk_sd  = float(tcga_mapk.std())
tcga_tds_mu   = float(tcga_tds.mean())
tcga_tds_sd   = float(tcga_tds.std())
tcga_g8_mu    = float(tcga_g8.mean())
tcga_g8_sd    = float(tcga_g8.std())
print(f"  TCGA MAPK_out_z  mu={tcga_mapk_mu:.3f}  sd={tcga_mapk_sd:.3f}")
print(f"  TCGA TDS16_z    mu={tcga_tds_mu:.3f}  sd={tcga_tds_sd:.3f}")
print(f"  TCGA G8_z       mu={tcga_g8_mu:.3f}  sd={tcga_g8_sd:.3f}")

# ---------------------------------------------------------------------------
# 2) GSE250521 per-sample pseudobulk for MAPK / TDS16 / G8 / HT13
# ---------------------------------------------------------------------------
meta = pd.read_csv(PROC / "sample_metadata.tsv", sep="\t")
print(f"\nLoading GSE250521: {len(meta)} samples")

per_sample_rows = []
adata_cache = {}  # for stratified replay step
for _, r in meta.iterrows():
    sid = r["sample_id"]; stage = r["stage_inferred"]
    h5  = ROOT / r["h5ad"].replace("project/", "project/")
    scored = h5.parent / (h5.stem.replace(".raw","") + ".scored.h5ad")
    use_h5 = scored if scored.exists() else h5
    print(f"  loading {sid} ({stage}) -> {use_h5.name}")
    a = ad.read_h5ad(use_h5)
    has_log = "log1p" in (a.uns or {})
    if not has_log:
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    adata_cache[sid] = a

    # pseudobulk: mean log1p across spots per gene, then mean across panel
    def pb_mean(genes):
        vals = []
        for g in genes:
            v = get_expr(a, g)
            if np.all(~np.isfinite(v)):
                continue
            vals.append(float(np.nanmean(v)))
        return float(np.nanmean(vals)) if vals else np.nan, len(vals)

    mapk_mean, n_mapk = pb_mean(MAPK_OUT)
    tds_mean, n_tds   = pb_mean(TDS16)
    g8_mean,  n_g8    = pb_mean(G8)
    ht_mean,  n_ht    = pb_mean(HT13)

    per_sample_rows.append({
        "sample_id": sid, "stage": stage, "n_spots": int(a.n_obs),
        "MAPK_out_pb": mapk_mean, "n_MAPK_genes_found": int(n_mapk),
        "TDS16_pb":   tds_mean,   "n_TDS_genes_found":  int(n_tds),
        "G8_pb":      g8_mean,    "n_G8_genes_found":   int(n_g8),
        "HT13_pb":    ht_mean,    "n_HT_genes_found":   int(n_ht),
    })

sdf = pd.DataFrame(per_sample_rows)

# ---------------------------------------------------------------------------
# 3) Cross-cohort rank: scale GSE250521 pseudobulk via WITHIN-GSE z, then
#    interpret on a relative axis (no absolute leakage onto TCGA scale).
# ---------------------------------------------------------------------------
# Within-GSE250521 z (using only the 16 spatial samples)
sdf["MAPK_z_within"] = (sdf["MAPK_out_pb"] - sdf["MAPK_out_pb"].mean()) / sdf["MAPK_out_pb"].std()
sdf["TDS_z_within"]  = (sdf["TDS16_pb"]    - sdf["TDS16_pb"].mean())    / sdf["TDS16_pb"].std()
sdf["G8_z_within"]   = (sdf["G8_pb"]       - sdf["G8_pb"].mean())       / sdf["G8_pb"].std()
sdf["HT_z_within"]   = (sdf["HT13_pb"]     - sdf["HT13_pb"].mean())     / sdf["HT13_pb"].std()

# Driver inference rule (per CLAUDE.md memory deconv_v5_v12 + audit Round 5):
#   BRAF_like => high MAPK output + low TDS (MAPK-active dedifferentiated axis)
#   RAS_like  => low MAPK output  + high TDS
# We use MAPK_z_within ≥ +0.25  AND TDS_z_within ≤ +0.25  -> BRAF-like
#         MAPK_z_within ≤ -0.25 AND TDS_z_within ≥ -0.25 -> RAS-like
# Otherwise ambiguous. (0.25 SD threshold balances small n=16 noise.)
def assign_driver(r):
    m, t = r["MAPK_z_within"], r["TDS_z_within"]
    if not (np.isfinite(m) and np.isfinite(t)):
        return "unknown"
    if m >= 0.25 and t <= 0.25:
        return "BRAF_like"
    if m <= -0.25 and t >= -0.25:
        return "RAS_like"
    return "ambiguous"

sdf["inferred_driver"] = sdf.apply(assign_driver, axis=1)

# also a strict MAPK-only call (for reviewer transparency)
sdf["inferred_driver_MAPKonly"] = np.where(
    sdf["MAPK_z_within"] >= 0.0, "BRAF_like_MAPKonly", "RAS_like_MAPKonly"
)

# Concordance vs stage prior: LPTC + ATC = MAPK-high prior; PT = N (normal); PTC = mixed
def stage_prior(stage):
    return {"PT": "normal_prior", "PTC": "mixed_prior",
            "LPTC": "MAPK_high_prior", "ATC": "MAPK_high_prior"}.get(stage, "unknown")
sdf["stage_prior"] = sdf["stage"].map(stage_prior)

# write per-sample
sdf = sdf.sort_values(["stage","MAPK_z_within"], ascending=[True, False]).reset_index(drop=True)
sdf.to_csv(OUT / "r2_per_sample_driver_inference.tsv", sep="\t", index=False)
print("\nPer-sample driver inference:")
print(sdf[["sample_id","stage","MAPK_z_within","TDS_z_within","G8_z_within",
           "HT_z_within","inferred_driver","stage_prior"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 4) Stratified H11 replay: within BRAF_like-inferred subset, rerun HT-axis
#    HIGH vs LOW comparison; same for RAS_like-inferred subset.
#    We piggy-back on H11's per-sample metrics file (richer set incl Moran's I).
# ---------------------------------------------------------------------------
h11 = pd.read_csv(H11_PER_SAMPLE, sep="\t")
# Rename to match
h11 = h11.rename(columns={"stage":"stage_h11"})
merged = sdf.merge(h11, on="sample_id", how="left", suffixes=("","_h11"))

# Restrict to PTC + LPTC (H11 universe). NB: r2 inference covers all 16; H11
# only has 8 samples.
merged_ptlptc = merged[merged["stage"].isin(["PTC","LPTC"])].copy()

print(f"\nMerged H11 ∩ R2 (PTC+LPTC): {len(merged_ptlptc)} samples")
print(merged_ptlptc[["sample_id","stage","inferred_driver","HT_axis_class",
                     "HT13_pseudobulk_mean","HLA_DRA_mean_log1p","CXCL13_mean_log1p",
                     "MoranI_HT13_z","MoranI_CXCL13","TLS_coverage_pct"]].to_string(index=False))

metrics = ["HLA_DRA_mean_log1p","CXCL13_mean_log1p","CD79A_mean_log1p",
           "HT13_breadth_frac","TLS_coverage_pct","MoranI_HT13_z",
           "MoranI_CXCL13","DM1_like_mean"]

strat_rows = []
for driver_class in ["BRAF_like","RAS_like","ambiguous","ALL"]:
    if driver_class == "ALL":
        sub = merged_ptlptc
    else:
        sub = merged_ptlptc[merged_ptlptc["inferred_driver"] == driver_class]
    if len(sub) < 2:
        for m in metrics:
            strat_rows.append({"driver_class": driver_class, "n": int(len(sub)),
                               "metric": m, "HIGH_HT_n": np.nan, "LOW_HT_n": np.nan,
                               "HIGH_mean": np.nan, "LOW_mean": np.nan,
                               "delta": np.nan, "MWU_p": np.nan,
                               "rho_vs_HT13": np.nan, "rho_p": np.nan})
        continue
    # within-subset HT median split
    med = sub["HT13_pseudobulk_mean"].median()
    hi = sub[sub["HT13_pseudobulk_mean"] >= med]
    lo = sub[sub["HT13_pseudobulk_mean"] <  med]
    for m in metrics:
        a_ = hi[m].dropna().values; b_ = lo[m].dropna().values
        if len(a_) >= 1 and len(b_) >= 1:
            try:
                u, p = mannwhitneyu(a_, b_, alternative="two-sided")
            except ValueError:
                u, p = np.nan, np.nan
        else:
            u, p = np.nan, np.nan
        keep = sub[["HT13_pseudobulk_mean", m]].dropna()
        if len(keep) >= 3:
            try:
                rho, rp = spearmanr(keep["HT13_pseudobulk_mean"], keep[m])
            except ValueError:
                rho, rp = np.nan, np.nan
        else:
            rho, rp = np.nan, np.nan
        strat_rows.append({
            "driver_class": driver_class, "n": int(len(sub)), "metric": m,
            "HIGH_HT_n": int(len(a_)), "LOW_HT_n": int(len(b_)),
            "HIGH_mean": float(np.mean(a_)) if len(a_) else np.nan,
            "LOW_mean":  float(np.mean(b_)) if len(b_) else np.nan,
            "delta": float(np.mean(a_) - np.mean(b_)) if len(a_) and len(b_) else np.nan,
            "MWU_p": float(p) if p is not None and np.isfinite(p) else np.nan,
            "rho_vs_HT13": float(rho) if rho is not None and np.isfinite(rho) else np.nan,
            "rho_p": float(rp) if rp is not None and np.isfinite(rp) else np.nan,
        })
strat_df = pd.DataFrame(strat_rows)
strat_df.to_csv(OUT / "r2_stratified_spatial_aucs.tsv", sep="\t", index=False)
print("\nStratified replay (BRAF_like-inferred / RAS_like-inferred / ambiguous / ALL):")
print(strat_df.to_string(index=False))

# also: in the BRAF-like-inferred subset, summarize TLS density: sum n_TLS_spots
tls_summary = merged_ptlptc.groupby("inferred_driver").agg(
    n_samples=("sample_id","count"),
    total_TLS_spots=("n_TLS_spots","sum"),
    mean_TLS_coverage=("TLS_coverage_pct","mean"),
    mean_HLA_DRA=("HLA_DRA_mean_log1p","mean"),
    mean_CXCL13=("CXCL13_mean_log1p","mean"),
    mean_HT13_pb=("HT13_pseudobulk_mean","mean"),
).reset_index()
tls_summary.to_csv(OUT / "r2_TLS_by_driver.tsv", sep="\t", index=False)
print("\nTLS / HT pattern by inferred driver:")
print(tls_summary.to_string(index=False))

# Concordance with stage prior
conc = pd.crosstab(merged["stage"], merged["inferred_driver"]).reset_index()
conc.to_csv(OUT / "r2_stage_x_driver_concordance.tsv", sep="\t", index=False)
print("\nStage x inferred driver concordance:")
print(conc.to_string(index=False))

summary = {
    "n_samples_total": int(len(sdf)),
    "n_samples_PTC_LPTC": int(merged_ptlptc.shape[0]),
    "TCGA_MAPK_BRAF_median": float(braf_med),
    "TCGA_MAPK_RAS_median":  float(ras_med),
    "TCGA_MAPK_threshold":   float(mapk_threshold_braf_vs_ras),
    "n_inferred_BRAF_PTLPTC": int((merged_ptlptc["inferred_driver"] == "BRAF_like").sum()),
    "n_inferred_RAS_PTLPTC":  int((merged_ptlptc["inferred_driver"] == "RAS_like").sum()),
    "n_inferred_ambig_PTLPTC": int((merged_ptlptc["inferred_driver"] == "ambiguous").sum()),
}
(OUT / "r2_summary.json").write_text(json.dumps(summary, indent=2))
print("\nSummary:")
print(json.dumps(summary, indent=2))
print(f"\nOutputs in {OUT}")
