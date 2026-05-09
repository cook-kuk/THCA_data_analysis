"""
H11 — spatial validation of HT/TLS axis in PTC / LPTC (GSE250521 Visium n=16).

GSE250521 has NO BRAF/RAS/HT annotation in metadata — only stage (N/PTC/LPTC/ATC).
We therefore restrict to PTC (n=4 primary) + LPTC (n=4 lateral metastasis) = 8 samples
and stratify by HT-13-axis bimodality on spatial pseudobulk (HIGH vs LOW HT median split).

Per-sample metrics:
    - HT-13 panel score (mean log1p z-scored across 13 genes)
    - DM1_like_score mean (already in scored.h5ad obs)
    - TLS-domain coverage = % spots with high CXCL13 + high HLA-DR + high CD79A
    - LR diversity (CXCL13<->CCR6, CXCL13<->CCR7) co-expression Spearman per sample
    - Moran's I of HT-axis score (k-nearest spatial weights, k=6)

Compare HIGH-HT vs LOW-HT samples within PTC+LPTC by Mann-Whitney.
Outputs:
    h11_per_sample_metrics.tsv
    h11_spatial_HTaxis_correlation.tsv
    H11_REPORT.md
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
from scipy.stats import spearmanr, mannwhitneyu, zscore
from libpysal.weights import KNN
from esda.moran import Moran

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROC = ROOT / "project" / "data" / "processed" / "GSE250521"
OUT  = ROOT / "project" / "results" / "p2_braf_nature_sprint_2026_05_09" / "h11_spatial_braf"
OUT.mkdir(parents=True, exist_ok=True)

HT13 = ["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
        "CD79A","CD79B","MS4A1","AICDA","CXCL13","CCR6","IFNG"]
TLS_DOMAIN_GENES = {"CXCL13": 0.50, "HLA-DRA": 0.50, "CD79A": 0.50}
# spot is "TLS-like" if all three genes exceed sample-internal q=0.50 (median) AND
# at least one gene exceeds sample-internal q=0.85 (high outlier).

# stage filter
SAMPLES_TARGET = ["PTC", "LPTC"]


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


def ht13_score_per_spot(adata: ad.AnnData) -> np.ndarray:
    """Mean of z-scored log1p expression across 13 panel genes (within-sample z).

    The within-sample z preserves spatial heterogeneity of the panel even when
    panel mean differs across samples; for cross-sample comparisons we use
    HT13_pseudobulk (mean log1p) and HT13_q90 (within-sample 90th percentile),
    both retained in the per-sample table.
    """
    cols = []
    for g in HT13:
        v = get_expr(adata, g)
        if np.all(~np.isfinite(v)) or np.nanstd(v) == 0:
            cols.append(np.zeros_like(v))
            continue
        mu = np.nanmean(v); sd = np.nanstd(v)
        z = (v - mu) / (sd + 1e-9)
        cols.append(z)
    arr = np.column_stack(cols)
    return arr.mean(axis=1)


def ht13_pseudobulk(adata: ad.AnnData) -> dict:
    """Cross-sample-comparable HT-13 panel scores.

    HT13_pseudobulk_mean: mean of (mean log1p expression per gene across spots),
                          averaged across the 13 genes — single scalar per sample.
    HT13_pseudobulk_pos:  fraction of spots where >= 6 of 13 genes are detected (>0).
    HT13_max_gene_d_q90 : 90th-percentile of within-sample mean-log1p expression.
    """
    g_means = []
    for g in HT13:
        v = get_expr(adata, g)
        if np.all(~np.isfinite(v)):
            g_means.append(np.nan); continue
        g_means.append(float(np.nanmean(v)))
    arr = np.array(g_means, dtype=float)
    panel_mean = float(np.nanmean(arr))
    # detection breadth per spot
    det = np.zeros(adata.n_obs, dtype=int)
    for g in HT13:
        v = get_expr(adata, g)
        det = det + (v > 0).astype(int)
    frac_breadth = float(np.mean(det >= 6))  # ≥6/13 panel genes detected at the spot
    return {
        "HT13_pseudobulk_mean": panel_mean,
        "HT13_breadth_frac": frac_breadth,
        "HT13_per_gene_means": arr.tolist(),
    }


def tls_domain_mask(adata: ad.AnnData) -> np.ndarray:
    """Spot is TLS-like if CXCL13 + HLA-DRA + CD79A all > sample-median AND
    at least one of (CXCL13, CD79A) > q85 within sample."""
    cxcl13 = get_expr(adata, "CXCL13")
    hladra = get_expr(adata, "HLA-DRA")
    cd79a  = get_expr(adata, "CD79A")
    g_ok   = (cxcl13 > np.nanmedian(cxcl13)) & (hladra > np.nanmedian(hladra)) & (cd79a > np.nanmedian(cd79a))
    high_outlier = (cxcl13 > np.nanquantile(cxcl13, 0.85)) | (cd79a > np.nanquantile(cd79a, 0.85))
    return (g_ok & high_outlier).astype(int)


def lr_coexpr(adata: ad.AnnData, lig: str, rec: str) -> float:
    """Spearman correlation of ligand and receptor over spots within sample.
    Positive => co-localized; near 0 => uncorrelated; negative => anti-correlated.
    Uses spots where either gene > 0 (sparser background)."""
    a = get_expr(adata, lig); b = get_expr(adata, rec)
    if np.all(~np.isfinite(a)) or np.all(~np.isfinite(b)):
        return np.nan
    keep = (a > 0) | (b > 0)
    if keep.sum() < 30:
        return np.nan
    r, _ = spearmanr(a[keep], b[keep], nan_policy="omit")
    return float(r) if r is not None and np.isfinite(r) else np.nan


def morans_i_score(adata: ad.AnnData, score: np.ndarray, k: int = 6) -> tuple[float, float]:
    """Moran's I via libpysal/esda using k-NN on (array_row, array_col) (Visium grid)."""
    coords = np.column_stack([
        adata.obs["array_row"].astype(float).to_numpy(),
        adata.obs["array_col"].astype(float).to_numpy(),
    ])
    finite = np.isfinite(score) & np.isfinite(coords).all(axis=1)
    if finite.sum() < 100:
        return np.nan, np.nan
    coords = coords[finite]; score = score[finite]
    try:
        w = KNN.from_array(coords, k=min(k, len(coords)-1))
        w.transform = "r"
        m = Moran(score, w, permutations=99)
        return float(m.I), float(m.p_sim)
    except Exception:
        return np.nan, np.nan


# ---------------------------------------------------------------------------
# 1) load metadata & target samples (PTC + LPTC)
# ---------------------------------------------------------------------------
meta = pd.read_csv(PROC / "sample_metadata.tsv", sep="\t")
meta = meta[meta["stage_inferred"].isin(SAMPLES_TARGET)].reset_index(drop=True)
print(f"Target samples (PTC+LPTC): {len(meta)}")

rows = []
spot_dump_for_corr = []
for _, r in meta.iterrows():
    sid = r["sample_id"]; stage = r["stage_inferred"]
    h5  = ROOT / r["h5ad"].replace("project/", "project/")
    # use scored.h5ad if available
    scored = h5.parent / (h5.stem.replace(".raw","") + ".scored.h5ad")
    use_h5 = scored if scored.exists() else h5
    print(f"  loading {sid} stage={stage}: {use_h5.name}")
    a = ad.read_h5ad(use_h5)
    # ensure log1p (scored has log1p applied)
    has_log = "log1p" in (a.uns or {})
    if not has_log:
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)

    ht13_z = ht13_score_per_spot(a)
    a.obs["HT13_z_within"] = ht13_z
    pb = ht13_pseudobulk(a)
    tls = tls_domain_mask(a)
    a.obs["TLS_like"] = tls

    # raw mean log1p of CXCL13 + CD79A as a TLS-strength scalar
    cxcl13_mean = float(np.nanmean(get_expr(a, "CXCL13")))
    cd79a_mean  = float(np.nanmean(get_expr(a, "CD79A")))
    hladra_mean = float(np.nanmean(get_expr(a, "HLA-DRA")))

    dm1_mean = float(np.nanmean(a.obs["DM1_like_score"])) if "DM1_like_score" in a.obs else np.nan
    rai_mean = float(np.nanmean(a.obs["RAI_8_score"])) if "RAI_8_score" in a.obs else np.nan
    epi_mean = float(np.nanmean(a.obs["Epithelial_score"])) if "Epithelial_score" in a.obs else np.nan

    cx_ccr6 = lr_coexpr(a, "CXCL13", "CCR6")
    cx_ccr7 = lr_coexpr(a, "CXCL13", "CCR7")
    cd79a_ms4a1 = lr_coexpr(a, "CD79A", "MS4A1")  # B-cell co-localization sanity

    # spatial autocorrelation of within-sample z-scored panel
    moran_I, moran_p = morans_i_score(a, ht13_z, k=6)
    # spatial autocorrelation of CXCL13 expression (raw log1p)
    moran_I_cxcl13, moran_p_cxcl13 = morans_i_score(a, get_expr(a, "CXCL13"), k=6)
    # cluster vs diffuse: Moran's I of TLS_like (binary), only if any TLS spots
    if tls.sum() >= 5:
        moran_I_tls, moran_p_tls = morans_i_score(a, tls.astype(float), k=6)
    else:
        moran_I_tls, moran_p_tls = np.nan, np.nan

    rows.append({
        "sample_id": sid,
        "stage": stage,
        "n_spots": int(a.n_obs),
        "HT13_pseudobulk_mean":  pb["HT13_pseudobulk_mean"],   # cross-sample comparable
        "HT13_breadth_frac":     pb["HT13_breadth_frac"],
        "CXCL13_mean_log1p":     cxcl13_mean,
        "CD79A_mean_log1p":      cd79a_mean,
        "HLA_DRA_mean_log1p":    hladra_mean,
        "DM1_like_mean": dm1_mean,
        "RAI_8_mean": rai_mean,
        "Epithelial_mean": epi_mean,
        "TLS_coverage_pct": float(100 * np.nanmean(tls)),
        "n_TLS_spots": int(np.nansum(tls)),
        "LR_CXCL13_CCR6_rho": cx_ccr6,
        "LR_CXCL13_CCR7_rho": cx_ccr7,
        "LR_CD79A_MS4A1_rho": cd79a_ms4a1,
        "MoranI_HT13_z": moran_I,
        "MoranI_HT13_z_p": moran_p,
        "MoranI_CXCL13": moran_I_cxcl13,
        "MoranI_CXCL13_p": moran_p_cxcl13,
        "MoranI_TLS_like": moran_I_tls,
        "MoranI_TLS_p": moran_p_tls,
    })

per_sample = pd.DataFrame(rows)
per_sample = per_sample.sort_values(["stage","HT13_pseudobulk_mean"], ascending=[True, False]).reset_index(drop=True)

# ---------------------------------------------------------------------------
# 2) HT-axis median split within PTC+LPTC pool (cross-sample comparable score)
# ---------------------------------------------------------------------------
median_ht = per_sample["HT13_pseudobulk_mean"].median()
per_sample["HT_axis_class"] = np.where(per_sample["HT13_pseudobulk_mean"] >= median_ht, "HIGH_HT", "LOW_HT")
# DM1 surrogate: bimodal split on DM1_like_mean (correlated proxy).
# Memory: TCGA HM450 mean_8g_beta DM1 < DM2 d=-1.75; sub-A MAPK-active reaches panel silencing.
# Higher DM1_like_score = MORE DM1 phenotype = BRAF/MAPK-active cluster in this scoring.
median_dm1 = per_sample["DM1_like_mean"].median()
per_sample["DM_surrogate"] = np.where(per_sample["DM1_like_mean"] >= median_dm1, "DM1_like", "DM2_like")

per_sample.to_csv(OUT / "h11_per_sample_metrics.tsv", sep="\t", index=False)
print("\nPer-sample table:")
print(per_sample[["sample_id","stage","HT13_pseudobulk_mean","CXCL13_mean_log1p","DM1_like_mean",
                  "TLS_coverage_pct","MoranI_HT13_z","MoranI_TLS_like","HT_axis_class","DM_surrogate"]].to_string(index=False))

# ---------------------------------------------------------------------------
# 3) HIGH vs LOW HT-axis spatial-pattern stats (Mann-Whitney)
# ---------------------------------------------------------------------------
metrics_to_test = [
    "TLS_coverage_pct",
    "n_TLS_spots",
    "CXCL13_mean_log1p",
    "CD79A_mean_log1p",
    "HLA_DRA_mean_log1p",
    "HT13_breadth_frac",
    "LR_CXCL13_CCR6_rho",
    "LR_CXCL13_CCR7_rho",
    "LR_CD79A_MS4A1_rho",
    "MoranI_HT13_z",
    "MoranI_CXCL13",
    "MoranI_TLS_like",
    "DM1_like_mean",
]
hi = per_sample[per_sample["HT_axis_class"] == "HIGH_HT"]
lo = per_sample[per_sample["HT_axis_class"] == "LOW_HT"]
corr_rows = []
for m in metrics_to_test:
    a_ = hi[m].dropna().values
    b_ = lo[m].dropna().values
    if len(a_) < 2 or len(b_) < 2:
        u, p = np.nan, np.nan
    else:
        try:
            u, p = mannwhitneyu(a_, b_, alternative="two-sided")
        except ValueError:
            u, p = np.nan, np.nan
    corr_rows.append({
        "metric": m,
        "HIGH_HT_n": int(len(a_)),
        "HIGH_HT_mean": float(np.mean(a_)) if len(a_) else np.nan,
        "HIGH_HT_median": float(np.median(a_)) if len(a_) else np.nan,
        "LOW_HT_n": int(len(b_)),
        "LOW_HT_mean": float(np.mean(b_)) if len(b_) else np.nan,
        "LOW_HT_median": float(np.median(b_)) if len(b_) else np.nan,
        "delta_HIGH_minus_LOW_mean": float(np.mean(a_) - np.mean(b_)) if len(a_) and len(b_) else np.nan,
        "MannWhitney_U": float(u) if u is not None and np.isfinite(u) else np.nan,
        "MannWhitney_p": float(p) if p is not None and np.isfinite(p) else np.nan,
    })
corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(OUT / "h11_spatial_HTaxis_correlation.tsv", sep="\t", index=False)
print("\nHIGH vs LOW HT-axis comparison:")
print(corr_df.to_string(index=False))

# also: Spearman correlation of HT13_pseudobulk_mean with each spatial metric across n=8 samples
sp_rows = []
for m in metrics_to_test:
    keep = per_sample[["HT13_pseudobulk_mean", m]].dropna()
    if len(keep) < 4:
        rho, p = np.nan, np.nan
    else:
        rho, p = spearmanr(keep["HT13_pseudobulk_mean"], keep[m])
    sp_rows.append({"metric": m, "n": int(len(keep)),
                    "Spearman_rho_vs_HT13": float(rho) if rho is not None and np.isfinite(rho) else np.nan,
                    "p": float(p) if p is not None and np.isfinite(p) else np.nan})
sp_df = pd.DataFrame(sp_rows)
sp_df.to_csv(OUT / "h11_spearman_vs_HT13.tsv", sep="\t", index=False)
print("\nSpearman with HT-13 axis (n=8 samples):")
print(sp_df.to_string(index=False))

summary = {
    "n_samples": int(len(per_sample)),
    "stages_present": sorted(per_sample["stage"].unique().tolist()),
    "HT13_median_split": float(median_ht),
    "DM1_median_split": float(median_dm1),
    "HIGH_HT_samples": hi["sample_id"].tolist(),
    "LOW_HT_samples":  lo["sample_id"].tolist(),
    "TLS_coverage_HIGH_mean": float(hi["TLS_coverage_pct"].mean()),
    "TLS_coverage_LOW_mean":  float(lo["TLS_coverage_pct"].mean()),
    "MoranI_HT13_z_HIGH_mean":  float(hi["MoranI_HT13_z"].mean()),
    "MoranI_HT13_z_LOW_mean":   float(lo["MoranI_HT13_z"].mean()),
    "MoranI_CXCL13_HIGH_mean":  float(hi["MoranI_CXCL13"].mean()),
    "MoranI_CXCL13_LOW_mean":   float(lo["MoranI_CXCL13"].mean()),
}
(OUT / "h11_summary.json").write_text(json.dumps(summary, indent=2))
print("\nSummary:", json.dumps(summary, indent=2))
print(f"\nOutputs in {OUT}")
