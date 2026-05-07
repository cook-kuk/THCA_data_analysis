#!/usr/bin/env python3
"""
Spatial transcriptomics full-package strengthening (A-F).

Datasets (already scored, on disk):
  A. GSE250521 (16 Visium slides; 4 PT + 4 PTC + 4 LPTC + 4 ATC) — `project/results/01_spatial_score/`
  B. GSE230424 (4 PTC+HT slides) — `project_external_st/results/scores/`
  C. GSE248205 (2 CONTROL + 3 HT + 3 GD; no cancer) — `project_external_st/results/scores/`

6 analyses:
  A. Per-spot DM1 × niche (Epithelial × Proliferation) 2D density
  B. Spot-level driver-orthogonal (DM1 spot scatter vs sample-level driver class) — sample-level proxy
  C. GSE230424 PTC+HT TLS / B-cell / IGHV niche scoring + Moran's I (h5ad-based)
  D. GSE248205 HT/GD vs CONTROL spot-level immune axis baseline
  E. TROP2 (TACSTD2) per-sample spot distribution + Moran's I niche test (h5ad-based)
  F. 3-dataset cross-cohort integration — all 28 samples × 8-gene + auxiliary axes

Outputs (project/results/spatial_full_2026_05_06/):
  spatial_full_summary.json
  spatial_A_DM1_niche_density.tsv
  spatial_B_driver_orthogonal_per_sample.tsv
  spatial_C_HT_TLS_spatial.tsv
  spatial_D_GSE248205_baseline.tsv
  spatial_E_TROP2_spot_distribution.tsv
  spatial_F_cross_cohort_integration.tsv

Figures (project/papers_hub_2026_05_04/assets/spatial_full/):
  S_F1_DM1_niche_2d_density.png      (analysis A)
  S_F2_driver_orthogonal_box.png     (analysis B)
  S_F3_HT_TLS_spatial_panel.png      (analysis C)
  S_F4_GSE248205_baseline.png        (analysis D)
  S_F5_TROP2_spot_distribution.png   (analysis E)
  S_F6_cross_cohort_integration.png  (analysis F)

Discipline: CPU-only. No new download. No causal claim. Paper 1 main story preserved.
TROP2 spatial diagnosis = honest "spot-level co-localisation" test, not enrichment claim.
"""
from __future__ import annotations
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, pearsonr, kendalltau
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize, LinearSegmentedColormap
import warnings
warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/results/spatial_full_2026_05_06"
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/spatial_full"
OUT.mkdir(parents=True, exist_ok=True)
ASSETS.mkdir(parents=True, exist_ok=True)

GSE250521_TSV = ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz"
EXTERNAL_SCORES_DIR = ROOT/"project_external_st/results/scores"
GSE250521_H5AD = ROOT/"project/data/processed/GSE250521"
EXTERNAL_H5AD = ROOT/"project_external_st/data/processed"

# Paper-hub color palette (light theme tone-matched)
COLORS = {
    "PT": "#3C6B4F", "PTC": "#34547A", "LPTC": "#B8893C", "ATC": "#7B1F2A",
    "PTC_HT": "#962E2E", "CONTROL": "#52525a", "HT": "#34547A", "GD": "#B8893C",
}

print("="*70)
print("SPATIAL FULL PACKAGE 2026-05-06 — start")
print("="*70)

# ============ load all per-spot data once ============
print("\n[load] GSE250521 all-spots...")
g521 = pd.read_csv(GSE250521_TSV, sep="\t")
print(f"  {len(g521)} spots, {g521['sample_id'].nunique()} samples")
print(f"  per-stage spots: {g521.groupby('stage').size().to_dict()}")

print("\n[load] external_st (GSE230424 + GSE248205)...")
ext_files = sorted(EXTERNAL_SCORES_DIR.glob("*_spot_scores.tsv.gz"))
print(f"  {len(ext_files)} files")
ext_dfs = []
for f in ext_files:
    df = pd.read_csv(f, sep="\t")
    ext_dfs.append(df)
ext = pd.concat(ext_dfs, axis=0, ignore_index=True)
print(f"  {len(ext)} spots, {ext['sample_id'].nunique()} samples")
print(f"  per-condition spots: {ext.groupby(['dataset','condition_inferred']).size().to_dict()}")

# Standardize column names — give g521 condition_inferred = stage
g521 = g521.rename(columns={"stage": "condition_inferred"})
g521["dataset"] = "GSE250521"
g521["disease_axis"] = "thyroid_cancer"

# Harmonize column subset for cross-cohort integration
HARMONIZED = ["sample_id","dataset","condition_inferred","disease_axis",
              "RAI_8_score","DM1_like_score","TDS_overlap_score" if "TDS_overlap_score" in ext.columns else "TDS_like_score",
              "Epithelial_score","Proliferation_score","total_counts","n_genes_by_counts","spot_id",
              "array_row","array_col"]

# Build harmonized cross-cohort table
g521_h = g521.copy()
g521_h = g521_h.rename(columns={"TDS_like_score":"TDS_overlap_score"})  # both refer to RAI∪NONOVERLAP combined
ext_h = ext.copy().rename(columns={
    "RAI_8_score_raw":"RAI_8_score","DM1_like_score_raw":"DM1_like_score",
    "TDS_overlap_score_raw":"TDS_overlap_score","Epithelial_score_raw":"Epithelial_score",
    "Proliferation_score_raw":"Proliferation_score","THYROID_NONOVERLAP_score_raw":"THYROID_NONOVERLAP_score",
})
COMMON = ["sample_id","dataset","condition_inferred","disease_axis","RAI_8_score","DM1_like_score",
          "TDS_overlap_score","Epithelial_score","Proliferation_score",
          "total_counts","n_genes_by_counts","spot_id","array_row","array_col"]
g521_h["disease_axis"] = "thyroid_cancer"
for c in COMMON:
    if c not in g521_h.columns: g521_h[c] = np.nan
    if c not in ext_h.columns: ext_h[c] = np.nan
all_spots = pd.concat([g521_h[COMMON], ext_h[COMMON]], axis=0, ignore_index=True)
print(f"\n[harmonized] {len(all_spots)} total spots across {all_spots['sample_id'].nunique()} samples")
print(f"  datasets: {all_spots['dataset'].value_counts().to_dict()}")

# ====================================================================
# A. Per-spot DM1 × niche (Epithelial × Proliferation) 2D density
# ====================================================================
print("\n" + "="*70)
print("ANALYSIS A — DM1 × niche 2D density (GSE250521 16-slide)")
print("="*70)

a_rows = []
for stage, sub in g521.groupby("condition_inferred"):
    # within stage: DM1 vs Epithelial + Proliferation
    rho_dm_epi, p_dm_epi = spearmanr(sub["DM1_like_score"], sub["Epithelial_score"], nan_policy="omit")
    rho_dm_prol, p_dm_prol = spearmanr(sub["DM1_like_score"], sub["Proliferation_score"], nan_policy="omit")
    rho_dm_tds, p_dm_tds = spearmanr(sub["DM1_like_score"], sub["TDS_like_score"] if "TDS_like_score" in sub else sub.get("TDS_overlap_score", pd.Series([np.nan]*len(sub))), nan_policy="omit")
    a_rows.append({
        "stage": stage, "n_spots": len(sub),
        "rho_DM1_vs_Epithelial": rho_dm_epi, "p_DM1_vs_Epi": p_dm_epi,
        "rho_DM1_vs_Proliferation": rho_dm_prol, "p_DM1_vs_Prol": p_dm_prol,
        "rho_DM1_vs_TDS": rho_dm_tds, "p_DM1_vs_TDS": p_dm_tds,
        "DM1_mean": sub["DM1_like_score"].mean(), "DM1_std": sub["DM1_like_score"].std(),
        "Epi_mean": sub["Epithelial_score"].mean(),
        "Prol_mean": sub["Proliferation_score"].mean(),
    })
df_A = pd.DataFrame(a_rows)
df_A.to_csv(OUT/"spatial_A_DM1_niche_density.tsv", sep="\t", index=False)
print(df_A.to_string(index=False))

# Figure A: 2D density per stage (4-panel) + scatter overlay
fig, axes = plt.subplots(1, 4, figsize=(15, 4), sharey=True)
stages = ["PT","PTC","LPTC","ATC"]
for ax, stage in zip(axes, stages):
    sub = g521[g521.condition_inferred == stage]
    x = sub["Epithelial_score"].values
    y = sub["DM1_like_score"].values
    h = ax.hexbin(x, y, gridsize=35, cmap="OrRd", mincnt=2, alpha=0.85)
    rho_a = df_A[df_A.stage == stage]["rho_DM1_vs_Epithelial"].iloc[0] if not df_A[df_A.stage == stage].empty else np.nan
    ax.set_title(f"{stage}  (n={len(sub):,} spots)\nρ(DM1, Epi) = {rho_a:.3f}", fontsize=10)
    ax.set_xlabel("Epithelial score (z)")
    if ax is axes[0]: ax.set_ylabel("DM1_like score (z)")
    ax.axhline(0, ls=":", lw=0.6, color="gray"); ax.axvline(0, ls=":", lw=0.6, color="gray")
fig.suptitle("S_F1  Per-spot DM1 × Epithelial niche density across stages — GSE250521 (16 slides)",
             fontsize=11, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F1_DM1_niche_2d_density.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ====================================================================
# B. Spot-level driver-orthogonal — sample-level driver-class proxy
# ====================================================================
# Approach: GSE250521 stages (PT/PTC/LPTC/ATC) used as driver-class proxy.
# Within each stage, characterize DM1 spot distribution: spread, extreme-tail fraction.
# Real driver mutation per spot is unobtainable from public ST; this is honest proxy.
print("\n" + "="*70)
print("ANALYSIS B — driver-orthogonal proxy (per-stage DM1 distribution)")
print("="*70)

b_rows = []
for sample, sub in g521.groupby("sample_id"):
    stage = sub["condition_inferred"].iloc[0]
    dm = sub["DM1_like_score"].dropna()
    rai = sub["RAI_8_score"].dropna()
    # extreme-tail
    high_dm = (dm > dm.quantile(0.75)).sum() / len(dm) if len(dm) else np.nan
    low_dm = (dm < dm.quantile(0.25)).sum() / len(dm) if len(dm) else np.nan
    b_rows.append({
        "sample_id": sample, "stage": stage, "n_spots": len(dm),
        "DM1_mean": dm.mean(), "DM1_std": dm.std(),
        "DM1_q75": dm.quantile(0.75) if len(dm) else np.nan,
        "DM1_q25": dm.quantile(0.25) if len(dm) else np.nan,
        "DM1_iqr": dm.quantile(0.75)-dm.quantile(0.25) if len(dm) else np.nan,
        "high_DM1_spot_frac": high_dm,
        "rho_DM1_RAI8_within_sample": spearmanr(dm, rai, nan_policy="omit")[0] if len(dm) > 5 else np.nan,
    })
df_B = pd.DataFrame(b_rows)
df_B.to_csv(OUT/"spatial_B_driver_orthogonal_per_sample.tsv", sep="\t", index=False)
print(df_B.to_string(index=False))

# Figure B: per-stage DM1 IQR + within-sample DM1↔RAI8 correlation
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
order = ["PT","PTC","LPTC","ATC"]
for i, stage in enumerate(order):
    sub = df_B[df_B.stage == stage]
    ax.scatter([i]*len(sub), sub["DM1_iqr"].values, s=80,
               color=COLORS.get(stage, "#666"), alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
ax.set_xticks(range(4)); ax.set_xticklabels(order)
ax.set_ylabel("Within-sample DM1 IQR (spread)")
ax.set_title("B1. Per-sample DM1 spread by stage\n(spread of within-sample DM1 across spots)", fontsize=10)
ax.grid(axis="y", alpha=0.3)

ax = axes[1]
for i, stage in enumerate(order):
    sub = df_B[df_B.stage == stage]
    ax.scatter([i]*len(sub), sub["rho_DM1_RAI8_within_sample"].values, s=80,
               color=COLORS.get(stage, "#666"), alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
ax.axhline(-1, ls=":", color="gray", lw=0.6)
ax.set_xticks(range(4)); ax.set_xticklabels(order)
ax.set_ylabel("ρ(DM1_like, RAI_8) within sample")
ax.set_title("B2. Within-sample DM1 ↔ RAI_8 mathematical identity check\n(should be exactly −1 if DM1 = −RAI_8)", fontsize=10)
ax.set_ylim(-1.05, 0.0)
ax.grid(axis="y", alpha=0.3)

fig.suptitle("S_F2  Spot-level driver-orthogonal proxy — DM1 distribution per sample × stage",
             fontsize=11, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F2_driver_orthogonal_box.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ====================================================================
# F. Cross-cohort integration (28 samples, 8-gene + axes)
# ====================================================================
# Run F BEFORE C/E (which are h5ad-based and slower). F uses the harmonized
# all_spots table.
print("\n" + "="*70)
print("ANALYSIS F — 3-dataset cross-cohort integration (28 samples)")
print("="*70)

f_rows = []
for ds in ["GSE250521","GSE230424","GSE248205"]:
    sub_ds = all_spots[all_spots.dataset == ds]
    for cond, sub in sub_ds.groupby("condition_inferred"):
        f_rows.append({
            "dataset": ds, "condition": cond,
            "n_samples": sub["sample_id"].nunique(), "n_spots": len(sub),
            "DM1_mean": sub["DM1_like_score"].mean(), "DM1_median": sub["DM1_like_score"].median(),
            "DM1_q25": sub["DM1_like_score"].quantile(0.25),
            "DM1_q75": sub["DM1_like_score"].quantile(0.75),
            "Epi_mean": sub["Epithelial_score"].mean(),
            "Prol_mean": sub["Proliferation_score"].mean(),
        })
df_F = pd.DataFrame(f_rows)
df_F.to_csv(OUT/"spatial_F_cross_cohort_integration.tsv", sep="\t", index=False)
print(df_F.to_string(index=False))

# Figure F: cross-cohort DM1 distribution comparison
fig, ax = plt.subplots(figsize=(13, 5.5))
groups = []
for ds in ["GSE250521","GSE230424","GSE248205"]:
    sub_ds = all_spots[all_spots.dataset == ds]
    for cond in sub_ds["condition_inferred"].dropna().unique():
        sub = sub_ds[sub_ds.condition_inferred == cond]
        groups.append((ds, cond, sub["DM1_like_score"].dropna().values, sub["sample_id"].nunique()))
positions = list(range(len(groups)))
data = [g[2] for g in groups]
bp = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True, showfliers=False)
for patch, (ds, cond, _, _) in zip(bp["boxes"], groups):
    patch.set_facecolor(COLORS.get(cond, "#999"))
    patch.set_alpha(0.7)
ax.set_xticks(positions)
ax.set_xticklabels([f"{ds}\n{cond}\n(n_samp={n})" for ds, cond, _, n in groups],
                   fontsize=8.5, rotation=0)
ax.axhline(0, ls=":", lw=0.6, color="gray")
ax.set_ylabel("DM1_like score (z, within-sample)")
ax.set_title("S_F6  Cross-cohort DM1 distribution — 28 samples × 3 datasets · all spots pooled per condition",
             fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F6_cross_cohort_integration.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ====================================================================
# D. GSE248205 HT/GD vs CONTROL baseline
# ====================================================================
print("\n" + "="*70)
print("ANALYSIS D — GSE248205 baseline (HT/GD vs CONTROL)")
print("="*70)

g248 = ext[ext.dataset == "GSE248205"]
d_rows = []
for cond, sub in g248.groupby("condition_inferred"):
    for axis in ["RAI_8_score_raw","DM1_like_score_raw","Epithelial_score_raw","Proliferation_score_raw","CAF_ECM_score_raw","Hypoxia_score_raw"]:
        d_rows.append({
            "condition": cond, "axis": axis,
            "n_samples": sub["sample_id"].nunique(), "n_spots": len(sub),
            "mean": sub[axis].mean(), "median": sub[axis].median(),
            "q25": sub[axis].quantile(0.25), "q75": sub[axis].quantile(0.75),
        })
df_D = pd.DataFrame(d_rows)

# Pairwise MW between conditions, per axis
mw_rows = []
conds = sorted(g248["condition_inferred"].dropna().unique())
for axis in ["RAI_8_score_raw","DM1_like_score_raw","Epithelial_score_raw","Proliferation_score_raw","CAF_ECM_score_raw"]:
    for i, c1 in enumerate(conds):
        for c2 in conds[i+1:]:
            x = g248[g248.condition_inferred == c1][axis].dropna().values
            y = g248[g248.condition_inferred == c2][axis].dropna().values
            if len(x) < 5 or len(y) < 5: continue
            u, p = mannwhitneyu(x, y, alternative="two-sided")
            mw_rows.append({"axis": axis, "comparison": f"{c1}_vs_{c2}",
                            "n1": len(x), "n2": len(y), "mean_diff": x.mean()-y.mean(),
                            "MW_U": u, "MW_p": p})
df_D_mw = pd.DataFrame(mw_rows)
df_D = pd.concat([df_D, df_D_mw], axis=0, ignore_index=True, sort=False)
df_D.to_csv(OUT/"spatial_D_GSE248205_baseline.tsv", sep="\t", index=False)
print("\n--- Pairwise MW (GSE248205) ---")
print(df_D_mw.to_string(index=False))

fig, ax = plt.subplots(figsize=(11, 5.5))
order_d = ["CONTROL","HT","GD"]
axes_d = ["RAI_8_score_raw","DM1_like_score_raw","Epithelial_score_raw","Proliferation_score_raw","CAF_ECM_score_raw","Hypoxia_score_raw"]
x_pos = np.arange(len(axes_d))
w = 0.25
for i, cond in enumerate(order_d):
    sub = g248[g248.condition_inferred == cond]
    means = [sub[a].mean() for a in axes_d]
    ax.bar(x_pos + (i-1)*w, means, w, label=f"{cond} (n_samp={sub['sample_id'].nunique()})",
           color=COLORS.get(cond, "#666"), alpha=0.85, edgecolor="black", linewidth=0.5)
ax.axhline(0, ls=":", lw=0.6, color="gray")
ax.set_xticks(x_pos); ax.set_xticklabels([a.replace("_score_raw","") for a in axes_d], fontsize=10)
ax.set_ylabel("Mean per-spot score (z, within-sample)")
ax.set_title("S_F4  GSE248205 baseline — autoimmune-only HT/GD vs CONTROL\n"
             "Negative control: cancer-side DM1/lineage axes 신호 없음 / 약함을 확인",
             fontsize=11, fontweight="bold")
ax.legend(loc="best", fontsize=9)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F4_GSE248205_baseline.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ====================================================================
# C. GSE230424 PTC+HT TLS / B-cell / IGHV niche (h5ad-based)
# E. TROP2 per-spot distribution + spatial autocorrelation
# ====================================================================
# These need h5ad raw counts to score additional gene sets.
print("\n" + "="*70)
print("ANALYSIS C + E — h5ad-based TLS / TROP2 niche scoring")
print("="*70)

import anndata as ad

# Gene sets for TLS / B-cell / TROP2
TROP2_GENES = ["TACSTD2"]
HLA_II_GENES = ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
B_CELL_GENES = ["CD19","MS4A1","CD79A","CD79B","SDC1","JCHAIN"]
TLS_GENES   = ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"]
IGHV_GENES  = ["AICDA","IGHM","IGHG1","IGHA1","IGKC","IGLC2"]

def score_h5ad(h5ad_path, gene_set):
    """Load h5ad, compute mean log1p of available gene_set per spot, return DataFrame with score column."""
    try:
        a = ad.read_h5ad(h5ad_path)
    except Exception as e:
        return None
    var_names = a.var_names.astype(str)
    avail = [g for g in gene_set if g in var_names.values]
    if len(avail) == 0:
        return None
    # Use raw if available else X
    X = a.raw.X if (a.raw is not None) else a.X
    # find indices in raw.var
    if a.raw is not None:
        rvar = a.raw.var_names.astype(str)
        idx = [list(rvar.values).index(g) for g in avail if g in rvar.values]
    else:
        idx = [list(var_names.values).index(g) for g in avail]
    if len(idx) == 0:
        return None
    sub = X[:, idx]
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    # log1p normalize per spot total counts
    tot = np.array(a.obs.get("total_counts", a.obs.get("n_counts", np.ones(a.n_obs)))).flatten()
    tot = np.where(tot == 0, 1, tot)
    norm = sub / tot[:, None] * 1e4
    log1p = np.log1p(norm)
    score = log1p.mean(axis=1)
    out = pd.DataFrame({
        "spot_id": a.obs_names.astype(str),
        "n_genes_used": len(idx),
        "score": score,
    })
    if "array_row" in a.obs.columns:
        out["array_row"] = a.obs["array_row"].values
        out["array_col"] = a.obs["array_col"].values
    return out

def morans_I(values, rows, cols, max_neighbors=4):
    """Cheap Moran's I via grid-adjacency (Visium hex approximation)."""
    n = len(values)
    if n < 30: return np.nan
    v = np.asarray(values, dtype=float)
    if np.nanstd(v) == 0: return np.nan
    v_centered = v - np.nanmean(v)
    rows = np.asarray(rows); cols = np.asarray(cols)
    # build neighbor index by row/col integer grid (coarse)
    pos = list(zip(rows, cols))
    pos_to_idx = {p: i for i, p in enumerate(pos)}
    W_num = 0.0; W_denom_v2 = float(np.nansum(v_centered**2))
    W_sum = 0
    for i, (r, c) in enumerate(pos):
        # 6-connected neighbors (Visium hex). approximate with 4 cardinal + 2 diagonal
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1)]:
            j = pos_to_idx.get((r+dr, c+dc))
            if j is not None and not np.isnan(v[i]) and not np.isnan(v[j]):
                W_num += v_centered[i] * v_centered[j]
                W_sum += 1
    if W_sum == 0 or W_denom_v2 == 0: return np.nan
    I = (n / W_sum) * (W_num / W_denom_v2)
    return I

# C: GSE230424 PTC+HT TLS spatial
print("\n[C] GSE230424 PTC+HT TLS / B-cell / IGHV scoring (4 slides)")
g230_h5ads = sorted((EXTERNAL_H5AD/"GSE230424").glob("*/GSM*.raw.h5ad"))
c_rows = []
for h in g230_h5ads:
    sample = h.parent.name
    for set_name, genes in [("HLA_II",HLA_II_GENES),("B_cell",B_CELL_GENES),("TLS",TLS_GENES),("IGHV_AICDA",IGHV_GENES)]:
        sc = score_h5ad(h, genes)
        if sc is None: continue
        # spatial autocorrelation via Moran's I
        if "array_row" in sc.columns:
            I = morans_I(sc["score"].values, sc["array_row"].values, sc["array_col"].values)
        else: I = np.nan
        c_rows.append({
            "sample_id": sample, "dataset": "GSE230424", "condition": "PTC_HT",
            "gene_set": set_name, "n_genes_used": sc["n_genes_used"].iloc[0],
            "n_spots": len(sc), "mean_score": sc["score"].mean(),
            "median_score": sc["score"].median(), "p90": sc["score"].quantile(0.9),
            "morans_I": I,
        })
df_C = pd.DataFrame(c_rows)
df_C.to_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t", index=False)
print(df_C.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
sets = ["HLA_II","B_cell","TLS","IGHV_AICDA"]
xs = np.arange(len(sets))
samples = sorted(df_C["sample_id"].unique())
for i, s in enumerate(samples):
    sub = df_C[df_C.sample_id == s]
    sub_means = [sub[sub.gene_set == g]["mean_score"].iloc[0] if not sub[sub.gene_set == g].empty else np.nan for g in sets]
    ax.plot(xs, sub_means, "o-", label=s, alpha=0.85, lw=1.4)
ax.set_xticks(xs); ax.set_xticklabels(sets, fontsize=9.5)
ax.set_ylabel("Mean log-normalized score (per spot)")
ax.set_title("C1. GSE230424 PTC+HT 4 slides × 4 immune gene sets\n(HLA-II, B-cell, TLS, IGHV/AICDA)",
             fontsize=10)
ax.grid(axis="y", alpha=0.3)
ax.legend(loc="best", fontsize=9)

ax = axes[1]
for i, s in enumerate(samples):
    sub = df_C[df_C.sample_id == s]
    moran = [sub[sub.gene_set == g]["morans_I"].iloc[0] if not sub[sub.gene_set == g].empty else np.nan for g in sets]
    ax.plot(xs, moran, "o-", label=s, alpha=0.85, lw=1.4)
ax.axhline(0, ls=":", lw=0.6, color="gray")
ax.set_xticks(xs); ax.set_xticklabels(sets, fontsize=9.5)
ax.set_ylabel("Moran's I (spatial autocorrelation)")
ax.set_title("C2. Spatial niche organization (Moran's I)\nHigher I → niche-organized; near 0 → scattered",
             fontsize=10)
ax.grid(axis="y", alpha=0.3)
ax.legend(loc="best", fontsize=9)

fig.suptitle("S_F3  GSE230424 PTC+HT — TLS/B-cell/IGHV niche identification",
             fontsize=11, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F3_HT_TLS_spatial_panel.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# E: TROP2 spot-level distribution across 28 samples
print("\n[E] TROP2 (TACSTD2) per-sample spot distribution + Moran's I")
all_h5ads = []
all_h5ads += [(p, "GSE250521") for p in sorted(GSE250521_H5AD.glob("*/GSM*.raw.h5ad"))]
all_h5ads += [(p, "GSE230424") for p in sorted((EXTERNAL_H5AD/"GSE230424").glob("*/GSM*.raw.h5ad"))]
all_h5ads += [(p, "GSE248205") for p in sorted((EXTERNAL_H5AD/"GSE248205").glob("*/GSM*.raw.h5ad"))]
print(f"  total h5ads: {len(all_h5ads)}")

e_rows = []
for path, ds in all_h5ads:
    sample = path.parent.name
    sc = score_h5ad(path, TROP2_GENES)
    if sc is None: continue
    # condition lookup
    if ds == "GSE250521":
        meta_row = g521[g521.sample_id == sample]
        cond = meta_row["condition_inferred"].iloc[0] if not meta_row.empty else "?"
    else:
        meta_row = ext[ext.sample_id == sample]
        cond = meta_row["condition_inferred"].iloc[0] if not meta_row.empty else "?"
    if "array_row" in sc.columns:
        I = morans_I(sc["score"].values, sc["array_row"].values, sc["array_col"].values)
    else: I = np.nan
    pos = (sc["score"] > sc["score"].quantile(0.9)).sum()
    e_rows.append({
        "sample_id": sample, "dataset": ds, "condition": cond,
        "n_spots": len(sc), "TROP2_mean": sc["score"].mean(),
        "TROP2_median": sc["score"].median(), "TROP2_p90": sc["score"].quantile(0.9),
        "TROP2_p99": sc["score"].quantile(0.99),
        "TROP2_high_spot_n": int(pos),
        "TROP2_high_spot_frac": pos / len(sc) if len(sc) else np.nan,
        "morans_I_TROP2": I,
    })
df_E = pd.DataFrame(e_rows)
df_E.to_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t", index=False)
print(df_E.to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
for ds in ["GSE250521","GSE230424","GSE248205"]:
    sub = df_E[df_E.dataset == ds]
    for cond in sub["condition"].dropna().unique():
        sub_c = sub[sub.condition == cond]
        x = [f"{ds[:7]}.{cond}" for _ in range(len(sub_c))]
        ax.scatter(x, sub_c["TROP2_mean"].values, s=80,
                   color=COLORS.get(cond, "#666"), alpha=0.85, edgecolor="black", linewidth=0.6)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=9)
ax.axhline(0, ls=":", lw=0.6, color="gray")
ax.set_ylabel("Mean TROP2 (TACSTD2) log-norm per spot")
ax.set_title("E1. TROP2 per-sample mean expression across 28 samples (3 cohorts)", fontsize=10)
ax.grid(axis="y", alpha=0.3)

ax = axes[1]
for ds in ["GSE250521","GSE230424","GSE248205"]:
    sub = df_E[df_E.dataset == ds]
    for cond in sub["condition"].dropna().unique():
        sub_c = sub[sub.condition == cond]
        x = [f"{ds[:7]}.{cond}" for _ in range(len(sub_c))]
        ax.scatter(x, sub_c["morans_I_TROP2"].values, s=80,
                   color=COLORS.get(cond, "#666"), alpha=0.85, edgecolor="black", linewidth=0.6)
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=9)
ax.axhline(0, ls=":", lw=0.6, color="gray")
ax.axhline(0.3, ls="--", lw=0.6, color="#7B1F2A", alpha=0.4)
ax.set_ylabel("Moran's I (TROP2 spatial autocorrelation)")
ax.set_title("E2. TROP2 spatial niche organization\nHigher = niche-clustered (TROP2-high regions); ~0 = scattered", fontsize=10)
ax.grid(axis="y", alpha=0.3)

fig.suptitle("S_F5  TROP2 (TACSTD2) per-spot distribution + spatial niche test — 28 samples × 3 cohorts",
             fontsize=11, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F5_TROP2_spot_distribution.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ====================================================================
# Summary JSON
# ====================================================================
summary = {
    "date": "2026-05-06",
    "datasets": {
        "GSE250521": {"slides": 16, "stages": ["PT","PTC","LPTC","ATC"], "n_spots": int(len(g521))},
        "GSE230424": {"slides": 4, "condition": "PTC_HT", "n_spots": int((ext.dataset=="GSE230424").sum())},
        "GSE248205": {"slides": 8, "conditions": ["CONTROL","HT","GD"], "n_spots": int((ext.dataset=="GSE248205").sum())},
        "total_samples": 28, "total_spots_pooled": int(len(g521) + len(ext)),
    },
    "analyses": {
        "A_DM1_niche_2d": {"n_stages": int(df_A.shape[0]), "table": "spatial_A_DM1_niche_density.tsv"},
        "B_driver_orth": {"n_samples": int(df_B.shape[0]), "table": "spatial_B_driver_orthogonal_per_sample.tsv"},
        "C_HT_TLS": {"n_rows": int(df_C.shape[0]), "table": "spatial_C_HT_TLS_spatial.tsv"},
        "D_GSE248205_baseline": {"n_rows": int(df_D.shape[0]), "table": "spatial_D_GSE248205_baseline.tsv"},
        "E_TROP2_spot": {"n_samples": int(df_E.shape[0]), "table": "spatial_E_TROP2_spot_distribution.tsv"},
        "F_cross_cohort": {"n_rows": int(df_F.shape[0]), "table": "spatial_F_cross_cohort_integration.tsv"},
    },
    "figures": [
        "S_F1_DM1_niche_2d_density.png",
        "S_F2_driver_orthogonal_box.png",
        "S_F3_HT_TLS_spatial_panel.png",
        "S_F4_GSE248205_baseline.png",
        "S_F5_TROP2_spot_distribution.png",
        "S_F6_cross_cohort_integration.png",
    ],
}
(OUT/"spatial_full_summary.json").write_text(json.dumps(summary, indent=2))
print(f"\n[save] {OUT}/spatial_full_summary.json")
print("\nDONE — all 6 analyses complete.")
