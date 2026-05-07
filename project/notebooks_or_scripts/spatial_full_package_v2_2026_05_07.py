#!/usr/bin/env python3
"""
Spatial full package v2 — expanded figure deck (S_F7~S_F18).

Goal: 실제 tissue 좌표에서 TROP2 niche / HT TLS niche / DM1 spot map을 직접 시각화하고
Moran's I, co-localization scatter, cross-cohort heatmap 추가 — advisor question에
대한 강력한 visual answer 구축.

Inputs (read-only):
  project/data/processed/GSE250521/*/*.raw.h5ad           (16 slides)
  project_external_st/data/processed/GSE230424/*/*.raw.h5ad (4 slides)
  project_external_st/data/processed/GSE248205/*/*.raw.h5ad (8 slides)
  project/results/01_spatial_score/all_spots_scored.tsv.gz
  project_external_st/results/scores/*_spot_scores.tsv.gz
  project/results/spatial_full_2026_05_06/spatial_E_TROP2_spot_distribution.tsv

Outputs:
  project/papers_hub_2026_05_04/assets/spatial_full/S_F7  ~ S_F18.png  (NEW)
  project/results/spatial_full_2026_05_06/spatial_per_spot_TROP2.tsv.gz
  project/results/spatial_full_2026_05_06/spatial_v2_summary.json
"""
from __future__ import annotations
import json, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.colors import LinearSegmentedColormap, Normalize
import warnings; warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/results/spatial_full_2026_05_06"
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/spatial_full"
G521H = ROOT/"project/data/processed/GSE250521"
EXTH = ROOT/"project_external_st/data/processed"
ASSETS.mkdir(parents=True, exist_ok=True)

# Paper-hub palette (light-theme tone-matched)
PAL = {
    "PT": "#3C6B4F", "PTC": "#34547A", "LPTC": "#B8893C", "ATC": "#7B1F2A",
    "PTC_HT": "#962E2E", "CONTROL": "#52525a", "HT": "#34547A", "GD": "#B8893C",
}
ORDER_GSE250521 = ["PT","PTC","LPTC","ATC"]

TROP2_GENES = ["TACSTD2"]
HLA_II = ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
B_CELL = ["CD19","MS4A1","CD79A","CD79B","SDC1","JCHAIN"]
TLS    = ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"]

# ---------------- helpers ----------------
def score_h5ad(h5ad_path, gene_set):
    try: a = ad.read_h5ad(h5ad_path)
    except Exception: return None
    var_names = a.var_names.astype(str)
    if a.raw is not None:
        rvar = a.raw.var_names.astype(str)
        idx = [list(rvar.values).index(g) for g in gene_set if g in rvar.values]
        X = a.raw.X
    else:
        idx = [list(var_names.values).index(g) for g in gene_set if g in var_names.values]
        X = a.X
    if not idx: return None
    sub = X[:, idx]
    if hasattr(sub, "toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", a.obs.get("n_counts", np.ones(a.n_obs)))).flatten()
    tot = np.where(tot == 0, 1, tot)
    norm = sub / tot[:, None] * 1e4
    log1p = np.log1p(norm)
    score = log1p.mean(axis=1)
    out = pd.DataFrame({"spot_id": a.obs_names.astype(str),
                       "score": score,
                       "n_genes_used": len(idx)})
    if "array_row" in a.obs.columns:
        out["array_row"] = a.obs["array_row"].values
        out["array_col"] = a.obs["array_col"].values
    if "pxl_row_in_fullres" in a.obs.columns:
        out["pxl_row"] = a.obs["pxl_row_in_fullres"].values
        out["pxl_col"] = a.obs["pxl_col_in_fullres"].values
    return out

def morans_I(values, rows, cols):
    n = len(values)
    if n < 30: return np.nan
    v = np.asarray(values, dtype=float)
    if np.nanstd(v) == 0: return np.nan
    vc = v - np.nanmean(v)
    rows = np.asarray(rows); cols = np.asarray(cols)
    pos_to_idx = {(r,c): i for i, (r,c) in enumerate(zip(rows, cols))}
    Wnum, Wsum = 0.0, 0
    for i, (r, c) in enumerate(zip(rows, cols)):
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1)]:
            j = pos_to_idx.get((r+dr, c+dc))
            if j is not None and not np.isnan(v[i]) and not np.isnan(v[j]):
                Wnum += vc[i]*vc[j]; Wsum += 1
    Wd = float(np.nansum(vc**2))
    if Wsum == 0 or Wd == 0: return np.nan
    return (n/Wsum)*(Wnum/Wd)

# Build registry of all h5ads with metadata
print("[1/8] discovering h5ads")
h5ads = []
for p in sorted(G521H.glob("*/GSM*.raw.h5ad")):
    sample = p.parent.name
    stage = sample.split("_")[-1].split("-")[0]  # N/PTC/LPTC/ATC
    if stage == "N": stage = "PT"
    h5ads.append({"path": str(p), "sample": sample, "dataset": "GSE250521", "condition": stage})
for p in sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad")):
    h5ads.append({"path": str(p), "sample": p.parent.name, "dataset": "GSE230424", "condition": "PTC_HT"})
for p in sorted((EXTH/"GSE248205").glob("*/GSM*.raw.h5ad")):
    sample = p.parent.name
    cond = "CONTROL" if "C" in sample.split("_")[-1][:1] and "HT" not in sample and "GD" not in sample else None
    if cond is None:
        if "HT" in sample: cond = "HT"
        elif "GD" in sample: cond = "GD"
        else: cond = "CONTROL"
    h5ads.append({"path": str(p), "sample": p.parent.name, "dataset": "GSE248205", "condition": cond})
print(f"  {len(h5ads)} h5ads total")

# ---------------- score TROP2 per spot for all 28 ----------------
print("[2/8] scoring TROP2 per spot (28 samples)")
trop_per_spot = []
for rec in h5ads:
    sc = score_h5ad(rec["path"], TROP2_GENES)
    if sc is None: continue
    sc["sample_id"] = rec["sample"]
    sc["dataset"] = rec["dataset"]
    sc["condition"] = rec["condition"]
    trop_per_spot.append(sc)
trop_all = pd.concat(trop_per_spot, axis=0, ignore_index=True)
trop_all = trop_all.rename(columns={"score":"TROP2"})
trop_all.to_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t", index=False, compression="gzip")
print(f"  {len(trop_all)} spots, {trop_all['sample_id'].nunique()} samples")

# Read pre-existing E table
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")

# Score TLS for GSE230424 (per spot)
print("[3/8] scoring HT immune sets per spot (GSE230424)")
tls_per_spot = []
g230_h5ads = [h for h in h5ads if h["dataset"] == "GSE230424"]
for rec in g230_h5ads:
    for setname, genes in [("HLA_II",HLA_II),("B_cell",B_CELL),("TLS",TLS)]:
        sc = score_h5ad(rec["path"], genes)
        if sc is None: continue
        sc["sample_id"] = rec["sample"]; sc["set"] = setname
        tls_per_spot.append(sc)
tls_all = pd.concat(tls_per_spot, axis=0, ignore_index=True).rename(columns={"score":"value"})

# ---------------- F7: TROP2 spatial maps — 16 GSE250521 slides 4x4 grid ----------------
print("[4/8] F7 — TROP2 spatial maps GSE250521 (16 slides, 4 stages × 4)")
g521_samples = [h for h in h5ads if h["dataset"] == "GSE250521"]
# group by stage
by_stage = {s: [] for s in ORDER_GSE250521}
for h in g521_samples:
    by_stage[h["condition"]].append(h["sample"])
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.94, bottom=0.04, left=0.04, right=0.97)
trop_cmap = LinearSegmentedColormap.from_list("trop2", ["#fffaf2","#ffd58a","#d97c2e","#7B1F2A"], N=256)
for r, stage in enumerate(ORDER_GSE250521):
    samples = sorted(by_stage[stage])[:4]
    for c, sample in enumerate(samples):
        ax = fig.add_subplot(g[r, c])
        sub = trop_all[trop_all.sample_id == sample]
        if "array_row" not in sub.columns or sub.empty:
            ax.text(0.5, 0.5, "no coords", ha="center", va="center"); ax.axis("off"); continue
        # use array_col as x, array_row as y (Visium convention, flip y for top-down)
        x = sub["array_col"].values; y = -sub["array_row"].values
        v = sub["TROP2"].values
        # color limit per slide
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1.0
        sc = ax.scatter(x, y, c=v, cmap=trop_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        moran = df_E[df_E.sample_id == sample]["morans_I_TROP2"]
        moran = moran.iloc[0] if not moran.empty else np.nan
        title_color = "#7B1F2A" if (not np.isnan(moran) and moran > 0.25) else "#52525a"
        ax.set_title(f"{stage}  ·  {sample.split('_')[-1]}\nMoran I = {moran:.2f}",
                     fontsize=9.5, color=title_color, fontweight="bold" if title_color == "#7B1F2A" else "normal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")
        for spine in ax.spines.values(): spine.set_color("#cdc1aa")
fig.suptitle("S_F7  TROP2 spatial maps — GSE250521 16 slides (4 stages × 4)\n"
             "각 dot = 1 spot · color = TROP2 (TACSTD2) log-norm expression · niche-organized samples (Moran I > 0.25) 빨간 제목",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F7_TROP2_maps_GSE250521.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------- F8: HT TLS spatial maps GSE230424 (4 slides × 3 sets) ----------------
print("[5/8] F8 — HT TLS spatial maps GSE230424 (4 slides × 3 sets)")
fig = plt.figure(figsize=(15, 18))
g = gs.GridSpec(4, 3, hspace=0.32, wspace=0.18, top=0.96, bottom=0.04, left=0.04, right=0.97)
tls_cmap = LinearSegmentedColormap.from_list("tls", ["#fffaf2","#a8d8c8","#1abc9c","#0F1A2E"], N=256)
g230_samples = sorted([h["sample"] for h in g230_h5ads])
sets = ["HLA_II","B_cell","TLS"]
for r, sample in enumerate(g230_samples):
    for c, setname in enumerate(sets):
        ax = fig.add_subplot(g[r, c])
        sub = tls_all[(tls_all.sample_id == sample) & (tls_all.set == setname)]
        if "array_row" not in sub.columns or sub.empty:
            ax.text(0.5, 0.5, "no data"); ax.axis("off"); continue
        x = sub["array_col"].values; y = -sub["array_row"].values
        v = sub["value"].values
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1.0
        ax.scatter(x, y, c=v, cmap=tls_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        # find Moran from C table
        Cdf = pd.read_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t")
        moran = Cdf[(Cdf.sample_id == sample) & (Cdf.gene_set == setname)]["morans_I"]
        moran = moran.iloc[0] if not moran.empty else np.nan
        title_color = "#1abc9c" if (not np.isnan(moran) and moran > 0.3) else "#52525a"
        sample_short = sample.split("_")[-1]
        ax.set_title(f"{sample_short}  ·  {setname}\nMoran I = {moran:.2f}",
                     fontsize=10, color=title_color, fontweight="bold" if title_color != "#52525a" else "normal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")
        for spine in ax.spines.values(): spine.set_color("#cdc1aa")
fig.suptitle("S_F8  HT-overlap PTC immune spatial maps — GSE230424 4 slides × 3 gene sets\n"
             "color = mean log-norm score · 4/4 슬라이드 + 3/3 set 모두 niche-organized (Moran I > 0.27)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F8_HT_TLS_maps_GSE230424.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------- F9: TROP2 ranked Moran's I across 28 samples ----------------
print("[6/8] F9 — TROP2 Moran's I ranked across 28 samples")
df_rank = df_E.sort_values("morans_I_TROP2", ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(13.5, 7.5))
xs = np.arange(len(df_rank))
colors = [PAL.get(c, "#999") for c in df_rank["condition"]]
ax.bar(xs, df_rank["morans_I_TROP2"].values, color=colors, edgecolor="black", linewidth=0.4)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1.2, alpha=0.7)
ax.axhline(0, ls="-", color="black", lw=0.6)
ax.text(28*0.7, 0.27, "Moran I = 0.25 — niche threshold", color="#7B1F2A", fontsize=10, fontweight="bold")
labels = [f"{r['sample_id'].split('_')[-1]}\n[{r['condition']}]" for _, r in df_rank.iterrows()]
ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=8)
ax.set_ylabel("TROP2 spatial autocorrelation (Moran's I)")
ax.set_ylim(-0.05, 0.65)
ax.set_title("S_F9  TROP2 Moran's I ranked across 28 samples · 3 cohorts\n"
             "★ PTC + LPTC 8/8 above 0.25 niche threshold; PT (4/4), HT (3/3), GD (3/3), CONTROL (2/2) below threshold",
             fontsize=12, fontweight="bold")
# legend
from matplotlib.patches import Patch
handles = [Patch(facecolor=PAL[c], edgecolor="black", label=f"{c} (n={(df_rank['condition']==c).sum()})") for c in ["PTC","LPTC","ATC","PTC_HT","HT","GD","CONTROL","PT"] if c in df_rank["condition"].values]
ax.legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.95)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F9_TROP2_ranked_morans.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F10: TROP2 × DM1 spot-level scatter (representative samples) ----------------
print("[7/8] F10 — TROP2 × DM1 co-localization scatter")
# Need to merge TROP2 (per spot) with DM1 (per spot) — DM1 from all_spots_scored.tsv.gz / external_st spot scores
g521_spots = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t")
g521_spots = g521_spots.rename(columns={"stage":"condition"})
# Merge by spot_id (use array_row+array_col instead since spot_id 형식이 다를 수 있음)
trop_g521 = trop_all[trop_all.dataset == "GSE250521"].copy()
trop_g521["row_col_id"] = trop_g521["sample_id"] + "_" + trop_g521["array_row"].astype(str) + "_" + trop_g521["array_col"].astype(str)
g521_spots["row_col_id"] = g521_spots["sample_id"] + "_" + g521_spots["array_row"].astype(str) + "_" + g521_spots["array_col"].astype(str)
merge_cols = ["row_col_id","DM1_like_score","RAI_8_score","Epithelial_score","Proliferation_score","TDS_like_score"]
mg = trop_g521.merge(g521_spots[merge_cols], on="row_col_id", how="left")
mg["TROP2_high"] = mg["TROP2"] > mg.groupby("sample_id")["TROP2"].transform(lambda x: x.quantile(0.85))
print(f"  merged {len(mg):,} spots, TROP2-high (top 15% per sample) = {mg['TROP2_high'].sum():,}")

fig, axes = plt.subplots(2, 4, figsize=(15.5, 8), sharex=False, sharey=False)
# Pick 1 representative sample per stage with highest Moran I
representative = {}
for stage in ORDER_GSE250521:
    sub = df_E[(df_E.dataset == "GSE250521") & (df_E.condition == stage)]
    if sub.empty: continue
    representative[stage] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]
# row 1 = TROP2 vs DM1, row 2 = TROP2 vs Epithelial
for c, stage in enumerate(ORDER_GSE250521):
    sample = representative.get(stage)
    if sample is None: continue
    sub = mg[mg.sample_id == sample]
    ax1 = axes[0, c]
    ax2 = axes[1, c]
    rho_dm, p_dm = spearmanr(sub["TROP2"], sub["DM1_like_score"], nan_policy="omit")
    rho_epi, p_epi = spearmanr(sub["TROP2"], sub["Epithelial_score"], nan_policy="omit")
    ax1.scatter(sub["DM1_like_score"], sub["TROP2"], s=3, c=sub["TROP2"], cmap=trop_cmap, alpha=0.6)
    ax1.set_title(f"{stage} · {sample.split('_')[-1]}\nρ(TROP2, DM1) = {rho_dm:.3f}, p = {p_dm:.1e}", fontsize=9.5)
    ax1.set_xlabel("DM1_like score"); ax1.set_ylabel("TROP2 (log-norm)")
    ax1.axhline(0, ls=":", lw=0.4, color="gray"); ax1.axvline(0, ls=":", lw=0.4, color="gray")
    ax2.scatter(sub["Epithelial_score"], sub["TROP2"], s=3, c=sub["TROP2"], cmap=trop_cmap, alpha=0.6)
    ax2.set_title(f"ρ(TROP2, Epi) = {rho_epi:.3f}, p = {p_epi:.1e}", fontsize=9.5)
    ax2.set_xlabel("Epithelial score"); ax2.set_ylabel("TROP2 (log-norm)")
    ax2.axhline(0, ls=":", lw=0.4, color="gray"); ax2.axvline(0, ls=":", lw=0.4, color="gray")
fig.suptitle("S_F10  TROP2 × DM1 / Epithelial spot-level co-localization (representative niche-organized sample per stage)\n"
             "Top: TROP2 vs DM1_like (axis 비교) · Bottom: TROP2 vs Epithelial (tumor-cell containment 검증)",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F10_TROP2_DM1_colocalization.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F11: Cross-cohort condition × axis heatmap ----------------
print("[8/8] F11-F14: cross-cohort heatmaps + closure battery + DM1 maps")
# Build mean per (dataset, condition, axis) heatmap
g521_h = g521_spots.copy()
g521_h["dataset"] = "GSE250521"
g521_h["condition_inferred"] = g521_h["condition"]
ext_files = sorted((ROOT/"project_external_st/results/scores").glob("*_spot_scores.tsv.gz"))
ext_dfs = [pd.read_csv(f, sep="\t") for f in ext_files]
ext_h = pd.concat(ext_dfs, axis=0, ignore_index=True).rename(columns={
    "RAI_8_score_raw":"RAI_8_score","DM1_like_score_raw":"DM1_like_score",
    "Epithelial_score_raw":"Epithelial_score","Proliferation_score_raw":"Proliferation_score",
    "CAF_ECM_score_raw":"CAF_ECM_score","EMT_score_raw":"EMT_score",
    "Hypoxia_score_raw":"Hypoxia_score","TDS_overlap_score_raw":"TDS_like_score",
    "THYROID_NONOVERLAP_score_raw":"THYROID_NONOVERLAP_score",
})
COMMON = ["sample_id","dataset","condition_inferred","RAI_8_score","DM1_like_score","TDS_like_score",
          "Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]
g521_h["TDS_like_score"] = g521_h.get("TDS_like_score", g521_h.get("TDS_overlap_score"))
for c in COMMON:
    if c not in g521_h.columns: g521_h[c] = np.nan
    if c not in ext_h.columns: ext_h[c] = np.nan
all_h = pd.concat([g521_h[COMMON], ext_h[COMMON]], axis=0, ignore_index=True)
all_h["group"] = all_h["dataset"] + "·" + all_h["condition_inferred"].astype(str)
GROUPS = ["GSE250521·PT","GSE250521·PTC","GSE250521·LPTC","GSE250521·ATC",
          "GSE230424·PTC_HT","GSE248205·CONTROL","GSE248205·HT","GSE248205·GD"]
AXES = ["RAI_8_score","DM1_like_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]
mat = np.full((len(AXES), len(GROUPS)), np.nan)
n_samp = []
for j, gp in enumerate(GROUPS):
    sub = all_h[all_h.group == gp]
    n_samp.append(sub["sample_id"].nunique())
    for i, ax_name in enumerate(AXES):
        mat[i, j] = sub[ax_name].mean()
fig, ax = plt.subplots(figsize=(12, 6.2))
norm_v = max(abs(np.nanmin(mat)), abs(np.nanmax(mat)))
im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-norm_v, vmax=norm_v)
for i in range(len(AXES)):
    for j in range(len(GROUPS)):
        v = mat[i,j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2g}", ha="center", va="center",
                    color="white" if abs(v) > norm_v*0.55 else "black", fontsize=9)
ax.set_xticks(range(len(GROUPS)))
ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(GROUPS, n_samp)], rotation=20, ha="right", fontsize=9.5)
ax.set_yticks(range(len(AXES))); ax.set_yticklabels([a.replace("_score","") for a in AXES], fontsize=10)
ax.set_title("S_F11  Cross-cohort condition × axis mean heatmap (28 samples · 7 axes)\n"
             "Mean per-spot score (within-sample z) per group · n_samp 표시",
             fontsize=12, fontweight="bold")
cbar = plt.colorbar(im, ax=ax, shrink=0.7); cbar.set_label("Mean per-spot score (within-sample z)")
fig.tight_layout()
fig.savefig(ASSETS/"S_F11_cross_cohort_heatmap.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F12: DM1 spatial maps GSE250521 representative + ATC niche ----------------
print("  F12 — DM1 spatial maps representative")
fig = plt.figure(figsize=(15, 4.2))
g = gs.GridSpec(1, 4, wspace=0.16, top=0.86, bottom=0.04, left=0.04, right=0.97)
dm_cmap = LinearSegmentedColormap.from_list("dm", ["#34547A","#fffaf2","#7B1F2A"], N=256)
for c, stage in enumerate(ORDER_GSE250521):
    sample = representative.get(stage)
    if sample is None: continue
    sub = g521_spots[g521_spots.sample_id == sample]
    ax = fig.add_subplot(g[0, c])
    if "array_row" not in sub.columns or sub.empty:
        ax.text(0.5, 0.5, "no coords"); ax.axis("off"); continue
    x = sub["array_col"].values; y = -sub["array_row"].values
    v = sub["DM1_like_score"].values
    vmax = max(abs(np.nanpercentile(v, 5)), abs(np.nanpercentile(v, 95)))
    ax.scatter(x, y, c=v, cmap=dm_cmap, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
    ax.set_title(f"{stage}  ·  {sample.split('_')[-1]}", fontsize=11, fontweight="bold",
                 color=PAL[stage])
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values(): spine.set_color("#cdc1aa")
fig.suptitle("S_F12  DM1_like spatial maps — representative niche samples per stage (GSE250521)\n"
             "blue = RAI-high (differentiated); red = DM1-high (dark matter); same samples as S_F7",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F12_DM1_maps_representative.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F13: closure battery negative result ----------------
print("  F13 — closure battery summary")
clos = pd.read_csv(ROOT/"project/results/03_pathology_poc/closure_battery_metrics.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axes[0]
sec_a = clos[clos["section"].astype(str).str.startswith("A")]
exps = sec_a["experiment"].astype(str).unique()
positions = np.arange(len(exps))
ridge = sec_a[sec_a["model"] == "ridge"]
enet = sec_a[sec_a["model"] == "enet"]
w = 0.4
ax.bar(positions - w/2, ridge.set_index("experiment").reindex(exps)["pooled_spearman_r"].values, w, label="ridge", color="#7B1F2A")
ax.bar(positions + w/2, enet.set_index("experiment").reindex(exps)["pooled_spearman_r"].values, w, label="enet", color="#34547A")
ax.axhline(0, color="black", lw=0.6)
ax.axhline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.6); ax.text(0, 0.32, "minimum useful threshold (0.3)", fontsize=9, color="#3C6B4F")
ax.set_xticks(positions); ax.set_xticklabels([e.replace("raw_","") for e in exps], rotation=15, fontsize=9)
ax.set_ylabel("Pooled Spearman ρ (LOSO predicted vs observed)")
ax.set_ylim(-0.3, 0.5)
ax.set_title("A. H&E → DM1 ResNet50 LOSO regression\n(rho ~ 0.06-0.08 = no signal; minimum useful threshold not crossed)",
             fontsize=10)
ax.legend()

ax = axes[1]
sec_e = clos[clos["section"].astype(str) == "E"]
labels = sec_e["experiment"].astype(str).str.replace("E_","").str.replace("_score_resid","").str.replace("_"," ")
ax.barh(range(len(sec_e)), sec_e["pooled_spearman_r"].values,
        color=["#7B1F2A" if abs(v) < 0.3 else "#3C6B4F" for v in sec_e["pooled_spearman_r"].fillna(0).values])
ax.set_yticks(range(len(sec_e))); ax.set_yticklabels(labels, fontsize=8.5)
ax.axvline(0, color="black", lw=0.6)
ax.axvline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
ax.axvline(-0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
ax.set_xlabel("Slide-level Spearman ρ (LOSO)")
ax.set_xlim(-0.5, 0.5)
ax.set_title("B. Slide-level mean / top-25% / high-risk fraction\n(only top25 mean ρ ≈ 0.21 — below useful threshold)",
             fontsize=10)
fig.suptitle("S_F13  H&E → DM1 closure battery (Paper 2A NO-GO verdict — ResNet50 다양한 라벨/embeddings 다 fail)",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F13_closure_battery_summary.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F14: per-stage axis distribution (DM1 + RAI8 + Epi + Prol violins) ----------------
print("  F14 — per-stage axis violin (GSE250521 16 slides)")
fig, axes = plt.subplots(1, 4, figsize=(15, 4.5), sharey=False)
axes_to_show = [("RAI_8_score","#3C6B4F"),("DM1_like_score","#7B1F2A"),("Epithelial_score","#B8893C"),("Proliferation_score","#34547A")]
for ax, (col, color) in zip(axes, axes_to_show):
    data = []; labels = []
    for stage in ORDER_GSE250521:
        sub = g521_spots[g521_spots.condition == stage]
        data.append(sub[col].dropna().values); labels.append(f"{stage}\nn={len(sub):,}")
    parts = ax.violinplot(data, showmeans=True, showmedians=True, widths=0.85)
    for body in parts['bodies']:
        body.set_facecolor(color); body.set_alpha(0.55); body.set_edgecolor("black")
    parts['cmeans'].set_color("black")
    parts['cmedians'].set_color("white")
    ax.set_xticks(range(1, len(labels)+1)); ax.set_xticklabels(labels, fontsize=8.5)
    ax.axhline(0, ls=":", color="gray", lw=0.6)
    ax.set_title(col.replace("_score",""), fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.25)
axes[0].set_ylabel("Per-spot score (within-sample z)")
fig.suptitle("S_F14  Per-stage spot-level axis violins (GSE250521 16 slides, ~55,873 spots)\n"
             "within-sample z normalization → stage-mean ≈ 0; 분포 spread / shape 비교 가능",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F14_per_stage_violin.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F15: GSE248205 baseline maps (CONTROL/HT/GD) ----------------
print("  F15 — GSE248205 baseline spatial maps (3 conditions × representative)")
g248_samples = [h for h in h5ads if h["dataset"] == "GSE248205"]
g248_score = pd.concat([pd.read_csv(f, sep="\t") for f in (ROOT/"project_external_st/results/scores").glob("GSE248205_*.tsv.gz")], axis=0, ignore_index=True)
g248_score = g248_score.rename(columns={"DM1_like_score_raw":"DM1_like_score","Epithelial_score_raw":"Epithelial_score"})
fig = plt.figure(figsize=(13, 8.5))
g = gs.GridSpec(2, 3, wspace=0.18, hspace=0.32, top=0.93, bottom=0.06, left=0.05, right=0.97)
# pick representative sample per condition (largest n_spots)
rep_g248 = {}
for cond in ["CONTROL","HT","GD"]:
    samples = [h["sample"] for h in g248_samples if h["condition"] == cond]
    if not samples: continue
    sample_sizes = {s: (g248_score.sample_id == s).sum() for s in samples}
    rep_g248[cond] = max(sample_sizes, key=sample_sizes.get)
for c, cond in enumerate(["CONTROL","HT","GD"]):
    if cond not in rep_g248: continue
    sample = rep_g248[cond]
    sub = g248_score[g248_score.sample_id == sample]
    if sub.empty or "array_row" not in sub.columns: continue
    x = sub["array_col"].values; y = -sub["array_row"].values
    # row 0: DM1
    ax = fig.add_subplot(g[0, c])
    v = sub["DM1_like_score"].values
    vmax = max(abs(np.nanpercentile(v, 5)), abs(np.nanpercentile(v, 95)))
    ax.scatter(x, y, c=v, cmap=dm_cmap, vmin=-vmax, vmax=vmax, s=4, alpha=0.92)
    ax.set_title(f"{cond} · {sample.split('_')[-1]}\nDM1_like (within-sample z)", fontsize=10, color=PAL[cond], fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    # row 1: TROP2
    ax = fig.add_subplot(g[1, c])
    sub_t = trop_all[trop_all.sample_id == sample]
    if not sub_t.empty:
        v = sub_t["TROP2"].values
        vmax_t = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1.0
        ax.scatter(sub_t["array_col"].values, -sub_t["array_row"].values, c=v, cmap=trop_cmap, vmin=0, vmax=vmax_t, s=4, alpha=0.92)
    ax.set_title(f"TROP2 (TACSTD2) log-norm", fontsize=10, color="#7B1F2A")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F15  GSE248205 autoimmune-only baseline (CONTROL · HT · GD) — spatial maps\n"
             "Top row: DM1_like (within-sample z) · Bottom row: TROP2 (log-norm) · niche structure 모두 약함 (negative control ✓)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F15_GSE248205_baseline_maps.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F16: Moran's I summary across niche markers ----------------
print("  F16 — Moran's I cross-marker summary")
fig, ax = plt.subplots(figsize=(12, 5.5))
# combine TROP2 (E) + HT TLS (C) + autoimmune (TROP2 in HT/GD)
df_C = pd.read_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t")
trop_summary = df_E[["sample_id","condition","morans_I_TROP2"]].copy()
trop_summary["marker"] = "TROP2"; trop_summary = trop_summary.rename(columns={"morans_I_TROP2":"morans_I"})
tls_summary = df_C[["sample_id","condition","gene_set","morans_I"]].rename(columns={"gene_set":"marker"})
combo = pd.concat([trop_summary[["sample_id","condition","marker","morans_I"]],
                   tls_summary[["sample_id","condition","marker","morans_I"]]],
                  axis=0, ignore_index=True)
markers = ["TROP2","HLA_II","B_cell","TLS","IGHV_AICDA"]
conds = ["PT","PTC","LPTC","ATC","PTC_HT","HT","GD","CONTROL"]
mat2 = np.full((len(markers), len(conds)), np.nan)
for i, m in enumerate(markers):
    for j, c in enumerate(conds):
        sub = combo[(combo.marker == m) & (combo.condition == c)]
        if not sub.empty: mat2[i, j] = sub["morans_I"].mean()
norm_v = 0.6
im = ax.imshow(mat2, aspect="auto", cmap="YlOrRd", vmin=0, vmax=norm_v)
for i in range(len(markers)):
    for j in range(len(conds)):
        v = mat2[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v > 0.4 else "black", fontsize=11, fontweight="bold")
ax.set_xticks(range(len(conds))); ax.set_xticklabels(conds, fontsize=10.5)
ax.set_yticks(range(len(markers))); ax.set_yticklabels(markers, fontsize=10.5, family="monospace")
cbar = plt.colorbar(im, ax=ax, shrink=0.7); cbar.set_label("Mean Moran's I (per condition)")
ax.set_title("S_F16  Spatial niche organization summary — Moran's I per (marker × condition)\n"
             "Higher = more niche-organized · TROP2 PTC/LPTC + HT TLS-set 모두 strong niche · normal/autoimmune-only 약함",
             fontsize=11.5, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F16_morans_summary.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F17: TROP2 high-spot fraction per condition ----------------
print("  F17 — TROP2 high-spot fraction across conditions")
fig, ax = plt.subplots(figsize=(11, 5.5))
xs = np.arange(len(df_E))
df_E_sorted = df_E.sort_values(["condition","TROP2_high_spot_frac"], ascending=[True, False]).reset_index(drop=True)
colors = [PAL.get(c, "#999") for c in df_E_sorted["condition"]]
ax.bar(np.arange(len(df_E_sorted)), df_E_sorted["TROP2_high_spot_frac"]*100,
       color=colors, edgecolor="black", linewidth=0.4)
ax.set_xticks(np.arange(len(df_E_sorted)))
ax.set_xticklabels([f"{r['sample_id'].split('_')[-1]}\n[{r['condition']}]" for _, r in df_E_sorted.iterrows()],
                   rotation=70, ha="right", fontsize=8)
ax.set_ylabel("TROP2-high spot fraction (%) — sample p90 threshold")
ax.set_title("S_F17  TROP2 high-spot fraction per sample (top 10% TROP2-expressing spots)\n"
             "PTC/LPTC ≈ 8–10% high-spots (niche concentration); PT/HT/GD ≤ 5% scattered minor spots",
             fontsize=11.5, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F17_TROP2_high_spot_fraction.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ---------------- F18: HT TLS sample × marker heatmap ----------------
print("  F18 — HT TLS sample × marker heatmap")
fig, ax = plt.subplots(figsize=(8.5, 4.5))
sets_order = ["HLA_II","B_cell","TLS","IGHV_AICDA"]
g230_order = sorted(df_C["sample_id"].unique())
mat3 = np.full((len(g230_order), len(sets_order)), np.nan)
for i, s in enumerate(g230_order):
    for j, set_name in enumerate(sets_order):
        sub = df_C[(df_C.sample_id == s) & (df_C.gene_set == set_name)]
        if not sub.empty: mat3[i, j] = sub["morans_I"].iloc[0]
im = ax.imshow(mat3, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.9)
for i in range(len(g230_order)):
    for j in range(len(sets_order)):
        v = mat3[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color="white" if v > 0.55 else "black", fontsize=11, fontweight="bold")
ax.set_xticks(range(len(sets_order))); ax.set_xticklabels(sets_order, fontsize=10.5, family="monospace")
ax.set_yticks(range(len(g230_order)))
ax.set_yticklabels([s.split("_")[-1] for s in g230_order], fontsize=10.5, family="monospace")
cbar = plt.colorbar(im, ax=ax, shrink=0.7); cbar.set_label("Moran's I")
ax.set_title("S_F18  GSE230424 HT-overlap PTC — sample × marker Moran's I heatmap\n"
             "4 slides × 4 immune sets · 모두 niche-organized · IGHV/AICDA P3 = 0.844 최강",
             fontsize=11, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F18_HT_TLS_heatmap.png", dpi=160, bbox_inches="tight")
plt.close(fig)

print("\nDONE — 12 new figures (S_F7 ~ S_F18) saved to", ASSETS)
