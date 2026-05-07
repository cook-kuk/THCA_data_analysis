#!/usr/bin/env python3
"""Spatial v4 — 6 more figures (S_F27 ~ S_F32)."""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr
import warnings; warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/results/spatial_full_2026_05_06"
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/spatial_full"
G521H = ROOT/"project/data/processed/GSE250521"
EXTH = ROOT/"project_external_st/data/processed"

PAL = {"PT":"#3C6B4F","PTC":"#34547A","LPTC":"#B8893C","ATC":"#7B1F2A",
       "PTC_HT":"#962E2E","CONTROL":"#52525a","HT":"#34547A","GD":"#B8893C"}
ORDER = ["PT","PTC","LPTC","ATC"]
trop_cmap = LinearSegmentedColormap.from_list("trop2",["#fffaf2","#ffd58a","#d97c2e","#7B1F2A"],N=256)
dm_cmap   = LinearSegmentedColormap.from_list("dm",["#34547A","#fffaf2","#7B1F2A"],N=256)

def score_h5ad(p, gs_list):
    a = ad.read_h5ad(p)
    rvar = a.raw.var_names.astype(str) if a.raw is not None else a.var_names.astype(str)
    X = a.raw.X if a.raw is not None else a.X
    idx = [list(rvar.values).index(g) for g in gs_list if g in rvar.values]
    if not idx: return None, []
    sub = X[:, idx]
    if hasattr(sub,"toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", np.ones(a.n_obs))).flatten()
    tot = np.where(tot==0,1,tot)
    norm_log = np.log1p(sub/tot[:,None]*1e4)
    avail = [g for g in gs_list if g in rvar.values]
    return norm_log, avail

def coords(p):
    a = ad.read_h5ad(p)
    return (a.obs["array_row"].values, a.obs["array_col"].values) if "array_row" in a.obs.columns else (None, None)

print("[load] data")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
trop_all = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")

# representative niche samples (highest Moran I per stage)
rep = {}
for s in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==s)]
    if not sub.empty:
        rep[s] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]

# ===== S_F27: 3-channel RGB overlay (TROP2 / DM1 / Epithelial) per stage =====
print("[F27] 3-channel RGB overlay per stage")
fig = plt.figure(figsize=(15, 5))
g = gs.GridSpec(1, 4, wspace=0.18, top=0.86, bottom=0.04, left=0.04, right=0.97)
for c, stage in enumerate(ORDER):
    sample = rep.get(stage)
    if sample is None: continue
    sub_t = trop_all[trop_all.sample_id == sample]
    sub_g = g521[g521.sample_id == sample]
    sub_g["row_col"] = sub_g["sample_id"]+"_"+sub_g["array_row"].astype(str)+"_"+sub_g["array_col"].astype(str)
    sub_t["row_col"] = sub_t["sample_id"]+"_"+sub_t["array_row"].astype(str)+"_"+sub_t["array_col"].astype(str)
    mg = sub_t.merge(sub_g[["row_col","DM1_like_score","Epithelial_score"]], on="row_col")
    if mg.empty: continue
    ax = fig.add_subplot(g[0, c])
    # normalize each channel 0-1
    def to01(x):
        x = np.array(x, dtype=float)
        x = x - np.nanmin(x); x = x/np.nanmax(x) if np.nanmax(x)>0 else x
        return np.nan_to_num(x, 0)
    R = to01(mg["TROP2"].values)        # red = TROP2
    G = to01(mg["Epithelial_score"].values)  # green = Epithelial
    B = to01(np.maximum(mg["DM1_like_score"].values, 0))  # blue = DM1+
    rgb = np.stack([R, G, B], axis=1)
    x = mg["array_col"].values; y = -mg["array_row"].values
    ax.scatter(x, y, c=rgb, s=4, alpha=0.92)
    ax.set_title(f"{stage} · {sample.split('_')[-1]}\nR=TROP2  G=Epithelial  B=DM1", fontsize=10, color=PAL[stage], fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F27  3-channel RGB overlay per stage — TROP2 (R) · Epithelial (G) · DM1 (B)\n"
             "Yellow = TROP2+Epi (tumor TROP2 niche) · Magenta = TROP2+DM1 · Cyan = Epi+DM1 · White = all 3",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F27_RGB_overlay.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# ===== S_F28: BRAF / HRAS / KRAS / NRAS spatial expression per stage representative =====
print("[F28] Driver gene spatial expression")
DRIVERS = ["BRAF","HRAS","KRAS","NRAS"]
fig = plt.figure(figsize=(15, 16))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.04, right=0.97)
g521_paths = {}
for p in sorted(G521H.glob("*/GSM*.raw.h5ad")):
    g521_paths[p.parent.name] = p
for r, gene in enumerate(DRIVERS):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if sample is None or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, avail = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0:
            ax.text(0.5,0.5,f"{gene} not found"); ax.axis("off"); continue
        rows, cols = coords(g521_paths[sample])
        if rows is None: ax.axis("off"); continue
        v = scores[:,0]
        x = cols; y = -np.array(rows)
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(x, y, c=v, cmap=trop_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=12, fontweight="bold", color="#7B1F2A", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F28  Driver gene RNA expression spatial maps (BRAF / HRAS / KRAS / NRAS) — stage representative\n"
             "RNA expression visualisation only · NOT mutation status (mutation requires DNA-seq)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F28_driver_genes_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F29: GSE230424 PTC+HT TROP2 spatial maps =====
print("[F29] GSE230424 PTC+HT TROP2 spatial maps")
g230_paths = sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad"))
fig = plt.figure(figsize=(15, 4.4))
g = gs.GridSpec(1, 4, wspace=0.18, top=0.85, bottom=0.04, left=0.04, right=0.97)
for c, p in enumerate(g230_paths):
    sample = p.parent.name
    ax = fig.add_subplot(g[0, c])
    sub = trop_all[trop_all.sample_id == sample]
    if sub.empty or "array_row" not in sub.columns: ax.axis("off"); continue
    x = sub["array_col"].values; y = -sub["array_row"].values
    v = sub["TROP2"].values
    vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
    ax.scatter(x, y, c=v, cmap=trop_cmap, vmin=0, vmax=vmax, s=4, alpha=0.92)
    moran = df_E[df_E.sample_id == sample]["morans_I_TROP2"]
    moran = moran.iloc[0] if not moran.empty else np.nan
    ax.set_title(f"PTC+HT · {sample.split('_')[-1]}\nTROP2 Moran I = {moran:.2f}", fontsize=11, color="#962E2E", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F29  GSE230424 PTC+HT TROP2 spatial maps — niche structure weak (GSE250521 PTC와 대비)\n"
             "HT 동반 시 TROP2 niche 약해질 가능성 (sample n=4 작음, future work)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F29_TROP2_maps_GSE230424.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F30: Closure battery LOSO per-slide breakdown =====
print("[F30] Closure battery per-slide LOSO predictions")
clos_pred = pd.read_csv(ROOT/"project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz", sep="\t")
print(f"  LOSO predictions: {len(clos_pred)} rows, columns: {list(clos_pred.columns)[:8]}")
fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), sharex=False, sharey=False)
# panel A: per-slide pred vs obs
if "DM1_like_score" in clos_pred.columns and "DM1_like_score_pred" in clos_pred.columns:
    ax = axes[0,0]
    slides = sorted(clos_pred["slide"].dropna().unique()) if "slide" in clos_pred.columns else []
    if not slides and "sample" in clos_pred.columns:
        slides = sorted(clos_pred["sample"].dropna().unique())
    samp_col = "slide" if "slide" in clos_pred.columns else ("sample" if "sample" in clos_pred.columns else None)
    if samp_col:
        for sl in slides[:16]:
            sub = clos_pred[clos_pred[samp_col] == sl]
            ax.scatter(sub["DM1_like_score"], sub["DM1_like_score_pred"], s=4, alpha=0.4, label=sl[:12])
    ax.plot([-3,3],[-3,3], "--", color="black", lw=0.8)
    ax.set_xlabel("Observed DM1_like (within-sample z)")
    ax.set_ylabel("Predicted DM1_like (LOSO ResNet50 ridge)")
    ax.set_title("A. Pred vs obs per slide — diagonal = perfect, scatter = poor", fontsize=10.5)
# panel B: residual vs proliferation (still no signal)
ax = axes[0,1]
if "DM1_like_score" in clos_pred.columns and "DM1_like_score_pred" in clos_pred.columns:
    resid = (clos_pred["DM1_like_score"] - clos_pred["DM1_like_score_pred"]).values
    ax.hist(resid, bins=60, color="#7B1F2A", alpha=0.85, edgecolor="black", linewidth=0.3)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("DM1_like prediction residual (obs - pred)")
    ax.set_ylabel("Spot count")
    ax.set_title(f"B. Residual distribution (mean={np.mean(resid):.2f}, std={np.std(resid):.2f})", fontsize=10.5)
# panel C: closure battery summary metrics (re-show)
ax = axes[1,0]
clos = pd.read_csv(ROOT/"project/results/03_pathology_poc/closure_battery_metrics.tsv", sep="\t")
sec_a = clos[clos["section"].astype(str).str.startswith("A")]
labels = sec_a["experiment"].astype(str).str.replace("raw_","")
ridge = sec_a[sec_a["model"]=="ridge"]
xs = np.arange(len(ridge))
ax.bar(xs, ridge["pooled_spearman_r"].values, color="#7B1F2A", edgecolor="black", label="ridge")
ax.axhline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.6)
ax.text(0, 0.32, "useful threshold (0.3)", color="#3C6B4F", fontsize=9)
ax.axhline(0, color="black", lw=0.6)
ax.set_xticks(xs); ax.set_xticklabels(labels.iloc[:len(xs)], rotation=15, fontsize=8.5, ha="right")
ax.set_ylabel("Pooled Spearman ρ")
ax.set_ylim(-0.2, 0.5)
ax.set_title("C. Section A pooled Spearman — all below threshold", fontsize=10.5)
# panel D: section E slide-level (re-show)
ax = axes[1,1]
sec_e = clos[clos["section"].astype(str)=="E"]
labels_e = sec_e["experiment"].astype(str).str.replace("E_","").str.replace("_score_resid","").str.replace("_"," ")
vals_e = sec_e["pooled_spearman_r"].fillna(0).values
colors_e = ["#7B1F2A" if abs(v)<0.3 else "#3C6B4F" for v in vals_e]
ax.barh(np.arange(len(sec_e)), vals_e, color=colors_e, edgecolor="black")
ax.set_yticks(np.arange(len(sec_e))); ax.set_yticklabels(labels_e, fontsize=8.5)
ax.axvline(0, color="black", lw=0.6); ax.axvline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5); ax.axvline(-0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
ax.set_xlim(-0.5, 0.5)
ax.set_xlabel("Slide-level LOSO ρ")
ax.set_title("D. Slide-level LOSO breakdown — top25 mean only marginal (0.21)", fontsize=10.5)
fig.suptitle("S_F30  Closure battery LOSO per-slide breakdown (Paper 2A NO-GO 직접 evidence)\n"
             "ResNet50 224·448, multiple labels, multiple models — 모두 useful threshold 미달",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F30_closure_LOSO_breakdown.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F31: TROP2 niche-niche distance distribution (PTC+LPTC niche samples) =====
print("[F31] TROP2 niche cluster spacing")
from scipy.ndimage import label, center_of_mass
fig, ax = plt.subplots(figsize=(11.5, 5.5))
all_dists = {"PTC": [], "LPTC": []}
for stage in ["PTC","LPTC"]:
    samples = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==stage)]["sample_id"].tolist()
    for sample in samples:
        sub = trop_all[trop_all.sample_id == sample]
        if "array_row" not in sub.columns or sub.empty: continue
        thresh = sub["TROP2"].quantile(0.85)
        high = sub[sub["TROP2"] > thresh]
        if high.empty: continue
        rows = high["array_row"].astype(int).values
        cols = high["array_col"].astype(int).values
        rmin, cmin = rows.min(), cols.min()
        grid = np.zeros((rows.max()-rmin+2, cols.max()-cmin+2), dtype=int)
        for r,c in zip(rows-rmin, cols-cmin): grid[r,c] = 1
        labeled, n = label(grid)
        if n < 2: continue
        centers = center_of_mass(grid, labeled, range(1, n+1))
        # pairwise distances
        for i in range(len(centers)):
            for j in range(i+1, len(centers)):
                d = np.sqrt((centers[i][0]-centers[j][0])**2 + (centers[i][1]-centers[j][1])**2)
                all_dists[stage].append(d)
for stage, color in [("PTC","#34547A"),("LPTC","#B8893C")]:
    if all_dists[stage]:
        ax.hist(all_dists[stage], bins=30, alpha=0.7, label=f"{stage} (n={len(all_dists[stage])} pairs)",
                color=color, edgecolor="black", linewidth=0.4)
ax.set_xlabel("TROP2 niche cluster pairwise distance (Visium spot units)")
ax.set_ylabel("Pair count")
ax.set_title("S_F31  TROP2 niche cluster spacing — pairwise distance distribution\n"
             "Multi-modal = multiple niches per slide; single peak = niche boundary defined",
             fontsize=12, fontweight="bold")
ax.legend(); ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F31_niche_cluster_spacing.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F32: ATC degradation evidence (transcriptional collapse) =====
print("[F32] ATC degradation — RAI/DM1 axis collapse")
fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
# panel A: per-stage RAI_8 raw + Epithelial mean (within-sample variance reduction)
ax = axes[0]
g521_pos = trop_all[trop_all.dataset == "GSE250521"].copy()
g521_pos["row_col"] = g521_pos["sample_id"]+"_"+g521_pos["array_row"].astype(str)+"_"+g521_pos["array_col"].astype(str)
g521["row_col"] = g521["sample_id"]+"_"+g521["array_row"].astype(str)+"_"+g521["array_col"].astype(str)
mg_all = g521_pos.merge(g521[["row_col","DM1_like_score","Epithelial_score","Proliferation_score"]], on="row_col")
trop_per_sample = mg_all.groupby(["sample_id","condition"])["TROP2"].agg(["mean","std"]).reset_index() if "TROP2" in mg_all.columns else None
positions = []
for i, stage in enumerate(ORDER):
    sub = mg_all[mg_all.condition == stage]
    samp_data = sub.groupby("sample_id")["TROP2"].agg(["mean","std"]).reset_index()
    ax.scatter([i]*len(samp_data), samp_data["mean"], s=80,
               color=PAL[stage], alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
ax.set_xticks(range(4)); ax.set_xticklabels(ORDER)
ax.set_ylabel("Per-sample mean TROP2 (log-norm)")
ax.set_title("A. ATC TROP2 expression collapse\n(per-sample mean - many ATC samples drop to ~0)", fontsize=10.5)
ax.grid(axis="y", alpha=0.3)

# panel B: per-stage TROP2 std (intra-tumor variance)
ax = axes[1]
for i, stage in enumerate(ORDER):
    sub = mg_all[mg_all.condition == stage]
    samp_data = sub.groupby("sample_id")["TROP2"].agg(["mean","std"]).reset_index()
    ax.scatter([i]*len(samp_data), samp_data["std"], s=80,
               color=PAL[stage], alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
ax.set_xticks(range(4)); ax.set_xticklabels(ORDER)
ax.set_ylabel("Per-sample TROP2 std (intra-tumor variance)")
ax.set_title("B. Intra-tumor TROP2 variance\n(PTC+LPTC = high variance niches; ATC = collapse)", fontsize=10.5)
ax.grid(axis="y", alpha=0.3)

# panel C: per-stage Moran I + per-stage % spots non-zero
ax = axes[2]
e = df_E[df_E.dataset == "GSE250521"]
xs = np.arange(len(ORDER))
moran_means = [e[e.condition==s]["morans_I_TROP2"].mean() for s in ORDER]
nonzero_pcts = []
for s in ORDER:
    samples_s = e[e.condition==s]["sample_id"].tolist()
    sub_t = trop_all[trop_all.sample_id.isin(samples_s)]
    nonzero_pcts.append(((sub_t["TROP2"] > 0.001).sum() / len(sub_t) * 100) if len(sub_t) else 0)
ax2 = ax.twinx()
ax.bar(xs - 0.18, moran_means, 0.36, color="#7B1F2A", edgecolor="black", label="mean Moran's I")
ax2.bar(xs + 0.18, nonzero_pcts, 0.36, color="#34547A", edgecolor="black", label="% non-zero TROP2 spots")
ax.set_xticks(xs); ax.set_xticklabels(ORDER)
ax.set_ylabel("Mean Moran's I", color="#7B1F2A")
ax2.set_ylabel("% non-zero TROP2 spots", color="#34547A")
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.5)
ax.set_title("C. ATC niche collapse + expression sparsity\n(Moran I + non-zero spot fraction)", fontsize=10.5)
ax.legend(loc="upper left", fontsize=9); ax2.legend(loc="upper right", fontsize=9)
fig.suptitle("S_F32  ATC dedifferentiation evidence — TROP2 niche structure collapse in 3/4 ATC samples\n"
             "PTC+LPTC = niche-organized + high expression; ATC = either collapsed expression or scattered remaining",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F32_ATC_degradation.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 6 new figures (S_F27 ~ S_F32) saved")
