#!/usr/bin/env python3
"""Spatial v9 — 8 more (S_F65 ~ S_F72)."""
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
from sklearn.decomposition import PCA
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
tls_cmap = LinearSegmentedColormap.from_list("tls",["#fffaf2","#a8d8c8","#1abc9c","#0F1A2E"],N=256)

def score_h5ad(p, gs_list):
    a = ad.read_h5ad(p)
    rvar = a.raw.var_names.astype(str) if a.raw is not None else a.var_names.astype(str)
    X = a.raw.X if a.raw is not None else a.X
    idx = [list(rvar.values).index(g) for g in gs_list if g in rvar.values]
    if not idx: return None, None
    sub = X[:, idx]
    if hasattr(sub,"toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", np.ones(a.n_obs))).flatten()
    tot = np.where(tot==0,1,tot)
    return np.log1p(sub/tot[:,None]*1e4), (a.obs["array_row"].values, a.obs["array_col"].values) if "array_row" in a.obs.columns else (None, None)

print("[load]")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
trop_all = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
df_C = pd.read_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t")
df_moran = pd.read_csv(OUT/"spatial_F40_per_sample_per_axis_morans.tsv", sep="\t")
df_combined = pd.read_csv(OUT/"spatial_F64_28sample_multiaxis_morans.tsv", sep="\t")

# ===== S_F65: TLS individual gene maps (CXCL13 / CCL19 / CCL21 / CXCR5) on P3 ★ =====
print("[F65] TLS individual gene maps on P3")
TLS_GENES_INDIV = ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"]
g230_paths = sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad"))
sample_p3 = "GSM7221917_P3"
p3_path = next((p for p in g230_paths if p.parent.name == sample_p3), None)
if p3_path:
    fig = plt.figure(figsize=(15, 9))
    g = gs.GridSpec(2, 3, hspace=0.35, wspace=0.18, top=0.92, bottom=0.04, left=0.04, right=0.97)
    for i, gene in enumerate(TLS_GENES_INDIV):
        ax = fig.add_subplot(g[i // 3, i % 3])
        scores, coords = score_h5ad(p3_path, [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows, cols = coords
        if rows is None: ax.axis("off"); continue
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols), -np.array(rows), c=v, cmap=tls_cmap, vmin=0, vmax=vmax, s=4, alpha=0.92)
        nonzero_pct = (v > 0.01).sum() / len(v) * 100
        ax.set_title(f"{gene}\nnon-zero in {nonzero_pct:.1f}% spots, max={v.max():.1f}", fontsize=10.5, color="#1abc9c", fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
    fig.suptitle(f"S_F65  TLS individual gene maps on {sample_p3.split('_')[-1]} (best TLS niche, Moran I=0.81)\n"
                 "TLS의 6 gene 모두 같은 region에 collocated → tertiary lymphoid structure 형성 직접 evidence",
                 fontsize=12, fontweight="bold")
    fig.savefig(ASSETS/"S_F65_TLS_individual_genes_P3.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F66: 28-sample similarity matrix (Moran I correlation) =====
print("[F66] 28-sample similarity matrix")
# Use available axes per sample, compute pairwise correlation of Moran I profiles
axes_for_sim = ["TROP2","DM1_like","Epithelial","Proliferation","CAF_ECM","EMT","Hypoxia"]
mat_for_sim = np.full((len(df_combined), len(axes_for_sim)), np.nan)
for i in range(len(df_combined)):
    for j, a in enumerate(axes_for_sim):
        v = df_combined.iloc[i].get(a, np.nan)
        if not pd.isna(v): mat_for_sim[i, j] = float(v)
sample_labels = [f"{r['condition']}·{r['sample'].split('_')[-1]}" for _, r in df_combined.iterrows()]
# Compute pairwise correlation
sim = np.full((len(df_combined), len(df_combined)), np.nan)
for i in range(len(df_combined)):
    for j in range(len(df_combined)):
        a, b = mat_for_sim[i], mat_for_sim[j]
        valid = ~(np.isnan(a) | np.isnan(b))
        if valid.sum() < 3: continue
        if a[valid].std() == 0 or b[valid].std() == 0: continue
        sim[i, j] = np.corrcoef(a[valid], b[valid])[0, 1]
fig, ax = plt.subplots(figsize=(13, 12))
im = ax.imshow(sim, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(df_combined))); ax.set_xticklabels(sample_labels, fontsize=7.5, rotation=70, ha="right")
ax.set_yticks(range(len(df_combined))); ax.set_yticklabels(sample_labels, fontsize=7.5)
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Moran's I profile correlation")
ax.set_title("S_F66  28-sample similarity matrix — pairwise correlation of axis-Moran's I profiles\n"
             "Block structure = samples with similar niche organization patterns",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F66_28_sample_similarity.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F67: PCA of samples in axis-Moran I space =====
print("[F67] Sample PCA in Moran I space")
# Use samples with all 7 axes available
valid_idx = ~np.isnan(mat_for_sim).any(axis=1)
if valid_idx.sum() >= 5:
    X_pca = mat_for_sim[valid_idx]
    pca = PCA(n_components=2)
    Y = pca.fit_transform(X_pca)
    fig, ax = plt.subplots(figsize=(10, 7))
    sub_combined = df_combined[valid_idx].reset_index(drop=True)
    for cond, color in PAL.items():
        mask = sub_combined["condition"] == cond
        if not mask.any(): continue
        ax.scatter(Y[mask, 0], Y[mask, 1], s=180, color=color, alpha=0.85,
                   edgecolor="black", linewidth=0.7, label=cond, zorder=3)
        for k, m in enumerate(mask):
            if m:
                ax.annotate(sub_combined.iloc[k]["sample"].split("_")[-1], (Y[k,0], Y[k,1]),
                            fontsize=8, xytext=(6,6), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
    ax.axhline(0, color="black", lw=0.5); ax.axvline(0, color="black", lw=0.5)
    ax.set_title("S_F67  Sample PCA in 7-axis Moran's I space (GSE250521 only, n=16)\n"
                 "PTC+LPTC samples cluster together (high TROP2/Epi/Prol niche); PT/ATC separate",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="best", fontsize=10, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ASSETS/"S_F67_sample_PCA_morans.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F68: Per-sample axis trajectory line plot (16 GSE250521) =====
print("[F68] Per-sample axis trajectory")
fig, ax = plt.subplots(figsize=(13, 6.5))
axes_traj = ["RAI_8_score","DM1_like_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score","TROP2"]
xs = np.arange(len(axes_traj))
for sample in df_moran["sample"].unique():
    sub = df_moran[df_moran["sample"] == sample]
    if sub.empty: continue
    cond = sub["condition"].iloc[0]
    vals = [sub.iloc[0].get(a, np.nan) for a in axes_traj]
    ax.plot(xs, vals, "o-", color=PAL.get(cond,"#999"), alpha=0.85, lw=1.4, ms=7,
            label=cond if sample == df_moran[df_moran.condition == cond]["sample"].iloc[0] else None)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1.2, alpha=0.6)
ax.set_xticks(xs); ax.set_xticklabels([a.replace("_score","") for a in axes_traj], fontsize=10, rotation=20, ha="right")
ax.set_ylabel("Moran's I per axis (per sample)")
ax.set_title("S_F68  Per-sample axis trajectory of Moran's I (16 GSE250521 samples × 8 axes)\n"
             "각 line = 1 sample; PTC+LPTC samples (blue+gold) consistently elevated across multiple axes",
             fontsize=12, fontweight="bold")
ax.legend(loc="upper right", fontsize=9, ncol=2)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F68_per_sample_axis_trajectory.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F69: Stage representative 3-panel (DM1 + Epi + Prol) =====
print("[F69] Stage representative 3-axis spatial")
rep = {}
for s in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==s)]
    if not sub.empty:
        rep[s] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]
fig = plt.figure(figsize=(15, 11))
g = gs.GridSpec(3, 4, hspace=0.32, wspace=0.18, top=0.94, bottom=0.04, left=0.05, right=0.97)
axes_view = [
    ("DM1_like_score","DM1_like","#7B1F2A"),
    ("Epithelial_score","Epithelial","#3C6B4F"),
    ("Proliferation_score","Proliferation","#34547A"),
]
for r, (col, label, color) in enumerate(axes_view):
    cmap_use = LinearSegmentedColormap.from_list("c",["#fffaf2","#cccccc",color],N=256)
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample: continue
        ax = fig.add_subplot(g[r, c])
        sub = g521[g521.sample_id == sample]
        if sub.empty or "array_row" not in sub.columns: ax.axis("off"); continue
        x = sub["array_col"].values; y = -sub["array_row"].values
        v = sub[col].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
        ax.scatter(x, y, c=v, cmap=cmap_use, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(label, fontsize=12, fontweight="bold", color=color, rotation=0, labelpad=35, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F69  Stage representative — DM1 / Epithelial / Proliferation 3-axis spatial maps",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F69_stage_3axis_spatial.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F70: TROP2 niche compactness (size vs spatial extent ratio) =====
print("[F70] TROP2 niche compactness")
def visium_hex_clusters_with_extent(rows, cols, mask):
    pos_to_idx = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols)) if mask[i]}
    parent = {p: p for p in pos_to_idx.values()}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (r, c), i in pos_to_idx.items():
        for dr, dc in [(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]:
            j = pos_to_idx.get((r+dr, c+dc))
            if j is not None:
                ra, rb = find(i), find(j)
                if ra != rb: parent[ra] = rb
    cluster_to_pts = {}
    for (r, c), i in pos_to_idx.items():
        root = find(i)
        cluster_to_pts.setdefault(root, []).append((r, c))
    return cluster_to_pts

niche_compact = []
for sample in df_E["sample_id"].unique():
    sub = trop_all[trop_all.sample_id == sample]
    cond = sub["condition"].iloc[0]
    if cond not in ("PTC","LPTC"): continue
    if sub.empty or "array_row" not in sub.columns: continue
    rows_arr = sub["array_row"].astype(int).values
    cols_arr = sub["array_col"].astype(int).values
    threshold = sub["TROP2"].quantile(0.85)
    mask = (sub["TROP2"] > threshold).values
    if mask.sum() < 5: continue
    clusters = visium_hex_clusters_with_extent(rows_arr, cols_arr, mask)
    for root, pts in clusters.items():
        if len(pts) < 3: continue  # too small for shape
        rs = [p[0] for p in pts]; cs = [p[1] for p in pts]
        bbox_h = (max(rs) - min(rs)) + 1
        bbox_w = (max(cs) - min(cs)) + 2  # account for visium hex 2-step col
        bbox_area = bbox_h * bbox_w / 2  # hex grid effective area
        compactness = len(pts) / bbox_area if bbox_area > 0 else 0
        niche_compact.append({"sample":sample, "stage":cond, "size":len(pts),
                              "bbox_area":bbox_area, "compactness":compactness})
df_compact = pd.DataFrame(niche_compact)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
ax = axes[0]
for stage, color in [("PTC","#34547A"),("LPTC","#B8893C")]:
    sub = df_compact[df_compact.stage == stage]
    if sub.empty: continue
    ax.scatter(sub["size"], sub["compactness"], s=80, color=color, alpha=0.7,
               edgecolor="black", linewidth=0.5, label=f"{stage} (n={len(sub)})")
ax.set_xlabel("Niche cluster size (n_spots)"); ax.set_ylabel("Compactness (n_spots / bbox area)")
ax.set_xscale("log")
ax.axhline(0.5, ls=":", color="gray", lw=0.6); ax.text(2, 0.55, "0.5 = compact threshold", fontsize=9, color="gray")
ax.set_title("A. Cluster compactness vs size", fontsize=11, fontweight="bold")
ax.legend(); ax.grid(alpha=0.3)
ax = axes[1]
for stage, color in [("PTC","#34547A"),("LPTC","#B8893C")]:
    sub = df_compact[df_compact.stage == stage]
    if sub.empty: continue
    ax.hist(sub["compactness"].values, bins=30, color=color, alpha=0.7, edgecolor="black", linewidth=0.4, label=stage)
ax.set_xlabel("Compactness"); ax.set_ylabel("Cluster count")
ax.set_title("B. Compactness distribution", fontsize=11, fontweight="bold")
ax.legend(); ax.grid(axis="y", alpha=0.3)
fig.suptitle("S_F70  TROP2 niche compactness analysis (Visium hex 6-neighbor clusters, PTC+LPTC niche samples)",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F70_niche_compactness.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F71: Slide × axis correlation heatmap (within slide cross-axis) =====
print("[F71] Within-slide cross-axis Spearman heatmap")
# For each GSE250521 slide, compute spot-level Spearman ρ between TROP2 and other axes
trop_g521 = trop_all[trop_all.dataset == "GSE250521"].copy()
trop_g521["row_col"] = trop_g521["sample_id"]+"_"+trop_g521["array_row"].astype(str)+"_"+trop_g521["array_col"].astype(str)
g521["row_col"] = g521["sample_id"]+"_"+g521["array_row"].astype(str)+"_"+g521["array_col"].astype(str)
mg = trop_g521.merge(g521[["row_col","DM1_like_score","RAI_8_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]], on="row_col")
axes_for = ["DM1_like_score","RAI_8_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]
samples_g521 = sorted(mg["sample_id"].unique())
mat_corr = np.full((len(samples_g521), len(axes_for)), np.nan)
for i, sample in enumerate(samples_g521):
    sub = mg[mg.sample_id == sample]
    for j, ax_name in enumerate(axes_for):
        if sub[ax_name].nunique() > 1:
            rho, _ = spearmanr(sub["TROP2"], sub[ax_name], nan_policy="omit")
            mat_corr[i, j] = rho
fig, ax = plt.subplots(figsize=(13, 7))
im = ax.imshow(mat_corr, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
for i in range(len(samples_g521)):
    for j in range(len(axes_for)):
        v = mat_corr[i,j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                    color="white" if abs(v)>0.4 else "black", fontsize=8.5, fontweight="bold")
ax.set_xticks(range(len(axes_for))); ax.set_xticklabels([a.replace("_score","") for a in axes_for], fontsize=10, family="monospace", rotation=20, ha="right")
ax.set_yticks(range(len(samples_g521)))
def stage_of(s):
    if "_N-" in s: return "PT"
    if "PTC-" in s and "L" not in s: return "PTC"
    if "LPTC-" in s: return "LPTC"
    if "ATC-" in s: return "ATC"
    return "?"
ax.set_yticklabels([f"{stage_of(s)} · {s.split('_')[-1]}" for s in samples_g521], fontsize=9, family="monospace")
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.04, label="Spot-level Spearman ρ (TROP2 vs axis)")
ax.set_title("S_F71  Within-slide cross-axis Spearman ρ — 16 GSE250521 slides × 7 axes (TROP2 reference)\n"
             "TROP2 × Epithelial 일관되게 +; TROP2 × DM1 stage-dependent; CAF/ECM 음의 상관 (stromal separation)",
             fontsize=11.5, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F71_within_slide_cross_axis_rho.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F72: ULTIMATE summary infographic =====
print("[F72] Ultimate summary infographic")
fig = plt.figure(figsize=(18, 12))
g = gs.GridSpec(3, 3, hspace=0.4, wspace=0.25, top=0.93, bottom=0.04, left=0.05, right=0.97)

# top-left: TROP2 niche per condition (boxplot)
ax = fig.add_subplot(g[0, 0])
data_box = []; labels_box = []
for c in ["PT","PTC","LPTC","ATC","PTC_HT","HT","GD","CONTROL"]:
    sub = df_E[df_E.condition == c]
    if not sub.empty: data_box.append(sub["morans_I_TROP2"].values); labels_box.append(c)
bp = ax.boxplot(data_box, labels=labels_box, patch_artist=True, showfliers=False)
for patch, c in zip(bp["boxes"], labels_box):
    patch.set_facecolor(PAL.get(c,"#999")); patch.set_alpha(0.8)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1.4)
ax.set_ylabel("TROP2 Moran's I")
ax.set_title("Q1. TROP2 niche tumor-specific", fontsize=11, color="#7B1F2A", fontweight="bold")
ax.tick_params(axis='x', rotation=30); ax.grid(axis="y", alpha=0.3)

# top-mid: HT TLS heatmap
ax = fig.add_subplot(g[0, 1])
samps = sorted(df_C["sample_id"].unique())
markers = ["HLA_II","B_cell","TLS","IGHV_AICDA"]
mat_ht = np.full((len(samps), len(markers)), np.nan)
for i, s in enumerate(samps):
    for j, m in enumerate(markers):
        sub = df_C[(df_C.sample_id==s) & (df_C.gene_set==m)]
        if not sub.empty: mat_ht[i,j] = sub["morans_I"].iloc[0]
im = ax.imshow(mat_ht, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.9)
for i in range(len(samps)):
    for j in range(len(markers)):
        ax.text(j, i, f"{mat_ht[i,j]:.2f}", ha="center", va="center",
                color="white" if mat_ht[i,j]>0.55 else "black", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(markers))); ax.set_xticklabels(markers, fontsize=9, family="monospace")
ax.set_yticks(range(len(samps))); ax.set_yticklabels([s.split("_")[-1] for s in samps], fontsize=9, family="monospace")
ax.set_title("Q2. HT TLS niche", fontsize=11, color="#1abc9c", fontweight="bold")

# top-right: stage trajectory
ax = fig.add_subplot(g[0, 2])
for ax_name, color in [("TROP2","#7B1F2A"),("DM1_like_score","#34547A"),("Epithelial_score","#3C6B4F"),("Proliferation_score","#B8893C")]:
    if ax_name not in df_moran.columns: continue
    means = []; stds = []
    for cond in ORDER:
        sub = df_moran[df_moran.condition == cond][ax_name].dropna()
        means.append(sub.mean() if len(sub) else np.nan)
        stds.append(sub.std() if len(sub) > 1 else 0)
    ax.errorbar(range(len(ORDER)), means, yerr=stds, marker="o", lw=2.4,
                color=color, label=ax_name.replace("_score",""), capsize=5, alpha=0.9, ms=9)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.5)
ax.set_ylabel("Moran's I")
ax.set_title("Q3. Stage trajectory", fontsize=11, color="#34547A", fontweight="bold")
ax.legend(fontsize=9, loc="best"); ax.grid(axis="y", alpha=0.3)

# middle-left: closure NO-GO
ax = fig.add_subplot(g[1, 0])
clos_metric = pd.read_csv(ROOT/"project/results/03_pathology_poc/closure_battery_metrics.tsv", sep="\t")
sec_a = clos_metric[(clos_metric["section"].astype(str).str.startswith("A")) & (clos_metric["model"]=="ridge")]
labels_c = sec_a["experiment"].astype(str).str.replace("raw_","")
ax.bar(range(len(sec_a)), sec_a["pooled_spearman_r"].values, color="#7B1F2A", edgecolor="black")
ax.axhline(0.3, ls="--", color="#3C6B4F", lw=1.5); ax.text(0, 0.32, "useful threshold", fontsize=9, color="#3C6B4F")
ax.axhline(0, color="black", lw=0.5)
ax.set_xticks(range(len(sec_a))); ax.set_xticklabels(labels_c, fontsize=9, rotation=15, ha="right")
ax.set_ylabel("LOSO ρ"); ax.set_ylim(-0.1, 0.5)
ax.set_title("Q4. H&E NO-GO", fontsize=11, color="#52525a", fontweight="bold")

# middle-mid: TROP2 vs DM1 Moran scatter (S_F51 condensed)
ax = fig.add_subplot(g[1, 1])
for cond, color in [("PT","#3C6B4F"),("PTC","#34547A"),("LPTC","#B8893C"),("ATC","#7B1F2A")]:
    sub = df_moran[df_moran.condition == cond]
    if sub.empty: continue
    ax.scatter(sub["DM1_like_score"], sub["TROP2"], s=120, color=color,
               alpha=0.85, edgecolor="black", linewidth=0.6, label=f"{cond}")
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.4)
ax.axvline(0.25, ls="--", color="#34547A", lw=1, alpha=0.4)
ax.axhline(0, color="black", lw=0.4); ax.axvline(0, color="black", lw=0.4)
ax.set_xlabel("DM1_like Moran's I"); ax.set_ylabel("TROP2 Moran's I")
ax.set_title("Q5. TROP2 × DM1 niche", fontsize=11, color="#7B1F2A", fontweight="bold")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# middle-right: cohort coverage
ax = fig.add_subplot(g[1, 2])
cohort_summary = df_E.groupby("dataset").agg(n_samp=("sample_id","count"), n_spots=("n_spots","sum")).reset_index()
ax.barh(range(len(cohort_summary)), cohort_summary["n_spots"].values,
        color=["#7B1F2A","#1abc9c","#3C6B4F"], edgecolor="black")
ax.set_yticks(range(len(cohort_summary)))
ax.set_yticklabels([f"{r['dataset']}\n(n_samp={int(r['n_samp'])})" for _,r in cohort_summary.iterrows()], fontsize=10)
ax.set_xlabel("Total spots")
ax.set_title("Cohort coverage", fontsize=11, fontweight="bold")
for i, v in enumerate(cohort_summary["n_spots"].values):
    ax.text(v + 1500, i, f"{v:,}", va="center", fontsize=10, fontweight="bold")

# bottom row: 3 large headline numbers
ax = fig.add_subplot(g[2, 0])
ax.axis("off")
ax.text(0.5, 0.7, "8/8", ha="center", va="center", fontsize=64, fontweight="bold", color="#7B1F2A", transform=ax.transAxes)
ax.text(0.5, 0.3, "PTC + LPTC samples\nTROP2 niche-organized", ha="center", va="center", fontsize=14, transform=ax.transAxes)

ax = fig.add_subplot(g[2, 1])
ax.axis("off")
ax.text(0.5, 0.7, "0.84", ha="center", va="center", fontsize=64, fontweight="bold", color="#1abc9c", transform=ax.transAxes)
ax.text(0.5, 0.3, "Strongest HT TLS Moran's I\n(P3 IGHV/AICDA)", ha="center", va="center", fontsize=14, transform=ax.transAxes)

ax = fig.add_subplot(g[2, 2])
ax.axis("off")
ax.text(0.5, 0.7, "ρ=0.06", ha="center", va="center", fontsize=64, fontweight="bold", color="#52525a", transform=ax.transAxes)
ax.text(0.5, 0.3, "H&E → DM1 LOSO\n(NO-GO, below 0.3)", ha="center", va="center", fontsize=14, transform=ax.transAxes)

fig.suptitle("S_F72  ULTIMATE SPATIAL SUMMARY — 5 panels + 3 headline numbers (Spatial full package mega graphical abstract)",
             fontsize=15, fontweight="bold", y=0.99)
fig.savefig(ASSETS/"S_F72_ultimate_summary.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 8 new figures (S_F65 ~ S_F72) saved")
