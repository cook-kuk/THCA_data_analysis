#!/usr/bin/env python3
"""Spatial full package v3 — 8 more figures (S_F19~S_F26)."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
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

PAL = {"PT":"#3C6B4F","PTC":"#34547A","LPTC":"#B8893C","ATC":"#7B1F2A",
       "PTC_HT":"#962E2E","CONTROL":"#52525a","HT":"#34547A","GD":"#B8893C"}
ORDER = ["PT","PTC","LPTC","ATC"]
trop_cmap = LinearSegmentedColormap.from_list("trop2",["#fffaf2","#ffd58a","#d97c2e","#7B1F2A"],N=256)
dm_cmap   = LinearSegmentedColormap.from_list("dm",["#34547A","#fffaf2","#7B1F2A"],N=256)
prol_cmap = LinearSegmentedColormap.from_list("prol",["#fffaf2","#a8c8e2","#34547A","#0F1A2E"],N=256)

print("[load] data tables")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
trop_all = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
df_C = pd.read_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t")
ext_files = sorted((ROOT/"project_external_st/results/scores").glob("*_spot_scores.tsv.gz"))
ext = pd.concat([pd.read_csv(f, sep="\t") for f in ext_files], axis=0, ignore_index=True)

# Build by-stage sample mapping
g521_samples = g521["sample_id"].unique()
by_stage = {s:[] for s in ORDER}
for sid in g521_samples:
    if "_N-" in sid: by_stage["PT"].append(sid)
    elif "PTC-" in sid and "L" not in sid: by_stage["PTC"].append(sid)
    elif "LPTC-" in sid: by_stage["LPTC"].append(sid)
    elif "ATC-" in sid: by_stage["ATC"].append(sid)

# representative per stage (max Moran I)
rep = {}
for s in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==s)]
    if not sub.empty:
        rep[s] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]

# ===== S_F19: GSE250521 16-slide DM1_like maps (parallel to S_F7) =====
print("[F19] GSE250521 16-slide DM1 maps")
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.94, bottom=0.04, left=0.04, right=0.97)
for r, stage in enumerate(ORDER):
    samples = sorted(by_stage[stage])[:4]
    for c, sample in enumerate(samples):
        ax = fig.add_subplot(g[r, c])
        sub = g521[g521.sample_id == sample]
        if "array_row" not in sub.columns or sub.empty:
            ax.text(0.5,0.5,"no coords"); ax.axis("off"); continue
        x = sub["array_col"].values; y = -sub["array_row"].values
        v = sub["DM1_like_score"].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95))) if not np.isnan(v).all() else 1
        ax.scatter(x, y, c=v, cmap=dm_cmap, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        title_color = PAL[stage]
        ax.set_title(f"{stage}  ·  {sample.split('_')[-1]}", fontsize=10, color=title_color, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal", adjustable="box")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F19  DM1_like spatial maps — GSE250521 16 slides (parallel to S_F7 TROP2)\n"
             "blue = RAI-high (differentiated); red = DM1-high (dark matter); within-sample z normalization",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F19_DM1_maps_GSE250521.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F20: GSE230424 4-slide DM1 maps (DM1 axis on PTC+HT) =====
print("[F20] GSE230424 4-slide DM1 maps")
g230 = ext[ext.dataset == "GSE230424"]
g230_samples = sorted(g230["sample_id"].unique())
fig = plt.figure(figsize=(15, 4.4))
g = gs.GridSpec(1, 4, wspace=0.18, top=0.86, bottom=0.04, left=0.04, right=0.97)
for c, sample in enumerate(g230_samples):
    ax = fig.add_subplot(g[0, c])
    sub = g230[g230.sample_id == sample]
    x = sub["array_col"].values; y = -sub["array_row"].values
    v = sub["DM1_like_score_raw"].values
    vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95))) if not np.isnan(v).all() else 1
    ax.scatter(x, y, c=v, cmap=dm_cmap, vmin=-vmax, vmax=vmax, s=4, alpha=0.92)
    ax.set_title(f"PTC+HT · {sample.split('_')[-1]}", fontsize=11, color="#962E2E", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal", adjustable="box")
fig.suptitle("S_F20  DM1_like spatial maps — GSE230424 PTC+HT 4 slides\n"
             "Same axis as S_F19 (Paper 1) on Hashimoto-overlap PTC slides — 같은 region에 differentiation/dark-matter polarization 보임",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F20_DM1_maps_GSE230424.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F21: TROP2 × microenvironment (Proliferation / Hypoxia / CAF/ECM) =====
print("[F21] TROP2 × microenvironment scatter")
trop_g521 = trop_all[trop_all.dataset == "GSE250521"].copy()
trop_g521["row_col_id"] = trop_g521["sample_id"]+"_"+trop_g521["array_row"].astype(str)+"_"+trop_g521["array_col"].astype(str)
g521["row_col_id"] = g521["sample_id"]+"_"+g521["array_row"].astype(str)+"_"+g521["array_col"].astype(str)
mg = trop_g521.merge(g521[["row_col_id","Proliferation_score","Hypoxia_score","CAF_ECM_score","EMT_score"]], on="row_col_id", how="left")
fig, axes = plt.subplots(1, 4, figsize=(16, 4.2), sharey=True)
ENV = [("Proliferation_score","Proliferation","#34547A"),
       ("Hypoxia_score","Hypoxia","#7B1F2A"),
       ("CAF_ECM_score","CAF/ECM","#3C6B4F"),
       ("EMT_score","EMT","#B8893C")]
# Use only PTC+LPTC niche samples (where TROP2 is meaningful)
niche_samples = list(by_stage["PTC"]) + list(by_stage["LPTC"])
mg_niche = mg[mg.sample_id.isin(niche_samples)]
for ax, (col, label, color) in zip(axes, ENV):
    if col not in mg_niche.columns or mg_niche[col].isna().all():
        ax.text(0.5,0.5,"no data"); ax.axis("off"); continue
    rho, p = spearmanr(mg_niche["TROP2"], mg_niche[col], nan_policy="omit")
    h = ax.hexbin(mg_niche[col], mg_niche["TROP2"], gridsize=40, cmap="OrRd", mincnt=2)
    ax.set_title(f"{label}\nrho = {rho:.3f}, p = {p:.1e}", fontsize=11, color=color, fontweight="bold")
    ax.set_xlabel(f"{label} score (z)")
    if ax is axes[0]: ax.set_ylabel("TROP2 (log-norm)")
    ax.axhline(0, ls=":", lw=0.5, color="gray"); ax.axvline(0, ls=":", lw=0.5, color="gray")
fig.suptitle("S_F21  TROP2 × tumor microenvironment axes (PTC+LPTC niche samples, n=8)\n"
             "TROP2-high spot이 어떤 microenvironment niche와 co-localize하는가",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F21_TROP2_microenv.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F22: Per-stage Moran I distribution boxplot =====
print("[F22] Per-stage Moran I boxplot")
fig, ax = plt.subplots(figsize=(11, 5.5))
groups = []
for stage in ["PT","PTC","LPTC","ATC","PTC_HT","CONTROL","HT","GD"]:
    sub = df_E[df_E.condition == stage]
    if not sub.empty:
        groups.append((stage, sub["morans_I_TROP2"].values))
positions = list(range(len(groups)))
data = [g[1] for g in groups]
bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
for patch, (st, _) in zip(bp["boxes"], groups):
    patch.set_facecolor(PAL.get(st,"#999")); patch.set_alpha(0.8); patch.set_edgecolor("black")
# overlay points
for i, (st, vals) in enumerate(groups):
    jitter = np.random.RandomState(42+i).uniform(-0.08, 0.08, len(vals))
    ax.scatter(np.full_like(vals, i)+jitter, vals, s=42, color="black", alpha=0.85, zorder=3, edgecolor="white", linewidth=0.5)
ax.set_xticks(positions)
ax.set_xticklabels([f"{st}\n(n={len(v)})" for st, v in groups], fontsize=10)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1.4, alpha=0.7)
ax.text(0.1, 0.27, "niche threshold (0.25)", color="#7B1F2A", fontsize=10, fontweight="bold")
ax.axhline(0, color="black", lw=0.6)
ax.set_ylabel("TROP2 spatial autocorrelation (Moran's I)")
ax.set_title("S_F22  TROP2 Moran's I per condition — boxplot + scatter\n"
             "PTC+LPTC clearly separated above niche threshold; 다른 모든 condition below",
             fontsize=12, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F22_morans_per_condition_box.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F23: TROP2 expression histogram per condition =====
print("[F23] TROP2 expression distribution")
fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey=False)
conds = ["PT","PTC","LPTC","ATC","PTC_HT","CONTROL","HT","GD"]
for ax, cond in zip(axes.flat, conds):
    sub = trop_all[trop_all.condition == cond]
    n_samples = sub["sample_id"].nunique()
    if sub.empty:
        ax.text(0.5,0.5,"no data"); ax.axis("off"); continue
    nonzero = sub[sub["TROP2"] > 0.001]["TROP2"].values
    pct_nonzero = len(nonzero) / len(sub) * 100 if len(sub) else 0
    if len(nonzero) > 0:
        ax.hist(nonzero, bins=40, color=PAL.get(cond,"#999"), alpha=0.85, edgecolor="black", linewidth=0.3)
    ax.set_title(f"{cond}  ·  n_samp={n_samples}  ·  n_spots={len(sub):,}\nNon-zero TROP2: {pct_nonzero:.1f}% of spots",
                 fontsize=10, color=PAL.get(cond,"#333"))
    ax.set_xlabel("TROP2 (log-norm, non-zero only)")
    if ax in axes[:,0]: ax.set_ylabel("Spot count")
    ax.set_xlim(0, 14)
    ax.grid(axis="y", alpha=0.3)
fig.suptitle("S_F23  TROP2 expression distribution per condition — non-zero spot histograms\n"
             "PTC/LPTC: 7-10% spots non-zero, distribution shifted higher · PT/HT/GD: < 5% non-zero",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F23_TROP2_distribution_hist.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F24: HT TLS pairwise co-localization heatmap =====
print("[F24] HT TLS pairwise correlation per slide (immune set co-organization)")
import anndata as ad
EXTH = ROOT/"project_external_st/data/processed"
HLA_II = ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
B_CELL = ["CD19","MS4A1","CD79A","CD79B","SDC1","JCHAIN"]
TLS    = ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"]
IGHV   = ["AICDA","IGHM","IGHG1","IGHA1","IGKC","IGLC2"]

def score_h5ad(p, gs_list):
    a = ad.read_h5ad(p)
    rvar = a.raw.var_names.astype(str) if a.raw is not None else a.var_names.astype(str)
    X = a.raw.X if a.raw is not None else a.X
    idx = [list(rvar.values).index(g) for g in gs_list if g in rvar.values]
    if not idx: return None
    sub = X[:, idx]
    if hasattr(sub,"toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", np.ones(a.n_obs))).flatten()
    tot = np.where(tot==0,1,tot)
    return np.log1p(sub/tot[:,None]*1e4).mean(axis=1)

g230_paths = sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad"))
fig = plt.figure(figsize=(15, 4.4))
g = gs.GridSpec(1, 4, wspace=0.32, top=0.86, bottom=0.05, left=0.05, right=0.97)
for c, p in enumerate(g230_paths):
    sample = p.parent.name
    ax = fig.add_subplot(g[0, c])
    scores = {}
    for nm, gs_list in [("HLA-II",HLA_II),("B-cell",B_CELL),("TLS",TLS),("IGHV/AICDA",IGHV)]:
        sc = score_h5ad(p, gs_list)
        if sc is not None: scores[nm] = sc
    if not scores:
        ax.text(0.5,0.5,"no data"); ax.axis("off"); continue
    df_sc = pd.DataFrame(scores)
    cor = df_sc.corr(method="spearman").values
    im = ax.imshow(cor, vmin=-1, vmax=1, cmap="RdBu_r")
    for i in range(cor.shape[0]):
        for j in range(cor.shape[1]):
            ax.text(j, i, f"{cor[i,j]:.2f}", ha="center", va="center",
                    color="white" if abs(cor[i,j])>0.5 else "black", fontsize=9.5, fontweight="bold")
    ax.set_xticks(range(4)); ax.set_xticklabels(list(scores.keys()), fontsize=9, rotation=20, ha="right")
    ax.set_yticks(range(4)); ax.set_yticklabels(list(scores.keys()), fontsize=9)
    ax.set_title(f"PTC+HT · {sample.split('_')[-1]}", fontsize=11, color="#962E2E", fontweight="bold")
fig.suptitle("S_F24  GSE230424 PTC+HT — HLA-II × B-cell × TLS × IGHV/AICDA spot-level co-localization\n"
             "각 cell = 같은 slide의 spot-level Spearman ρ — 모두 양의 강한 상관 = 같은 spot에 immune signature collocated",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F24_HT_immune_colocalization.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F25: TROP2 niche size distribution =====
print("[F25] TROP2 niche cluster size distribution")
fig, ax = plt.subplots(figsize=(11, 5.5))
# For each PTC+LPTC niche sample, count contiguous clusters of TROP2-high spots
from scipy.ndimage import label
niche_size_data = []
for sample in sorted(set(by_stage["PTC"]+by_stage["LPTC"])):
    sub = trop_all[trop_all.sample_id == sample]
    if "array_row" not in sub.columns or sub.empty: continue
    thresh = sub["TROP2"].quantile(0.85)
    high = sub[sub["TROP2"] > thresh]
    if high.empty: continue
    rows = high["array_row"].astype(int).values
    cols = high["array_col"].astype(int).values
    if len(rows) == 0: continue
    rmin, cmin = rows.min(), cols.min()
    grid = np.zeros((rows.max()-rmin+2, cols.max()-cmin+2), dtype=int)
    for r,c in zip(rows-rmin, cols-cmin): grid[r,c] = 1
    labeled, n = label(grid)
    sizes = [int((labeled == i).sum()) for i in range(1, n+1)]
    for sz in sizes:
        niche_size_data.append({"sample":sample, "stage": "PTC" if sample in by_stage["PTC"] else "LPTC", "cluster_size": sz})
df_niche = pd.DataFrame(niche_size_data)
df_niche.to_csv(OUT/"spatial_F25_TROP2_niche_clusters.tsv", sep="\t", index=False)
# Plot size histograms per stage
for stage, color in [("PTC","#34547A"),("LPTC","#B8893C")]:
    sizes = df_niche[df_niche.stage == stage]["cluster_size"].values
    if len(sizes)==0: continue
    ax.hist(sizes, bins=np.logspace(0, 3.5, 30), alpha=0.7, label=f"{stage} (n_clusters={len(sizes)}, max={sizes.max()})",
            color=color, edgecolor="black", linewidth=0.4)
ax.set_xscale("log")
ax.set_xlabel("Niche cluster size (n_spots, log scale)")
ax.set_ylabel("Cluster count")
ax.set_title("S_F25  TROP2-high niche cluster size distribution (PTC + LPTC, top 15% TROP2 per sample)\n"
             "Size > 100 contiguous spots = large niche; many small clusters + few large = niche organization",
             fontsize=12, fontweight="bold")
ax.legend(loc="upper right", fontsize=10)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F25_TROP2_niche_cluster_sizes.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F26: Mega summary 6-panel =====
print("[F26] Mega summary 6-panel (one-glance overview)")
fig = plt.figure(figsize=(16, 10))
g = gs.GridSpec(2, 3, hspace=0.35, wspace=0.3, top=0.94, bottom=0.06, left=0.06, right=0.96)

# panel 1: TROP2 niche % per condition
ax = fig.add_subplot(g[0, 0])
cond_summary = df_E.groupby("condition").agg(
    mean_morans=("morans_I_TROP2","mean"),
    n_samp=("sample_id","count")
).reset_index().sort_values("mean_morans", ascending=False)
colors = [PAL.get(c,"#999") for c in cond_summary["condition"]]
bars = ax.barh(np.arange(len(cond_summary)), cond_summary["mean_morans"].values,
               color=colors, edgecolor="black", linewidth=0.5)
ax.set_yticks(np.arange(len(cond_summary)))
ax.set_yticklabels([f"{r['condition']} (n={int(r['n_samp'])})" for _,r in cond_summary.iterrows()], fontsize=10)
ax.invert_yaxis()
ax.axvline(0.25, ls="--", color="#7B1F2A", lw=1.2)
ax.set_xlabel("Mean TROP2 Moran's I"); ax.set_title("Panel 1. TROP2 niche per condition", fontsize=11, fontweight="bold")
ax.grid(axis="x", alpha=0.3)

# panel 2: HT TLS Moran summary
ax = fig.add_subplot(g[0, 1])
ht_means = df_C.groupby("gene_set")["morans_I"].mean().sort_values(ascending=False)
colors_ht = ["#1abc9c","#7B1F2A","#34547A","#B8893C"]
ax.bar(np.arange(len(ht_means)), ht_means.values, color=colors_ht[:len(ht_means)], edgecolor="black")
ax.set_xticks(np.arange(len(ht_means))); ax.set_xticklabels(ht_means.index, fontsize=10, family="monospace")
ax.set_ylabel("Mean Moran's I (4 PTC+HT slides)")
ax.axhline(0.3, ls="--", color="#3C6B4F", lw=1)
ax.set_title("Panel 2. HT immune niche organization", fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)

# panel 3: TROP2 vs DM1 spot ρ across stages (from S_F10 data)
ax = fig.add_subplot(g[0, 2])
trop_dm_rhos = []
for stage in ORDER:
    sample = rep.get(stage)
    if sample is None: continue
    sub_t = trop_all[trop_all.sample_id == sample]
    sub_t["row_col_id"] = sub_t["sample_id"]+"_"+sub_t["array_row"].astype(str)+"_"+sub_t["array_col"].astype(str)
    mg2 = sub_t.merge(g521[["row_col_id","DM1_like_score","Epithelial_score"]], on="row_col_id", how="left")
    rho_dm = spearmanr(mg2["TROP2"], mg2["DM1_like_score"], nan_policy="omit")[0]
    rho_epi = spearmanr(mg2["TROP2"], mg2["Epithelial_score"], nan_policy="omit")[0]
    trop_dm_rhos.append({"stage":stage, "TROP2_DM1":rho_dm, "TROP2_Epi":rho_epi})
df_rho = pd.DataFrame(trop_dm_rhos)
xs = np.arange(len(df_rho))
ax.bar(xs-0.18, df_rho["TROP2_DM1"], 0.36, color="#7B1F2A", edgecolor="black", label="TROP2 ↔ DM1")
ax.bar(xs+0.18, df_rho["TROP2_Epi"], 0.36, color="#B8893C", edgecolor="black", label="TROP2 ↔ Epithelial")
ax.set_xticks(xs); ax.set_xticklabels(df_rho["stage"], fontsize=10)
ax.axhline(0, color="black", lw=0.6)
ax.set_ylabel("Spot-level Spearman ρ"); ax.set_title("Panel 3. TROP2 axis specificity (rep niche per stage)", fontsize=11, fontweight="bold")
ax.legend(fontsize=9); ax.grid(axis="y", alpha=0.3)

# panel 4: Per-stage DM1 IQR (heterogeneity proxy)
ax = fig.add_subplot(g[1, 0])
df_B = pd.read_csv(OUT/"spatial_B_driver_orthogonal_per_sample.tsv", sep="\t")
positions = []; labels = []; data_iqr = []
for i, stage in enumerate(ORDER):
    sub = df_B[df_B.stage == stage]
    if not sub.empty:
        positions.append(i); labels.append(f"{stage}\nn={len(sub)}")
        data_iqr.append(sub["DM1_iqr"].values)
bp2 = ax.boxplot(data_iqr, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
for patch, st in zip(bp2["boxes"], ORDER):
    patch.set_facecolor(PAL[st]); patch.set_alpha(0.7)
for i, vals in enumerate(data_iqr):
    jit = np.random.RandomState(i).uniform(-0.08, 0.08, len(vals))
    ax.scatter(np.full_like(vals, positions[i])+jit, vals, s=40, color="black", zorder=3, edgecolor="white")
ax.set_xticks(positions); ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("DM1 IQR (within-sample spread)"); ax.set_title("Panel 4. DM1 spatial heterogeneity per stage", fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)

# panel 5: HT TLS Moran I per slide × 4 markers
ax = fig.add_subplot(g[1, 1])
samps = sorted(df_C["sample_id"].unique())
markers = ["HLA_II","B_cell","TLS","IGHV_AICDA"]
mat = np.full((len(samps), len(markers)), np.nan)
for i, s in enumerate(samps):
    for j, m in enumerate(markers):
        sub = df_C[(df_C.sample_id==s) & (df_C.gene_set==m)]
        if not sub.empty: mat[i,j] = sub["morans_I"].iloc[0]
im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.9)
for i in range(len(samps)):
    for j in range(len(markers)):
        v = mat[i,j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                color="white" if v>0.55 else "black", fontsize=9, fontweight="bold")
ax.set_xticks(range(len(markers))); ax.set_xticklabels(markers, fontsize=9.5, family="monospace")
ax.set_yticks(range(len(samps))); ax.set_yticklabels([s.split("_")[-1] for s in samps], fontsize=9.5, family="monospace")
ax.set_title("Panel 5. HT immune niche heatmap (4 slides × 4 sets)", fontsize=11, fontweight="bold")

# panel 6: closure battery ρ summary
ax = fig.add_subplot(g[1, 2])
clos = pd.read_csv(ROOT/"project/results/03_pathology_poc/closure_battery_metrics.tsv", sep="\t")
sec_e = clos[clos["section"].astype(str)=="E"]
labels_e = sec_e["experiment"].astype(str).str.replace("E_","").str.replace("_score_resid","").str.replace("_"," ")
vals_e = sec_e["pooled_spearman_r"].fillna(0).values
colors_e = ["#7B1F2A" if abs(v) < 0.3 else "#3C6B4F" for v in vals_e]
ax.barh(np.arange(len(sec_e)), vals_e, color=colors_e, edgecolor="black")
ax.set_yticks(np.arange(len(sec_e))); ax.set_yticklabels(labels_e, fontsize=8.5)
ax.axvline(0, color="black", lw=0.6)
ax.axvline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5); ax.axvline(-0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
ax.set_xlim(-0.5, 0.5)
ax.set_xlabel("LOSO ρ"); ax.set_title("Panel 6. H&E → DM1 closure (NO-GO)", fontsize=11, fontweight="bold")

fig.suptitle("S_F26  Spatial full-package mega summary — 6 key results in one frame",
             fontsize=14, fontweight="bold", y=0.99)
fig.savefig(ASSETS/"S_F26_mega_summary.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 8 new figures (S_F19-S_F26) saved")
