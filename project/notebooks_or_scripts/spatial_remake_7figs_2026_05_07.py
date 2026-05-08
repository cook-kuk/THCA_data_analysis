#!/usr/bin/env python3
"""Spatial REMAKE — 7 figures with corrected Visium hex 6-neighbor adjacency.

Targets (from SPATIAL_HEX_ADJACENCY_AUDIT_2026_05_07.md):
  S_F9   TROP2 ranked Moran's I per sample
  S_F11  Cross-cohort Moran I heatmap (MAIN M2)
  S_F18  HT TLS Moran I heatmap (MAIN M4)
  S_F22  Per-condition Moran I box
  S_F31  TROP2 niche cluster spacing (Bug A — scipy.ndimage.label default 4-conn)
  S_F40  Moran axis correlation
  S_F47  Per-axis Moran summary

Output:
  - project/results/spatial_full_2026_05_06/spatial_morans_v7plus_2026_05_07.tsv
  - project/papers_hub_2026_05_04/assets/spatial_full/{S_F9,11,18,22,31,40,47}_*.png (REMAKE)
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.colors import LinearSegmentedColormap
import warnings; warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT/"project/results/spatial_full_2026_05_06"
ASSETS = ROOT/"project/papers_hub_2026_05_04/assets/spatial_full"
G521H = ROOT/"project/data/processed/GSE250521"
EXTH = ROOT/"project_external_st/data/processed"
EXT_SC = ROOT/"project_external_st/results/scores"

# CORRECT Visium hex 6-neighbor offsets (post-fix)
HEX_OFFSETS = [(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]

PAL = {"PT":"#3C6B4F","PTC":"#34547A","LPTC":"#B8893C","ATC":"#7B1F2A",
       "PTC_HT":"#962E2E","CONTROL":"#52525a","HT":"#34547A","GD":"#B8893C"}
COND_ORDER = ["PT","PTC","LPTC","ATC","PTC_HT","CONTROL","HT","GD"]
DATASET_ORDER = ["GSE250521","GSE230424","GSE248205"]

def morans_hex(values, rows, cols):
    """Moran's I with Visium odd-r hex 6-neighbor adjacency."""
    n = len(values); v = np.asarray(values, dtype=float)
    if n < 30 or np.nanstd(v) == 0: return np.nan
    vc = v - np.nanmean(v)
    rows = np.asarray(rows); cols = np.asarray(cols)
    pos = {(int(r),int(c)): i for i,(r,c) in enumerate(zip(rows, cols))}
    Wnum, Wsum = 0.0, 0
    for i,(r,c) in enumerate(zip(rows, cols)):
        for dr, dc in HEX_OFFSETS:
            j = pos.get((int(r)+dr, int(c)+dc))
            if j is not None and not np.isnan(v[i]) and not np.isnan(v[j]):
                Wnum += vc[i]*vc[j]; Wsum += 1
    Wd = float(np.nansum(vc**2))
    if Wsum == 0 or Wd == 0: return np.nan
    return (n/Wsum)*(Wnum/Wd)

def visium_hex_clusters(rows, cols, mask):
    """Connected components on Visium hex grid using union-find with hex 6-nbr."""
    rows = np.asarray(rows); cols = np.asarray(cols); mask = np.asarray(mask, dtype=bool)
    pos = {(int(r),int(c)): i for i,(r,c) in enumerate(zip(rows, cols))}
    parent = list(range(len(rows)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    idx_mask = np.where(mask)[0]
    set_mask = set(idx_mask.tolist())
    for i in idx_mask:
        r, c = int(rows[i]), int(cols[i])
        for dr, dc in HEX_OFFSETS:
            j = pos.get((r+dr, c+dc))
            if j is not None and j in set_mask:
                union(i, j)
    labels = -np.ones(len(rows), dtype=int)
    for i in idx_mask:
        labels[i] = find(i)
    # remap labels to 1..K
    uniq = sorted(set(labels[idx_mask].tolist()))
    remap = {u: k+1 for k, u in enumerate(uniq)}
    labels_out = labels.copy()
    for i in idx_mask: labels_out[i] = remap[labels[i]]
    return labels_out  # -1 if not in mask, else 1..K

def score_h5ad(p, gene_set):
    a = ad.read_h5ad(p)
    rvar = a.raw.var_names.astype(str) if a.raw is not None else a.var_names.astype(str)
    X = a.raw.X if a.raw is not None else a.X
    avail = [g for g in gene_set if g in rvar.values]
    if not avail: return None, None, None
    idx = [list(rvar.values).index(g) for g in avail]
    sub = X[:, idx]
    if hasattr(sub,"toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", np.ones(a.n_obs))).flatten()
    tot = np.where(tot==0,1,tot)
    norm_log = np.log1p(sub/tot[:,None]*1e4).mean(axis=1)
    coords = (a.obs["array_row"].values, a.obs["array_col"].values) if "array_row" in a.obs.columns else (None, None)
    return norm_log, coords[0], coords[1]

print("[load] master per-spot tables")
trop = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t")
ext = pd.read_csv(ROOT/"project_external_st/results/scores/all_external_spots_scored.tsv.gz", sep="\t")

# Sample → (dataset, condition) mapping from trop file
sample_meta = trop[["sample_id","dataset","condition"]].drop_duplicates().set_index("sample_id")
print(f"  28 sample meta rows: {len(sample_meta)}")

# ---------- Compute Moran for each (sample, axis) ----------
print("[compute] Moran I (hex 6-nbr) for all 28 samples × axes")
records = []

# 1) TROP2 across all 28 samples
for s, row in sample_meta.iterrows():
    sub = trop[trop.sample_id == s]
    I = morans_hex(sub["TROP2"].values, sub["array_row"].values, sub["array_col"].values)
    records.append({"sample": s, "dataset": row["dataset"], "condition": row["condition"],
                    "axis": "TROP2", "morans_I": I, "n_spots": len(sub)})

# 2) 8 axes for GSE250521 (16 samples)
G521_AXES = ["DM1_like_score","RAI_8_score","TDS_like_score","Epithelial_score",
             "Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]
for s in g521["sample_id"].unique():
    sub = g521[g521.sample_id == s]
    cond_full = sample_meta.loc[s, "condition"]
    for axis in G521_AXES:
        I = morans_hex(sub[axis].values, sub["array_row"].values, sub["array_col"].values)
        records.append({"sample": s, "dataset": "GSE250521", "condition": cond_full,
                        "axis": axis, "morans_I": I, "n_spots": len(sub)})

# 3) Same 8 axes for external (GSE230424 + GSE248205) — using *_score_raw columns
EXT_AXES = [("RAI_8_score","RAI_8_score_raw"),
            ("TDS_like_score","TDS_overlap_score_raw"),
            ("Epithelial_score","Epithelial_score_raw"),
            ("Proliferation_score","Proliferation_score_raw"),
            ("CAF_ECM_score","CAF_ECM_score_raw"),
            ("EMT_score","EMT_score_raw"),
            ("Hypoxia_score","Hypoxia_score_raw"),
            ("DM1_like_score","DM1_like_score_raw")]
for s in ext["sample_id"].unique():
    sub = ext[ext.sample_id == s]
    if s not in sample_meta.index: continue
    cond_full = sample_meta.loc[s, "condition"]
    ds = sample_meta.loc[s, "dataset"]
    for axis_canon, col_raw in EXT_AXES:
        if col_raw not in sub.columns: continue
        I = morans_hex(sub[col_raw].values, sub["array_row"].values, sub["array_col"].values)
        records.append({"sample": s, "dataset": ds, "condition": cond_full,
                        "axis": axis_canon, "morans_I": I, "n_spots": len(sub)})

# 4) HT immune axes (HLA_II, B_cell, TLS, IGHV_AICDA) — score from h5ads
HLA_II_GENES = ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
B_CELL_GENES = ["CD19","MS4A1","CD79A","CD79B","SDC1","JCHAIN"]
TLS_GENES   = ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"]
IGHV_GENES  = ["AICDA","IGHM","IGHG1","IGHA1","IGKC","IGLC2"]
HT_SETS = [("HLA_II", HLA_II_GENES), ("B_cell", B_CELL_GENES),
           ("TLS", TLS_GENES), ("IGHV_AICDA", IGHV_GENES)]
ALL_H5AD = {}
for ds_dir in [G521H, EXTH/"GSE230424", EXTH/"GSE248205"]:
    for p in sorted(ds_dir.glob("*/GSM*.raw.h5ad")):
        ALL_H5AD[p.parent.name] = p
print(f"[h5ad] discovered {len(ALL_H5AD)} h5ads")
for s in sample_meta.index:
    if s not in ALL_H5AD: continue
    cond_full = sample_meta.loc[s, "condition"]
    ds = sample_meta.loc[s, "dataset"]
    for set_name, genes in HT_SETS:
        v, rows, cols = score_h5ad(ALL_H5AD[s], genes)
        if v is None or rows is None: continue
        I = morans_hex(v, rows, cols)
        records.append({"sample": s, "dataset": ds, "condition": cond_full,
                        "axis": set_name, "morans_I": I, "n_spots": len(v)})

df = pd.DataFrame(records)
df.to_csv(OUT/"spatial_morans_v7plus_2026_05_07.tsv", sep="\t", index=False)
print(f"[write] {OUT/'spatial_morans_v7plus_2026_05_07.tsv'}  rows={len(df)}")
piv = df.pivot_table(index=["sample","dataset","condition"], columns="axis",
                     values="morans_I", aggfunc="first")
print(piv.round(2).to_string())

# =============== FIGURES ===============

cmap_div = LinearSegmentedColormap.from_list("div", ["#34547A","#fffaf2","#7B1F2A"], N=256)
cmap_seq = LinearSegmentedColormap.from_list("seq", ["#fffaf2","#ffd58a","#7B1F2A","#0F1A2E"], N=256)

# ----- S_F9 — TROP2 ranked Moran's I per sample (28 bars) -----
print("[F9] TROP2 ranked Moran's I per sample")
trop_morans = df[df.axis == "TROP2"].copy().sort_values("morans_I", ascending=False)
fig = plt.figure(figsize=(15, 5.6))
ax = fig.add_subplot(111)
xs = np.arange(len(trop_morans))
colors = [PAL.get(c, "#888") for c in trop_morans["condition"].values]
ax.bar(xs, trop_morans["morans_I"].values, color=colors, edgecolor="black", linewidth=0.4)
labels = [s.split("_",1)[1] for s in trop_morans["sample"].values]
ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=8.5)
ax.axhline(0.13, color="#888", lw=0.7, ls="--", label="niche threshold (Moran's I = 0.13)")
ax.set_ylabel("Moran's I (TROP2, hex 6-nbr)", fontsize=11)
ax.set_title("S_F9  TROP2 ranked Moran's I per sample (28 slides × 3 cohorts) — REMADE 2026-05-07 with Visium hex 6-neighbor adjacency",
             fontsize=11.5, fontweight="bold")
import matplotlib.patches as mp
legend_handles = [mp.Patch(facecolor=PAL[c], label=c) for c in COND_ORDER if c in trop_morans["condition"].values]
ax.legend(handles=legend_handles + [plt.Line2D([0],[0], color="#888", ls="--", label="threshold 0.13")],
          loc="upper right", fontsize=9, frameon=False)
ax.grid(True, ls=":", alpha=0.4, axis="y")
fig.tight_layout()
fig.savefig(ASSETS/"S_F9_TROP2_ranked_morans.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ----- S_F11 — Cross-cohort Moran I heatmap (MAIN M2) -----
print("[F11] Cross-cohort Moran I heatmap (MAIN M2)")
PIVOT_AXES = ["TROP2","DM1_like_score","RAI_8_score","Epithelial_score","Proliferation_score",
              "CAF_ECM_score","EMT_score","Hypoxia_score","HLA_II","B_cell","TLS","IGHV_AICDA"]
piv_full = df.pivot_table(index="sample", columns="axis", values="morans_I", aggfunc="first")
# order samples by dataset then condition
sample_meta_sorted = sample_meta.copy()
sample_meta_sorted["cond_rank"] = sample_meta_sorted["condition"].map({c:i for i,c in enumerate(COND_ORDER)})
sample_meta_sorted["ds_rank"] = sample_meta_sorted["dataset"].map({d:i for i,d in enumerate(DATASET_ORDER)})
sample_order = sample_meta_sorted.sort_values(["ds_rank","cond_rank"]).index.tolist()
piv_full = piv_full.reindex(sample_order)
piv_full = piv_full[[a for a in PIVOT_AXES if a in piv_full.columns]]

fig = plt.figure(figsize=(13, 11))
ax = fig.add_subplot(111)
M = piv_full.values
im = ax.imshow(M, cmap=cmap_div, vmin=-0.5, vmax=0.85, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=7.5, color="white" if abs(v) > 0.5 else "black")
ax.set_yticks(range(len(piv_full.index)))
y_labels = [f"{s.split('_',1)[1]}  [{sample_meta.loc[s,'condition']}]" for s in piv_full.index]
ax.set_yticklabels(y_labels, fontsize=8.5)
ax.set_xticks(range(len(piv_full.columns)))
ax.set_xticklabels(piv_full.columns, rotation=45, ha="right", fontsize=10)
# color y-tick labels by condition
for tick, s in zip(ax.get_yticklabels(), piv_full.index):
    tick.set_color(PAL.get(sample_meta.loc[s, "condition"], "#0F1A2E"))
ax.set_title("S_F11  Cross-cohort Moran's I heatmap (MAIN M2) — 28 samples × 12 axes  ·  REMADE 2026-05-07 with Visium hex 6-neighbor adjacency\n"
             "TROP2 niche tumor-specific (PTC+LPTC 8/8 high; PT/autoimmune ~0); HT-overlap PTC TLS niche directly visible",
             fontsize=11.5, fontweight="bold")
plt.colorbar(im, ax=ax, fraction=0.018, pad=0.02, label="Moran's I (hex 6-nbr)")
fig.tight_layout()
fig.savefig(ASSETS/"S_F11_cross_cohort_heatmap.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ----- S_F18 — HT TLS Moran I heatmap (MAIN M4) -----
print("[F18] HT TLS Moran I heatmap (MAIN M4)")
ht_axes = ["HLA_II","B_cell","TLS","IGHV_AICDA"]
ht_samples = sample_meta[sample_meta["dataset"] == "GSE230424"].index.tolist()
ht_pivot = df[(df["sample"].isin(ht_samples)) & (df["axis"].isin(ht_axes))].pivot_table(
    index="sample", columns="axis", values="morans_I", aggfunc="first")
ht_pivot = ht_pivot.reindex(ht_samples)[ht_axes]

fig = plt.figure(figsize=(8, 4.6))
ax = fig.add_subplot(111)
M = ht_pivot.values
im = ax.imshow(M, cmap=cmap_seq, vmin=0, vmax=0.85, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=11, fontweight="bold",
                    color="white" if v > 0.5 else "black")
ax.set_yticks(range(len(ht_samples)))
ax.set_yticklabels([s.split("_",1)[1] for s in ht_samples], fontsize=11, color="#962E2E", fontweight="bold")
ax.set_xticks(range(len(ht_axes)))
ax.set_xticklabels(ht_axes, rotation=0, fontsize=11)
ax.set_title("S_F18  HT-overlap PTC immune Moran's I heatmap (MAIN M4)\n"
             "GSE230424 4/4 PTC+HT slides — REMADE 2026-05-07 with Visium hex 6-nbr",
             fontsize=11.5, fontweight="bold")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Moran's I (hex 6-nbr)")
fig.tight_layout()
fig.savefig(ASSETS/"S_F18_HT_TLS_heatmap.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ----- S_F22 — Per-condition Moran I (TROP2) box -----
print("[F22] Per-condition Moran I (TROP2) box")
fig = plt.figure(figsize=(11, 5.2))
ax = fig.add_subplot(111)
trop_per_cond = []
labels_cond = []
positions = []
for i, c in enumerate(COND_ORDER):
    vals = df[(df.axis=="TROP2") & (df.condition==c)]["morans_I"].dropna().values
    if len(vals) == 0: continue
    trop_per_cond.append(vals); labels_cond.append(c); positions.append(i)
bp = ax.boxplot(trop_per_cond, positions=positions, widths=0.55,
                patch_artist=True, showfliers=False,
                medianprops=dict(color="#0F1A2E", lw=1.6))
for box, c in zip(bp["boxes"], labels_cond):
    box.set_facecolor(PAL[c]); box.set_edgecolor("#0F1A2E"); box.set_alpha(0.85)
# scatter overlays
for i, (vals, c) in enumerate(zip(trop_per_cond, labels_cond)):
    jitter = (np.random.RandomState(i).rand(len(vals)) - 0.5) * 0.3
    ax.scatter(np.full(len(vals), positions[i]) + jitter, vals,
               color=PAL[c], edgecolor="black", lw=0.5, s=58, zorder=3)
ax.set_xticks(positions); ax.set_xticklabels(labels_cond, fontsize=11, fontweight="bold")
for tick, c in zip(ax.get_xticklabels(), labels_cond):
    tick.set_color(PAL[c])
ax.set_ylabel("Moran's I (TROP2, hex 6-nbr)", fontsize=11)
ax.axhline(0.13, color="#888", lw=0.7, ls="--", label="niche threshold")
ax.set_title("S_F22  TROP2 Moran's I distribution per condition (28 samples × 3 cohorts) — REMADE 2026-05-07",
             fontsize=11.5, fontweight="bold")
ax.grid(True, ls=":", alpha=0.4)
ax.legend(loc="upper right", fontsize=9, frameon=False)
fig.tight_layout()
fig.savefig(ASSETS/"S_F22_morans_per_condition_box.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ----- S_F31 — TROP2 niche cluster spacing (Bug A fix) -----
print("[F31] TROP2 niche cluster spacing (hex 6-nbr fix)")
ptc_lptc = sample_meta[sample_meta["condition"].isin(["PTC","LPTC"])].index.tolist()
spacings = []
size_recs = []
for s in ptc_lptc:
    sub = trop[trop.sample_id == s].sort_values("spot_id")
    if len(sub) < 100: continue
    rows_arr = sub["array_row"].values
    cols_arr = sub["array_col"].values
    pxr = sub["pxl_row"].values
    pxc = sub["pxl_col"].values
    v = sub["TROP2"].values
    if np.nanstd(v) == 0: continue
    thr = np.nanpercentile(v, 85)  # top 15%
    mask = v >= thr
    labels = visium_hex_clusters(rows_arr, cols_arr, mask)
    K = labels.max()
    if K < 2: continue
    # cluster centroids in pixel space
    centroids = []
    for k in range(1, K+1):
        m_k = labels == k
        sz = int(m_k.sum())
        size_recs.append({"sample": s, "cluster_id": k, "size": sz})
        if sz < 3: continue
        centroids.append((np.mean(pxr[m_k]), np.mean(pxc[m_k])))
    if len(centroids) < 2: continue
    centroids = np.array(centroids)
    # pairwise distances (Visium spot units approx via px / median spot pitch)
    d = np.sqrt(((centroids[:,None,:] - centroids[None,:,:])**2).sum(-1))
    iu = np.triu_indices(len(centroids), k=1)
    spacings.extend(list(d[iu]))

sizes_df = pd.DataFrame(size_recs)
sizes_df.to_csv(OUT/"spatial_F31_niche_clusters_v7plus.tsv", sep="\t", index=False)

fig = plt.figure(figsize=(11, 5.2)); ax = fig.add_subplot(111)
if spacings:
    sp_arr = np.asarray(spacings)
    bins = np.logspace(np.log10(max(1, sp_arr.min())), np.log10(sp_arr.max()+1), 40)
    ax.hist(sp_arr, bins=bins, color="#7B1F2A", edgecolor="black", lw=0.4, alpha=0.85)
    ax.set_xscale("log")
    ax.set_xlabel("Pairwise distance between TROP2-high niche cluster centroids (pixel units)", fontsize=11)
    ax.set_ylabel("Cluster pair count", fontsize=11)
    ax.set_title(f"S_F31  TROP2 niche cluster spacing — PTC+LPTC pairwise centroid distances\n"
                 f"REMADE 2026-05-07 with Visium hex 6-neighbor union-find  ·  n={len(spacings)} cluster pairs across {len(set([r['sample'] for r in size_recs]))} slides",
                 fontsize=11.5, fontweight="bold")
else:
    ax.text(0.5, 0.5, "No clusters detected (post-fix)", ha="center", va="center", fontsize=14)
ax.grid(True, ls=":", alpha=0.4)
fig.tight_layout()
fig.savefig(ASSETS/"S_F31_niche_cluster_spacing.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ----- S_F40 — Moran axis correlation (across samples) -----
print("[F40] Moran axis correlation")
piv_for_corr = df.pivot_table(index="sample", columns="axis", values="morans_I", aggfunc="first")
piv_for_corr = piv_for_corr[[a for a in PIVOT_AXES if a in piv_for_corr.columns]]
corr = piv_for_corr.corr(method="spearman")
fig = plt.figure(figsize=(10, 9)); ax = fig.add_subplot(111)
M = corr.values
im = ax.imshow(M, cmap=cmap_div, vmin=-1, vmax=1)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=8, color="white" if abs(v) > 0.6 else "black")
ax.set_xticks(range(len(corr.columns))); ax.set_yticks(range(len(corr.index)))
ax.set_xticklabels(corr.columns, rotation=45, ha="right", fontsize=9.5)
ax.set_yticklabels(corr.index, fontsize=9.5)
ax.set_title("S_F40  Cross-axis Spearman correlation of Moran's I across 28 samples\n"
             "REMADE 2026-05-07 with Visium hex 6-neighbor",
             fontsize=11.5, fontweight="bold")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman ρ")
fig.tight_layout()
fig.savefig(ASSETS/"S_F40_morans_axis_correlation.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ----- S_F47 — Per-axis Moran summary (mean ± SD per condition) -----
print("[F47] Per-axis Moran summary (per condition)")
summary = df.groupby(["axis","condition"])["morans_I"].agg(["mean","std","count"]).reset_index()
summary.to_csv(OUT/"spatial_F47_per_axis_moran_summary_v7plus.tsv", sep="\t", index=False)

fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(111)
axes_order = [a for a in PIVOT_AXES if a in summary["axis"].values]
n_ax = len(axes_order); n_cond = len(COND_ORDER)
xs = np.arange(n_ax)
bar_w = 0.10
for i, c in enumerate(COND_ORDER):
    means, sds = [], []
    for a in axes_order:
        sub_s = summary[(summary.axis==a)&(summary.condition==c)]
        if sub_s.empty:
            means.append(np.nan); sds.append(0)
        else:
            means.append(sub_s["mean"].iloc[0])
            sds.append(sub_s["std"].iloc[0] if not np.isnan(sub_s["std"].iloc[0]) else 0)
    pos = xs + (i - n_cond/2) * bar_w + bar_w/2
    ax.bar(pos, means, width=bar_w, color=PAL[c], edgecolor="black", lw=0.4, label=c, alpha=0.88)
    ax.errorbar(pos, means, yerr=sds, fmt="none", ecolor="#0F1A2E", elinewidth=0.7, capsize=2)
ax.set_xticks(xs); ax.set_xticklabels(axes_order, rotation=35, ha="right", fontsize=10)
ax.set_ylabel("Mean Moran's I per condition (± SD)", fontsize=11)
ax.axhline(0, color="#888", lw=0.6)
ax.legend(loc="upper right", fontsize=9, frameon=False, ncol=2)
ax.set_title("S_F47  Per-axis Moran's I summary by condition — REMADE 2026-05-07 with Visium hex 6-neighbor (28 samples × 3 cohorts)",
             fontsize=11.5, fontweight="bold")
ax.grid(True, ls=":", alpha=0.4, axis="y")
fig.tight_layout()
fig.savefig(ASSETS/"S_F47_per_axis_moran_summary.png", dpi=160, bbox_inches="tight")
plt.close(fig)

print("\n=== DONE — 7 figures REMADE ===")
print(f"  master TSV:  {OUT/'spatial_morans_v7plus_2026_05_07.tsv'}")
print(f"  cluster TSV: {OUT/'spatial_F31_niche_clusters_v7plus.tsv'}")
print(f"  summary TSV: {OUT/'spatial_F47_per_axis_moran_summary_v7plus.tsv'}")
print(f"  PNGs:        S_F9, S_F11, S_F18, S_F22, S_F31, S_F40, S_F47")
print("\n[verify] TROP2 by condition (post-fix):")
chk = df[df.axis=="TROP2"].groupby("condition")["morans_I"].agg(["min","mean","max","count"])
print(chk.round(3))
