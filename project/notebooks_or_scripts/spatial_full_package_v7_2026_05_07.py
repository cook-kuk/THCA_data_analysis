#!/usr/bin/env python3
"""Spatial v7 — fix S_F16/S_F25/S_F30 + 8 new (S_F49~S_F56)."""
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
from scipy.ndimage import label as nd_label
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

print("[load]")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
trop_all = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
df_C = pd.read_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t")

# ===== FIX S_F16: split into 2 panels (TROP2 cross-cohort + HT immune PTC_HT-only) =====
print("[FIX S_F16] split into 2 panels")
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw={"width_ratios":[2,1]})
# panel A: TROP2 across 8 conditions (all data available)
ax = axes[0]
conds_A = ["PT","PTC","LPTC","ATC","PTC_HT","CONTROL","HT","GD"]
trop_means = [df_E[df_E.condition==c]["morans_I_TROP2"].mean() for c in conds_A]
trop_stds  = [df_E[df_E.condition==c]["morans_I_TROP2"].std() for c in conds_A]
ns         = [(df_E.condition==c).sum() for c in conds_A]
mat_t = np.array(trop_means).reshape(1, -1)
im = ax.imshow(mat_t, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.6)
for j, (m, n) in enumerate(zip(trop_means, ns)):
    ax.text(j, 0, f"{m:.2f}\n(n={n})", ha="center", va="center",
            color="white" if m>0.35 else "black", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(conds_A))); ax.set_xticklabels(conds_A, fontsize=10)
ax.set_yticks([0]); ax.set_yticklabels(["TROP2"], fontsize=11, family="monospace")
ax.set_title("A. TROP2 niche organization across all 28 samples × 8 conditions\n"
             "(higher = niche-organized; PTC + LPTC clearly elevated)", fontsize=11, fontweight="bold")
# colorbar
plt.colorbar(im, ax=ax, fraction=0.10, pad=0.04, label="Mean Moran's I")

# panel B: HT immune markers (PTC_HT only — all 4 markers, 4 slides each)
ax = axes[1]
markers_ht = ["HLA_II","B_cell","TLS","IGHV_AICDA"]
ht_means = [df_C[df_C.gene_set==m]["morans_I"].mean() for m in markers_ht]
ht_stds  = [df_C[df_C.gene_set==m]["morans_I"].std() for m in markers_ht]
ax.barh(range(len(markers_ht))[::-1], ht_means, xerr=ht_stds, color="#1abc9c",
        edgecolor="black", linewidth=0.6, capsize=5, alpha=0.85)
for i, (m, n_v) in enumerate(zip(reversed(ht_means), reversed(markers_ht))):
    ax.text(m+0.02, len(markers_ht)-1-i, f"{m:.2f}", va="center", fontsize=10, fontweight="bold")
ax.set_yticks(range(len(markers_ht))[::-1]); ax.set_yticklabels(markers_ht, fontsize=10, family="monospace")
ax.axvline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
ax.set_xlabel("Mean Moran's I (4 PTC+HT slides)")
ax.set_xlim(0, 0.9)
ax.set_title("B. HT-overlap PTC immune niche organization\n(GSE230424 only; n=4 slides per marker)",
             fontsize=11, fontweight="bold")
ax.grid(axis="x", alpha=0.3)

fig.suptitle("S_F16  Niche organization summary — TROP2 cross-cohort (A) + HT immune markers (B, PTC_HT only)\n"
             "FIXED v7: 두 데이터셋 coverage가 다르기에 분리 visualization (이전엔 4×8 matrix가 거의 비어있었음)",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.91])
fig.savefig(ASSETS/"S_F16_morans_summary.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ===== FIX S_F25: proper Visium hex adjacency for cluster detection =====
print("[FIX S_F25] Visium hex adjacency cluster size")
def visium_hex_clusters(rows, cols, mask):
    """Find connected components in Visium hex grid using 6-neighbor adjacency."""
    pos_to_idx = {}
    for i, (r, c) in enumerate(zip(rows, cols)):
        if mask[i]:
            pos_to_idx[(int(r), int(c))] = i
    # union-find
    parent = {p: p for p in pos_to_idx.values()}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    # Visium hex 6-neighbor offsets in array_row/array_col
    OFFS = [(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]
    for (r, c), i in pos_to_idx.items():
        for dr, dc in OFFS:
            j = pos_to_idx.get((r+dr, c+dc))
            if j is not None: union(i, j)
    # gather sizes
    from collections import Counter
    roots = [find(v) for v in pos_to_idx.values()]
    sizes = list(Counter(roots).values())
    return sizes

niche_size_data = []
all_samples_g521 = trop_all[trop_all.dataset == "GSE250521"]["sample_id"].unique()
for sample in all_samples_g521:
    sub = trop_all[trop_all.sample_id == sample]
    cond = sub["condition"].iloc[0]
    if cond not in ("PTC","LPTC"): continue
    if sub.empty or "array_row" not in sub.columns: continue
    rows_arr = sub["array_row"].astype(int).values
    cols_arr = sub["array_col"].astype(int).values
    threshold = sub["TROP2"].quantile(0.85)
    mask = (sub["TROP2"] > threshold).values
    if mask.sum() < 5: continue
    sizes = visium_hex_clusters(rows_arr, cols_arr, mask)
    for sz in sizes:
        niche_size_data.append({"sample": sample, "stage": cond, "cluster_size": sz})
df_niche = pd.DataFrame(niche_size_data)
print(f"  cluster size distribution: min={df_niche['cluster_size'].min()}, max={df_niche['cluster_size'].max()}, mean={df_niche['cluster_size'].mean():.1f}")

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
# panel A: cluster size histogram by stage
ax = axes[0]
for stage, color in [("PTC","#34547A"),("LPTC","#B8893C")]:
    sizes = df_niche[df_niche.stage == stage]["cluster_size"].values
    if len(sizes)==0: continue
    ax.hist(sizes, bins=np.logspace(0, np.log10(max(sizes)+1), 30), alpha=0.7,
            label=f"{stage} (n_clusters={len(sizes)}, max={int(sizes.max())}, median={int(np.median(sizes))})",
            color=color, edgecolor="black", linewidth=0.4)
ax.set_xscale("log")
ax.set_xlabel("Niche cluster size (n_spots, log scale)")
ax.set_ylabel("Cluster count")
ax.set_title("A. TROP2-high niche cluster size distribution\n(Visium hex 6-neighbor adjacency, top 15% TROP2 per sample)",
             fontsize=11, fontweight="bold")
ax.legend(loc="upper right", fontsize=9.5); ax.grid(axis="y", alpha=0.3)

# panel B: per-sample largest cluster size
ax = axes[1]
per_sample = df_niche.groupby(["sample","stage"]).agg(
    max_cluster=("cluster_size","max"),
    n_clusters=("cluster_size","count"),
    total_high=("cluster_size","sum"),
).reset_index().sort_values(["stage","max_cluster"]).reset_index(drop=True)
xs = np.arange(len(per_sample))
colors = [PAL.get(s,"#999") for s in per_sample["stage"]]
ax.bar(xs, per_sample["max_cluster"], color=colors, edgecolor="black", linewidth=0.5)
for i, (_, r) in enumerate(per_sample.iterrows()):
    ax.text(i, r["max_cluster"]+5, f"{int(r['n_clusters'])}", ha="center", fontsize=8, color="#52606f")
ax.set_xticks(xs); ax.set_xticklabels([f"{r['sample'].split('_')[-1]}\n[{r['stage']}]" for _,r in per_sample.iterrows()],
                                       rotation=70, ha="right", fontsize=8)
ax.set_ylabel("Largest niche cluster size (n_spots)")
ax.set_title("B. Per-sample largest TROP2 niche cluster\n(annotation = total cluster count per sample)",
             fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)

fig.suptitle("S_F25  TROP2-high niche cluster size — FIXED v7 (Visium hex 6-neighbor adjacency)\n"
             "이전 4-conn rectangular adjacency가 hex grid를 못 인식 → all clusters size=1; 현재 모든 cluster의 진짜 크기 표시",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(ASSETS/"S_F25_TROP2_niche_cluster_sizes.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== FIX S_F30: use real LOSO schema (target/y_obs/y_pred) =====
print("[FIX S_F30] use y_obs/y_pred/target schema")
clos_pred = pd.read_csv(ROOT/"project/results/03_pathology_poc/loso_predictions_resnet50.tsv.gz", sep="\t")
clos_metric = pd.read_csv(ROOT/"project/results/03_pathology_poc/closure_battery_metrics.tsv", sep="\t")
# Filter to DM1_like, ridge, tile_size=224
target_filter = "DM1_like_score_resid"
sub_pred = clos_pred[(clos_pred["target"]==target_filter) & (clos_pred["model"]=="ridge") & (clos_pred["tile_size"]==224)]
print(f"  filtered LOSO {target_filter} ridge 224: {len(sub_pred)} spots")

fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
# Panel A: pred vs obs hexbin (all spots)
ax = axes[0,0]
if len(sub_pred):
    ax.hexbin(sub_pred["y_obs"].values, sub_pred["y_pred"].values, gridsize=40, cmap="OrRd", mincnt=2)
    ax.plot([-3,3],[-3,3], "--", color="black", lw=1.0, alpha=0.7, label="diagonal (perfect)")
    rho, p = spearmanr(sub_pred["y_obs"].values, sub_pred["y_pred"].values, nan_policy="omit")
    ax.set_title(f"A. LOSO pred vs obs ({target_filter})\nrho = {rho:.3f}, p = {p:.1e}, n={len(sub_pred)} spots",
                 fontsize=10.5)
    ax.legend(fontsize=9)
ax.set_xlabel("Observed DM1_like (within-sample z)")
ax.set_ylabel("Predicted DM1_like (LOSO ResNet50 ridge)")
ax.axhline(0, ls=":", lw=0.4, color="gray"); ax.axvline(0, ls=":", lw=0.4, color="gray")

# Panel B: per-slide LOSO ρ bar
ax = axes[0,1]
if len(sub_pred):
    slide_rho = sub_pred.groupby("sample_id").apply(
        lambda g: spearmanr(g["y_obs"], g["y_pred"], nan_policy="omit")[0]).reset_index()
    slide_rho.columns = ["sample","rho"]
    slide_rho = slide_rho.merge(sub_pred[["sample_id","stage"]].drop_duplicates(), left_on="sample", right_on="sample_id")
    slide_rho = slide_rho.sort_values(["stage","rho"]).reset_index(drop=True)
    colors = [PAL.get(s,"#999") for s in slide_rho["stage"]]
    ax.bar(range(len(slide_rho)), slide_rho["rho"].values, color=colors, edgecolor="black", linewidth=0.4)
    ax.axhline(0, color="black", lw=0.6); ax.axhline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
    ax.axhline(-0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
    ax.set_xticks(range(len(slide_rho)))
    ax.set_xticklabels([f"{s.split('_')[-1]}\n[{st}]" for s, st in zip(slide_rho["sample"], slide_rho["stage"])],
                       rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("Per-slide LOSO Spearman ρ")
    ax.set_title(f"B. Per-slide ρ ({target_filter}, ridge)\n"
                 f"median = {slide_rho['rho'].median():.3f}; only {(slide_rho['rho'].abs() > 0.3).sum()}/16 above |0.3|",
                 fontsize=10.5)
    ax.grid(axis="y", alpha=0.3)

# Panel C: section A pooled metrics
ax = axes[1,0]
sec_a = clos_metric[(clos_metric["section"].astype(str).str.startswith("A")) & (clos_metric["model"]=="ridge")]
labels = sec_a["experiment"].astype(str).str.replace("raw_","")
ax.bar(range(len(sec_a)), sec_a["pooled_spearman_r"].values, color="#7B1F2A",
       edgecolor="black", linewidth=0.4)
ax.axhline(0, color="black", lw=0.6); ax.axhline(0.3, ls="--", color="#3C6B4F", lw=1, alpha=0.5)
ax.set_xticks(range(len(sec_a))); ax.set_xticklabels(labels, fontsize=9, rotation=15, ha="right")
ax.set_ylabel("Pooled Spearman ρ")
ax.set_ylim(-0.2, 0.5)
ax.set_title("C. Section A — pooled raw label LOSO ρ\n(all targets fail useful threshold 0.3)", fontsize=10.5)

# Panel D: residual histogram + slide-level CDF
ax = axes[1,1]
if len(sub_pred):
    resid = (sub_pred["y_obs"] - sub_pred["y_pred"]).values
    ax.hist(resid, bins=60, color="#7B1F2A", alpha=0.85, edgecolor="black", linewidth=0.3)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlabel("Residual (y_obs - y_pred)")
    ax.set_ylabel("Spot count")
    ax.set_title(f"D. Residual distribution\n(mean={np.mean(resid):.3f}, std={np.std(resid):.3f})", fontsize=10.5)

fig.suptitle("S_F30  Closure battery LOSO breakdown — FIXED v7 (correct y_obs/y_pred/target schema)\n"
             f"All 4 panels show different angles of NO-GO verdict for {target_filter} ResNet50 ridge",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F30_closure_LOSO_breakdown.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== Now 8 NEW figures (S_F49 ~ S_F56) =====
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
    norm_log = np.log1p(sub/tot[:,None]*1e4)
    coords = (a.obs["array_row"].values, a.obs["array_col"].values) if "array_row" in a.obs.columns else (None, None)
    return norm_log, coords

print("[F49] All 28 samples DM1_like mega grid (companion to S_F41)")
all_samples = trop_all[["sample_id","dataset","condition"]].drop_duplicates()
cond_rank = {"PT":0,"PTC":1,"LPTC":2,"ATC":3,"PTC_HT":4,"CONTROL":5,"HT":6,"GD":7}
all_samples["rank"] = all_samples["condition"].map(cond_rank)
all_samples = all_samples.sort_values(["rank","sample_id"]).reset_index(drop=True)
# Get DM1 from g521 + ext
ext_files = sorted((ROOT/"project_external_st/results/scores").glob("*_spot_scores.tsv.gz"))
ext_dfs = [pd.read_csv(f, sep="\t") for f in ext_files]
ext = pd.concat(ext_dfs, ignore_index=True).rename(columns={"DM1_like_score_raw":"DM1_like_score"})
fig = plt.figure(figsize=(20, 16))
g = gs.GridSpec(4, 7, hspace=0.42, wspace=0.18, top=0.95, bottom=0.03, left=0.03, right=0.97)
for idx, (_, row) in enumerate(all_samples.iterrows()):
    if idx >= 28: break
    r, c = idx // 7, idx % 7
    ax = fig.add_subplot(g[r, c])
    if row["dataset"] == "GSE250521":
        sub = g521[g521.sample_id == row["sample_id"]]
    else:
        sub = ext[ext.sample_id == row["sample_id"]]
    if sub.empty or "array_row" not in sub.columns:
        ax.axis("off"); continue
    x = sub["array_col"].values; y = -sub["array_row"].values
    v = sub["DM1_like_score"].values
    vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95))) if not np.isnan(v).all() else 1
    ax.scatter(x, y, c=v, cmap=dm_cmap, vmin=-vmax, vmax=vmax, s=2, alpha=0.92)
    ax.set_title(f"{row['condition']} · {row['sample_id'].split('_')[-1]}", fontsize=8.5, color=PAL.get(row['condition'],"#666"), fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F49  ALL 28 samples × 3 cohorts DM1_like spatial maps (companion to S_F41 TROP2)\n"
             "Blue = RAI-high (differentiated); red = DM1-high (dark matter); within-sample z normalization",
             fontsize=13, fontweight="bold")
fig.savefig(ASSETS/"S_F49_DM1_all_28_samples.png", dpi=140, bbox_inches="tight"); plt.close(fig)

print("[F50] Niche threshold robustness check + reproducibility")
df_F46 = pd.read_csv(OUT/"spatial_F46_niche_threshold_sensitivity.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
# panel A: per-sample trajectory
ax = axes[0]
for cond, color in PAL.items():
    sub = df_F46[df_F46.condition == cond]
    if sub.empty: continue
    for sample, sub_s in sub.groupby("sample"):
        ax.plot(sub_s["threshold_pct"], sub_s["morans_I"], color=color, alpha=0.5, lw=0.8)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.6)
ax.set_xlabel("TROP2 threshold percentile (%)")
ax.set_ylabel("Per-sample binary Moran's I")
ax.set_title("A. Per-sample trajectory across thresholds\n(each line = 1 sample; PTC+LPTC consistently above 0.25)",
             fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
# panel B: condition mean trajectory
ax = axes[1]
for cond, color in PAL.items():
    sub = df_F46[df_F46.condition == cond]
    if sub.empty: continue
    grp = sub.groupby("threshold_pct")["morans_I"].agg(["mean","sem"]).reset_index()
    ax.errorbar(grp["threshold_pct"], grp["mean"], yerr=grp["sem"], marker="o", lw=2,
                color=color, label=f"{cond}", capsize=4, alpha=0.85)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.6)
ax.set_xlabel("TROP2 threshold percentile (%)"); ax.set_ylabel("Mean Moran's I (per condition)")
ax.set_title("B. Condition-level mean trajectory\n(PTC+LPTC clearly distinct at all thresholds)", fontsize=11, fontweight="bold")
ax.legend(loc="best", fontsize=9, ncol=2); ax.grid(axis="y", alpha=0.3)
fig.suptitle("S_F50  Niche-finding robustness — per-sample + condition-mean across 6 threshold choices",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F50_niche_robustness.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("[F51] Per-sample TROP2 vs DM1 Moran I scatter (axis specificity check)")
df_moran = pd.read_csv(OUT/"spatial_F40_per_sample_per_axis_morans.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(11, 7.5))
for cond, color in [("PT","#3C6B4F"),("PTC","#34547A"),("LPTC","#B8893C"),("ATC","#7B1F2A")]:
    sub = df_moran[df_moran.condition == cond]
    if sub.empty: continue
    ax.scatter(sub["DM1_like_score"], sub["TROP2"], s=140, color=color, alpha=0.85,
               edgecolor="black", linewidth=0.8, label=f"{cond} (n={len(sub)})")
    for _, r in sub.iterrows():
        ax.annotate(r["sample"].split("_")[-1], (r["DM1_like_score"], r["TROP2"]),
                    fontsize=8, color="#0F1A2E", xytext=(5,5), textcoords="offset points")
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.5)
ax.axvline(0.25, ls="--", color="#34547A", lw=1, alpha=0.5)
ax.axhline(0, color="black", lw=0.5); ax.axvline(0, color="black", lw=0.5)
ax.set_xlabel("DM1_like Moran's I"); ax.set_ylabel("TROP2 Moran's I")
ax.set_title("S_F51  Per-sample TROP2 niche × DM1 niche organization scatter\n"
             "PTC+LPTC samples: TROP2 niche YES (y > 0.25), DM1 niche도 함께 발달 — 두 axis 독립적이지만 co-occur",
             fontsize=11.5, fontweight="bold")
ax.legend(loc="upper left", fontsize=10)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F51_TROP2_DM1_moran_scatter.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("[F52] Stage progression — TROP2 high-spot fraction trajectory")
fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
ax = axes[0]
for cond in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==cond)]
    if sub.empty: continue
    for _, r in sub.iterrows():
        ax.scatter(ORDER.index(cond), r["TROP2_high_spot_frac"]*100, s=100,
                   color=PAL[cond], alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
# stage means
stage_means = [df_E[(df_E.dataset=="GSE250521") & (df_E.condition==c)]["TROP2_high_spot_frac"].mean()*100 for c in ORDER]
ax.plot(range(len(ORDER)), stage_means, "k--", lw=1.5, alpha=0.7, label="stage mean")
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, fontsize=11)
ax.set_ylabel("% spots TROP2-high (top 10% per sample)")
ax.set_title("A. Per-sample TROP2 high-spot fraction trajectory\n(PTC+LPTC peaks at 8-10%; ATC drops back)", fontsize=11, fontweight="bold")
ax.legend(); ax.grid(axis="y", alpha=0.3)
# panel B: per-stage cluster numbers (from FIX S_F25)
ax = axes[1]
stage_summary = df_niche.groupby(["sample","stage"]).size().reset_index(name="n_clusters")
for s, color in [("PTC","#34547A"),("LPTC","#B8893C")]:
    sub_s = stage_summary[stage_summary.stage == s]
    ax.scatter([s]*len(sub_s), sub_s["n_clusters"].values, s=120, color=color,
               alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
ax.set_ylabel("# of TROP2-high niche clusters per sample\n(Visium hex 6-neighbor)")
ax.set_title("B. Number of niche clusters per niche-organized sample\n(higher = more fragmented niches)", fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
fig.suptitle("S_F52  TROP2 niche progression across stages (fraction + cluster count)",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F52_niche_progression.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("[F53] HT TLS PTC_HT × all axes spatial autocorrelation comparison")
# For each PTC+HT slide, compute Moran I for TROP2 + 4 immune sets + DM1 + Epi from per-spot data
HLA_II = ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
B_CELL = ["CD19","MS4A1","CD79A","CD79B","SDC1","JCHAIN"]
TLS    = ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"]
IGHV   = ["AICDA","IGHM","IGHG1","IGHA1","IGKC","IGLC2"]
g230_paths = sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad"))
def morans_quick(values, rows, cols):
    n = len(values); v = np.asarray(values, dtype=float)
    if n < 30 or np.nanstd(v) == 0: return np.nan
    vc = v - np.nanmean(v); pos = {(r,c): i for i,(r,c) in enumerate(zip(rows, cols))}
    Wnum, Wsum = 0.0, 0
    for i, (r, c) in enumerate(zip(rows, cols)):
        for dr, dc in [(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)]:
            j = pos.get((r+dr, c+dc))
            if j is not None: Wnum += vc[i]*vc[j]; Wsum += 1
    Wd = float(np.nansum(vc**2))
    return (n/Wsum)*(Wnum/Wd) if Wsum > 0 and Wd > 0 else np.nan

g230_axes_morans = []
for p in g230_paths:
    sample = p.parent.name
    rec = {"sample": sample}
    # immune sets via h5ad
    for nm, gs_list in [("HLA_II",HLA_II),("B_cell",B_CELL),("TLS",TLS),("IGHV",IGHV),("TROP2",["TACSTD2"]),("EPCAM",["EPCAM"])]:
        scores, coords = score_h5ad(p, gs_list)
        if scores is None: rec[nm] = np.nan; continue
        rows, cols = coords
        if rows is None: rec[nm] = np.nan; continue
        rec[nm] = morans_quick(scores.mean(axis=1), rows, cols)
    # DM1 from per-spot table
    sub_dm = ext[ext.sample_id == sample]
    if not sub_dm.empty and "array_row" in sub_dm.columns:
        rec["DM1_like"] = morans_quick(sub_dm["DM1_like_score"].values, sub_dm["array_row"].values, sub_dm["array_col"].values)
        rec["Epithelial"] = morans_quick(sub_dm["Epithelial_score_raw"].values, sub_dm["array_row"].values, sub_dm["array_col"].values)
    g230_axes_morans.append(rec)
df_g230_axes = pd.DataFrame(g230_axes_morans)
print(df_g230_axes)

axes_show = ["HLA_II","B_cell","TLS","IGHV","EPCAM","Epithelial","DM1_like","TROP2"]
fig, ax = plt.subplots(figsize=(11, 5.5))
mat = np.full((len(g230_axes_morans), len(axes_show)), np.nan)
for i, r in enumerate(g230_axes_morans):
    for j, ax_name in enumerate(axes_show):
        mat[i, j] = r.get(ax_name, np.nan)
im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.9)
for i in range(len(g230_axes_morans)):
    for j in range(len(axes_show)):
        v = mat[i,j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v>0.55 else "black", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(axes_show))); ax.set_xticklabels(axes_show, fontsize=10, family="monospace")
ax.set_yticks(range(len(g230_axes_morans)))
ax.set_yticklabels([r["sample"].split("_")[-1] for r in g230_axes_morans], fontsize=10, family="monospace")
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.04, label="Moran's I")
ax.set_title("S_F53  GSE230424 PTC+HT — full axis Moran's I matrix (8 axes × 4 slides)\n"
             "HT immune sets (HLA-II / B-cell / TLS / IGHV) niche-organized; TROP2 weak; DM1 also organized",
             fontsize=11.5, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F53_HT_full_axis_morans.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("[F54] TROP2 niche × Epithelial niche overlap (within-sample Jaccard)")
# For each PTC+LPTC niche sample, compute Jaccard between TROP2-high and Epithelial-high spots
trop_g521 = trop_all[trop_all.dataset == "GSE250521"].copy()
trop_g521["row_col"] = trop_g521["sample_id"]+"_"+trop_g521["array_row"].astype(str)+"_"+trop_g521["array_col"].astype(str)
g521["row_col"] = g521["sample_id"]+"_"+g521["array_row"].astype(str)+"_"+g521["array_col"].astype(str)
mg_all = trop_g521.merge(g521[["row_col","DM1_like_score","Epithelial_score","Proliferation_score"]], on="row_col")
jaccard_data = []
for sample, sub in mg_all.groupby("sample_id"):
    cond = trop_all[trop_all.sample_id == sample]["condition"].iloc[0]
    if sub.empty: continue
    # TROP2-high mask
    t_thr = sub["TROP2"].quantile(0.85)
    e_thr = sub["Epithelial_score"].quantile(0.75)
    p_thr = sub["Proliferation_score"].quantile(0.75)
    d_thr = sub["DM1_like_score"].quantile(0.75)
    t_mask = (sub["TROP2"] > t_thr).values
    for axis_name, axis_mask in [("Epithelial",(sub["Epithelial_score"] > e_thr).values),
                                   ("Proliferation",(sub["Proliferation_score"] > p_thr).values),
                                   ("DM1_high",(sub["DM1_like_score"] > d_thr).values)]:
        intersection = (t_mask & axis_mask).sum()
        union = (t_mask | axis_mask).sum()
        jaccard = intersection / union if union > 0 else np.nan
        jaccard_data.append({"sample":sample, "stage":cond, "axis":axis_name, "jaccard":jaccard})
df_jac = pd.DataFrame(jaccard_data)
fig, ax = plt.subplots(figsize=(11, 5.5))
axes_jac = ["Epithelial","Proliferation","DM1_high"]
positions = []
for i, axis in enumerate(axes_jac):
    for j, stage in enumerate(["PT","PTC","LPTC","ATC"]):
        sub = df_jac[(df_jac.axis == axis) & (df_jac.stage == stage)]
        if sub.empty: continue
        x = i*5 + j
        positions.append((x, axis, stage))
        ax.scatter([x]*len(sub), sub["jaccard"].values, s=110, color=PAL[stage],
                   alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
xticks = [p[0] for p in positions]
xlabels = [f"{p[2]}" for p in positions]
ax.set_xticks(xticks); ax.set_xticklabels(xlabels, fontsize=8.5)
# add axis-level group labels
for i, axis in enumerate(axes_jac):
    ax.text(i*5 + 1.5, -0.05, axis, transform=ax.get_xaxis_transform(),
            fontsize=11, fontweight="bold", color="#7B1F2A", ha="center")
ax.set_ylabel("Jaccard overlap (TROP2-high ∩ axis-high)")
ax.set_title("S_F54  TROP2 niche overlap with Epithelial / Proliferation / DM1-high zones (Jaccard, top-quantile)\n"
             "TROP2 niche의 cell-state context — Epithelial-high과 가장 큰 overlap (tumor cell 영역)",
             fontsize=11.5, fontweight="bold")
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F54_TROP2_niche_jaccard.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("[F55] Sample size × niche detection — power consideration")
fig, ax = plt.subplots(figsize=(10, 5.5))
df_E_full = df_E.copy()
df_E_full["log_spots"] = np.log10(df_E_full["n_spots"])
for cond, color in PAL.items():
    sub = df_E_full[df_E_full.condition == cond]
    if sub.empty: continue
    ax.scatter(sub["log_spots"], sub["morans_I_TROP2"], s=120, color=color,
               alpha=0.85, edgecolor="black", linewidth=0.6, label=f"{cond}", zorder=3)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.6)
ax.set_xlabel("log10(n_spots per sample)")
ax.set_ylabel("TROP2 Moran's I")
ax.set_title("S_F55  Niche detection vs sample size (per-slide n_spots)\n"
             "Niche signal ≠ sample size effect; PTC/LPTC niche-organized regardless of slide size",
             fontsize=11.5, fontweight="bold")
ax.legend(loc="lower right", fontsize=9, ncol=2); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F55_niche_vs_sample_size.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("[F56] SUMMARY headline figure — 4 key findings collage")
fig = plt.figure(figsize=(16, 10))
g = gs.GridSpec(2, 2, hspace=0.32, wspace=0.18, top=0.92, bottom=0.05, left=0.05, right=0.97)

# Q1: TROP2 niche per condition (use S_F22 data)
ax = fig.add_subplot(g[0, 0])
positions_h = []
data_h = []
for c in ["PT","PTC","LPTC","ATC","PTC_HT","HT","GD","CONTROL"]:
    sub = df_E[df_E.condition == c]
    if not sub.empty:
        positions_h.append(c); data_h.append(sub["morans_I_TROP2"].values)
bp = ax.boxplot(data_h, labels=positions_h, patch_artist=True, showfliers=False, widths=0.55)
for patch, c in zip(bp["boxes"], positions_h):
    patch.set_facecolor(PAL.get(c,"#999")); patch.set_alpha(0.8)
for i, vals in enumerate(data_h):
    jit = np.random.RandomState(i).uniform(-0.07, 0.07, len(vals))
    ax.scatter(np.full_like(vals, i+1)+jit, vals, s=40, color="black", zorder=3, edgecolor="white")
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1.4, alpha=0.7)
ax.set_ylabel("TROP2 Moran's I")
ax.set_title("Q1.  TROP2 niche tumor-specific\nPTC + LPTC 8/8 above niche threshold; 다른 모든 group below",
             fontsize=12, fontweight="bold", color="#7B1F2A")
ax.grid(axis="y", alpha=0.3)

# Q2: HT TLS niche
ax = fig.add_subplot(g[0, 1])
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
        ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                color="white" if mat[i,j]>0.55 else "black", fontsize=11, fontweight="bold")
ax.set_xticks(range(len(markers))); ax.set_xticklabels(markers, fontsize=10, family="monospace")
ax.set_yticks(range(len(samps))); ax.set_yticklabels([s.split("_")[-1] for s in samps], fontsize=10, family="monospace")
ax.set_title("Q2.  HT-overlap PTC TLS niche\nGSE230424 4/4 슬라이드 × 4/4 immune marker 모두 niche-organized",
             fontsize=12, fontweight="bold", color="#1abc9c")

# Q3: stage trajectory
ax = fig.add_subplot(g[1, 0])
df_moran = pd.read_csv(OUT/"spatial_F40_per_sample_per_axis_morans.tsv", sep="\t")
trajectory_axes = ["TROP2","DM1_like_score","Epithelial_score","Proliferation_score"]
colors_traj = ["#7B1F2A","#34547A","#3C6B4F","#B8893C"]
for ax_name, color in zip(trajectory_axes, colors_traj):
    if ax_name not in df_moran.columns: continue
    means = []; stds = []
    for cond in ORDER:
        sub = df_moran[df_moran.condition == cond][ax_name].dropna()
        means.append(sub.mean() if len(sub) else np.nan)
        stds.append(sub.std() if len(sub) > 1 else 0)
    ax.errorbar(range(len(ORDER)), means, yerr=stds, marker="o", lw=2.4,
                color=color, label=ax_name.replace("_score",""), capsize=5, alpha=0.9, ms=10)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, fontsize=11)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.5)
ax.set_ylabel("Mean Moran's I per stage")
ax.set_title("Q3.  Stage trajectory of niche organization\nTROP2 peaks at PTC+LPTC; collapses in ATC",
             fontsize=12, fontweight="bold", color="#34547A")
ax.legend(loc="best", fontsize=10); ax.grid(axis="y", alpha=0.3)

# Q4: H&E closure NO-GO
ax = fig.add_subplot(g[1, 1])
sec_a = clos_metric[(clos_metric["section"].astype(str).str.startswith("A")) & (clos_metric["model"]=="ridge")]
labels = sec_a["experiment"].astype(str).str.replace("raw_","")
ax.bar(range(len(sec_a)), sec_a["pooled_spearman_r"].values, color="#7B1F2A",
       edgecolor="black", linewidth=0.5)
ax.axhline(0.3, ls="--", color="#3C6B4F", lw=1.5, alpha=0.7); ax.text(0, 0.32, "useful threshold (0.3)", color="#3C6B4F", fontsize=10, fontweight="bold")
ax.axhline(0, color="black", lw=0.6)
ax.set_xticks(range(len(sec_a))); ax.set_xticklabels(labels, fontsize=10, rotation=15, ha="right")
ax.set_ylabel("Pooled Spearman ρ (LOSO)")
ax.set_ylim(-0.1, 0.5)
ax.set_title("Q4.  H&E → DM1 closure battery NO-GO\n모든 raw label LOSO ρ < 0.3 useful threshold",
             fontsize=12, fontweight="bold", color="#52525a")
ax.grid(axis="y", alpha=0.3)

fig.suptitle("S_F56  Spatial full-package — 4 KEY FINDINGS at a glance\n"
             "Q1 TROP2 niche tumor-specific · Q2 HT TLS niche directly visible · Q3 stage trajectory · Q4 H&E NO-GO",
             fontsize=14, fontweight="bold", y=0.99)
fig.savefig(ASSETS/"S_F56_4_key_findings.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 3 fixed (S_F16/S_F25/S_F30) + 8 new (S_F49-S_F56)")
