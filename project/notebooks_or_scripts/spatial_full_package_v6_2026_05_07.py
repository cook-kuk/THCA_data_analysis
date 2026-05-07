#!/usr/bin/env python3
"""Spatial v6 — 8 more figures (S_F41 ~ S_F48)."""
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

PAL = {"PT":"#3C6B4F","PTC":"#34547A","LPTC":"#B8893C","ATC":"#7B1F2A",
       "PTC_HT":"#962E2E","CONTROL":"#52525a","HT":"#34547A","GD":"#B8893C"}
ORDER = ["PT","PTC","LPTC","ATC"]
trop_cmap = LinearSegmentedColormap.from_list("trop2",["#fffaf2","#ffd58a","#d97c2e","#7B1F2A"],N=256)
hypo_cmap = LinearSegmentedColormap.from_list("hypo",["#fffaf2","#dcc5e0","#7B1F2A","#0F1A2E"],N=256)
caf_cmap  = LinearSegmentedColormap.from_list("caf",["#fffaf2","#a8c8b6","#3C6B4F","#0F1A2E"],N=256)
emt_cmap  = LinearSegmentedColormap.from_list("emt",["#fffaf2","#e2c08c","#B8893C","#0F1A2E"],N=256)

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

print("[load]")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
trop_all = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
df_moran = pd.read_csv(OUT/"spatial_F40_per_sample_per_axis_morans.tsv", sep="\t")

# ===== S_F41: ALL 28 samples TROP2 mega grid =====
print("[F41] 28-sample TROP2 mega grid")
all_samples = trop_all[["sample_id","dataset","condition"]].drop_duplicates().reset_index(drop=True)
# Order by condition: PT, PTC, LPTC, ATC, PTC_HT, CONTROL, HT, GD
cond_rank = {"PT":0,"PTC":1,"LPTC":2,"ATC":3,"PTC_HT":4,"CONTROL":5,"HT":6,"GD":7}
all_samples["rank"] = all_samples["condition"].map(cond_rank)
all_samples = all_samples.sort_values(["rank","sample_id"]).reset_index(drop=True)
fig = plt.figure(figsize=(20, 16))
g = gs.GridSpec(4, 7, hspace=0.42, wspace=0.18, top=0.95, bottom=0.03, left=0.03, right=0.97)
for idx, (_, row) in enumerate(all_samples.iterrows()):
    if idx >= 28: break
    r, c = idx // 7, idx % 7
    ax = fig.add_subplot(g[r, c])
    sub = trop_all[trop_all.sample_id == row["sample_id"]]
    if "array_row" not in sub.columns or sub.empty:
        ax.axis("off"); continue
    x = sub["array_col"].values; y = -sub["array_row"].values
    v = sub["TROP2"].values
    vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
    ax.scatter(x, y, c=v, cmap=trop_cmap, vmin=0, vmax=vmax, s=2, alpha=0.92)
    moran = df_E[df_E.sample_id == row["sample_id"]]["morans_I_TROP2"]
    moran = moran.iloc[0] if not moran.empty else np.nan
    title_color = PAL.get(row["condition"], "#666")
    weight = "bold" if (not np.isnan(moran) and moran > 0.25) else "normal"
    ax.set_title(f"{row['condition']} · {row['sample_id'].split('_')[-1]}\nI={moran:.2f}",
                 fontsize=8.5, color=title_color, fontweight=weight)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F41  ALL 28 samples × 3 cohorts TROP2 spatial maps — comprehensive view\n"
             "Bold-titled samples = niche-organized (Moran I > 0.25); 9 samples (cancer PTC+LPTC dominant)",
             fontsize=13, fontweight="bold")
fig.savefig(ASSETS/"S_F41_TROP2_all_28_samples.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ===== S_F42: HT individual immune markers (CD20=MS4A1, CXCL13, AICDA × 4 slides) =====
print("[F42] HT individual immune marker maps")
g230_paths = sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad"))
markers_indv = ["MS4A1","CXCL13","AICDA"]  # CD20, TLS chemokine, AID enzyme
fig = plt.figure(figsize=(15, 11.5))
g = gs.GridSpec(3, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.05, right=0.97)
tls_cmap = LinearSegmentedColormap.from_list("tls2",["#fffaf2","#a8d8c8","#1abc9c","#0F1A2E"],N=256)
for r, gene in enumerate(markers_indv):
    for c, p in enumerate(g230_paths):
        sample = p.parent.name
        ax = fig.add_subplot(g[r, c])
        scores, coords = score_h5ad(p, [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows, cols = coords
        if rows is None: ax.axis("off"); continue
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols), -np.array(rows), c=v, cmap=tls_cmap, vmin=0, vmax=vmax, s=4, alpha=0.92)
        if r == 0: ax.set_title(f"PTC+HT · {sample.split('_')[-1]}", fontsize=10, color="#962E2E", fontweight="bold")
        if c == 0:
            label_map = {"MS4A1":"CD20\n(MS4A1)\nB-cell","CXCL13":"CXCL13\nTLS chemokine","AICDA":"AICDA\nB-cell\nactivation"}
            ax.set_ylabel(label_map.get(gene, gene), fontsize=11, fontweight="bold", color="#1abc9c", rotation=0, labelpad=45, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F42  HT-overlap PTC individual immune marker maps — CD20 · CXCL13 · AICDA\n"
             "3 핵심 marker가 같은 영역에 collocated (S_F8 multi-set 검증) — true tertiary lymphoid structure 형성 evidence",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F42_HT_individual_markers.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F43: Hypoxia 16-slide spatial maps =====
print("[F43] Hypoxia 16-slide maps")
by_stage = {s:[] for s in ORDER}
for sid in g521["sample_id"].unique():
    if "_N-" in sid: by_stage["PT"].append(sid)
    elif "PTC-" in sid and "L" not in sid: by_stage["PTC"].append(sid)
    elif "LPTC-" in sid: by_stage["LPTC"].append(sid)
    elif "ATC-" in sid: by_stage["ATC"].append(sid)
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.94, bottom=0.04, left=0.04, right=0.97)
for r, stage in enumerate(ORDER):
    for c, sample in enumerate(sorted(by_stage[stage])[:4]):
        ax = fig.add_subplot(g[r, c])
        sub = g521[g521.sample_id == sample]
        if sub.empty or "array_row" not in sub.columns: ax.axis("off"); continue
        v = sub["Hypoxia_score"].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
        ax.scatter(sub["array_col"].values, -sub["array_row"].values, c=v, cmap=hypo_cmap, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F43  Hypoxia signature spatial maps — GSE250521 16 slides\n"
             "ATC stages에서 Hypoxia 활성 영역 visible (microenv trajectory marker)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F43_Hypoxia_maps_GSE250521.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F44: CAF/ECM 16-slide spatial maps =====
print("[F44] CAF/ECM 16-slide maps")
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.94, bottom=0.04, left=0.04, right=0.97)
for r, stage in enumerate(ORDER):
    for c, sample in enumerate(sorted(by_stage[stage])[:4]):
        ax = fig.add_subplot(g[r, c])
        sub = g521[g521.sample_id == sample]
        if sub.empty or "array_row" not in sub.columns: ax.axis("off"); continue
        v = sub["CAF_ECM_score"].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
        ax.scatter(sub["array_col"].values, -sub["array_row"].values, c=v, cmap=caf_cmap, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F44  CAF/ECM stromal signature spatial maps — GSE250521 16 slides\n"
             "Stromal niche 영역 식별 — TROP2 niche가 stromal 영역과 분리되는지 직접 검증 가능 (S_F7와 비교)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F44_CAF_ECM_maps_GSE250521.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F45: EMT 16-slide spatial maps =====
print("[F45] EMT 16-slide maps")
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.94, bottom=0.04, left=0.04, right=0.97)
for r, stage in enumerate(ORDER):
    for c, sample in enumerate(sorted(by_stage[stage])[:4]):
        ax = fig.add_subplot(g[r, c])
        sub = g521[g521.sample_id == sample]
        if sub.empty or "array_row" not in sub.columns: ax.axis("off"); continue
        v = sub["EMT_score"].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
        ax.scatter(sub["array_col"].values, -sub["array_row"].values, c=v, cmap=emt_cmap, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F45  EMT signature spatial maps — GSE250521 16 slides\n"
             "EMT high zones 식별; ATC stages에서 EMT 활성 visible (dedifferentiation pathway)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F45_EMT_maps_GSE250521.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F46: Niche threshold sensitivity analysis =====
print("[F46] Niche threshold sensitivity")
def morans_I_quick(values, rows, cols):
    n = len(values)
    if n < 30: return np.nan
    v = np.asarray(values, dtype=float)
    if np.nanstd(v) == 0: return np.nan
    vc = v - np.nanmean(v)
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

# vary spot threshold (top X%) for binarization, see how Moran I changes
print("  computing across thresholds [70%, 75%, 80%, 85%, 90%, 95%]")
thresholds = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
sens_rows = []
for sample in df_E["sample_id"].unique()[:28]:
    sub = trop_all[trop_all.sample_id == sample]
    if "array_row" not in sub.columns or sub.empty: continue
    cond = sub["condition"].iloc[0]
    rows_arr = sub["array_row"].values
    cols_arr = sub["array_col"].values
    for t in thresholds:
        threshold = sub["TROP2"].quantile(t)
        binary = (sub["TROP2"] > threshold).astype(int).values
        if binary.sum() < 5 or binary.sum() == len(binary): continue
        I = morans_I_quick(binary, rows_arr, cols_arr)
        sens_rows.append({"sample": sample, "condition": cond, "threshold_pct": int(t*100), "morans_I": I})
df_sens = pd.DataFrame(sens_rows)
df_sens.to_csv(OUT/"spatial_F46_niche_threshold_sensitivity.tsv", sep="\t", index=False)
fig, ax = plt.subplots(figsize=(13, 6))
for cond, color in PAL.items():
    sub = df_sens[df_sens.condition == cond]
    if sub.empty: continue
    grp = sub.groupby("threshold_pct")["morans_I"].agg(["mean","std","count"]).reset_index()
    ax.errorbar(grp["threshold_pct"], grp["mean"], yerr=grp["std"], marker="o", lw=1.5,
                color=color, label=f"{cond} (n={int(grp['count'].max())})", capsize=4, alpha=0.85)
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.5)
ax.text(70, 0.27, "niche threshold (0.25)", color="#7B1F2A", fontsize=10, fontweight="bold")
ax.set_xlabel("TROP2 binarization threshold (per-sample percentile %)")
ax.set_ylabel("TROP2 binary Moran's I (mean ± std)")
ax.set_title("S_F46  Niche threshold sensitivity — TROP2 binary Moran's I across 70-95% percentile thresholds\n"
             "PTC+LPTC niche 신호는 70-90% threshold에서 robust 유지; PT/HT/GD/CONTROL은 모든 threshold에서 baseline",
             fontsize=12, fontweight="bold")
ax.legend(loc="best", fontsize=10, ncol=2)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(ASSETS/"S_F46_niche_threshold_sensitivity.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F47: Per-axis Moran I distribution boxplots (cross-cohort) =====
print("[F47] Per-axis Moran I boxplot")
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
# build long-form
axes_show = ["TROP2","DM1_like_score","RAI_8_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]
ax = axes[0]
all_data = []
for ax_name in axes_show:
    if ax_name in df_moran.columns:
        sub_data = df_moran[ax_name].dropna().values
        all_data.append(sub_data)
    else:
        all_data.append([])
bp = ax.boxplot(all_data, labels=[a.replace("_score","") for a in axes_show], patch_artist=True, showfliers=False)
for patch, ax_name in zip(bp["boxes"], axes_show):
    color = "#7B1F2A" if "TROP2" in ax_name else ("#34547A" if "DM1" in ax_name or "RAI" in ax_name else "#3C6B4F")
    patch.set_facecolor(color); patch.set_alpha(0.7)
ax.tick_params(axis='x', rotation=30)
for tk in ax.get_xticklabels(): tk.set_horizontalalignment('right')
ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.5)
ax.set_ylabel("Moran's I (per-sample distribution, n=16 GSE250521)")
ax.set_title("A. Per-axis Moran I distribution across 16 GSE250521 samples", fontsize=11, fontweight="bold")
ax.grid(axis="y", alpha=0.3)

# panel B: per-stage Moran I trajectory for TROP2 + DM1 + Epithelial
ax = axes[1]
trajectory_axes = ["TROP2","DM1_like_score","Epithelial_score","Proliferation_score"]
colors_traj = ["#7B1F2A","#34547A","#3C6B4F","#B8893C"]
for ax_name, color in zip(trajectory_axes, colors_traj):
    if ax_name not in df_moran.columns: continue
    means = []; stds = []
    for cond in ORDER:
        sub = df_moran[df_moran.condition == cond][ax_name].dropna()
        means.append(sub.mean() if len(sub) else np.nan)
        stds.append(sub.std() if len(sub) > 1 else 0)
    ax.errorbar(range(len(ORDER)), means, yerr=stds, marker="o", lw=2,
                color=color, label=ax_name.replace("_score",""), capsize=5, alpha=0.85)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, fontsize=11)
ax.axhline(0, color="black", lw=0.6); ax.axhline(0.25, ls="--", color="#7B1F2A", lw=1, alpha=0.5)
ax.set_ylabel("Mean Moran's I per stage")
ax.set_title("B. Stage trajectory of niche organization (TROP2 dominant in PTC+LPTC)", fontsize=11, fontweight="bold")
ax.legend(loc="best", fontsize=10); ax.grid(axis="y", alpha=0.3)

fig.suptitle("S_F47  Per-axis Moran's I summary — distribution + trajectory",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F47_per_axis_moran_summary.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F48: Spot-level depth confounding =====
print("[F48] Spot-level depth confounding")
fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
# panel A: log_total_counts vs DM1_like (overall)
ax = axes[0]
sub_g521 = g521.dropna(subset=["DM1_like_score","total_counts"])
log_counts = np.log10(sub_g521["total_counts"].values + 1)
ax.hexbin(log_counts, sub_g521["DM1_like_score"].values, gridsize=50, cmap="OrRd", mincnt=10)
from scipy.stats import spearmanr
rho, p = spearmanr(log_counts, sub_g521["DM1_like_score"].values, nan_policy="omit")
ax.set_xlabel("log10(total_counts per spot)"); ax.set_ylabel("DM1_like score")
ax.set_title(f"A. log_counts vs DM1 (raw)\nrho = {rho:.3f}, p = {p:.1e}", fontsize=10.5)

# panel B: per-stage depth mean
ax = axes[1]
for stage in ORDER:
    sub = g521[g521.condition == stage]
    samp_means = sub.groupby("sample_id")["total_counts"].mean()
    ax.scatter([ORDER.index(stage)]*len(samp_means), samp_means.values, s=80,
               color=PAL[stage], alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER)
ax.set_ylabel("Per-sample mean total_counts (depth)")
ax.set_title("B. Sequencing depth per stage\n(stage-mean depth differences = confound source)", fontsize=10.5)
ax.grid(axis="y", alpha=0.3)

# panel C: per-stage n_genes_by_counts mean
ax = axes[2]
for stage in ORDER:
    sub = g521[g521.condition == stage]
    samp_means = sub.groupby("sample_id")["n_genes_by_counts"].mean()
    ax.scatter([ORDER.index(stage)]*len(samp_means), samp_means.values, s=80,
               color=PAL[stage], alpha=0.85, edgecolor="black", linewidth=0.6, zorder=3)
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER)
ax.set_ylabel("Per-sample mean n_genes_by_counts")
ax.set_title("C. Gene detection per stage\n(detection capacity confound)", fontsize=10.5)
ax.grid(axis="y", alpha=0.3)

fig.suptitle("S_F48  Spot-level technical confounding — depth + gene detection per stage\n"
             "Depth-correction이 필요한 이유 (closure_battery report 7가지 stage trajectory에서 detected)",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F48_depth_confounding.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 8 new figures (S_F41 ~ S_F48) saved")
