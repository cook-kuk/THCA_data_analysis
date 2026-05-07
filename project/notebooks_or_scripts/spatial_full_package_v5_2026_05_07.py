#!/usr/bin/env python3
"""Spatial v5 — 8 more figures (S_F33 ~ S_F40)."""
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
epi_cmap  = LinearSegmentedColormap.from_list("epi",["#fffaf2","#a8c8b6","#3C6B4F","#0F1A2E"],N=256)
prol_cmap = LinearSegmentedColormap.from_list("prol",["#fffaf2","#a8c8e2","#34547A","#0F1A2E"],N=256)
hypo_cmap = LinearSegmentedColormap.from_list("hypo",["#fffaf2","#ddc5e0","#7B1F2A","#0F1A2E"],N=256)

def score_h5ad(p, gs_list):
    a = ad.read_h5ad(p)
    rvar = a.raw.var_names.astype(str) if a.raw is not None else a.var_names.astype(str)
    X = a.raw.X if a.raw is not None else a.X
    idx = [list(rvar.values).index(g) for g in gs_list if g in rvar.values]
    if not idx: return None, [], None
    sub = X[:, idx]
    if hasattr(sub,"toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", np.ones(a.n_obs))).flatten()
    tot = np.where(tot==0,1,tot)
    norm_log = np.log1p(sub/tot[:,None]*1e4)
    avail = [g for g in gs_list if g in rvar.values]
    coords = (a.obs["array_row"].values, a.obs["array_col"].values) if "array_row" in a.obs.columns else (None, None)
    return norm_log, avail, coords

print("[load]")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
trop_all = pd.read_csv(OUT/"spatial_per_spot_TROP2.tsv.gz", sep="\t")
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
df_C = pd.read_csv(OUT/"spatial_C_HT_TLS_spatial.tsv", sep="\t")

# representative per stage
rep = {}
for s in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==s)]
    if not sub.empty:
        rep[s] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]

g521_paths = {p.parent.name: p for p in sorted(G521H.glob("*/GSM*.raw.h5ad"))}

# ===== S_F33: RAI_8 individual gene spatial maps (4 stages x 4 representative genes) =====
print("[F33] RAI_8 individual gene maps")
RAI8_REP = ["TG","TPO","SLC5A5","FOXE1"]
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.05, right=0.97)
for r, gene in enumerate(RAI8_REP):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, avail, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None: ax.axis("off"); continue
        rows, cols = coords
        if rows is None: ax.axis("off"); continue
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols), -np.array(rows), c=v, cmap=epi_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#3C6B4F", rotation=0, labelpad=35, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F33  RAI_8 individual gene RNA spatial maps (TG · TPO · SLC5A5 · FOXE1)\n"
             "PT에서 thyroid 분화 마커 모두 high; ATC에서 거의 사라짐 — RAI_8 axis loss visible per gene",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F33_RAI8_individual_gene_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F34: GSE250521 16-slide Epithelial maps =====
print("[F34] Epithelial 16-slide maps")
by_stage = {s:[] for s in ORDER}
for sid in g521["sample_id"].unique():
    if "_N-" in sid: by_stage["PT"].append(sid)
    elif "PTC-" in sid and "L" not in sid: by_stage["PTC"].append(sid)
    elif "LPTC-" in sid: by_stage["LPTC"].append(sid)
    elif "ATC-" in sid: by_stage["ATC"].append(sid)

fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.94, bottom=0.04, left=0.04, right=0.97)
for r, stage in enumerate(ORDER):
    samples = sorted(by_stage[stage])[:4]
    for c, sample in enumerate(samples):
        ax = fig.add_subplot(g[r, c])
        sub = g521[g521.sample_id == sample]
        if "array_row" not in sub.columns or sub.empty: ax.axis("off"); continue
        x = sub["array_col"].values; y = -sub["array_row"].values
        v = sub["Epithelial_score"].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
        ax.scatter(x, y, c=v, cmap=epi_cmap, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F34  Epithelial signature spatial maps — GSE250521 16 slides (S_F7 TROP2 / S_F19 DM1 parallel)\n"
             "Epithelial-high region이 tumor cell 영역 정의 — TROP2 niche가 epithelial 영역 안에 있는지 직접 비교 가능",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F34_Epithelial_maps_GSE250521.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F35: GSE250521 16-slide Proliferation maps =====
print("[F35] Proliferation 16-slide maps")
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.94, bottom=0.04, left=0.04, right=0.97)
for r, stage in enumerate(ORDER):
    samples = sorted(by_stage[stage])[:4]
    for c, sample in enumerate(samples):
        ax = fig.add_subplot(g[r, c])
        sub = g521[g521.sample_id == sample]
        if "array_row" not in sub.columns or sub.empty: ax.axis("off"); continue
        x = sub["array_col"].values; y = -sub["array_row"].values
        v = sub["Proliferation_score"].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
        ax.scatter(x, y, c=v, cmap=prol_cmap, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F35  Proliferation spatial maps — GSE250521 16 slides\n"
             "Proliferation gradient는 closure battery에서 depth-corrected ρ=+0.58 (positive control 통과한 axis)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F35_Proliferation_maps_GSE250521.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F36: GSE248205 full 8-slide TROP2 + DM1 maps =====
print("[F36] GSE248205 full 8-slide maps")
fig = plt.figure(figsize=(20, 8))
g = gs.GridSpec(2, 8, hspace=0.32, wspace=0.18, top=0.92, bottom=0.04, left=0.03, right=0.99)
g248 = ext = pd.concat([pd.read_csv(f, sep="\t") for f in (ROOT/"project_external_st/results/scores").glob("GSE248205_*.tsv.gz")], axis=0, ignore_index=True).rename(columns={"DM1_like_score_raw":"DM1_like_score"})
g248_samples = sorted(g248["sample_id"].unique())
for c, sample in enumerate(g248_samples):
    sub = g248[g248.sample_id == sample]
    cond = sub["condition_inferred"].iloc[0] if not sub.empty else "?"
    if "array_row" not in sub.columns or sub.empty: continue
    # row 0: DM1
    ax = fig.add_subplot(g[0, c])
    x = sub["array_col"].values; y = -sub["array_row"].values
    v = sub["DM1_like_score"].values
    vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
    ax.scatter(x, y, c=v, cmap=dm_cmap, vmin=-vmax, vmax=vmax, s=4, alpha=0.92)
    ax.set_title(f"{cond} · {sample.split('_')[-1]}", fontsize=10, color=PAL.get(cond,"#666"), fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    if c == 0: ax.set_ylabel("DM1_like", fontsize=10, color="#7B1F2A", fontweight="bold")
    # row 1: TROP2
    ax = fig.add_subplot(g[1, c])
    sub_t = trop_all[trop_all.sample_id == sample]
    if not sub_t.empty:
        v = sub_t["TROP2"].values
        vmax_t = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(sub_t["array_col"].values, -sub_t["array_row"].values, c=v, cmap=trop_cmap, vmin=0, vmax=vmax_t, s=4, alpha=0.92)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    if c == 0: ax.set_ylabel("TROP2", fontsize=10, color="#7B1F2A", fontweight="bold")
fig.suptitle("S_F36  GSE248205 autoimmune-only full 8 슬라이드 — DM1_like (top) · TROP2 (bottom) — comprehensive negative control\n"
             "8/8 슬라이드 모두 niche structure 거의 없음 → TROP2 niche specificity = tumor-only",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F36_GSE248205_full_maps.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ===== S_F37: HT combined immune density map per slide =====
print("[F37] HT combined immune niche density")
HLA_II = ["HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"]
B_CELL = ["CD19","MS4A1","CD79A","CD79B","SDC1","JCHAIN"]
TLS    = ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL"]
IGHV   = ["AICDA","IGHM","IGHG1","IGHA1","IGKC","IGLC2"]
g230_paths = sorted((EXTH/"GSE230424").glob("*/GSM*.raw.h5ad"))
fig = plt.figure(figsize=(15, 4.4))
g = gs.GridSpec(1, 4, wspace=0.18, top=0.85, bottom=0.04, left=0.04, right=0.97)
combined_cmap = LinearSegmentedColormap.from_list("combined",["#fffaf2","#a8d8c8","#1abc9c","#0F1A2E"],N=256)
for c, p in enumerate(g230_paths):
    sample = p.parent.name
    ax = fig.add_subplot(g[0, c])
    combined = []
    for genes in [HLA_II, B_CELL, TLS, IGHV]:
        scores, avail, coords = score_h5ad(p, genes)
        if scores is None: continue
        v = scores.mean(axis=1)
        # z-score within sample
        if v.std() > 0: v = (v - v.mean()) / v.std()
        combined.append(v)
    if not combined: continue
    rows, cols = coords
    avg = np.mean(np.stack(combined, axis=0), axis=0)
    x = np.array(cols); y = -np.array(rows)
    vmax = np.nanpercentile(avg, 99); vmin = np.nanpercentile(avg, 1)
    ax.scatter(x, y, c=avg, cmap=combined_cmap, vmin=vmin, vmax=vmax, s=4, alpha=0.92)
    ax.set_title(f"PTC+HT · {sample.split('_')[-1]}\nCombined immune z (4 sets avg)", fontsize=10, color="#962E2E", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F37  Combined HT immune niche density — HLA-II + B-cell + TLS + IGHV/AICDA averaged per spot (z-scored within sample)\n"
             "다크 영역 = 4 immune signature 모두 공동 활성 = TLS-organized region",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F37_HT_combined_immune_density.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F38: per-sample QC summary =====
print("[F38] Per-sample QC summary")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
qc_data = []
for sample in trop_all["sample_id"].unique():
    sub = trop_all[trop_all.sample_id == sample]
    cond = sub["condition"].iloc[0]
    ds = sub["dataset"].iloc[0]
    qc_data.append({
        "sample": sample, "dataset": ds, "condition": cond,
        "n_spots": len(sub),
    })
df_qc = pd.DataFrame(qc_data)
# panel A: spot count per sample
ax = axes[0]
df_qc_sorted = df_qc.sort_values(["dataset","condition","n_spots"]).reset_index(drop=True)
xs = np.arange(len(df_qc_sorted))
ax.bar(xs, df_qc_sorted["n_spots"].values,
       color=[PAL.get(c,"#999") for c in df_qc_sorted["condition"]], edgecolor="black", linewidth=0.4)
ax.set_xticks(xs)
ax.set_xticklabels([f"{r['sample'].split('_')[-1]}\n[{r['condition']}]" for _,r in df_qc_sorted.iterrows()],
                   rotation=70, ha="right", fontsize=7.5)
ax.set_ylabel("n_spots per sample"); ax.set_title("A. Spot count per sample", fontsize=10.5)
ax.grid(axis="y", alpha=0.3)
# panel B: per-condition mean spot count
ax = axes[1]
cond_avg = df_qc.groupby(["dataset","condition"])["n_spots"].mean().reset_index()
labels = [f"{r['dataset'][:7]}\n{r['condition']}" for _, r in cond_avg.iterrows()]
ax.bar(range(len(cond_avg)), cond_avg["n_spots"].values,
       color=[PAL.get(r["condition"],"#999") for _,r in cond_avg.iterrows()], edgecolor="black", linewidth=0.4)
ax.set_xticks(range(len(cond_avg))); ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
ax.set_ylabel("Mean n_spots per condition"); ax.set_title("B. Mean spots per condition", fontsize=10.5)
ax.grid(axis="y", alpha=0.3)
# panel C: TROP2 gene coverage per sample
ax = axes[2]
df_qc["TROP2_nonzero_pct"] = [
    (trop_all[(trop_all.sample_id == s) & (trop_all["TROP2"] > 0.001)].shape[0] / df_qc[df_qc["sample"] == s]["n_spots"].iloc[0] * 100)
    for s in df_qc["sample"]
]
df_qc_sorted2 = df_qc.sort_values(["condition","TROP2_nonzero_pct"]).reset_index(drop=True)
ax.bar(np.arange(len(df_qc_sorted2)), df_qc_sorted2["TROP2_nonzero_pct"].values,
       color=[PAL.get(c,"#999") for c in df_qc_sorted2["condition"]], edgecolor="black", linewidth=0.4)
ax.set_xticks(np.arange(len(df_qc_sorted2)))
ax.set_xticklabels([f"{r['sample'].split('_')[-1]}\n[{r['condition']}]" for _,r in df_qc_sorted2.iterrows()],
                   rotation=70, ha="right", fontsize=7.5)
ax.set_ylabel("% spots with TROP2 > 0"); ax.set_title("C. TROP2 detection coverage per sample", fontsize=10.5)
ax.grid(axis="y", alpha=0.3)
fig.suptitle("S_F38  Per-sample QC summary — spot count + condition averages + TROP2 detection coverage",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(ASSETS/"S_F38_per_sample_QC.png", dpi=160, bbox_inches="tight"); plt.close(fig)

# ===== S_F39: Microenvironment trajectory (CAF/ECM + Hypoxia + EMT per stage) =====
print("[F39] Microenvironment trajectory")
fig, axes = plt.subplots(3, 4, figsize=(15, 11), sharex=True, sharey=True)
ENV_AXES = [("CAF_ECM_score","CAF/ECM","#3C6B4F"),
            ("Hypoxia_score","Hypoxia","#7B1F2A"),
            ("EMT_score","EMT","#B8893C")]
for r, (col, label, color) in enumerate(ENV_AXES):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample: continue
        ax = axes[r, c]
        sub = g521[g521.sample_id == sample]
        if "array_row" not in sub.columns: ax.axis("off"); continue
        x = sub["array_col"].values; y = -sub["array_row"].values
        v = sub[col].values
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
        cmap_use = LinearSegmentedColormap.from_list("c",["#fffaf2","#cccccc",color],N=256)
        ax.scatter(x, y, c=v, cmap=cmap_use, vmin=-vmax, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage}  ·  {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(label, fontsize=12, fontweight="bold", color=color, rotation=0, labelpad=35, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F39  Tumor microenvironment trajectory — CAF/ECM · Hypoxia · EMT stage representative\n"
             "ATC stage에서 EMT/Hypoxia +; LPTC stage에서 CAF/ECM stromal 활성화 visible",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F39_microenv_trajectory.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F40: Cross-axis Moran I correlation matrix =====
print("[F40] Cross-axis Moran's I correlation")
# For each GSE250521 sample, compute Moran I for each axis (DM1_like, RAI_8, Epithelial, Proliferation, CAF_ECM, EMT, Hypoxia)
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

print("  computing per-sample Moran I across 7 axes (16 GSE250521 samples)")
moran_rows = []
for sample in sorted(by_stage["PT"]+by_stage["PTC"]+by_stage["LPTC"]+by_stage["ATC"]):
    sub = g521[g521.sample_id == sample]
    if "array_row" not in sub.columns or sub.empty: continue
    rows_arr = sub["array_row"].values; cols_arr = sub["array_col"].values
    rec = {"sample": sample}
    # condition
    if "_N-" in sample: rec["condition"] = "PT"
    elif "PTC-" in sample and "L" not in sample: rec["condition"] = "PTC"
    elif "LPTC-" in sample: rec["condition"] = "LPTC"
    elif "ATC-" in sample: rec["condition"] = "ATC"
    for axis in ["RAI_8_score","DM1_like_score","TDS_like_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]:
        rec[axis] = morans_I_quick(sub[axis].values, rows_arr, cols_arr)
    # TROP2
    sub_t = trop_all[trop_all.sample_id == sample]
    if not sub_t.empty and "array_row" in sub_t.columns:
        rec["TROP2"] = morans_I_quick(sub_t["TROP2"].values, sub_t["array_row"].values, sub_t["array_col"].values)
    moran_rows.append(rec)
df_moran = pd.DataFrame(moran_rows)
df_moran.to_csv(OUT/"spatial_F40_per_sample_per_axis_morans.tsv", sep="\t", index=False)

axes_to_show = ["TROP2","DM1_like_score","RAI_8_score","Epithelial_score","Proliferation_score","CAF_ECM_score","EMT_score","Hypoxia_score"]
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
# panel A: per-condition mean Moran I per axis
ax = axes[0]
mat = np.full((len(axes_to_show), len(ORDER)), np.nan)
for j, cond in enumerate(ORDER):
    sub = df_moran[df_moran.condition == cond]
    for i, ax_name in enumerate(axes_to_show):
        if ax_name in sub.columns: mat[i, j] = sub[ax_name].mean()
im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=-0.1, vmax=0.7)
for i in range(len(axes_to_show)):
    for j in range(len(ORDER)):
        v = mat[i,j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v>0.4 else "black", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, fontsize=11)
ax.set_yticks(range(len(axes_to_show))); ax.set_yticklabels([a.replace("_score","") for a in axes_to_show], fontsize=10)
ax.set_title("A. Mean Moran's I per axis × stage (GSE250521 16 samples)", fontsize=11, fontweight="bold")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Mean Moran's I")
# panel B: cross-axis Moran I correlation (which axes co-vary in their niche organization)
ax = axes[1]
sub_data = df_moran[axes_to_show].dropna()
cor = sub_data.corr(method="spearman").values
im = ax.imshow(cor, vmin=-1, vmax=1, cmap="RdBu_r")
for i in range(len(axes_to_show)):
    for j in range(len(axes_to_show)):
        ax.text(j, i, f"{cor[i,j]:.2f}", ha="center", va="center",
                color="white" if abs(cor[i,j])>0.6 else "black", fontsize=8.5, fontweight="bold")
ax.set_xticks(range(len(axes_to_show))); ax.set_xticklabels([a.replace("_score","") for a in axes_to_show], fontsize=8.5, rotation=45, ha="right")
ax.set_yticks(range(len(axes_to_show))); ax.set_yticklabels([a.replace("_score","") for a in axes_to_show], fontsize=8.5)
ax.set_title("B. Cross-axis Moran's I correlation (same-sample organization coupling)", fontsize=11, fontweight="bold")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman ρ")
fig.suptitle("S_F40  Per-axis Moran's I — stage trajectory + cross-axis coupling (GSE250521 16 samples × 8 axes)",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(ASSETS/"S_F40_morans_axis_correlation.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 8 new figures (S_F33 ~ S_F40) saved")
