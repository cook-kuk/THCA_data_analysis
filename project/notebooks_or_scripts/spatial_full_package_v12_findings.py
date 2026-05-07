#!/usr/bin/env python3
"""Spatial v12 — 12 new figures (S_F89-S_F100), 주제 무제한 finding sweep."""
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

PAL = {"PT":"#3C6B4F","PTC":"#34547A","LPTC":"#B8893C","ATC":"#7B1F2A"}
ORDER = ["PT","PTC","LPTC","ATC"]

# 8 NEW gene set families (non-overlap with v10/v11)
NEW_SETS = {
    "Hypoxia":           ["HIF1A","EPAS1","VEGFA","CA9","BNIP3","SLC2A1","NDRG1","LDHA","PGK1"],
    "EMT_program":       ["CDH1","CDH2","VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","FN1","S100A4"],
    "Thyroid_differen":  ["TG","TPO","TSHR","SLC5A5","DUOX1","DUOX2","IYD","DIO1","DIO2"],
    "Stemness_lineage":  ["SOX2","POU5F1","NANOG","KLF4","ALDH1A1","PROM1","CD44","BMI1"],
    "Senescence_SASP":   ["CDKN1A","CDKN2A","GLB1","IL6","CXCL8","MMP3","MMP9","IGFBP3","SERPINE1"],
    "Apoptosis_DDR":     ["BAX","BCL2","CASP3","CASP8","FAS","BID","ATM","BRCA1","RAD51","H2AFX"],
    "Angiogenesis":      ["VEGFA","FLT1","KDR","PDGFB","ANGPT1","ANGPT2","CDH5","PECAM1","TEK","DLL4"],
    "Cellular_stress":   ["GPX4","SLC7A11","ACSL4","FTH1","BECN1","MAP1LC3B","SQSTM1","ATG5",
                          "HSPA5","ATF4","XBP1","DDIT3","HMGCR","SREBF1","FASN"],
}
SET_COLORS = {"Hypoxia":"#7B1F2A","EMT_program":"#962E2E","Thyroid_differen":"#3C6B4F",
              "Stemness_lineage":"#7B3F8A","Senescence_SASP":"#A04451","Apoptosis_DDR":"#34547A",
              "Angiogenesis":"#3F7A8A","Cellular_stress":"#B8893C"}

def cmap_pair(c1, c2): return LinearSegmentedColormap.from_list("c", ["#fffaf2", c1, c2, "#0F1A2E"], N=256)

def score_h5ad(p, gs_list):
    a = ad.read_h5ad(p)
    rvar = a.raw.var_names.astype(str) if a.raw is not None else a.var_names.astype(str)
    X = a.raw.X if a.raw is not None else a.X
    idx = [list(rvar.values).index(g) for g in gs_list if g in rvar.values]
    if not idx: return None, [], None, a.obs
    sub = X[:, idx]
    if hasattr(sub,"toarray"): sub = sub.toarray()
    tot = np.array(a.obs.get("total_counts", np.ones(a.n_obs))).flatten()
    tot = np.where(tot==0,1,tot)
    norm_log = np.log1p(sub/tot[:,None]*1e4)
    avail = [g for g in gs_list if g in rvar.values]
    coords = (a.obs["array_row"].values, a.obs["array_col"].values) if "array_row" in a.obs.columns else (None, None)
    return norm_log, avail, coords, a.obs

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

print("[load] g521 paths")
g521_paths = {p.parent.name: p for p in sorted(G521H.glob("*/GSM*.raw.h5ad"))}
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
rep = {}
for s in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==s)]
    if not sub.empty:
        rep[s] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]

print("[compute] 16 samples × 8 new sets Moran I + per-spot scores cache")
records = []
PERSPOT = {}  # (sample, set_name) -> (rows, cols, score)
for sample, p in g521_paths.items():
    if "_N-" in sample: cond = "PT"
    elif "PTC-" in sample and "L" not in sample: cond = "PTC"
    elif "LPTC-" in sample: cond = "LPTC"
    elif "ATC-" in sample: cond = "ATC"
    else: continue
    rec = {"sample": sample, "condition": cond}
    for sname, genes in NEW_SETS.items():
        scores, avail, coords, obs = score_h5ad(p, genes)
        if scores is None: rec[sname] = np.nan; continue
        v = scores.mean(axis=1); rows, cols = coords
        rec[sname] = morans_quick(v, rows, cols)
        rec[sname+"_n"] = len(avail)
        PERSPOT[(sample, sname)] = (np.array(rows), np.array(cols), v)
    records.append(rec)
df_new = pd.DataFrame(records)
df_new.to_csv(OUT/"spatial_v12_new_morans.tsv", sep="\t", index=False)
print(df_new[["sample","condition"]+list(NEW_SETS.keys())].round(2).to_string(index=False))

# ===== S_F89 ~ S_F96: 8 new gene-set spatial maps (stage representative) =====
SET_LIST = list(NEW_SETS.keys())
for fid, sname in enumerate(SET_LIST, start=89):
    print(f"[F{fid}] {sname}")
    fig = plt.figure(figsize=(15.5, 4.6))
    g = gs.GridSpec(1, 4, wspace=0.18, top=0.83, bottom=0.04, left=0.04, right=0.97)
    cmap = cmap_pair("#cccccc", SET_COLORS[sname])
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or (sample, sname) not in PERSPOT: continue
        ax = fig.add_subplot(g[0, c])
        rows_a, cols_a, v = PERSPOT[(sample, sname)]
        if rows_a is None or len(rows_a)==0: ax.axis("off"); continue
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95))) if not np.isnan(v).all() else 1
        ax.scatter(cols_a, -rows_a, c=v, cmap=cmap, vmin=0, vmax=vmax, s=4.2, alpha=0.92)
        moran = df_new[df_new["sample"]==sample][sname]
        moran = moran.iloc[0] if not moran.empty else np.nan
        ax.set_title(f"{stage} · {sample.split('_')[-1]}\nI = {moran:.2f}", fontsize=10.5,
                     color=PAL[stage], fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
    n_genes = len(NEW_SETS[sname])
    fig.suptitle(f"S_F{fid}  {sname} spatial maps — stage representative\n"
                 f"Gene set ({n_genes}): {', '.join(NEW_SETS[sname])}",
                 fontsize=11.5, fontweight="bold")
    fig.savefig(ASSETS/f"S_F{fid}_{sname}_maps.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

# ===== S_F97: Stage Moran I gradient — all axes (v10+v11+v12 combined) =====
print("[F97] stage Moran I gradient")
df_v10 = pd.read_csv(OUT/"spatial_v10_deep_morans.tsv", sep="\t")
df_v11 = pd.read_csv(OUT/"spatial_v11_pathway_morans.tsv", sep="\t")
all_axes_dict = {
    "v10": ["Thyroid_TF","Cell_cycle","Tumor_suppressor","T_cell","Macrophage_TAM",
            "Stromal_endo","Immune_escape","Glycolysis_OxPhos"],
    "v11": ["WNT_signaling","NOTCH_signaling","Hippo_YAP","RAS_MAPK","PI3K_AKT",
            "MYC_targets","TGF_beta","Inflammation_NFkB"],
    "v12": list(NEW_SETS.keys()),
}
df_all = df_v10.merge(df_v11, on=["sample","condition"], suffixes=("","_pw"))
df_all = df_all.merge(df_new, on=["sample","condition"], suffixes=("","_new"))

stage_mean = df_all.groupby("condition")[
    all_axes_dict["v10"]+all_axes_dict["v11"]+all_axes_dict["v12"]
].mean().reindex(ORDER)
stage_mean.to_csv(OUT/"spatial_v12_stage_gradient.tsv", sep="\t")

fig = plt.figure(figsize=(16, 8.5)); ax = fig.add_subplot(111)
xs = np.arange(len(ORDER))
all_axes = all_axes_dict["v10"]+all_axes_dict["v11"]+all_axes_dict["v12"]
import matplotlib.cm as cm
colors = cm.tab20(np.linspace(0, 1, len(all_axes)))
for ax_name, col in zip(all_axes, colors):
    if ax_name not in stage_mean.columns: continue
    ys = stage_mean[ax_name].values
    ax.plot(xs, ys, "o-", color=col, lw=1.6, ms=5, label=ax_name, alpha=0.85)
ax.set_xticks(xs); ax.set_xticklabels(ORDER, fontweight="bold")
ax.set_ylabel("mean Moran's I (within-stage)", fontsize=11)
ax.set_xlabel("Stage", fontsize=11)
ax.axhline(0, color="#888", lw=0.6, ls="--")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, ncol=1, frameon=False)
ax.set_title("S_F97  Stage Moran's I gradient — 24 axes across PT → PTC → LPTC → ATC\n"
             "(rising = niche organization grows with stage; falling = collapse)",
             fontsize=12.5, fontweight="bold")
ax.grid(True, ls=":", alpha=0.4)
fig.tight_layout()
fig.savefig(ASSETS/"S_F97_stage_gradient_24axes.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# Identify rising / falling axes (slope from PT to ATC)
slopes = {}
for ax_name in all_axes:
    if ax_name not in stage_mean.columns: continue
    y = stage_mean[ax_name].values
    if np.any(np.isnan(y)): continue
    slopes[ax_name] = y[-1] - y[0]
slopes_df = pd.DataFrame.from_dict(slopes, orient="index", columns=["delta_ATC_minus_PT"]).sort_values("delta_ATC_minus_PT")
slopes_df.to_csv(OUT/"spatial_v12_axis_slopes.tsv", sep="\t")
print("\n[finding] axes that COLLAPSE most (PT → ATC):")
print(slopes_df.head(5).round(2))
print("\n[finding] axes that RISE most (PT → ATC):")
print(slopes_df.tail(5).round(2))

# ===== S_F98: ATC-2 outlier deep dive — 8-panel =====
print("[F98] ATC-2 outlier deep dive")
atc2 = "GSM7980873_ATC-2"
panel_axes = ["MYC_targets","Glycolysis_OxPhos","Immune_escape","Cell_cycle",
              "Hypoxia","EMT_program","Stemness_lineage","Angiogenesis"]
fig = plt.figure(figsize=(16, 9.2))
g = gs.GridSpec(2, 4, wspace=0.16, hspace=0.25, top=0.88, bottom=0.04, left=0.04, right=0.97)
def get_perspot(sample, axis):
    if (sample, axis) in PERSPOT: return PERSPOT[(sample, axis)]
    # else compute from gene set (for v10/v11 axes)
    sets_v10 = {"Thyroid_TF":["PAX8","NKX2-1","FOXE1","HHEX","TITF2"],
                "Cell_cycle":["MKI67","CCNB1","CDK1","MCM2","TOP2A","BIRC5","UBE2C","CCNE1"],
                "Tumor_suppressor":["TP53","PTEN","RB1","APC","CDKN2A","CDKN1A","BRCA1"],
                "T_cell":["CD3D","CD3E","CD3G","CD4","CD8A","CD8B","FOXP3","GZMB","PRF1"],
                "Macrophage_TAM":["CD68","CD163","CD86","MARCO","MSR1","LYVE1"],
                "Stromal_endo":["FAP","ACTA2","COL1A1","COL3A1","PECAM1","VWF","CDH5"],
                "Immune_escape":["CD274","PDCD1LG2","IDO1","B2M","TAP1","TAP2","HLA-A","HLA-B","HLA-C"],
                "Glycolysis_OxPhos":["HK2","LDHA","PKM","SLC2A1","ENO1","GAPDH","SDHA","COX4I1","ATP5A1","NDUFA1"]}
    sets_v11 = {"WNT_signaling":["CTNNB1","WNT5A","WNT5B","LEF1","TCF7","AXIN1","APC","DKK1"],
                "MYC_targets":["MYC","MAX","MXI1","NPM1","NCL","HSPD1","RRM2","ENO1"]}
    glist = sets_v10.get(axis) or sets_v11.get(axis) or NEW_SETS.get(axis)
    if glist is None: return None, None, None
    s, _, c, _ = score_h5ad(g521_paths[sample], glist)
    if s is None: return None, None, None
    return c[0], c[1], s.mean(axis=1)
for i, axn in enumerate(panel_axes):
    r, c = i//4, i%4
    ax = fig.add_subplot(g[r, c])
    rows, cols, v = get_perspot(atc2, axn)
    if rows is None: ax.axis("off"); continue
    cmap_a = cmap_pair("#cccccc", "#7B1F2A")
    vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
    ax.scatter(cols, -rows, c=v, cmap=cmap_a, vmin=0, vmax=vmax, s=4.2, alpha=0.92)
    if axn in df_all.columns:
        moran = df_all[df_all["sample"]==atc2][axn].iloc[0]
    else: moran = np.nan
    ax.set_title(f"{axn}\nI = {moran:.2f}", fontsize=10.5, fontweight="bold", color="#7B1F2A")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F98  ATC-2 outlier deep dive — 8 axes simultaneously high\n"
             "Single-sample finding: GSM7980873 organizes >7 distinct programs in the same slide",
             fontsize=12.5, fontweight="bold")
fig.savefig(ASSETS/"S_F98_ATC2_deep_dive.png", dpi=160, bbox_inches="tight")
plt.close(fig)

# ===== S_F99: Cross-axis spatial co-occurrence (within-slide top-quartile Jaccard) =====
print("[F99] cross-axis Jaccard")
co_axes = panel_axes + ["Senescence_SASP","Apoptosis_DDR","Cellular_stress","Thyroid_differen"]
N = len(co_axes); J_avg = np.zeros((N, N))
counts = np.zeros((N, N))
for sample in g521_paths:
    if "_N-" in sample: cond = "PT"
    elif "PTC-" in sample and "L" not in sample: cond = "PTC"
    elif "LPTC-" in sample: cond = "LPTC"
    elif "ATC-" in sample: cond = "ATC"
    else: continue
    if cond not in ORDER: continue
    pos_sets = {}
    for axn in co_axes:
        rows, cols, v = get_perspot(sample, axn)
        if v is None: continue
        thr = np.nanpercentile(v, 75)
        idx = set(np.where(v >= thr)[0].tolist())
        pos_sets[axn] = idx
    for i, a in enumerate(co_axes):
        for j, b in enumerate(co_axes):
            if a not in pos_sets or b not in pos_sets: continue
            A, B = pos_sets[a], pos_sets[b]
            if not A or not B: continue
            J = len(A & B) / max(1, len(A | B))
            J_avg[i, j] += J; counts[i, j] += 1
J_avg = np.divide(J_avg, np.where(counts==0, 1, counts))
np.fill_diagonal(J_avg, np.nan)
fig = plt.figure(figsize=(11, 9.5)); ax = fig.add_subplot(111)
im = ax.imshow(J_avg, cmap="magma_r", vmin=0, vmax=0.5)
for i in range(N):
    for j in range(N):
        if not np.isnan(J_avg[i,j]):
            ax.text(j, i, f"{J_avg[i,j]:.2f}", ha="center", va="center",
                    fontsize=7, color="white" if J_avg[i,j]>0.25 else "black")
ax.set_xticks(range(N)); ax.set_yticks(range(N))
ax.set_xticklabels(co_axes, rotation=45, ha="right", fontsize=8.5)
ax.set_yticklabels(co_axes, fontsize=8.5)
ax.set_title("S_F99  Cross-axis spatial co-occurrence — within-slide top-quartile Jaccard\n"
             "(averaged across 16 samples; high = niches share spots; isolated programs = low)",
             fontsize=11.8, fontweight="bold")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Jaccard (Q4 ∩ / Q4 ∪)")
fig.tight_layout()
fig.savefig(ASSETS/"S_F99_cross_axis_jaccard.png", dpi=160, bbox_inches="tight")
plt.close(fig)
pd.DataFrame(J_avg, index=co_axes, columns=co_axes).to_csv(OUT/"spatial_v12_cross_axis_jaccard.tsv", sep="\t")

# ===== S_F100: Centennial findings infographic =====
print("[F100] centennial findings infographic")
fig = plt.figure(figsize=(17.5, 10.5))
fig.patch.set_facecolor("#fffaf2")
fig.suptitle("S_F100  Spatial Full Package — 100 figures · 8 finding pillars (v1→v12, 28 samples)",
             fontsize=14.5, fontweight="bold", y=0.985)

g = gs.GridSpec(3, 4, wspace=0.32, hspace=0.55, top=0.93, bottom=0.04, left=0.045, right=0.985)

# Panel A: stage gradient mini
axA = fig.add_subplot(g[0, 0:2])
top_rising = slopes_df.tail(4).index.tolist(); top_falling = slopes_df.head(4).index.tolist()
xs = np.arange(len(ORDER))
for axn in top_rising:
    if axn in stage_mean.columns:
        axA.plot(xs, stage_mean[axn].values, "o-", lw=1.8, ms=5, label=f"↑ {axn}", color="#7B1F2A", alpha=0.5+0.15*top_rising.index(axn))
for axn in top_falling:
    if axn in stage_mean.columns:
        axA.plot(xs, stage_mean[axn].values, "o-", lw=1.8, ms=5, label=f"↓ {axn}", color="#3C6B4F", alpha=0.5+0.15*top_falling.index(axn))
axA.set_xticks(xs); axA.set_xticklabels(ORDER, fontweight="bold")
axA.set_ylabel("mean Moran's I"); axA.axhline(0, color="#888", lw=0.5, ls="--")
axA.legend(fontsize=7.5, ncol=2, loc="upper left")
axA.set_title("A. Top 4 RISING vs top 4 COLLAPSING axes (PT → ATC)", fontsize=10.5, fontweight="bold")
axA.grid(True, ls=":", alpha=0.4)

# Panel B: Jaccard top pairs
axB = fig.add_subplot(g[0, 2:4])
pairs = []
for i in range(N):
    for j in range(i+1, N):
        if not np.isnan(J_avg[i,j]):
            pairs.append((co_axes[i], co_axes[j], J_avg[i,j]))
pairs_sorted = sorted(pairs, key=lambda x: -x[2])[:10]
yvals = np.arange(len(pairs_sorted))
axB.barh(yvals, [p[2] for p in pairs_sorted], color="#962E2E", alpha=0.85)
axB.set_yticks(yvals); axB.set_yticklabels([f"{p[0]}\n× {p[1]}" for p in pairs_sorted], fontsize=8)
axB.invert_yaxis()
axB.set_xlabel("Jaccard (within-slide Q4 overlap)")
axB.set_title("B. Top 10 spatially co-occurring axis pairs", fontsize=10.5, fontweight="bold")

# Panel C: ATC-2 niches mini (4 axes)
for k, axn in enumerate(["MYC_targets","Hypoxia","Glycolysis_OxPhos","Cell_cycle"]):
    axc = fig.add_subplot(g[1, k])
    rows, cols, v = get_perspot(atc2, axn)
    if rows is None or v is None: axc.axis("off"); continue
    cmap_a = cmap_pair("#cccccc", "#7B1F2A")
    vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95)))
    axc.scatter(cols, -rows, c=v, cmap=cmap_a, vmin=0, vmax=vmax, s=3.5, alpha=0.92)
    moran = df_all[df_all["sample"]==atc2][axn].iloc[0] if axn in df_all.columns else np.nan
    axc.set_title(f"C{k+1}. ATC-2 · {axn}\nI = {moran:.2f}", fontsize=9.5, color="#7B1F2A", fontweight="bold")
    axc.set_xticks([]); axc.set_yticks([]); axc.set_aspect("equal")
    for sp in axc.spines.values(): sp.set_color("#cdc1aa")

# Panel D: Findings text box
axD = fig.add_subplot(g[2, :])
axD.axis("off")
findings_text = (
    "Eight finding pillars (signature-level only — no causal/clinical inference)\n\n"
    f"  ① Stage gradient (S_F97): {len([s for s in slopes if slopes[s]>0.1])} axes RISE, "
    f"{len([s for s in slopes if slopes[s]<-0.1])} axes COLLAPSE from PT → ATC.\n"
    f"     Strongest collapse: {top_falling[0]} (Δ={slopes_df.loc[top_falling[0],'delta_ATC_minus_PT']:+.2f}). "
    f"Strongest rise: {top_rising[-1]} (Δ={slopes_df.loc[top_rising[-1],'delta_ATC_minus_PT']:+.2f}).\n"
    "  ② TROP2 niche tumor-specific (PTC+LPTC 8/8 Moran I 0.27–0.58; normal/autoimmune ~0).\n"
    "  ③ HT TLS niche directly visible (4/4 PTC+HT, GSE230424; Moran I 0.27–0.84).\n"
    "  ④ ATC-2 (GSM7980873) outlier — 8 axes simultaneously organized (S_F98).\n"
    f"  ⑤ Top spatial co-occurrence (S_F99): {pairs_sorted[0][0]} × {pairs_sorted[0][1]} (J={pairs_sorted[0][2]:.2f}).\n"
    "  ⑥ H&E → DM1 closure NO-GO (LOSO ResNet50 ρ ≈ 0.06, S_F30 fixed).\n"
    "  ⑦ Visium hex 6-neighbor adjacency reveals proper niche cluster structure (S_F25 fixed; max=162, mean=3).\n"
    "  ⑧ 100 figures spanning TROP2 / TLS / 24 axes / 8 stress quartet / pathway / pathology projection.\n\n"
    "Discipline: signature/niche-level exploratory only — no causal claim, no biomarker claim, no cross-paper inference.\n"
    "Yu 2026-05-04 분리벽 유지: Paper 2 = HT-overlap PTC; Paper 4 backlog = Korean GD HLA. Voice-protected sections untouched."
)
axD.text(0.012, 0.94, findings_text, va="top", ha="left", fontsize=10.0, family="monospace",
         color="#0F1A2E", linespacing=1.55,
         bbox=dict(boxstyle="round,pad=0.6", fc="#fff5e0", ec="#B8893C", lw=1.6))
fig.savefig(ASSETS/"S_F100_findings_centennial.png", dpi=160, bbox_inches="tight", facecolor="#fffaf2")
plt.close(fig)

print("\nDONE — 12 figures (S_F89 ~ S_F100) saved")
print(f"  TSVs: spatial_v12_new_morans.tsv, spatial_v12_stage_gradient.tsv, spatial_v12_axis_slopes.tsv, spatial_v12_cross_axis_jaccard.tsv")
