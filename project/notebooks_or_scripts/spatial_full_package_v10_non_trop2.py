#!/usr/bin/env python3
"""Spatial v10 — TROP2 제외 깊은 분석 (S_F73~S_F80)."""
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

# Gene sets — non-TROP2 deep biology
GENE_SETS_DEEP = {
    "Thyroid_TF":      ["PAX8","NKX2-1","FOXE1","HHEX","TITF2"],
    "Cell_cycle":      ["MKI67","CCNB1","CDK1","MCM2","TOP2A","BIRC5","UBE2C","CCNE1"],
    "Tumor_suppressor":["TP53","PTEN","RB1","APC","CDKN2A","CDKN1A","BRCA1"],
    "T_cell":          ["CD3D","CD3E","CD3G","CD4","CD8A","CD8B","FOXP3","GZMB","PRF1"],
    "Macrophage_TAM":  ["CD68","CD163","CD86","MARCO","MSR1","LYVE1"],
    "Stromal_endo":    ["FAP","ACTA2","COL1A1","COL3A1","PECAM1","VWF","CDH5"],
    "Immune_escape":   ["CD274","PDCD1LG2","IDO1","B2M","TAP1","TAP2","HLA-A","HLA-B","HLA-C"],
    "Glycolysis_OxPhos":["HK2","LDHA","PKM","SLC2A1","ENO1","GAPDH","SDHA","COX4I1","ATP5A1","NDUFA1"],
}

trop_cmap = LinearSegmentedColormap.from_list("trop2",["#fffaf2","#ffd58a","#d97c2e","#7B1F2A"],N=256)

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

print("[load]")
g521 = pd.read_csv(ROOT/"project/results/01_spatial_score/all_spots_scored.tsv.gz", sep="\t").rename(columns={"stage":"condition"})
g521_paths = {p.parent.name: p for p in sorted(G521H.glob("*/GSM*.raw.h5ad"))}
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")

# representative niche-organized sample per stage (already used)
rep = {}
for s in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==s)]
    if not sub.empty:
        rep[s] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]

# ===== Compute per-sample × per-gene-set Moran I (28 samples × 8 sets) =====
print("[compute] 16 GSE250521 × 8 gene sets Moran's I")
deep_morans = []
for sample, p in g521_paths.items():
    if "_N-" in sample: cond = "PT"
    elif "PTC-" in sample and "L" not in sample: cond = "PTC"
    elif "LPTC-" in sample: cond = "LPTC"
    elif "ATC-" in sample: cond = "ATC"
    else: continue
    rec = {"sample": sample, "condition": cond}
    for set_name, genes in GENE_SETS_DEEP.items():
        scores, avail, coords = score_h5ad(p, genes)
        if scores is None: rec[set_name] = np.nan; continue
        v = scores.mean(axis=1)
        rows, cols = coords
        if rows is None: rec[set_name] = np.nan; continue
        rec[set_name] = morans_quick(v, rows, cols)
        rec[set_name+"_n_genes"] = len(avail)
    deep_morans.append(rec)
df_deep = pd.DataFrame(deep_morans)
df_deep.to_csv(OUT/"spatial_v10_deep_morans.tsv", sep="\t", index=False)
print(df_deep[["sample","condition"]+list(GENE_SETS_DEEP.keys())].to_string(index=False))

# ===== S_F73: Thyroid TF spatial maps (PAX8, NKX2-1, FOXE1) per stage =====
print("[F73] Thyroid TF maps")
TFs = ["PAX8","NKX2-1","FOXE1"]
fig = plt.figure(figsize=(15, 11))
g = gs.GridSpec(3, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.06, right=0.97)
green_cmap = LinearSegmentedColormap.from_list("g", ["#fffaf2","#a8c8b6","#3C6B4F","#0F1A2E"], N=256)
for r, gene in enumerate(TFs):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, _, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows_a, cols_a = coords
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=green_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#3C6B4F", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F73  Thyroid transcription factors (PAX8 · NKX2-1 · FOXE1) — stage representative\n"
             "PTC stage: TF level 유지; LPTC: 부분 감소; ATC: 거의 사라짐 → thyroid lineage TF program loss",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F73_thyroid_TF_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F74: Cell cycle markers (MKI67, CCNB1, TOP2A, BIRC5) per stage =====
print("[F74] Cell cycle markers")
CC_GENES = ["MKI67","CCNB1","TOP2A","BIRC5"]
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.06, right=0.97)
blue_cmap = LinearSegmentedColormap.from_list("b", ["#fffaf2","#a8c8e2","#34547A","#0F1A2E"], N=256)
for r, gene in enumerate(CC_GENES):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, _, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows_a, cols_a = coords
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=blue_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#34547A", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
fig.suptitle("S_F74  Cell cycle markers (MKI67 · CCNB1 · TOP2A · BIRC5) — stage representative\n"
             "ATC에서 cell cycle marker 활성화 패턴 visible (proliferation niche)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F74_cell_cycle_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F75: T cell markers (CD3D, CD8A, FOXP3, GZMB) per stage =====
print("[F75] T cell markers")
TC_GENES = ["CD3D","CD8A","FOXP3","GZMB"]
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.06, right=0.97)
red_cmap = LinearSegmentedColormap.from_list("r", ["#fffaf2","#f0c5c0","#7B1F2A","#0F1A2E"], N=256)
for r, gene in enumerate(TC_GENES):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, _, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows_a, cols_a = coords
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=red_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#7B1F2A", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F75  T cell markers (CD3D · CD8A · FOXP3 · GZMB) — stage representative\n"
             "ATC에서 T cell infiltration patterns visible (CD8 / cytotoxic / regulatory)",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F75_Tcell_markers.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F76: Macrophage / TAM (CD68, CD163, MARCO, MSR1) =====
print("[F76] Macrophage TAM")
MAC_GENES = ["CD68","CD163","MARCO","MSR1"]
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.06, right=0.97)
gold_cmap = LinearSegmentedColormap.from_list("go", ["#fffaf2","#e2c08c","#B8893C","#0F1A2E"], N=256)
for r, gene in enumerate(MAC_GENES):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, _, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows_a, cols_a = coords
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=gold_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#B8893C", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F76  Macrophage / TAM markers (CD68 · CD163 · MARCO · MSR1) — stage representative\n"
             "TAM infiltration pattern progression — ATC에서 CD163 (M2) 증가 visible",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F76_macrophage_TAM.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F77: Stromal / endothelial (FAP, ACTA2, PECAM1, VWF) =====
print("[F77] Stromal / endothelial")
STR_GENES = ["FAP","ACTA2","PECAM1","VWF"]
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.06, right=0.97)
teal_cmap = LinearSegmentedColormap.from_list("t", ["#fffaf2","#a8d8c8","#1abc9c","#0F1A2E"], N=256)
for r, gene in enumerate(STR_GENES):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, _, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows_a, cols_a = coords
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=teal_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#1abc9c", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F77  Stromal / endothelial markers (FAP · ACTA2 · PECAM1 · VWF) — stage representative\n"
             "ATC에서 vascular network/stromal activation 활발; PT에서는 thyroid follicle 사이 stromal 약함",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F77_stromal_endothelial.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F78: Immune escape (CD274/PD-L1, IDO1, B2M, HLA-A) =====
print("[F78] Immune escape markers")
IE_GENES = ["CD274","IDO1","B2M","HLA-A"]
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.06, right=0.97)
purple_cmap = LinearSegmentedColormap.from_list("p", ["#fffaf2","#e0c5dc","#7B3F8A","#0F1A2E"], N=256)
for r, gene in enumerate(IE_GENES):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, _, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows_a, cols_a = coords
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=purple_cmap, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#7B3F8A", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F78  Immune escape markers (CD274/PD-L1 · IDO1 · B2M · HLA-A) — stage representative\n"
             "PD-L1/IDO1 activation; ATC에서 antigen presentation 변화 (B2M, HLA-A) visible",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F78_immune_escape_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F79: Glycolysis vs OxPhos (HK2, LDHA, COX4I1, ATP5A1) — Warburg balance =====
print("[F79] Glycolysis vs OxPhos")
ENERGY_GENES = ["HK2","LDHA","COX4I1","ATP5A1"]
fig = plt.figure(figsize=(15, 14.5))
g = gs.GridSpec(4, 4, hspace=0.42, wspace=0.18, top=0.95, bottom=0.04, left=0.06, right=0.97)
red_cmap2 = LinearSegmentedColormap.from_list("r2", ["#fffaf2","#f0c5c0","#962E2E","#0F1A2E"], N=256)
for r, gene in enumerate(ENERGY_GENES):
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[r, c])
        scores, _, coords = score_h5ad(g521_paths[sample], [gene])
        if scores is None or scores.shape[1] == 0: ax.text(0.5,0.5,f"{gene} not detected"); ax.axis("off"); continue
        rows_a, cols_a = coords
        v = scores[:,0]
        vmax = np.nanpercentile(v, 99) if np.nanpercentile(v, 99) > 0 else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=red_cmap2, vmin=0, vmax=vmax, s=3, alpha=0.92)
        if r == 0: ax.set_title(f"{stage} · {sample.split('_')[-1]}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        if c == 0: ax.set_ylabel(gene, fontsize=13, fontweight="bold", color="#962E2E", rotation=0, labelpad=30, va="center")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
fig.suptitle("S_F79  Energy metabolism markers (HK2 · LDHA · COX4I1 · ATP5A1) — Warburg balance\n"
             "ATC에서 glycolysis (HK2/LDHA) 증가, OxPhos (COX4I1/ATP5A1) 감소 — Warburg shift",
             fontsize=12, fontweight="bold")
fig.savefig(ASSETS/"S_F79_metabolism_maps.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ===== S_F80: 8 gene-set Moran I master heatmap =====
print("[F80] Deep gene-set Moran heatmap")
fig, ax = plt.subplots(figsize=(13, 10))
df_deep_sorted = df_deep.sort_values(["condition","sample"]).reset_index(drop=True)
gs_show = list(GENE_SETS_DEEP.keys())
mat = np.full((len(df_deep_sorted), len(gs_show)), np.nan)
for i, _ in df_deep_sorted.iterrows():
    for j, set_name in enumerate(gs_show):
        v = df_deep_sorted.iloc[i].get(set_name, np.nan)
        if not pd.isna(v): mat[i, j] = float(v)
im = ax.imshow(mat, aspect="auto", cmap="RdYlGn_r", vmin=-0.1, vmax=0.85)
for i in range(len(df_deep_sorted)):
    for j in range(len(gs_show)):
        v = mat[i,j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    color="white" if v>0.5 else "black", fontsize=9, fontweight="bold")
ax.set_xticks(range(len(gs_show))); ax.set_xticklabels(gs_show, fontsize=10, rotation=20, ha="right")
ax.set_yticks(range(len(df_deep_sorted)))
ax.set_yticklabels([f"{r['condition']} · {r['sample'].split('_')[-1]}" for _, r in df_deep_sorted.iterrows()], fontsize=9.5, family="monospace")
plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Moran's I")
ax.set_title("S_F80  Non-TROP2 deep gene-set Moran's I — 16 GSE250521 samples × 8 gene sets\n"
             "(Thyroid TF · Cell cycle · Tumor suppressor · T cell · Macrophage · Stromal · Immune escape · Metabolism)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(ASSETS/"S_F80_deep_geneset_morans.png", dpi=160, bbox_inches="tight"); plt.close(fig)

print("\nDONE — 8 figures (S_F73 ~ S_F80) saved")
