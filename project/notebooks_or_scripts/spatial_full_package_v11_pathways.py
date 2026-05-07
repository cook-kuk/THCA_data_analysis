#!/usr/bin/env python3
"""Spatial v11 — 8 additional non-TROP2 pathway analyses (S_F81~S_F88)."""
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

# Cancer pathway gene sets — non-TROP2
PATHWAYS = {
    "WNT_signaling":     ["CTNNB1","WNT5A","WNT5B","LEF1","TCF7","AXIN1","APC","DKK1"],
    "NOTCH_signaling":   ["NOTCH1","NOTCH2","NOTCH3","JAG1","JAG2","DLL1","HES1","HEY1"],
    "Hippo_YAP":         ["YAP1","WWTR1","TEAD1","TEAD2","CTGF","CYR61","AMOTL1","LATS1"],
    "RAS_MAPK":          ["KRAS","HRAS","NRAS","BRAF","MAP2K1","MAPK1","MAPK3","DUSP6","SPRY2"],
    "PI3K_AKT":          ["PIK3CA","AKT1","AKT2","AKT3","PTEN","MTOR","RPS6","GSK3B"],
    "MYC_targets":       ["MYC","MAX","MXI1","NPM1","NCL","HSPD1","RRM2","ENO1"],
    "TGF_beta":          ["TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SMAD7"],
    "Inflammation_NFkB": ["NFKB1","NFKB2","RELA","RELB","TNF","IL6","IL1B","CXCL8","NFKBIA"],
}

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

print("[load] g521 paths")
g521_paths = {p.parent.name: p for p in sorted(G521H.glob("*/GSM*.raw.h5ad"))}
df_E = pd.read_csv(OUT/"spatial_E_TROP2_spot_distribution.tsv", sep="\t")
rep = {}
for s in ORDER:
    sub = df_E[(df_E.dataset=="GSE250521") & (df_E.condition==s)]
    if not sub.empty:
        rep[s] = sub.sort_values("morans_I_TROP2", ascending=False).iloc[0]["sample_id"]

print("[compute] 16 samples × 8 pathways Moran I")
pathway_morans = []
pathway_scores = {}
for sample, p in g521_paths.items():
    if "_N-" in sample: cond = "PT"
    elif "PTC-" in sample and "L" not in sample: cond = "PTC"
    elif "LPTC-" in sample: cond = "LPTC"
    elif "ATC-" in sample: cond = "ATC"
    else: continue
    rec = {"sample": sample, "condition": cond}
    for pw_name, genes in PATHWAYS.items():
        scores, avail, coords = score_h5ad(p, genes)
        if scores is None: rec[pw_name] = np.nan; continue
        v = scores.mean(axis=1)
        rows, cols = coords
        rec[pw_name] = morans_quick(v, rows, cols)
        rec[pw_name+"_n"] = len(avail)
    pathway_morans.append(rec)
df_pw = pd.DataFrame(pathway_morans)
df_pw.to_csv(OUT/"spatial_v11_pathway_morans.tsv", sep="\t", index=False)
print(df_pw[["sample","condition"]+list(PATHWAYS.keys())].to_string(index=False))

# ===== Build cmaps for each pathway =====
def cmap_pair(c1, c2): return LinearSegmentedColormap.from_list("c", ["#fffaf2", c1, c2, "#0F1A2E"], N=256)
PW_COLORS = {
    "WNT_signaling":     "#962E2E",
    "NOTCH_signaling":   "#3F7A8A",
    "Hippo_YAP":         "#34547A",
    "RAS_MAPK":          "#B8893C",
    "PI3K_AKT":          "#7B3F8A",
    "MYC_targets":       "#7B1F2A",
    "TGF_beta":          "#3C6B4F",
    "Inflammation_NFkB": "#A04451",
}

# ===== S_F81~S_F88: 8 pathway × stage representative spatial maps =====
PW_LIST = list(PATHWAYS.keys())
for fid, pw_name in enumerate(PW_LIST, start=81):
    print(f"[F{fid}] {pw_name}")
    fig = plt.figure(figsize=(15, 4.4))
    g = gs.GridSpec(1, 4, wspace=0.18, top=0.85, bottom=0.04, left=0.04, right=0.97)
    cmap_pw = cmap_pair("#cccccc", PW_COLORS[pw_name])
    for c, stage in enumerate(ORDER):
        sample = rep.get(stage)
        if not sample or sample not in g521_paths: continue
        ax = fig.add_subplot(g[0, c])
        scores, avail, coords = score_h5ad(g521_paths[sample], PATHWAYS[pw_name])
        if scores is None: ax.axis("off"); continue
        rows_a, cols_a = coords
        if rows_a is None: ax.axis("off"); continue
        v = scores.mean(axis=1)
        vmax = max(abs(np.nanpercentile(v,5)), abs(np.nanpercentile(v,95))) if not np.isnan(v).all() else 1
        ax.scatter(np.array(cols_a), -np.array(rows_a), c=v, cmap=cmap_pw, vmin=0, vmax=vmax, s=4, alpha=0.92)
        moran = df_pw[df_pw["sample"]==sample][pw_name]
        moran = moran.iloc[0] if not moran.empty else np.nan
        ax.set_title(f"{stage} · {sample.split('_')[-1]}\nI = {moran:.2f}", fontsize=10.5, color=PAL[stage], fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_color("#cdc1aa")
    n_genes = df_pw[df_pw["condition"]=="PTC"][pw_name+"_n"].iloc[0] if pw_name+"_n" in df_pw.columns else "?"
    fig.suptitle(f"S_F{fid}  {pw_name} pathway spatial maps — stage representative\n"
                 f"Gene set ({len(PATHWAYS[pw_name])} genes): {', '.join(PATHWAYS[pw_name])}",
                 fontsize=11.5, fontweight="bold")
    fig.savefig(ASSETS/f"S_F{fid}_{pw_name}_maps.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

print("\nDONE — 8 figures (S_F81 ~ S_F88) saved")
