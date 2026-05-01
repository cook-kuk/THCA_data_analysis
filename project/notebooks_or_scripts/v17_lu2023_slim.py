#!/usr/bin/env python3
"""GSE193581 (Lu 2023 JCI) — v17p35 8-gene panel transfer, MEMORY-SLIM version.

Skips full HVG/scale/PCA/UMAP (which densified 67k×38k matrix → swap thrashing).
Computes 8-panel z-score directly on log-normalized sparse expression. Outputs:
- per-cell DM_score TSV
- per-sample summary TSV
- histology × DM-class crosstab
- per-sample box plot (panel_z by histology, scRNA cells)
- summary JSON + markdown
"""
from __future__ import annotations
import os, sys, json, time, logging
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats, sparse
from sklearn.metrics import adjusted_rand_score
import plotly.graph_objects as go

RAW_DIR  = Path("/data/thca/v17_lu2023_GSE193581/raw")
ANN_FILE = Path("/data/thca/v17_lu2023_GSE193581/celltype_annotation.txt.gz")
OUT      = Path("/opt/thyroid-dash/project/results/v17_lu2023");          OUT.mkdir(parents=True, exist_ok=True)
FIG      = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17"); FIG.mkdir(parents=True, exist_ok=True)
RPT      = Path("/opt/thyroid-dash/project/reports/v17p35");              RPT.mkdir(parents=True, exist_ok=True)
LOG      = Path("/opt/thyroid-dash/project/logs/v17_lu2023_slim.log")

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.FileHandler(LOG, mode="w"), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
log.info("=== v17_lu2023 SLIM panel transfer START ===")

BG, INK = "#0b0e12", "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font=dict(color=INK, size=14))
COLOR = {"NORM":"#2ECC71","PTC":"#F5A623","ATC":"#c24c4c"}
SEED = 42; np.random.seed(SEED)

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
log.info(f"panel: {PANEL}")

# ---------- annotation ----------
ann = pd.read_csv(ANN_FILE, sep="\t")
ann.columns = ["sample_id","celltype"]
log.info(f"annotation rows: {len(ann)}, samples in ann: {ann['sample_id'].nunique()}, celltypes: {ann['celltype'].nunique()}")
ann_by_sample = {s: g["celltype"].values for s,g in ann.groupby("sample_id")}

def histo(s):
    if s.startswith("PTC"): return "PTC"
    if s.startswith("ATC"): return "ATC"
    if s.startswith("NORM"): return "NORM"
    return "OTHER"

# Process each sample one at a time, only keeping the 8 panel rows + metadata
# so the merged matrix is 8 × 67k (~2 MB), not 38k × 67k (~10 GB).
samples_data = []   # list of (sample_id, histology, cell_barcodes, panel_X_log_normed [8 × n_cells], celltypes)

t_start = time.time()
for fp in sorted(RAW_DIR.glob("*_UMI.txt.gz")):
    sname = fp.name.split("_")[1]
    ann_key = sname if sname in ann_by_sample else (sname + "T")
    if ann_key not in ann_by_sample:
        log.warning(f"  {sname}: no annotation, SKIP")
        continue
    t0 = time.time()
    df = pd.read_csv(fp, sep="\t", index_col=0)  # genes × cells
    cells = df.columns.tolist()
    expected = ann_by_sample[ann_key]
    n = min(len(cells), len(expected))
    df = df.iloc[:, :n]
    cells = cells[:n]
    ctypes = expected[:n]

    # QC: filter cells with reasonable n_genes & low mito
    cell_n_genes = (df > 0).sum(axis=0).values
    mt_mask = df.index.str.startswith("MT-")
    cell_pct_mt = 100 * df.loc[mt_mask].sum(axis=0).values / np.maximum(df.sum(axis=0).values, 1)
    keep_cells = (cell_n_genes > 200) & (cell_n_genes < 6000) & (cell_pct_mt < 20)
    df = df.loc[:, keep_cells]
    cells = [c for c,k in zip(cells, keep_cells) if k]
    ctypes_kept = [c for c,k in zip(ctypes, keep_cells) if k]
    if df.shape[1] == 0:
        log.warning(f"  {sname}: 0 cells post-QC, skip"); continue

    # Library-size normalize per cell + log1p (in-place, only on dense values for this sample's small matrix)
    libsize = df.sum(axis=0).values
    # CPM-like: counts per 10000
    df_norm = df.div(libsize, axis=1) * 1e4
    df_log = np.log1p(df_norm.values).astype(np.float32)

    # Pull only 8 panel rows
    present_genes = [g for g in PANEL if g in df.index]
    if len(present_genes) < 8:
        log.warning(f"  {sname}: only {len(present_genes)}/8 panel genes — {set(PANEL)-set(present_genes)} missing")
    panel_X = pd.DataFrame(df_log[df.index.get_indexer(present_genes), :],
                           index=present_genes, columns=cells)
    samples_data.append({
        "sample": sname,
        "histology": histo(sname),
        "cells": cells,
        "panel_X_log": panel_X,        # 8 × n_cells, log-normed
        "celltypes": np.array(ctypes_kept),
    })
    log.info(f"  {sname} ({histo(sname)}): {df.shape[1]} cells post-QC ({time.time()-t0:.1f}s)")

log.info(f"loaded {len(samples_data)} samples in {time.time()-t_start:.1f}s")

# Merge all panel X into one wide DataFrame
all_panel = pd.concat([sd["panel_X_log"] for sd in samples_data], axis=1)
all_meta = pd.DataFrame({
    "cell": np.concatenate([sd["cells"] for sd in samples_data]),
    "sample": np.concatenate([[sd["sample"]]*len(sd["cells"]) for sd in samples_data]),
    "histology": np.concatenate([[sd["histology"]]*len(sd["cells"]) for sd in samples_data]),
    "author_celltype": np.concatenate([sd["celltypes"] for sd in samples_data]),
}).set_index("cell")
all_panel = all_panel.loc[:, all_meta.index]
log.info(f"merged: panel matrix {all_panel.shape}; meta {all_meta.shape}")
log.info(f"histology counts (all cells): {all_meta['histology'].value_counts().to_dict()}")
log.info(f"author_celltype top 10: {all_meta['author_celltype'].value_counts().head(10).to_dict()}")

# ---------- filter to malignant/epithelial cells ----------
malignant_ct = ["Malignant cell","Epithelial cell","Thyroid follicular cell","Thyrocyte","Tumor cell","Cancer cell","Malignant"]
mal_mask = all_meta["author_celltype"].isin(malignant_ct)
if mal_mask.sum() == 0:
    # author labels are different — guess by marker presence (TG/TPO high)
    log.warning(f"no cells matched malignant labels; using TG-high cells as proxy")
    tg_expr = all_panel.loc["TG"] if "TG" in all_panel.index else None
    if tg_expr is not None:
        thresh = tg_expr.quantile(0.5)
        mal_mask = (tg_expr >= thresh)
    else:
        mal_mask = pd.Series(True, index=all_meta.index)

mal_meta = all_meta[mal_mask].copy()
mal_panel = all_panel.loc[:, mal_meta.index]
log.info(f"malignant/epithelial cells: {mal_mask.sum()} ({mal_meta['histology'].value_counts().to_dict()})")

# ---------- panel z-score per cell ----------
# z-score per gene across all malignant cells (cohort-level normalization, like bulk panel transfer)
gene_mean = mal_panel.mean(axis=1)
gene_std  = mal_panel.std(axis=1) + 1e-9
z = mal_panel.sub(gene_mean, axis=0).div(gene_std, axis=0)
DM_score = z.mean(axis=0)  # mean z across 8 (or fewer) genes per cell
mal_meta["DM_score"] = DM_score.values
log.info(f"DM_score: mean={mal_meta['DM_score'].mean():.3f} sd={mal_meta['DM_score'].std():.3f} median={mal_meta['DM_score'].median():.3f}")

# DM1 = high DM_score (well-differentiated), DM2 = low (less-differentiated)
med = mal_meta["DM_score"].median()
mal_meta["DM_class"] = np.where(mal_meta["DM_score"] >= med, "DM1_high", "DM2_low")
log.info(f"DM_class (median split @ {med:.3f}): {mal_meta['DM_class'].value_counts().to_dict()}")

# ---------- cross-tab + statistics ----------
ct = pd.crosstab(mal_meta["histology"], mal_meta["DM_class"])
log.info(f"\ncross-tab histology × DM_class:\n{ct}")
chi2, p_chi, dof, _ = stats.chi2_contingency(ct)
log.info(f"chi-square: chi2={chi2:.2f} dof={dof} p={p_chi:.3e}")

# ARI
hist_codes = pd.Categorical(mal_meta["histology"]).codes
dm_codes   = pd.Categorical(mal_meta["DM_class"]).codes
ari = adjusted_rand_score(hist_codes, dm_codes)
log.info(f"ARI(histology, DM_class) = {ari:.4f}")

# Mann-Whitney PTC vs ATC DM_score
dm_ptc = mal_meta.loc[mal_meta["histology"]=="PTC", "DM_score"].values
dm_atc = mal_meta.loc[mal_meta["histology"]=="ATC", "DM_score"].values
dm_norm = mal_meta.loc[mal_meta["histology"]=="NORM", "DM_score"].values
u, p_mw = stats.mannwhitneyu(dm_ptc, dm_atc, alternative="two-sided")
log.info(f"MW PTC vs ATC DM_score: U={u:.0f} p={p_mw:.3e}; PTC median={np.median(dm_ptc):.3f}, ATC median={np.median(dm_atc):.3f}")
# Kruskal-Wallis 3-group
groups = []
hlabels = []
for h in ["NORM","PTC","ATC"]:
    g = mal_meta.loc[mal_meta["histology"]==h, "DM_score"].values
    if len(g) > 0:
        groups.append(g); hlabels.append(h)
H, p_kw = stats.kruskal(*groups)
log.info(f"Kruskal-Wallis 3-group ({hlabels}): H={H:.2f} p={p_kw:.3e}")

# Per-sample
per_sample = mal_meta.groupby("sample").agg(
    histology=("histology","first"),
    n_cells=("DM_score","count"),
    median_DM=("DM_score","median"),
    mean_DM=("DM_score","mean"),
    pct_DM1=("DM_class", lambda s: (s=="DM1_high").mean()),
).sort_values(["histology","median_DM"])
log.info(f"\nper-sample summary:\n{per_sample}")

# ---------- save ----------
mal_meta[["sample","histology","author_celltype","DM_score","DM_class"]].to_csv(
    OUT/"GSE193581_cell_panel_score.tsv", sep="\t")
per_sample.to_csv(OUT/"GSE193581_per_sample_summary.tsv", sep="\t")
ct.to_csv(OUT/"GSE193581_crosstab_histology_DM.tsv", sep="\t")

summary = {
    "cohort": "GSE193581 (Lu et al. 2023 JCI)",
    "study": "Single-cell map of differentiated → anaplastic thyroid transformation",
    "n_samples": int(all_meta["sample"].nunique()),
    "n_cells_total": int(len(all_meta)),
    "n_malignant_cells": int(mal_meta.shape[0]),
    "panel_genes_present": list(mal_panel.index),
    "panel_genes_missing": list(set(PANEL) - set(mal_panel.index)),
    "histology_counts_all": all_meta["histology"].value_counts().to_dict(),
    "histology_counts_malignant": mal_meta["histology"].value_counts().to_dict(),
    "DM_score_median_all": float(med),
    "DM_class_counts": mal_meta["DM_class"].value_counts().to_dict(),
    "crosstab_histology_DM": ct.to_dict(),
    "chi_square_p": float(p_chi),
    "chi_square_stat": float(chi2),
    "ARI_histology_DM": float(ari),
    "MW_PTC_vs_ATC_DM_score_p": float(p_mw),
    "DM_score_PTC_median": float(np.median(dm_ptc)),
    "DM_score_ATC_median": float(np.median(dm_atc)),
    "DM_score_NORM_median": float(np.median(dm_norm)) if len(dm_norm) else None,
    "kruskal_wallis_3grp_p": float(p_kw),
    "kruskal_wallis_3grp_H": float(H),
    "per_sample_summary": per_sample.reset_index().to_dict("records"),
    "seed": SEED,
}
(OUT/"GSE193581_summary.json").write_text(json.dumps(summary, indent=2, default=str))
log.info(f"saved {OUT}/GSE193581_summary.json")

# ---------- figures ----------
# Fig A: per-cell DM_score by histology (box+strip)
fig = go.Figure()
order = ["NORM","PTC","ATC"]
for h in order:
    sub = mal_meta[mal_meta["histology"]==h]
    if len(sub) == 0: continue
    fig.add_trace(go.Box(y=sub["DM_score"], x=[h]*len(sub),
                         name=f"{h} (n={len(sub):,} cells)",
                         marker=dict(color=COLOR.get(h,"#888"), size=2, opacity=0.4),
                         boxpoints="outliers", line=dict(width=2)))
fig.update_layout(title=dict(
    text=f"<b>GSE193581 (Lu 2023 JCI) — single-cell 8-gene panel DM_score by histology</b><br>"
         f"<sub style='color:#7ccfcd'>{int(mal_meta.shape[0]):,} malignant/epithelial cells from {len(samples_data)} samples; "
         f"NORM median={np.median(dm_norm):+.2f}, PTC median={np.median(dm_ptc):+.2f}, ATC median={np.median(dm_atc):+.2f}; "
         f"Kruskal-Wallis p={p_kw:.2e}, MW PTC vs ATC p={p_mw:.2e}</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="DM_score (mean z across 8 panel genes per cell)",
               gridcolor="rgba(255,255,255,0.08)", zeroline=True, zerolinecolor="rgba(255,255,255,0.2)"),
    xaxis=dict(tickfont=dict(size=12, color=INK)),
    showlegend=False, height=520, margin=dict(l=80,r=20,t=120,b=60), **DARK)
fig.write_html(FIG/"v17_lu2023_dmscore_by_histology.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_dmscore_by_histology.html")

# Fig B: per-sample DM_score box, ordered by median, colored by histology
fig = go.Figure()
sample_order = per_sample.sort_values(["histology","median_DM"], ascending=[True, False]).index.tolist()
for s in sample_order:
    sub = mal_meta[mal_meta["sample"]==s]
    h = sub["histology"].iloc[0]
    color = COLOR.get(h, "#888")
    fig.add_trace(go.Box(y=sub["DM_score"], x=[s]*len(sub),
                         name=f"{s} ({h}, n={len(sub):,})",
                         marker=dict(color=color, size=2, opacity=0.5),
                         boxpoints="outliers", line=dict(width=1.5)))
fig.update_layout(title=dict(
    text=f"<b>GSE193581 — per-sample DM_score (malignant cells, 23 samples sorted by median)</b><br>"
         f"<sub style='color:#7ccfcd'>green=NORM ({(per_sample.histology=='NORM').sum()}), "
         f"orange=PTC ({(per_sample.histology=='PTC').sum()}), "
         f"red=ATC ({(per_sample.histology=='ATC').sum()})</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="DM_score", gridcolor="rgba(255,255,255,0.08)",
               zeroline=True, zerolinecolor="rgba(255,255,255,0.2)"),
    xaxis=dict(tickangle=-45, tickfont=dict(size=10)),
    showlegend=False, height=560, width=1100,
    margin=dict(l=80,r=20,t=110,b=140), **DARK)
fig.write_html(FIG/"v17_lu2023_per_sample_box.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_per_sample_box.html")

# Fig C: cross-tab heatmap (% of row, with counts annotated)
ct_norm = ct.div(ct.sum(axis=1), axis=0) * 100
fig = go.Figure(go.Heatmap(z=ct_norm.values, x=ct_norm.columns.tolist(),
    y=ct_norm.index.tolist(), colorscale="RdBu_r", zmid=50, zmin=0, zmax=100,
    text=[[f"<b>{int(ct.loc[r,c]):,}</b><br>({ct_norm.loc[r,c]:.1f}%)"
           for c in ct.columns] for r in ct.index],
    texttemplate="%{text}",
    textfont=dict(size=14, color=INK),
    colorbar=dict(title=dict(text="% of row", font=dict(color=INK)), tickfont=dict(color=INK))))
fig.update_layout(title=dict(
    text=f"<b>GSE193581 — histology × DM-class cross-tab (malignant cells)</b><br>"
         f"<sub style='color:#7ccfcd'>chi²={chi2:.0f}, dof={dof}, p={p_chi:.2e}; ARI={ari:.3f}</sub>",
    font=dict(size=14, color=INK)),
    xaxis=dict(title="DM class (median split)", tickfont=dict(color=INK)),
    yaxis=dict(title="Histology (author)", tickfont=dict(color=INK)),
    height=420, width=720, margin=dict(l=120,r=40,t=110,b=80), **DARK)
fig.write_html(FIG/"v17_lu2023_crosstab.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_crosstab.html")

# Fig D: per-sample DM1 percentage bar
fig = go.Figure()
for s in sample_order:
    h = per_sample.loc[s,"histology"]
    pct = per_sample.loc[s,"pct_DM1"] * 100
    fig.add_trace(go.Bar(x=[s], y=[pct], marker=dict(color=COLOR.get(h,"#888")),
                         name=f"{s} ({h})", showlegend=False,
                         text=[f"{pct:.0f}%"], textposition="outside"))
fig.update_layout(title=dict(
    text=f"<b>GSE193581 — per-sample DM1-high cell %</b><br>"
         f"<sub style='color:#7ccfcd'>NORM samples DM1% expected high → ATC samples DM1% expected low</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="% cells DM1 (high differentiation)", range=[0, 110]),
    xaxis=dict(tickangle=-45),
    height=460, width=1100, margin=dict(l=80,r=20,t=100,b=140), **DARK)
fig.write_html(FIG/"v17_lu2023_dm1_pct_bar.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_dm1_pct_bar.html")

# ============================================================
# UMAP layer — memory-safe (HVG subset BEFORE scale)
# ============================================================
log.info("=== UMAP layer (HVG 2000 subset) ===")
try:
    import scanpy as sc
    sc.settings.verbosity = 1
    # Build full sparse matrix from ALL samples (cells × genes), reuse memory carefully
    # Only need this for UMAP visualization
    log.info("loading full sparse matrices for UMAP (one sample at a time)...")
    adatas = []
    for fp in sorted(RAW_DIR.glob("*_UMI.txt.gz")):
        sname = fp.name.split("_")[1]
        ann_key = sname if sname in ann_by_sample else (sname + "T")
        if ann_key not in ann_by_sample: continue
        df = pd.read_csv(fp, sep="\t", index_col=0)
        cells = df.columns.tolist()
        expected = ann_by_sample[ann_key]
        n = min(len(cells), len(expected))
        df = df.iloc[:, :n]
        cells = cells[:n]
        ctypes = expected[:n]
        # cells × genes sparse
        X = sparse.csr_matrix(df.values.T.astype(np.float32))
        a = sc.AnnData(X=X)
        a.obs_names = cells
        a.var_names = df.index.tolist()
        a.obs["sample"] = sname
        a.obs["histology"] = histo(sname)
        a.obs["author_celltype"] = pd.Categorical(ctypes)
        adatas.append(a)
        del df
    adata = sc.concat(adatas, axis=0, join="outer", merge="same", index_unique=None)
    del adatas
    adata.var_names_make_unique()
    log.info(f"merged adata for UMAP: {adata.shape}")

    # QC
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    mask = (adata.obs["n_genes_by_counts"] > 200) & (adata.obs["n_genes_by_counts"] < 6000) \
           & (adata.obs["pct_counts_mt"] < 20)
    adata = adata[mask].copy()
    log.info(f"post-QC: {adata.shape}")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    log.info(f"HVG flagged: {adata.var.highly_variable.sum()}")

    # CRITICAL: subset to HVG BEFORE scale → keeps dense matrix small
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    log.info(f"HVG subset: {adata_hvg.shape} (dense ~ {adata_hvg.shape[0]*adata_hvg.shape[1]*4/1e9:.2f} GB float32)")
    sc.pp.scale(adata_hvg, max_value=10)
    log.info("scale done")

    sc.tl.pca(adata_hvg, n_comps=50, random_state=SEED)
    log.info("PCA done")
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=50, random_state=SEED)
    log.info("neighbors done")
    sc.tl.umap(adata_hvg, random_state=SEED)
    log.info(f"UMAP done: {adata_hvg.obsm['X_umap'].shape}")

    # propagate DM_score back to UMAP cells
    cell_dm = mal_meta["DM_score"]   # only malignant cells have DM_score
    adata_hvg.obs["DM_score"] = np.nan
    adata_hvg.obs.loc[cell_dm.index.intersection(adata_hvg.obs_names), "DM_score"] = \
        cell_dm.loc[cell_dm.index.intersection(adata_hvg.obs_names)].values

    um = pd.DataFrame(adata_hvg.obsm["X_umap"], index=adata_hvg.obs_names, columns=["UMAP1","UMAP2"])
    um["author_celltype"] = adata_hvg.obs["author_celltype"].astype(str).values
    um["histology"]       = adata_hvg.obs["histology"].values
    um["sample"]          = adata_hvg.obs["sample"].values
    um["DM_score"]        = adata_hvg.obs["DM_score"].values

    # Fig E: UMAP by histology
    fig = go.Figure()
    for h in ["NORM","PTC","ATC"]:
        sub = um[um["histology"]==h]
        if len(sub) == 0: continue
        fig.add_trace(go.Scattergl(x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
            marker=dict(size=2.5, color=COLOR.get(h,"#888"), opacity=0.6),
            name=f"{h} (n={len(sub):,})", hoverinfo="skip"))
    fig.update_layout(title=dict(
        text=f"<b>GSE193581 — UMAP colored by histology (n={adata_hvg.n_obs:,} cells, 23 samples)</b><br>"
             f"<sub style='color:#7ccfcd'>green=NORM (n={int((um['histology']=='NORM').sum()):,}), "
             f"orange=PTC ({int((um['histology']=='PTC').sum()):,}), "
             f"red=ATC ({int((um['histology']=='ATC').sum()):,})</sub>",
        font=dict(size=14, color=INK)),
        xaxis=dict(title="UMAP1", showgrid=False, zeroline=False),
        yaxis=dict(title="UMAP2", showgrid=False, zeroline=False),
        height=620, width=900, margin=dict(l=60,r=20,t=90,b=60), **DARK,
        legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)"))
    fig.write_html(FIG/"v17_lu2023_umap_histology.html", include_plotlyjs="cdn", full_html=True)
    log.info(f"wrote {FIG}/v17_lu2023_umap_histology.html")

    # Fig F: UMAP by author cell type
    ctype_palette = {"T cell":"#7ccfcd","Malignant cell":"#c24c4c","Myeloid cell":"#F5A623",
                     "B cell":"#9b59b6","NK cell":"#1abc9c","Fibroblast":"#e67e22",
                     "Endothelial cell":"#3498db","Epithelial cell":"#F2F2F2",
                     "Stromal cell":"#888","Plasma cell":"#bb9bd6"}
    fig = go.Figure()
    cts_seen = um["author_celltype"].value_counts()
    for ct_label, n_cells in cts_seen.head(12).items():
        sub = um[um["author_celltype"]==ct_label]
        color = ctype_palette.get(ct_label, "#aaa")
        fig.add_trace(go.Scattergl(x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
            marker=dict(size=2, color=color, opacity=0.6),
            name=f"{ct_label} (n={int(n_cells):,})", hoverinfo="skip"))
    fig.update_layout(title=dict(
        text=f"<b>GSE193581 — UMAP colored by author cell-type label</b><br>"
             f"<sub style='color:#7ccfcd'>{int(cts_seen.size)} cell types in dataset; top 12 shown</sub>",
        font=dict(size=14, color=INK)),
        xaxis=dict(title="UMAP1", showgrid=False), yaxis=dict(title="UMAP2", showgrid=False),
        height=620, width=950, margin=dict(l=60,r=20,t=90,b=60), **DARK,
        legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)"))
    fig.write_html(FIG/"v17_lu2023_umap_celltype.html", include_plotlyjs="cdn", full_html=True)
    log.info(f"wrote {FIG}/v17_lu2023_umap_celltype.html")

    # Fig G: UMAP by DM_score (malignant cells overlay, others gray)
    um_mal = um.dropna(subset=["DM_score"])
    um_oth = um[um["DM_score"].isna()]
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=um_oth["UMAP1"], y=um_oth["UMAP2"], mode="markers",
        marker=dict(size=2, color="rgba(120,120,120,0.20)"),
        name=f"non-malignant (n={len(um_oth):,})", hoverinfo="skip"))
    fig.add_trace(go.Scattergl(x=um_mal["UMAP1"], y=um_mal["UMAP2"], mode="markers",
        marker=dict(size=3, color=um_mal["DM_score"], colorscale="RdBu",
                    cmid=0, cmin=-1.5, cmax=1.5, showscale=True,
                    colorbar=dict(title=dict(text="DM_score", font=dict(color=INK)),
                                  tickfont=dict(color=INK))),
        name=f"malignant (n={len(um_mal):,})", hoverinfo="skip"))
    fig.update_layout(title=dict(
        text=f"<b>GSE193581 — UMAP colored by 8-gene panel DM_score (malignant only)</b><br>"
             f"<sub style='color:#7ccfcd'>panel: SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1; "
             f"NORM med={np.median(dm_norm):+.2f}, PTC med={np.median(dm_ptc):+.2f}, ATC med={np.median(dm_atc):+.2f}; "
             f"KW p={p_kw:.2e}</sub>",
        font=dict(size=14, color=INK)),
        xaxis=dict(title="UMAP1", showgrid=False), yaxis=dict(title="UMAP2", showgrid=False),
        height=620, width=950, margin=dict(l=60,r=20,t=90,b=60), **DARK)
    fig.write_html(FIG/"v17_lu2023_umap_dmscore.html", include_plotlyjs="cdn", full_html=True)
    log.info(f"wrote {FIG}/v17_lu2023_umap_dmscore.html")

    # Save adata_hvg to disk for reuse
    adata_hvg.write_h5ad(OUT/"GSE193581_hvg_adata.h5ad")
    log.info(f"wrote {OUT}/GSE193581_hvg_adata.h5ad")

except Exception as e:
    log.exception(f"UMAP layer failed: {e}")
    log.info("Continuing without UMAP — slim panel z-score figures already written.")

# ---------- markdown summary ----------
md = f"""# Lu et al. 2023 (GSE193581) — v17p35 8-gene panel transfer (SLIM)

**Generated:** 2026-04-28
**Cohort:** GSE193581 — Lu et al. *J Clin Invest* 2023;133:e169653
**Design:** 10x Genomics 3' scRNA-seq, 23 samples (NORM + PTC + ATC)
**Strategy:** memory-slim — extract only 8 panel rows post-normalization, no full HVG/scale/PCA/UMAP

---

## Summary

| metric                              | value                                    |
|-------------------------------------|------------------------------------------|
| Samples used                        | {len(samples_data)} ({mal_meta.groupby('histology')['sample'].nunique().to_dict()}) |
| Total cells post-QC                 | {len(all_meta):,}                        |
| Malignant cells used                | {mal_meta.shape[0]:,}                    |
| Malignant by histology              | NORM={int(mal_meta['histology'].value_counts().get('NORM',0)):,}, PTC={int(mal_meta['histology'].value_counts().get('PTC',0)):,}, ATC={int(mal_meta['histology'].value_counts().get('ATC',0)):,} |
| Panel genes found                   | {len(mal_panel.index)}/8                |
| Panel genes missing                 | {set(PANEL) - set(mal_panel.index) or 'none'} |

## Cross-cohort statistics

* Kruskal-Wallis 3-group ({hlabels}) DM_score: H=**{H:.1f}**, p=**{p_kw:.2e}**
* Mann-Whitney PTC vs ATC: U={u:.0f}, p=**{p_mw:.2e}**
* chi-square histology × DM_class: chi²={chi2:.1f}, p=**{p_chi:.2e}**
* ARI(histology, DM_class) = **{ari:.3f}**

NORM median DM_score = {np.median(dm_norm):+.3f}; PTC median = {np.median(dm_ptc):+.3f}; ATC median = {np.median(dm_atc):+.3f}.

## Per-sample
```
{per_sample.round(3).to_string()}
```

## Files

* `/opt/thyroid-dash/project/results/v17_lu2023/GSE193581_*.tsv` + `.json`
* `/opt/thyroid-dash/project/reports/html/figs_interactive/v17/v17_lu2023_*.html` (4 figures)
* `/opt/thyroid-dash/project/logs/v17_lu2023_slim.log`

*author: Seungho Cook · seed=42*
"""
(RPT/"Lu2023_summary.md").write_text(md)
log.info(f"wrote {RPT}/Lu2023_summary.md")
log.info("=== DONE ===")
print(json.dumps({k:v for k,v in summary.items() if k != "per_sample_summary"}, indent=2, default=str))
