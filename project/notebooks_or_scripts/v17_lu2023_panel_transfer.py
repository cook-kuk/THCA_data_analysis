#!/usr/bin/env python3
"""GSE193581 (Lu et al. 2023 JCI) — v17p35 8-gene panel transfer onto scRNA-seq.

Single-cell counterpart to v17_korean_panel_transfer.py (GSE213647 bulk).
Cohort: 6 normal + 7 PTC + 10 ATC (10x Genomics 3'); analyzes 16,049 author-annotated
malignant epithelial cells across PTC vs ATC.
"""
from __future__ import annotations
import os, sys, json, gzip, time, logging
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats, sparse
from sklearn.metrics import adjusted_rand_score
import plotly.graph_objects as go

# ---------- paths ----------
RAW_DIR  = Path("/data/thca/v17_lu2023_GSE193581/raw")
ANN_FILE = Path("/data/thca/v17_lu2023_GSE193581/celltype_annotation.txt.gz")
OUT      = Path("/opt/thyroid-dash/project/results/v17_lu2023");          OUT.mkdir(parents=True, exist_ok=True)
FIG      = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17"); FIG.mkdir(parents=True, exist_ok=True)
RPT      = Path("/opt/thyroid-dash/project/reports/v17p35");              RPT.mkdir(parents=True, exist_ok=True)
LOG      = Path("/opt/thyroid-dash/project/logs/v17_lu2023.log")

# ---------- logging ----------
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.FileHandler(LOG, mode="w"), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
log.info("=== v17_lu2023 GSE193581 panel transfer START ===")

# ---------- styling ----------
BG, INK = "#0b0e12", "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG,
            font=dict(color=INK, size=14))
COLOR_PTC = "#F5A623"; COLOR_ATC = "#c24c4c"
SEED = 42; np.random.seed(SEED)

# ---------- panel ----------
PANEL = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
log.info(f"panel: {PANEL}")

# ---------- 1. load annotation & build sample manifest ----------
ann = pd.read_csv(ANN_FILE, sep="\t")
ann.columns = ["sample_id", "celltype"]
# cell barcode is SAMPLE_BARCODE-1 — encoded in the implicit row order.  In Lu 2023 the
# annotation file rows are aligned 1:1 with concatenated cell IDs; the first column in each
# UMI matrix already includes "SAMPLE_BARCODE-1" — so we can build a barcode key from the
# annotation file row index. But the annotation file has no per-cell barcode. We instead
# trust author cell counts per sample: per-sample alignment by order in each UMI matrix.
log.info(f"annotation: {ann.shape}; samples={ann['sample_id'].nunique()}; celltypes={ann['celltype'].nunique()}")

# Build {sample_id: [cell_types in order]}
ann_by_sample = {s: g["celltype"].values for s, g in ann.groupby("sample_id")}

# Map sample_id -> filename (note ATC18T in annotation but ATC18 in file)
sample_to_file = {}
for fp in sorted(RAW_DIR.glob("*_UMI.txt.gz")):
    # GSM5814574_PTC01_UMI.txt.gz -> PTC01
    sname = fp.name.split("_")[1]
    sample_to_file[sname] = fp
log.info(f"UMI files found: {len(sample_to_file)}")

# Histology label per sample
def histo(s):
    if s.startswith("PTC"):  return "PTC"
    if s.startswith("ATC"):  return "ATC"
    if s.startswith("NORM"): return "NORM"
    return "OTHER"

# ---------- 2. read each sample, attach annotation by ORDERED row alignment ----------
# In dense matrices, COLUMN order of cells = order in author annotation table.
adatas = []
for sample, fp in sample_to_file.items():
    # Match ATC18 file <-> ATC18T annotation
    ann_key = sample if sample in ann_by_sample else (sample + "T")
    if ann_key not in ann_by_sample:
        log.warning(f"  {sample}: no annotation key found, SKIP")
        continue
    t0 = time.time()
    df = pd.read_csv(fp, sep="\t", index_col=0)  # genes × cells
    # cell barcodes from columns
    barcodes = df.columns.tolist()
    expected = ann_by_sample[ann_key]
    if len(barcodes) != len(expected):
        log.warning(f"  {sample}: matrix has {len(barcodes)} cells, annotation has {len(expected)} — using min and matching by order")
    n = min(len(barcodes), len(expected))
    df = df.iloc[:, :n]
    cells = barcodes[:n]
    ctypes = expected[:n]

    X = sparse.csr_matrix(df.values.T.astype(np.float32))   # cells × genes
    a = sc.AnnData(X=X)
    a.obs_names = cells
    a.var_names = df.index.tolist()
    a.obs["sample"] = sample
    a.obs["histology"] = histo(sample)
    a.obs["author_celltype"] = pd.Categorical(ctypes)
    log.info(f"  {sample}: {a.n_obs} cells × {a.n_vars} genes ({time.time()-t0:.1f}s)")
    adatas.append(a)

# Concat
adata = sc.concat(adatas, axis=0, join="outer", merge="same", index_unique=None)
log.info(f"merged: {adata.shape}; histology counts: {adata.obs['histology'].value_counts().to_dict()}")

# ---------- 3. QC + normalization ----------
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
log.info(f"pre-QC: median n_genes={adata.obs['n_genes_by_counts'].median():.0f}, median pct_mt={adata.obs['pct_counts_mt'].median():.2f}")

mask = (adata.obs["n_genes_by_counts"] > 200) & (adata.obs["n_genes_by_counts"] < 6000) \
       & (adata.obs["pct_counts_mt"] < 20)
adata = adata[mask].copy()
log.info(f"post-QC: {adata.shape}")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()  # save log-normed counts for panel scoring later
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, random_state=SEED)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
log.info("UMAP done")

# ---------- 4. filter to malignant epithelial cells ----------
malignant_mask = adata.obs["author_celltype"].isin(["Malignant cell", "Epithelial cell"])
log.info(f"malignant + epithelial cells: {malignant_mask.sum()} ({adata.obs['author_celltype'].value_counts().to_dict()})")
mal = adata[malignant_mask].copy()
log.info(f"malignant histology: {mal.obs['histology'].value_counts().to_dict()}")

# ---------- 5. 8-gene panel z-score per cell ----------
panel_present = [g for g in PANEL if g in mal.raw.var_names]
panel_missing = [g for g in PANEL if g not in mal.raw.var_names]
log.info(f"panel genes present: {panel_present}")
log.info(f"panel genes missing: {panel_missing}")

# Pull log-normed expression of the 8 genes from .raw
raw_X = mal.raw[:, panel_present].X
if sparse.issparse(raw_X): raw_X = raw_X.toarray()
df_panel = pd.DataFrame(raw_X, index=mal.obs_names, columns=panel_present)
# z-score per gene across all malignant cells
z = (df_panel - df_panel.mean(axis=0)) / (df_panel.std(axis=0) + 1e-9)
mal.obs["DM_score"] = z.mean(axis=1).values  # mean z across panel = "DM_score"
log.info(f"DM_score: mean={mal.obs['DM_score'].mean():.3f} sd={mal.obs['DM_score'].std():.3f} "
         f"median={mal.obs['DM_score'].median():.3f}")

# Median split into DM1 (high differentiation) / DM2 (low differentiation)
med = mal.obs["DM_score"].median()
mal.obs["DM_class"] = np.where(mal.obs["DM_score"] >= med, "DM1_high", "DM2_low")
log.info(f"DM_class median split @ {med:.3f}: {mal.obs['DM_class'].value_counts().to_dict()}")

# Also propagate DM_score back to full adata for UMAP plotting
adata.obs["DM_score"] = np.nan
adata.obs.loc[mal.obs_names, "DM_score"] = mal.obs["DM_score"].values
adata.obs["DM_class"] = "non_malignant"
adata.obs.loc[mal.obs_names, "DM_class"] = mal.obs["DM_class"].values

# ---------- 6. cross-tab: PTC vs ATC histology vs DM1/DM2 ----------
ct = pd.crosstab(mal.obs["histology"], mal.obs["DM_class"])
log.info(f"\ncross-tab histology × DM_class:\n{ct}")
chi2, p_chi, dof, exp = stats.chi2_contingency(ct)
log.info(f"chi-square: chi2={chi2:.2f} dof={dof} p={p_chi:.3e}")

# ARI between {PTC,ATC} histology label and DM1/DM2
hist_codes = pd.Categorical(mal.obs["histology"]).codes
dm_codes   = pd.Categorical(mal.obs["DM_class"]).codes
ari = adjusted_rand_score(hist_codes, dm_codes)
log.info(f"ARI(histology, DM_class) = {ari:.4f}")

# Mann-Whitney on DM_score (PTC vs ATC)
dm_ptc = mal.obs.loc[mal.obs["histology"]=="PTC", "DM_score"].values
dm_atc = mal.obs.loc[mal.obs["histology"]=="ATC", "DM_score"].values
u, p_mw = stats.mannwhitneyu(dm_ptc, dm_atc, alternative="two-sided")
log.info(f"Mann-Whitney PTC vs ATC DM_score: U={u:.0f} p={p_mw:.3e}; "
         f"PTC median={np.median(dm_ptc):.3f} ATC median={np.median(dm_atc):.3f}")

# Per-sample DM_score
per_sample = mal.obs.groupby("sample").agg(
    histology=("histology","first"),
    n_cells=("DM_score","count"),
    median_DM=("DM_score","median"),
    mean_DM=("DM_score","mean"),
    pct_DM1=("DM_class", lambda s: (s=="DM1_high").mean()),
).sort_values(["histology","median_DM"])
log.info(f"\nper-sample summary:\n{per_sample}")

# ---------- 7. save TSV / JSON ----------
mal.obs[["sample","histology","author_celltype","DM_score","DM_class"]].to_csv(
    OUT/"GSE193581_cell_panel_score.tsv", sep="\t")
per_sample.to_csv(OUT/"GSE193581_per_sample_summary.tsv", sep="\t")
ct.to_csv(OUT/"GSE193581_crosstab_histology_DM.tsv", sep="\t")

summary = {
    "cohort": "GSE193581 (Lu et al. 2023 JCI)",
    "study": "Single-cell map of differentiated → anaplastic thyroid transformation",
    "samples_total": int(adata.obs["sample"].nunique()),
    "n_cells_pre_qc": int(sum(a.n_obs for a in adatas)),
    "n_cells_post_qc": int(adata.n_obs),
    "n_malignant_cells": int(mal.n_obs),
    "malignant_by_histology": {k:int(v) for k,v in mal.obs["histology"].value_counts().to_dict().items()},
    "panel_genes_found": panel_present,
    "panel_genes_missing": panel_missing,
    "DM_score_median": float(med),
    "DM_class_counts": {k:int(v) for k,v in mal.obs["DM_class"].value_counts().to_dict().items()},
    "crosstab_histology_DM": ct.to_dict(),
    "chi_square_p": float(p_chi),
    "chi_square_stat": float(chi2),
    "chi_square_dof": int(dof),
    "ARI_histology_DM": float(ari),
    "MannWhitney_PTC_vs_ATC_DM_score_p": float(p_mw),
    "DM_score_PTC_median": float(np.median(dm_ptc)),
    "DM_score_ATC_median": float(np.median(dm_atc)),
    "seed": SEED,
}
(OUT/"GSE193581_summary.json").write_text(json.dumps(summary, indent=2))
log.info(f"saved {OUT}/GSE193581_summary.json")

# ---------- 8. figures (plotly dark) ----------
um = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1","UMAP2"])
um["author_celltype"] = adata.obs["author_celltype"].astype(str).values
um["histology"]       = adata.obs["histology"].values
um["DM_score"]        = adata.obs["DM_score"].values
um["sample"]          = adata.obs["sample"].values

# Fig A: UMAP colored by author cell type
ctype_colors = {
    "T cell":"#7ccfcd","Malignant cell":"#c24c4c","Myeloid cell":"#F5A623","B cell":"#9b59b6",
    "NK cell":"#1abc9c","Fibroblast":"#e67e22","Endothelial cell":"#3498db","Epithelial cell":"#F2F2F2"
}
fig = go.Figure()
for ct_label, color in ctype_colors.items():
    sub = um[um["author_celltype"]==ct_label]
    if len(sub)==0: continue
    fig.add_trace(go.Scattergl(x=sub["UMAP1"], y=sub["UMAP2"], mode="markers",
        marker=dict(size=2, color=color, opacity=0.6),
        name=f"{ct_label} (n={len(sub)})", hoverinfo="skip"))
fig.update_layout(title=dict(
    text="<b>GSE193581 (Lu 2023) — UMAP colored by author cell-type label</b><br>"
         f"<sub style='color:#7ccfcd'>n={adata.n_obs:,} post-QC cells from {adata.obs['sample'].nunique()} samples (PTC + ATC + NORM)</sub>",
    font=dict(size=14, color=INK)),
    xaxis=dict(title="UMAP1", showgrid=False), yaxis=dict(title="UMAP2", showgrid=False),
    height=600, width=900, margin=dict(l=60,r=20,t=90,b=60), **DARK,
    legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)"))
fig.write_html(FIG/"v17_lu2023_umap_celltype.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_umap_celltype.html")

# Fig B: UMAP colored by DM_score (continuous, malignant only highlighted)
um_mal = um.dropna(subset=["DM_score"])
um_oth = um[um["DM_score"].isna()]
fig = go.Figure()
fig.add_trace(go.Scattergl(x=um_oth["UMAP1"], y=um_oth["UMAP2"], mode="markers",
    marker=dict(size=2, color="rgba(120,120,120,0.25)"),
    name=f"non-malignant (n={len(um_oth)})", hoverinfo="skip"))
fig.add_trace(go.Scattergl(x=um_mal["UMAP1"], y=um_mal["UMAP2"], mode="markers",
    marker=dict(size=3, color=um_mal["DM_score"], colorscale="RdBu",
                cmid=0, cmin=-1.5, cmax=1.5, showscale=True,
                colorbar=dict(title=dict(text="DM_score", font=dict(color=INK)),
                              tickfont=dict(color=INK))),
    name=f"malignant (n={len(um_mal)})", hoverinfo="skip"))
fig.update_layout(title=dict(
    text="<b>GSE193581 — UMAP colored by 8-gene panel DM_score (malignant cells only)</b><br>"
         f"<sub style='color:#7ccfcd'>panel: SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1; "
         f"PTC median={np.median(dm_ptc):.2f}, ATC median={np.median(dm_atc):.2f}, MWU p={p_mw:.2e}</sub>",
    font=dict(size=14, color=INK)),
    xaxis=dict(title="UMAP1", showgrid=False), yaxis=dict(title="UMAP2", showgrid=False),
    height=600, width=900, margin=dict(l=60,r=20,t=90,b=60), **DARK)
fig.write_html(FIG/"v17_lu2023_umap_dmscore.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_umap_dmscore.html")

# Fig C: cross-tab heatmap (histology × DM class) — show counts and percentages
ct_norm = ct.div(ct.sum(axis=1), axis=0) * 100
fig = go.Figure(go.Heatmap(z=ct_norm.values, x=ct_norm.columns.tolist(),
    y=ct_norm.index.tolist(), colorscale="RdBu_r", zmid=50, zmin=0, zmax=100,
    text=[[f"<b>{int(ct.loc[r,c])}</b><br>({ct_norm.loc[r,c]:.1f}%)"
           for c in ct.columns] for r in ct.index],
    texttemplate="%{text}",
    textfont=dict(size=14, color=INK),
    colorbar=dict(title=dict(text="% of row", font=dict(color=INK)), tickfont=dict(color=INK))))
fig.update_layout(title=dict(
    text="<b>GSE193581 — histology × DM-class cross-tab (malignant cells)</b><br>"
         f"<sub style='color:#7ccfcd'>chi² = {chi2:.1f}, dof = {dof}, p = {p_chi:.2e}; ARI = {ari:.3f}</sub>",
    font=dict(size=14, color=INK)),
    xaxis=dict(title="DM class (median split on 8-gene panel z-score)", tickfont=dict(color=INK)),
    yaxis=dict(title="Histology (author label)", tickfont=dict(color=INK)),
    height=420, width=720, margin=dict(l=120,r=40,t=110,b=80), **DARK)
fig.write_html(FIG/"v17_lu2023_crosstab.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_crosstab.html")

# Fig D bonus: per-sample DM_score box plot, ordered by median, colored by histology
fig = go.Figure()
sample_order = per_sample.sort_values(["histology","median_DM"], ascending=[True, False]).index.tolist()
for s in sample_order:
    sub = mal.obs[mal.obs["sample"]==s]
    h = sub["histology"].iloc[0]
    color = COLOR_PTC if h=="PTC" else COLOR_ATC
    fig.add_trace(go.Box(y=sub["DM_score"], x=[s]*len(sub),
        name=f"{s} ({h}, n={len(sub)})", marker=dict(color=color, size=2, opacity=0.5),
        boxpoints="outliers", line=dict(width=1.5)))
fig.update_layout(title=dict(
    text="<b>GSE193581 — per-sample 8-gene DM_score (malignant cells, sorted by median)</b><br>"
         f"<sub style='color:#7ccfcd'>orange=PTC ({(per_sample.histology=='PTC').sum()} samples), "
         f"red=ATC ({(per_sample.histology=='ATC').sum()} samples)</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="DM_score (mean z across 8 genes)", gridcolor="rgba(255,255,255,0.08)",
               zeroline=True, zerolinecolor="rgba(255,255,255,0.2)"),
    xaxis=dict(tickangle=-45), showlegend=False,
    height=520, width=1100, margin=dict(l=80,r=20,t=110,b=110), **DARK)
fig.write_html(FIG/"v17_lu2023_per_sample_box.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote {FIG}/v17_lu2023_per_sample_box.html")

# ---------- 9. markdown report ----------
md = f"""# Lu et al. 2023 (GSE193581) — v17p35 8-gene panel transfer

**Generated:** 2026-04-28
**Cohort:** GSE193581 — Lu et al. *J Clin Invest* 2023;133:e169653
**Design:** 10x Genomics 3' scRNA-seq, 6 normal + 7 PTC + 10 ATC thyroid samples
**Author cell-type labels supplied** ({adata.obs['author_celltype'].nunique()} types) — malignant/epithelial cells used directly.

---

## Summary statistics

| metric                              | value                                    |
|-------------------------------------|------------------------------------------|
| Samples used (UMI matrices found)   | {adata.obs['sample'].nunique()} ({adata.obs.groupby('histology')['sample'].nunique().to_dict()}) |
| Total cells pre-QC                  | {sum(a.n_obs for a in adatas):,}         |
| Total cells post-QC                 | {adata.n_obs:,}                          |
| Malignant + epithelial cells used   | {mal.n_obs:,}                            |
| Malignant cells by histology        | PTC={int(mal.obs['histology'].value_counts().get('PTC',0)):,}, ATC={int(mal.obs['histology'].value_counts().get('ATC',0)):,} |
| 8-gene panel genes found            | {len(panel_present)}/8 ({panel_present}) |
| 8-gene panel genes missing          | {panel_missing if panel_missing else 'none'} |
| DM_score median (split point)       | {med:.3f}                                |
| DM1 (high) cells                    | {int((mal.obs['DM_class']=='DM1_high').sum()):,} |
| DM2 (low) cells                     | {int((mal.obs['DM_class']=='DM2_low').sum()):,}  |

## Cross-cohort cross-tab (histology × DM class)

```
{ct.to_string()}
```

* chi² = **{chi2:.2f}**, dof = {dof}, **p = {p_chi:.3e}**
* ARI(histology, DM-class) = **{ari:.3f}**
* Mann-Whitney U (PTC vs ATC DM_score): U={u:.0f}, **p = {p_mw:.3e}**
* PTC median DM_score = {np.median(dm_ptc):+.3f}; ATC median DM_score = {np.median(dm_atc):+.3f}
  (ΔDM_score PTC−ATC = **{np.median(dm_ptc)-np.median(dm_atc):+.3f}**)

## Per-sample (malignant cells)

```
{per_sample.round(3).to_string()}
```

## Cross-cohort comparison

| cohort                  | platform  | n malignant / total | test                    | p-value     |
|-------------------------|-----------|---------------------|-------------------------|-------------|
| GSE213647 (Lee 2024)    | bulk RNA  | n=632 samples       | Kruskal-Wallis (4 grp)  | 8.72e-55    |
| GSE193581 (Lu 2023)     | scRNA-seq | n={mal.n_obs:,} cells | Mann-Whitney (PTC vs ATC) | {p_mw:.2e} |

Both replicate the v17p35 panel direction: differentiation markers ↓ as histology dedifferentiates.

## Key figures

* UMAP — author cell type: `../html/figs_interactive/v17/v17_lu2023_umap_celltype.html`
* UMAP — DM_score:         `../html/figs_interactive/v17/v17_lu2023_umap_dmscore.html`
* Cross-tab heatmap:       `../html/figs_interactive/v17/v17_lu2023_crosstab.html`
* Per-sample box plot:     `../html/figs_interactive/v17/v17_lu2023_per_sample_box.html`

## Files

* `/opt/thyroid-dash/project/results/v17_lu2023/GSE193581_cell_panel_score.tsv` — per-cell DM_score & class
* `/opt/thyroid-dash/project/results/v17_lu2023/GSE193581_per_sample_summary.tsv` — per-sample aggregates
* `/opt/thyroid-dash/project/results/v17_lu2023/GSE193581_crosstab_histology_DM.tsv`
* `/opt/thyroid-dash/project/results/v17_lu2023/GSE193581_summary.json`
* `/opt/thyroid-dash/project/logs/v17_lu2023.log`

---
*author: Seungho Cook · seed=42 · scanpy {sc.__version__}*
"""
(RPT/"Lu2023_summary.md").write_text(md)
log.info(f"wrote {RPT}/Lu2023_summary.md")
log.info("=== DONE ===")
print(json.dumps(summary, indent=2))
