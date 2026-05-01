#!/usr/bin/env python
"""v17p35 - GSE241184 single-cell projection of 8-gene DM panel.

Author: Seungho Cook (2026-04-28)

NOTE / COHORT MISMATCH (recorded by author):
  Task brief described GSE241184 as "50 tumor + 14 normal Chinese single-cell
  thyroid atlas". GEO record (PubMed 38061122, Chen W & Zhong S, Nanjing 2023)
  shows GSE241184 is in fact a 3-sample case study from a single 17-year-old
  female PTC patient: 1 thyroid tumor (TT), 1 adjacent normal (NT), 1 lymph-
  node metastasis (LN). We proceed with single-cell DM-axis projection and
  per-tissue (rather than per-patient) statistics, and flag this in the report.
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import mannwhitneyu, kruskal, spearmanr

DATA = Path("/data/thca/v17_gse241184")
OUT  = Path("/opt/thyroid-dash/project/results/v17_gse241184")
FIG  = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
LOGD = Path("/opt/thyroid-dash/project/logs")
for p in (OUT, FIG, LOGD): p.mkdir(parents=True, exist_ok=True)

LOG = LOGD / "v17_gse241184.log"
logging.basicConfig(
    filename=LOG, level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger("gse241184")

RNG = 42
np.random.seed(RNG)
sc.settings.seed = RNG
sc.settings.verbosity = 1

PANEL = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]

SAMPLES = {
    "TT": ("Thyroid_tumor",   "Tumor"),
    "NT": ("Normal_thyroid",  "Normal"),
    "LN": ("Lymph_node",      "LN_Met"),
}

# Marker panels for cell-type labelling
MARKERS = {
    "Epithelial": ["EPCAM", "KRT18", "KRT8", "KRT19"],
    "Thyrocyte" : ["TG", "TPO", "PAX8", "TSHR"],
    "T_cell"    : ["CD3D", "CD3E", "CD3G", "CD2"],
    "B_cell"    : ["CD79A", "CD79B", "MS4A1", "IGHM"],
    "NK_cell"   : ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "Myeloid"   : ["LYZ", "CD68", "CD14", "C1QA", "C1QB"],
    "Endothelial":["PECAM1", "VWF", "CDH5"],
    "Fibroblast": ["COL1A1", "COL3A1", "DCN", "LUM"],
}

# Plotly dark style ----------------------------------------------------------
PLOTLY_BG   = "#0b0e12"
PLOTLY_FONT = "#F2F2F2"

def style(fig, title=""):
    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor=PLOTLY_BG, plot_bgcolor=PLOTLY_BG,
        font=dict(color=PLOTLY_FONT, family="Inter, system-ui, sans-serif"),
        title=dict(text=title, x=0.02, xanchor="left",
                   font=dict(color=PLOTLY_FONT, size=18)),
        margin=dict(l=60, r=30, t=60, b=50),
    )
    return fig

# 1. Load 10x ----------------------------------------------------------------
def load_sample(key):
    folder, lbl = SAMPLES[key]
    path = DATA / folder
    log.info("loading %s from %s", key, path)
    a = sc.read_mtx(path / "matrix.mtx.gz").T  # cells x genes
    barcodes = pd.read_csv(path / "barcodes.tsv.gz", header=None, sep="\t")[0].values
    feats    = pd.read_csv(path / "features.tsv.gz", header=None, sep="\t")
    a.obs_names = [f"{key}_{b}" for b in barcodes]
    a.var_names = feats[1].values  # symbols
    a.var["ensembl"] = feats[0].values
    a.var_names_make_unique()
    a.X = sparse.csr_matrix(a.X)
    a.obs["sample"]    = key
    a.obs["tissue"]    = lbl
    return a

def main():
    log.info("=== v17p35 GSE241184 single-cell DM panel ===")
    adatas = {k: load_sample(k) for k in SAMPLES}
    for k, a in adatas.items():
        log.info("  %s shape %s", k, a.shape)

    # Concatenate
    A = ad.concat(adatas, join="outer", index_unique=None)
    A.obs["sample"] = A.obs["sample"].astype("category")
    A.obs["tissue"] = A.obs["tissue"].astype("category")
    log.info("concat shape %s, dtype %s", A.shape, type(A.X).__name__)

    # 2. QC ------------------------------------------------------------------
    A.var["mt"] = A.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None,
                               log1p=False)
    pre = A.shape
    A = A[(A.obs["pct_counts_mt"] < 20) &
          (A.obs["n_genes_by_counts"] >= 200) &
          (A.obs["n_genes_by_counts"] <= 6000)].copy()
    log.info("after QC %s -> %s", pre, A.shape)
    sc.pp.filter_genes(A, min_cells=3)
    log.info("after gene filter %s", A.shape)

    raw_counts = A.X.copy()  # keep counts for marker scoring

    # 3. Normalize -----------------------------------------------------------
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.layers["lognorm"] = A.X.copy()
    sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat",
                                batch_key="sample")
    A_hvg = A[:, A.var.highly_variable].copy()
    sc.pp.scale(A_hvg, max_value=10)
    sc.tl.pca(A_hvg, n_comps=50, random_state=RNG)

    # Harmony integration ----------------------------------------------------
    try:
        sc.external.pp.harmony_integrate(A_hvg, key="sample",
                                         basis="X_pca",
                                         adjusted_basis="X_pca_harmony",
                                         random_state=RNG)
        emb_key = "X_pca_harmony"
    except Exception as e:
        log.warning("harmony failed (%s) - falling back to raw PCA", e)
        emb_key = "X_pca"

    sc.pp.neighbors(A_hvg, use_rep=emb_key, n_neighbors=15, random_state=RNG)
    sc.tl.umap(A_hvg, random_state=RNG)
    sc.tl.leiden(A_hvg, resolution=0.6, random_state=RNG)
    A.obsm["X_umap"]    = A_hvg.obsm["X_umap"]
    A.obs["leiden"]     = A_hvg.obs["leiden"].values

    # 4. Cell-type annotation via marker scores ------------------------------
    for ct, genes in MARKERS.items():
        present = [g for g in genes if g in A.var_names]
        if not present:
            A.obs[f"score_{ct}"] = 0.0
            continue
        sc.tl.score_genes(A, gene_list=present, score_name=f"score_{ct}",
                          random_state=RNG)
    score_cols = [f"score_{ct}" for ct in MARKERS]
    A.obs["celltype"] = A.obs[score_cols].idxmax(axis=1).str.replace("score_", "")
    # promote thyrocyte over generic epithelial when both elevated
    thy_mask = (A.obs["score_Thyrocyte"] > 0.25) & (A.obs["celltype"]=="Epithelial")
    A.obs.loc[thy_mask, "celltype"] = "Thyrocyte"

    log.info("cell-type counts:\n%s", A.obs["celltype"].value_counts().to_string())

    # 5. 8-gene DM_score per cell -------------------------------------------
    present_panel = [g for g in PANEL if g in A.var_names]
    missing = [g for g in PANEL if g not in A.var_names]
    log.info("panel present: %s; missing: %s", present_panel, missing)
    # z-score per gene across cells
    expr = A[:, present_panel].layers["lognorm"]
    if sparse.issparse(expr): expr = expr.toarray()
    expr_df = pd.DataFrame(expr, columns=present_panel, index=A.obs_names)
    z = (expr_df - expr_df.mean()) / expr_df.std(ddof=0).replace(0, np.nan)
    z = z.fillna(0.0)
    A.obs["DM_score"] = z.mean(axis=1).values

    # restrict DM1/DM2 split to malignant epithelial (Thyrocyte) cells
    thy = A.obs["celltype"]=="Thyrocyte"
    n_thy = int(thy.sum())
    log.info("Thyrocyte cells: %d", n_thy)
    if n_thy >= 50:
        thy_med = float(A.obs.loc[thy, "DM_score"].median())
    else:
        thy_med = float(A.obs["DM_score"].median())
    A.obs["DM_class"] = np.where(A.obs["DM_score"] >= thy_med, "DM1", "DM2")
    A.obs.loc[~thy, "DM_class"] = "non_thyrocyte"

    # 6. Per-sample summary --------------------------------------------------
    rows = []
    for k, (folder, lbl) in SAMPLES.items():
        sub = A.obs[A.obs["sample"]==k]
        n_total = len(sub)
        n_thy   = int((sub["celltype"]=="Thyrocyte").sum())
        thy_sub = sub[sub["celltype"]=="Thyrocyte"]
        n_dm1   = int((thy_sub["DM_class"]=="DM1").sum())
        n_dm2   = int((thy_sub["DM_class"]=="DM2").sum())
        immune_ct = ["T_cell","B_cell","NK_cell","Myeloid"]
        n_imm   = int(sub["celltype"].isin(immune_ct).sum())
        rows.append({
            "sample":k, "tissue":lbl, "n_cells":n_total,
            "n_thyrocyte":n_thy, "thy_DM1":n_dm1, "thy_DM2":n_dm2,
            "DM1_pct_of_thy": (n_dm1/max(n_thy,1))*100,
            "immune_n": n_imm,
            "immune_pct": (n_imm/max(n_total,1))*100,
            "DM_score_mean_thy": float(thy_sub["DM_score"].mean()) if n_thy else np.nan,
            "DM_score_median_thy": float(thy_sub["DM_score"].median()) if n_thy else np.nan,
        })
    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(OUT/"per_sample_summary.tsv", sep="\t", index=False)
    log.info("per-sample summary:\n%s", per_sample.to_string(index=False))

    # 7. Cross-tissue test on DM_score (Thyrocytes only) ---------------------
    thy_obs = A.obs[A.obs["celltype"]=="Thyrocyte"]
    groups  = [thy_obs.loc[thy_obs["sample"]==k, "DM_score"].values
               for k in SAMPLES if (thy_obs["sample"]==k).any()]
    kw_p = float("nan")
    if len(groups) >= 2 and all(len(g)>1 for g in groups):
        try:
            kw_stat, kw_p = kruskal(*groups)
        except ValueError as e:
            log.warning("kruskal: %s", e); kw_stat = float("nan")
    # NT vs (TT+LN) Mann-Whitney
    mw_p = float("nan"); mw_stat = float("nan")
    nt = thy_obs.loc[thy_obs["sample"]=="NT","DM_score"].values
    tu = thy_obs.loc[thy_obs["sample"].isin(["TT","LN"]),"DM_score"].values
    if len(nt) > 5 and len(tu) > 5:
        mw_stat, mw_p = mannwhitneyu(nt, tu, alternative="greater")

    # 8. Patient-level dominance hypothesis (n=1 patient -> just report) ------
    dom_pct = per_sample["DM1_pct_of_thy"].max()
    dominant_label = per_sample.loc[per_sample["DM1_pct_of_thy"].idxmax(), "sample"]
    # Spearman immune vs DM1% across the 3 samples (under-powered, reported only)
    try:
        sp_rho, sp_p = spearmanr(per_sample["DM1_pct_of_thy"], per_sample["immune_pct"])
    except Exception:
        sp_rho, sp_p = float("nan"), float("nan")

    summary = {
        "geo": "GSE241184",
        "pubmed": "38061122",
        "n_samples": int(len(SAMPLES)),
        "n_patients": 1,
        "cohort_mismatch_note":
            "Task brief expected '50 tumor + 14 normal'; GEO record is single "
            "17-year-old PTC patient (TT/NT/LN). Analysis adapted accordingly.",
        "tissues": list(per_sample["tissue"]),
        "n_cells_total_post_qc": int(A.n_obs),
        "n_genes_post_filter": int(A.n_vars),
        "panel_genes_present": present_panel,
        "panel_genes_missing": missing,
        "n_thyrocyte_total": int((A.obs["celltype"]=="Thyrocyte").sum()),
        "DM_split_median_used": float(thy_med),
        "kruskal_DM_score_across_samples_p": (None if np.isnan(kw_p) else float(kw_p)),
        "MWU_normal_vs_tumor_DM_p": (None if np.isnan(mw_p) else float(mw_p)),
        "MWU_alt": "greater (Normal > Tumor expected)",
        "patient_dominance_pct_max": float(dom_pct),
        "patient_dominance_sample": str(dominant_label),
        "patient_dominance_threshold_70pct": bool(dom_pct >= 70),
        "spearman_DM1pct_vs_immunepct_rho": (None if np.isnan(sp_rho) else float(sp_rho)),
        "spearman_DM1pct_vs_immunepct_p":   (None if np.isnan(sp_p)   else float(sp_p)),
        "spearman_n": int(len(per_sample)),
        "spearman_caveat":
            "n=3 samples from 1 patient -> Spearman is descriptive only, not inferential.",
        "random_state": RNG,
    }
    with open(OUT/"summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    log.info("summary: %s", json.dumps(summary, indent=2))

    # Per-cell metadata dump
    keep_cols = ["sample","tissue","leiden","celltype","DM_score","DM_class",
                 "n_genes_by_counts","total_counts","pct_counts_mt"]
    A.obs[keep_cols].to_csv(OUT/"per_cell_metadata.tsv.gz",
                            sep="\t", compression="gzip")

    # 9. Plotly figures ------------------------------------------------------
    import plotly.express as px
    import plotly.graph_objects as go

    umap = pd.DataFrame(A.obsm["X_umap"], columns=["UMAP1","UMAP2"],
                        index=A.obs_names)
    umap = umap.join(A.obs[keep_cols])

    # F1: UMAP by tissue
    fig1 = px.scatter(umap, x="UMAP1", y="UMAP2", color="tissue",
                      color_discrete_map={"Normal":"#3aa6ff","Tumor":"#ff5277",
                                          "LN_Met":"#f4c542"},
                      opacity=0.6, render_mode="webgl",
                      hover_data=["celltype","DM_score"])
    fig1.update_traces(marker=dict(size=3, line=dict(width=0)))
    style(fig1, "GSE241184 UMAP - tissue (Harmony-integrated)")
    fig1.write_html(FIG/"v17_gse241184_umap_tissue.html",
                    include_plotlyjs="cdn", full_html=True)

    # F2: UMAP by DM_score
    fig2 = px.scatter(umap, x="UMAP1", y="UMAP2", color="DM_score",
                      color_continuous_scale="RdBu_r",
                      range_color=[-2,2],
                      opacity=0.7, render_mode="webgl",
                      hover_data=["sample","celltype"])
    fig2.update_traces(marker=dict(size=3, line=dict(width=0)))
    style(fig2, "GSE241184 UMAP - 8-gene DM_score (per-cell)")
    fig2.write_html(FIG/"v17_gse241184_umap_dmscore.html",
                    include_plotlyjs="cdn", full_html=True)

    # F2b: UMAP by celltype
    fig2b = px.scatter(umap, x="UMAP1", y="UMAP2", color="celltype",
                       opacity=0.6, render_mode="webgl",
                       hover_data=["sample"])
    fig2b.update_traces(marker=dict(size=3, line=dict(width=0)))
    style(fig2b, "GSE241184 UMAP - cell type (marker-score)")
    fig2b.write_html(FIG/"v17_gse241184_umap_celltype.html",
                     include_plotlyjs="cdn", full_html=True)

    # F3: per-sample DM1% bar
    bar = per_sample.copy().sort_values("DM1_pct_of_thy", ascending=False)
    fig3 = go.Figure()
    fig3.add_bar(x=bar["sample"]+" ("+bar["tissue"]+")",
                 y=bar["DM1_pct_of_thy"],
                 marker_color=["#3aa6ff" if t=="Normal" else
                               "#ff5277" if t=="Tumor" else "#f4c542"
                               for t in bar["tissue"]],
                 text=[f"{v:.1f}%" for v in bar["DM1_pct_of_thy"]],
                 textposition="outside")
    fig3.add_hline(y=50, line_dash="dot", line_color="#888",
                   annotation_text="50% (median)", annotation_position="top right")
    fig3.add_hline(y=70, line_dash="dash", line_color="#f4c542",
                   annotation_text="70% dominance", annotation_position="bottom right")
    fig3.update_yaxes(title="DM1 % of Thyrocytes", range=[0,110])
    style(fig3, "GSE241184 per-sample DM1 dominance (Thyrocytes)")
    fig3.write_html(FIG/"v17_gse241184_dm1_per_sample.html",
                    include_plotlyjs="cdn", full_html=True)

    # F4: immune vs DM1% scatter (n=3)
    fig4 = go.Figure()
    fig4.add_scatter(
        x=per_sample["DM1_pct_of_thy"], y=per_sample["immune_pct"],
        mode="markers+text", text=per_sample["sample"]+"/"+per_sample["tissue"],
        textposition="top center",
        marker=dict(size=14,
                    color=["#3aa6ff" if t=="Normal" else
                           "#ff5277" if t=="Tumor" else "#f4c542"
                           for t in per_sample["tissue"]],
                    line=dict(width=1, color="#222")))
    fig4.update_xaxes(title="DM1 % of Thyrocytes")
    fig4.update_yaxes(title="Immune cell fraction (%)")
    rho_txt = ("Spearman rho=NA" if sp_rho is None or np.isnan(sp_rho)
               else f"Spearman rho={sp_rho:.3f}, p={sp_p:.3g} (n=3, descriptive)")
    style(fig4, f"GSE241184 immune-DM1 (n=3) | {rho_txt}")
    fig4.write_html(FIG/"v17_gse241184_immune_vs_dm1.html",
                    include_plotlyjs="cdn", full_html=True)

    # F5: violin DM_score per tissue (Thyrocytes)
    thy_df = umap[umap["celltype"]=="Thyrocyte"].copy()
    fig5 = go.Figure()
    for t, col in [("Normal","#3aa6ff"),("Tumor","#ff5277"),("LN_Met","#f4c542")]:
        sub = thy_df[thy_df["tissue"]==t]["DM_score"]
        if len(sub) == 0: continue
        fig5.add_violin(y=sub, name=f"{t} (n={len(sub)})",
                        box_visible=True, meanline_visible=True,
                        line_color=col, fillcolor=col, opacity=0.55)
    p_str = "p=NA" if mw_p is None or np.isnan(mw_p) else f"MWU(N>T) p={mw_p:.3g}"
    kw_str = "KW=NA" if kw_p is None or np.isnan(kw_p) else f"KW p={kw_p:.3g}"
    style(fig5, f"GSE241184 DM_score by tissue (Thyrocytes) | {kw_str} | {p_str}")
    fig5.update_yaxes(title="DM_score (per cell)")
    fig5.write_html(FIG/"v17_gse241184_dmscore_violin.html",
                    include_plotlyjs="cdn", full_html=True)

    log.info("figures written to %s", FIG)
    log.info("DONE")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
