"""
Track 12 — scRNA HLA-I/II module per cell type in thyroid cohorts.

Boundary: HLA gene-expression module ONLY, not allele genotype.
See project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md (Section 1.1).

Inputs (already processed locally):
  - /data/thca/scrna/processed/classical_baseline.h5ad        (GSE184362, Pu et al. PTC scRNA, leiden-clustered)
  - /data/thca/scrna/processed/classical_baseline_F12.h5ad    (GSE193581, Lu 2023 PTC+ATC scRNA, BRAF/RAS-stratified)
  - /home/seungho/personal/THCA_data_analysis/project/data/processed/GSE250521/<gsm>/<gsm>.scored.h5ad
                                                              (16 Visium spots × stages = N/PTC/LPTC/ATC, full var, has DM1_like_score)
  - /home/seungho/personal/THCA_data_analysis/project/results/pantheonos_demo/<sample>/visium_with_mapped.h5ad
                                                              (mapped_celltype from MOSCOT scRNA reference)

Outputs land in:
  /home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track12_scrna_hla_celltype/
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# ---------- paths ----------
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track12_scrna_hla_celltype"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

PU_PATH = Path("/data/thca/scrna/processed/classical_baseline.h5ad")  # GSE184362
LU_PATH = Path("/data/thca/scrna/processed/classical_baseline_F12.h5ad")  # GSE193581
VISIUM_DIR = ROOT / "project/data/processed/GSE250521"  # full-gene scored
PANTHEON_DIR = ROOT / "project/results/pantheonos_demo"  # mapped_celltype

# ---------- gene lists ----------
HLA_I = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP",
         "NLRC5", "IRF1", "PSMB8", "PSMB9", "ERAP1", "ERAP2",
         "HLA-E", "HLA-F", "HLA-G", "CALR", "CANX", "PDIA3"]
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
          "HLA-DQB1", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
          "CIITA", "CD74"]

DM1_GENES = ["TG", "TPO", "TSHR", "SLC5A5", "DUOX1", "DUOX2", "DIO1", "DIO2"]
TLS_GENES = ["CXCL13", "CCR6", "CCR7", "MS4A1", "CD79A"]

# canonical lineage markers (single-cell-grade, conservative)
LINEAGE_MARKERS = {
    "Thyrocyte": ["TG", "TPO", "TSHR", "PAX8", "SLC5A5", "FOXE1", "NKX2-1"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2"],
    "B_cell": ["MS4A1", "CD19", "CD79A", "CD79B", "BANK1"],
    "Plasma": ["MZB1", "JCHAIN", "IGHG1", "IGHA1", "XBP1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "AIF1"],
    "DC": ["CLEC9A", "CLEC10A", "FCER1A", "CD1C", "LAMP3"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "ENG"],
    "Fibroblast": ["DCN", "COL1A1", "COL1A2", "LUM", "PDGFRA"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
}


def _ensure_symbol_index(a):
    """Some processed h5ads (Pu) use Ensembl IDs with a gene_symbol column. Make
    a copy with symbols on var_names, deduplicated by sum, lazily."""
    if "gene_symbol" not in a.var.columns:
        return a
    if any(v.startswith("HLA-") or v == "B2M" for v in a.var_names[:5000]):
        return a
    sym = a.var["gene_symbol"].astype(str).values
    # deduplicate symbols (keep first occurrence; unique enough for our targets)
    seen = {}
    keep_idx = []
    new_names = []
    for i, s in enumerate(sym):
        if s == "" or s == "nan":
            continue
        if s in seen:
            continue
        seen[s] = i
        keep_idx.append(i)
        new_names.append(s)
    a2 = a[:, keep_idx].copy()
    a2.var_names = new_names
    a2.var.index = pd.Index(new_names)
    return a2


def _genes_in(a, glist):
    s = set(a.var_names.tolist())
    return [g for g in glist if g in s]


def score_modules(a, label, hla_i=HLA_I, hla_ii=HLA_II):
    """Add hla_I_score, hla_II_score, plus per-lineage module scores via sc.tl.score_genes."""
    h1 = _genes_in(a, hla_i)
    h2 = _genes_in(a, hla_ii)
    sc.tl.score_genes(a, gene_list=h1, score_name="HLA_I_score", use_raw=False)
    sc.tl.score_genes(a, gene_list=h2, score_name="HLA_II_score", use_raw=False)
    info = {"dataset": label, "n_HLA_I_in_var": len(h1), "n_HLA_II_in_var": len(h2)}
    for lname, gs in LINEAGE_MARKERS.items():
        gs_p = _genes_in(a, gs)
        if len(gs_p) >= 2:
            sc.tl.score_genes(a, gene_list=gs_p, score_name=f"score_{lname}",
                              use_raw=False)
        else:
            a.obs[f"score_{lname}"] = np.nan
        info[f"n_{lname}_markers_in_var"] = len(gs_p)
    return info


def assign_celltype(a):
    """Assign coarse cell type per cell by argmax of lineage scores; require >= 0 and margin."""
    score_cols = [f"score_{k}" for k in LINEAGE_MARKERS.keys()]
    M = a.obs[score_cols].values
    best_ix = np.argmax(M, axis=1)
    best = M[np.arange(M.shape[0]), best_ix]
    # second-best
    M2 = M.copy()
    M2[np.arange(M.shape[0]), best_ix] = -np.inf
    second = M2.max(axis=1)
    margin = best - second
    labels = np.array(list(LINEAGE_MARKERS.keys()))[best_ix]
    # require best > 0 and margin > 0.05; else 'Unassigned'
    bad = (best < 0) | (margin < 0.05)
    labels = labels.astype(object)
    labels[bad] = "Unassigned"
    a.obs["celltype"] = pd.Categorical(labels)
    return a


def dotplot_celltype_genes(a, label, ds_short):
    """Custom dotplot: cell types × HLA-I + HLA-II genes (mean expr & pct expressing)."""
    genes = [g for g in HLA_I + HLA_II if g in a.var_names]
    cts = sorted([c for c in a.obs["celltype"].unique() if c != "Unassigned"])
    X = a[:, genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    df_rows = []
    for ct in cts:
        idx = (a.obs["celltype"] == ct).values
        if idx.sum() < 5:
            continue
        sub = X[idx]
        mean = sub.mean(axis=0)
        pct = (sub > 0).mean(axis=0) * 100.0
        for j, g in enumerate(genes):
            df_rows.append({"celltype": ct, "gene": g,
                            "mean_expr": float(mean[j]),
                            "pct_pos": float(pct[j]),
                            "module": "HLA-I" if g in HLA_I else "HLA-II"})
    df = pd.DataFrame(df_rows)
    df.to_csv(TAB / f"T2_{ds_short}_celltype_x_HLA_gene_dotplot.csv", index=False)

    # plot
    fig, ax = plt.subplots(figsize=(0.45 * len(genes) + 2.5, 0.55 * len(cts) + 2))
    # normalize mean per gene to 0..1 for color
    g_means = df.pivot_table(index="celltype", columns="gene", values="mean_expr").reindex(columns=genes)
    g_pct = df.pivot_table(index="celltype", columns="gene", values="pct_pos").reindex(columns=genes)
    cells = g_means.index.tolist()
    norm = (g_means - g_means.min()) / (g_means.max() - g_means.min() + 1e-12)
    for i, ct in enumerate(cells):
        for j, g in enumerate(genes):
            mv = norm.loc[ct, g] if not pd.isna(norm.loc[ct, g]) else 0
            sz = g_pct.loc[ct, g] if not pd.isna(g_pct.loc[ct, g]) else 0
            color = plt.cm.Reds(mv)
            ax.scatter(j, i, s=10 + sz * 4, c=[color], edgecolors="grey", linewidths=0.3)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(cells)))
    ax.set_yticklabels(cells, fontsize=8)
    # vertical separator between HLA-I and HLA-II
    nI = sum(1 for g in genes if g in HLA_I)
    ax.axvline(nI - 0.5, ls="--", color="grey", lw=0.5)
    ax.text(nI / 2, len(cells), "HLA-I", ha="center", fontsize=8, fontweight="bold")
    ax.text(nI + (len(genes) - nI) / 2, len(cells), "HLA-II", ha="center", fontsize=8, fontweight="bold")
    ax.set_title(f"{label} — cell type x HLA-I/II gene\n(mean expr color, pct expressing size)\n"
                 "HLA gene-expression module — not allele genotype.", fontsize=9)
    plt.tight_layout()
    plt.savefig(FIG / f"F2_{ds_short}_celltype_x_HLA_gene_dotplot.png", dpi=180, bbox_inches="tight")
    plt.close()


def violin_module_per_celltype(a, label, ds_short):
    df = a.obs[["celltype", "HLA_I_score", "HLA_II_score"]].copy()
    df = df[df["celltype"] != "Unassigned"]
    long = df.melt(id_vars="celltype", value_vars=["HLA_I_score", "HLA_II_score"],
                   var_name="module", value_name="score")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)
    for ax, mod in zip(axes, ["HLA_I_score", "HLA_II_score"]):
        sub = long[long["module"] == mod]
        order = (df.groupby("celltype")[mod].median().sort_values(ascending=False).index.tolist())
        sns.violinplot(data=sub, x="celltype", y="score", order=order, ax=ax,
                       inner="quartile", cut=0, linewidth=0.6, palette="Set2")
        ax.set_title(f"{label} — {mod} per cell type")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle("HLA gene-expression module — not allele genotype.", fontsize=9, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG / f"F3_{ds_short}_HLA_module_violin.png", dpi=180, bbox_inches="tight")
    plt.close()

    # table
    tbl = (df.groupby("celltype")[["HLA_I_score", "HLA_II_score"]]
              .agg(["mean", "median", "std", "count"]))
    tbl.columns = ["_".join(c) for c in tbl.columns]
    tbl.to_csv(TAB / f"T3_{ds_short}_HLA_module_per_celltype.csv")
    return tbl


def condition_stratified(a, label, ds_short, cond_col):
    """Boxplot per cell type per condition for HLA-I and HLA-II module."""
    if cond_col not in a.obs.columns:
        return None
    df = a.obs[["celltype", cond_col, "HLA_I_score", "HLA_II_score"]].copy()
    df = df[df["celltype"] != "Unassigned"]
    df.columns = ["celltype", "condition", "HLA_I_score", "HLA_II_score"]
    if df["condition"].nunique() < 2:
        return None
    # Compute per (celltype, condition) summary
    summ = (df.groupby(["celltype", "condition"])
              [["HLA_I_score", "HLA_II_score"]]
              .agg(["mean", "median", "count"]))
    summ.columns = ["_".join(c) for c in summ.columns]
    summ.to_csv(TAB / f"T4_{ds_short}_celltype_condition_HLA.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    for ax, mod in zip(axes, ["HLA_I_score", "HLA_II_score"]):
        sns.boxplot(data=df, x="celltype", y=mod, hue="condition",
                    ax=ax, fliersize=1, linewidth=0.6, palette="Set1")
        ax.set_title(f"{label} — {mod} by cell type and {cond_col}")
        ax.tick_params(axis="x", rotation=45)
        ax.set_xlabel("")
    fig.suptitle("HLA gene-expression module — not allele genotype.", fontsize=9, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG / f"F4_{ds_short}_celltype_condition_box.png", dpi=180, bbox_inches="tight")
    plt.close()
    return summ


def thyrocyte_pseudobulk_dm1_corr(a, label, ds_short, sample_col="sample"):
    """
    For each sample, take only Thyrocyte cells, average normalized expression →
    pseudobulk vector. Compute DM1 score (mean of DM1 genes) and HLA-I/II module
    means from this sample-level thyrocyte-only vector. Correlate across samples.
    """
    df = a.obs.copy()
    mask = (df["celltype"] == "Thyrocyte").values
    if mask.sum() < 50:
        return None
    samples = df.loc[mask, sample_col].astype(str).unique().tolist()
    rows = []
    Xnorm = a.X
    if hasattr(Xnorm, "toarray"):
        # we will subset with numpy via sparse later
        pass
    var = a.var_names
    dm1_g = [g for g in DM1_GENES if g in var]
    h1_g = [g for g in HLA_I if g in var]
    h2_g = [g for g in HLA_II if g in var]
    for s in samples:
        sel = (df[sample_col].astype(str) == s) & (df["celltype"] == "Thyrocyte")
        ix = np.where(sel.values)[0]
        if len(ix) < 10:
            continue
        sub = a.X[ix, :]
        if hasattr(sub, "toarray"):
            sub = sub.toarray()
        mean_expr = pd.Series(sub.mean(axis=0), index=var)
        rows.append({"sample": s, "n_cells": int(len(ix)),
                     "DM1_score": float(mean_expr[dm1_g].mean()),
                     "HLA_I_score": float(mean_expr[h1_g].mean()),
                     "HLA_II_score": float(mean_expr[h2_g].mean()),
                     })
    if not rows:
        return None
    pb = pd.DataFrame(rows)
    pb.to_csv(TAB / f"T5_{ds_short}_thyrocyte_pseudobulk.csv", index=False)
    out = {}
    for mod in ["HLA_I_score", "HLA_II_score"]:
        if pb["DM1_score"].std() < 1e-6 or pb[mod].std() < 1e-6 or len(pb) < 4:
            rho, p = np.nan, np.nan
        else:
            rho, p = stats.spearmanr(pb["DM1_score"], pb[mod])
        out[f"DM1_x_{mod}"] = {"rho": float(rho) if rho == rho else None,
                                "p": float(p) if p == p else None,
                                "n_samples": int(len(pb))}
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, mod in zip(axes, ["HLA_I_score", "HLA_II_score"]):
        sns.regplot(data=pb, x="DM1_score", y=mod, ax=ax,
                    scatter_kws={"s": 70, "color": "#cc4444"},
                    line_kws={"color": "grey"})
        rr = out[f"DM1_x_{mod}"]
        ax.set_title(f"{label} thyrocyte pseudobulk\n{mod}: rho={rr['rho']}, p={rr['p']}, n={rr['n_samples']}",
                     fontsize=9)
    fig.suptitle("Thyrocyte-only pseudobulk DM1 x HLA module — gene-expression module, not allele.",
                 fontsize=9, y=1.02)
    plt.tight_layout()
    plt.savefig(FIG / f"F5_{ds_short}_thyrocyte_pseudobulk_DM1xHLA.png", dpi=180, bbox_inches="tight")
    plt.close()
    return out


def decomposition_table(a, label, ds_short):
    """Approximate decomposition: per-sample bulk-pseudobulk DM1 x HLA-I r,
    compared to thyrocyte-only and immune-only pseudobulk. Fraction = how much
    each subset reproduces the bulk r."""
    df = a.obs.copy()
    if "sample" not in df.columns:
        return None
    var = a.var_names
    dm1_g = [g for g in DM1_GENES if g in var]
    h1_g = [g for g in HLA_I if g in var]
    h2_g = [g for g in HLA_II if g in var]

    def _pb(mask):
        rows = []
        for s in df.loc[mask, "sample"].astype(str).unique().tolist():
            sel = mask & (df["sample"].astype(str) == s).values
            ix = np.where(sel)[0]
            if len(ix) < 10:
                continue
            sub = a.X[ix, :]
            if hasattr(sub, "toarray"):
                sub = sub.toarray()
            me = pd.Series(sub.mean(axis=0), index=var)
            rows.append({"sample": s, "n": int(len(ix)),
                         "DM1": float(me[dm1_g].mean()),
                         "HLA_I": float(me[h1_g].mean()),
                         "HLA_II": float(me[h2_g].mean())})
        return pd.DataFrame(rows)

    immune_set = {"T_cell", "B_cell", "Plasma", "Myeloid", "DC", "NK", "Mast"}
    stromal_set = {"Fibroblast", "Endothelial"}
    masks = {
        "bulk_all_cells": np.ones(len(df), bool),
        "thyrocyte_only": (df["celltype"] == "Thyrocyte").values,
        "immune_only": df["celltype"].isin(immune_set).values,
        "stromal_only": df["celltype"].isin(stromal_set).values,
    }
    rows = []
    for name, m in masks.items():
        pb = _pb(m)
        if len(pb) < 4:
            rows.append({"subset": name, "n_samples": len(pb),
                         "rho_DM1_HLA_I": None, "rho_DM1_HLA_II": None})
            continue
        rho1, p1 = stats.spearmanr(pb["DM1"], pb["HLA_I"]) if pb["HLA_I"].std() > 1e-6 else (np.nan, np.nan)
        rho2, p2 = stats.spearmanr(pb["DM1"], pb["HLA_II"]) if pb["HLA_II"].std() > 1e-6 else (np.nan, np.nan)
        rows.append({"subset": name, "n_samples": int(len(pb)),
                     "rho_DM1_HLA_I": float(rho1) if rho1 == rho1 else None,
                     "p_DM1_HLA_I": float(p1) if p1 == p1 else None,
                     "rho_DM1_HLA_II": float(rho2) if rho2 == rho2 else None,
                     "p_DM1_HLA_II": float(p2) if p2 == p2 else None})
    decomp = pd.DataFrame(rows)
    decomp.to_csv(TAB / f"T6_{ds_short}_decomposition.csv", index=False)
    # bar plot
    fig, ax = plt.subplots(figsize=(8, 4))
    width = 0.35
    xs = np.arange(len(decomp))
    ax.bar(xs - width / 2, [r if r is not None else 0 for r in decomp["rho_DM1_HLA_I"]],
           width, label="DM1 x HLA-I", color="#cc4444")
    ax.bar(xs + width / 2, [r if r is not None else 0 for r in decomp["rho_DM1_HLA_II"]],
           width, label="DM1 x HLA-II", color="#4477aa")
    ax.set_xticks(xs)
    ax.set_xticklabels(decomp["subset"], rotation=20)
    ax.axhline(0.318, ls="--", color="#cc4444", lw=0.7, alpha=0.5,
               label="bulk T5 ref rho HLA-I=0.318")
    ax.axhline(0.393, ls="--", color="#4477aa", lw=0.7, alpha=0.5,
               label="bulk T5 ref rho HLA-II=0.393")
    ax.set_ylabel("Spearman rho (sample-level)")
    ax.set_title(f"{label} — decomposition vs bulk T5 (TCGA n=527)\n"
                 "HLA gene-expression module — not allele genotype.", fontsize=9)
    ax.legend(fontsize=7, loc="best")
    plt.tight_layout()
    plt.savefig(FIG / f"F6_{ds_short}_decomposition.png", dpi=180, bbox_inches="tight")
    plt.close()
    return decomp


# ---------- driver: run on Pu and Lu cohorts ----------
def run_scrna_dataset(path, label, ds_short, cond_col=None):
    print(f"[{label}] loading {path}")
    a = ad.read_h5ad(path)
    a = _ensure_symbol_index(a)
    print(f"  shape {a.shape}, sample col '{[c for c in a.obs.columns if c=='sample']}'")
    info = score_modules(a, label)
    a = assign_celltype(a)
    ct_counts = a.obs["celltype"].value_counts().to_dict()
    print(f"  celltype counts: {ct_counts}")

    dotplot_celltype_genes(a, label, ds_short)
    tbl3 = violin_module_per_celltype(a, label, ds_short)
    cond_summary = None
    if cond_col is not None:
        cond_summary = condition_stratified(a, label, ds_short, cond_col)
    pb_corr = thyrocyte_pseudobulk_dm1_corr(a, label, ds_short)
    decomp = decomposition_table(a, label, ds_short)
    return {
        "label": label,
        "ds_short": ds_short,
        "n_cells": int(a.shape[0]),
        "n_genes": int(a.shape[1]),
        "n_samples": int(a.obs["sample"].nunique()),
        "celltype_counts": {k: int(v) for k, v in ct_counts.items()},
        "module_info": info,
        "thyrocyte_pseudobulk_corr": pb_corr,
        "decomposition": decomp.to_dict(orient="records") if decomp is not None else None,
        "condition_present": bool(cond_col),
    }


def run_visium_spatial():
    """Per-spot HLA-I/II module on Visium across stages, mapped_celltype overlay."""
    rows = []
    overlay_rows = []
    for samp_dir in sorted(VISIUM_DIR.iterdir()):
        if not samp_dir.is_dir() or not samp_dir.name.startswith("GSM"):
            continue
        gsm_sample = samp_dir.name  # e.g. GSM7980864_PTC-1
        sample_short = gsm_sample.split("_", 1)[1]  # PTC-1
        scored = samp_dir / f"{gsm_sample}.scored.h5ad"
        if not scored.exists():
            continue
        a = ad.read_h5ad(scored)
        # add HLA modules
        sc.tl.score_genes(a, gene_list=_genes_in(a, HLA_I), score_name="HLA_I_score", use_raw=False)
        sc.tl.score_genes(a, gene_list=_genes_in(a, HLA_II), score_name="HLA_II_score", use_raw=False)
        sc.tl.score_genes(a, gene_list=_genes_in(a, TLS_GENES), score_name="TLS_score", use_raw=False)
        # HLA-G individual gene
        if "HLA-G" in a.var_names:
            x = a[:, "HLA-G"].X
            if hasattr(x, "toarray"):
                x = x.toarray()
            a.obs["HLA_G_expr"] = np.asarray(x).ravel()
        else:
            a.obs["HLA_G_expr"] = 0.0
        # add mapped_celltype from pantheonos_demo if present
        pdir = PANTHEON_DIR / sample_short / "visium_with_mapped.h5ad"
        if pdir.exists():
            ap = ad.read_h5ad(pdir, backed="r")
            mp = ap.obs[["mapped_celltype"]] if "mapped_celltype" in ap.obs.columns else None
            if mp is not None:
                a.obs = a.obs.join(mp, how="left")
        else:
            a.obs["mapped_celltype"] = "Unmapped"
        # spot-level summary
        stage = a.obs["stage"].iloc[0]
        gsm = a.obs["gsm"].iloc[0]
        for ct, sub in a.obs.groupby("mapped_celltype"):
            rows.append({"sample": sample_short, "stage": str(stage), "gsm": str(gsm),
                         "celltype": str(ct), "n_spots": int(len(sub)),
                         "HLA_I_mean": float(sub["HLA_I_score"].mean()),
                         "HLA_II_mean": float(sub["HLA_II_score"].mean()),
                         "TLS_mean": float(sub["TLS_score"].mean()),
                         "HLA_G_mean": float(sub["HLA_G_expr"].mean()),
                         "DM1_like_mean": float(sub.get("DM1_like_score", pd.Series(np.nan)).mean()),
                         })
        # spot-level table for this sample (TLS overlap)
        spot = a.obs[["mapped_celltype", "HLA_I_score", "HLA_II_score",
                      "TLS_score", "HLA_G_expr", "DM1_like_score"]].copy()
        spot["sample"] = sample_short
        spot["stage"] = str(stage)
        overlay_rows.append(spot)
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "T7_visium_celltype_HLA_per_sample.csv", index=False)
    big = pd.concat(overlay_rows, axis=0)
    big.to_csv(TAB / "T8_visium_spotlevel.csv.gz", index=False, compression="gzip")

    # F7: Visium HLA module by stage and mapped_celltype (heatmap)
    pivI = df.groupby(["stage", "celltype"])["HLA_I_mean"].mean().unstack()
    pivII = df.groupby(["stage", "celltype"])["HLA_II_mean"].mean().unstack()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    sns.heatmap(pivI, annot=True, fmt=".2f", cmap="Reds", ax=axes[0])
    axes[0].set_title("HLA-I module (mean across samples)")
    sns.heatmap(pivII, annot=True, fmt=".2f", cmap="Blues", ax=axes[1])
    axes[1].set_title("HLA-II module")
    fig.suptitle("Visium GSE250521 — HLA module by stage x mapped_celltype\n"
                 "HLA gene-expression module — not allele genotype.", y=1.05)
    plt.tight_layout()
    plt.savefig(FIG / "F7_visium_HLA_stage_celltype_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close()

    # F8: TLS x HLA-II spot-level overlap, plus HLA-G distribution
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    s = big.dropna(subset=["TLS_score", "HLA_II_score"])
    sns.scatterplot(data=s.sample(min(20000, len(s)), random_state=0),
                    x="TLS_score", y="HLA_II_score", hue="stage",
                    s=4, alpha=0.4, ax=axes[0])
    rho_tls, p_tls = stats.spearmanr(s["TLS_score"], s["HLA_II_score"])
    axes[0].set_title(f"TLS x HLA-II spot-level\nrho={rho_tls:.3f}, p={p_tls:.2e}, n={len(s):,}",
                      fontsize=9)
    # HLA-G by celltype + stage
    hg = big.dropna(subset=["HLA_G_expr", "mapped_celltype"])
    if hg["HLA_G_expr"].sum() > 0:
        sns.boxplot(data=hg, x="mapped_celltype", y="HLA_G_expr", hue="stage",
                    ax=axes[1], fliersize=1, linewidth=0.5)
        axes[1].set_title("HLA-G expression by mapped_celltype x stage")
        axes[1].tick_params(axis="x", rotation=30)
    else:
        axes[1].text(0.5, 0.5, "HLA-G: no detectable expression", ha="center")
    fig.suptitle("Spatial HLA-II / TLS / HLA-G — gene-expression module only.", y=1.05)
    plt.tight_layout()
    plt.savefig(FIG / "F8_visium_TLS_HLAG.png", dpi=180, bbox_inches="tight")
    plt.close()

    # spatial pseudobulk (per Visium sample, all spots): DM1_like x HLA module, similar to bulk
    pb = (big.groupby("sample")
              .agg(DM1=("DM1_like_score", "mean"),
                   HLA_I=("HLA_I_score", "mean"),
                   HLA_II=("HLA_II_score", "mean"),
                   stage=("stage", "first"))
              .reset_index())
    pb.to_csv(TAB / "T9_visium_sample_pseudobulk.csv", index=False)
    # decomposition: thyrocyte-only Visium (mapped_celltype == "Malignant cell" approx Thyrocyte)
    rows = []
    for label_subset, mask in [
        ("bulk_all_spots", big["DM1_like_score"].notna()),
        ("malignant_only", (big["mapped_celltype"] == "Malignant cell") & big["DM1_like_score"].notna()),
        ("immune_only", (big["mapped_celltype"].isin(["T cell", "B cell", "Myeloid cell"]))
                          & big["DM1_like_score"].notna()),
        ("stromal_only", (big["mapped_celltype"] == "Fibroblast") & big["DM1_like_score"].notna()),
    ]:
        sub = big[mask]
        if sub.empty:
            continue
        ag = (sub.groupby("sample")
                  .agg(DM1=("DM1_like_score", "mean"),
                       HLA_I=("HLA_I_score", "mean"),
                       HLA_II=("HLA_II_score", "mean"))
                  .reset_index())
        if len(ag) >= 4 and ag["HLA_I"].std() > 1e-6:
            r1, p1 = stats.spearmanr(ag["DM1"], ag["HLA_I"])
        else:
            r1, p1 = np.nan, np.nan
        if len(ag) >= 4 and ag["HLA_II"].std() > 1e-6:
            r2, p2 = stats.spearmanr(ag["DM1"], ag["HLA_II"])
        else:
            r2, p2 = np.nan, np.nan
        rows.append({"subset": label_subset, "n_samples": int(len(ag)),
                     "rho_DM1_HLA_I": float(r1) if r1 == r1 else None,
                     "p_DM1_HLA_I": float(p1) if p1 == p1 else None,
                     "rho_DM1_HLA_II": float(r2) if r2 == r2 else None,
                     "p_DM1_HLA_II": float(p2) if p2 == p2 else None})
    visium_decomp = pd.DataFrame(rows)
    visium_decomp.to_csv(TAB / "T10_visium_decomposition.csv", index=False)

    return {"per_sample_celltype": df.head(40).to_dict(orient="records"),
            "n_visium_samples": int(df["sample"].nunique()),
            "tls_x_hla_ii_spot_level": {"rho": float(rho_tls), "p": float(p_tls), "n": int(len(s))},
            "decomposition": visium_decomp.to_dict(orient="records")}


def cohort_manifest():
    rows = []
    if PU_PATH.exists():
        a = ad.read_h5ad(PU_PATH, backed="r")
        rows.append({"dataset": "GSE184362 (Pu et al. PTC scRNA)", "path": str(PU_PATH),
                     "n_cells": int(a.shape[0]), "n_samples": int(a.obs["sample"].nunique()),
                     "conditions": "PTC tumor + normal/adjacent (Pu)",
                     "celltype_label": "leiden + marker score (this track)",
                     "role": "primary scRNA #1"})
    if LU_PATH.exists():
        a = ad.read_h5ad(LU_PATH, backed="r")
        rows.append({"dataset": "GSE193581 (Lu 2023 PTC+ATC scRNA)", "path": str(LU_PATH),
                     "n_cells": int(a.shape[0]), "n_samples": int(a.obs["sample"].nunique()),
                     "conditions": "PTC vs ATC; BRAF / RAS / WT",
                     "celltype_label": "leiden + marker score (this track)",
                     "role": "primary scRNA #2"})
    n_vis = sum(1 for p in VISIUM_DIR.iterdir() if p.is_dir() and p.name.startswith("GSM"))
    rows.append({"dataset": "GSE250521 (Visium spatial atlas N/PTC/LPTC/ATC)",
                 "path": str(VISIUM_DIR),
                 "n_cells": "spots, ~2k/sample", "n_samples": n_vis,
                 "conditions": "Normal, PTC, LPTC, ATC",
                 "celltype_label": "MOSCOT-mapped (pantheonos_demo)",
                 "role": "spatial overlay"})
    rows.append({"dataset": "GSE163203 (PTC + HT scRNA, bridge)", "path": "not_downloaded",
                 "n_cells": "n/a", "n_samples": "n/a", "conditions": "PTC+HT",
                 "celltype_label": "n/a",
                 "role": "BRIDGE — not loaded locally; flagged in registry as not_downloaded"})
    rows.append({"dataset": "GSE191288 (PTC scRNA)", "path": "not_loaded",
                 "n_cells": "n/a", "n_samples": "n/a", "conditions": "PTC",
                 "celltype_label": "n/a", "role": "not loaded locally"})
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "T1_cohort_manifest.csv", index=False)
    return df


def main():
    summary = {
        "track": "Track 12 - scRNA HLA-I/II module per cell type (thyroid)",
        "boundary": ("HLA gene-expression module ONLY, not allele genotype. "
                     "Per HLA_CANCER_SEPARATION_RULES Section 1.1 — Paper 1 transcriptomic immune-context. "
                     "Bridge zone hypothesis-only: GSE163203 PTC+HT was not loaded locally."),
        "bulk_T5_reference": {"DM1_x_HLA_I_rho": 0.318, "DM1_x_HLA_II_rho": 0.393, "n": 527},
    }
    cm = cohort_manifest()
    summary["cohorts"] = cm.to_dict(orient="records")

    pu_summary = run_scrna_dataset(PU_PATH, "GSE184362 (Pu PTC)", "pu", cond_col=None)
    summary["GSE184362_pu"] = pu_summary

    lu_summary = run_scrna_dataset(LU_PATH, "GSE193581 (Lu PTC+ATC)", "lu", cond_col="cancer_type")
    summary["GSE193581_lu"] = lu_summary

    visium_summary = run_visium_spatial()
    summary["GSE250521_visium"] = visium_summary

    with open(OUT / "track12_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print("done; summary written to", OUT / "track12_summary.json")
    return summary


if __name__ == "__main__":
    main()
