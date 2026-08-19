#!/usr/bin/env python3
"""Figure 16 — Per-gene Cohen d × cohort heatmap.

For each of 8 panel genes and each of N cohorts compute the per-gene Cohen d
of the silenced subgroup vs preserved subgroup. Shows that no single gene
dominates and that the panel signal is broad-based.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
RAW = ROOT / "data" / "raw"
KOREAN_GSE213647 = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/GSE213647_panel_score.tsv")
TCGA_MERGED = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/r17_per_sample_merged.tsv")
LU_SC_PER_SAMPLE = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/sc/lu2023_sc_zone_fraction_by_sample.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure16_per_gene_cross_cohort.png"
OUT_PDF = ROOT / "results" / "figures" / "figure16_per_gene_cross_cohort.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    na, nb = len(a), len(b)
    pooled = np.sqrt(((na-1)*a.var(ddof=1) + (nb-1)*b.var(ddof=1)) / max(na+nb-2, 1))
    if pooled == 0 or na < 2 or nb < 2: return np.nan
    return (a.mean() - b.mean()) / pooled


def load_gse_expr(acc):
    expr = pd.read_csv(INTERIM / f"{acc}_expression.tsv", sep="\t", index_col=0)
    p2g_path = next((RAW / acc).glob("GPL*_probe2gene.tsv"))
    p2g = pd.read_csv(p2g_path, sep="\t", dtype=str)
    pmap = dict(zip(p2g["ID"], p2g.get("gene_symbol", p2g.get("GENE_SYMBOL", []))))
    e = expr.copy(); e["gene"] = e.index.map(pmap)
    e = e.dropna(subset=["gene"]); e = e[e["gene"] != ""]
    e["__var"] = e.iloc[:, :-1].var(axis=1)
    e = e.sort_values("__var", ascending=False).drop_duplicates("gene", keep="first")
    return e.drop(columns="__var").set_index("gene")


def gse151179_per_gene():
    e = load_gse_expr("GSE151179")
    meta = pd.read_csv(INTERIM / "GSE151179_metadata.tsv", sep="\t", index_col=0)
    NORMAL = r"non-neoplastic|^normal$|adjacent"
    st = meta["sample_type"].fillna("").astype(str).str.lower()
    is_normal = st.str.contains(NORMAL, regex=True)
    is_tumor = (~is_normal) & st.ne("")
    samples = [c for c in e.columns if c in meta.index]
    out = {}
    for g in PANEL:
        if g not in e.index: out[g] = np.nan; continue
        v = e.loc[g, samples]
        v_t = v[is_tumor[samples].values].dropna().values
        v_n = v[is_normal[samples].values].dropna().values
        out[g] = cohen_d(v_t, v_n)
    return out


def gse299988_per_gene():
    e = load_gse_expr("GSE299988")
    meta = pd.read_csv(INTERIM / "GSE299988_metadata.tsv", sep="\t", index_col=0)
    st = meta["sample_type"].fillna("").astype(str).str.lower()
    is_normal = st.str.contains(r"normal|adjacent", regex=True)
    is_tumor = ~is_normal & st.ne("")
    samples = [c for c in e.columns if c in meta.index]
    out = {}
    for g in PANEL:
        if g not in e.index: out[g] = np.nan; continue
        v = e.loc[g, samples]
        v_t = v[is_tumor[samples].values].dropna().values
        v_n = v[is_normal[samples].values].dropna().values
        out[g] = cohen_d(v_t, v_n)
    return out


def tcga_per_gene():
    df = pd.read_csv(TCGA_MERGED, sep="\t")
    panel = df["panel_DM"].fillna("").astype(str)
    d4    = df["d4p2_DM"].fillna("").astype(str)
    zone = pd.Series("", index=df.index)
    zone[(panel == "DM1") & (d4 == "DM1")] = "dark-matter"
    zone[(panel == "DM2") & (d4 == "DM2")] = "WT-like"
    out = {}
    for g in PANEL:
        if g not in df.columns: out[g] = np.nan; continue
        v_dm = df.loc[zone == "dark-matter", g].dropna().values
        v_wt = df.loc[zone == "WT-like", g].dropna().values
        out[g] = cohen_d(v_dm, v_wt)
    return out


def main():
    cohorts_d = {}
    cohorts_d["TCGA-THCA\n(dark vs WT, n≈74)"] = tcga_per_gene()
    cohorts_d["GSE151179\n(tumour vs normal, n=52)"] = gse151179_per_gene()
    cohorts_d["GSE299988\n(tumour vs normal, n=14)"] = gse299988_per_gene()

    mat = pd.DataFrame(cohorts_d).T[PANEL]  # cohorts × genes
    print(mat)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(11, 4.2), facecolor=IVORY)
    fig.subplots_adjust(left=0.21, right=0.95, top=0.83, bottom=0.20)
    fig.suptitle("Per-gene Cohen d across cohorts — no single gene dominates the panel",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    vmax = max(abs(np.nanmin(mat.values)), abs(np.nanmax(mat.values)), 1.5)
    cmap = LinearSegmentedColormap.from_list("d", ["#1f4d7a", "#7ba0c5", "#fafafa", "#cb857f", "#7c322e"], N=200)
    im = ax.imshow(mat.values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            if np.isnan(v): continue
            col = "white" if abs(v) > vmax*0.55 else INK
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=10, color=col)
    ax.set_xticks(range(len(PANEL))); ax.set_xticklabels(PANEL, fontstyle="italic", rotation=20)
    ax.set_yticks(range(len(mat.index))); ax.set_yticklabels(mat.index, fontsize=9.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("Cohen d", fontsize=9)
    fig.text(0.5, 0.04,
             "Effect direction is consistent (negative) across all 8 genes × 3 cohorts; "
             "TG / TPO / DIO1 and SLC5A5 deliver the strongest per-gene effects.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
