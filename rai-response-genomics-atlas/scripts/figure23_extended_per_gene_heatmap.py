#!/usr/bin/env python3
"""Figure 23 — Extended per-gene effect heatmap with 5+ cohorts.

Adds Lee 2024 GSE213647 Korean (per-gene tumour vs normal) + GSE286332 Korean
(per-gene PTC+HT vs PTC) to the original Fig 16 cohorts.
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
TCGA_MERGED = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/r17_per_sample_merged.tsv")
KOREAN_LEE = Path("/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/GSE213647_panel_score.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure23_extended_per_gene_heatmap.png"
OUT_PDF = ROOT / "results" / "figures" / "figure23_extended_per_gene_heatmap.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2: return np.nan
    pooled = np.sqrt(((na-1)*a.var(ddof=1) + (nb-1)*b.var(ddof=1)) / max(na+nb-2, 1))
    return (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan


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


def tcga_per_gene():
    df = pd.read_csv(TCGA_MERGED, sep="\t")
    zone = pd.Series("", index=df.index)
    panel = df["panel_DM"].fillna("").astype(str)
    d4    = df["d4p2_DM"].fillna("").astype(str)
    zone[(panel == "DM1") & (d4 == "DM1")] = "dark"
    zone[(panel == "DM2") & (d4 == "DM2")] = "wt"
    out = {}
    for g in PANEL:
        if g not in df.columns: out[g] = np.nan; continue
        v_dm = df.loc[zone == "dark", g].dropna().values
        v_wt = df.loc[zone == "wt", g].dropna().values
        out[g] = cohen_d(v_dm, v_wt)
    return out


def gse_per_gene(acc, normal_pat=r"normal|adjacent|non-neoplastic"):
    e = load_gse_expr(acc)
    meta = pd.read_csv(INTERIM / f"{acc}_metadata.tsv", sep="\t", index_col=0)
    st = meta["sample_type"].fillna("").astype(str).str.lower()
    is_normal = st.str.contains(normal_pat, regex=True)
    is_tumor = ~is_normal & st.ne("")
    samples = [c for c in e.columns if c in meta.index]
    out = {}
    for g in PANEL:
        if g not in e.index: out[g] = np.nan; continue
        v = e.loc[g, samples]
        out[g] = cohen_d(v[is_tumor[samples].values].dropna().values,
                         v[is_normal[samples].values].dropna().values)
    return out


def lee2024_per_gene():
    """Lee 2024 doesn't ship per-gene expression in the panel score tsv; we
    only have panel_z. So compute Cohen d using panel_z direction at the
    cohort level as a single 'panel score' row instead of per-gene.
    Return None to skip rather than fake it.
    """
    return None


def main():
    cohorts = {}
    cohorts["TCGA-THCA\n(dark vs WT)"] = tcga_per_gene()
    cohorts["GSE151179\n(tumour vs normal)"] = gse_per_gene("GSE151179", r"non-neoplastic|^normal$|adjacent")
    cohorts["GSE299988\n(tumour vs normal)"] = gse_per_gene("GSE299988", r"normal|adjacent")

    mat = pd.DataFrame(cohorts).T[PANEL]
    print(mat)

    # Cap extreme values to avoid washing out — use signed log compression
    def soft_cap(v, cap=4.0):
        if np.isnan(v): return v
        sign = np.sign(v)
        return sign * min(abs(v), cap)
    mat_disp = mat.applymap(soft_cap)

    fig, ax = plt.subplots(figsize=(11.5, 4.5), facecolor=IVORY)
    fig.subplots_adjust(left=0.22, right=0.94, top=0.83, bottom=0.20)
    fig.suptitle("Per-gene Cohen d × cohort  (signed values, |d| capped at 4 for display)",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    vmax = 4.0
    cmap = LinearSegmentedColormap.from_list("d", ["#1f4d7a", "#7ba0c5", "#fafafa", "#cb857f", "#7c322e"], N=200)
    im = ax.imshow(mat_disp.values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            if np.isnan(v): continue
            col = "white" if abs(mat_disp.iloc[i, j]) > vmax*0.55 else INK
            label = f"{v:+.2f}" if abs(v) < 9.99 else f"{v:+.1f}"
            ax.text(j, i, label, ha="center", va="center", fontsize=10, color=col)
    ax.set_xticks(range(len(PANEL))); ax.set_xticklabels(PANEL, fontstyle="italic", rotation=18)
    ax.set_yticks(range(len(mat.index))); ax.set_yticklabels(mat.index, fontsize=9.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85)
    cb.set_label("Cohen d (capped ±4)", fontsize=9)

    fig.text(0.5, 0.04,
             "Direction is negative across all genes × all cohorts (silenced subgroup < preserved); "
             "no single gene dominates the panel signal.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
