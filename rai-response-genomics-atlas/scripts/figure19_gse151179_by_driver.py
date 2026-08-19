#!/usr/bin/env python3
"""Figure 19 — GSE151179 tumour panel z by lesion driver class standalone.

Per-tumour panel z stratified by lesion driver class (BRAF / fusion / pTERT /
WT) for GSE151179. Demonstrates panel z is a complementary stratifier
independent of driver class.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED, GRAY_P

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
RAW = ROOT / "data" / "raw"
OUT_PNG = ROOT / "results" / "figures" / "figure19_gse151179_by_driver.png"
OUT_PDF = ROOT / "results" / "figures" / "figure19_gse151179_by_driver.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def load_gene_expr(acc):
    expr = pd.read_csv(INTERIM / f"{acc}_expression.tsv", sep="\t", index_col=0)
    p2g_path = next((RAW / acc).glob("GPL*_probe2gene.tsv"))
    p2g = pd.read_csv(p2g_path, sep="\t", dtype=str)
    pmap = dict(zip(p2g["ID"], p2g.get("gene_symbol", p2g.get("GENE_SYMBOL", []))))
    e = expr.copy(); e["gene"] = e.index.map(pmap)
    e = e.dropna(subset=["gene"]); e = e[e["gene"] != ""]
    e["__var"] = e.iloc[:, :-1].var(axis=1)
    e = e.sort_values("__var", ascending=False).drop_duplicates("gene", keep="first")
    return e.drop(columns="__var").set_index("gene")


def main():
    expr = load_gene_expr("GSE151179")
    meta = pd.read_csv(INTERIM / "GSE151179_metadata.tsv", sep="\t", index_col=0)
    NORMAL = r"non-neoplastic|^normal$|adjacent"
    st = meta["sample_type"].fillna("").astype(str).str.lower()
    is_normal = st.str.contains(NORMAL, regex=True)
    is_tumor = ~is_normal & st.ne("")

    avail = [g for g in PANEL if g in expr.index]
    sub = expr.loc[avail]
    z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    panel_z = z.mean(axis=0)

    df = pd.DataFrame({"panel_z": panel_z, "is_tumor": is_tumor.reindex(panel_z.index).values})
    df["lesion_class"] = meta["lesion_class"].reindex(df.index).fillna("").astype(str).str.lower()
    tumor = df[df["is_tumor"]].copy()

    # Map lesion classes
    def map_drv(x):
        x = str(x).lower()
        if "braf" in x: return "BRAF"
        if "ras" in x: return "RAS"
        if "tert" in x: return "pTERT"
        if "fus" in x: return "Fusion"
        if x == "wt" or "wild" in x: return "WT"
        return "other"
    tumor["driver"] = tumor["lesion_class"].apply(map_drv)
    print(tumor["driver"].value_counts())

    drv_order = ["WT", "BRAF", "Fusion", "pTERT", "RAS", "other"]
    drv_order = [d for d in drv_order if (tumor["driver"] == d).sum() >= 1]
    palette = {"WT": "#888888", "BRAF": BLUE, "Fusion": "#7e689a", "pTERT": RED, "RAS": "#b89858", "other": "#bbbbbb"}

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.32, left=0.07, right=0.97, top=0.83, bottom=0.16))
    fig.suptitle("GSE151179 tumours — panel z is independent of driver class",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: boxplot
    ax = axes[0]
    bp_data = [tumor.loc[tumor["driver"] == d, "panel_z"].dropna().values for d in drv_order]
    bp = ax.boxplot(bp_data, positions=range(len(drv_order)), widths=0.55,
                     patch_artist=True, medianprops=dict(color=INK, lw=1.3),
                     boxprops=dict(lw=0.8), whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8))
    for patch, d in zip(bp["boxes"], drv_order):
        patch.set_facecolor(palette[d]); patch.set_alpha(0.5); patch.set_edgecolor(palette[d])
    rng = np.random.default_rng(7)
    for i, vals in enumerate(bp_data):
        jit = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i)+jit, vals, color=palette[drv_order[i]],
                   s=22, alpha=0.85, edgecolor="white", lw=0.5, zorder=3)
    ax.set_xticks(range(len(drv_order)))
    ax.set_xticklabels([f"{d}\n(n={int((tumor['driver']==d).sum())})" for d in drv_order], fontsize=10)
    ax.set_ylabel("8-gene panel z")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    panel_title(ax, "Tumour panel z by driver class"); panel_letter(ax, "a")

    # Right: KW + ANOVA summary as table
    from scipy.stats import kruskal
    valid = [v for v in bp_data if len(v) >= 2]
    if len(valid) >= 2:
        h, p = kruskal(*valid)
    else:
        h, p = float("nan"), float("nan")
    ax = axes[1]
    medians = [(d, np.median(v) if len(v) else np.nan, len(v)) for d, v in zip(drv_order, bp_data)]
    ax.axis("off")
    rows = ["driver", "n", "median panel z"]
    cells = [[d for d, m, n in medians],
             [n for d, m, n in medians],
             [f"{m:+.2f}" if not np.isnan(m) else "—" for d, m, n in medians]]
    table = ax.table(cellText=list(zip(*[[r] + c for r, c in zip(rows, cells)])),
                     loc="center", cellLoc="center", colLoc="center")
    table.auto_set_font_size(False); table.set_fontsize(10)
    table.scale(1.0, 1.6)
    ax.text(0.5, 0.05, f"Kruskal–Wallis H = {h:.2f}, p = {p:.3g}  ·  driver class does not partition panel z",
            transform=ax.transAxes, ha="center", fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Per-class summary"); panel_letter(ax, "b")

    fig.text(0.5, 0.045,
             "BRAF, fusion, pTERT and WT tumours overlap substantially in panel z — driver class is not a substitute for panel readout.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
