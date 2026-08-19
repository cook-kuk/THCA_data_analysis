#!/usr/bin/env python3
"""Figure 36 — Pooled R17 zone composition across cohorts.

Side-by-side stacked bars showing zone composition across TCGA-THCA, Lu 2023
sc-RNA (PTC and ATC pools), Mun 2025 proteome (PTC, PDTC, ATC), GSE286332
Korean PTC. Demonstrates dark-matter zone is consistently present.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
TCGA = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
LU_HIST = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/sc/lu2023_sc_zone_fraction_by_histology.tsv")
MUN = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/proteome_mun2025/mun2025_zone_fraction_by_group.tsv")
GSE286 = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/cross_ethnic/gse286332_zone_fraction_by_group.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure36_pooled_zone_composition.png"
OUT_PDF = ROOT / "results" / "figures" / "figure36_pooled_zone_composition.pdf"

ZONE_ORDER = ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]
ZONE_COLOR = {"WT-like": "#cdc6c0", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}
ZONE_ALIASES = {
    "wild-type-like (preserved, no HT)": "WT-like",
    "WT-like (preserved,no HT)": "WT-like",
    "RAS-like zone (preserved + HT)": "RAS-like",
    "RAS-like (preserved+HT)": "RAS-like",
    "BRAF-like zone (silenced, no HT)": "BRAF-like",
    "dark-matter zone (silenced + HT)": "dark-matter",
    "dark-matter (silenced+HT)": "dark-matter",
}


def norm_zone(c):
    return ZONE_ALIASES.get(c, c)


def main():
    rows = []

    # TCGA
    tcga = pd.read_csv(TCGA, sep="\t").dropna(subset=["zone"])
    fr = tcga["zone"].value_counts(normalize=True).reindex(ZONE_ORDER).fillna(0)
    rows.append(("TCGA-THCA bulk\n(n = " + f"{len(tcga):,})", fr.values))

    # Lu 2023 sc PTC
    lu = pd.read_csv(LU_HIST, sep="\t").set_index("histology")
    lu.columns = [norm_zone(c) for c in lu.columns]
    if lu.values.max() > 1.5: lu = lu / 100.0
    for hist in ["PTC", "ATC"]:
        if hist in lu.index:
            fr = lu.loc[hist].reindex(ZONE_ORDER).fillna(0)
            rows.append((f"Lu 2023 sc-RNA {hist}", fr.values))

    # Mun 2025 proteome
    mun = pd.read_csv(MUN, sep="\t")
    mun = mun.set_index(mun.columns[0])
    if mun.values.max() > 1.5: mun = mun / 100.0
    for grp in ["PTC", "PDTC", "ATC"]:
        if grp in mun.index:
            fr = mun.loc[grp].reindex(ZONE_ORDER).fillna(0)
            rows.append((f"Mun 2025 proteome {grp}", fr.values))

    # GSE286332 Korean
    gse286 = pd.read_csv(GSE286, sep="\t")
    if "zone" in gse286.columns and "fraction" in gse286.columns:
        for grp, sub in gse286.groupby(gse286.columns[0]):
            fr = sub.set_index("zone")["fraction"].reindex(ZONE_ORDER).fillna(0)
            if fr.values.max() > 1.5: fr = fr / 100.0
            rows.append((f"GSE286332 {grp}", fr.values))
    else:
        gse286 = gse286.set_index(gse286.columns[0])
        gse286.columns = [norm_zone(c) for c in gse286.columns]
        if gse286.values.max() > 1.5: gse286 = gse286 / 100.0
        for grp in gse286.index:
            fr = gse286.loc[grp].reindex(ZONE_ORDER).fillna(0)
            rows.append((f"GSE286332 {grp}", fr.values))

    print(f"# cohort/group rows = {len(rows)}")

    fig, ax = plt.subplots(figsize=(15, 7), facecolor=IVORY)
    fig.subplots_adjust(left=0.05, right=0.97, top=0.88, bottom=0.32)
    fig.suptitle("Pooled R17 zone composition across cohorts and histologies",
                 fontsize=14, fontweight="bold", color=INK, y=0.96)

    labels = [r[0] for r in rows]
    vals = np.array([r[1] for r in rows])
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for j, z in enumerate(ZONE_ORDER):
        ax.bar(x, vals[:, j], bottom=bottom, color=ZONE_COLOR[z], label=z, width=0.65,
               edgecolor="white", lw=0.5)
        for xi, v, b in zip(x, vals[:, j], bottom):
            if v > 0.06:
                ax.text(xi, b + v/2, f"{v*100:.0f}%", ha="center", va="center",
                        fontsize=9, color="white", fontweight="bold")
        bottom += vals[:, j]

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9, rotation=30, ha="right")
    ax.set_ylabel("R17 zone fraction")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=4, frameon=False)

    fig.text(0.5, 0.05,
             "Dark-matter and BRAF-like zones consistently appear across bulk · single-cell · proteome · Korean cohorts — "
             "the silencing axis is platform-independent.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
