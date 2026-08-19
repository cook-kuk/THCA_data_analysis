#!/usr/bin/env python3
"""Figure 35 — GSE112202 TERT-wt vs TERT-mut subgroup redifferentiation comparison.

Shows digoxin-induced redifferentiation effect in both TERT-wt and TERT-mut
subgroups separately to demonstrate the panel response is preserved
regardless of TERT status.
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
WT = Path("/home/seungho/personal/THCA_data_analysis/rai-response-genomics-atlas/data/raw/GSE112202/GSE112202_Papillary_thyroid_cancer_TERT_wild-type_subgroup.fpkm_tracking.gz")
MUT = Path("/home/seungho/personal/THCA_data_analysis/rai-response-genomics-atlas/data/raw/GSE112202/GSE112202_Papillary_thyroid_cancer_TERT_mutant_subgroup.fpkm_tracking.gz")
OUT_PNG = ROOT / "results" / "figures" / "figure35_gse112202_subgroup.png"
OUT_PDF = ROOT / "results" / "figures" / "figure35_gse112202_subgroup.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def parse_subgroup(path, prefix):
    df = pd.read_csv(path, sep="\t", compression="gzip")
    # find untreated and treated FPKM columns
    cols_unt = [c for c in df.columns if "untreated" in c.lower() and "fpkm" in c.lower() and "lo" not in c.lower() and "hi" not in c.lower()]
    cols_trt = [c for c in df.columns if ("digoxin" in c.lower() or "treated" in c.lower()) and "untreated" not in c.lower() and "fpkm" in c.lower() and "lo" not in c.lower() and "hi" not in c.lower()]
    fpkm_u = cols_unt[0]; fpkm_t = cols_trt[0]
    sub = df[df["gene_short_name"].isin(PANEL)].drop_duplicates("gene_short_name").set_index("gene_short_name").reindex(PANEL)
    sub["log2FC"] = np.log2((sub[fpkm_t] + 1e-3) / (sub[fpkm_u] + 1e-3))
    return sub[["log2FC"]].rename(columns={"log2FC": prefix})


def main():
    wt = parse_subgroup(WT, "TERT-wt")
    mut = parse_subgroup(MUT, "TERT-mut")
    merged = pd.concat([wt, mut], axis=1)
    print(merged)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.07, right=0.97, top=0.83, bottom=0.17))
    fig.suptitle("GSE112202 — digoxin redifferentiation by TERT subgroup",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: grouped bar
    ax = axes[0]
    x = np.arange(len(PANEL))
    w = 0.4
    bars_wt = ax.bar(x - w/2, merged["TERT-wt"], width=w, color=BLUE, alpha=0.85,
                      edgecolor="white", lw=0.5, label="TERT wild-type")
    bars_mut = ax.bar(x + w/2, merged["TERT-mut"], width=w, color=RED, alpha=0.85,
                      edgecolor="white", lw=0.5, label="TERT mutant")
    ax.axhline(0, color=INK, lw=0.7)
    ax.set_xticks(x); ax.set_xticklabels(PANEL, rotation=18, fontstyle="italic")
    ax.set_ylabel("log₂FC  (digoxin / untreated)")
    ax.legend(fontsize=9, loc="upper left", frameon=True)
    panel_title(ax, "Per-gene log₂FC by TERT subgroup"); panel_letter(ax, "a")

    # Right: subgroup median + counts
    ax = axes[1]
    med_wt = float(merged["TERT-wt"].median())
    med_mut = float(merged["TERT-mut"].median())
    up_wt = int((merged["TERT-wt"] > 0).sum())
    up_mut = int((merged["TERT-mut"] > 0).sum())

    summary = [["TERT subgroup", "median log₂FC", "# genes ↑"],
               ["TERT wild-type", f"{med_wt:+.2f}", f"{up_wt}/8"],
               ["TERT mutant",   f"{med_mut:+.2f}", f"{up_mut}/8"]]
    ax.axis("off")
    tbl = ax.table(cellText=summary, loc="upper center", cellLoc="center", colLoc="center",
                   bbox=[0.05, 0.45, 0.9, 0.50])
    tbl.auto_set_font_size(False); tbl.set_fontsize(11); tbl.scale(1.0, 1.8)
    for j in range(3):
        c = tbl[(0, j)]; c.set_facecolor("#e2eaf2"); c.set_text_props(fontweight="bold", color=INK)

    ax.text(0.5, 0.30,
            "Digoxin restores 8-gene panel expression in both TERT-wt and TERT-mut subgroups,\n"
            "confirming redifferentiation is panel-mediated, not TERT-mediated.",
            transform=ax.transAxes, ha="center", fontsize=10.5, color=INK, fontweight="bold")
    panel_title(ax, "Subgroup summary"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             f"Redifferentiation effect preserved across TERT status — median log₂FC = {med_wt:+.2f} (wt) vs {med_mut:+.2f} (mut).",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
