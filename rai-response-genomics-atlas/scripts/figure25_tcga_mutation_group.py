#!/usr/bin/env python3
"""Figure 25 — TCGA-THCA panel z by mutation group.

Standalone view of panel z stratified by mutation group (BRAF V600E only,
RAS only, BRAF+RAS double, WT no driver). Demonstrates panel z is not a
linear function of driver class.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kruskal

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
TCGA = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure25_tcga_mutation_group.png"
OUT_PDF = ROOT / "results" / "figures" / "figure25_tcga_mutation_group.pdf"


def main():
    df = pd.read_csv(TCGA, sep="\t").dropna(subset=["RAI_8"])
    # Build mutation group
    def assign(row):
        b = bool(row["has_braf_v600e"])
        r = bool(row["has_ras_mut"])
        if b and r: return "BRAF + RAS"
        if b: return "BRAF V600E only"
        if r: return "RAS only"
        # Check for fusion / other drivers via 'driver' field
        drv = str(row.get("driver", "")).lower()
        if "fus" in drv or "ret" in drv or "alk" in drv: return "Fusion / other"
        return "BRAF·RAS-neg"
    df["mut_grp"] = df.apply(assign, axis=1)
    print(df["mut_grp"].value_counts())

    order = ["BRAF V600E only", "RAS only", "BRAF + RAS", "Fusion / other", "BRAF·RAS-neg"]
    order = [g for g in order if (df["mut_grp"] == g).sum() >= 1]
    palette = {
        "BRAF V600E only": BLUE,
        "RAS only": RED,
        "BRAF + RAS": "#7e689a",
        "Fusion / other": "#b89858",
        "BRAF·RAS-neg": "#888888",
    }

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.07, right=0.97, top=0.83, bottom=0.18))
    fig.suptitle("TCGA-THCA panel z by mutation group (n ≈ 500)",
                 fontsize=13, fontweight="bold", color=INK, y=0.97)

    # Left: boxplot
    ax = axes[0]
    data = [df.loc[df["mut_grp"] == g, "RAI_8"].dropna().values for g in order]
    bp = ax.boxplot(data, positions=range(len(order)), widths=0.55,
                     patch_artist=True, medianprops=dict(color=INK, lw=1.3),
                     boxprops=dict(lw=0.7), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                     flierprops=dict(marker="o", markersize=2, alpha=0.4))
    for patch, g in zip(bp["boxes"], order):
        patch.set_facecolor(palette[g]); patch.set_alpha(0.55); patch.set_edgecolor(palette[g])
    rng = np.random.default_rng(13)
    for i, vals in enumerate(data):
        jit = rng.uniform(-0.10, 0.10, len(vals))
        ax.scatter(np.full(len(vals), i)+jit, vals, color=palette[order[i]], s=10,
                   alpha=0.55, edgecolor="none", zorder=2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{g}\n(n={int((df['mut_grp']==g).sum())})" for g in order],
                       fontsize=9, rotation=12)
    ax.set_ylabel("panel z (RAI_8)")
    ax.axhline(0, color=MUTED, ls="--", lw=0.7)
    valid = [v for v in data if len(v) >= 2]
    h, p = kruskal(*valid) if len(valid) >= 2 else (float("nan"), float("nan"))
    ax.text(0.02, 0.04, f"Kruskal–Wallis H = {h:.1f}, p = {p:.3g}",
            transform=ax.transAxes, fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Panel z by mutation group"); panel_letter(ax, "a")

    # Right: median + IQR summary table
    ax = axes[1]
    ax.axis("off")
    rows = []
    for g, vals in zip(order, data):
        if len(vals) == 0:
            rows.append([g, "—", "—", "—"]); continue
        med = np.median(vals); q1, q3 = np.percentile(vals, [25, 75])
        rows.append([g, f"{len(vals)}", f"{med:+.2f}", f"[{q1:+.2f}, {q3:+.2f}]"])
    table_data = [["mutation group", "n", "median", "IQR"]] + rows
    tbl = ax.table(cellText=table_data, loc="center", cellLoc="center", colLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1.0, 1.6)
    # color header
    for j in range(4):
        c = tbl[(0, j)]; c.set_facecolor("#e2eaf2"); c.set_text_props(fontweight="bold", color=INK)
    panel_title(ax, "Per-group panel z summary"); panel_letter(ax, "b")

    fig.text(0.5, 0.04,
             "Driver mutation modifies but does not determine panel z — BRAF V600E tumours occupy both silenced and preserved ranges.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
