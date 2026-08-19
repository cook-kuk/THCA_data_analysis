#!/usr/bin/env python3
"""Figure 3 v2 — Nature-style label-anchored validation panel.

3 sub-panels: GSE151179 boxplot, GSE299988 boxplot with caveat, Mu 2024 4-class.
Restyled with shared Nature Medicine palette + typography.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import (setup_rc, panel_letter, panel_title, style_blank,
                            IVORY, INK, MUTED, BLUE, BLUE_L, RED, RED_L, GRAY_P, GRAY_L,
                            BEIGE, BEIGE_L)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures"; FIG.mkdir(parents=True, exist_ok=True)


def boxplot_panel(ax, df, label_col, order, colors, edge_colors, title, score_col="panel_z"):
    data = [df.loc[df[label_col] == g, score_col].dropna().values for g in order]
    bp = ax.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True, widths=0.5,
                    medianprops={"color": "#222", "linewidth": 1.2},
                    whiskerprops={"color": INK, "linewidth": 0.7},
                    capprops={"color": INK, "linewidth": 0.7})
    for patch, c, e in zip(bp["boxes"], colors, edge_colors):
        patch.set_facecolor(c); patch.set_alpha(0.85); patch.set_edgecolor(e); patch.set_linewidth(0.8)
    rng = np.random.default_rng(7)
    for i, vals in enumerate(data):
        x = rng.uniform(-0.10, 0.10, size=len(vals)) + (i + 1)
        ax.scatter(x, vals, s=22, alpha=0.7, c=edge_colors[i], edgecolors="white",
                   linewidths=0.5, zorder=3)
        ax.text(i + 1, 0.04, f"n = {len(vals)}", ha="center", va="bottom", fontsize=8.5,
                color=MUTED, transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#999", lw=0.5, ls="--", alpha=0.7)
    ax.set_ylabel("8-gene panel z (within-cohort)")
    panel_title(ax, title)


def panel_a_gse151179(ax):
    df = pd.read_csv(PROCESSED / "GSE151179_label_joined.tsv", sep="\t", index_col=0)
    df["group"] = "tumour · refractory"
    df.loc[df["sample_type"].astype(str).str.contains("non-neoplastic", case=False, na=False), "group"] = "normal"
    rai = df["rai_response"].astype(str).str.lower()
    df.loc[(df["group"] != "normal") & rai.eq("avid"), "group"] = "tumour · avid"
    order = ["normal", "tumour · avid", "tumour · refractory"]
    boxplot_panel(ax, df, "group", order,
                   colors=[BEIGE_L, BLUE_L, RED_L],
                   edge_colors=["#7a6038", BLUE, RED],
                   title=f"GSE151179  ·  Tier 2  ·  n = 13 normal + 4 avid + 35 refractory  ·  tumour vs normal AUC = 0.96")


def panel_b_gse299988(ax):
    df = pd.read_csv(PROCESSED / "GSE299988_label_joined.tsv", sep="\t", index_col=0)
    NORMAL = r"non-neoplastic|^normal$|adjacent"
    df["group"] = "tumour · refractive"
    df.loc[df["sample_type"].astype(str).str.contains(NORMAL, case=False, na=False, regex=True), "group"] = "normal"
    rai = df["rai_response"].astype(str).str.lower().replace({"refractive": "refractive"})
    df.loc[(df["group"] != "normal") & rai.eq("avid"), "group"] = "tumour · avid"
    order = ["normal", "tumour · avid", "tumour · refractive"]
    boxplot_panel(ax, df, "group", order,
                   colors=[BEIGE_L, BLUE_L, RED_L],
                   edge_colors=["#7a6038", BLUE, RED],
                   title="GSE299988  ·  Tier 2 supportive  ·  AUC tumour vs normal = 1.00  ·  binary direction reversed (caveat)")


def panel_c_mu2024(ax):
    style_blank(ax, xlim=(0, 100), ylim=(0, 100))
    panel_title(ax, "Mu 2024 (HRA004166)  ·  Chinese  ·  n = 214 metastatic DTC  ·  4-class RAI uptake patterns")
    classes = [
        ("I-RAIR",  "initially RAI-refractory",  80, RED,    RED_L,    "BRAF V600E 61%  ·  TERT 51%"),
        ("G-RAIR",  "gradually RAI-refractory",  19, GRAY_P, GRAY_L,   "late-hit accrual"),
        ("P-RAIR",  "partly RAI-refractory",     10, GRAY_P, GRAY_L,   "mixed-lesion avidity"),
        ("C-RAIA",  "continually RAI-avid",      48, BLUE,   BLUE_L,   "RAS ↑   late-hit 27 vs 50%"),
    ]
    total_n = sum(c[2] for c in classes)
    by = 50; bh = 16
    x0 = 6
    for code, name, n, edge, face, note in classes:
        w = 86 * (n / total_n)
        ax.add_patch(FancyBboxPatch((x0, by - bh/2), w, bh, boxstyle="round,pad=0,rounding_size=1.2",
                                    facecolor=face, edgecolor=edge, lw=1.0))
        ax.text(x0 + w/2, by + 2.5, code, ha="center", fontsize=11.5, color=edge, fontweight="bold")
        ax.text(x0 + w/2, by - 3.5, f"n = {n}", ha="center", fontsize=9, color=INK)
        ax.text(x0 + w/2, by - bh/2 - 3.5, name, ha="center", fontsize=8, color=MUTED, style="italic")
        ax.text(x0 + w/2, by - bh/2 - 7.5, note, ha="center", fontsize=8, color=INK)
        x0 += w
    # Gray-zone bracket
    g_start = 6 + 86 * (80 / total_n)
    g_end   = 6 + 86 * ((80 + 19 + 10) / total_n)
    ax.annotate("", xy=(g_end, 78), xytext=(g_start, 78),
                arrowprops=dict(arrowstyle="-", color=GRAY_P, lw=2.2))
    ax.plot([g_start, g_start], [76, 80], color=GRAY_P, lw=1.5)
    ax.plot([g_end, g_end], [76, 80], color=GRAY_P, lw=1.5)
    ax.text((g_start + g_end) / 2, 84, "molecular gray zone  (G-RAIR + P-RAIR = 29 / 214 = 14 %)",
            ha="center", fontsize=10.5, color=GRAY_P, fontweight="bold")
    ax.text(50, 20,
            "Driver mutations alone do not partition the gray-zone classes  →  expression-based panel score required.",
            ha="center", fontsize=9.5, color=INK, style="italic")


def main():
    fig = plt.figure(figsize=(16, 10), facecolor=IVORY)
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.22,
                          left=0.06, right=0.97, top=0.91, bottom=0.06,
                          height_ratios=[1.0, 0.85])
    fig.suptitle("Label-anchored validation of the 8-gene differentiation-silencing axis",
                 fontsize=14.5, fontweight="bold", color=INK, y=0.965)
    fig.text(0.5, 0.93,
             "Tier 2 RAI-avidity cohorts (GSE151179, GSE299988) plus Mu 2024 4-class gray-zone framework",
             ha="center", fontsize=10, color=MUTED, style="italic")

    ax_a = fig.add_subplot(gs[0, 0])
    panel_a_gse151179(ax_a); panel_letter(ax_a, "a")
    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_gse299988(ax_b); panel_letter(ax_b, "b")
    ax_c = fig.add_subplot(gs[1, :])
    panel_c_mu2024(ax_c); panel_letter(ax_c, "c")

    out_png = FIG / "figure3_label_anchored_validation.png"
    out_pdf = FIG / "figure3_label_anchored_validation.pdf"
    fig.savefig(out_png, dpi=210, facecolor=IVORY)
    fig.savefig(out_pdf, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
