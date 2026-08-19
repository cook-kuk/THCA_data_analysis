#!/usr/bin/env python3
"""Figure 3 — Label-anchored validation: GSE151179 + GSE299988 + Mu 2024 4-class summary.

Combines the two completed Tier-2 validations plus a pictorial summary of Mu 2024's
4-class uptake patterns (driver frequencies fetched from PMC11031230, 2026-05-21).

Output: results/figures/figure3_label_anchored_validation.{png,pdf}
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
import numpy as np
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED = ROOT / "data" / "processed"

INK = "#222"
BLUE = "#3a6ea8"; BLUE_L = "#a9c2dc"
RED  = "#b04a44"; RED_L  = "#d8a8a4"
GRAY_P = "#8e7aa5"; GRAY_L = "#c4baca"
BEIGE = "#d8c8a8"
MUTED = "#555"


def panel_a_boxplot_gse151179(ax):
    df = pd.read_csv(PROCESSED / "GSE151179_label_joined.tsv", sep="\t", index_col=0)
    # Two groups: tumor (refractory+avid+anything) vs non-neoplastic
    df["group"] = "tumor"
    df.loc[df["sample_type"].astype(str).str.contains("non-neoplastic", case=False, na=False), "group"] = "normal"
    # Within tumor, split avid vs refractory
    rai = df["rai_response"].astype(str).str.lower()
    df["plot_label"] = df["group"]
    df.loc[(df["group"] == "tumor") & rai.eq("avid"), "plot_label"] = "tumor · avid"
    df.loc[(df["group"] == "tumor") & rai.eq("refractory"), "plot_label"] = "tumor · refractory"
    order = ["normal", "tumor · avid", "tumor · refractory"]
    data = [df.loc[df["plot_label"] == g, "panel_z"].dropna().values for g in order]
    colors = [BEIGE, BLUE_L, RED_L]
    edge_colors = ["#7a6038", BLUE, RED]
    bp = ax.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True, widths=0.55)
    for patch, c, e in zip(bp["boxes"], colors, edge_colors):
        patch.set_facecolor(c); patch.set_alpha(0.85); patch.set_edgecolor(e)
    for i, vals in enumerate(data):
        x = np.random.uniform(-0.12, 0.12, size=len(vals)) + (i + 1)
        ax.scatter(x, vals, s=22, alpha=0.8, c=edge_colors[i], edgecolors="white", linewidths=0.6)
        ax.text(i + 1, 0.04, f"n = {len(vals)}", ha="center", va="bottom", fontsize=8.5,
                color="#444", transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#444", lw=0.6, ls="--")
    ax.set_ylabel("8-gene panel z (within-cohort)", fontsize=9.5)
    ax.set_title("GSE151179 · Tier 2  ·  n = 13 normal · 4 avid · 35 refractory  ·  tumor vs normal AUC = 0.96", fontsize=10, color=INK)
    ax.tick_params(axis='x', labelsize=9)


def panel_b_boxplot_gse299988(ax):
    df = pd.read_csv(PROCESSED / "GSE299988_label_joined.tsv", sep="\t", index_col=0)
    df["group"] = "tumor"
    NORMAL = r"non-neoplastic|^normal$|adjacent"
    df.loc[df["sample_type"].astype(str).str.contains(NORMAL, case=False, na=False, regex=True), "group"] = "normal"
    rai = df["rai_response"].astype(str).str.lower().replace({"refractive": "refractory"})
    df["plot_label"] = df["group"]
    df.loc[(df["group"] == "tumor") & rai.eq("avid"), "plot_label"] = "tumor · avid"
    df.loc[(df["group"] == "tumor") & rai.eq("refractory"), "plot_label"] = "tumor · refractory"
    order = ["normal", "tumor · avid", "tumor · refractory"]
    data = [df.loc[df["plot_label"] == g, "panel_z"].dropna().values for g in order]
    colors = [BEIGE, BLUE_L, RED_L]; edge_colors = ["#7a6038", BLUE, RED]
    bp = ax.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True, widths=0.55)
    for patch, c, e in zip(bp["boxes"], colors, edge_colors):
        patch.set_facecolor(c); patch.set_alpha(0.85); patch.set_edgecolor(e)
    for i, vals in enumerate(data):
        x = np.random.uniform(-0.12, 0.12, size=len(vals)) + (i + 1)
        ax.scatter(x, vals, s=30, alpha=0.85, c=edge_colors[i], edgecolors="white", linewidths=0.7)
        ax.text(i + 1, 0.04, f"n = {len(vals)}", ha="center", va="bottom", fontsize=8.5,
                color="#444", transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#444", lw=0.6, ls="--")
    ax.set_ylabel("8-gene panel z (within-cohort)", fontsize=9.5)
    ax.set_title("GSE299988 · Tier 2 supportive  ·  n = 4 normal · 5 avid · 5 refractive  ·  tumor vs normal AUC = 1.00\n"
                 "Direction reversed in 5 vs 5 binary (caveat — LN+/LN− selection confound likely)",
                 fontsize=9.8, color=INK)
    ax.tick_params(axis='x', labelsize=9)


def panel_c_mu2024_summary(ax):
    """Pictorial summary of Mu 2024 4-class uptake patterns + driver frequencies."""
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])

    ax.text(5.0, 9.55, "Mu 2024 JCEM · HRA004166 · n = 214 metastatic DTC · 4-class uptake patterns",
            ha="center", fontsize=10, color=INK, fontweight="bold")

    classes = [
        ("I-RAIR",  "initially RAI-refractory",  80, RED,    RED_L,    "BRAF V600E 61% / TERT 51% / TP53↑"),
        ("G-RAIR",  "gradually RAI-refractory",  19, GRAY_P, GRAY_L,   "late-hit ↑ over time"),
        ("P-RAIR",  "partly RAI-refractory",     10, GRAY_P, GRAY_L,   "mixed-lesion avidity"),
        ("C-RAIA",  "continually RAI-avid",      48, BLUE,   BLUE_L,   "RAS ↑ / late-hit ↓ (27% vs 50%)"),
    ]
    by = 6.5; bh = 1.0
    total_n = sum(c[2] for c in classes)
    x0 = 0.6
    for i, (code, name, n, edge, face, note) in enumerate(classes):
        w = (8.8) * (n / total_n)
        b = FancyBboxPatch((x0, by - bh/2), w, bh, boxstyle="round,pad=0,rounding_size=0.10",
                            facecolor=face, edgecolor=edge, lw=0.9)
        ax.add_patch(b)
        ax.text(x0 + w/2, by + 0.05, code, ha="center", fontsize=10, color=edge, fontweight="bold")
        ax.text(x0 + w/2, by - 0.35, f"n = {n}", ha="center", fontsize=8.5, color=INK)
        ax.text(x0 + w/2, by - bh/2 - 0.45, name, ha="center", fontsize=8, color=MUTED)
        ax.text(x0 + w/2, by - bh/2 - 0.95, note, ha="center", fontsize=7.5, color=MUTED, style="italic")
        x0 += w

    # gray zone bracket
    gray_start = 0.6 + (8.8) * (80 / total_n)
    gray_end = 0.6 + (8.8) * ((80 + 19 + 10) / total_n)
    ax.annotate("", xy=(gray_end, by + bh/2 + 0.55), xytext=(gray_start, by + bh/2 + 0.55),
                arrowprops=dict(arrowstyle="-", color=GRAY_P, lw=2.0))
    ax.text((gray_start + gray_end) / 2, by + bh/2 + 0.95,
            "molecular gray zone (n = 29 of 214)",
            ha="center", fontsize=9.5, color=GRAY_P, fontweight="bold")

    # interpretation
    ax.text(5.0, 3.5,
            "RAI uptake is NOT binary — 4 classes with intermediate gray zone (P-RAIR + G-RAIR, 14 %).",
            ha="center", fontsize=9.5, color=INK)
    ax.text(5.0, 2.9,
            "Driver mutations alone do not partition the classes cleanly:",
            ha="center", fontsize=9, color=MUTED, style="italic")
    ax.text(5.0, 2.4,
            "BRAF + TERT enriched in I-RAIR · RAS enriched in C-RAIA · late-hit (TERT/TP53/PIK3CA) 50 % I-RAIR vs 27 % I-RAIA",
            ha="center", fontsize=8.5, color=MUTED)
    ax.text(5.0, 1.3,
            "→ This justifies the manuscript's gray-zone framing for an expression-based panel score.",
            ha="center", fontsize=9, color=INK, fontweight="bold")
    ax.text(5.0, 0.7,
            "Access: NGDC HRA004166 (controlled). Frequencies usable from supplement; per-patient raw on request.",
            ha="center", fontsize=7.5, color=MUTED, style="italic")


def main():
    fig = plt.figure(figsize=(15.5, 12), facecolor="white")
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.18,
                          left=0.05, right=0.97, top=0.92, bottom=0.05,
                          height_ratios=[1.0, 1.1])

    fig.suptitle("Label-anchored validation of the 8-gene differentiation-silencing axis",
                 fontsize=14, fontweight="bold", color="#0e2a4a", y=0.985, family="DejaVu Sans")

    ax_a = fig.add_subplot(gs[0, 0])
    panel_a_boxplot_gse151179(ax_a)
    ax_a.text(-0.07, 1.10, "a", transform=ax_a.transAxes, fontsize=15, fontweight="bold",
              family="DejaVu Sans", va="top")

    ax_b = fig.add_subplot(gs[0, 1])
    panel_b_boxplot_gse299988(ax_b)
    ax_b.text(-0.07, 1.10, "b", transform=ax_b.transAxes, fontsize=15, fontweight="bold",
              family="DejaVu Sans", va="top")

    ax_c = fig.add_subplot(gs[1, :])
    panel_c_mu2024_summary(ax_c)
    ax_c.text(0.0, 1.05, "c", transform=ax_c.transAxes, fontsize=15, fontweight="bold",
              family="DejaVu Sans", va="top")

    out_png = FIG_DIR / "figure3_label_anchored_validation.png"
    out_pdf = FIG_DIR / "figure3_label_anchored_validation.pdf"
    fig.savefig(out_png, dpi=180, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
