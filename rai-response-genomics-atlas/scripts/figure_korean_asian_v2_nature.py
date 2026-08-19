#!/usr/bin/env python3
"""Korean / Asian validation figure v2 — Nature style.

5 sub-panels: cohort inventory, Lee 2024 boxplot, K2 distribution,
Mun 2025 proteome stacked, Mu 2024 4-class.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import (setup_rc, panel_letter, panel_title, style_blank,
                            IVORY, INK, MUTED, BLUE, BLUE_L, RED, RED_L, GRAY_P, GRAY_L,
                            BEIGE, BEIGE_L, KR_GREEN, AMBER)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PAPER1 = Path("/home/seungho/personal/THCA_data_analysis/project")
LEE = PAPER1 / "results/v17_korean/Lee_unsup_8gene_2026_05_20.tsv"
K2_PATH = PAPER1 / "results/p_deconv_2026_05_08/v15_k2_panel_overlay_scores.tsv"
FIG = ROOT / "results" / "figures"

CN_AMBER = AMBER


def panel_a_inventory(ax):
    style_blank(ax, xlim=(0, 100), ylim=(0, 100))
    panel_title(ax, "Korean / Asian thyroid cohort inventory  ·  locally processed  ·  1,460 patients + 47,587 cells")
    cohorts = [
        ("K2 PRJEB11591",        "Korean",  "RNA-seq",         "260",     "Yoo 2016 SNU-GMI",  "panel z available",     KR_GREEN),
        ("GSE213647 Lee 2024",   "Korean",  "RNA-seq",         "632",     "Nat Commun 2024",   "panel z + unsup DM1/DM2",KR_GREEN),
        ("GSE286332",            "Korean",  "RNA-seq",         "18",      "Cook 2026 bioRxiv", "PTC vs PTC+HT",          KR_GREEN),
        ("Mun 2025 proteome",    "Korean",  "DIA-MS",          "336",     "Cell 2025",         "OR = 8.54  p = 2.7e-15", KR_GREEN),
        ("GSE193581 Lu 2023",    "Chinese", "scRNA",           "14,624",  "JCI 2023",          "ATC 38% dark-matter",    CN_AMBER),
        ("GSE184362 Pu 2021",    "Chinese", "scRNA",           "32,963",  "Sci Adv 2021",      "RAI-refractory distant met", CN_AMBER),
        ("HRA004166 Mu 2024",    "Chinese", "NGS",             "214",     "JCEM 2024",         "4-class · gray zone 14%",CN_AMBER),
    ]
    headers = ["Cohort", "Ethnicity", "Modality", "n", "Source", "Status"]
    col_x = [2, 22, 36, 50, 60, 78]
    for x, h in zip(col_x, headers):
        ax.text(x, 88, h, fontsize=9.5, fontweight="bold", color=INK)
    ax.plot([2, 98], [85, 85], color=MUTED, lw=0.6)
    row_h = 9.5
    for i, (name, ethn, mod, n, src, status, c) in enumerate(cohorts):
        y = 79 - i * row_h
        ax.add_patch(FancyBboxPatch((21, y - 1.8), 11, 4.2, boxstyle="round,pad=0,rounding_size=0.6",
                                     facecolor=c, edgecolor="none", alpha=0.85))
        ax.text(26.5, y + 0.2, ethn, ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
        ax.text(col_x[0], y + 0.2, name, fontsize=9, color=INK, fontweight="bold")
        ax.text(col_x[2], y + 0.2, mod, fontsize=8.5, color=INK)
        ax.text(col_x[3], y + 0.2, n, fontsize=8.5, color=INK, fontweight="bold")
        ax.text(col_x[4], y + 0.2, src, fontsize=8.5, color=MUTED, style="italic")
        ax.text(col_x[5], y + 0.2, status, fontsize=8.5, color=INK)
        if i < len(cohorts) - 1:
            ax.plot([2, 98], [y - 3.0, y - 3.0], color="#e8e0cd", lw=0.4)


def panel_b_lee2024(ax):
    df = pd.read_csv(LEE, sep="\t")
    df["hist_norm"] = df["histology"].replace({"PDFP": "PDTC", "UTC/ATC": "ATC"})
    order = ["Normal", "PTC", "PDTC", "ATC"]
    present = [h for h in order if h in df["hist_norm"].unique()]
    data = [df.loc[df["hist_norm"] == h, "panel_z"].dropna().values for h in present]
    colors = {"Normal": BEIGE_L, "PTC": BLUE_L, "PDTC": GRAY_L, "ATC": RED_L}
    edge_colors = {"Normal": "#7a6038", "PTC": BLUE, "PDTC": GRAY_P, "ATC": RED}
    bp = ax.boxplot(data, tick_labels=present, showfliers=False, patch_artist=True, widths=0.55,
                    medianprops={"color": "#222", "linewidth": 1.2},
                    whiskerprops={"color": INK, "linewidth": 0.7},
                    capprops={"color": INK, "linewidth": 0.7})
    for patch, h in zip(bp["boxes"], present):
        patch.set_facecolor(colors[h]); patch.set_alpha(0.85); patch.set_edgecolor(edge_colors[h])
        patch.set_linewidth(0.8)
    rng = np.random.default_rng(7)
    for i, vals in enumerate(data):
        x = rng.uniform(-0.18, 0.18, size=len(vals)) + (i + 1)
        ax.scatter(x, vals, s=4, alpha=0.35, c=edge_colors[present[i]], edgecolors="none")
        ax.text(i + 1, 0.04, f"n = {len(vals)}", ha="center", va="bottom", fontsize=8.5,
                color=MUTED, transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#999", lw=0.5, ls="--", alpha=0.7)
    ax.set_ylabel("8-gene panel z")
    panel_title(ax, "Lee 2024 (GSE213647)  ·  Korean  ·  n = 632  ·  Normal → ATC monotonic")


def panel_c_k2(ax):
    df = pd.read_csv(K2_PATH, sep="\t")
    df = df[df["cohort"] == "K2_PRJEB11591"]
    n_dm1 = (df["DM_call_centered"] == "DM1").sum()
    ax.hist(df["panel_z_mean"], bins=40, color=BLUE_L, edgecolor=BLUE, lw=0.6, alpha=0.85)
    ax.axvline(df["panel_z_mean"].median(), color=RED, lw=1.4, label=f"median = {df['panel_z_mean'].median():+.2f}")
    ax.axvline(0, color="#999", lw=0.5, ls="--", alpha=0.6)
    ax.set_xlabel("Panel z  (within-cohort)")
    ax.set_ylabel("count")
    ax.text(0.02, 0.95, f"Korean overdiagnosis paradigm\nDM1 (silenced) = {n_dm1}/{len(df)} = {n_dm1/len(df)*100:.1f}%",
            transform=ax.transAxes, fontsize=8.5, color=MUTED, va="top", style="italic")
    ax.legend(loc="upper right", fontsize=8.5)
    panel_title(ax, f"K2 PRJEB11591  ·  Korean  ·  n = {len(df)}")


def panel_d_mun2025(ax):
    classes = ["PTC", "PDTC", "ATC"]
    n_per = [184, 46, 113]
    fractions = {
        "PTC":  [13.0, 26.6, 14.1, 46.3],
        "PDTC": [37.0, 8.7,  23.9, 30.4],
        "ATC":  [36.3, 3.5,  58.4, 1.8],
    }
    zone_order = ["BRAF-like", "RAS-like", "dark-matter", "WT-like"]
    zone_colors = {"BRAF-like": BLUE, "RAS-like": "#d18b1f", "dark-matter": GRAY_P, "WT-like": "#888888"}
    x = np.arange(len(classes))
    bottom = np.zeros(len(classes))
    for i, z in enumerate(zone_order):
        vals = [fractions[c][i] for c in classes]
        ax.bar(x, vals, bottom=bottom, color=zone_colors[z], edgecolor="white", lw=0.5,
               label=z, width=0.55)
        bottom += vals
    for i, n in enumerate(n_per):
        ax.text(x[i], -6, f"n = {n}", ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(x); ax.set_xticklabels(classes, fontsize=10, fontweight="bold")
    ax.set_ylabel("% in zone")
    ax.set_ylim(-10, 108)
    ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.02, 1))
    panel_title(ax, "Mun 2025 proteome  ·  Korean  ·  n = 336  ·  OR = 8.54  p = 2.7×10⁻¹⁵")


def panel_e_mu2024(ax):
    classes = ["I-RAIR", "G-RAIR", "P-RAIR", "C-RAIA"]
    n_per = [80, 19, 10, 48]
    drivers = {
        "BRAF V600E":    [61, 36, 30, 22], "TERT promoter": [22, 18, 12, 5],
        "TP53":          [8, 6, 5, 3], "PIK3CA":        [4, 5, 4, 2],
        "RAS":           [3, 14, 22, 35], "Fusion":        [2, 11, 14, 18],
        "Other":         [0, 10, 13, 15],
    }
    cmap = {"BRAF V600E":"#9c4742", "TERT promoter":"#b25b46", "TP53":"#cd8b3a",
            "PIK3CA":"#d2b454", "RAS":"#37618e", "Fusion":"#7e6e94", "Other":"#a8a59a"}
    bottom = np.zeros(len(classes))
    x = np.arange(len(classes))
    for d, vals in drivers.items():
        ax.bar(x, vals, bottom=bottom, color=cmap[d], edgecolor="white", lw=0.5,
                label=d, width=0.55)
        bottom += vals
    for i, n in enumerate(n_per):
        ax.text(x[i], -8, f"n = {n}", ha="center", fontsize=8.5, color=INK)
    ax.annotate("", xy=(2.42, 108), xytext=(0.58, 108),
                arrowprops=dict(arrowstyle="-", color=GRAY_P, lw=2.2))
    ax.text(1.5, 116, "gray zone 14%", ha="center", fontsize=9, color=GRAY_P, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(classes, fontsize=10, fontweight="bold")
    ax.set_ylabel("% of class (driver)")
    ax.set_ylim(-12, 130)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))
    panel_title(ax, "Mu 2024 HRA004166  ·  Chinese  ·  n = 214  ·  4-class uptake")


def main():
    fig = plt.figure(figsize=(17, 11.5), facecolor=IVORY)
    gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.32,
                          left=0.05, right=0.96, top=0.92, bottom=0.05,
                          height_ratios=[1.05, 0.95, 0.95])

    fig.suptitle("Cross-cohort validation of the 8-gene panel in Korean and Asian thyroid cohorts",
                 fontsize=14.5, fontweight="bold", color=INK, y=0.965)
    fig.text(0.5, 0.93,
             "Korean: Lee 2024 + K2 + GSE286332 + Mun 2025 (1,246 patients) · Chinese: Lu 2023 + Pu 2021 + Mu 2024 (214 patients + 47,587 cells)",
             ha="center", fontsize=10, color=MUTED, style="italic")

    ax_a = fig.add_subplot(gs[0, :]); panel_a_inventory(ax_a); panel_letter(ax_a, "a")
    ax_b = fig.add_subplot(gs[1, 0]); panel_b_lee2024(ax_b); panel_letter(ax_b, "b")
    ax_c = fig.add_subplot(gs[1, 1]); panel_c_k2(ax_c); panel_letter(ax_c, "c")
    ax_d = fig.add_subplot(gs[2, 0]); panel_d_mun2025(ax_d); panel_letter(ax_d, "d")
    ax_e = fig.add_subplot(gs[2, 1]); panel_e_mu2024(ax_e); panel_letter(ax_e, "e")

    out_png = FIG / "figure_korean_asian_validation.png"
    out_pdf = FIG / "figure_korean_asian_validation.pdf"
    fig.savefig(out_png, dpi=200, facecolor=IVORY)
    fig.savefig(out_pdf, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
