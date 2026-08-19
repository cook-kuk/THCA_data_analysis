#!/usr/bin/env python3
"""Figure 4 v2 — Nature-style molecular gray-zone framework.

4 sub-panels: TCGA driver × zone heatmap, Mu 2024 4-class composition,
GSE151179 panel z by driver, GSE112202 redifferentiation.
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
from matplotlib.colors import LinearSegmentedColormap

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PAPER1_R17 = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/r17_per_sample_merged.tsv")
GSE151179 = ROOT / "data" / "processed" / "GSE151179_label_joined.tsv"
FIG = ROOT / "results" / "figures"


def panel_a_tcga_heatmap(ax):
    df = pd.read_csv(PAPER1_R17, sep="\t")
    df["panel_silenced"] = df["DM_call"] == "DM1"
    df["ht_high"] = df["sig_score"] > 0
    def zone_of(row):
        if row["panel_silenced"] and row["ht_high"]: return "dark-matter"
        if row["panel_silenced"]: return "BRAF-like"
        if row["ht_high"]: return "RAS-like"
        return "WT-like"
    df["zone"] = df.apply(zone_of, axis=1)
    drivers = ["BRAF V600E", "RAS", "BRAF·RAS-neg"]
    zones   = ["BRAF-like", "RAS-like", "dark-matter", "WT-like"]
    ct = pd.crosstab(df["driver"], df["zone"]).reindex(index=drivers, columns=zones, fill_value=0)
    pct = ct.div(ct.sum(axis=1), axis=0) * 100

    cmap = LinearSegmentedColormap.from_list("nat_purp", ["#fbf8f1", "#dbcae0", "#9e84b2", "#5e4a7a"], N=200)
    im = ax.imshow(pct.values, cmap=cmap, vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(zones))); ax.set_xticklabels(zones, rotation=12, fontsize=9.5)
    ax.set_yticks(range(len(drivers)))
    ax.set_yticklabels([f"{d}\n(n = {int(ct.loc[d].sum())})" for d in drivers], fontsize=9.5)
    for i, drv in enumerate(drivers):
        for j, zone in enumerate(zones):
            v = pct.iloc[i, j]; n = ct.iloc[i, j]
            ax.text(j, i, f"{v:.0f}%\nn = {n}", ha="center", va="center", fontsize=8.5,
                    color="white" if v > 55 else INK, fontweight="bold" if v > 30 else "normal")
    cb = plt.colorbar(im, ax=ax, fraction=0.044, pad=0.04)
    cb.set_label("% of driver in zone", fontsize=8.5)
    cb.outline.set_linewidth(0.5)
    panel_title(ax, "TCGA-THCA  ·  driver × panel zone  ·  R17 L1  ·  n = 500")
    ax.set_xlabel("")


def panel_b_mu2024(ax):
    classes = ["I-RAIR", "G-RAIR", "P-RAIR", "C-RAIA"]
    n_per = [80, 19, 10, 48]
    drivers = {
        "BRAF V600E":    [61, 36, 30, 22],
        "TERT promoter": [22, 18, 12, 5],
        "TP53":          [8, 6, 5, 3],
        "PIK3CA":        [4, 5, 4, 2],
        "RAS":           [3, 14, 22, 35],
        "Fusion":        [2, 11, 14, 18],
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
        ax.text(x[i], -7, f"n = {n}", ha="center", fontsize=8.5, color=INK)
    ax.annotate("", xy=(2.42, 108), xytext=(0.58, 108),
                arrowprops=dict(arrowstyle="-", color=GRAY_P, lw=2.2))
    ax.plot([0.58, 0.58], [105, 111], color=GRAY_P, lw=1.5)
    ax.plot([2.42, 2.42], [105, 111], color=GRAY_P, lw=1.5)
    ax.text(1.5, 116, "gray zone  (G-RAIR + P-RAIR = 29 / 214)",
            ha="center", fontsize=9.2, color=GRAY_P, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(classes, fontsize=10, fontweight="bold")
    ax.set_ylabel("% of class  (driver composition)")
    ax.set_ylim(-12, 130)
    ax.legend(fontsize=7.5, loc="upper left", ncol=1, framealpha=0.92, bbox_to_anchor=(1.02, 1))
    panel_title(ax, "Mu 2024 HRA004166  ·  Chinese  ·  n = 214  ·  4-class RAI uptake")


def panel_c_gse151179_driver(ax):
    df = pd.read_csv(GSE151179, sep="\t", index_col=0)
    df = df[~df["sample_type"].astype(str).str.contains("non-neoplastic", case=False, na=False)]
    df["lc"] = df["lesion_class"].astype(str).str.lower()
    orders = ["brafv600e", "fusion", "ptert", "wt"]
    labels = ["BRAF V600E", "Fusion", "pTERT", "WT"]
    data = [df.loc[df["lc"] == d, "panel_z"].dropna().values for d in orders]
    colors = ["#9c4742", "#7e6e94", "#cd8b3a", "#37618e"]
    bp = ax.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True, widths=0.5,
                    medianprops={"color": "#222", "linewidth": 1.2},
                    whiskerprops={"color": INK, "linewidth": 0.7},
                    capprops={"color": INK, "linewidth": 0.7})
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor(c); patch.set_linewidth(0.8)
    rng = np.random.default_rng(7)
    for i, vals in enumerate(data):
        x = rng.uniform(-0.10, 0.10, size=len(vals)) + (i + 1)
        ax.scatter(x, vals, s=22, alpha=0.85, c=colors[i], edgecolors="white", linewidths=0.5)
        ax.text(i + 1, 0.04, f"n = {len(vals)}", ha="center", va="bottom", fontsize=8.5,
                color=MUTED, transform=ax.get_xaxis_transform())
    ax.axhline(0, color="#999", lw=0.5, ls="--", alpha=0.7)
    ax.set_ylabel("8-gene panel z")
    panel_title(ax, "GSE151179 tumour-only  ·  panel z by driver class  ·  Tier 2")


def panel_d_gse112202_redifferentiation(ax):
    """Recompute from local GSE112202 fpkm to ensure styling consistency."""
    import gzip
    src = ROOT / "data" / "raw" / "GSE112202" / "GSE112202_overall_comparison.fpkm_tracking.gz"
    df = pd.read_csv(src, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    untr = next(c for c in df.columns if "untreated" in c.lower() and c.endswith("_FPKM"))
    dig  = next(c for c in df.columns if "digoxin"   in c.lower() and c.endswith("_FPKM"))
    sym = "gene_short_name"
    df = df[[sym, untr, dig]].copy()
    df[untr] = pd.to_numeric(df[untr], errors="coerce")
    df[dig]  = pd.to_numeric(df[dig], errors="coerce")
    df = df.dropna().drop_duplicates(subset=sym, keep="first")
    df["log2fc"] = np.log2((df[dig] + 0.1) / (df[untr] + 0.1))

    PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
    TDS_extra = ["DIO2", "DUOX1", "DUOX2", "GLIS3", "SLC26A4", "SLC5A8", "THRA", "THRB"]
    sub = df[df[sym].isin(PANEL_8 + TDS_extra)].copy()
    sub["set"] = sub[sym].apply(lambda g: "8-gene panel" if g in PANEL_8 else "TDS-16 extras")
    sub = sub.sort_values("log2fc")

    colors = {"8-gene panel": BLUE, "TDS-16 extras": GRAY_P}
    edges  = {"8-gene panel": "#22456b", "TDS-16 extras": "#5b4d70"}
    y = np.arange(len(sub))
    for s_name in ["8-gene panel", "TDS-16 extras"]:
        idx = sub[sub["set"] == s_name].index
        plot_pos = [list(sub.index).index(i) for i in idx]
        vals = sub.loc[idx, "log2fc"]
        ax.barh(plot_pos, vals, color=colors[s_name], edgecolor=edges[s_name], lw=0.5,
                alpha=0.85, label=f"{s_name} ({len(idx)})", height=0.65)
    ax.axvline(0, color="#444", lw=0.5, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(sub[sym], fontsize=8.5)
    ax.set_xlabel("log₂ FC  (digoxin / untreated)")
    n_panel = (sub["set"] == "8-gene panel").sum()
    n_up = ((sub["set"] == "8-gene panel") & (sub["log2fc"] > 0)).sum()
    median_panel = sub.loc[sub["set"] == "8-gene panel", "log2fc"].median()
    ax.text(0.98, 0.04,
            f"8-gene panel up: {n_up} / {n_panel}  ·  median log₂FC = {median_panel:+.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=BLUE,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=BLUE, linewidth=0.5))
    ax.legend(fontsize=8.5, loc="lower right", bbox_to_anchor=(1, 0.16))
    panel_title(ax, "GSE112202 redifferentiation  ·  digoxin n = 11 vs untreated n = 11  ·  Tier 5")


def main():
    fig = plt.figure(figsize=(15.5, 11.5), facecolor=IVORY)
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.30,
                          left=0.07, right=0.96, top=0.92, bottom=0.06)
    fig.suptitle("Molecular gray-zone framework  ·  driver mutations alone do not partition RAI response classes",
                 fontsize=14, fontweight="bold", color=INK, y=0.97)
    fig.text(0.5, 0.937,
             "Tier 4 discovery (TCGA-THCA R17)  ·  Tier 2 gray-zone anchor (Mu 2024)  ·  Tier 2 driver subset (GSE151179)  ·  Tier 5 redifferentiation (GSE112202)",
             ha="center", fontsize=10, color=MUTED, style="italic")

    ax_a = fig.add_subplot(gs[0, 0]); panel_a_tcga_heatmap(ax_a); panel_letter(ax_a, "a")
    ax_b = fig.add_subplot(gs[0, 1]); panel_b_mu2024(ax_b); panel_letter(ax_b, "b")
    ax_c = fig.add_subplot(gs[1, 0]); panel_c_gse151179_driver(ax_c); panel_letter(ax_c, "c")
    ax_d = fig.add_subplot(gs[1, 1]); panel_d_gse112202_redifferentiation(ax_d); panel_letter(ax_d, "d")

    out_png = FIG / "figure4_driver_grayzone.png"
    out_pdf = FIG / "figure4_driver_grayzone.pdf"
    fig.savefig(out_png, dpi=210, facecolor=IVORY)
    fig.savefig(out_pdf, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}\nwrote {out_pdf}")


if __name__ == "__main__":
    main()
