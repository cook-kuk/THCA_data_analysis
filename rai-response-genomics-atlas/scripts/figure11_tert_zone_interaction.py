#!/usr/bin/env python3
"""Figure 11 — TERT × R17 zone interaction (PFI event rate stratified)."""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED, GRAY_P, BEIGE

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PFI = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_pfi_by_zone_tert.tsv")
FISHER = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_zone_fisher.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure11_tert_zone_interaction.png"
OUT_PDF = ROOT / "results" / "figures" / "figure11_tert_zone_interaction.pdf"

ZONE_ORDER = ["WT-like", "RAS-like", "BRAF-like", "dark-matter"]
ZONE_COLOR = {"WT-like": "#b5b5b5", "RAS-like": RED, "BRAF-like": BLUE, "dark-matter": "#5a4470"}


def main():
    df = pd.read_csv(PFI, sep="\t")
    fi = pd.read_csv(FISHER, sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.30, left=0.075, right=0.97, top=0.83, bottom=0.16))
    fig.suptitle("TERT promoter × R17 zone interaction on TCGA-THCA",
                 fontsize=13.5, fontweight="bold", color=INK, y=0.97)

    # Left: PFI event rate by zone × TERT
    ax = axes[0]
    zones = ZONE_ORDER
    x = np.arange(len(zones))
    w = 0.36
    tneg_rates = []; tpos_rates = []; tneg_n = []; tpos_n = []
    for z in zones:
        rn = df[(df["zone"] == z) & (df["tert_bool"] == False)]
        rp = df[(df["zone"] == z) & (df["tert_bool"] == True)]
        tneg_rates.append(float(rn["event_rate_pct"].iloc[0]) if len(rn) else 0)
        tpos_rates.append(float(rp["event_rate_pct"].iloc[0]) if len(rp) else 0)
        tneg_n.append(int(rn["n"].iloc[0]) if len(rn) else 0)
        tpos_n.append(int(rp["n"].iloc[0]) if len(rp) else 0)

    b1 = ax.bar(x - w/2, tneg_rates, w, color=[ZONE_COLOR[z] for z in zones], alpha=0.6,
                edgecolor=INK, lw=0.6, label="TERT promoter wild-type")
    b2 = ax.bar(x + w/2, tpos_rates, w, color=[ZONE_COLOR[z] for z in zones], alpha=1.0,
                edgecolor=INK, lw=0.8, label="TERT promoter mutant", hatch="///")

    for xi, n_neg, n_pos, r_neg, r_pos in zip(x, tneg_n, tpos_n, tneg_rates, tpos_rates):
        ax.text(xi - w/2, r_neg + 1.0, f"n = {n_neg}", ha="center", fontsize=8, color=MUTED)
        ax.text(xi + w/2, r_pos + 1.0, f"n = {n_pos}", ha="center", fontsize=8, color=MUTED)
        ax.text(xi - w/2, r_neg/2, f"{r_neg:.1f}%", ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
        if r_pos > 5:
            ax.text(xi + w/2, r_pos/2, f"{r_pos:.1f}%", ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
        else:
            ax.text(xi + w/2, r_pos + 1.6, f"{r_pos:.1f}%", ha="center", fontsize=8.5, color=INK, fontweight="bold")

    ax.set_xticks(x); ax.set_xticklabels(zones, fontsize=9.5)
    ax.set_ylabel("PFI event rate (%)")
    ax.set_ylim(0, max(max(tpos_rates), max(tneg_rates))*1.30 + 5)
    ax.legend(fontsize=8.5, loc="upper left", frameon=True)
    panel_title(ax, "PFI event rate by zone × TERT promoter"); panel_letter(ax, "a")

    # Right: Fisher OR for TERT enrichment per zone
    ax = axes[1]
    fi_sorted = fi.set_index("zone").loc[zones].reset_index()
    ors = fi_sorted["OR"].values
    ps = fi_sorted["p"].values
    # Replace OR=0 with small floor for log axis
    ors_plot = np.where(ors == 0, 0.05, ors)
    colors = [ZONE_COLOR[z] for z in fi_sorted["zone"]]
    bars = ax.barh(x, ors_plot, color=colors, alpha=0.85, edgecolor=INK, lw=0.6)
    ax.axvline(1, color=MUTED, ls="--", lw=0.9)
    ax.set_yticks(x); ax.set_yticklabels(fi_sorted["zone"], fontsize=9.5)
    ax.set_xscale("log")
    ax.set_xlim(0.04, 10)
    ax.set_xticks([0.05, 0.1, 0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["≤0.05", "0.1", "0.5", "1", "2", "5", "10"])
    ax.set_xlabel("Fisher odds ratio  (TERT+ enrichment vs rest)")
    for xi, or_v, p_v, or_plot in zip(x, ors, ps, ors_plot):
        label = f"OR = {or_v:.2f}, p = {p_v:.3f}" if or_v > 0 else f"OR = 0, p = {p_v:.3f}"
        ax.text(or_plot*1.06, xi, label, va="center", fontsize=8.8, color=INK)
    panel_title(ax, "TERT promoter Fisher enrichment by zone"); panel_letter(ax, "b")

    fig.text(0.5, 0.045,
             "TERT promoter mutation is uniquely enriched in the dark-matter zone (OR = 2.34, p = 0.016) "
             "and amplifies PFI event rate (13.6% → 27.8%).",
             ha="center", fontsize=9.5, color=INK, style="italic")

    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
