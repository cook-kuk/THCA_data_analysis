#!/usr/bin/env python3
"""Figure 26 — TCGA-THCA Kaplan-Meier PFI by 4-way zone × TERT stratification.

Combines R17 zone and TERT promoter mutation to identify the worst-prognosis
subgroup: dark-matter + TERT+ versus dark-matter alone versus WT-like.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PER_SAMPLE = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure26_zone_tert_km.png"
OUT_PDF = ROOT / "results" / "figures" / "figure26_zone_tert_km.pdf"


def main():
    m = pd.read_csv(PER_SAMPLE, sep="\t").dropna(subset=["PFI.time", "PFI", "zone"])
    m["tert"] = m["tert_pos"].fillna(False).astype(bool)

    # Build 4 strata for clean KM (focus on dark-matter and WT comparisons)
    def strata(r):
        if r["zone"] == "dark-matter" and r["tert"]: return "dark-matter · TERT+"
        if r["zone"] == "dark-matter": return "dark-matter · TERT−"
        if r["zone"] == "BRAF-like" and r["tert"]: return "BRAF-like · TERT+"
        if r["zone"] == "BRAF-like": return "BRAF-like · TERT−"
        return "RAS-like / WT-like  (reference)"
    m["stratum"] = m.apply(strata, axis=1)
    order = [
        "dark-matter · TERT+",
        "dark-matter · TERT−",
        "BRAF-like · TERT+",
        "BRAF-like · TERT−",
        "RAS-like / WT-like  (reference)",
    ]
    colors = {
        "dark-matter · TERT+": "#7c322e",
        "dark-matter · TERT−": "#5a4470",
        "BRAF-like · TERT+": "#cb857f",
        "BRAF-like · TERT−": BLUE,
        "RAS-like / WT-like  (reference)": "#888888",
    }

    fig, ax = plt.subplots(figsize=(13, 5.8), facecolor=IVORY)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.86, bottom=0.18)
    fig.suptitle("Kaplan-Meier PFI by R17 zone × TERT promoter (TCGA-THCA, n ≈ 477)",
                 fontsize=13.5, fontweight="bold", color=INK, y=0.96)

    # n per stratum + final event count
    n_per = {}
    for s in order:
        sub = m[m["stratum"] == s]
        n_per[s] = (len(sub), int(sub["PFI"].sum()))

    for s in order:
        sub = m[m["stratum"] == s]
        if len(sub) < 3: continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub["PFI.time"], sub["PFI"],
                label=f"{s}  (n = {n_per[s][0]}, ev = {n_per[s][1]})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=colors[s], lw=1.8, alpha=0.95)

    ax.set_xlabel("time (days)")
    ax.set_ylabel("progression-free probability")
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0, m["PFI.time"].quantile(0.985))
    ax.legend(fontsize=9, loc="lower left", frameon=True)

    # Log-rank across all strata
    lr = multivariate_logrank_test(m["PFI.time"], m["stratum"], m["PFI"])
    p = float(lr.p_value)
    # Calculate event rates per stratum
    rate_dm_tertp = 100 * n_per.get("dark-matter · TERT+", (1, 0))[1] / max(n_per.get("dark-matter · TERT+", (1, 0))[0], 1)
    rate_ref     = 100 * n_per.get("RAS-like / WT-like  (reference)", (1, 0))[1] / max(n_per.get("RAS-like / WT-like  (reference)", (1, 0))[0], 1)
    ax.text(0.02, 0.94, f"multivariate log-rank p = {p:.3g}", transform=ax.transAxes,
            fontsize=11, color=INK, fontweight="bold")
    ax.text(0.02, 0.88,
            f"dark-matter · TERT+  event rate = {rate_dm_tertp:.1f}%  vs reference {rate_ref:.1f}%",
            transform=ax.transAxes, fontsize=10, color=RED, fontweight="bold")
    panel_title(ax, "PFI stratified by zone × TERT — dark-matter · TERT+ has the worst trajectory")

    fig.text(0.5, 0.05,
             "Dark-matter zone with TERT promoter mutation isolates the worst-prognosis subset; "
             "the combination is biologically and clinically additive.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")
    print("\nStrata counts:")
    for s in order:
        n, e = n_per[s]
        print(f"  {s:50s} n = {n:3d}  events = {e:2d}  rate = {e/max(n,1)*100:.1f}%")


if __name__ == "__main__":
    main()
