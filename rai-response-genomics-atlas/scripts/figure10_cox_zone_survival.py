#!/usr/bin/env python3
"""Figure 10 — TCGA-THCA per-zone Cox survival forest (PFI / OS / DSS).

Uses pre-computed r17_zone_cox_results.tsv from the parent project's
zone-survival module. Plots univariate + adjusted HR with 95% CI for
BRAF-like, dark-matter, RAS-like zones vs wild-type-like (reference).
"""
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
COX = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/zone_survival/r17_zone_cox_results.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure10_zone_cox_forest.png"
OUT_PDF = ROOT / "results" / "figures" / "figure10_zone_cox_forest.pdf"


ZONE_ORDER = ["zone_BRAF", "zone_dark", "zone_RAS"]
ZONE_LABEL = {"zone_BRAF": "BRAF-like", "zone_dark": "dark-matter", "zone_RAS": "RAS-like"}
ZONE_COLOR = {"zone_BRAF": BLUE, "zone_dark": "#5a4470", "zone_RAS": RED}


def main():
    df = pd.read_csv(COX, sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.45, left=0.10, right=0.985, top=0.85, bottom=0.18))
    fig.suptitle("Per-zone Cox proportional hazards on TCGA-THCA (n = 552–563)",
                 fontsize=13.5, fontweight="bold", color=INK, y=0.97)

    outcomes = ["PFI", "OS", "DSS"]
    panel_letters = ["a", "b", "c"]

    for ax, oc, ltr in zip(axes, outcomes, panel_letters):
        sub = df[df["outcome"] == oc].copy()
        rows = []
        for model in ["univariate", "adjusted_age_sex_stage"]:
            for zone in ZONE_ORDER:
                r = sub[(sub["model"] == model) & (sub["term"] == zone)]
                if len(r) == 0: continue
                rows.append({
                    "model": "univariate" if model == "univariate" else "adjusted",
                    "zone": zone,
                    "HR": float(r["HR"].iloc[0]),
                    "lo": float(r["HR_lower_95"].iloc[0]),
                    "hi": float(r["HR_upper_95"].iloc[0]),
                    "p":  float(r["p"].iloc[0]),
                    "n":  int(r["n"].iloc[0]),
                    "ev": int(r["events"].iloc[0]),
                })
        plot = pd.DataFrame(rows)

        # Y positions: alternating uni/adj for each zone
        y_pos = []
        labels = []
        colors = []
        for i, zone in enumerate(ZONE_ORDER):
            for j, model in enumerate(["univariate", "adjusted"]):
                row = plot[(plot["zone"] == zone) & (plot["model"] == model)]
                if len(row) == 0: continue
                y = i*2.5 + (0 if model == "univariate" else 0.95)
                y_pos.append(y)
                r = row.iloc[0]
                lab = f"{ZONE_LABEL[zone]} · {model}\n   HR {r['HR']:.2f}  p {r['p']:.3g}"
                labels.append(lab)
                colors.append(ZONE_COLOR[zone] if model == "univariate" else "#888888")
                ax.errorbar(r["HR"], y,
                            xerr=[[max(r["HR"] - r["lo"], 0)], [max(r["hi"] - r["HR"], 0)]],
                            fmt="o", color=colors[-1], ecolor=colors[-1],
                            markersize=7, capsize=3.5, capthick=1.2, lw=1.2, alpha=0.92,
                            mec="white", mew=0.8)

        ax.axvline(1.0, color=MUTED, ls="--", lw=0.9, alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8.5)
        ax.set_xlabel("hazard ratio  (vs wild-type-like)")
        ax.set_xscale("log")
        ax.set_xlim(0.1, 6)
        ax.set_xticks([0.1, 0.3, 1, 3, 6])
        ax.set_xticklabels(["0.1", "0.3", "1", "3", "6"])
        n_ev_row = plot.iloc[0] if len(plot) > 0 else None
        if n_ev_row is not None:
            panel_title(ax, f"{oc} · n = {n_ev_row['n']:,} (events = {n_ev_row['ev']})")
        panel_letter(ax, ltr)

    fig.text(0.5, 0.04,
             "Dark-matter and BRAF-like zones show consistent PFI hazard elevation; OS / DSS underpowered (events ≤ 14).",
             ha="center", fontsize=9.5, color=INK, style="italic")

    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
