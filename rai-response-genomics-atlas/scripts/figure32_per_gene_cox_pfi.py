#!/usr/bin/env python3
"""Figure 32 — Per-gene continuous Cox PFI forest.

Each of 8 panel genes is fit as a continuous PFI predictor in TCGA-THCA
(n=477). Reports per-gene HR (95% CI) and p. Shows the panel is a
collaboration, not a single dominant gene.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
PER = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv")
OUT_PNG = ROOT / "results" / "figures" / "figure32_per_gene_cox_pfi.png"
OUT_PDF = ROOT / "results" / "figures" / "figure32_per_gene_cox_pfi.pdf"

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]


def main():
    df = pd.read_csv(PER, sep="\t").dropna(subset=["PFI.time", "PFI"])
    df = df[df["PFI.time"] > 0]
    print(f"n = {len(df)}, events = {int(df['PFI'].sum())}")

    rows = []
    for g in PANEL:
        if g not in df.columns: continue
        sub = df[[g, "PFI.time", "PFI"]].dropna().rename(columns={g: "x", "PFI.time": "duration", "PFI": "event"})
        if len(sub) < 30: continue
        # log10 + 1 to stabilize, then z-score for comparability across genes
        sub["x"] = np.log10(sub["x"].clip(lower=1e-3) + 1.0)
        sub["x"] = (sub["x"] - sub["x"].mean()) / sub["x"].std()
        cph = CoxPHFitter(penalizer=0.01)
        cph.fit(sub[["duration", "event", "x"]], duration_col="duration", event_col="event")
        hr = float(cph.hazard_ratios_["x"])
        ci = cph.confidence_intervals_.loc["x"].values
        rows.append({"gene": g, "HR": hr, "lo": float(np.exp(ci[0])), "hi": float(np.exp(ci[1])),
                     "p": float(cph.summary.loc["x", "p"])})

    res = pd.DataFrame(rows).sort_values("HR")
    print(res)

    fig, ax = plt.subplots(figsize=(11, 5.5), facecolor=IVORY)
    fig.subplots_adjust(left=0.20, right=0.85, top=0.85, bottom=0.16)
    fig.suptitle("Per-gene continuous Cox PFI on TCGA-THCA (n = " + f"{len(df):,}, events = {int(df['PFI'].sum())})",
                 fontsize=13, fontweight="bold", color=INK, y=0.96)

    y = np.arange(len(res))[::-1]
    for yi, r in zip(y, res.itertuples()):
        col = BLUE if r.HR < 1 else RED
        ax.errorbar(r.HR, yi, xerr=[[r.HR - r.lo], [r.hi - r.HR]],
                    fmt="o", color=col, ecolor=col, markersize=8,
                    capsize=4, capthick=1.2, lw=1.2, mec="white", mew=0.8, alpha=0.92)
        ax.text(1.02, yi, f"HR = {r.HR:.2f}  [{r.lo:.2f}, {r.hi:.2f}]  p = {r.p:.3g}",
                transform=ax.get_yaxis_transform(), fontsize=9.5, va="center", color=INK)

    ax.axvline(1, color=MUTED, lw=0.9, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(res["gene"].values, fontstyle="italic")
    ax.set_xscale("log")
    ax.set_xlim(0.4, 1.6)
    ax.set_xticks([0.5, 0.7, 1, 1.4])
    ax.set_xticklabels(["0.5", "0.7", "1", "1.4"])
    ax.set_xlabel("hazard ratio  (per 1 SD of log10 expression)")
    panel_title(ax, "Each gene independently predicts PFI direction-of-effect")

    n_sig = (res["p"] < 0.05).sum()
    n_protective = (res["HR"] < 1).sum()
    fig.text(0.5, 0.04,
             f"{n_protective}/{len(res)} genes show protective direction (HR < 1); "
             f"{n_sig} of {len(res)} are individually significant at p < 0.05 — "
             "the panel aggregates multiple independent signals.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
