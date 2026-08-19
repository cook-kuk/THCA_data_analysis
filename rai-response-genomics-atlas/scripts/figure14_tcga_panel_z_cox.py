#!/usr/bin/env python3
"""Figure 14 — TCGA-THCA continuous panel z Cox model + KM by tertile.

Treats panel z as a continuous predictor of PFI and shows Kaplan-Meier curves
stratified by panel z tertile, with Cox HR and log-rank p annotations.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from _nature_style import setup_rc, panel_letter, panel_title, IVORY, INK, MUTED, BLUE, RED, GRAY_P

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

setup_rc()

ROOT = Path(__file__).resolve().parent.parent
TCGA_MERGED = Path("/home/seungho/personal/THCA_data_analysis/project/results/r17_tcga_panel_d4p2_reconciliation/r17_per_sample_merged.tsv")
TCGA_CLIN_CANDIDATES = list(Path("/home/seungho/personal/THCA_data_analysis/project/data").rglob("tcga*clinical*"))
OUT_PNG = ROOT / "results" / "figures" / "figure14_panel_z_continuous_cox.png"
OUT_PDF = ROOT / "results" / "figures" / "figure14_panel_z_continuous_cox.pdf"


CLIN = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv")


def main():
    df = pd.read_csv(TCGA_MERGED, sep="\t").dropna(subset=["RAI_8"])
    clin = pd.read_csv(CLIN, sep="\t")
    df["pat"] = df["sampleId"].str.slice(0, 12)
    merged = df.merge(clin, left_on="pat", right_on="tcga_short", how="inner")
    pfi_t, pfi_e = "PFI.time", "PFI"
    print(f"Using PFI; merged n = {merged[pfi_t].notna().sum()}")

    sub = merged[[pfi_t, pfi_e, "RAI_8"]].dropna().rename(columns={pfi_t:"T", pfi_e:"E", "RAI_8":"panel_z"})
    sub["T"] = sub["T"].astype(float); sub["E"] = sub["E"].astype(int)
    sub = sub[sub["T"] > 0]

    # Continuous Cox PH
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(sub.rename(columns={"T":"duration","E":"event"})[["duration","event","panel_z"]],
            duration_col="duration", event_col="event")
    hr = float(cph.hazard_ratios_["panel_z"])
    ci = cph.confidence_intervals_.loc["panel_z"].values
    hr_lo, hr_hi = float(np.exp(ci[0])), float(np.exp(ci[1]))
    p_cox = float(cph.summary.loc["panel_z","p"])
    print(f"Continuous panel_z HR = {hr:.3f} [{hr_lo:.3f}, {hr_hi:.3f}], p = {p_cox:.3g}")

    # Tertile KM
    sub["tertile"] = pd.qcut(sub["panel_z"], 3, labels=["low (T1)","mid (T2)","high (T3)"])
    lr = multivariate_logrank_test(sub["T"], sub["tertile"], sub["E"])
    p_lr = float(lr.p_value)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), facecolor=IVORY,
                              gridspec_kw=dict(wspace=0.28, left=0.08, right=0.97, top=0.84, bottom=0.16))
    fig.suptitle("Continuous panel z is a graded predictor of PFI in TCGA-THCA",
                 fontsize=13.5, fontweight="bold", color=INK, y=0.97)

    # KM curves
    ax = axes[0]
    colors_t = {"low (T1)": RED, "mid (T2)": "#7e689a", "high (T3)": BLUE}
    for label in ["low (T1)", "mid (T2)", "high (T3)"]:
        kmf = KaplanMeierFitter()
        msk = sub["tertile"] == label
        kmf.fit(sub.loc[msk, "T"], sub.loc[msk, "E"], label=f"{label}  (n = {int(msk.sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors_t[label], lw=1.6, alpha=0.92)
    ax.set_xlabel("time (days)")
    ax.set_ylabel("progression-free probability")
    ax.set_ylim(0.5, 1.02)
    ax.text(0.02, 0.06, f"log-rank p = {p_lr:.3g}", transform=ax.transAxes,
            fontsize=10, color=INK, fontweight="bold")
    panel_title(ax, "Kaplan-Meier by panel z tertile"); panel_letter(ax, "a")

    # Right: HR forest of continuous panel z
    ax = axes[1]
    ax.errorbar(hr, 1, xerr=[[hr - hr_lo], [hr_hi - hr]],
                fmt="D", color=BLUE, ecolor=BLUE, markersize=12, capsize=4, capthick=1.4, lw=1.4,
                mec="white", mew=1.0)
    ax.text(hr_hi*1.05, 1, f"HR = {hr:.2f}  [{hr_lo:.2f}, {hr_hi:.2f}]\np = {p_cox:.3g}",
            fontsize=11, va="center", color=INK)
    ax.axvline(1, color=MUTED, ls="--", lw=0.9)
    ax.set_yticks([1]); ax.set_yticklabels(["continuous panel z\n(per unit increase)"], fontsize=10)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 5)
    ax.set_xticks([0.1, 0.3, 1, 3, 5])
    ax.set_xticklabels(["0.1","0.3","1","3","5"])
    ax.set_xlabel("hazard ratio  (continuous panel z)")
    ax.set_ylim(0, 2.2)
    panel_title(ax, "Continuous Cox PH model"); panel_letter(ax, "b")

    fig.text(0.5, 0.025,
             f"Each unit increase in panel z reduces PFI hazard by {(1-hr)*100:.1f}%; the lowest tertile carries the highest hazard.",
             ha="center", fontsize=9.5, color=INK, style="italic")
    fig.savefig(OUT_PNG, dpi=220, facecolor=IVORY)
    fig.savefig(OUT_PDF, facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {OUT_PNG}\nwrote {OUT_PDF}")


if __name__ == "__main__":
    main()
