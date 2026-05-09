"""Wave 5A — Build comparison figures.

Generates `fig_wave5a_uplift.{png,pdf}`:
  - Bar chart of AUROC on ITSNdb_no_overlap with bootstrap CI95:
      Wave 1 baseline (ESM2-Bayesian, 0.411)
      Wave 5A ablation (no_distill, biophys-only)
      Wave 5A distill (full distillation)
      MHCflurry teacher (0.668)
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W5A = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave5a")


def load_row(tsv_path, testset):
    if not tsv_path.exists():
        return None
    df = pd.read_csv(tsv_path, sep="\t")
    sub = df[df["testset"] == testset]
    if len(sub) == 0:
        return None
    return sub.iloc[0].to_dict()


def main():
    rows = []
    # Wave 1 baseline (hard-coded from BAYESIAN_NEO_REPORT.md)
    rows.append({
        "label": "Wave 1\nESM2-Bayes\n(no distill)",
        "AUROC": 0.411, "lo": None, "hi": None,
        "color": "#888888",
    })
    # Wave 5A ablation
    abl = load_row(W5A / "wave5a_results_ablation.tsv", "ITSNdb_no_overlap")
    if abl:
        rows.append({
            "label": "Wave 5A\nablation\n(no distill)",
            "AUROC": abl["AUROC"], "lo": abl["AUROC_lo95"], "hi": abl["AUROC_hi95"],
            "color": "#cc8855",
        })
    # Wave 5A distill
    dist = load_row(W5A / "wave5a_results_distill.tsv", "ITSNdb_no_overlap")
    if dist:
        rows.append({
            "label": "Wave 5A\ndistill\n(MHCflurry)",
            "AUROC": dist["AUROC"], "lo": dist["AUROC_lo95"], "hi": dist["AUROC_hi95"],
            "color": "#3366aa",
        })
    # MHCflurry teacher (hard-coded — its no_overlap is reported as 0.668 in user's note)
    rows.append({
        "label": "MHCflurry\nteacher",
        "AUROC": 0.668, "lo": None, "hi": None,
        "color": "#229944",
    })

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    xs = np.arange(len(rows))
    bars = ax.bar(xs, [r["AUROC"] for r in rows],
                  color=[r["color"] for r in rows], edgecolor="black", linewidth=0.6)
    for i, r in enumerate(rows):
        if r["lo"] is not None and r["hi"] is not None and not pd.isna(r["lo"]):
            ax.errorbar(i, r["AUROC"], yerr=[[r["AUROC"] - r["lo"]], [r["hi"] - r["AUROC"]]],
                        fmt="none", ecolor="black", capsize=5)
        ax.text(i, r["AUROC"] + 0.018, f"{r['AUROC']:.3f}",
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.axhline(0.5, color="grey", lw=0.5, ls="--", alpha=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([r["label"] for r in rows], fontsize=9)
    ax.set_ylabel("AUROC (ITSNdb_no_overlap)")
    ax.set_ylim(0, 0.85)
    ax.set_title("Wave 5A — MHCflurry distillation uplift on truly-external test set\n"
                 "(ITSNdb_no_overlap, n=106, peptides absent from master training data)",
                 fontsize=10)
    plt.tight_layout()
    fig.savefig(W5A / "fig_wave5a_uplift.png", dpi=180)
    fig.savefig(W5A / "fig_wave5a_uplift.pdf")
    print(f"saved: fig_wave5a_uplift.png/pdf")


if __name__ == "__main__":
    main()
