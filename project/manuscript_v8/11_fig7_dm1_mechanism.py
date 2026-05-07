"""
Figure 7 — DM1 mechanism (paradigm) Cell Press style 5-panel figure.
Paper 1 manuscript v8.

Panels:
  A — Sub-A vs Sub-B silhouette + score profile (silhouette 0.584)
  B — Fusion enrichment 76.8% vs 30.9% (Fisher OR 7.41)
  C — Fusion partner stack (RET/NTRK/ALK/BRAF)
  D — Sub-A vs Sub-B phenotype panel (age, stage, immune)
  E — DM1 captures 81.8% TCGA RET+ (reflex algorithm)

Inputs:
  - results/audit_2026_04_30/round3/cbio_sv_thca.tsv          (TCGA fusion landscape)
  - results/audit_2026_04_30/round4/r4_2_fusion_phenotype.tsv  (Sub-A vs Sub-B)
  - results/audit_2026_04_30/round4/r4_3_reflex_capture.tsv    (RET capture)
  - results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv      (DM1 sub-A/B)
  - results/tables/tcga_thca_clinical_extended.tsv             (clinical)

Output: project/manuscript_v8/figures/Fig7_dm1_mechanism.{png,pdf}

Usage: python 11_fig7_dm1_mechanism.py
"""

from __future__ import annotations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parents[1]  # project/ (was parents[2] — bugfix 2026-05-07)
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_DIR = PROJECT_ROOT / "manuscript_v8" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
})

# Cell Press color palette
COLOR_DM1 = "#C44E52"    # red
COLOR_DM2 = "#4C72B0"    # blue
COLOR_FUSION = "#55A868"  # green
COLOR_NOFUSION = "#CCCCCC"


# -----------------------------------------------------------------------------
# Data loaders (stubs — fill from existing TSV files)
# -----------------------------------------------------------------------------
def load_subcluster_labels() -> pd.DataFrame:
    """DM1 sub-A vs sub-B per-sample labels (n=140)."""
    return pd.read_csv(
        RESULTS_DIR / "d6p7_dm1_subcluster" / "dm1_subcluster_labels.tsv",
        sep="\t",
    )


def load_subcluster_scores() -> pd.DataFrame:
    """Sub-A/B silhouette + 8-gene score profile."""
    return pd.read_csv(
        RESULTS_DIR / "d6p7_dm1_subcluster" / "subcluster_scores.tsv",
        sep="\t",
    )


def load_sv_data() -> pd.DataFrame:
    """TCGA fusion landscape (RET/NTRK/ALK/BRAF) — n=542 SV-tested."""
    return pd.read_csv(
        RESULTS_DIR / "audit_2026_04_30" / "round3" / "cbio_sv_thca.tsv",
        sep="\t",
    )


# -----------------------------------------------------------------------------
# Panel A — Sub-A vs Sub-B silhouette
# -----------------------------------------------------------------------------
def panel_A(ax: plt.Axes) -> None:
    """A — Silhouette plot of DM1 sub-A (n=72) vs sub-B (n=19), score 0.584."""
    # TODO populate from subcluster_scores.tsv
    sub_a_scores = np.random.normal(0.7, 0.15, 72)  # placeholder
    sub_b_scores = np.random.normal(0.3, 0.12, 19)  # placeholder

    ax.barh(np.arange(len(sub_a_scores)), sorted(sub_a_scores, reverse=True),
            color=COLOR_DM1, alpha=0.8, label="Sub-A (n=72)")
    ax.barh(np.arange(len(sub_a_scores), len(sub_a_scores) + len(sub_b_scores)),
            sorted(sub_b_scores, reverse=True),
            color="#E07B7B", alpha=0.8, label="Sub-B (n=19)")
    ax.axvline(0.584, color="black", linestyle="--", linewidth=0.5, alpha=0.5)
    ax.set_xlabel("Silhouette coefficient")
    ax.set_ylabel("Sample (sorted)")
    ax.set_title("(A) DM1 sub-A vs sub-B silhouette\n(score = 0.584)")
    ax.legend(loc="lower right")


# -----------------------------------------------------------------------------
# Panel B — Fusion enrichment OR 7.41
# -----------------------------------------------------------------------------
def panel_B(ax: plt.Axes) -> None:
    """B — DM1 76.8% vs DM2 30.9% fusion+ stacked bar."""
    clusters = ["DM1\n(n=82)", "DM2\n(n=460)"]
    fusion_pos_pct = [76.8, 30.9]
    fusion_neg_pct = [100 - p for p in fusion_pos_pct]

    ax.bar(clusters, fusion_pos_pct, color=COLOR_FUSION, label="Fusion+")
    ax.bar(clusters, fusion_neg_pct, bottom=fusion_pos_pct,
           color=COLOR_NOFUSION, label="Fusion−")
    for i, p in enumerate(fusion_pos_pct):
        ax.text(i, p / 2, f"{p:.1f}%", ha="center", va="center", fontweight="bold")
    ax.set_ylabel("% of DM cluster")
    ax.set_ylim(0, 100)
    ax.set_title("(B) Fusion enrichment\nFisher OR 7.41 (CI 4.4–12.6), p=1.9e-13")
    ax.legend(loc="upper right")


# -----------------------------------------------------------------------------
# Panel C — Fusion partner stack
# -----------------------------------------------------------------------------
def panel_C(ax: plt.Axes) -> None:
    """C — Fusion partner spectrum within DM1 fusion+ (n=63)."""
    partners = {
        "RET (CCDC6)": 17,
        "RET (NCOA4)": 3,
        "RET (other)": 13,
        "NTRK1/3": 10,
        "ALK": 4,
        "BRAF fusion": 5,
        "Other": 11,
    }
    labels = list(partners.keys())
    counts = list(partners.values())
    colors = ["#C44E52", "#E07B7B", "#F4A0A0", "#55A868", "#8172B2", "#CCB974", "#999999"]

    bars = ax.barh(labels, counts, color=colors)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                f"n={count}", va="center", fontsize=7)
    ax.set_xlabel("# DM1 fusion+ cases")
    ax.set_title("(C) Fusion partner spectrum\n(DM1 fusion+, n=63)")
    ax.invert_yaxis()


# -----------------------------------------------------------------------------
# Panel D — Sub-A vs Sub-B phenotype
# -----------------------------------------------------------------------------
def panel_D(ax: plt.Axes) -> None:
    """D — Phenotype contrast: age, stage III/IV, fusion%, immune signature."""
    metrics = ["Age (yr)", "Stage III/IV (%)", "Fusion+ (%)", "CD8/IFN-γ z"]
    sub_a_vals = [37.3, 15.3, 84.7, -0.5]
    sub_b_vals = [51.3, 44.4, 57.9, 0.8]
    x = np.arange(len(metrics))
    width = 0.35

    ax.bar(x - width / 2, sub_a_vals, width, color=COLOR_DM1, label="Sub-A")
    ax.bar(x + width / 2, sub_b_vals, width, color="#E07B7B", label="Sub-B")
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=20, ha="right")
    ax.set_title("(D) Sub-A vs sub-B phenotype")
    ax.legend()
    ax.axhline(0, color="black", linewidth=0.5)


# -----------------------------------------------------------------------------
# Panel E — DM1 captures 81.8% TCGA RET+
# -----------------------------------------------------------------------------
def panel_E(ax: plt.Axes) -> None:
    """E — DM1 reflex algorithm: 27/33 = 81.8% RET+ capture."""
    # Donut chart of TCGA RET-fusion-positive tumors (n=33) by DM cluster
    captured = 27
    missed = 6
    sizes = [captured, missed]
    labels = [f"DM1+\n(captured)\n{captured}/33", f"DM2\n(missed)\n{missed}/33"]
    colors = [COLOR_DM1, COLOR_DM2]

    wedges, texts = ax.pie(sizes, labels=labels, colors=colors,
                            wedgeprops={"width": 0.4, "edgecolor": "white"},
                            startangle=90)
    ax.text(0, 0, f"{captured/(captured+missed)*100:.1f}%\nRET+\ncapture",
            ha="center", va="center", fontsize=10, fontweight="bold")
    ax.set_title("(E) DM1 captures 81.8% of TCGA RET+\n(reflex testing algorithm)")


# -----------------------------------------------------------------------------
# Main figure assembly
# -----------------------------------------------------------------------------
def build_figure() -> None:
    fig = plt.figure(figsize=(8.5, 7.5))
    gs = GridSpec(3, 2, figure=fig, hspace=0.5, wspace=0.35,
                  height_ratios=[1, 1, 1])
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, :])

    panel_A(axA)
    panel_B(axB)
    panel_C(axC)
    panel_D(axD)
    panel_E(axE)

    fig.suptitle(
        "Figure 7. DM1 is a fusion-driven dark matter subtype with mechanistic heterogeneity.",
        fontsize=10, fontweight="bold", y=0.995
    )
    fig.savefig(OUT_DIR / "Fig7_dm1_mechanism.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Fig7_dm1_mechanism.pdf", bbox_inches="tight")
    print(f"Figure 7 written: {OUT_DIR / 'Fig7_dm1_mechanism.png'}")


if __name__ == "__main__":
    build_figure()
