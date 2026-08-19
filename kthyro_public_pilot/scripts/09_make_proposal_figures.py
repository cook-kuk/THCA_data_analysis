#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pilot_utils import ensure_standard_dirs, pilot_root, savefig, setup_logging, zscore_series


COLORS = {
    "RAI": "#2aa7ff",
    "HLA/APM": "#8a63d2",
    "CD8/T cell": "#20b486",
    "Myeloid/CAF": "#f28e2b",
    "Drug delivery": "#f4c542",
    "Aggressive": "#e84a5f",
}


def load(name: str) -> pd.DataFrame:
    p = pilot_root() / "results" / "tables" / name
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def darken(fig, axes):
    fig.patch.set_facecolor("#111318")
    for ax in np.ravel(axes):
        ax.set_facecolor("#111318")
        ax.tick_params(colors="white")
        ax.xaxis.label.set_color("white")
        ax.yaxis.label.set_color("white")
        ax.title.set_color("white")
        for spine in ax.spines.values():
            spine.set_color("#777777")


def draw_flow(dark=False):
    fig, ax = plt.subplots(figsize=(13, 6))
    if dark:
        fig.patch.set_facecolor("#111318")
        ax.set_facecolor("#111318")
        txt = "white"
    else:
        txt = "black"
    ax.axis("off")
    boxes = [
        ("TCGA bulk", 0.08, 0.65, COLORS["RAI"]),
        ("Public spatial", 0.08, 0.40, COLORS["HLA/APM"]),
        ("scRNA reference", 0.08, 0.15, COLORS["CD8/T cell"]),
        ("DepMap/PRISM", 0.38, 0.15, COLORS["Drug delivery"]),
        ("STVS score table", 0.38, 0.58, "#4c78a8"),
        ("Spatial niche map", 0.63, 0.58, "#59a14f"),
        ("Candidate therapy axes", 0.63, 0.22, "#e15759"),
        ("Samsung validation design", 0.82, 0.40, "#b07aa1"),
    ]
    for text, x, y, color in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.18, 0.12, color=color, alpha=0.9))
        ax.text(x + 0.09, y + 0.06, text, ha="center", va="center", color="white", fontsize=12, weight="bold")
    arrows = [((0.26, 0.71), (0.38, 0.64)), ((0.26, 0.46), (0.63, 0.64)), ((0.26, 0.21), (0.63, 0.28)), ((0.56, 0.21), (0.63, 0.28)), ((0.56, 0.64), (0.82, 0.46)), ((0.81, 0.64), (0.82, 0.46)), ((0.81, 0.28), (0.82, 0.46))]
    for a, b in arrows:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", color=txt, lw=1.8))
    ax.text(0.5, 0.92, "Figure 1. Public Pilot Data Flow", ha="center", color=txt, fontsize=18, weight="bold")
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description="Create proposal-quality summary figures.")
    parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("09_make_proposal_figures")
    fig_dir = pilot_root() / "results" / "figures"
    prop = fig_dir / "proposal"
    prop.mkdir(parents=True, exist_ok=True)
    tcga = load("tcga_thca_patient_vulnerability_scores.tsv")
    spatial = load("spatial_slide_niche_summary.tsv")
    spots = load("spatial_spot_vulnerability_scores.tsv")
    coherence = load("spatial_coherence_statistics.tsv")
    drug = load("drug_pilot_candidate_rankings.tsv")

    for dark in [False, True]:
        fig = draw_flow(dark=dark)
        savefig(fig, prop / f"figure1_public_pilot_data_flow{'_dark' if dark else ''}.png")
        plt.close(fig)

    if not tcga.empty:
        axes = ["rai_differentiation_score", "immune_visibility_score", "cd8_exclusion_proxy", "drug_delivery_failure_proxy", "aggressive_dedifferentiation_score"]
        axes = [c for c in axes if c in tcga]
        fig, ax = plt.subplots(figsize=(10, 7))
        heat = tcga[axes].apply(zscore_series).fillna(0).sample(min(250, len(tcga)), random_state=1)
        sns.heatmap(heat.T, cmap="vlag", center=0, cbar_kws={"label": "z-score"}, ax=ax)
        ax.set_title("Figure 2. TCGA Therapeutic Vulnerability Landscape")
        ax.set_xlabel("Representative TCGA-THCA tumors")
        ax.set_ylabel("Therapeutic axis")
        savefig(fig, prop / "figure2_tcga_therapeutic_vulnerability_landscape.png")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 6))
        sc = ax.scatter(tcga["rai_differentiation_score"], tcga["immune_visibility_score"], c=tcga.get("myeloid_caf_barrier_score", 0), s=35, cmap="Oranges", alpha=0.8)
        ax.set_xlabel("RAI differentiation")
        ax.set_ylabel("Immune visibility")
        ax.set_title("Figure 5. Therapeutic Quadrants")
        fig.colorbar(sc, ax=ax, label="Myeloid/CAF barrier")
        savefig(fig, prop / "figure5_therapeutic_quadrants_tcga.png")
        plt.close(fig)

    if not spots.empty:
        rep = spots["sample_id"].value_counts().index[0]
        if not spatial.empty and "same_niche_z" in spatial:
            rep = spatial.sort_values("same_niche_z", ascending=False).iloc[0]["sample_id"]
        s = spots[spots["sample_id"].eq(rep)]
        fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
        panels = [("hla_i_apm_score", "HLA/APM"), ("rai_differentiation_score", "RAI"), ("spatial_myeloid_caf_barrier_score", "CAF/myeloid"), ("spatial_drug_delivery_failure_proxy", "Delivery failure")]
        for ax, (col, title) in zip(axes, panels):
            ax.scatter(s["x"], s["y"], c=s[col], s=6, cmap="viridis", linewidths=0)
            ax.invert_yaxis()
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set_title(title)
        fig.suptitle(f"Figure 3. Spatial Therapeutic Niche Maps ({rep})")
        savefig(fig, prop / "figure3_spatial_therapeutic_niche_maps.png")
        plt.close(fig)

    if not spatial.empty:
        frac_cols = [c for c in spatial.columns if c.startswith("fraction_")]
        long = spatial.melt(id_vars=["sample_id", "condition"], value_vars=frac_cols, var_name="niche", value_name="fraction")
        long["niche"] = long["niche"].str.replace("fraction_", "", regex=False)
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.barplot(data=long, x="condition", y="fraction", hue="niche", estimator=np.mean, errorbar="se", ax=ax)
        ax.set_title("Figure 4. Spatial Niche Fractions by Disease Condition")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
        savefig(fig, prop / "figure4_spatial_coherence_statistics.png")
        plt.close(fig)

    if not drug.empty and "evidence_category" in drug:
        top = drug[drug["evidence_category"].astype(str).str.startswith(("A", "B", "C"))].head(18)
        if top.empty:
            top = drug.head(18)
        fig, ax = plt.subplots(figsize=(11, 5))
        plot = top[["drug_label", "mechanism_axis", "rank_score"]].copy()
        sns.barplot(data=plot, y="drug_label", x="rank_score", hue="mechanism_axis", dodge=False, ax=ax)
        ax.set_title("Figure 6. Drug/Perturbation Candidate Pilot")
        ax.set_xlabel("Exploratory rank score")
        ax.set_ylabel("")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, frameon=False)
        savefig(fig, prop / "figure6_drug_perturbation_candidate_pilot.png")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 5))
    ax.axis("off")
    steps = ["Computational inference", "FFPE mIHC/IF", "GeoMx ROI", "Fresh slice/organoid", "Functional assay"]
    markers = "PAX8/TG/TPO/NIS | HLA-I/B2M/TAP1 | CD8/GZMB | CD68/CD163 | ACTA2/FAP/COL1A1 | PD-L1 | CA9/VEGFA"
    for i, step in enumerate(steps):
        x = 0.05 + i * 0.19
        ax.add_patch(plt.Rectangle((x, 0.55), 0.15, 0.16, color=list(COLORS.values())[i % len(COLORS)]))
        ax.text(x + 0.075, 0.63, step, ha="center", va="center", color="white", weight="bold", fontsize=10)
        if i < len(steps) - 1:
            ax.annotate("", xy=(x + 0.18, 0.63), xytext=(x + 0.15, 0.63), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.35, markers, ha="center", fontsize=12)
    ax.text(0.5, 0.88, "Figure 7. Samsung Experimental Validation Design", ha="center", fontsize=18, weight="bold")
    savefig(fig, prop / "figure7_samsung_experimental_validation_design.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.axis("off")
    bullets = [
        "Public data support separable therapeutic vulnerability axes.",
        "Spatial data are needed because bulk average hides niche-level treatment logic.",
        "Experimental validation should focus on RAI restoration, HLA/APM visibility, CD8 exclusion, and drug-delivery failure.",
    ]
    ax.text(0.06, 0.82, "Figure 8. Pilot Conclusion", fontsize=24, weight="bold")
    for i, b in enumerate(bullets):
        ax.text(0.09, 0.62 - i * 0.18, f"{i+1}. {b}", fontsize=18)
    savefig(fig, prop / "figure8_pilot_conclusion_slide.png")
    plt.close(fig)
    logger.info("Proposal figures written to %s", prop)


if __name__ == "__main__":
    main()
