#!/usr/bin/env python3
"""Generate ultra-priority visualization panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
FIG = OUT / "figures"
ULTRA = OUT / "ultra_priority"


def save(fig, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df = pd.read_csv(ULTRA / "ultra_wetlab_priority_candidates.tsv", sep="\t")
    top = df.head(12).copy()
    top["label"] = top["peptide"].astype(str) + "\n" + top["hla_4digit"].astype(str)

    fig, ax = plt.subplots(figsize=(10.5, 6))
    colors = top["recommendation_tier"].map(
        {
            "TIER_A_EXPERIMENT_NOW_WITH_CONTROLS": "#2fbf71",
            "TIER_A_PENDING_MD_COMPLETION_THEN_EXPERIMENT": "#f2c14e",
            "TIER_B_EXPERIMENT_AFTER_CONTROL_MD": "#4d96ff",
            "TIER_C_MD_OR_STRUCTURE_SCREEN_FIRST": "#9b8cff",
            "TIER_D_TCR_OR_WT_CURATION_FIRST": "#8a9aa9",
            "HOLD_FOR_NOW": "#56616f",
        }
    ).fillna("#56616f")
    ax.barh(top["label"][::-1], top["ultra_priority_score"][::-1], color=list(colors[::-1]))
    ax.set_xlabel("Ultra priority score (recommendation, not immunogenicity probability)")
    ax.set_title("CROSS-Neo-TCR-MD wetlab prioritization")
    ax.set_xlim(0, max(0.7, float(top["ultra_priority_score"].max()) + 0.08))
    ax.grid(axis="x", alpha=0.25)
    save(fig, "fig_md14_ultra_priority_ranking")

    blocks = [
        "model_evidence_score",
        "tcr_resource_score",
        "md_structural_score",
        "control_readiness_score",
    ]
    show = df[df["recommendation_tier"].ne("HOLD_FOR_NOW")].head(10).copy()
    if show.empty:
        show = top.head(10).copy()
    fig, ax = plt.subplots(figsize=(11, 6.2))
    bottom = pd.Series([0.0] * len(show), index=show.index)
    palette = {
        "model_evidence_score": "#4d96ff",
        "tcr_resource_score": "#2fbf71",
        "md_structural_score": "#f2c14e",
        "control_readiness_score": "#e76f51",
    }
    labels = show["peptide"].astype(str) + "\n" + show["hla_4digit"].astype(str)
    for col in blocks:
        vals = pd.to_numeric(show[col], errors="coerce").fillna(0) / len(blocks)
        ax.bar(labels, vals, bottom=bottom, label=col.replace("_score", "").replace("_", " "), color=palette[col])
        bottom += vals
    ax.scatter(labels, show["ultra_priority_score"], color="black", s=35, label="final score")
    ax.set_ylabel("Evidence block contribution")
    ax.set_title("Evidence composition for actionable candidates")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(frameon=False, ncols=2)
    ax.grid(axis="y", alpha=0.25)
    save(fig, "fig_md15_ultra_priority_evidence_blocks")

    (ULTRA / "ultra_priority_figure_report.md").write_text(
        "# Ultra Priority Figure Report\n\n"
        "- `fig_md14_ultra_priority_ranking.png/pdf`\n"
        "- `fig_md15_ultra_priority_evidence_blocks.png/pdf`\n"
    )
    print("[ultra-figures] generated fig_md14, fig_md15")


if __name__ == "__main__":
    main()
