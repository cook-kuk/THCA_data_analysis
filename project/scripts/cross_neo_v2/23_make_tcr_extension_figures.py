#!/usr/bin/env python3
"""Generate CROSS-Neo-TCR extension figures as PNG/PDF."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

from common import OUT


TCR_OUT = OUT / "tcr_extension"
FIG = TCR_OUT / "figures"

BG = "#fbfaf7"
INK = "#1f2933"
MUTED = "#637381"
BLUE = "#2f80ed"
TEAL = "#0f9f8f"
GOLD = "#c9952f"
RED = "#c43c39"
GREEN = "#2e7d32"
PURPLE = "#6a4c93"
GRAY = "#d9dee7"


def save(fig: plt.Figure, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def box(ax, xy, w, h, text, fc="white", ec=INK, color=INK, fs=10, lw=1.2, radius=0.045):
    patch = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle=f"round,pad=0.018,rounding_size={radius}",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", color=color, fontsize=fs, wrap=True)
    return patch


def arrow(ax, a, b, color=MUTED, lw=1.5, style="-|>", rad=0.0):
    ax.add_patch(
        FancyArrowPatch(
            a,
            b,
            arrowstyle=style,
            mutation_scale=13,
            linewidth=lw,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
        )
    )


def load_tables():
    registry_counts = pd.read_csv(TCR_OUT / "tcr_registry_source_counts.tsv", sep="\t")
    linkage = pd.read_csv(TCR_OUT / "tcr_neo_linkage_summary.tsv", sep="\t")
    metrics = pd.read_csv(TCR_OUT / "model_discovery/external_tcr_expert_benchmark/external_tcr_expert_benchmark_metrics.tsv", sep="\t")
    source_metrics = pd.read_csv(TCR_OUT / "model_discovery/external_tcr_expert_benchmark/external_tcr_expert_source_heldout_metrics.tsv", sep="\t")
    wetlab = pd.read_csv(TCR_OUT / "model_discovery/wetlab_external_tcr_expert_scores/wetlab_candidates_external_tcr_expert_ranked.tsv", sep="\t")
    jobs = pd.read_csv(TCR_OUT / "structure_jobs/tcr_structure_job_summary.tsv", sep="\t")
    tool_avail = pd.read_csv(TCR_OUT / "structure_jobs/tcr_structure_tool_availability.tsv", sep="\t")
    return registry_counts, linkage, metrics, source_metrics, wetlab, jobs, tool_avail


def fig_tcr1(registry_counts, linkage):
    total = int(registry_counts["n"].sum())
    paired = int(registry_counts["paired_tcr"].sum())
    peptide_hla = int(registry_counts["peptide_hla"].sum())
    exact = int(linkage.loc[linkage["tcr_link_category"].eq("exact_tcr_pmhc_match"), "n_neo_rows"].sum())
    no_match = int(linkage.loc[linkage["tcr_link_category"].eq("no_tcr_match"), "n_neo_rows"].sum())

    fig, ax = plt.subplots(figsize=(10, 5.7), facecolor=BG)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.05, 0.93, "Presentation is not recognition", fontsize=20, fontweight="bold", color=INK)
    ax.text(0.05, 0.88, "CROSS-Neo keeps pMHC ranking as the main model and adds TCR evidence only when available.", fontsize=10.5, color=MUTED)

    box(ax, (0.06, 0.55), 0.26, 0.15, "Peptide\nmutation-derived antigen", fc="#eef4ff", ec=BLUE, fs=11)
    box(ax, (0.39, 0.55), 0.26, 0.15, "HLA binding /\npresentation layer", fc="#edf7f4", ec=TEAL, fs=11)
    box(ax, (0.72, 0.55), 0.22, 0.15, "pMHC score\nmain ranker", fc="#fff6df", ec=GOLD, fs=11)
    arrow(ax, (0.32, 0.625), (0.39, 0.625), BLUE)
    arrow(ax, (0.65, 0.625), (0.72, 0.625), TEAL)

    box(ax, (0.39, 0.25), 0.26, 0.15, "TCR alpha/beta\nrecognition layer", fc="#f2edff", ec=PURPLE, fs=11)
    box(ax, (0.72, 0.25), 0.22, 0.15, "Optional expert\nwhen TCR exists", fc="#fef2f2", ec=RED, fs=11)
    arrow(ax, (0.52, 0.55), (0.52, 0.40), PURPLE)
    arrow(ax, (0.65, 0.325), (0.72, 0.325), PURPLE)

    stats = [
        ("TCR registry rows", f"{total:,}"),
        ("peptide-HLA rows", f"{peptide_hla:,}"),
        ("paired alpha/beta", f"{paired:,}"),
        ("CROSS-Neo exact paired links", f"{exact:,}"),
        ("CROSS-Neo no TCR match", f"{no_match:,}"),
    ]
    x0 = 0.06
    for i, (k, v) in enumerate(stats):
        x = x0 + i * 0.18
        ax.text(x, 0.12, v, fontsize=17, fontweight="bold", color=INK, ha="left")
        ax.text(x, 0.075, k, fontsize=8.2, color=MUTED, ha="left")
    save(fig, "fig_tcr1_recognition_gap")


def fig_tcr2():
    fig, ax = plt.subplots(figsize=(11, 6), facecolor=BG)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.05, 0.93, "CROSS-Neo-TCR architecture", fontsize=20, fontweight="bold", color=INK)
    ax.text(0.05, 0.88, "TCR-aware signals are gated by availability and confidence; missing TCR data does not create synthetic signal.", fontsize=10.5, color=MUTED)

    box(ax, (0.06, 0.66), 0.19, 0.12, "pMHC branch\npresentation score", fc="#eef4ff", ec=BLUE)
    box(ax, (0.06, 0.46), 0.19, 0.12, "Counterfactual branch\nmutant vs WT", fc="#edf7f4", ec=TEAL)
    box(ax, (0.06, 0.26), 0.19, 0.12, "Main immunogenicity\nranker", fc="#fff6df", ec=GOLD)
    box(ax, (0.35, 0.62), 0.21, 0.12, "TCR sequence branch\npMTnet + TEPCAM", fc="#f2edff", ec=PURPLE)
    box(ax, (0.35, 0.38), 0.21, 0.12, "Structure branch\ncontacts / geometry", fc="#f7f7f7", ec=MUTED)
    box(ax, (0.35, 0.17), 0.21, 0.12, "Missingness branch\nno TCR / beta-only / no structure", fc="#fef2f2", ec=RED)
    box(ax, (0.68, 0.42), 0.17, 0.16, "MoE gate\ninner-fold only", fc="white", ec=INK, fs=11)
    box(ax, (0.72, 0.16), 0.20, 0.12, "Wetlab priority\nscore + diagnostics", fc="#e9f7ef", ec=GREEN, fs=11)
    for y in [0.72, 0.52, 0.32]:
        arrow(ax, (0.25, y), (0.68, 0.50), BLUE if y == 0.72 else MUTED, rad=0.08)
    for y in [0.68, 0.44, 0.23]:
        arrow(ax, (0.56, y), (0.68, 0.50), PURPLE if y == 0.68 else MUTED, rad=-0.06)
    arrow(ax, (0.765, 0.42), (0.80, 0.28), GREEN)
    ax.text(0.68, 0.68, "TCR expert score only when TCR data are present", fontsize=9, color=PURPLE)
    save(fig, "fig_tcr2_architecture")


def fig_tcr3(jobs, tool_avail):
    fig, ax = plt.subplots(figsize=(11, 5.8), facecolor=BG)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.05, 0.93, "TCR-pMHC structure pipeline", fontsize=20, fontweight="bold", color=INK)
    ax.text(0.05, 0.88, "Current bottleneck is sequence completeness: full TCR and MHC chains are required before structure claims.", fontsize=10.5, color=MUTED)

    total_jobs = int(jobs["n_jobs"].sum()) if "n_jobs" in jobs else len(jobs)
    p0 = int(jobs.loc[jobs.get("priority", pd.Series(dtype=str)).astype(str).eq("P0"), "n_jobs"].sum()) if "priority" in jobs and "n_jobs" in jobs else 0
    tools = {}
    if {"tool", "available"}.issubset(tool_avail.columns):
        for r in tool_avail.itertuples(index=False):
            tools[str(r.tool)] = bool(r.available)

    stages = [
        ("Candidate\nPMHC/TCR", 0.08, 0.55, BLUE),
        ("Recover full\nTCR/MHC chains", 0.27, 0.55, TEAL),
        ("Run structure tool\nTCRdock / AF / Boltz / Chai", 0.48, 0.55, PURPLE),
        ("Parse contacts\nconfidence / geometry", 0.70, 0.55, GOLD),
        ("Diagnostic\nstructure features", 0.86, 0.55, GREEN),
    ]
    for text, x, y, color in stages:
        box(ax, (x - 0.07, y - 0.08), 0.14, 0.16, text, fc="white", ec=color, fs=9.5)
    for (_, x1, y1, _), (_, x2, y2, _) in zip(stages[:-1], stages[1:]):
        arrow(ax, (x1 + 0.07, y1), (x2 - 0.07, y2), MUTED)

    ax.text(0.09, 0.27, f"manifest jobs\n{total_jobs:,}", fontsize=16, fontweight="bold", color=INK, ha="center")
    ax.text(0.27, 0.27, f"P0 jobs\n{p0:,}", fontsize=16, fontweight="bold", color=INK, ha="center")
    available_count = sum(v for v in tools.values())
    ax.text(0.48, 0.27, f"available tools\n{available_count}", fontsize=16, fontweight="bold", color=INK, ha="center")
    ax.text(0.70, 0.27, "claim status\nQC/diagnostic", fontsize=16, fontweight="bold", color=INK, ha="center")
    ax.text(0.86, 0.27, "structure claim\nnot yet", fontsize=16, fontweight="bold", color=RED, ha="center")
    save(fig, "fig_tcr3_structure_pipeline")


def fig_tcr4():
    fig, ax = plt.subplots(figsize=(10, 5.7), facecolor=BG)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.05, 0.93, "Mutant vs WT interface delta", fontsize=20, fontweight="bold", color=INK)
    ax.text(0.05, 0.88, "Planned counterfactual structure readout: prioritize mutant recognition with low WT cross-reactivity risk.", fontsize=10.5, color=MUTED)

    for x, title, pep_color in [(0.24, "Mutant peptide", RED), (0.72, "Wildtype peptide", TEAL)]:
        ax.add_patch(Rectangle((x - 0.18, 0.22), 0.36, 0.08, facecolor="#e8edf3", edgecolor=MUTED, linewidth=1))
        ax.text(x, 0.26, "MHC groove", ha="center", va="center", fontsize=9, color=INK)
        ax.add_patch(Rectangle((x - 0.14, 0.34), 0.28, 0.045, facecolor=pep_color, edgecolor=pep_color))
        ax.text(x, 0.405, title, ha="center", fontsize=11, color=pep_color, fontweight="bold")
        ax.add_patch(Circle((x - 0.08, 0.64), 0.075, facecolor="#efe7ff", edgecolor=PURPLE, linewidth=1.2))
        ax.add_patch(Circle((x + 0.08, 0.64), 0.075, facecolor="#efe7ff", edgecolor=PURPLE, linewidth=1.2))
        ax.text(x - 0.08, 0.64, "TCRα", ha="center", va="center", fontsize=9, color=INK)
        ax.text(x + 0.08, 0.64, "TCRβ", ha="center", va="center", fontsize=9, color=INK)
        contacts = 7 if title.startswith("Mutant") else 3
        for i in range(contacts):
            xx = x - 0.11 + i * (0.22 / max(1, contacts - 1))
            ax.plot([xx, xx + 0.02 * math.sin(i)], [0.55, 0.39], color=GOLD if title.startswith("Mutant") else MUTED, linewidth=1.2, alpha=0.85)
        ax.text(x, 0.13, f"contact proxy: {contacts}\nconfidence-gated", ha="center", fontsize=9, color=MUTED)
    arrow(ax, (0.43, 0.52), (0.53, 0.52), INK, lw=1.4)
    ax.text(0.48, 0.57, "delta", ha="center", fontsize=11, color=INK, fontweight="bold")
    save(fig, "fig_tcr4_mutant_wt_interface_delta")


def fig_tcr5(metrics, source_metrics):
    pooled = metrics[metrics["panel"].eq("pooled")].copy()
    order = ["pMTnet", "TEPCAM", "mean_pMTnet_TEPCAM", "logreg_groupcv_pMTnet_TEPCAM"]
    pooled["model"] = pd.Categorical(pooled["model"], categories=order, ordered=True)
    pooled = pooled.sort_values("model")
    labels = ["pMTnet", "TEPCAM", "Mean\nensemble", "GroupCV\nensemble"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), facecolor=BG, gridspec_kw={"width_ratios": [1.1, 0.9]})
    ax = axes[0]
    x = np.arange(len(pooled))
    width = 0.36
    ax.bar(x - width / 2, pooled["auprc"], width, label="AUPRC", color=BLUE)
    ax.bar(x + width / 2, pooled["auroc"], width, label="AUROC", color=TEAL)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0.45, 0.86)
    ax.set_ylabel("Performance")
    ax.set_title("TCR expert benchmark, pooled", loc="left", fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax = axes[1]
    if not source_metrics.empty:
        y = np.arange(len(source_metrics))
        ax.barh(y - 0.18, source_metrics["auprc"], 0.34, color=GOLD, label="AUPRC")
        ax.barh(y + 0.18, source_metrics["auroc"], 0.34, color=PURPLE, label="AUROC")
        ax.set_yticks(y)
        ax.set_yticklabels(source_metrics["heldout_source"], fontsize=9)
        ax.invert_yaxis()
        ax.set_xlim(0.45, 0.90)
        ax.set_title("Source-heldout ensemble", loc="left", fontweight="bold")
        ax.legend(frameon=False)
        ax.grid(axis="x", alpha=0.2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.suptitle("TCR-available subset performance", x=0.04, y=1.03, ha="left", fontsize=20, fontweight="bold", color=INK)
    fig.text(0.04, 0.95, "Synthetic decoy benchmark for model selection, not an external SOTA claim.", color=MUTED, fontsize=10)
    save(fig, "fig_tcr5_tcr_available_subset_performance")


def fig_tcr6(wetlab):
    top = wetlab.head(10).copy()
    top["label"] = top["peptide"].astype(str) + "\n" + top["hla_4digit"].astype(str)
    fig, ax = plt.subplots(figsize=(11, 6), facecolor=BG)
    y = np.arange(len(top))
    old = pd.to_numeric(top["wetlab_priority_score_evidence_adjusted"], errors="coerce").fillna(0)
    new = pd.to_numeric(top["wetlab_priority_score_external_adjusted"], errors="coerce").fillna(old)
    ext = pd.to_numeric(top["external_tcr_expert_mean"], errors="coerce")
    ax.barh(y, old, color=GRAY, label="Evidence-adjusted priority")
    ax.barh(y, new, color=GREEN, alpha=0.78, label="External TCR-adjusted priority")
    ax.scatter(new, y, s=70, c=ext.fillna(0), cmap="viridis", vmin=0, vmax=max(0.8, ext.max(skipna=True) if ext.notna().any() else 0.8), edgecolor="white", linewidth=0.8, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(top["label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Wetlab priority score")
    ax.set_title("Case studies and wetlab candidate rescoring", loc="left", fontsize=18, fontweight="bold")
    ax.text(0.0, -1.05, "Color encodes mean pMTnet/TEPCAM score where exact TCR rows were modelable.", fontsize=9.5, color=MUTED)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", alpha=0.2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    save(fig, "fig_tcr6_case_studies")


def fig_tcr7():
    fig, ax = plt.subplots(figsize=(11, 6), facecolor=BG)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.05, 0.93, "Claim boundary", fontsize=20, fontweight="bold", color=INK)
    ax.text(0.05, 0.88, "TCR-aware evidence is useful, but incomplete paired TCR coverage limits the main claim.", fontsize=10.5, color=MUTED)
    allowed = [
        "Optional TCR expert for TCR-available rows",
        "Diagnostic explanation of selected false positives",
        "Wetlab prioritization evidence",
        "Missing paired TCR data is a benchmark limitation",
    ]
    forbidden = [
        "Universal TCR-aware neoantigen prediction",
        "Clinical utility",
        "Structure predictions are always correct",
        "TCR recognition from peptide-HLA alone",
    ]
    box(ax, (0.08, 0.18), 0.38, 0.58, "Allowed claims", fc="#ecf8ef", ec=GREEN, fs=14)
    box(ax, (0.54, 0.18), 0.38, 0.58, "Not allowed yet", fc="#fff1f0", ec=RED, fs=14)
    for i, text in enumerate(allowed):
        ax.text(0.12, 0.66 - i * 0.10, f"+ {text}", fontsize=10, color=INK, ha="left")
    for i, text in enumerate(forbidden):
        ax.text(0.58, 0.66 - i * 0.10, f"- {text}", fontsize=10, color=INK, ha="left")
    ax.text(0.50, 0.08, "Decision: diagnostic / wetlab tool now; main-method claim only after strict external paired-TCR benchmarks.", ha="center", fontsize=11, color=INK, fontweight="bold")
    save(fig, "fig_tcr7_claim_boundary")


def write_report() -> None:
    names = [
        "fig_tcr1_recognition_gap",
        "fig_tcr2_architecture",
        "fig_tcr3_structure_pipeline",
        "fig_tcr4_mutant_wt_interface_delta",
        "fig_tcr5_tcr_available_subset_performance",
        "fig_tcr6_case_studies",
        "fig_tcr7_claim_boundary",
    ]
    lines = ["# TCR Extension Figures", "", "| Figure | PNG | PDF |", "|---|---|---|"]
    for n in names:
        lines.append(f"| {n} | `{FIG / (n + '.png')}` | `{FIG / (n + '.pdf')}` |")
    (TCR_OUT / "tcr_extension_figure_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    registry_counts, linkage, metrics, source_metrics, wetlab, jobs, tool_avail = load_tables()
    fig_tcr1(registry_counts, linkage)
    fig_tcr2()
    fig_tcr3(jobs, tool_avail)
    fig_tcr4()
    fig_tcr5(metrics, source_metrics)
    fig_tcr6(wetlab)
    fig_tcr7()
    write_report()
    print(f"[tcr-figures] wrote 7 png/pdf figure pairs to {FIG}")


if __name__ == "__main__":
    main()
