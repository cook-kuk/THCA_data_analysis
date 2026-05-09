#!/usr/bin/env python3
"""Publication-grade first-pass figures for CROSS-Neo 2.0."""

from __future__ import annotations

import shutil
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import OUT, ensure_dirs


def savefig(name: str) -> None:
    path = OUT / "figures" / name
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=240)
    plt.savefig(path.with_suffix(".pdf"))
    plt.close()


def box(ax, xy, text, fc="#f4f7fb"):
    ax.text(xy[0], xy[1], "\n".join(textwrap.wrap(text, 24)), ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.35", fc=fc, ec="#34495e", lw=1.2), fontsize=9)


def fig1_design() -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axis("off")
    ax.set_title("Internal locked contamination/source-shift audit design", fontsize=13)
    labels = [
        ("Canonical registry", (0.16, 0.72), "#e8f4f8"),
        ("Exact/near peptide-HLA overlap audit", (0.42, 0.72), "#fff4e6"),
        ("Strict locked splits", (0.68, 0.72), "#f2f0ff"),
        ("Counterfactual + PLM + structure features", (0.25, 0.35), "#eef7ee"),
        ("Nested train-only model selection", (0.55, 0.35), "#f8eeee"),
        ("Source-heldout failure/rescue report", (0.82, 0.35), "#fff7cc"),
    ]
    for txt, xy, color in labels:
        box(ax, xy, txt, color)
    for a, b in [((0.25, 0.72), (0.34, 0.72)), ((0.51, 0.72), (0.60, 0.72)), ((0.68, 0.64), (0.58, 0.43)), ((0.42, 0.64), (0.32, 0.43)), ((0.63, 0.35), (0.73, 0.35))]:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", lw=1.4))
    savefig("fig1_problem_and_audit_design")


def fig2_architecture() -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axis("off")
    ax.set_title("CROSS-Neo 2.0 internal architecture", fontsize=13)
    for txt, xy, color in [
        ("Mutant-WT counterfactual encoding", (0.18, 0.68), "#d8ecff"),
        ("Frozen PLM fallback embeddings", (0.18, 0.42), "#d8ecff"),
        ("Structure/processing/missingness", (0.18, 0.16), "#d8ecff"),
        ("Source-robust rankers", (0.50, 0.55), "#e8f5e9"),
        ("Rule MoE + abstention", (0.50, 0.25), "#e8f5e9"),
        ("Locked split report + case audit", (0.82, 0.40), "#fff3cd"),
    ]:
        box(ax, xy, txt, color)
    for y in [0.68, 0.42, 0.16]:
        ax.annotate("", xy=(0.40, 0.50 if y > 0.3 else 0.28), xytext=(0.28, y), arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(0.72, 0.40), xytext=(0.58, 0.43), arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(0.72, 0.40), xytext=(0.58, 0.28), arrowprops=dict(arrowstyle="->"))
    savefig("fig2_cross_neo_v2_architecture")


def fig3_locked(metrics: pd.DataFrame) -> None:
    primary = ["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"]
    sub = metrics[metrics["split_name"].isin(primary) & metrics["claim_status"].eq("reviewer_safe_internal_locked")].copy()
    if sub.empty:
        return
    best = sub.sort_values(["split_name", "AUPRC"], ascending=[True, False]).groupby("split_name").head(6)
    fig, ax = plt.subplots(figsize=(10, 5))
    xlabels, vals, colors = [], [], []
    for _, r in best.iterrows():
        xlabels.append(r["split_name"].replace("_", "\n") + "\n" + str(r["model_name"])[:18])
        vals.append(r["AUPRC"])
        colors.append("#2f6f9f" if str(r["model_name"]).startswith("v2") else "#7a7a7a")
    ax.bar(range(len(vals)), vals, color=colors)
    ax.set_ylabel("AUPRC")
    ax.set_title("Internal locked split performance: AUPRC")
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(xlabels, rotation=70, ha="right", fontsize=7)
    savefig("fig3_locked_split_performance")


def fig4_source(metrics: pd.DataFrame) -> None:
    sub = metrics[metrics["split_name"].astype(str).str.startswith("source_heldout_")].copy()
    if sub.empty:
        return
    show_models = ["anchor_rf", "v2_source_balanced_cf_rf", "v2_groupdro_proxy_cf_lr", "v2_rule_moe_cf_plm_structure"]
    sub = sub[sub["model_name"].isin(show_models)]
    piv = sub.pivot_table(index="split_name", columns="model_name", values="top20_precision", aggfunc="max").fillna(0)
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    im = ax.imshow(piv.values, aspect="auto", cmap="viridis", vmin=0, vmax=max(0.01, np.nanmax(piv.values)))
    ax.set_yticks(range(len(piv.index)))
    ax.set_yticklabels([x.replace("source_heldout_", "") for x in piv.index])
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels(piv.columns, rotation=35, ha="right", fontsize=8)
    ax.set_title("Source-heldout internal top20 recovery (descriptive small-n)")
    plt.colorbar(im, ax=ax, label="top20 precision")
    savefig("fig4_source_heldout_rescue")


def fig5_overlap() -> None:
    src = OUT / "overlap_heatmap.png"
    dst_png = OUT / "figures/fig5_public_overlap_audit.png"
    dst_pdf = OUT / "figures/fig5_public_overlap_audit.pdf"
    if src.exists():
        shutil.copyfile(src, dst_png)
        img = plt.imread(src)
        plt.figure(figsize=(8, 5))
        plt.imshow(img)
        plt.axis("off")
        plt.title("Internal/public overlap audit status")
        plt.savefig(dst_pdf)
        plt.close()
        return
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.axis("off")
    ax.set_title("Internal/public overlap audit status")
    ax.text(0.5, 0.5, "Overlap heatmap unavailable", ha="center")
    savefig("fig5_public_overlap_audit")


def fig6_abstention() -> None:
    path = OUT / "metrics/source_abstention_or_rescue_curve.tsv"
    if not path.exists():
        return
    curve = pd.read_csv(path, sep="\t")
    show = curve[curve["model_name"].isin(["anchor_rf", "v2_source_balanced_cf_rf", "v2_rule_moe_cf_plm_structure", "v2_cf_plm_rf"])]
    if show.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for model, g in show.groupby("model_name"):
        d = g.groupby("coverage")["top20_precision_after_score_abstention"].mean().sort_index()
        ax.plot(d.index, d.values, marker="o", label=model)
    ax.set_xlabel("coverage retained")
    ax.set_ylabel("source-heldout top20 precision")
    ax.set_title("Source-heldout internal abstention curve")
    ax.legend(fontsize=7)
    savefig("fig6_abstention_curve")


def fig7_case() -> None:
    rescue = pd.read_csv(OUT / "case_audit_rescued_positives.tsv", sep="\t") if (OUT / "case_audit_rescued_positives.tsv").exists() else pd.DataFrame()
    harm = pd.read_csv(OUT / "case_audit_harmed_negatives.tsv", sep="\t") if (OUT / "case_audit_harmed_negatives.tsv").exists() else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(["rescued positives", "harmed negatives"], [len(rescue), len(harm)], color=["#1f77b4", "#d62728"])
    ax.set_ylabel("case count")
    ax.set_title("Internal case-level rescue/harm audit")
    savefig("fig7_case_level_rescue_harm")


def fig8_claim_boundary() -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.axis("off")
    ax.set_title("Locked claim boundary decision tree", fontsize=13)
    items = [
        ("Public overlap resolved?", (0.20, 0.70), "#fff3cd"),
        ("Source-heldout rescued?", (0.50, 0.70), "#fff3cd"),
        ("Improves locked AUPRC/top-k?", (0.35, 0.38), "#d8ecff"),
        ("Allowed: robust internal method / benchmark", (0.70, 0.42), "#e8f5e9"),
        ("Forbidden now: external validation, clinical utility, quantum advantage", (0.50, 0.13), "#f8d7da"),
    ]
    for txt, xy, color in items:
        box(ax, xy, txt, color)
    for a, b in [((0.28, 0.70), (0.42, 0.70)), ((0.50, 0.62), (0.39, 0.46)), ((0.44, 0.38), (0.60, 0.42)), ((0.50, 0.30), (0.50, 0.20))]:
        ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->"))
    savefig("fig8_claim_boundary_decision_tree")


def main() -> None:
    ensure_dirs()
    metrics = pd.read_csv(OUT / "metrics/all_model_all_split_metrics.tsv", sep="\t")
    fig1_design()
    fig2_architecture()
    fig3_locked(metrics)
    fig4_source(metrics)
    fig5_overlap()
    fig6_abstention()
    fig7_case()
    fig8_claim_boundary()
    print("[v2-figures] generated figures 1-8 png/pdf where data available")


if __name__ == "__main__":
    main()
