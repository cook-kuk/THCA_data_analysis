#!/usr/bin/env python3
"""Generate CROSS-Neo v1 lockdown figures."""

from __future__ import annotations

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cross_neo_v1_lockdown_common import FIG, OUT, ensure_dirs


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIG / name, dpi=180)
    plt.close(fig)


def fig1_anchor_vs_fusion(comp: pd.DataFrame) -> None:
    keep_methods = ["anchor_rf", "anchor_lr", "qk_no_anchor", "qk_quantum_only", "prespecified_rf_qk_no_anchor_w0.5", "prespecified_rf_qk_quantum_only_w0.5", "rule_gate_rf_qk_fallback_train_selected"]
    d = comp[comp["method"].isin(keep_methods) & comp["split_name"].isin(["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"])].copy()
    fig, ax = plt.subplots(figsize=(11, 5))
    if len(d):
        piv = d.pivot_table(index="split_name", columns="method", values="AUPRC")
        piv.plot(kind="bar", ax=ax, width=0.85)
        for i, split in enumerate(piv.index):
            prev = float(d[d["split_name"].eq(split)]["prevalence"].iloc[0])
            ax.hlines(prev, i - 0.45, i + 0.45, colors="black", linestyles="dotted", linewidth=1)
    ax.set_title("Internal locked anchor vs fold-safe fusion AUPRC")
    ax.set_ylabel("AUPRC")
    ax.set_xlabel("")
    ax.legend(fontsize=7, ncol=2)
    save(fig, "fig1_anchor_vs_fusion_auprc.png")


def fig2_topk(comp: pd.DataFrame) -> None:
    d = comp[comp["split_name"].isin(["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"])].copy()
    d = d[d["method"].isin(["anchor_rf", "anchor_lr", "qk_no_anchor", "qk_quantum_only", "prespecified_rf_qk_no_anchor_w0.5", "prespecified_rf_qk_quantum_only_w0.5"])]
    fig, ax = plt.subplots(figsize=(11, 5))
    if len(d):
        d.pivot_table(index="split_name", columns="method", values="top10_precision").plot(kind="bar", ax=ax, width=0.85)
        for i, split in enumerate(sorted(d["split_name"].unique())):
            prev = float(d[d["split_name"].eq(split)]["prevalence"].iloc[0])
            ax.hlines(prev, i - 0.45, i + 0.45, colors="black", linestyles="dotted", linewidth=1)
    ax.set_title("Internal locked top-10 precision by split")
    ax.set_ylabel("Top-10 precision")
    ax.set_xlabel("")
    ax.legend(fontsize=7, ncol=2)
    save(fig, "fig2_topk_precision_by_split.png")


def fig3_rescue_harm() -> None:
    path = OUT / "qk_rescue_harm_summary.tsv"
    d = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(10, 5))
    if len(d):
        dd = d[d["event"].isin(["rescued_positive", "harmed_negative"])].copy()
        if len(dd):
            dd["key"] = dd["split_name"].str.replace("_holdout", "", regex=False) + "\n" + dd["comparator"]
            piv = dd.pivot_table(index="key", columns="event", values="n", aggfunc="sum", fill_value=0)
            piv.plot(kind="bar", ax=ax)
    ax.set_title("Internal locked QK rescue/harm counts")
    ax.set_ylabel("Case count")
    ax.set_xlabel("")
    save(fig, "fig3_qk_rescue_harm.png")


def fig4_source_ranks() -> None:
    path = OUT / "source_positive_rank_positions.tsv"
    d = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(10, 5))
    if len(d):
        use = d[d["method"].eq("sourceheld_prespecified_late_fusion_w0.5")]
        data = [use[use["heldout_study"].eq(s)]["rank_pct"].dropna().to_numpy() for s in sorted(use["heldout_study"].unique())]
        labels = sorted(use["heldout_study"].unique())
        if data:
            ax.boxplot(data, labels=labels, showfliers=False)
            ax.axhline(10 / max(1, use.groupby("heldout_study").size().median()), color="black", linestyle="dotted", linewidth=1)
    ax.set_title("Source-heldout descriptive positive rank positions")
    ax.set_ylabel("Positive rank percentile, lower is better")
    save(fig, "fig4_source_collapse_positive_ranks.png")


def fig5_shift() -> None:
    path = OUT / "source_feature_shift.tsv"
    d = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(9, 5))
    if len(d):
        dd = d[d["shift_metric"].isin(["categorical_total_variation", "standardized_mean_delta"])].copy()
        if len(dd):
            piv = dd.pivot_table(index="heldout_study", columns="feature", values="value", aggfunc="max").fillna(0)
            im = ax.imshow(piv.to_numpy(float), aspect="auto", cmap="viridis")
            ax.set_xticks(range(len(piv.columns)))
            ax.set_xticklabels(piv.columns, rotation=45, ha="right")
            ax.set_yticks(range(len(piv.index)))
            ax.set_yticklabels(piv.index)
            fig.colorbar(im, ax=ax, fraction=0.04, label="Shift magnitude")
    ax.set_title("Source-heldout descriptive feature shift heatmap")
    save(fig, "fig5_source_shift_heatmap.png")


def fig6_overlap() -> None:
    path = OUT / "public_overlap_manifest.tsv"
    d = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(9, 4))
    if len(d):
        counts = d["status"].value_counts()
        ax.bar(counts.index, counts.values, color=["#4c78a8", "#f58518", "#e45756"][: len(counts)])
        ax.tick_params(axis="x", rotation=30)
    ax.set_title("Locked public overlap audit status")
    ax.set_ylabel("Sources")
    save(fig, "fig6_overlap_audit_status.png")


def fig7_rescue_curve(comp: pd.DataFrame) -> None:
    d = comp[comp["method"].str.contains("abstention|rescue|diversification|rank_normalization|qk_compact", na=False)].copy()
    fig, ax = plt.subplots(figsize=(10, 5))
    if len(d):
        source = d[d["split_name"].str.startswith("source_heldout_")]
        for split, g in source.groupby("split_name"):
            g = g.sort_values("top10_precision", ascending=False).head(6)
            ax.plot(range(len(g)), g["top10_precision"], marker="o", label=split.replace("source_heldout_", ""))
            if len(g):
                ax.axhline(float(g["prevalence"].iloc[0]), color="gray", linestyle="dotted", linewidth=0.8)
        ax.legend(fontsize=8)
    ax.set_title("Source-heldout locked rescue/abstention top-10 curve")
    ax.set_ylabel("Top-10 precision")
    ax.set_xlabel("Top rescue variants")
    save(fig, "fig7_abstention_or_rescue_curve.png")


def main() -> None:
    ensure_dirs()
    comp = pd.read_csv(OUT / "v1_lockdown_comparison.tsv", sep="\t")
    fig1_anchor_vs_fusion(comp)
    fig2_topk(comp)
    fig3_rescue_harm()
    fig4_source_ranks()
    fig5_shift()
    fig6_overlap()
    fig7_rescue_curve(comp)
    print(f"[v1-figures] wrote {FIG}")


if __name__ == "__main__":
    main()
