#!/usr/bin/env python3
"""Generate CROSS-Neo v0 figures."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from cross_neo_v0_common import FIG, OUT, ensure_dirs


def save_bar(df, x, y, path, title, ylabel):
    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(df))), constrained_layout=True)
    ax.barh(df[x], df[y], color="#3b6ea8")
    ax.set_xlabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    metrics = pd.read_csv(OUT / "metrics_by_split.tsv", sep="\t")
    ab = metrics[metrics["split_name"] == "repeated_stratified_5x5_internal"].copy()
    ab["name"] = ab["feature_group"] + " / " + ab["model"]
    ab = ab.sort_values("AUPRC").tail(15)
    save_bar(
        ab,
        "name",
        "AUPRC",
        FIG / "figure_ablation_auprc.png",
        "CROSS-Neo v0 Internal Repeated 5x5 AUPRC Ablation",
        "AUPRC",
    )
    save_bar(
        ab.sort_values("top10_precision").tail(15),
        "name",
        "top10_precision",
        FIG / "figure_topk_precision.png",
        "CROSS-Neo v0 Internal Top-10 Precision",
        "Top-10 precision",
    )

    heat = metrics.pivot_table(index="feature_group", columns="split_name", values="AUPRC", aggfunc="max")
    fig, ax = plt.subplots(figsize=(12, max(5, 0.35 * len(heat))), constrained_layout=True)
    im = ax.imshow(heat.fillna(0).to_numpy(), aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(heat.index)), heat.index)
    ax.set_xticks(range(len(heat.columns)), heat.columns, rotation=35, ha="right")
    ax.set_title("CROSS-Neo v0 Locked/Internal Split Robustness Heatmap (AUPRC)")
    fig.colorbar(im, ax=ax, label="AUPRC")
    fig.savefig(FIG / "figure_split_robustness_heatmap.png", dpi=180)
    plt.close(fig)

    audit = pd.read_csv(OUT / "retrieval_leakage_audit.tsv", sep="\t")
    audit["name"] = audit["split_name"] + " / " + audit["retrieval_leakage_risk"]
    save_bar(
        audit.sort_values("n"),
        "name",
        "n",
        FIG / "figure_retrieval_leakage_flags.png",
        "CROSS-Neo v0 Train-Only Retrieval Leakage Flags",
        "Rows",
    )

    abst = pd.read_csv(OUT / "abstention_metrics.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(abst["coverage"], abst["top10_precision"], marker="o", label="top10 precision")
    ax.plot(abst["coverage"], abst["AUPRC"], marker="s", label="AUPRC")
    ax.set_xlabel("Coverage retained")
    ax.set_ylabel("Metric")
    ax.set_title("CROSS-Neo v0 Internal Calibration/Abstention Curve")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "figure_calibration_abstention.png", dpi=180)
    fig.savefig(FIG / "abstention_precision_curve.png", dpi=180)
    plt.close(fig)
    print(f"[figures] wrote {FIG}")


if __name__ == "__main__":
    main()
