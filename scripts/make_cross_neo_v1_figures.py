#!/usr/bin/env python3
"""Generate CROSS-Neo v1 figures."""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cross_neo_v1_common import FIG, V1, ensure_v1_dirs, metrics


def barh(df: pd.DataFrame, label: str, value: str, path, title: str, xlabel: str, baseline: float | None = None) -> None:
    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(df))), constrained_layout=True)
    ax.barh(df[label], df[value], color="#356c9b")
    if baseline is not None and np.isfinite(baseline):
        ax.axvline(baseline, color="#9b3d35", ls="--", lw=1.3, label="prevalence baseline")
        ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    ensure_v1_dirs()
    comp = pd.read_csv(V1 / "v1_model_comparison.tsv", sep="\t")
    hla = comp[(comp["split_name"] == "hla_stratified_group_5fold") & (comp["subset"] == "all_rows")].copy()
    hla["name"] = hla["family"] + " / " + hla["model"]
    top = hla.sort_values("AUPRC").tail(18)
    baseline = float(hla["prevalence"].dropna().iloc[0]) if len(hla) else None
    barh(top, "name", "AUPRC", FIG / "figure_v1_vs_v0_auprc.png", "CROSS-Neo v1 vs v0 Locked HLA-Stratified AUPRC", "AUPRC", baseline)
    barh(hla.sort_values("top10_precision").tail(18), "name", "top10_precision", FIG / "figure_v1_topk_precision.png", "CROSS-Neo v1 Locked HLA-Stratified Top-k Precision", "Top-10 precision", baseline)

    if (V1 / "branch_complementarity.tsv").exists():
        bc = pd.read_csv(V1 / "branch_complementarity.tsv", sep="\t")
        own = bc[~bc["expert"].astype(str).str.startswith("overlap::")].groupby("expert")["top10_precision"].mean().reset_index().sort_values("top10_precision")
        barh(own.tail(18), "expert", "top10_precision", FIG / "figure_branch_complementarity.png", "CROSS-Neo v1 Internal Branch Complementarity", "Mean top-10 precision", baseline)

    qk_rescue = pd.read_csv(V1 / "qk_rescue_cases.tsv", sep="\t") if (V1 / "qk_rescue_cases.tsv").exists() else pd.DataFrame()
    qk_harm = pd.read_csv(V1 / "qk_harm_cases.tsv", sep="\t") if (V1 / "qk_harm_cases.tsv").exists() else pd.DataFrame()
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    ax.bar(["QK rescue", "QK harm"], [len(qk_rescue), len(qk_harm)], color=["#397a54", "#a34b43"])
    ax.set_ylabel("Rows")
    ax.set_title("CROSS-Neo v1 Internal QK Rescue vs Harm")
    fig.savefig(FIG / "figure_qk_rescue_vs_harm.png", dpi=180)
    plt.close(fig)

    if (V1 / "source_shift_diagnostics.tsv").exists():
        sd = pd.read_csv(V1 / "source_shift_diagnostics.tsv", sep="\t")
        fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
        ax.bar(sd["source"], sd["prevalence"], color="#7a6a39")
        ax.set_ylabel("Label prevalence")
        ax.set_xticklabels(sd["source"], rotation=25, ha="right")
        ax.set_title("CROSS-Neo v1 Source-Heldout Prevalence Shift")
        fig.savefig(FIG / "figure_source_shift.png", dpi=180)
        plt.close(fig)

    if (V1 / "gated_moe_fold_weights.tsv").exists():
        w = pd.read_csv(V1 / "gated_moe_fold_weights.tsv", sep="\t")
        w = w[w["variant"].isin(["nested_selected_weights", "rule_gated_mean_weights"])].copy()
        rows = []
        for _, r in w.iterrows():
            try:
                ww = json.loads(r["weights"])
            except Exception:
                continue
            for k, v in ww.items():
                rows.append({"variant": r["variant"], "expert": k, "weight": float(v)})
        wd = pd.DataFrame(rows)
        if len(wd):
            piv = wd.groupby(["variant", "expert"])["weight"].mean().reset_index()
            fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
            for i, variant in enumerate(sorted(piv["variant"].unique())):
                sub = piv[piv["variant"] == variant]
                ax.bar(np.arange(len(sub)) + i * 0.38, sub["weight"], width=0.36, label=variant)
                labels = sub["expert"].tolist()
            ax.set_xticks(np.arange(len(labels)) + 0.18, labels, rotation=25, ha="right")
            ax.set_ylabel("Mean fold weight")
            ax.set_title("CROSS-Neo v1 Locked Gated MoE Weights")
            ax.legend()
            fig.savefig(FIG / "figure_gated_moe_weights.png", dpi=180)
            plt.close(fig)

    pred = pd.read_csv(V1 / "v1_all_model_predictions.tsv", sep="\t")
    best = hla[hla["family"].str.startswith("v1")].sort_values("AUPRC", ascending=False).head(1)
    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    if len(best):
        r = best.iloc[0]
        p = pred[(pred["family"] == r["family"]) & (pred["model"] == r["model"]) & (pred["split_name"] == "hla_stratified_group_5fold")].copy()
        p["confidence"] = np.abs(p["score"] - 0.5)
        xs, ys = [], []
        for cov in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]:
            kept = p.sort_values("confidence", ascending=False).head(max(1, int(round(len(p) * cov))))
            xs.append(cov)
            ys.append(metrics(kept["label"].to_numpy(), kept["score"].to_numpy())["top10_precision"])
        ax.plot(xs, ys, marker="o", label="top10 precision")
        ax.axhline(float(p["label"].mean()), color="#9b3d35", ls="--", label="prevalence")
    ax.set_xlabel("Coverage retained")
    ax.set_ylabel("Precision")
    ax.set_title("CROSS-Neo v1 Locked Abstention Coverage/Precision")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "figure_abstention_coverage_precision.png", dpi=180)
    plt.close(fig)
    print(f"[v1-figures] wrote {FIG}")


if __name__ == "__main__":
    main()
