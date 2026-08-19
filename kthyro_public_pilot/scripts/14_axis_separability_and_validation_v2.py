#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pilot_utils import pilot_root, setup_logging


AXES = {
    "RAI differentiation": ("rai_differentiation_score", "rai_differentiation_score"),
    "HLA/APM visibility": ("hla_i_apm_score", "hla_i_apm_score"),
    "CD8 cytotoxicity": ("cytotoxic_t_score", "cytotoxic_t_score"),
    "Myeloid/CAF barrier": ("myeloid_caf_barrier_score", "spatial_myeloid_caf_barrier_score"),
    "Drug-delivery failure proxy": ("drug_delivery_failure_proxy", "spatial_drug_delivery_failure_proxy"),
    "Aggressive dedifferentiation": ("aggressive_dedifferentiation_score", "spatial_aggressive_dedifferentiation_score"),
}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def axis_frame(df: pd.DataFrame, spatial: bool = False) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for label, (bulk_col, spatial_col) in AXES.items():
        col = spatial_col if spatial else bulk_col
        if col in df.columns:
            out[label] = pd.to_numeric(df[col], errors="coerce")
    return out.dropna(axis=1, how="all")


def corr_long(frame: pd.DataFrame, dataset: str, unit: str, n_units: int, sample_group: str = "") -> pd.DataFrame:
    if frame.shape[0] < 4 or frame.shape[1] < 2:
        return pd.DataFrame()
    corr = frame.corr(method="spearman", min_periods=4)
    rows = []
    cols = list(corr.columns)
    for i, a in enumerate(cols):
        for b in cols[i + 1 :]:
            val = corr.loc[a, b]
            rows.append(
                {
                    "dataset": dataset,
                    "sample_group": sample_group,
                    "unit": unit,
                    "n_units": n_units,
                    "axis_a": a,
                    "axis_b": b,
                    "spearman_rho": val,
                    "abs_rho": abs(val) if np.isfinite(val) else np.nan,
                    "interpretation": "partly_overlapping" if np.isfinite(val) and abs(val) >= 0.75 else "separable_or_moderately_related",
                    "claim_boundary": "Correlation checks axis separability; it does not prove causal independence or treatment response.",
                }
            )
    return pd.DataFrame(rows)


def compute_axis_separability(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    tables = root / "results" / "tables"
    expansion = root / "results" / "public_validation_expansion" / "tables"
    pieces = []
    heatmaps: dict[str, pd.DataFrame] = {}

    tcga = read_tsv(tables / "tcga_thca_patient_vulnerability_scores.tsv")
    if not tcga.empty:
        frame = axis_frame(tcga, spatial=False)
        pieces.append(corr_long(frame, "TCGA-THCA", "patient", len(frame), "primary_tumor"))
        heatmaps["TCGA patients"] = frame.corr(method="spearman", min_periods=10)

    external = read_tsv(expansion / "exact_kthyro_external_bulk_scores.tsv")
    if not external.empty:
        frame = axis_frame(external, spatial=False)
        pieces.append(corr_long(frame, "External GPL570 cohorts", "sample", len(frame), "all"))
        heatmaps["External bulk"] = frame.corr(method="spearman", min_periods=10)
        for dataset, ddf in external.groupby("dataset", sort=False):
            dframe = axis_frame(ddf, spatial=False)
            pieces.append(corr_long(dframe, dataset, "sample", len(dframe), "external_bulk"))

    spatial = read_tsv(tables / "spatial_spot_vulnerability_scores.tsv")
    spatial_corrs = []
    if not spatial.empty:
        for sample_id, sdf in spatial.groupby("sample_id", sort=False):
            sframe = axis_frame(sdf, spatial=True)
            one = corr_long(sframe, "GSE250521", "spot_within_slide", len(sframe), sample_id)
            if not one.empty:
                one["condition"] = str(sdf["condition"].iloc[0]) if "condition" in sdf.columns else "unknown"
                spatial_corrs.append(one)
                pieces.append(one)
        if spatial_corrs:
            all_sp = pd.concat(spatial_corrs, ignore_index=True)
            pivot = all_sp.pivot_table(index=["axis_a", "axis_b"], values="spearman_rho", aggfunc="median").reset_index()
            labels = list(AXES.keys())
            mat = pd.DataFrame(np.eye(len(labels)), index=labels, columns=labels)
            for _, r in pivot.iterrows():
                mat.loc[r["axis_a"], r["axis_b"]] = r["spearman_rho"]
                mat.loc[r["axis_b"], r["axis_a"]] = r["spearman_rho"]
            heatmaps["Spatial slide median"] = mat

    all_corr = pd.concat([p for p in pieces if p is not None and not p.empty], ignore_index=True) if pieces else pd.DataFrame()
    summary_rows = []
    if not all_corr.empty:
        for (dataset, unit, sample_group), ddf in all_corr.groupby(["dataset", "unit", "sample_group"], dropna=False, sort=False):
            summary_rows.append(
                {
                    "dataset": dataset,
                    "unit": unit,
                    "sample_group": sample_group,
                    "n_units": int(ddf["n_units"].max()),
                    "n_axis_pairs": int(ddf.shape[0]),
                    "median_abs_rho": float(ddf["abs_rho"].median()),
                    "max_abs_rho": float(ddf["abs_rho"].max()),
                    "n_pairs_abs_rho_ge_0_75": int((ddf["abs_rho"] >= 0.75).sum()),
                    "n_pairs_abs_rho_ge_0_90": int((ddf["abs_rho"] >= 0.90).sum()),
                    "interpretation": "axes show partial overlap but are not a single collapsed score"
                    if float(ddf["abs_rho"].median()) < 0.75
                    else "axes are highly coupled in this layer; interpret cautiously",
                }
            )
    return all_corr, pd.DataFrame(summary_rows), heatmaps


def plot_heatmaps(heatmaps: dict[str, pd.DataFrame], fig_dir: Path) -> None:
    if not heatmaps:
        return
    fig, axes = plt.subplots(1, len(heatmaps), figsize=(6.2 * len(heatmaps), 5.8), squeeze=False)
    for ax, (title, mat) in zip(axes.ravel(), heatmaps.items()):
        sns.heatmap(
            mat,
            cmap="vlag",
            center=0,
            vmin=-1,
            vmax=1,
            square=True,
            linewidths=0.4,
            linecolor="#222222",
            cbar=title == list(heatmaps.keys())[-1],
            cbar_kws={"label": "Spearman rho"} if title == list(heatmaps.keys())[-1] else None,
            ax=ax,
        )
        ax.set_title(title, fontsize=13, weight="bold")
        ax.tick_params(axis="x", rotation=45, labelsize=9)
        ax.tick_params(axis="y", rotation=0, labelsize=9)
    fig.suptitle("K-Thyro therapeutic axes are related but not collapsed into one signal", fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(fig_dir / "axis_separability_correlation_heatmaps.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "axis_separability_correlation_heatmaps.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def strength(effect: float, q: float, direction_match: str) -> str:
    if not np.isfinite(effect):
        return "not_tested"
    if direction_match == "opposite":
        return "opposite"
    ae = abs(effect)
    if ae >= 1.0 and np.isfinite(q) and q < 0.05:
        return "strong"
    if ae >= 0.75 and np.isfinite(q) and q < 0.10:
        return "moderate"
    if ae >= 0.5:
        return "weak"
    return "weak_or_none"


def build_validation_v2(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expansion = root / "results" / "public_validation_expansion" / "tables"
    base = read_tsv(expansion / "external_public_validation_signal_matrix.tsv")
    exact = read_tsv(expansion / "exact_kthyro_external_bulk_contrasts.tsv")
    exact_rows = []
    if not exact.empty:
        for _, r in exact.iterrows():
            exact_rows.append(
                {
                    "dataset": r["dataset"],
                    "validation_layer": "exact_KThyro_external_GPL570_bulk",
                    "axis": r["axis_label"],
                    "endpoint": r["contrast"],
                    "n_group_a": r["n_group_a"],
                    "n_group_b": r["n_group_b"],
                    "effect_size_cohens_d_or_reported": r["effect_size_cohens_d"],
                    "p_value": r["p_value"],
                    "direction": f"expected={r['expected_direction']}; observed={r['direction_match']}",
                    "support_strength": strength(r["effect_size_cohens_d"], r["fdr_q_value"], r["direction_match"]),
                    "caveat": r["claim_boundary"],
                    "source_file": "results/public_validation_expansion/tables/exact_kthyro_external_bulk_contrasts.tsv",
                }
            )
    exact_df = pd.DataFrame(exact_rows)
    if base.empty:
        v2 = exact_df
    elif exact_df.empty:
        v2 = base
    else:
        v2 = pd.concat([base, exact_df], ignore_index=True, sort=False)
    if v2.empty:
        return v2, pd.DataFrame()
    summary = (
        v2.groupby(["dataset", "validation_layer", "support_strength"], dropna=False)
        .size()
        .reset_index(name="n_tests")
        .sort_values(["validation_layer", "dataset", "support_strength"])
    )
    return v2, summary


def plot_support_counts(summary: pd.DataFrame, fig_dir: Path) -> None:
    if summary.empty:
        return
    order = ["strong", "moderate", "weak", "weak_or_none", "opposite", "not_tested"]
    df = summary.copy()
    df["layer_dataset"] = df["validation_layer"] + "\n" + df["dataset"]
    pivot = df.pivot_table(index="layer_dataset", columns="support_strength", values="n_tests", aggfunc="sum", fill_value=0)
    pivot = pivot[[c for c in order if c in pivot.columns]]
    fig_h = max(6, 0.28 * len(pivot) + 2)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    colors = {
        "strong": "#2ca25f",
        "moderate": "#99d8c9",
        "weak": "#fed976",
        "weak_or_none": "#bdbdbd",
        "opposite": "#de2d26",
        "not_tested": "#636363",
    }
    left = np.zeros(len(pivot))
    y = np.arange(len(pivot))
    for col in pivot.columns:
        ax.barh(y, pivot[col].values, left=left, label=col, color=colors.get(col, "#cccccc"))
        left += pivot[col].values
    ax.set_yticks(y)
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Number of validation tests")
    ax.set_title("Expanded public validation matrix v2", fontsize=15, weight="bold")
    ax.legend(frameon=False, ncol=3, loc="lower right")
    fig.tight_layout()
    fig.savefig(fig_dir / "public_validation_expansion_v2_support_counts.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "public_validation_expansion_v2_support_counts.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize axis separability and combine validation expansion v2.")
    parser.add_argument("--root", type=Path, default=pilot_root())
    args = parser.parse_args()
    logger = setup_logging("14_axis_separability_and_validation_v2")
    out = args.root / "results" / "public_validation_expansion"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    corr, corr_summary, heatmaps = compute_axis_separability(args.root)
    corr.to_csv(tables / "axis_separability_pairwise_correlations.tsv", sep="\t", index=False)
    corr_summary.to_csv(tables / "axis_separability_dataset_summary.tsv", sep="\t", index=False)
    plot_heatmaps(heatmaps, figs)

    v2, v2_summary = build_validation_v2(args.root)
    v2.to_csv(tables / "external_public_validation_signal_matrix_v2.tsv", sep="\t", index=False)
    v2_summary.to_csv(tables / "external_public_validation_layer_summary_v2.tsv", sep="\t", index=False)
    plot_support_counts(v2_summary, figs)

    logger.info("Axis separability pairs: %d", len(corr))
    logger.info("Validation v2 tests: %d", len(v2))
    logger.info("Output: %s", out)


if __name__ == "__main__":
    main()
