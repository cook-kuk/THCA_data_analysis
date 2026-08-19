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
from scipy.stats import chi2_contingency

from pilot_utils import pilot_root, setup_logging


AXES = [
    "rai_differentiation_score",
    "hla_i_apm_score",
    "immune_visibility_score",
    "cd8_exclusion_proxy",
    "myeloid_caf_barrier_score",
    "drug_delivery_failure_proxy",
    "aggressive_dedifferentiation_score",
    "proliferation_score",
]

AXIS_LABELS = {
    "rai_differentiation_score": "RAI differentiation",
    "hla_i_apm_score": "HLA-I/APM",
    "immune_visibility_score": "Immune visibility",
    "cd8_exclusion_proxy": "CD8 exclusion",
    "myeloid_caf_barrier_score": "Myeloid/CAF barrier",
    "drug_delivery_failure_proxy": "Drug-delivery failure proxy",
    "aggressive_dedifferentiation_score": "Aggressive dedifferentiation",
    "proliferation_score": "Proliferation",
}


def shannon_entropy(counts: pd.Series) -> float:
    p = counts[counts > 0] / counts.sum()
    if len(p) == 0:
        return np.nan
    return float(-(p * np.log2(p)).sum())


def cramer_v(table: pd.DataFrame) -> tuple[float, float]:
    chi2, p, _, _ = chi2_contingency(table)
    n = table.to_numpy().sum()
    r, k = table.shape
    denom = n * (min(k - 1, r - 1))
    return (float(np.sqrt(chi2 / denom)) if denom > 0 else np.nan, float(p))


def eta_squared(df: pd.DataFrame, group_col: str, value_col: str) -> float:
    sub = df[[group_col, value_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if sub[group_col].nunique() < 2:
        return np.nan
    grand = sub[value_col].mean()
    ss_total = ((sub[value_col] - grand) ** 2).sum()
    ss_between = sub.groupby(group_col)[value_col].apply(lambda x: len(x) * (x.mean() - grand) ** 2).sum()
    return float(ss_between / ss_total) if ss_total > 0 else np.nan


def main() -> None:
    parser = argparse.ArgumentParser(description="Quantify why TCGA driver mutation alone is insufficient for K-Thyro vulnerability labels.")
    parser.add_argument("--root", type=Path, default=pilot_root())
    args = parser.parse_args()
    logger = setup_logging("15_tcga_mutation_not_enough_validation")

    tables = args.root / "results" / "public_validation_expansion" / "tables"
    figs = args.root / "results" / "public_validation_expansion" / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    path = args.root / "results" / "tables" / "tcga_thca_patient_vulnerability_scores.tsv"
    df = pd.read_csv(path, sep="\t")
    if "driver_anchor" not in df.columns or "primary_vulnerability_label" not in df.columns:
        raise RuntimeError("Required TCGA driver/label columns missing.")
    work = df.copy()
    work["driver_anchor_filled"] = work["driver_anchor"].fillna("No BRAF/RAS call")
    work["primary_vulnerability_label"] = work["primary_vulnerability_label"].fillna("Unlabeled")

    count_table = pd.crosstab(work["driver_anchor_filled"], work["primary_vulnerability_label"])
    frac_table = count_table.div(count_table.sum(axis=1), axis=0)
    count_table.to_csv(tables / "tcga_driver_by_vulnerability_label_counts.tsv", sep="\t")
    frac_table.to_csv(tables / "tcga_driver_by_vulnerability_label_fractions.tsv", sep="\t")

    cv, p_chi = cramer_v(count_table)
    diversity_rows = []
    for driver, counts in count_table.iterrows():
        n = int(counts.sum())
        largest = int(counts.max())
        diversity_rows.append(
            {
                "driver_group": driver,
                "n_patients": n,
                "n_nonzero_vulnerability_labels": int((counts > 0).sum()),
                "largest_label": counts.idxmax(),
                "largest_label_count": largest,
                "largest_label_fraction": largest / n if n else np.nan,
                "mutation_only_ambiguity_fraction": 1 - largest / n if n else np.nan,
                "shannon_entropy_bits": shannon_entropy(counts),
                "driver_label_cramers_v_overall": cv,
                "driver_label_chi_square_p_overall": p_chi,
                "claim_boundary": "TCGA mutation groups are public bulk molecular annotations; this quantifies heterogeneity within driver groups, not treatment-response prediction.",
            }
        )
    diversity = pd.DataFrame(diversity_rows)
    diversity.to_csv(tables / "tcga_driver_within_group_vulnerability_diversity.tsv", sep="\t", index=False)

    eta_rows = []
    for axis in AXES:
        if axis not in work.columns:
            continue
        eta = eta_squared(work, "driver_anchor_filled", axis)
        eta_known = eta_squared(work[work["driver_anchor"].notna()].copy(), "driver_anchor", axis)
        eta_rows.append(
            {
                "axis": axis,
                "axis_label": AXIS_LABELS.get(axis, axis),
                "eta_squared_driver_with_no_call_group": eta,
                "eta_squared_braf_vs_ras_known_only": eta_known,
                "unexplained_fraction_1_minus_eta_with_no_call_group": 1 - eta if np.isfinite(eta) else np.nan,
                "claim_boundary": "Eta-squared is descriptive; unmeasured fusions/CNV/methylation and purity can contribute.",
            }
        )
    eta_df = pd.DataFrame(eta_rows).sort_values("eta_squared_driver_with_no_call_group", ascending=False)
    eta_df.to_csv(tables / "tcga_driver_axis_variance_explained.tsv", sep="\t", index=False)

    plot_order = count_table.sum(axis=1).sort_values(ascending=False).index.tolist()
    label_order = count_table.sum(axis=0).sort_values(ascending=False).index.tolist()
    plot_frac = frac_table.loc[plot_order, label_order]
    fig, ax = plt.subplots(figsize=(13.5, 6.2))
    left = np.zeros(len(plot_frac))
    palette = sns.color_palette("tab10", n_colors=len(label_order))
    for color, label in zip(palette, label_order):
        vals = plot_frac[label].values
        ax.barh(plot_frac.index, vals, left=left, color=color, label=label)
        left += vals
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction within driver group")
    ax.set_title("TCGA driver groups split into multiple K-Thyro vulnerability labels", fontsize=15, weight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figs / "tcga_mutation_not_enough_driver_label_fraction.png", dpi=300, bbox_inches="tight")
    fig.savefig(figs / "tcga_mutation_not_enough_driver_label_fraction.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    sns.barplot(data=eta_df, x="eta_squared_driver_with_no_call_group", y="axis_label", color="#7bccc4", ax=ax)
    ax.set_xlim(0, max(0.35, eta_df["eta_squared_driver_with_no_call_group"].max() * 1.15))
    ax.set_xlabel("Eta-squared explained by BRAF/RAS/no-call group")
    ax.set_ylabel("")
    ax.set_title("Driver group explains only part of therapeutic-axis variance", fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(figs / "tcga_mutation_not_enough_axis_variance.png", dpi=300, bbox_inches="tight")
    fig.savefig(figs / "tcga_mutation_not_enough_axis_variance.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

    logger.info("Driver-label Cramer's V: %.3f, chi-square p: %.3g", cv, p_chi)
    logger.info("Output: %s", tables)


if __name__ == "__main__":
    main()
