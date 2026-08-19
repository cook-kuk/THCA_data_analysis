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
from scipy.stats import mannwhitneyu

from pilot_utils import (
    add_derived_axes,
    clean_gene_symbol,
    compute_module_scores,
    fdr_bh,
    gene_sets,
    pilot_root,
    setup_logging,
)


KEY_AXES = [
    "rai_differentiation_score",
    "rai_restorable_score",
    "mapk_score",
    "hla_i_apm_score",
    "immune_visibility_score",
    "hla_low_invisible_score",
    "cytotoxic_t_score",
    "cd8_exclusion_proxy",
    "myeloid_caf_barrier_score",
    "drug_delivery_failure_proxy",
    "hypoxia_score",
    "proliferation_score",
    "aggressive_dedifferentiation_score",
    "tumor_epithelial_score",
]

AXIS_LABELS = {
    "rai_differentiation_score": "RAI differentiation",
    "rai_restorable_score": "RAI-restorable proxy",
    "mapk_score": "MAPK stress",
    "hla_i_apm_score": "HLA-I/APM",
    "immune_visibility_score": "Immune visibility",
    "hla_low_invisible_score": "HLA-low invisible",
    "cytotoxic_t_score": "CD8/cytotoxic",
    "cd8_exclusion_proxy": "CD8 exclusion proxy",
    "myeloid_caf_barrier_score": "Myeloid/CAF barrier",
    "drug_delivery_failure_proxy": "Drug-delivery failure proxy",
    "hypoxia_score": "Hypoxia",
    "proliferation_score": "Proliferation",
    "aggressive_dedifferentiation_score": "Aggressive dedifferentiation",
    "tumor_epithelial_score": "Tumor epithelial/thyroid lineage",
}

CONTRASTS = [
    ("ATC", "PTC"),
    ("ATC", "normal"),
    ("PTC", "normal"),
    ("PDTC", "PTC"),
    ("PDTC", "normal"),
    ("FVPTC", "PTC"),
    ("FTC", "normal"),
]

EXPECTED_ADVANCED = {
    "ATC_vs_PTC": {
        "rai_differentiation_score": "down",
        "rai_restorable_score": "down",
        "aggressive_dedifferentiation_score": "up",
        "proliferation_score": "up",
        "hypoxia_score": "up",
        "drug_delivery_failure_proxy": "up",
    },
    "ATC_vs_normal": {
        "rai_differentiation_score": "down",
        "rai_restorable_score": "down",
        "aggressive_dedifferentiation_score": "up",
        "proliferation_score": "up",
        "hypoxia_score": "up",
        "drug_delivery_failure_proxy": "up",
    },
    "PTC_vs_normal": {
        "rai_differentiation_score": "down",
        "rai_restorable_score": "down",
    },
    "PDTC_vs_PTC": {
        "rai_differentiation_score": "down",
        "rai_restorable_score": "down",
        "aggressive_dedifferentiation_score": "up",
        "proliferation_score": "up",
    },
    "PDTC_vs_normal": {
        "rai_differentiation_score": "down",
        "rai_restorable_score": "down",
        "aggressive_dedifferentiation_score": "up",
        "proliferation_score": "up",
    },
}


def cohens_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    denom = len(a) + len(b) - 2
    if denom <= 0:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / denom)
    if not np.isfinite(pooled) or pooled == 0:
        return np.nan
    return float((a.mean() - b.mean()) / pooled)


def mw_p(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    if a.nunique() < 2 and b.nunique() < 2:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def direction_match(effect: float, expected: str | None) -> str:
    if expected is None or not np.isfinite(effect):
        return "not_prespecified"
    if expected == "up":
        return "match" if effect > 0 else "opposite"
    if expected == "down":
        return "match" if effect < 0 else "opposite"
    return "not_prespecified"


def support_strength(effect: float, p: float, q: float, match: str) -> str:
    if not np.isfinite(effect):
        return "not_tested"
    if match == "opposite":
        return "opposite"
    p_use = q if np.isfinite(q) else p
    ae = abs(effect)
    if ae >= 1.0 and np.isfinite(p_use) and p_use < 0.05:
        return "strong"
    if ae >= 0.75 and np.isfinite(p_use) and p_use < 0.10:
        return "moderate"
    if ae >= 0.5:
        return "weak"
    return "weak_or_none"


def read_expression(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", compression="gzip")
    if "gene_symbol" not in raw.columns:
        raise ValueError(f"gene_symbol column missing in {path}")
    raw["gene_symbol"] = raw["gene_symbol"].map(clean_gene_symbol)
    raw = raw[raw["gene_symbol"].ne("")]
    value_cols = [c for c in raw.columns if c != "gene_symbol"]
    values = raw[value_cols].apply(pd.to_numeric, errors="coerce")
    row_mean = values.mean(axis=1)
    best_idx = row_mean.groupby(raw["gene_symbol"], sort=False).idxmax()
    raw = raw.loc[best_idx].copy()
    expr = raw.set_index("gene_symbol").T
    expr.index.name = "sample_id"
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr


def score_external_datasets(input_dir: Path, meta: pd.DataFrame, logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_scores = []
    all_coverage = []
    gs = gene_sets()
    for expr_path in sorted(input_dir.glob("GSE*_expression_gene_log.tsv.gz")):
        dataset = expr_path.name.split("_", 1)[0]
        logger.info("Scoring %s from %s", dataset, expr_path)
        expr = read_expression(expr_path)
        dmeta = meta[meta["dataset"].eq(dataset)].copy()
        keep = [s for s in expr.index if s in set(dmeta["sample_id"])]
        if not keep:
            logger.warning("No metadata-matched samples for %s", dataset)
            continue
        expr = expr.loc[keep]
        scores, coverage = compute_module_scores(expr, gs, method="mean_z")
        scores = add_derived_axes(scores)
        scores.index.name = None
        scores.insert(0, "sample_id", scores.index)
        scores.insert(0, "dataset", dataset)
        merged = dmeta.merge(scores, on=["dataset", "sample_id"], how="inner")
        all_scores.append(merged)
        coverage.insert(0, "dataset", dataset)
        all_coverage.append(coverage)
    if not all_scores:
        return pd.DataFrame(), pd.DataFrame()
    return pd.concat(all_scores, ignore_index=True), pd.concat(all_coverage, ignore_index=True)


def make_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, ddf in scores.groupby("dataset", sort=False):
        for group_a, group_b in CONTRASTS:
            a = ddf[ddf["histology_clean"].eq(group_a)]
            b = ddf[ddf["histology_clean"].eq(group_b)]
            if len(a) < 2 or len(b) < 2:
                continue
            contrast = f"{group_a}_vs_{group_b}"
            expected_map = EXPECTED_ADVANCED.get(contrast, {})
            for axis in KEY_AXES:
                if axis not in ddf.columns:
                    continue
                effect = cohens_d(a[axis], b[axis])
                p = mw_p(a[axis], b[axis])
                expected = expected_map.get(axis)
                rows.append(
                    {
                        "dataset": dataset,
                        "contrast": contrast,
                        "group_a": group_a,
                        "group_b": group_b,
                        "n_group_a": len(a),
                        "n_group_b": len(b),
                        "axis": axis,
                        "axis_label": AXIS_LABELS.get(axis, axis),
                        "mean_group_a": pd.to_numeric(a[axis], errors="coerce").mean(),
                        "mean_group_b": pd.to_numeric(b[axis], errors="coerce").mean(),
                        "effect_size_cohens_d": effect,
                        "p_value": p,
                        "expected_direction": expected or "not_prespecified",
                        "direction_match": direction_match(effect, expected),
                        "claim_boundary": "External GPL570 bulk expression validates axis directionality/separability, not clinical treatment response or spatial niche biology.",
                    }
                )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr_q_value"] = fdr_bh(out["p_value"])
    out["support_strength"] = [
        support_strength(r.effect_size_cohens_d, r.p_value, r.fdr_q_value, r.direction_match)
        for r in out.itertuples(index=False)
    ]
    return out.sort_values(["dataset", "contrast", "axis"]).reset_index(drop=True)


def write_summary(scores: pd.DataFrame, contrasts: pd.DataFrame, out_path: Path) -> None:
    rows = []
    for dataset, ddf in scores.groupby("dataset", sort=False):
        row = {
            "dataset": dataset,
            "n_samples": len(ddf),
            "histology_counts": "; ".join(f"{k}:{v}" for k, v in ddf["histology_clean"].value_counts().items()),
        }
        sub = contrasts[contrasts["dataset"].eq(dataset)] if not contrasts.empty else pd.DataFrame()
        row["n_contrasts_tested"] = int(sub[["contrast", "axis"]].drop_duplicates().shape[0]) if not sub.empty else 0
        row["n_strong"] = int(sub["support_strength"].eq("strong").sum()) if not sub.empty else 0
        row["n_moderate"] = int(sub["support_strength"].eq("moderate").sum()) if not sub.empty else 0
        row["n_opposite"] = int(sub["support_strength"].eq("opposite").sum()) if not sub.empty else 0
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def plot_heatmap(contrasts: pd.DataFrame, fig_dir: Path) -> None:
    if contrasts.empty:
        return
    selected = contrasts[contrasts["axis"].isin(KEY_AXES)].copy()
    selected["row"] = selected["dataset"] + " | " + selected["contrast"]
    pivot = selected.pivot_table(index="row", columns="axis_label", values="effect_size_cohens_d", aggfunc="mean")
    ordered_cols = [AXIS_LABELS[a] for a in KEY_AXES if AXIS_LABELS.get(a, a) in pivot.columns]
    pivot = pivot.loc[:, ordered_cols]
    fig_h = max(5.5, 0.38 * len(pivot) + 1.6)
    fig, ax = plt.subplots(figsize=(15.5, fig_h))
    sns.heatmap(
        pivot,
        cmap="vlag",
        center=0,
        vmin=-3,
        vmax=3,
        linewidths=0.35,
        linecolor="#222222",
        cbar_kws={"label": "Cohen's d (group A - group B)"},
        ax=ax,
    )
    ax.set_title("Exact K-Thyro axes in external GPL570 thyroid cohorts", fontsize=15, weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "exact_kthyro_external_bulk_axis_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "exact_kthyro_external_bulk_axis_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_key_boxplots(scores: pd.DataFrame, fig_dir: Path) -> None:
    if scores.empty:
        return
    axes = [
        "rai_differentiation_score",
        "immune_visibility_score",
        "myeloid_caf_barrier_score",
        "drug_delivery_failure_proxy",
        "aggressive_dedifferentiation_score",
    ]
    plot_df = scores.melt(
        id_vars=["dataset", "sample_id", "histology_clean"],
        value_vars=[a for a in axes if a in scores.columns],
        var_name="axis",
        value_name="score",
    )
    plot_df["axis"] = plot_df["axis"].map(lambda x: AXIS_LABELS.get(x, x))
    hist_order = ["normal", "PTC", "FVPTC", "FTC", "PDTC", "ATC"]
    fig, axes_arr = plt.subplots(2, 3, figsize=(17.5, 9.5), sharey=False)
    axes_arr = axes_arr.ravel()
    for ax, axis_name in zip(axes_arr, [AXIS_LABELS[a] for a in axes]):
        sub = plot_df[plot_df["axis"].eq(axis_name)]
        sns.boxplot(data=sub, x="histology_clean", y="score", order=hist_order, color="#9ecae1", fliersize=2, ax=ax)
        sns.stripplot(data=sub, x="histology_clean", y="score", order=hist_order, hue="dataset", dodge=True, size=2.5, alpha=0.65, ax=ax)
        ax.set_title(axis_name, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("mean-z module score")
        ax.tick_params(axis="x", rotation=35)
        if ax.legend_:
            ax.legend_.remove()
    axes_arr[-1].axis("off")
    handles, labels = axes_arr[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="lower right", frameon=False, ncol=2)
    fig.suptitle("External public cohorts reproduce K-Thyro axis separability", fontsize=16, weight="bold")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(fig_dir / "exact_kthyro_external_bulk_boxplots.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "exact_kthyro_external_bulk_boxplots.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Re-score local external GPL570 thyroid cohorts with exact K-Thyro gene sets.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=pilot_root().parent / "project" / "results" / "p_external_expression_validation",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=pilot_root() / "results" / "public_validation_expansion",
    )
    args = parser.parse_args()

    logger = setup_logging("13_exact_kthyro_external_bulk_validation")
    tables = args.output_dir / "tables"
    figs = args.output_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    meta_path = args.input_dir / "sample_metadata.tsv"
    if not meta_path.exists():
        raise FileNotFoundError(f"Missing metadata: {meta_path}")
    meta = pd.read_csv(meta_path, sep="\t")
    scores, coverage = score_external_datasets(args.input_dir, meta, logger)
    if scores.empty:
        raise RuntimeError("No external cohorts could be scored.")

    contrasts = make_contrasts(scores)
    scores.to_csv(tables / "exact_kthyro_external_bulk_scores.tsv", sep="\t", index=False)
    coverage.to_csv(tables / "exact_kthyro_external_bulk_gene_coverage.tsv", sep="\t", index=False)
    contrasts.to_csv(tables / "exact_kthyro_external_bulk_contrasts.tsv", sep="\t", index=False)
    write_summary(scores, contrasts, tables / "exact_kthyro_external_bulk_dataset_summary.tsv")
    plot_heatmap(contrasts, figs)
    plot_key_boxplots(scores, figs)

    logger.info("Wrote %d scored samples and %d contrast tests.", len(scores), len(contrasts))
    logger.info("Output: %s", args.output_dir)


if __name__ == "__main__":
    main()
