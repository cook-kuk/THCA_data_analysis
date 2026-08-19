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
from scipy.stats import mannwhitneyu, wilcoxon

from pilot_utils import fdr_bh, pilot_root, setup_logging


RAI_SCORES = ["rai8_z_mean", "rai6_z_mean", "SLC5A5", "TG", "TPO", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]

MODULES = [
    "rai_differentiation_score",
    "hla_ii_apc_score",
    "cytotoxic_t_score",
    "ifn_activation_score",
    "checkpoint_suppression_score",
    "myeloid_tam_score",
    "caf_ecm_tgfb_score",
    "vascular_delivery_proxy_score",
    "hypoxia_score",
    "proliferation_score",
    "tumor_epithelial_score",
    "dm1_dark_lineage_score",
]


def cohens_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled) if np.isfinite(pooled) and pooled > 0 else np.nan


def rai_group_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    contexts = {
        "all_tumor": df[df["is_tumor"].astype(str).str.lower().isin(["true", "1"])].copy(),
        "primary_pre_rai_tumor": df[
            df["is_tumor"].astype(str).str.lower().isin(["true", "1"])
            & df["is_primary"].astype(str).str.lower().isin(["true", "1"])
            & df["is_pre_rai"].astype(str).str.lower().isin(["true", "1"])
        ].copy(),
        "metastatic_or_nonprimary_tumor": df[
            df["is_tumor"].astype(str).str.lower().isin(["true", "1"])
            & ~df["is_primary"].astype(str).str.lower().isin(["true", "1"])
        ].copy(),
    }
    comparisons = [
        ("response_refractory", 1, 0, "refractory_vs_avid_response"),
        ("uptake_no", 1, 0, "no_uptake_vs_uptake"),
    ]
    for context, cdf in contexts.items():
        for col, a_val, b_val, name in comparisons:
            if col not in cdf.columns:
                continue
            a_mask = pd.to_numeric(cdf[col], errors="coerce").eq(a_val)
            b_mask = pd.to_numeric(cdf[col], errors="coerce").eq(b_val)
            for score in [s for s in RAI_SCORES if s in cdf.columns]:
                a = pd.to_numeric(cdf.loc[a_mask, score], errors="coerce").dropna()
                b = pd.to_numeric(cdf.loc[b_mask, score], errors="coerce").dropna()
                p = mannwhitneyu(a, b, alternative="two-sided").pvalue if len(a) >= 2 and len(b) >= 2 else np.nan
                rows.append(
                    {
                        "context": context,
                        "comparison": name,
                        "score": score,
                        "n_group_a": len(a),
                        "n_group_b": len(b),
                        "mean_group_a": a.mean() if len(a) else np.nan,
                        "mean_group_b": b.mean() if len(b) else np.nan,
                        "effect_cohens_d_group_a_minus_b": cohens_d(a, b),
                        "mannwhitney_p": p,
                        "expected_direction": "lower in refractory/no-uptake expected for RAI-readiness scores",
                        "claim_boundary": "GSE151179 public RAI labels are noisy/small for direct response validation; use to define prospective assay need.",
                    }
                )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_q"] = fdr_bh(out["mannwhitney_p"])
        out["support"] = np.select(
            [
                out["effect_cohens_d_group_a_minus_b"].lt(-0.8) & out["fdr_q"].lt(0.10),
                out["effect_cohens_d_group_a_minus_b"].lt(-0.5),
                out["effect_cohens_d_group_a_minus_b"].ge(0.3),
            ],
            ["moderate_or_strong_expected", "weak_expected", "opposite_or_noisy"],
            default="weak_or_none",
        )
    return out


def rai_paired_tests(df: pd.DataFrame) -> pd.DataFrame:
    tumor = df[df["is_tumor"].astype(str).str.lower().isin(["true", "1"])].copy()
    rows = []
    if "patient_id" not in tumor.columns or "collection_before_after_rai" not in tumor.columns:
        return pd.DataFrame()
    tumor["collection_before_after_rai"] = tumor["collection_before_after_rai"].astype(str)
    for score in [s for s in RAI_SCORES if s in tumor.columns]:
        piv = tumor.pivot_table(index="patient_id", columns="collection_before_after_rai", values=score, aggfunc="mean")
        before_cols = [c for c in piv.columns if str(c).lower().startswith("before")]
        after_cols = [c for c in piv.columns if str(c).lower().startswith("after")]
        if not before_cols or not after_cols:
            continue
        paired = pd.DataFrame({"before": piv[before_cols[0]], "after": piv[after_cols[0]]}).dropna()
        if len(paired) < 3:
            p = np.nan
        else:
            try:
                p = wilcoxon(paired["after"], paired["before"]).pvalue
            except Exception:
                p = np.nan
        diff = paired["after"] - paired["before"]
        rows.append(
            {
                "score": score,
                "n_paired_patients": len(paired),
                "mean_before": paired["before"].mean() if len(paired) else np.nan,
                "mean_after": paired["after"].mean() if len(paired) else np.nan,
                "mean_after_minus_before": diff.mean() if len(diff) else np.nan,
                "wilcoxon_p": p,
                "claim_boundary": "Before/after RAI public samples are observational and not a controlled perturbation experiment.",
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_q"] = fdr_bh(out["wilcoxon_p"])
    return out


def scrna_attribution(scrna: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module in MODULES:
        col = f"{module}_mean"
        count_col = f"{module}_count"
        if col not in scrna.columns:
            continue
        sub = scrna[["cell_type", col, count_col if count_col in scrna.columns else "cell_type"]].copy()
        sub[col] = pd.to_numeric(sub[col], errors="coerce")
        if count_col in scrna.columns:
            sub[count_col] = pd.to_numeric(sub[count_col], errors="coerce")
        ranked = sub.dropna(subset=[col]).sort_values(col, ascending=False)
        if ranked.empty:
            rows.append({"module": module, "status": "not_scored", "reason": "all NaN or zero coverage"})
            continue
        top = ranked.iloc[0]
        second = ranked.iloc[1] if len(ranked) > 1 else None
        rows.append(
            {
                "module": module,
                "top_cell_type": top["cell_type"],
                "top_mean_score": top[col],
                "second_cell_type": second["cell_type"] if second is not None else np.nan,
                "second_mean_score": second[col] if second is not None else np.nan,
                "top_minus_second": top[col] - second[col] if second is not None else np.nan,
                "n_cells_top": top[count_col] if count_col in sub.columns else np.nan,
                "status": "scored",
                "claim_boundary": "scRNA cell-type attribution supports marker context; it does not deconvolve Visium spots or prove tumor-cell intrinsic expression.",
            }
        )
    return pd.DataFrame(rows)


def plot_rai(df: pd.DataFrame, fig_dir: Path) -> None:
    tumor = df[df["is_tumor"].astype(str).str.lower().isin(["true", "1"])].copy()
    tumor["response_label"] = np.where(pd.to_numeric(tumor["response_refractory"], errors="coerce").eq(1), "Refractory", "Avid/other")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    sns.boxplot(data=tumor, x="response_label", y="rai8_z_mean", color="#80b1d3", ax=axes[0])
    sns.stripplot(data=tumor, x="response_label", y="rai8_z_mean", color="#111111", size=3, ax=axes[0])
    axes[0].set_title("GSE151179 RAI module by response label", weight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("RAI8 z-mean")
    if "uptake_no" in tumor.columns:
        tumor["uptake_label"] = np.where(pd.to_numeric(tumor["uptake_no"], errors="coerce").eq(1), "No uptake", "Uptake/other")
        sns.boxplot(data=tumor, x="uptake_label", y="SLC5A5", color="#b3de69", ax=axes[1])
        sns.stripplot(data=tumor, x="uptake_label", y="SLC5A5", color="#111111", size=3, ax=axes[1])
        axes[1].set_title("SLC5A5 by metastatic uptake label", weight="bold")
        axes[1].set_xlabel("")
    fig.suptitle("Direct RAI public label reanalysis (exploratory/noisy)", fontsize=14, weight="bold")
    fig.tight_layout()
    fig.savefig(fig_dir / "gse151179_rai_direct_label_reanalysis.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "gse151179_rai_direct_label_reanalysis.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_scrna(scrna: pd.DataFrame, fig_dir: Path) -> None:
    cols = [f"{m}_mean" for m in MODULES if f"{m}_mean" in scrna.columns]
    mat = scrna.set_index("cell_type")[cols].copy()
    mat.columns = [c.replace("_score_mean", "").replace("_mean", "") for c in mat.columns]
    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    sns.heatmap(mat, cmap="vlag", center=0, linewidths=0.35, linecolor="#222", ax=ax)
    ax.set_title("scRNA cell-type attribution of K-Thyro modules", fontsize=14, weight="bold")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    fig.savefig(fig_dir / "scrna_module_celltype_attribution_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "scrna_module_celltype_attribution_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=pilot_root())
    parser.add_argument("--rai-path", type=Path, default=pilot_root().parent / "project/results/aggressive_sprint_2026_05_06/gse151179_rai_scores.tsv")
    args = parser.parse_args()
    logger = setup_logging("18_rai_scrna_deepdive")
    out = args.root / "results" / "extra_analyses"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    rai = pd.read_csv(args.rai_path, sep="\t") if args.rai_path.exists() else pd.DataFrame()
    if not rai.empty:
        group = rai_group_tests(rai)
        paired = rai_paired_tests(rai)
        group.to_csv(tables / "gse151179_rai_direct_label_reanalysis.tsv", sep="\t", index=False)
        paired.to_csv(tables / "gse151179_rai_before_after_paired_reanalysis.tsv", sep="\t", index=False)
        plot_rai(rai, figs)
        logger.info("Wrote GSE151179 RAI reanalysis: %d group rows, %d paired rows.", len(group), len(paired))
    else:
        pd.DataFrame([{"status": "not_run", "reason": f"missing {args.rai_path}"}]).to_csv(tables / "gse151179_rai_direct_label_reanalysis.tsv", sep="\t", index=False)

    scrna_path = args.root / "results" / "tables" / "scrna_celltype_module_scores.tsv"
    if scrna_path.exists():
        scrna = pd.read_csv(scrna_path, sep="\t")
        attr = scrna_attribution(scrna)
        attr.to_csv(tables / "scrna_module_celltype_attribution.tsv", sep="\t", index=False)
        plot_scrna(scrna, figs)
        logger.info("Wrote scRNA attribution: %d rows.", len(attr))


if __name__ == "__main__":
    main()
