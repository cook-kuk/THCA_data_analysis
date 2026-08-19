#!/usr/bin/env python3
"""Hotspot concordance for Path2Space-inspired image-to-spatial-RNA evidence."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
GSE250 = P2 / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"
GSE230 = P2 / "analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"
GSE230_RESID = GSE230 / "residual_target_controls"
OUT = P2 / "analysis_supp/path2space_hotspot_concordance_2026_05_09"
QUANTILES = [0.80, 0.90]
N_PERM = 1000
RANDOM_SEED = 250230425


def top_mask_within_group(df: pd.DataFrame, value_col: str, group_col: str, q: float) -> np.ndarray:
    mask = np.zeros(len(df), dtype=bool)
    values = df[value_col].to_numpy(dtype=float)
    groups = df[group_col].astype(str).to_numpy()
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        if keep.sum() < 5:
            continue
        valid_idx = idx[keep]
        n_top = max(1, int(np.ceil((1.0 - q) * len(valid_idx))))
        # Exact top-k avoids inflated hotspot rates when smoothed labels have ties.
        order = np.argsort(vals[keep])[::-1]
        mask[valid_idx[order[:n_top]]] = True
    return mask


def safe_fisher(tp: int, fp: int, fn: int, tn: int) -> tuple[float, float]:
    if min(tp, fp, fn, tn) < 0:
        return np.nan, np.nan
    odds, p = stats.fisher_exact([[tp, fp], [fn, tn]], alternative="greater")
    return float(odds), float(p)


def permute_overlap(pred_mask: np.ndarray, obs_mask: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> float:
    observed_tp = int(np.sum(pred_mask & obs_mask))
    null = np.zeros(N_PERM, dtype=int)
    for b in range(N_PERM):
        perm = pred_mask.copy()
        for group in pd.unique(groups):
            idx = np.where(groups == group)[0]
            perm[idx] = rng.permutation(perm[idx])
        null[b] = int(np.sum(perm & obs_mask))
    return float((1 + np.sum(null >= observed_tp)) / (1 + N_PERM))


def evaluate(df: pd.DataFrame, evidence: str, group_col: str, obs_col: str, pred_col: str, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    sample_rows = []
    groups = df[group_col].astype(str).to_numpy()
    for q in QUANTILES:
        obs_hot = top_mask_within_group(df, obs_col, group_col, q)
        pred_hot = top_mask_within_group(df, pred_col, group_col, q)
        finite = np.isfinite(df[obs_col].to_numpy(dtype=float)) & np.isfinite(df[pred_col].to_numpy(dtype=float))
        obs_hot &= finite
        pred_hot &= finite

        tp = int(np.sum(obs_hot & pred_hot))
        fp = int(np.sum(~obs_hot & pred_hot & finite))
        fn = int(np.sum(obs_hot & ~pred_hot & finite))
        tn = int(np.sum(~obs_hot & ~pred_hot & finite))
        n = int(np.sum(finite))
        odds, fisher_p = safe_fisher(tp, fp, fn, tn)
        precision = tp / max(tp + fp, 1)
        recall = tp / max(tp + fn, 1)
        obs_rate = (tp + fn) / max(n, 1)
        pred_rate = (tp + fp) / max(n, 1)
        expected_tp = obs_rate * pred_rate * n
        perm_p = permute_overlap(pred_hot[finite], obs_hot[finite], groups[finite], rng)
        rows.append(
            {
                "evidence": evidence,
                "quantile": q,
                "top_percent": int(round((1 - q) * 100)),
                "n": n,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "tn": tn,
                "obs_hotspot_rate": obs_rate,
                "pred_hotspot_rate": pred_rate,
                "precision": precision,
                "recall": recall,
                "precision_lift": precision / obs_rate if obs_rate > 0 else np.nan,
                "jaccard": tp / max(tp + fp + fn, 1),
                "expected_tp_independent": expected_tp,
                "tp_over_expected": tp / expected_tp if expected_tp > 0 else np.nan,
                "fisher_odds_ratio": odds,
                "fisher_p_greater": fisher_p,
                "permutation_p_greater_equal_tp": perm_p,
            }
        )

        for group in pd.unique(groups):
            idx = np.where((groups == group) & finite)[0]
            if len(idx) == 0:
                continue
            g_obs = obs_hot[idx]
            g_pred = pred_hot[idx]
            g_tp = int(np.sum(g_obs & g_pred))
            g_fp = int(np.sum(~g_obs & g_pred))
            g_fn = int(np.sum(g_obs & ~g_pred))
            g_tn = int(np.sum(~g_obs & ~g_pred))
            g_odds, g_p = safe_fisher(g_tp, g_fp, g_fn, g_tn)
            sample_rows.append(
                {
                    "evidence": evidence,
                    "quantile": q,
                    "group": group,
                    "n": int(len(idx)),
                    "tp": g_tp,
                    "fp": g_fp,
                    "fn": g_fn,
                    "tn": g_tn,
                    "precision": g_tp / max(g_tp + g_fp, 1),
                    "recall": g_tp / max(g_tp + g_fn, 1),
                    "jaccard": g_tp / max(g_tp + g_fp + g_fn, 1),
                    "fisher_odds_ratio": g_odds,
                    "fisher_p_greater": g_p,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(sample_rows)


def load_evidence() -> list[tuple[str, pd.DataFrame, str, str, str]]:
    g250 = pd.read_csv(GSE250 / "path2space_spot_predictions.tsv.gz", sep="\t")
    g230 = pd.read_csv(GSE230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    g230_resid = pd.read_csv(GSE230_RESID / "gse230424_residual_target_predictions.tsv.gz", sep="\t")
    return [
        (
            "GSE250521 UNI DM1/RAI smoothed",
            g250,
            "sample_id",
            "obs_DM1_like_score_smooth8",
            "pred_DM1_like_score_smooth8",
        ),
        (
            "GSE230424 H&E DM1/low-RAI smoothed",
            g230,
            "sample",
            "obs_DM1_low_RAI_score_smooth8",
            "pred_DM1_low_RAI_score_HE_tile_features",
        ),
        (
            "GSE230424 H&E DM1/low-RAI coord+QC residual",
            g230_resid,
            "sample",
            "obs_resid_DM1_low_RAI_score_coord_qc",
            "pred_resid_DM1_low_RAI_score_HE_tile_features",
        ),
    ]


def write_figure(summary: pd.DataFrame, by_group: pd.DataFrame) -> None:
    q90 = summary[summary["quantile"] == 0.90].copy()
    q80 = summary[summary["quantile"] == 0.80].copy()
    order = q90.sort_values("precision_lift", ascending=True)["evidence"].tolist()

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5))

    for ax, metric, title, xlabel in [
        (axes[0, 0], "precision_lift", "A. Top-decile hotspot precision lift", "Precision / expected rate"),
        (axes[0, 1], "tp_over_expected", "B. Top-decile overlap over independence", "Observed TP / expected TP"),
    ]:
        vals = q90.set_index("evidence").loc[order][metric]
        ax.barh(np.arange(len(order)), vals, color="#266b73")
        ax.axvline(1, color="#333", lw=0.8)
        ax.set_yticks(np.arange(len(order)))
        ax.set_yticklabels([x.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for x in order])
        ax.set_xlabel(xlabel)
        ax.set_title(title)

    wide = summary.pivot(index="evidence", columns="quantile", values="precision")
    wide = wide.loc[order]
    x = np.arange(len(wide))
    axes[1, 0].bar(x - 0.18, wide[0.80], width=0.35, color="#8aa6a3", label="top 20%")
    axes[1, 0].bar(x + 0.18, wide[0.90], width=0.35, color="#d5a84d", label="top 10%")
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels([idx.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for idx in wide.index], rotation=25, ha="right")
    axes[1, 0].set_ylabel("Precision")
    axes[1, 0].set_title("C. Hotspot precision by threshold")
    axes[1, 0].legend(frameon=False)

    sub = by_group[by_group["quantile"] == 0.90].copy()
    evidence_to_x = {e: i for i, e in enumerate(order)}
    for evidence, part in sub.groupby("evidence", sort=False):
        xs = np.full(len(part), evidence_to_x[evidence], dtype=float)
        jitter = np.linspace(-0.13, 0.13, len(part)) if len(part) > 1 else np.zeros(len(part))
        axes[1, 1].scatter(xs + jitter, part["precision"], s=42, alpha=0.8, label=evidence)
    axes[1, 1].set_xticks(np.arange(len(order)))
    axes[1, 1].set_xticklabels([idx.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for idx in order], rotation=25, ha="right")
    axes[1, 1].set_ylabel("Per-slide/sample top-decile precision")
    axes[1, 1].set_title("D. Within-slide/sample consistency")

    fig.suptitle("Path2Space-inspired hotspot concordance: predicted high-DM1/RAI regions recover observed high-DM1/RAI regions", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_hotspot_concordance.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_hotspot_concordance.pdf")
    plt.close(fig)


def write_report(summary: pd.DataFrame, by_group: pd.DataFrame) -> dict:
    q90 = summary[summary["quantile"] == 0.90].set_index("evidence")
    result = {
        "n_permutations": N_PERM,
        "gse250521_top10_precision": float(q90.loc["GSE250521 UNI DM1/RAI smoothed", "precision"]),
        "gse250521_top10_lift": float(q90.loc["GSE250521 UNI DM1/RAI smoothed", "precision_lift"]),
        "gse250521_top10_or": float(q90.loc["GSE250521 UNI DM1/RAI smoothed", "fisher_odds_ratio"]),
        "gse250521_top10_perm_p": float(q90.loc["GSE250521 UNI DM1/RAI smoothed", "permutation_p_greater_equal_tp"]),
        "gse230424_top10_precision": float(q90.loc["GSE230424 H&E DM1/low-RAI smoothed", "precision"]),
        "gse230424_top10_lift": float(q90.loc["GSE230424 H&E DM1/low-RAI smoothed", "precision_lift"]),
        "gse230424_top10_or": float(q90.loc["GSE230424 H&E DM1/low-RAI smoothed", "fisher_odds_ratio"]),
        "gse230424_top10_perm_p": float(q90.loc["GSE230424 H&E DM1/low-RAI smoothed", "permutation_p_greater_equal_tp"]),
        "gse230424_resid_top10_precision": float(q90.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "precision"]),
        "gse230424_resid_top10_lift": float(q90.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "precision_lift"]),
        "gse230424_resid_top10_or": float(q90.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "fisher_odds_ratio"]),
        "gse230424_resid_top10_perm_p": float(q90.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "permutation_p_greater_equal_tp"]),
    }
    lines = [
        "# Path2Space-inspired hotspot concordance",
        "",
        "## Verdict",
        "",
        f"- GSE250521 top-decile hotspot precision: **{result['gse250521_top10_precision']:.3f}**; lift **{result['gse250521_top10_lift']:.2f}x**; Fisher OR **{result['gse250521_top10_or']:.2f}**; permutation p = **{result['gse250521_top10_perm_p']:.4f}**.",
        f"- GSE230424 top-decile hotspot precision: **{result['gse230424_top10_precision']:.3f}**; lift **{result['gse230424_top10_lift']:.2f}x**; Fisher OR **{result['gse230424_top10_or']:.2f}**; permutation p = **{result['gse230424_top10_perm_p']:.4f}**.",
        f"- GSE230424 coord+QC residual top-decile precision: **{result['gse230424_resid_top10_precision']:.3f}**; lift **{result['gse230424_resid_top10_lift']:.2f}x**; Fisher OR **{result['gse230424_resid_top10_or']:.2f}**; permutation p = **{result['gse230424_resid_top10_perm_p']:.4f}**.",
        "",
        "## Interpretation",
        "",
        "The hotspot test asks whether predicted high-DM1/RAI regions recover observed high-DM1/RAI regions within each slide/sample, not merely whether continuous scores are correlated. Top-decile enrichment supports spatial localization. The residual-target version is expectedly weaker but is the strictest localization test after coordinate+QC removal.",
        "",
        "## Summary",
        "",
        summary.round(4).to_markdown(index=False),
        "",
        "## Per-Slide/Sample Detail",
        "",
        by_group.round(4).to_markdown(index=False),
    ]
    (OUT / "PATH2SPACE_HOTSPOT_CONCORDANCE_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "PATH2SPACE_HOTSPOT_CONCORDANCE_SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    summaries = []
    group_rows = []
    for evidence, df, group_col, obs_col, pred_col in load_evidence():
        summary, by_group = evaluate(df, evidence, group_col, obs_col, pred_col, rng)
        summaries.append(summary)
        group_rows.append(by_group)
    summary_df = pd.concat(summaries, ignore_index=True)
    group_df = pd.concat(group_rows, ignore_index=True)
    summary_df.to_csv(OUT / "path2space_hotspot_concordance_summary.tsv", sep="\t", index=False, na_rep="NA")
    group_df.to_csv(OUT / "path2space_hotspot_concordance_by_group.tsv", sep="\t", index=False, na_rep="NA")
    write_figure(summary_df, group_df)
    result = write_report(summary_df, group_df)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
