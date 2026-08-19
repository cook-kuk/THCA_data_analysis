#!/usr/bin/env python3
"""Spatial-block controls for Path2Space-inspired hotspot evidence.

These controls ask whether image-to-DM1/RAI concordance remains after removing
broad coordinate-domain structure within each slide/sample.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
GSE250 = P2 / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"
GSE230 = P2 / "analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"
GSE230_RESID = GSE230 / "residual_target_controls"
OUT = P2 / "analysis_supp/path2space_spatial_block_controls_2026_05_09"
N_PERM = 1000
TOP_Q = 0.90
RANDOM_SEED = 250230426


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def center_by_group(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    out = np.full(len(values), np.nan)
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        if keep.sum() == 0:
            continue
        out[idx[keep]] = vals[keep] - np.nanmean(vals[keep])
    return out


def assign_domains(df: pd.DataFrame, group_col: str, row_col: str, col_col: str) -> pd.DataFrame:
    out = df.copy()
    out["spatial_block"] = "NA"
    for group, idx_raw in out.groupby(group_col).groups.items():
        idx = np.asarray(list(idx_raw))
        coords = out.loc[idx, [row_col, col_col]].to_numpy(dtype=float)
        keep = np.isfinite(coords).all(axis=1)
        if keep.sum() < 25:
            out.loc[idx[keep], "spatial_block"] = [f"{group}:B1"] * int(keep.sum())
            continue
        n_domains = int(np.clip(round(keep.sum() / 350), 4, 12))
        n_domains = min(n_domains, keep.sum())
        labs = KMeans(n_clusters=n_domains, n_init=50, random_state=RANDOM_SEED).fit_predict(
            StandardScaler().fit_transform(coords[keep])
        )
        block_names = np.array([f"{group}:B{lab + 1:02d}" for lab in labs], dtype=object)
        out.loc[idx[keep], "spatial_block"] = block_names
    return out


def top_k_mask(values: np.ndarray, groups: np.ndarray, q: float = TOP_Q) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    mask = np.zeros(len(values), dtype=bool)
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        if keep.sum() < 5:
            continue
        valid_idx = idx[keep]
        n_top = max(1, int(np.ceil((1.0 - q) * len(valid_idx))))
        order = np.argsort(vals[keep])[::-1]
        mask[valid_idx[order[:n_top]]] = True
    return mask


def hotspot_stats(obs_hot: np.ndarray, pred_hot: np.ndarray, finite: np.ndarray) -> dict:
    tp = int(np.sum(obs_hot & pred_hot & finite))
    fp = int(np.sum(~obs_hot & pred_hot & finite))
    fn = int(np.sum(obs_hot & ~pred_hot & finite))
    tn = int(np.sum(~obs_hot & ~pred_hot & finite))
    odds, fisher_p = stats.fisher_exact([[tp, fp], [fn, tn]], alternative="greater")
    n = int(np.sum(finite))
    obs_rate = (tp + fn) / max(n, 1)
    pred_rate = (tp + fp) / max(n, 1)
    precision = tp / max(tp + fp, 1)
    expected = obs_rate * pred_rate * n
    return {
        "n": n,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "obs_hotspot_rate": obs_rate,
        "pred_hotspot_rate": pred_rate,
        "precision": precision,
        "recall": tp / max(tp + fn, 1),
        "precision_lift": precision / obs_rate if obs_rate > 0 else np.nan,
        "jaccard": tp / max(tp + fp + fn, 1),
        "expected_tp_independent": expected,
        "tp_over_expected": tp / expected if expected > 0 else np.nan,
        "fisher_odds_ratio": float(odds),
        "fisher_p_greater": float(fisher_p),
    }


def domain_stratified_perm_p(
    obs_hot: np.ndarray,
    pred_hot: np.ndarray,
    finite: np.ndarray,
    blocks: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    observed_tp = int(np.sum(obs_hot & pred_hot & finite))
    null = np.zeros(N_PERM, dtype=int)
    finite_idx = np.where(finite)[0]
    block_vals = blocks.astype(str)
    for b in range(N_PERM):
        perm = pred_hot.copy()
        for block in pd.unique(block_vals[finite_idx]):
            idx = np.where(finite & (block_vals == block))[0]
            if len(idx) > 1:
                perm[idx] = rng.permutation(perm[idx])
        null[b] = int(np.sum(obs_hot & perm & finite))
    p = float((1 + np.sum(null >= observed_tp)) / (1 + N_PERM))
    return p, float(np.mean(null)), float(np.quantile(null, 0.95))


def evaluate(
    evidence: str,
    df: pd.DataFrame,
    group_col: str,
    row_col: str,
    col_col: str,
    obs_col: str,
    pred_col: str,
    rng: np.random.Generator,
) -> tuple[dict, pd.DataFrame]:
    work = assign_domains(df, group_col, row_col, col_col)
    obs = work[obs_col].to_numpy(dtype=float)
    pred = work[pred_col].to_numpy(dtype=float)
    groups = work[group_col].astype(str).to_numpy()
    blocks = work["spatial_block"].astype(str).to_numpy()
    finite = np.isfinite(obs) & np.isfinite(pred) & (blocks != "NA")

    sample_rho, sample_p, n = safe_spearman(center_by_group(obs, groups), center_by_group(pred, groups))
    block_rho, block_p, _ = safe_spearman(center_by_group(obs, blocks), center_by_group(pred, blocks))

    obs_hot = top_k_mask(obs, groups, TOP_Q)
    pred_hot = top_k_mask(pred, groups, TOP_Q)
    hot = hotspot_stats(obs_hot, pred_hot, finite)
    domain_p, null_mean, null_p95 = domain_stratified_perm_p(obs_hot, pred_hot, finite, blocks, rng)

    per_group_rows = []
    for group in pd.unique(groups):
        idx = finite & (groups == group)
        if idx.sum() < 10:
            continue
        g_sample_rho, _, _ = safe_spearman(obs[idx], pred[idx])
        g_block_rho, _, _ = safe_spearman(center_by_group(obs[idx], blocks[idx]), center_by_group(pred[idx], blocks[idx]))
        g_hot = hotspot_stats(obs_hot, pred_hot, idx)
        per_group_rows.append(
            {
                "evidence": evidence,
                "group": group,
                "n": int(idx.sum()),
                "n_blocks": int(pd.Series(blocks[idx]).nunique()),
                "sample_rho": g_sample_rho,
                "block_centered_rho": g_block_rho,
                "top10_precision": g_hot["precision"],
                "top10_precision_lift": g_hot["precision_lift"],
                "top10_or": g_hot["fisher_odds_ratio"],
            }
        )

    row = {
        "evidence": evidence,
        "n": n,
        "n_groups": int(pd.Series(groups[finite]).nunique()),
        "n_spatial_blocks": int(pd.Series(blocks[finite]).nunique()),
        "sample_centered_rho": sample_rho,
        "sample_centered_p": sample_p,
        "spatial_block_centered_rho": block_rho,
        "spatial_block_centered_p": block_p,
        "top10_precision": hot["precision"],
        "top10_precision_lift": hot["precision_lift"],
        "top10_tp_over_expected": hot["tp_over_expected"],
        "top10_or": hot["fisher_odds_ratio"],
        "top10_fisher_p": hot["fisher_p_greater"],
        "top10_domain_stratified_perm_p": domain_p,
        "domain_stratified_null_tp_mean": null_mean,
        "domain_stratified_null_tp_p95": null_p95,
        "observed_tp": hot["tp"],
    }
    return row, pd.DataFrame(per_group_rows)


def load_evidence() -> list[tuple[str, pd.DataFrame, str, str, str, str, str]]:
    g250 = pd.read_csv(GSE250 / "path2space_spot_predictions.tsv.gz", sep="\t")
    g230 = pd.read_csv(GSE230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    g230_resid = pd.read_csv(GSE230_RESID / "gse230424_residual_target_predictions.tsv.gz", sep="\t")
    return [
        (
            "GSE250521 UNI DM1/RAI smoothed",
            g250,
            "sample_id",
            "array_row",
            "array_col",
            "obs_DM1_like_score_smooth8",
            "pred_DM1_like_score_smooth8",
        ),
        (
            "GSE230424 H&E DM1/low-RAI smoothed",
            g230,
            "sample",
            "array_row",
            "array_col",
            "obs_DM1_low_RAI_score_smooth8",
            "pred_DM1_low_RAI_score_HE_tile_features",
        ),
        (
            "GSE230424 H&E DM1/low-RAI coord+QC residual",
            g230_resid,
            "sample",
            "array_row",
            "array_col",
            "obs_resid_DM1_low_RAI_score_coord_qc",
            "pred_resid_DM1_low_RAI_score_HE_tile_features",
        ),
    ]


def write_figure(summary: pd.DataFrame, by_group: pd.DataFrame) -> None:
    order = summary.sort_values("spatial_block_centered_rho", ascending=True)["evidence"].tolist()
    labels = [x.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for x in order]

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.2))
    y = np.arange(len(order))
    plot = summary.set_index("evidence").loc[order]

    axes[0, 0].barh(y - 0.18, plot["sample_centered_rho"], height=0.35, color="#8aa6a3", label="sample-centered")
    axes[0, 0].barh(y + 0.18, plot["spatial_block_centered_rho"], height=0.35, color="#d5a84d", label="block-centered")
    axes[0, 0].axvline(0, color="#333", lw=0.8)
    axes[0, 0].set_yticks(y)
    axes[0, 0].set_yticklabels(labels)
    axes[0, 0].set_xlabel("Spearman rho")
    axes[0, 0].set_title("A. Broad spatial-block removal")
    axes[0, 0].legend(frameon=False)

    axes[0, 1].barh(y, plot["top10_domain_stratified_perm_p"], color="#266b73")
    axes[0, 1].axvline(0.05, color="#b44d4d", lw=1.2)
    axes[0, 1].set_yticks(y)
    axes[0, 1].set_yticklabels(labels)
    axes[0, 1].set_xlabel("Domain-stratified permutation p")
    axes[0, 1].set_title("B. Hotspot enrichment after block-preserving null")

    axes[1, 0].barh(y, plot["top10_precision_lift"], color="#266b73")
    axes[1, 0].axvline(1, color="#333", lw=0.8)
    axes[1, 0].set_yticks(y)
    axes[1, 0].set_yticklabels(labels)
    axes[1, 0].set_xlabel("Top-decile precision lift")
    axes[1, 0].set_title("C. Exact top-decile hotspot lift")

    for evidence, part in by_group.groupby("evidence", sort=False):
        x0 = order.index(evidence)
        jitter = np.linspace(-0.13, 0.13, len(part)) if len(part) > 1 else np.zeros(len(part))
        axes[1, 1].scatter(np.full(len(part), x0) + jitter, part["block_centered_rho"], s=42, alpha=0.8)
    axes[1, 1].axhline(0, color="#333", lw=0.8)
    axes[1, 1].set_xticks(np.arange(len(order)))
    axes[1, 1].set_xticklabels(labels, rotation=25, ha="right")
    axes[1, 1].set_ylabel("Per-slide/sample block-centered rho")
    axes[1, 1].set_title("D. Per-slide/sample block-centered consistency")

    fig.suptitle("Spatial-block controls: image-to-DM1/RAI signal beyond broad coordinate domains", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_spatial_block_controls.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_spatial_block_controls.pdf")
    plt.close(fig)


def write_report(summary: pd.DataFrame, by_group: pd.DataFrame) -> dict:
    idx = summary.set_index("evidence")
    result = {
        "gse250521_block_centered_rho": float(idx.loc["GSE250521 UNI DM1/RAI smoothed", "spatial_block_centered_rho"]),
        "gse250521_block_hotspot_perm_p": float(idx.loc["GSE250521 UNI DM1/RAI smoothed", "top10_domain_stratified_perm_p"]),
        "gse230424_block_centered_rho": float(idx.loc["GSE230424 H&E DM1/low-RAI smoothed", "spatial_block_centered_rho"]),
        "gse230424_block_hotspot_perm_p": float(idx.loc["GSE230424 H&E DM1/low-RAI smoothed", "top10_domain_stratified_perm_p"]),
        "gse230424_resid_block_centered_rho": float(idx.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "spatial_block_centered_rho"]),
        "gse230424_resid_block_hotspot_perm_p": float(idx.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "top10_domain_stratified_perm_p"]),
        "n_permutations": N_PERM,
    }
    lines = [
        "# Path2Space-inspired spatial-block controls",
        "",
        "## Verdict",
        "",
        f"- GSE250521 block-centered rho: **{result['gse250521_block_centered_rho']:.3f}**; top-decile block-stratified permutation p = **{result['gse250521_block_hotspot_perm_p']:.4f}**.",
        f"- GSE230424 block-centered rho: **{result['gse230424_block_centered_rho']:.3f}**; top-decile block-stratified permutation p = **{result['gse230424_block_hotspot_perm_p']:.4f}**.",
        f"- GSE230424 coord+QC residual block-centered rho: **{result['gse230424_resid_block_centered_rho']:.3f}**; top-decile block-stratified permutation p = **{result['gse230424_resid_block_hotspot_perm_p']:.4f}**.",
        "",
        "## Interpretation",
        "",
        "The block-centered test removes broad coordinate-domain means before computing correlation. The block-stratified hotspot null shuffles predicted hotspot labels only within spatial blocks, preserving broad predicted hotspot density by region. Residual signal under these controls argues against a purely broad-domain artifact.",
        "",
        "## Summary",
        "",
        summary.round(4).to_markdown(index=False),
        "",
        "## Per-Slide/Sample Detail",
        "",
        by_group.round(4).to_markdown(index=False),
    ]
    (OUT / "PATH2SPACE_SPATIAL_BLOCK_CONTROLS_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "PATH2SPACE_SPATIAL_BLOCK_CONTROLS_SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    group_dfs = []
    for evidence in load_evidence():
        row, group_df = evaluate(*evidence, rng=rng)
        rows.append(row)
        group_dfs.append(group_df)
    summary = pd.DataFrame(rows)
    by_group = pd.concat(group_dfs, ignore_index=True)
    summary.to_csv(OUT / "path2space_spatial_block_controls_summary.tsv", sep="\t", index=False, na_rep="NA")
    by_group.to_csv(OUT / "path2space_spatial_block_controls_by_group.tsv", sep="\t", index=False, na_rep="NA")
    write_figure(summary, by_group)
    result = write_report(summary, by_group)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
