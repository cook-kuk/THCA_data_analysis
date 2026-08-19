#!/usr/bin/env python3
"""Nearest-hotspot distance controls for Path2Space-inspired evidence.

This asks whether predicted high-DM1/RAI spots are spatially closer to observed
high-DM1/RAI spots than expected under sample- or spatial-block-preserving nulls.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2 = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
GSE250 = P2 / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"
GSE230 = P2 / "analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"
GSE230_RESID = GSE230 / "residual_target_controls"
OUT = P2 / "analysis_supp/path2space_hotspot_distance_controls_2026_05_09"
TOP_Q = 0.90
N_PERM = 1000
RANDOM_SEED = 250230427


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
        out.loc[idx[keep], "spatial_block"] = np.array([f"{group}:B{lab + 1:02d}" for lab in labs], dtype=object)
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


def median_grid_spacing(coords: np.ndarray) -> float:
    coords = np.asarray(coords, dtype=float)
    keep = np.isfinite(coords).all(axis=1)
    coords = coords[keep]
    if len(coords) < 3:
        return 1.0
    dist, _ = cKDTree(coords).query(coords, k=2)
    spacing = float(np.nanmedian(dist[:, 1]))
    return spacing if np.isfinite(spacing) and spacing > 0 else 1.0


def nearest_distances(
    df: pd.DataFrame,
    group_col: str,
    row_col: str,
    col_col: str,
    pred_hot: np.ndarray,
    obs_hot: np.ndarray,
    finite: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    coords = df[[row_col, col_col]].to_numpy(dtype=float)
    groups = df[group_col].astype(str).to_numpy()
    all_dists = []
    rows = []
    for group in pd.unique(groups):
        idx = np.where((groups == group) & finite)[0]
        if len(idx) < 5:
            continue
        spacing = median_grid_spacing(coords[idx])
        pred_idx = idx[pred_hot[idx]]
        obs_idx = idx[obs_hot[idx]]
        if len(pred_idx) == 0 or len(obs_idx) == 0:
            continue
        d, _ = cKDTree(coords[obs_idx]).query(coords[pred_idx], k=1)
        d_norm = d / spacing
        all_dists.append(d_norm)
        rows.append(
            {
                "group": group,
                "n": int(len(idx)),
                "n_pred_hot": int(len(pred_idx)),
                "n_obs_hot": int(len(obs_idx)),
                "grid_spacing": spacing,
                "median_nearest_norm_dist": float(np.nanmedian(d_norm)),
                "mean_nearest_norm_dist": float(np.nanmean(d_norm)),
                "frac_within_1_grid": float(np.mean(d_norm <= 1.0)),
                "frac_within_2_grid": float(np.mean(d_norm <= 2.0)),
            }
        )
    if all_dists:
        return np.concatenate(all_dists), pd.DataFrame(rows)
    return np.array([], dtype=float), pd.DataFrame(rows)


def permute_mask(mask: np.ndarray, finite: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = mask.copy()
    strata = strata.astype(str)
    for stratum in pd.unique(strata[finite]):
        idx = np.where(finite & (strata == stratum))[0]
        if len(idx) > 1:
            out[idx] = rng.permutation(out[idx])
    return out


def distance_metric(dist: np.ndarray) -> tuple[float, float]:
    if len(dist) == 0:
        return np.nan, np.nan
    return float(np.nanmedian(dist)), float(np.mean(dist <= 2.0))


def permutation_summary(
    df: pd.DataFrame,
    group_col: str,
    row_col: str,
    col_col: str,
    pred_hot: np.ndarray,
    obs_hot: np.ndarray,
    finite: np.ndarray,
    strata: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    median_null = np.full(N_PERM, np.nan)
    frac2_null = np.full(N_PERM, np.nan)
    for i in range(N_PERM):
        perm_obs = permute_mask(obs_hot, finite, strata, rng)
        dist, _ = nearest_distances(df, group_col, row_col, col_col, pred_hot, perm_obs, finite)
        median_null[i], frac2_null[i] = distance_metric(dist)
    return median_null, frac2_null


def evaluate(
    evidence: str,
    df: pd.DataFrame,
    group_col: str,
    row_col: str,
    col_col: str,
    obs_col: str,
    pred_col: str,
    rng: np.random.Generator,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    work = assign_domains(df, group_col, row_col, col_col)
    obs = work[obs_col].to_numpy(dtype=float)
    pred = work[pred_col].to_numpy(dtype=float)
    groups = work[group_col].astype(str).to_numpy()
    blocks = work["spatial_block"].astype(str).to_numpy()
    finite = np.isfinite(obs) & np.isfinite(pred) & np.isfinite(work[[row_col, col_col]].to_numpy(dtype=float)).all(axis=1)

    obs_hot = top_k_mask(obs, groups, TOP_Q)
    pred_hot = top_k_mask(pred, groups, TOP_Q)
    dist, group_detail = nearest_distances(work, group_col, row_col, col_col, pred_hot, obs_hot, finite)
    obs_median, obs_frac2 = distance_metric(dist)

    sample_null_median, sample_null_frac2 = permutation_summary(
        work, group_col, row_col, col_col, pred_hot, obs_hot, finite, groups, rng
    )
    block_null_median, block_null_frac2 = permutation_summary(
        work, group_col, row_col, col_col, pred_hot, obs_hot, finite, blocks, rng
    )
    p_sample = float((1 + np.sum(sample_null_median <= obs_median)) / (1 + np.isfinite(sample_null_median).sum()))
    p_block = float((1 + np.sum(block_null_median <= obs_median)) / (1 + np.isfinite(block_null_median).sum()))
    p_frac2_sample = float((1 + np.sum(sample_null_frac2 >= obs_frac2)) / (1 + np.isfinite(sample_null_frac2).sum()))
    p_frac2_block = float((1 + np.sum(block_null_frac2 >= obs_frac2)) / (1 + np.isfinite(block_null_frac2).sum()))

    row = {
        "evidence": evidence,
        "n": int(np.sum(finite)),
        "n_groups": int(pd.Series(groups[finite]).nunique()),
        "n_spatial_blocks": int(pd.Series(blocks[finite]).nunique()),
        "n_pred_hotspots": int(np.sum(pred_hot & finite)),
        "n_obs_hotspots": int(np.sum(obs_hot & finite)),
        "median_nearest_norm_dist": obs_median,
        "frac_pred_hotspots_within_2_grid": obs_frac2,
        "sample_null_median_dist_mean": float(np.nanmean(sample_null_median)),
        "sample_null_median_dist_p05": float(np.nanquantile(sample_null_median, 0.05)),
        "sample_null_frac2_mean": float(np.nanmean(sample_null_frac2)),
        "sample_stratified_perm_p_lower_dist": p_sample,
        "sample_stratified_perm_p_higher_frac2": p_frac2_sample,
        "block_null_median_dist_mean": float(np.nanmean(block_null_median)),
        "block_null_median_dist_p05": float(np.nanquantile(block_null_median, 0.05)),
        "block_null_frac2_mean": float(np.nanmean(block_null_frac2)),
        "block_stratified_perm_p_lower_dist": p_block,
        "block_stratified_perm_p_higher_frac2": p_frac2_block,
    }
    group_detail.insert(0, "evidence", evidence)
    null_df = pd.DataFrame(
        {
            "evidence": evidence,
            "perm_idx": np.arange(N_PERM),
            "sample_null_median_dist": sample_null_median,
            "sample_null_frac2": sample_null_frac2,
            "block_null_median_dist": block_null_median,
            "block_null_frac2": block_null_frac2,
        }
    )
    return row, group_detail, null_df


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


def write_figure(summary: pd.DataFrame, group_detail: pd.DataFrame) -> None:
    order = summary.sort_values("median_nearest_norm_dist", ascending=False)["evidence"].tolist()
    labels = [x.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for x in order]
    plot = summary.set_index("evidence").loc[order]

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.2))
    y = np.arange(len(order))

    axes[0, 0].barh(y - 0.18, plot["median_nearest_norm_dist"], height=0.35, color="#266b73", label="observed")
    axes[0, 0].barh(y + 0.18, plot["block_null_median_dist_mean"], height=0.35, color="#9aa4ad", label="block null mean")
    axes[0, 0].set_yticks(y)
    axes[0, 0].set_yticklabels(labels)
    axes[0, 0].set_xlabel("Median nearest observed-hotspot distance (grid-normalized)")
    axes[0, 0].set_title("A. Predicted hotspots are closer to observed hotspots")
    axes[0, 0].legend(frameon=False)

    axes[0, 1].barh(y - 0.18, plot["frac_pred_hotspots_within_2_grid"], height=0.35, color="#266b73", label="observed")
    axes[0, 1].barh(y + 0.18, plot["block_null_frac2_mean"], height=0.35, color="#9aa4ad", label="block null mean")
    axes[0, 1].set_yticks(y)
    axes[0, 1].set_yticklabels(labels)
    axes[0, 1].set_xlabel("Fraction within 2 grid spacings")
    axes[0, 1].set_title("B. Nearby-hotspot recovery")
    axes[0, 1].legend(frameon=False)

    axes[1, 0].barh(y, plot["block_stratified_perm_p_lower_dist"], color="#d5a84d")
    axes[1, 0].axvline(0.05, color="#b44d4d", lw=1.2)
    axes[1, 0].set_yticks(y)
    axes[1, 0].set_yticklabels(labels)
    axes[1, 0].set_xlabel("Block-stratified permutation p (lower distance)")
    axes[1, 0].set_title("C. Distance enrichment p-value")

    for evidence, part in group_detail.groupby("evidence", sort=False):
        x0 = order.index(evidence)
        jitter = np.linspace(-0.13, 0.13, len(part)) if len(part) > 1 else np.zeros(len(part))
        axes[1, 1].scatter(np.full(len(part), x0) + jitter, part["median_nearest_norm_dist"], s=42, alpha=0.8)
    axes[1, 1].set_xticks(np.arange(len(order)))
    axes[1, 1].set_xticklabels(labels, rotation=25, ha="right")
    axes[1, 1].set_ylabel("Per-slide/sample median nearest distance")
    axes[1, 1].set_title("D. Per-slide/sample proximity")

    fig.suptitle("Nearest-hotspot distance controls: strong in GSE250521, caveated in dense GSE230424", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_hotspot_distance_controls.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_hotspot_distance_controls.pdf")
    plt.close(fig)


def write_report(summary: pd.DataFrame, group_detail: pd.DataFrame) -> dict:
    idx = summary.set_index("evidence")
    result = {
        "n_permutations": N_PERM,
        "gse250521_median_distance": float(idx.loc["GSE250521 UNI DM1/RAI smoothed", "median_nearest_norm_dist"]),
        "gse250521_block_null_median_distance": float(idx.loc["GSE250521 UNI DM1/RAI smoothed", "block_null_median_dist_mean"]),
        "gse250521_block_perm_p": float(idx.loc["GSE250521 UNI DM1/RAI smoothed", "block_stratified_perm_p_lower_dist"]),
        "gse230424_median_distance": float(idx.loc["GSE230424 H&E DM1/low-RAI smoothed", "median_nearest_norm_dist"]),
        "gse230424_block_null_median_distance": float(idx.loc["GSE230424 H&E DM1/low-RAI smoothed", "block_null_median_dist_mean"]),
        "gse230424_block_perm_p": float(idx.loc["GSE230424 H&E DM1/low-RAI smoothed", "block_stratified_perm_p_lower_dist"]),
        "gse230424_resid_median_distance": float(idx.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "median_nearest_norm_dist"]),
        "gse230424_resid_block_null_median_distance": float(idx.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "block_null_median_dist_mean"]),
        "gse230424_resid_block_perm_p": float(idx.loc["GSE230424 H&E DM1/low-RAI coord+QC residual", "block_stratified_perm_p_lower_dist"]),
    }
    lines = [
        "# Path2Space-inspired nearest-hotspot distance controls",
        "",
        "## Verdict",
        "",
        f"- GSE250521 predicted top-decile hotspots median distance to observed hotspots: **{result['gse250521_median_distance']:.3f} grid spacings** vs block-null mean **{result['gse250521_block_null_median_distance']:.3f}**; p = **{result['gse250521_block_perm_p']:.4f}**.",
        f"- GSE230424 predicted top-decile hotspots median distance: **{result['gse230424_median_distance']:.3f}** vs block-null mean **{result['gse230424_block_null_median_distance']:.3f}**; p = **{result['gse230424_block_perm_p']:.4f}**.",
        f"- GSE230424 coord+QC residual predicted top-decile hotspots median distance: **{result['gse230424_resid_median_distance']:.3f}** vs block-null mean **{result['gse230424_resid_block_null_median_distance']:.3f}**; p = **{result['gse230424_resid_block_perm_p']:.4f}**.",
        "",
        "## Interpretation",
        "",
        "Overlap asks whether the same spots are recovered. The nearest-hotspot distance test asks whether predicted hotspots are spatially proximal to observed hotspots, allowing small localization error. Block-stratified nulls preserve broad regional hotspot density before testing whether predicted hotspots are unusually close to observed hotspots. GSE250521 passes this proximity test. In dense GSE230424 arrays, the exact hotspot-overlap and block-centered correlation layers remain positive, but nearest-distance is dominated by spatial block density and should be treated as a caveat rather than support.",
        "",
        "## Summary",
        "",
        summary.round(4).to_markdown(index=False),
        "",
        "## Per-Slide/Sample Detail",
        "",
        group_detail.round(4).to_markdown(index=False),
    ]
    (OUT / "PATH2SPACE_HOTSPOT_DISTANCE_CONTROLS_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "PATH2SPACE_HOTSPOT_DISTANCE_CONTROLS_SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    group_details = []
    null_dfs = []
    for evidence in load_evidence():
        row, group_detail, null_df = evaluate(*evidence, rng=rng)
        rows.append(row)
        group_details.append(group_detail)
        null_dfs.append(null_df)
    summary = pd.DataFrame(rows)
    group_detail = pd.concat(group_details, ignore_index=True)
    null_df = pd.concat(null_dfs, ignore_index=True)
    summary.to_csv(OUT / "path2space_hotspot_distance_summary.tsv", sep="\t", index=False, na_rep="NA")
    group_detail.to_csv(OUT / "path2space_hotspot_distance_by_group.tsv", sep="\t", index=False, na_rep="NA")
    null_df.to_csv(OUT / "path2space_hotspot_distance_null.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    write_figure(summary, group_detail)
    result = write_report(summary, group_detail)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
