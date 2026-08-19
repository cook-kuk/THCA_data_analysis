#!/usr/bin/env python3
"""Spatial-autocorrelation-adjusted specificity for GSE230424 random modules.

The random-module control showed that DM1/RAI sits in the upper tail of
matched random modules. This follow-up asks whether that tail position is
explained only by the DM1/RAI target being more spatially autocorrelated.
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
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
BASE = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"
RANDOM_DIR = BASE / "random_module_specificity"
OUT = BASE / "spatial_autocorr_specificity"
K_NEIGH = 6


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def residualize_within_sample_many(y: np.ndarray, covars: np.ndarray, samples: np.ndarray) -> np.ndarray:
    out = np.full(y.shape, np.nan, dtype=float)
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        x = covars[idx]
        keep_x = np.isfinite(x).all(axis=1)
        if keep_x.sum() < 30:
            continue
        for j in range(y.shape[1]):
            vals = y[idx, j]
            keep = keep_x & np.isfinite(vals)
            if keep.sum() < 30:
                continue
            xs = StandardScaler().fit_transform(x[keep])
            fit = LinearRegression().fit(xs, vals[keep])
            out[idx[keep], j] = vals[keep] - fit.predict(xs)
    return out


def neighbor_autocorr(values: np.ndarray, coords: np.ndarray, samples: np.ndarray) -> tuple[float, float]:
    per_sample = []
    weights = []
    neighbor_vals_all = np.full(len(values), np.nan)
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        vals = values[idx]
        xy = coords[idx]
        keep = np.isfinite(vals) & np.isfinite(xy).all(axis=1)
        if keep.sum() < K_NEIGH + 3:
            continue
        valid_idx = idx[keep]
        tree = cKDTree(xy[keep])
        kk = min(K_NEIGH + 1, keep.sum())
        _, neigh = tree.query(xy[keep], k=kk)
        if neigh.ndim == 1:
            neigh = neigh[:, None]
        neigh = neigh[:, 1:] if neigh.shape[1] > 1 else neigh
        neigh_mean = np.nanmean(vals[keep][neigh], axis=1)
        neighbor_vals_all[valid_idx] = neigh_mean
        rho, _, n = safe_spearman(vals[keep], neigh_mean)
        if np.isfinite(rho):
            per_sample.append(rho)
            weights.append(n)
    pooled, _, _ = safe_spearman(values, neighbor_vals_all)
    if weights:
        weighted = float(np.average(per_sample, weights=weights))
    else:
        weighted = np.nan
    return pooled, weighted


def compute_spatial_blocks(samples: np.ndarray, coords: np.ndarray) -> np.ndarray:
    blocks = np.full(len(samples), "NA", dtype=object)
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        xy = coords[idx]
        keep = np.isfinite(xy).all(axis=1)
        if keep.sum() < 30:
            blocks[idx[keep]] = [f"{sample}:B1"] * int(keep.sum())
            continue
        n_domains = int(np.clip(round(keep.sum() / 350), 4, 12))
        from sklearn.cluster import KMeans

        labs = KMeans(n_clusters=n_domains, n_init=30, random_state=4242).fit_predict(
            StandardScaler().fit_transform(xy[keep])
        )
        valid_idx = idx[keep]
        blocks[valid_idx] = np.array([f"{sample}:B{lab + 1:02d}" for lab in labs], dtype=object)
    return blocks


def block_variance_fraction(values: np.ndarray, samples: np.ndarray, block_labels: np.ndarray) -> float:
    numer = 0.0
    denom = 0.0
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        vals = values[idx]
        blocks = block_labels[idx].astype(str)
        keep = np.isfinite(vals) & (blocks != "NA")
        if keep.sum() < 30:
            continue
        centered = vals[keep] - np.nanmean(vals[keep])
        denom += float(np.nansum(centered**2))
        means = np.zeros_like(centered)
        for lab in pd.unique(blocks[keep]):
            lab_mask = blocks[keep] == lab
            means[lab_mask] = np.nanmean(vals[keep][lab_mask])
        fitted = means - np.nanmean(vals[keep])
        numer += float(np.nansum(fitted**2))
    return numer / denom if denom > 0 else np.nan


def empirical_upper_p(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    return float((1.0 + np.sum(null >= observed)) / (1.0 + len(null)))


def percentile(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    return float(100.0 * np.mean(null <= observed))


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.read_csv(RANDOM_DIR / "gse230424_random_module_spot_scores.tsv.gz", sep="\t")
    results = pd.read_csv(RANDOM_DIR / "gse230424_random_module_specificity_results.tsv", sep="\t")
    return scores, results


def compute_metrics(scores: pd.DataFrame) -> pd.DataFrame:
    id_cols = ["sample", "barcode", "array_row", "array_col", "total_counts", "n_genes_by_counts", "pct_counts_mt"]
    target_cols = [c for c in scores.columns if c.endswith("_smooth8")]
    target_names = [c.replace("_smooth8", "") for c in target_cols]
    samples = scores["sample"].astype(str).to_numpy()
    coords = scores[["array_row", "array_col"]].to_numpy(dtype=float)
    block_labels = compute_spatial_blocks(samples, coords)
    y_raw = scores[target_cols].to_numpy(dtype=float)
    covars = np.column_stack(
        [
            scores[["array_row", "array_col"]].to_numpy(dtype=float),
            np.log1p(scores["total_counts"].to_numpy(dtype=float)),
            np.log1p(scores["n_genes_by_counts"].to_numpy(dtype=float)),
            scores["pct_counts_mt"].to_numpy(dtype=float),
        ]
    )
    y_resid = residualize_within_sample_many(y_raw, covars, samples)

    rows = []
    for j, target in enumerate(target_names):
        for mode, y in [("raw_smoothed", y_raw[:, j]), ("coord_qc_residual_target", y_resid[:, j])]:
            pooled, weighted = neighbor_autocorr(y, coords, samples)
            rows.append(
                {
                    "target": target,
                    "mode": mode,
                    "neighbor_autocorr_pooled": pooled,
                    "neighbor_autocorr_weighted_sample": weighted,
                    "block_variance_fraction": block_variance_fraction(y, samples, block_labels),
                    "target_sd": float(np.nanstd(y)),
                }
            )
    return pd.DataFrame(rows)


def classify(name: str) -> str:
    if name == "actual_DM1_low_RAI":
        return "actual_dm1"
    if name.startswith("loo_drop_"):
        return "panel_loo"
    return "random_matched"


def fit_adjustment(df: pd.DataFrame, mode: str) -> tuple[pd.DataFrame, dict]:
    sub = df[df["mode"] == mode].copy()
    random = sub[sub["class"] == "random_matched"].copy()
    actual = sub[sub["class"] == "actual_dm1"].iloc[0]
    features = ["neighbor_autocorr_weighted_sample", "block_variance_fraction", "target_sd"]
    keep = random[features + ["sample_centered_rho"]].replace([np.inf, -np.inf], np.nan).dropna()
    x = keep[features].to_numpy(dtype=float)
    y = keep["sample_centered_rho"].to_numpy(dtype=float)
    scaler = StandardScaler()
    model = LinearRegression().fit(scaler.fit_transform(x), y)
    sub["smoothness_expected_rho"] = np.nan
    sub["smoothness_adjusted_residual"] = np.nan
    ok = sub[features].replace([np.inf, -np.inf], np.nan).notna().all(axis=1)
    sub.loc[ok, "smoothness_expected_rho"] = model.predict(scaler.transform(sub.loc[ok, features].to_numpy(dtype=float)))
    sub.loc[ok, "smoothness_adjusted_residual"] = sub.loc[ok, "sample_centered_rho"] - sub.loc[ok, "smoothness_expected_rho"]
    random_resid = sub.loc[sub["class"] == "random_matched", "smoothness_adjusted_residual"].to_numpy(dtype=float)
    actual_resid = float(sub.loc[sub["class"] == "actual_dm1", "smoothness_adjusted_residual"].iloc[0])
    actual_expected = float(sub.loc[sub["class"] == "actual_dm1", "smoothness_expected_rho"].iloc[0])
    actual_rho = float(actual["sample_centered_rho"])
    summary = {
        "mode": mode,
        "actual_rho": actual_rho,
        "actual_expected_from_smoothness": actual_expected,
        "actual_smoothness_adjusted_residual": actual_resid,
        "residual_percentile_vs_random": percentile(actual_resid, random_resid),
        "residual_empirical_upper_p": empirical_upper_p(actual_resid, random_resid),
        "random_residual_p95": float(np.nanquantile(random_resid, 0.95)),
        "random_residual_mean": float(np.nanmean(random_resid)),
        "model_coef_neighbor_autocorr": float(model.coef_[0]),
        "model_coef_block_variance": float(model.coef_[1]),
        "model_coef_target_sd": float(model.coef_[2]),
        "model_r2_random": float(model.score(scaler.transform(x), y)),
    }
    return sub, summary


def write_figure(adjusted: pd.DataFrame, summary: dict) -> None:
    raw = adjusted[adjusted["mode"] == "raw_smoothed"].copy()
    resid = adjusted[adjusted["mode"] == "coord_qc_residual_target"].copy()

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5))
    for ax, sub, title in [
        (axes[0, 0], raw, "A. Raw target predictability vs spatial autocorrelation"),
        (axes[0, 1], resid, "B. Residual target predictability vs spatial autocorrelation"),
    ]:
        random = sub[sub["class"] == "random_matched"]
        actual = sub[sub["class"] == "actual_dm1"]
        ax.scatter(random["neighbor_autocorr_weighted_sample"], random["sample_centered_rho"], s=18, alpha=0.55, color="#888")
        ax.scatter(actual["neighbor_autocorr_weighted_sample"], actual["sample_centered_rho"], s=95, color="#b44d4d", label="DM1/RAI", zorder=5)
        ax.set_xlabel("Neighbor autocorrelation")
        ax.set_ylabel("H&E sample-centered rho")
        ax.set_title(title)
        ax.legend(frameon=False)

    for ax, sub, title in [
        (axes[1, 0], raw, "C. Raw smoothness-adjusted residual tail"),
        (axes[1, 1], resid, "D. Residual-target smoothness-adjusted residual tail"),
    ]:
        random = sub[sub["class"] == "random_matched"]
        actual_val = float(sub.loc[sub["class"] == "actual_dm1", "smoothness_adjusted_residual"].iloc[0])
        ax.hist(random["smoothness_adjusted_residual"], bins=32, color="#9aa4ad", edgecolor="white")
        ax.axvline(actual_val, color="#b44d4d", lw=2.5, label="DM1/RAI")
        ax.axvline(0, color="#333", lw=0.8)
        ax.set_xlabel("Observed rho - expected rho from smoothness")
        ax.set_ylabel("Matched random modules")
        ax.set_title(title)
        ax.legend(frameon=False)

    fig.suptitle("GSE230424 spatial-autocorrelation-adjusted random-module specificity", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse230424_spatial_autocorr_specificity.png", dpi=220)
    fig.savefig(OUT / "fig_gse230424_spatial_autocorr_specificity.pdf")
    plt.close(fig)


def write_report(adjusted: pd.DataFrame, summaries: list[dict]) -> dict:
    raw = next(s for s in summaries if s["mode"] == "raw_smoothed")
    resid = next(s for s in summaries if s["mode"] == "coord_qc_residual_target")
    summary = {
        "raw_actual_rho": raw["actual_rho"],
        "raw_expected_from_smoothness": raw["actual_expected_from_smoothness"],
        "raw_smoothness_adjusted_residual": raw["actual_smoothness_adjusted_residual"],
        "raw_residual_percentile": raw["residual_percentile_vs_random"],
        "raw_residual_empirical_p": raw["residual_empirical_upper_p"],
        "resid_actual_rho": resid["actual_rho"],
        "resid_expected_from_smoothness": resid["actual_expected_from_smoothness"],
        "resid_smoothness_adjusted_residual": resid["actual_smoothness_adjusted_residual"],
        "resid_residual_percentile": resid["residual_percentile_vs_random"],
        "resid_residual_empirical_p": resid["residual_empirical_upper_p"],
    }
    random_top = adjusted[adjusted["class"] == "random_matched"].sort_values(
        ["mode", "smoothness_adjusted_residual"], ascending=[True, False]
    ).groupby("mode").head(10)
    lines = [
        "# GSE230424 spatial-autocorrelation-adjusted specificity",
        "",
        "## Verdict",
        "",
        f"- Raw DM1/RAI H&E rho: **{raw['actual_rho']:.3f}**; expected from random-module smoothness model **{raw['actual_expected_from_smoothness']:.3f}**; adjusted residual **{raw['actual_smoothness_adjusted_residual']:.3f}**; percentile **{raw['residual_percentile_vs_random']:.1f}%**, empirical p = **{raw['residual_empirical_upper_p']:.4f}**.",
        f"- Coord+QC residual-target DM1/RAI H&E rho: **{resid['actual_rho']:.3f}**; expected from random-module smoothness model **{resid['actual_expected_from_smoothness']:.3f}**; adjusted residual **{resid['actual_smoothness_adjusted_residual']:.3f}**; percentile **{resid['residual_percentile_vs_random']:.1f}%**, empirical p = **{resid['residual_empirical_upper_p']:.4f}**.",
        "",
        "## Interpretation",
        "",
        "This control regresses random-module H&E predictability on spatial neighbor autocorrelation, coordinate-block variance fraction, and target standard deviation. DM1/RAI is then compared against the smoothness-adjusted random residual distribution. Passing this control argues that the DM1/RAI advantage is not explained only by being a more spatially smooth target.",
        "",
        "In this run, DM1/RAI does not pass that stricter control. The raw and coord+QC residual-target H&E effects remain positive, but their apparent random-module specificity is largely explained by target spatial smoothness and block structure.",
        "",
        "## Smoothness Model Summary",
        "",
        pd.DataFrame(summaries).round(4).to_markdown(index=False),
        "",
        "## Top Random Smoothness-Adjusted Residuals",
        "",
        random_top.round(4).to_markdown(index=False),
    ]
    (OUT / "GSE230424_SPATIAL_AUTOCORR_SPECIFICITY_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "GSE230424_SPATIAL_AUTOCORR_SPECIFICITY_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scores, results = load_inputs()
    metrics = compute_metrics(scores)
    results = results.copy()
    results["class"] = results["target"].map(classify)
    merged = results.merge(metrics, on=["target", "mode"], how="left")
    adjusted_parts = []
    summaries = []
    for mode in ["raw_smoothed", "coord_qc_residual_target"]:
        part, summary = fit_adjustment(merged, mode)
        adjusted_parts.append(part)
        summaries.append(summary)
    adjusted = pd.concat(adjusted_parts, ignore_index=True)
    metrics.to_csv(OUT / "gse230424_spatial_autocorr_target_metrics.tsv", sep="\t", index=False, na_rep="NA")
    adjusted.to_csv(OUT / "gse230424_spatial_autocorr_adjusted_specificity.tsv", sep="\t", index=False, na_rep="NA")
    write_figure(adjusted, {})
    summary = write_report(adjusted, summaries)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
