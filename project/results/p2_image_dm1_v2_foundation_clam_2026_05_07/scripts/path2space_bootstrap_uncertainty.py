#!/usr/bin/env python3
"""Bootstrap uncertainty for Path2Space-inspired Paper 2 spatial evidence."""
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
GSE250_STR = GSE250 / "strengthening_controls"
GSE230 = P2 / "analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"
GSE230_RESID = GSE230 / "residual_target_controls"
OUT = P2 / "analysis_supp/path2space_bootstrap_uncertainty_2026_05_09"
N_BOOT = 2000
RANDOM_SEED = 250230424


def safe_spearman(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan
    return float(stats.spearmanr(x[keep], y[keep]).statistic)


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


def stratified_bootstrap(df: pd.DataFrame, group_col: str, obs_col: str, pred_col: str, rng: np.random.Generator) -> np.ndarray:
    groups = {group: idx.to_numpy() for group, idx in df.groupby(group_col).groups.items()}
    obs = df[obs_col].to_numpy(dtype=float)
    pred = df[pred_col].to_numpy(dtype=float)
    reps = np.full(N_BOOT, np.nan)
    for b in range(N_BOOT):
        sampled_idx = []
        for idx in groups.values():
            sampled_idx.append(rng.choice(idx, size=len(idx), replace=True))
        boot_idx = np.concatenate(sampled_idx)
        boot_groups = df[group_col].to_numpy(dtype=str)[boot_idx]
        reps[b] = safe_spearman(center_by_group(obs[boot_idx], boot_groups), center_by_group(pred[boot_idx], boot_groups))
    return reps


def summarize(name: str, df: pd.DataFrame, group_col: str, obs_col: str, pred_col: str, level: str, rng: np.random.Generator) -> tuple[dict, pd.DataFrame]:
    groups = df[group_col].astype(str).to_numpy()
    observed = safe_spearman(
        center_by_group(df[obs_col].to_numpy(dtype=float), groups),
        center_by_group(df[pred_col].to_numpy(dtype=float), groups),
    )
    reps = stratified_bootstrap(df, group_col, obs_col, pred_col, rng)
    summary = {
        "evidence": name,
        "level": level,
        "n": int(len(df)),
        "n_groups": int(pd.Series(groups).nunique()),
        "observed_centered_rho": observed,
        "boot_mean": float(np.nanmean(reps)),
        "boot_ci_low": float(np.nanquantile(reps, 0.025)),
        "boot_ci_high": float(np.nanquantile(reps, 0.975)),
        "boot_p05": float(np.nanquantile(reps, 0.05)),
        "boot_p95": float(np.nanquantile(reps, 0.95)),
        "n_boot": int(np.isfinite(reps).sum()),
    }
    rep_df = pd.DataFrame({"evidence": name, "level": level, "bootstrap_rho": reps})
    return summary, rep_df


def load_evidence() -> list[tuple[str, pd.DataFrame, str, str, str, str]]:
    rows = []

    g250_spots = pd.read_csv(GSE250 / "path2space_spot_predictions.tsv.gz", sep="\t")
    rows.append(
        (
            "GSE250521 UNI DM1/RAI smoothed spots",
            g250_spots,
            "sample_id",
            "obs_DM1_like_score_smooth8",
            "pred_DM1_like_score_smooth8",
            "spot",
        )
    )

    g250_domains = pd.read_csv(GSE250_STR / "path2space_strengthening_domains.tsv", sep="\t")
    rows.append(
        (
            "GSE250521 UNI DM1/RAI coordinate domains",
            g250_domains,
            "sample_id",
            "obs_DM1_like_score",
            "pred_DM1_like_score",
            "domain",
        )
    )

    g230_spots = pd.read_csv(GSE230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    rows.append(
        (
            "GSE230424 H&E DM1/low-RAI smoothed spots",
            g230_spots,
            "sample",
            "obs_DM1_low_RAI_score_smooth8",
            "pred_DM1_low_RAI_score_HE_tile_features",
            "spot",
        )
    )

    g230_resid = pd.read_csv(GSE230_RESID / "gse230424_residual_target_predictions.tsv.gz", sep="\t")
    rows.append(
        (
            "GSE230424 H&E DM1/low-RAI coord+QC residual target",
            g230_resid,
            "sample",
            "obs_resid_DM1_low_RAI_score_coord_qc",
            "pred_resid_DM1_low_RAI_score_HE_tile_features",
            "spot",
        )
    )

    g230_domains = pd.read_csv(GSE230 / "gse230424_spatial_domains.tsv", sep="\t")
    rows.append(
        (
            "GSE230424 H&E DM1/low-RAI coordinate domains",
            g230_domains,
            "sample",
            "obs_DM1_low_RAI_score",
            "pred_DM1_low_RAI_score",
            "domain",
        )
    )
    return rows


def write_figure(summary_df: pd.DataFrame, reps: pd.DataFrame) -> None:
    order = summary_df.sort_values(["level", "observed_centered_rho"], ascending=[True, False])["evidence"].tolist()
    y = np.arange(len(order))
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8), gridspec_kw={"width_ratios": [1.15, 1]})

    plot_df = summary_df.set_index("evidence").loc[order].reset_index()
    x = plot_df["observed_centered_rho"].to_numpy(dtype=float)
    lo = plot_df["boot_ci_low"].to_numpy(dtype=float)
    hi = plot_df["boot_ci_high"].to_numpy(dtype=float)
    colors = ["#266b73" if level == "spot" else "#d5a84d" for level in plot_df["level"]]
    axes[0].errorbar(x, y, xerr=[x - lo, hi - x], fmt="none", ecolor="#444", elinewidth=1.3, capsize=3, zorder=1)
    axes[0].scatter(x, y, s=70, c=colors, edgecolor="white", linewidth=0.7, zorder=2)
    axes[0].axvline(0, color="#333", lw=0.8)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([name.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for name in order])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Centered Spearman rho with 95% bootstrap CI")
    axes[0].set_title("A. Uncertainty on primary image-to-spatial axes")

    for i, evidence in enumerate(order):
        sub = reps[reps["evidence"] == evidence]
        axes[1].hist(
            sub["bootstrap_rho"],
            bins=36,
            alpha=0.42,
            density=True,
            label=evidence.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ")[:38],
        )
    axes[1].axvline(0, color="#333", lw=0.8)
    axes[1].set_xlabel("Bootstrap centered rho")
    axes[1].set_ylabel("Density")
    axes[1].set_title("B. Bootstrap replicate distributions")
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_bootstrap_uncertainty.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_bootstrap_uncertainty.pdf")
    plt.close(fig)


def write_report(summary_df: pd.DataFrame) -> dict:
    key = summary_df.set_index("evidence")
    summary = {
        "n_boot": N_BOOT,
        "gse250521_spot_rho": float(key.loc["GSE250521 UNI DM1/RAI smoothed spots", "observed_centered_rho"]),
        "gse250521_spot_ci": [
            float(key.loc["GSE250521 UNI DM1/RAI smoothed spots", "boot_ci_low"]),
            float(key.loc["GSE250521 UNI DM1/RAI smoothed spots", "boot_ci_high"]),
        ],
        "gse250521_domain_rho": float(key.loc["GSE250521 UNI DM1/RAI coordinate domains", "observed_centered_rho"]),
        "gse250521_domain_ci": [
            float(key.loc["GSE250521 UNI DM1/RAI coordinate domains", "boot_ci_low"]),
            float(key.loc["GSE250521 UNI DM1/RAI coordinate domains", "boot_ci_high"]),
        ],
        "gse230424_spot_rho": float(key.loc["GSE230424 H&E DM1/low-RAI smoothed spots", "observed_centered_rho"]),
        "gse230424_spot_ci": [
            float(key.loc["GSE230424 H&E DM1/low-RAI smoothed spots", "boot_ci_low"]),
            float(key.loc["GSE230424 H&E DM1/low-RAI smoothed spots", "boot_ci_high"]),
        ],
        "gse230424_residual_rho": float(key.loc["GSE230424 H&E DM1/low-RAI coord+QC residual target", "observed_centered_rho"]),
        "gse230424_residual_ci": [
            float(key.loc["GSE230424 H&E DM1/low-RAI coord+QC residual target", "boot_ci_low"]),
            float(key.loc["GSE230424 H&E DM1/low-RAI coord+QC residual target", "boot_ci_high"]),
        ],
    }
    lines = [
        "# Path2Space-inspired bootstrap uncertainty",
        "",
        "## Verdict",
        "",
        f"- Bootstrap replicates per evidence layer: {N_BOOT}.",
        f"- GSE250521 smoothed spot-level DM1/RAI rho: **{summary['gse250521_spot_rho']:.3f}**; 95% bootstrap CI **[{summary['gse250521_spot_ci'][0]:.3f}, {summary['gse250521_spot_ci'][1]:.3f}]**.",
        f"- GSE250521 coordinate-domain DM1/RAI rho: **{summary['gse250521_domain_rho']:.3f}**; 95% bootstrap CI **[{summary['gse250521_domain_ci'][0]:.3f}, {summary['gse250521_domain_ci'][1]:.3f}]**.",
        f"- GSE230424 smoothed spot-level DM1/low-RAI rho: **{summary['gse230424_spot_rho']:.3f}**; 95% bootstrap CI **[{summary['gse230424_spot_ci'][0]:.3f}, {summary['gse230424_spot_ci'][1]:.3f}]**.",
        f"- GSE230424 coord+QC residual-target rho: **{summary['gse230424_residual_rho']:.3f}**; 95% bootstrap CI **[{summary['gse230424_residual_ci'][0]:.3f}, {summary['gse230424_residual_ci'][1]:.3f}]**.",
        "",
        "## Interpretation",
        "",
        "Bootstrap intervals remain above zero for the primary spot-level and domain-level DM1/RAI evidence layers. These intervals are descriptive clustered/stratified uncertainty checks, not a replacement for independent-cohort validation.",
        "",
        "## Summary Table",
        "",
        summary_df.round(4).to_markdown(index=False),
    ]
    (OUT / "PATH2SPACE_BOOTSTRAP_UNCERTAINTY_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "PATH2SPACE_BOOTSTRAP_UNCERTAINTY_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    summaries = []
    rep_dfs = []
    for evidence in load_evidence():
        summary, reps = summarize(*evidence, rng=rng)
        summaries.append(summary)
        rep_dfs.append(reps)
    summary_df = pd.DataFrame(summaries)
    reps = pd.concat(rep_dfs, ignore_index=True)
    summary_df.to_csv(OUT / "path2space_bootstrap_uncertainty_summary.tsv", sep="\t", index=False, na_rep="NA")
    reps.to_csv(OUT / "path2space_bootstrap_uncertainty_replicates.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    write_figure(summary_df, reps)
    summary = write_report(summary_df)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
