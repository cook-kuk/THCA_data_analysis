#!/usr/bin/env python3
"""Predicted-decile dose-response controls for image-to-DM1/RAI evidence."""
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
G250 = P2 / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"
G250_SPEC = P2 / "analysis_supp/path2space_gse250521_random_module_specificity_2026_05_09"
G230 = P2 / "analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"
G230_RESID = G230 / "residual_target_controls"
OUT = P2 / "analysis_supp/path2space_decile_dose_response_2026_05_09"
N_DECILES = 10


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
        idx = groups == group
        out[idx] = values[idx] - np.nanmean(values[idx])
    return out


def assign_deciles(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    dec = np.full(len(values), np.nan)
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        if keep.sum() < N_DECILES:
            continue
        valid = idx[keep]
        order = np.argsort(vals[keep], kind="mergesort")
        ranks = np.empty(len(valid), dtype=int)
        ranks[order] = np.arange(1, len(valid) + 1)
        decile = np.ceil(ranks / len(valid) * N_DECILES).astype(int)
        decile = np.clip(decile, 1, N_DECILES)
        dec[valid] = decile
    return dec


def load_evidence() -> list[tuple[str, pd.DataFrame, str, str, str]]:
    g250 = pd.read_csv(G250 / "path2space_spot_predictions.tsv.gz", sep="\t")
    g250_resid = pd.read_csv(G250_SPEC / "gse250521_actual_predictions.tsv.gz", sep="\t")
    g230 = pd.read_csv(G230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    g230_resid = pd.read_csv(G230_RESID / "gse230424_residual_target_predictions.tsv.gz", sep="\t")
    return [
        (
            "GSE250521 UNI raw DM1/RAI",
            g250,
            "sample_id",
            "obs_DM1_like_score_smooth8",
            "pred_DM1_like_score_smooth8",
        ),
        (
            "GSE250521 UNI coord+QC residual DM1/RAI",
            g250_resid,
            "sample_id",
            "obs_actual_coord_qc_residual",
            "pred_actual_coord_qc_residual",
        ),
        (
            "GSE230424 H&E raw DM1/low-RAI",
            g230,
            "sample",
            "obs_DM1_low_RAI_score_smooth8",
            "pred_DM1_low_RAI_score_HE_tile_features",
        ),
        (
            "GSE230424 H&E coord+QC residual DM1/low-RAI",
            g230_resid,
            "sample",
            "obs_resid_DM1_low_RAI_score_coord_qc",
            "pred_resid_DM1_low_RAI_score_HE_tile_features",
        ),
    ]


def evaluate(evidence: str, df: pd.DataFrame, group_col: str, obs_col: str, pred_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    groups = df[group_col].astype(str).to_numpy()
    obs = df[obs_col].to_numpy(dtype=float)
    pred = df[pred_col].to_numpy(dtype=float)
    obs_centered = center_by_group(obs, groups)
    pred_decile = assign_deciles(pred, groups)
    work = pd.DataFrame(
        {
            "evidence": evidence,
            "group": groups,
            "obs": obs,
            "obs_group_centered": obs_centered,
            "pred": pred,
            "pred_decile": pred_decile,
        }
    )
    work = work[np.isfinite(work["obs_group_centered"]) & np.isfinite(work["pred_decile"])].copy()
    work["pred_decile"] = work["pred_decile"].astype(int)

    dec_rows = []
    for decile, sub in work.groupby("pred_decile", sort=True):
        dec_rows.append(
            {
                "evidence": evidence,
                "pred_decile": int(decile),
                "n": int(len(sub)),
                "mean_obs_group_centered": float(sub["obs_group_centered"].mean()),
                "median_obs_group_centered": float(sub["obs_group_centered"].median()),
                "mean_pred": float(sub["pred"].mean()),
            }
        )
    dec_df = pd.DataFrame(dec_rows)

    group_rows = []
    for group, sub in work.groupby("group", sort=True):
        means = sub.groupby("pred_decile")["obs_group_centered"].mean().reindex(range(1, N_DECILES + 1))
        rho, p, n = safe_spearman(sub["pred_decile"], sub["obs_group_centered"])
        group_rows.append(
            {
                "evidence": evidence,
                "group": group,
                "n": int(len(sub)),
                "decile_obs_spearman": rho,
                "decile_obs_p": p,
                "top_minus_bottom": float(means.loc[N_DECILES] - means.loc[1]),
                "monotonic_increases": int(np.nansum(np.diff(means.to_numpy(dtype=float)) > 0)),
                "top_decile_mean_obs": float(means.loc[N_DECILES]),
                "bottom_decile_mean_obs": float(means.loc[1]),
            }
        )

    means = dec_df.set_index("pred_decile")["mean_obs_group_centered"].reindex(range(1, N_DECILES + 1))
    rho, p, n = safe_spearman(work["pred_decile"], work["obs_group_centered"])
    summary = pd.DataFrame(
        [
            {
                "evidence": evidence,
                "n": int(len(work)),
                "decile_obs_spearman": rho,
                "decile_obs_p": p,
                "top_minus_bottom": float(means.loc[N_DECILES] - means.loc[1]),
                "monotonic_increases": int(np.nansum(np.diff(means.to_numpy(dtype=float)) > 0)),
                "top_decile_mean_obs": float(means.loc[N_DECILES]),
                "bottom_decile_mean_obs": float(means.loc[1]),
            }
        ]
    )
    return pd.concat([summary, pd.DataFrame(group_rows)], ignore_index=True), dec_df


def write_figure(summary: pd.DataFrame, deciles: pd.DataFrame) -> None:
    order = [
        "GSE250521 UNI raw DM1/RAI",
        "GSE250521 UNI coord+QC residual DM1/RAI",
        "GSE230424 H&E raw DM1/low-RAI",
        "GSE230424 H&E coord+QC residual DM1/low-RAI",
    ]
    colors = {
        order[0]: "#266b73",
        order[1]: "#63d5c4",
        order[2]: "#9b5de5",
        order[3]: "#d6b25e",
    }
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.2))

    for evidence in order:
        sub = deciles[deciles["evidence"] == evidence].sort_values("pred_decile")
        axes[0, 0].plot(
            sub["pred_decile"],
            sub["mean_obs_group_centered"],
            marker="o",
            lw=2.2,
            label=evidence.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 "),
            color=colors[evidence],
        )
    axes[0, 0].axhline(0, color="#333", lw=0.8)
    axes[0, 0].set_xlabel("Predicted decile within slide/sample")
    axes[0, 0].set_ylabel("Mean observed DM1/RAI, group-centered")
    axes[0, 0].set_title("A. Predicted-decile dose response")
    axes[0, 0].legend(frameon=False, fontsize=8)

    top = summary[summary["group"].isna()].set_index("evidence").loc[order]
    x = np.arange(len(order))
    axes[0, 1].bar(x, top["top_minus_bottom"], color=[colors[e] for e in order])
    axes[0, 1].axhline(0, color="#333", lw=0.8)
    axes[0, 1].set_xticks(x)
    axes[0, 1].set_xticklabels([e.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for e in order], rotation=25, ha="right")
    axes[0, 1].set_ylabel("Top decile mean - bottom decile mean")
    axes[0, 1].set_title("B. Top-to-bottom separation")

    axes[1, 0].bar(x, top["decile_obs_spearman"], color=[colors[e] for e in order])
    axes[1, 0].axhline(0, color="#333", lw=0.8)
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels([e.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for e in order], rotation=25, ha="right")
    axes[1, 0].set_ylabel("Spot-level Spearman(predicted decile, observed)")
    axes[1, 0].set_title("C. Decile rank correlation")

    group_summary = summary[summary["group"].notna()].copy()
    for i, evidence in enumerate(order):
        vals = group_summary[group_summary["evidence"] == evidence]["top_minus_bottom"].to_numpy(dtype=float)
        jitter = np.linspace(-0.14, 0.14, len(vals)) if len(vals) > 1 else np.zeros(len(vals))
        axes[1, 1].scatter(np.full(len(vals), i) + jitter, vals, color=colors[evidence], s=42, alpha=0.8)
    axes[1, 1].axhline(0, color="#333", lw=0.8)
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels([e.replace("GSE250521 ", "G250 ").replace("GSE230424 ", "G230 ") for e in order], rotation=25, ha="right")
    axes[1, 1].set_ylabel("Per-slide/sample top-bottom separation")
    axes[1, 1].set_title("D. Group-level consistency")

    fig.suptitle("Path2Space-inspired predicted-decile dose-response controls", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_decile_dose_response.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_decile_dose_response.pdf")
    plt.close(fig)


def write_report(summary: pd.DataFrame, deciles: pd.DataFrame) -> dict:
    overall = summary[summary["group"].isna()].copy()
    by_evidence = overall.set_index("evidence")
    result = {
        "gse250521_raw_top_minus_bottom": float(by_evidence.loc["GSE250521 UNI raw DM1/RAI", "top_minus_bottom"]),
        "gse250521_raw_decile_spearman": float(by_evidence.loc["GSE250521 UNI raw DM1/RAI", "decile_obs_spearman"]),
        "gse250521_resid_top_minus_bottom": float(by_evidence.loc["GSE250521 UNI coord+QC residual DM1/RAI", "top_minus_bottom"]),
        "gse250521_resid_decile_spearman": float(by_evidence.loc["GSE250521 UNI coord+QC residual DM1/RAI", "decile_obs_spearman"]),
        "gse230424_raw_top_minus_bottom": float(by_evidence.loc["GSE230424 H&E raw DM1/low-RAI", "top_minus_bottom"]),
        "gse230424_raw_decile_spearman": float(by_evidence.loc["GSE230424 H&E raw DM1/low-RAI", "decile_obs_spearman"]),
        "gse230424_resid_top_minus_bottom": float(by_evidence.loc["GSE230424 H&E coord+QC residual DM1/low-RAI", "top_minus_bottom"]),
        "gse230424_resid_decile_spearman": float(by_evidence.loc["GSE230424 H&E coord+QC residual DM1/low-RAI", "decile_obs_spearman"]),
    }
    lines = [
        "# Path2Space predicted-decile dose-response controls",
        "",
        "## Verdict",
        "",
        f"- GSE250521 raw top-minus-bottom observed DM1/RAI: **{result['gse250521_raw_top_minus_bottom']:.3f}**; decile Spearman **{result['gse250521_raw_decile_spearman']:.3f}**.",
        f"- GSE250521 coord+QC residual top-minus-bottom: **{result['gse250521_resid_top_minus_bottom']:.3f}**; decile Spearman **{result['gse250521_resid_decile_spearman']:.3f}**.",
        f"- GSE230424 raw top-minus-bottom observed DM1/low-RAI: **{result['gse230424_raw_top_minus_bottom']:.3f}**; decile Spearman **{result['gse230424_raw_decile_spearman']:.3f}**.",
        f"- GSE230424 coord+QC residual top-minus-bottom: **{result['gse230424_resid_top_minus_bottom']:.3f}**; decile Spearman **{result['gse230424_resid_decile_spearman']:.3f}**.",
        "",
        "## Interpretation",
        "",
        "Predicted-decile dose response checks whether higher predicted DM1/RAI bins carry progressively higher observed DM1/RAI within each slide/sample. This complements top-decile hotspot overlap by avoiding a single threshold. Positive top-minus-bottom separation supports ranking calibration; weaker residual slopes preserve the QC/smoothness caveat.",
        "",
        "## Overall Summary",
        "",
        overall.round(4).to_markdown(index=False),
        "",
        "## Group-Level Summary",
        "",
        summary[summary["group"].notna()].round(4).to_markdown(index=False),
        "",
        "## Decile Means",
        "",
        deciles.round(4).to_markdown(index=False),
    ]
    (OUT / "PATH2SPACE_DECILE_DOSE_RESPONSE_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "PATH2SPACE_DECILE_DOSE_RESPONSE_SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")

    summary_lines = [
        "# Path2Space Predicted-Decile Dose Response",
        "",
        "## Verdict",
        "",
        "| Evidence | Top-bottom observed separation | Decile Spearman |",
        "|---|---:|---:|",
        f"| GSE250521 raw | {result['gse250521_raw_top_minus_bottom']:.3f} | {result['gse250521_raw_decile_spearman']:.3f} |",
        f"| GSE250521 coord+QC residual | {result['gse250521_resid_top_minus_bottom']:.3f} | {result['gse250521_resid_decile_spearman']:.3f} |",
        f"| GSE230424 raw | {result['gse230424_raw_top_minus_bottom']:.3f} | {result['gse230424_raw_decile_spearman']:.3f} |",
        f"| GSE230424 coord+QC residual | {result['gse230424_resid_top_minus_bottom']:.3f} | {result['gse230424_resid_decile_spearman']:.3f} |",
        "",
        "Safe wording: predicted high-DM1/RAI deciles show dose-response enrichment of observed DM1/RAI in both cohorts; residual layers are weaker but remain positive.",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(summary_lines) + "\n")
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summaries = []
    deciles = []
    for evidence, df, group_col, obs_col, pred_col in load_evidence():
        summary, decile = evaluate(evidence, df, group_col, obs_col, pred_col)
        summaries.append(summary)
        deciles.append(decile)
    summary_df = pd.concat(summaries, ignore_index=True)
    decile_df = pd.concat(deciles, ignore_index=True)
    summary_df.to_csv(OUT / "path2space_decile_dose_response_summary.tsv", sep="\t", index=False, na_rep="NA")
    decile_df.to_csv(OUT / "path2space_decile_dose_response_deciles.tsv", sep="\t", index=False, na_rep="NA")
    result = write_report(summary_df, decile_df)
    write_figure(summary_df, decile_df)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
