#!/usr/bin/env python3
"""Spot-budget downsampling controls for Path2Space-inspired evidence."""
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
OUT = P2 / "analysis_supp/path2space_downsample_robustness_2026_05_09"

N_PER_GROUP = 200
N_REPS = 1000
TOP_FRAC = 0.10
RANDOM_SEED = 250509


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


def exact_top_mask(values: np.ndarray, groups: np.ndarray, frac: float = TOP_FRAC) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    mask = np.zeros(len(values), dtype=bool)
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        valid = idx[keep]
        if len(valid) < 5:
            continue
        k = max(1, int(np.ceil(len(valid) * frac)))
        order = np.argsort(vals[keep], kind="mergesort")
        top = valid[order[-k:]]
        mask[top] = True
    return mask


def assign_deciles(values: np.ndarray, groups: np.ndarray, n_deciles: int = 10) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    dec = np.full(len(values), np.nan)
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        valid = idx[keep]
        if len(valid) < n_deciles:
            continue
        order = np.argsort(vals[keep], kind="mergesort")
        ranks = np.empty(len(valid), dtype=int)
        ranks[order] = np.arange(1, len(valid) + 1)
        d = np.ceil(ranks / len(valid) * n_deciles).astype(int)
        dec[valid] = np.clip(d, 1, n_deciles)
    return dec


def load_evidence() -> list[dict[str, object]]:
    g250 = pd.read_csv(G250 / "path2space_spot_predictions.tsv.gz", sep="\t")
    g250_resid = pd.read_csv(G250_SPEC / "gse250521_actual_predictions.tsv.gz", sep="\t")
    g230 = pd.read_csv(G230 / "gse230424_pathology_predictions.tsv.gz", sep="\t")
    g230_resid = pd.read_csv(G230_RESID / "gse230424_residual_target_predictions.tsv.gz", sep="\t")
    return [
        {
            "evidence": "GSE250521 UNI raw DM1/RAI",
            "short": "G250 raw",
            "df": g250,
            "group": "sample_id",
            "obs": "obs_DM1_like_score_smooth8",
            "pred": "pred_DM1_like_score_smooth8",
        },
        {
            "evidence": "GSE250521 UNI coord+QC residual DM1/RAI",
            "short": "G250 residual",
            "df": g250_resid,
            "group": "sample_id",
            "obs": "obs_actual_coord_qc_residual",
            "pred": "pred_actual_coord_qc_residual",
        },
        {
            "evidence": "GSE230424 H&E raw DM1/low-RAI",
            "short": "G230 raw",
            "df": g230,
            "group": "sample",
            "obs": "obs_DM1_low_RAI_score_smooth8",
            "pred": "pred_DM1_low_RAI_score_HE_tile_features",
        },
        {
            "evidence": "GSE230424 H&E coord+QC residual DM1/low-RAI",
            "short": "G230 residual",
            "df": g230_resid,
            "group": "sample",
            "obs": "obs_resid_DM1_low_RAI_score_coord_qc",
            "pred": "pred_resid_DM1_low_RAI_score_HE_tile_features",
        },
    ]


def prepare_work(df: pd.DataFrame, group_col: str, obs_col: str, pred_col: str) -> pd.DataFrame:
    work = pd.DataFrame(
        {
            "group": df[group_col].astype(str).to_numpy(),
            "obs": df[obs_col].to_numpy(dtype=float),
            "pred": df[pred_col].to_numpy(dtype=float),
        }
    )
    work = work[np.isfinite(work["obs"]) & np.isfinite(work["pred"])].copy()
    return work.reset_index(drop=True)


def evaluate_subset(work: pd.DataFrame, idx: np.ndarray) -> dict[str, float]:
    sub = work.iloc[idx].reset_index(drop=True)
    groups = sub["group"].to_numpy(dtype=str)
    obs_centered = center_by_group(sub["obs"].to_numpy(dtype=float), groups)
    pred_centered = center_by_group(sub["pred"].to_numpy(dtype=float), groups)
    rho, p, n = safe_spearman(pred_centered, obs_centered)

    pred_top = exact_top_mask(sub["pred"].to_numpy(dtype=float), groups)
    obs_top = exact_top_mask(sub["obs"].to_numpy(dtype=float), groups)
    hits = int((pred_top & obs_top).sum())
    pred_total = int(pred_top.sum())
    obs_total = int(obs_top.sum())
    precision = hits / pred_total if pred_total else np.nan
    base = obs_total / len(sub) if len(sub) else np.nan
    lift = precision / base if base and np.isfinite(base) else np.nan

    deciles = assign_deciles(sub["pred"].to_numpy(dtype=float), groups)
    dec_df = pd.DataFrame({"decile": deciles, "obs_centered": obs_centered})
    means = dec_df.groupby("decile")["obs_centered"].mean()
    top_bottom = float(means.get(10, np.nan) - means.get(1, np.nan))
    dec_rho, _, _ = safe_spearman(deciles, obs_centered)
    return {
        "n": int(len(sub)),
        "sample_centered_rho": rho,
        "sample_centered_p": p,
        "top_decile_precision": precision,
        "top_decile_lift": lift,
        "top_minus_bottom": top_bottom,
        "decile_spearman": dec_rho,
    }


def sample_indices(work: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    chunks = []
    for _, sub in work.groupby("group", sort=True):
        idx = sub.index.to_numpy()
        replace = len(idx) <= N_PER_GROUP
        chunks.append(rng.choice(idx, size=N_PER_GROUP, replace=replace))
    return np.concatenate(chunks)


def summarize_reps(rep_df: pd.DataFrame, metric: str) -> dict[str, float]:
    vals = rep_df[metric].to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    return {
        f"{metric}_median": float(np.nanmedian(vals)),
        f"{metric}_ci_low": float(np.nanpercentile(vals, 2.5)),
        f"{metric}_ci_high": float(np.nanpercentile(vals, 97.5)),
        f"{metric}_positive_fraction": float(np.mean(vals > 0)),
    }


def run_downsample() -> tuple[pd.DataFrame, pd.DataFrame]:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    full_rows = []
    rep_rows = []
    for spec in load_evidence():
        evidence = str(spec["evidence"])
        short = str(spec["short"])
        work = prepare_work(spec["df"], str(spec["group"]), str(spec["obs"]), str(spec["pred"]))
        full = evaluate_subset(work, work.index.to_numpy())
        full.update(
            {
                "evidence": evidence,
                "short": short,
                "n_groups": int(work["group"].nunique()),
                "n_total": int(len(work)),
                "spot_budget_per_group": N_PER_GROUP,
            }
        )
        full_rows.append(full)
        for rep in range(N_REPS):
            idx = sample_indices(work, rng)
            row = evaluate_subset(work, idx)
            row.update({"evidence": evidence, "short": short, "replicate": rep + 1})
            rep_rows.append(row)
    full_df = pd.DataFrame(full_rows)
    reps_df = pd.DataFrame(rep_rows)
    summary_rows = []
    for _, full in full_df.iterrows():
        rep_df = reps_df[reps_df["evidence"] == full["evidence"]]
        row = {
            "evidence": full["evidence"],
            "short": full["short"],
            "n_groups": int(full["n_groups"]),
            "n_total": int(full["n_total"]),
            "spot_budget_per_group": N_PER_GROUP,
            "n_reps": N_REPS,
            "full_sample_centered_rho": float(full["sample_centered_rho"]),
            "full_top_minus_bottom": float(full["top_minus_bottom"]),
            "full_top_decile_lift": float(full["top_decile_lift"]),
        }
        for metric in ["sample_centered_rho", "top_minus_bottom", "top_decile_lift", "decile_spearman"]:
            row.update(summarize_reps(rep_df, metric))
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "path2space_downsample_robustness_summary.tsv", sep="\t", index=False, na_rep="NA")
    reps_df.to_csv(OUT / "path2space_downsample_robustness_replicates.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    return summary, reps_df


def write_figure(summary: pd.DataFrame, reps: pd.DataFrame) -> None:
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
    labels = {
        order[0]: "G250 raw",
        order[1]: "G250 residual",
        order[2]: "G230 raw",
        order[3]: "G230 residual",
    }
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.4))
    metrics = [
        ("sample_centered_rho", "A. Downsampled centered rho", "Sample-centered rho"),
        ("top_minus_bottom", "B. Decile top-bottom separation", "Top decile - bottom decile"),
        ("top_decile_lift", "C. Top-decile hotspot lift", "Lift over 10% baseline"),
        ("decile_spearman", "D. Decile rank correlation", "Spearman(predicted decile, observed)"),
    ]
    x = np.arange(len(order))
    for ax, (metric, title, ylabel) in zip(axes.ravel(), metrics):
        med = summary.set_index("evidence").loc[order, f"{metric}_median"].to_numpy(dtype=float)
        lo = summary.set_index("evidence").loc[order, f"{metric}_ci_low"].to_numpy(dtype=float)
        hi = summary.set_index("evidence").loc[order, f"{metric}_ci_high"].to_numpy(dtype=float)
        yerr = np.vstack([med - lo, hi - med])
        ax.bar(x, med, yerr=yerr, color=[colors[e] for e in order], capsize=4, alpha=0.92)
        if metric == "top_decile_lift":
            ax.axhline(1, color="#333", lw=0.9, ls="--")
        else:
            ax.axhline(0, color="#333", lw=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels([labels[e] for e in order], rotation=25, ha="right")
        ax.set_title(title)
        ax.set_ylabel(ylabel)
    fig.suptitle(f"Path2Space spot-budget robustness ({N_REPS} reps, {N_PER_GROUP} spots per slide/sample)", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_downsample_robustness.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_downsample_robustness.pdf")
    plt.close(fig)


def write_report(summary: pd.DataFrame) -> dict[str, float]:
    key = summary.set_index("evidence")
    result = {
        "spot_budget_per_group": N_PER_GROUP,
        "n_reps": N_REPS,
        "gse250521_raw_rho_median": float(key.loc["GSE250521 UNI raw DM1/RAI", "sample_centered_rho_median"]),
        "gse250521_resid_rho_median": float(key.loc["GSE250521 UNI coord+QC residual DM1/RAI", "sample_centered_rho_median"]),
        "gse230424_raw_rho_median": float(key.loc["GSE230424 H&E raw DM1/low-RAI", "sample_centered_rho_median"]),
        "gse230424_resid_rho_median": float(key.loc["GSE230424 H&E coord+QC residual DM1/low-RAI", "sample_centered_rho_median"]),
        "gse230424_raw_top_minus_bottom_median": float(key.loc["GSE230424 H&E raw DM1/low-RAI", "top_minus_bottom_median"]),
        "gse230424_resid_top_minus_bottom_median": float(key.loc["GSE230424 H&E coord+QC residual DM1/low-RAI", "top_minus_bottom_median"]),
    }
    lines = [
        "# Path2Space spot-budget downsampling robustness",
        "",
        "## Verdict",
        "",
        f"- Repeated downsampling/bootstrapping used **{N_REPS}** replicates with **{N_PER_GROUP} spots per slide/sample**.",
        f"- GSE250521 raw median centered rho: **{result['gse250521_raw_rho_median']:.3f}**.",
        f"- GSE250521 coord+QC residual median centered rho: **{result['gse250521_resid_rho_median']:.3f}**.",
        f"- GSE230424 raw median centered rho after sample-size equalization: **{result['gse230424_raw_rho_median']:.3f}**.",
        f"- GSE230424 coord+QC residual median centered rho after sample-size equalization: **{result['gse230424_resid_rho_median']:.3f}**.",
        "",
        "## Interpretation",
        "",
        "This control asks whether the stronger GSE230424 result is mainly a large-spot-count artifact. Sampling the external cohort down to the same 200-spot-per-sample budget preserves a strong raw H&E-to-DM1/low-RAI signal and a weaker but positive residual signal.",
        "",
        "## Summary",
        "",
        summary.round(4).to_markdown(index=False),
    ]
    (OUT / "PATH2SPACE_DOWNSAMPLE_ROBUSTNESS_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "PATH2SPACE_DOWNSAMPLE_ROBUSTNESS_SUMMARY.json").write_text(json.dumps(result, indent=2) + "\n")
    (OUT / "SUMMARY.md").write_text(
        "# Path2Space Downsample Robustness\n\n"
        "## Verdict\n\n"
        f"- Spot budget: {N_PER_GROUP} per slide/sample; replicates: {N_REPS}.\n"
        f"- GSE230424 raw median rho after equalization: {result['gse230424_raw_rho_median']:.3f}.\n"
        f"- GSE230424 residual median rho after equalization: {result['gse230424_resid_rho_median']:.3f}.\n"
        "- Use as a sample-size robustness support layer; keep QC/smoothness caveats.\n"
    )
    return result


def main() -> None:
    summary, reps = run_downsample()
    write_figure(summary, reps)
    result = write_report(summary)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
