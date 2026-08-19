#!/usr/bin/env python3
"""GSE230424 DM1/RAI image signal specificity controls.

This extends the GSE230424 Path2Space-style analysis with two reviewer-facing
controls:

1. Panel leave-one-gene-out: the image-to-DM1/RAI signal should not collapse
   when any one of the eight thyroid-lineage genes is removed.
2. Expression/detection-matched random 8-gene modules: the DM1/RAI signal is
   compared against a null distribution of similarly expressed gene modules.
"""
from __future__ import annotations

import gzip
import json
import math
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler

from gse230424_pathology_thyroid_axis import (
    OUT as PARENT_OUT,
    PANEL_8,
    RANDOM_SEED,
    RAW,
    RIDGE_ALPHA,
    SMOOTH_K,
    center_by_group,
    read_barcodes,
    read_features,
    safe_spearman,
    sample_prefixes,
    spatial_smooth,
)


OUT = PARENT_OUT / "random_module_specificity"
N_RANDOM_MODULES = 250
MODULE_SIZE = 8
N_PCS = 32
RANDOM_POOL_NEAREST = 250

warnings.filterwarnings("ignore", category=PerformanceWarning)


def logit_clip(values: np.ndarray | pd.Series | float) -> np.ndarray:
    clipped = np.clip(np.asarray(values, dtype=float), 1e-4, 1 - 1e-4)
    return np.log(clipped / (1.0 - clipped))


def sample_from_prefix(prefix: str) -> str:
    return prefix.split("_")[-1]


def expression_metrics() -> pd.DataFrame:
    """Stream matrix-market files once to get gene-level expression metrics."""
    prefixes = sample_prefixes()
    first = read_features(prefixes[0])
    genes = first["gene"].astype(str).to_numpy()
    gene_to_idx: dict[str, int] = {}
    for i, gene in enumerate(genes):
        gene_to_idx.setdefault(gene, i)

    total_counts = np.zeros(len(genes), dtype=np.float64)
    detected_spots = np.zeros(len(genes), dtype=np.int64)
    n_spots = 0

    for prefix in prefixes:
        features = read_features(prefix)
        row_to_idx = {}
        for row_idx, gene in enumerate(features["gene"].astype(str), start=1):
            idx = gene_to_idx.get(gene)
            if idx is not None:
                row_to_idx[row_idx] = idx
        with gzip.open(RAW / f"{prefix}_matrix.mtx.gz", "rt") as handle:
            dims_seen = False
            for line in handle:
                if line.startswith("%"):
                    continue
                parts = line.strip().split()
                if not dims_seen:
                    _, n_cols, _ = map(int, parts)
                    n_spots += n_cols
                    dims_seen = True
                    continue
                row = int(parts[0])
                val = float(parts[2])
                idx = row_to_idx.get(row)
                if idx is None:
                    continue
                total_counts[idx] += val
                if val > 0:
                    detected_spots[idx] += 1

    metrics = pd.DataFrame(
        {
            "gene": genes,
            "total_counts": total_counts,
            "detected_spots": detected_spots,
        }
    )
    metrics = metrics.drop_duplicates("gene", keep="first")
    metrics["n_spots"] = n_spots
    metrics["detect_frac"] = metrics["detected_spots"] / max(n_spots, 1)
    metrics["mean_count"] = metrics["total_counts"] / max(n_spots, 1)
    metrics["log1p_mean_count"] = np.log1p(metrics["mean_count"])
    metrics.to_csv(OUT / "gse230424_gene_expression_metrics.tsv", sep="\t", index=False, na_rep="NA")
    return metrics


def eligible_universe(metrics: pd.DataFrame) -> pd.DataFrame:
    excluded_prefixes = ("MT-", "RPS", "RPL")
    excluded_genes = set(PANEL_8)
    universe = metrics.copy()
    universe = universe[
        (universe["detect_frac"] >= 0.01)
        & (universe["detect_frac"] <= 0.85)
        & (universe["total_counts"] >= 50)
        & (~universe["gene"].astype(str).str.startswith(excluded_prefixes))
        & (~universe["gene"].isin(excluded_genes))
    ].copy()
    universe = universe.reset_index(drop=True)
    if len(universe) < 500:
        raise RuntimeError(f"Eligible random-module universe too small: {len(universe)}")
    return universe


def make_matched_random_modules(metrics: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Build random modules matched to the eight panel genes by expression/detection."""
    universe = eligible_universe(metrics)
    panel_metrics = metrics.set_index("gene").loc[PANEL_8].reset_index()

    u_feat = np.column_stack(
        [
            np.log1p(universe["mean_count"].to_numpy(dtype=float)),
            logit_clip(universe["detect_frac"].to_numpy(dtype=float)),
        ]
    )
    scaler = StandardScaler()
    u_scaled = scaler.fit_transform(u_feat)

    pools: dict[str, list[str]] = {}
    for _, row in panel_metrics.iterrows():
        p_feat = np.array(
            [[
                math.log1p(float(row["mean_count"])),
                float(logit_clip(float(row["detect_frac"]))),
            ]]
        )
        p_scaled = scaler.transform(p_feat)
        dist = np.sqrt(((u_scaled - p_scaled) ** 2).sum(axis=1))
        nearest = np.argsort(dist)[: min(RANDOM_POOL_NEAREST, len(universe))]
        pools[str(row["gene"])] = universe.iloc[nearest]["gene"].astype(str).tolist()

    modules: dict[str, list[str]] = {}
    attempts = 0
    while len(modules) < N_RANDOM_MODULES and attempts < N_RANDOM_MODULES * 50:
        attempts += 1
        chosen: list[str] = []
        for panel_gene in PANEL_8:
            pool = [gene for gene in pools[panel_gene] if gene not in chosen]
            if not pool:
                break
            chosen.append(str(rng.choice(pool)))
        if len(set(chosen)) != MODULE_SIZE:
            continue
        module_name = f"random_matched_{len(modules) + 1:03d}"
        modules[module_name] = chosen

    if len(modules) < N_RANDOM_MODULES:
        raise RuntimeError(f"Only built {len(modules)} random modules")

    rows = []
    for module, genes in modules.items():
        for slot, gene in zip(PANEL_8, genes):
            rows.append({"module": module, "matched_panel_gene": slot, "random_gene": gene})
    module_df = pd.DataFrame(rows)
    module_df.to_csv(OUT / "gse230424_random_matched_modules.tsv", sep="\t", index=False)
    return module_df, modules


def stream_expression(selected_genes: list[str]) -> pd.DataFrame:
    selected = list(dict.fromkeys(selected_genes))
    selected_set = set(selected)
    gene_to_col = {gene: i for i, gene in enumerate(selected)}
    rows = []
    for prefix in sample_prefixes():
        sample = sample_from_prefix(prefix)
        features = read_features(prefix)
        barcodes = read_barcodes(prefix).to_numpy()
        row_to_col = {}
        for row_idx, gene in enumerate(features["gene"].astype(str), start=1):
            if gene in selected_set:
                row_to_col[row_idx] = gene_to_col[gene]

        counts = None
        total = None
        with gzip.open(RAW / f"{prefix}_matrix.mtx.gz", "rt") as handle:
            dims_seen = False
            for line in handle:
                if line.startswith("%"):
                    continue
                parts = line.strip().split()
                if not dims_seen:
                    _, n_cols, _ = map(int, parts)
                    counts = np.zeros((n_cols, len(selected)), dtype=np.float32)
                    total = np.zeros(n_cols, dtype=np.float64)
                    dims_seen = True
                    continue
                row = int(parts[0])
                col = int(parts[1]) - 1
                val = float(parts[2])
                total[col] += val
                target_col = row_to_col.get(row)
                if target_col is not None:
                    counts[col, target_col] += val

        assert counts is not None and total is not None
        denom = np.maximum(total, 1.0)
        log_cp10k = np.log1p(counts / denom[:, None] * 10000.0)
        part = pd.DataFrame(log_cp10k, columns=[f"g_{gene}" for gene in selected])
        part.insert(0, "barcode", barcodes)
        part.insert(0, "sample", sample)
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def zscore_matrix(expr: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, dict[str, int]]:
    cols = [f"g_{gene}" for gene in genes]
    x = expr[cols].to_numpy(dtype=np.float32)
    mean = np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    z = (x - mean[None, :]) / sd[None, :]
    return z, {gene: i for i, gene in enumerate(genes)}


def build_score_table(meta: pd.DataFrame, expr: pd.DataFrame, modules: dict[str, list[str]]) -> tuple[pd.DataFrame, list[str], list[str]]:
    meta_cols = [
        "sample",
        "disease_group",
        "geo_accession",
        "barcode",
        "array_row",
        "array_col",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ]
    merged = meta[meta_cols].merge(expr, on=["sample", "barcode"], how="inner")
    selected_genes = [c.replace("g_", "", 1) for c in expr.columns if c.startswith("g_")]
    z, gene_idx = zscore_matrix(merged, selected_genes)

    scores = merged[[
        "sample",
        "disease_group",
        "geo_accession",
        "barcode",
        "array_row",
        "array_col",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ]].copy()

    panel_cols = [gene_idx[gene] for gene in PANEL_8]
    rai8 = np.nanmean(z[:, panel_cols], axis=1)
    scores["actual_DM1_low_RAI"] = -rai8

    loo_names = []
    for removed in PANEL_8:
        keep = [gene_idx[gene] for gene in PANEL_8 if gene != removed]
        name = f"loo_drop_{removed}"
        scores[name] = -np.nanmean(z[:, keep], axis=1)
        loo_names.append(name)

    random_names = []
    for module, genes in modules.items():
        cols = [gene_idx[gene] for gene in genes]
        name = f"null_{module}"
        scores[name] = -np.nanmean(z[:, cols], axis=1)
        random_names.append(name)

    score_cols = ["actual_DM1_low_RAI"] + loo_names + random_names
    scores = spatial_smooth(scores, score_cols, k=SMOOTH_K)
    keep_cols = list(scores.columns[:9]) + [f"{col}_smooth{SMOOTH_K}" for col in score_cols]
    out = scores[keep_cols].copy()
    out.to_csv(OUT / "gse230424_random_module_spot_scores.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    return out, loo_names, random_names


def residualize_targets(y: np.ndarray, covars: np.ndarray, samples: np.ndarray) -> np.ndarray:
    out = np.full_like(y, np.nan, dtype=np.float64)
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        x = np.asarray(covars, dtype=float)[idx]
        keep_x = np.isfinite(x).all(axis=1)
        if keep_x.sum() < 30:
            continue
        scaler = StandardScaler()
        xs = scaler.fit_transform(x[keep_x])
        for j in range(y.shape[1]):
            vals = y[idx, j]
            keep = keep_x & np.isfinite(vals)
            if keep.sum() < 30:
                continue
            xs_j = scaler.fit_transform(x[keep])
            fit = LinearRegression().fit(xs_j, vals[keep])
            out[idx[keep], j] = vals[keep] - fit.predict(xs_j)
    return out


def loso_predict_many(x: np.ndarray, y: np.ndarray, samples: np.ndarray, use_pca: bool = True) -> np.ndarray:
    pred = np.full(y.shape, np.nan, dtype=np.float64)
    finite_x = np.isfinite(x).all(axis=1)
    for held in sorted(pd.unique(samples)):
        train = (samples != held) & finite_x & np.isfinite(y).all(axis=1)
        test = (samples == held) & finite_x
        if train.sum() < 30 or test.sum() == 0:
            continue
        scaler = StandardScaler()
        x_train = scaler.fit_transform(x[train])
        x_test = scaler.transform(x[test])
        if use_pca:
            n_comp = min(N_PCS, x_train.shape[0] - 1, x_train.shape[1])
            pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
            x_train = pca.fit_transform(x_train)
            x_test = pca.transform(x_test)
        mean = np.nanmean(y[train], axis=0)
        sd = np.nanstd(y[train], axis=0)
        sd[~np.isfinite(sd) | (sd == 0)] = 1.0
        model = Ridge(alpha=RIDGE_ALPHA)
        model.fit(x_train, (y[train] - mean[None, :]) / sd[None, :])
        pred[test] = model.predict(x_test)
    return pred


def summarize_predictions(names: list[str], y: np.ndarray, pred: np.ndarray, samples: np.ndarray, mode: str) -> pd.DataFrame:
    y_center = np.column_stack([center_by_group(y[:, j], samples) for j in range(y.shape[1])])
    pred_center = np.column_stack([center_by_group(pred[:, j], samples) for j in range(pred.shape[1])])
    rows = []
    for j, name in enumerate(names):
        rho, p, n = safe_spearman(y[:, j], pred[:, j])
        crho, cp, cn = safe_spearman(y_center[:, j], pred_center[:, j])
        per_sample = []
        for sample in pd.unique(samples):
            idx = np.where(samples == sample)[0]
            srho, _, _ = safe_spearman(y[idx, j], pred[idx, j])
            per_sample.append(srho)
        rows.append(
            {
                "target": name,
                "mode": mode,
                "n": n,
                "pooled_rho": rho,
                "sample_centered_rho": crho,
                "sample_centered_p": cp,
                "median_sample_rho": float(np.nanmedian(per_sample)),
                "samples_rho_gt_0_2": int(np.nansum(np.asarray(per_sample) > 0.2)),
            }
        )
    return pd.DataFrame(rows)


def empirical_upper_p(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    return float((1.0 + np.sum(null >= observed)) / (1.0 + len(null)))


def percentile(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    return float(100.0 * np.mean(null <= observed))


def write_figures(result_df: pd.DataFrame, loo_df: pd.DataFrame, summary: dict) -> None:
    random_raw = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "raw_smoothed")]
    random_res = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "coord_qc_residual_target")]
    actual = result_df[result_df["class"] == "actual_dm1"].set_index("mode")

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5))
    axes[0, 0].hist(random_raw["sample_centered_rho"], bins=32, color="#9aa4ad", edgecolor="white")
    axes[0, 0].axvline(actual.loc["raw_smoothed", "sample_centered_rho"], color="#b44d4d", lw=2.5, label="DM1/RAI")
    axes[0, 0].set_title("A. Raw H&E predictability vs matched random modules")
    axes[0, 0].set_xlabel("Sample-centered LOSO rho")
    axes[0, 0].set_ylabel("Random modules")
    axes[0, 0].legend(frameon=False)

    axes[0, 1].hist(random_res["sample_centered_rho"], bins=32, color="#9aa4ad", edgecolor="white")
    axes[0, 1].axvline(actual.loc["coord_qc_residual_target", "sample_centered_rho"], color="#b44d4d", lw=2.5, label="DM1/RAI")
    axes[0, 1].set_title("B. Coord+QC residual target vs matched random modules")
    axes[0, 1].set_xlabel("Sample-centered LOSO rho")
    axes[0, 1].legend(frameon=False)

    loo_plot = loo_df.pivot(index="target", columns="mode", values="sample_centered_rho")
    loo_plot = loo_plot.sort_values("raw_smoothed", ascending=True)
    y = np.arange(len(loo_plot))
    axes[1, 0].barh(y - 0.18, loo_plot["raw_smoothed"], height=0.35, color="#266b73", label="raw")
    axes[1, 0].barh(y + 0.18, loo_plot["coord_qc_residual_target"], height=0.35, color="#d5a84d", label="coord+QC residual")
    axes[1, 0].axvline(0, color="#333", lw=0.8)
    axes[1, 0].set_yticks(y)
    axes[1, 0].set_yticklabels([idx.replace("loo_drop_", "drop ") for idx in loo_plot.index])
    axes[1, 0].set_xlabel("Sample-centered LOSO rho")
    axes[1, 0].set_title("C. Panel leave-one-gene-out stability")
    axes[1, 0].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, -0.14), ncol=2)

    paired = random_raw[["target", "sample_centered_rho"]].rename(columns={"sample_centered_rho": "raw"}).merge(
        random_res[["target", "sample_centered_rho"]].rename(columns={"sample_centered_rho": "residual"}),
        on="target",
    )
    axes[1, 1].scatter(paired["raw"], paired["residual"], s=18, alpha=0.55, color="#777")
    axes[1, 1].scatter(
        [actual.loc["raw_smoothed", "sample_centered_rho"]],
        [actual.loc["coord_qc_residual_target", "sample_centered_rho"]],
        s=90,
        color="#b44d4d",
        label="DM1/RAI",
        zorder=5,
    )
    axes[1, 1].axhline(0, color="#333", lw=0.8)
    axes[1, 1].axvline(0, color="#333", lw=0.8)
    axes[1, 1].set_xlabel("Raw rho")
    axes[1, 1].set_ylabel("Coord+QC residual-target rho")
    axes[1, 1].set_title("D. Raw vs residual specificity")
    axes[1, 1].legend(frameon=False)

    fig.suptitle(
        "GSE230424 specificity controls: DM1/RAI image signal vs matched random modules",
        y=0.995,
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse230424_random_module_specificity.png", dpi=220)
    fig.savefig(OUT / "fig_gse230424_random_module_specificity.pdf")
    plt.close(fig)


def write_report(result_df: pd.DataFrame, module_df: pd.DataFrame, summary: dict) -> None:
    random_raw = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "raw_smoothed")]
    random_res = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "coord_qc_residual_target")]
    actual = result_df[result_df["class"] == "actual_dm1"].set_index("mode")
    loo = result_df[result_df["class"] == "panel_loo"].copy()
    top_raw = random_raw.sort_values("sample_centered_rho", ascending=False).head(12)
    top_res = random_res.sort_values("sample_centered_rho", ascending=False).head(12)

    lines = [
        "# GSE230424 random-module specificity controls",
        "",
        "## Verdict",
        "",
        f"- Matched random modules: {summary['n_random_modules']} expression/detection-matched 8-gene modules.",
        f"- Actual DM1/RAI raw H&E rho: **{summary['actual_raw_rho']:.3f}**; random-module percentile **{summary['actual_raw_percentile']:.1f}%**; empirical upper-tail p = **{summary['actual_raw_empirical_p']:.4f}**.",
        f"- Actual DM1/RAI coord+QC residual-target H&E rho: **{summary['actual_residual_rho']:.3f}**; random-module percentile **{summary['actual_residual_percentile']:.1f}%**; empirical upper-tail p = **{summary['actual_residual_empirical_p']:.4f}**.",
        f"- Panel leave-one-out raw rho range: **{summary['loo_raw_min']:.3f} to {summary['loo_raw_max']:.3f}**.",
        f"- Panel leave-one-out residual-target rho range: **{summary['loo_residual_min']:.3f} to {summary['loo_residual_max']:.3f}**.",
        "",
        "## Interpretation",
        "",
        "The control asks whether the H&E-to-DM1/RAI result is just another spatially smooth random gene module. A high percentile supports specificity; a lower percentile would force a broader tissue-state framing. Leave-one-out stability tests whether one extreme thyroid-lineage gene drives the result.",
        "",
        "## Actual And Leave-One-Out Results",
        "",
        result_df[result_df["class"].isin(["actual_dm1", "panel_loo"])].round(4).to_markdown(index=False),
        "",
        "## Random Module Null Summary",
        "",
        pd.DataFrame(
            [
                {
                    "mode": "raw_smoothed",
                    "null_mean": random_raw["sample_centered_rho"].mean(),
                    "null_median": random_raw["sample_centered_rho"].median(),
                    "null_p90": random_raw["sample_centered_rho"].quantile(0.90),
                    "null_p95": random_raw["sample_centered_rho"].quantile(0.95),
                    "actual": actual.loc["raw_smoothed", "sample_centered_rho"],
                },
                {
                    "mode": "coord_qc_residual_target",
                    "null_mean": random_res["sample_centered_rho"].mean(),
                    "null_median": random_res["sample_centered_rho"].median(),
                    "null_p90": random_res["sample_centered_rho"].quantile(0.90),
                    "null_p95": random_res["sample_centered_rho"].quantile(0.95),
                    "actual": actual.loc["coord_qc_residual_target", "sample_centered_rho"],
                },
            ]
        ).round(4).to_markdown(index=False),
        "",
        "## Top Random Raw Modules",
        "",
        top_raw.round(4).to_markdown(index=False),
        "",
        "## Top Random Residual Modules",
        "",
        top_res.round(4).to_markdown(index=False),
        "",
        "## Random Module Gene Map",
        "",
        "See `gse230424_random_matched_modules.tsv` for the full matched gene map.",
    ]
    (OUT / "GSE230424_RANDOM_MODULE_SPECIFICITY_REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED + 57)

    meta = pd.read_csv(PARENT_OUT / "gse230424_spot_module_scores.tsv.gz", sep="\t")
    feats = pd.read_csv(PARENT_OUT / "gse230424_he_tile_features.tsv.gz", sep="\t")

    metrics_path = OUT / "gse230424_gene_expression_metrics.tsv"
    if metrics_path.exists():
        metrics = pd.read_csv(metrics_path, sep="\t")
    else:
        metrics = expression_metrics()

    module_map_path = OUT / "gse230424_random_matched_modules.tsv"
    if module_map_path.exists():
        module_df = pd.read_csv(module_map_path, sep="\t")
        modules = {name: sub["random_gene"].astype(str).tolist() for name, sub in module_df.groupby("module", sort=True)}
    else:
        module_df, modules = make_matched_random_modules(metrics, rng)

    selected_genes = list(PANEL_8)
    for genes in modules.values():
        selected_genes.extend(genes)
    selected_genes = list(dict.fromkeys(selected_genes))

    expr_path = OUT / "gse230424_random_module_expression.tsv.gz"
    if expr_path.exists():
        expr = pd.read_csv(expr_path, sep="\t")
    else:
        expr = stream_expression(selected_genes)
        expr.to_csv(expr_path, sep="\t", index=False, compression="gzip", na_rep="NA")

    scores, loo_names, random_names = build_score_table(meta, expr, modules)
    score_cols = ["actual_DM1_low_RAI"] + loo_names + random_names
    smooth_cols = [f"{col}_smooth{SMOOTH_K}" for col in score_cols]

    merged = scores.merge(feats, on=["sample", "barcode", "array_row", "array_col"], how="inner")
    samples = merged["sample"].astype(str).to_numpy()
    he_cols = [c for c in feats.columns if c.startswith("r")]
    x_he = merged[he_cols].to_numpy(dtype=float)
    y_raw = merged[smooth_cols].to_numpy(dtype=float)

    coord_qc = np.column_stack(
        [
            merged[["array_row", "array_col"]].to_numpy(dtype=float),
            np.log1p(merged["total_counts"].to_numpy(dtype=float)),
            np.log1p(merged["n_genes_by_counts"].to_numpy(dtype=float)),
            merged["pct_counts_mt"].to_numpy(dtype=float),
        ]
    )
    y_resid = residualize_targets(y_raw, coord_qc, samples)

    pred_raw = loso_predict_many(x_he, y_raw, samples, use_pca=True)
    pred_resid = loso_predict_many(x_he, y_resid, samples, use_pca=True)

    names = score_cols
    raw_df = summarize_predictions(names, y_raw, pred_raw, samples, "raw_smoothed")
    res_df = summarize_predictions(names, y_resid, pred_resid, samples, "coord_qc_residual_target")
    result_df = pd.concat([raw_df, res_df], ignore_index=True)

    def classify(name: str) -> str:
        if name == "actual_DM1_low_RAI":
            return "actual_dm1"
        if name.startswith("loo_drop_"):
            return "panel_loo"
        return "random_matched"

    result_df["class"] = result_df["target"].map(classify)
    result_df.to_csv(OUT / "gse230424_random_module_specificity_results.tsv", sep="\t", index=False, na_rep="NA")

    pred_keep = merged[["sample", "barcode", "array_row", "array_col"]].copy()
    for name, obs, pred in [
        ("actual_raw", y_raw[:, 0], pred_raw[:, 0]),
        ("actual_coord_qc_residual", y_resid[:, 0], pred_resid[:, 0]),
    ]:
        pred_keep[f"obs_{name}"] = obs
        pred_keep[f"pred_{name}"] = pred
    pred_keep.to_csv(OUT / "gse230424_random_module_actual_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    random_raw = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "raw_smoothed")]["sample_centered_rho"].to_numpy(dtype=float)
    random_res = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "coord_qc_residual_target")]["sample_centered_rho"].to_numpy(dtype=float)
    actual_raw = float(result_df[(result_df["class"] == "actual_dm1") & (result_df["mode"] == "raw_smoothed")]["sample_centered_rho"].iloc[0])
    actual_res = float(result_df[(result_df["class"] == "actual_dm1") & (result_df["mode"] == "coord_qc_residual_target")]["sample_centered_rho"].iloc[0])

    loo_raw = result_df[(result_df["class"] == "panel_loo") & (result_df["mode"] == "raw_smoothed")]["sample_centered_rho"]
    loo_res = result_df[(result_df["class"] == "panel_loo") & (result_df["mode"] == "coord_qc_residual_target")]["sample_centered_rho"]

    summary = {
        "n_spots": int(len(merged)),
        "n_samples": int(pd.Series(samples).nunique()),
        "n_random_modules": int(N_RANDOM_MODULES),
        "random_module_size": int(MODULE_SIZE),
        "actual_raw_rho": actual_raw,
        "actual_raw_percentile": percentile(actual_raw, random_raw),
        "actual_raw_empirical_p": empirical_upper_p(actual_raw, random_raw),
        "actual_residual_rho": actual_res,
        "actual_residual_percentile": percentile(actual_res, random_res),
        "actual_residual_empirical_p": empirical_upper_p(actual_res, random_res),
        "random_raw_mean": float(np.nanmean(random_raw)),
        "random_raw_p95": float(np.nanquantile(random_raw, 0.95)),
        "random_residual_mean": float(np.nanmean(random_res)),
        "random_residual_p95": float(np.nanquantile(random_res, 0.95)),
        "loo_raw_min": float(loo_raw.min()),
        "loo_raw_max": float(loo_raw.max()),
        "loo_residual_min": float(loo_res.min()),
        "loo_residual_max": float(loo_res.max()),
    }
    (OUT / "GSE230424_RANDOM_MODULE_SPECIFICITY_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")

    loo_df = result_df[result_df["class"] == "panel_loo"].copy()
    write_figures(result_df, loo_df, summary)
    write_report(result_df, module_df, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
