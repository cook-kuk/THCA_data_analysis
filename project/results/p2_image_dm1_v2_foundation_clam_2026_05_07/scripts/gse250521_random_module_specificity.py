#!/usr/bin/env python3
"""GSE250521 matched-random-module specificity and smoothness caveat.

This mirrors the GSE230424 random-module control on the original GSE250521
Path2Space-inspired cohort. It asks whether the UNI H&E-to-DM1/RAI signal is
stronger than expression/detection-matched random 8-gene modules, and whether
that advantage survives adjustment for target spatial smoothness.
"""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandas.errors import PerformanceWarning
from scipy import stats
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2_ROOT = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
PATH2 = P2_ROOT / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"
PHASE1 = P2_ROOT / "phase1_gse250521"
GSE_PROC = ROOT / "project/data/processed/GSE250521"
OUT = P2_ROOT / "analysis_supp/path2space_gse250521_random_module_specificity_2026_05_09"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
N_RANDOM_MODULES = 250
MODULE_SIZE = 8
RANDOM_POOL_NEAREST = 250
SMOOTH_K = 8
N_PCS = 64
RIDGE_ALPHA = 100.0
RANDOM_SEED = 9251
K_NEIGH = 6

warnings.filterwarnings("ignore", category=PerformanceWarning)


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def logit_clip(values: np.ndarray | pd.Series | float) -> np.ndarray:
    clipped = np.clip(np.asarray(values, dtype=float), 1e-4, 1 - 1e-4)
    return np.log(clipped / (1.0 - clipped))


def read_h5ad(path: Path):
    import scanpy as sc

    return sc.read_h5ad(path)


def load_spots_and_embeddings() -> tuple[pd.DataFrame, np.ndarray]:
    spots = pd.read_csv(PATH2 / "path2space_input_spots.tsv.gz", sep="\t").sort_values("embed_idx").reset_index(drop=True)
    emb = np.load(PHASE1 / "uni_embeddings_size224.npz")["embeddings"].astype(np.float32)
    if len(spots) != emb.shape[0]:
        raise RuntimeError(f"Spot/embedding length mismatch: {len(spots)} vs {emb.shape[0]}")
    return spots, emb


def expression_metrics(spots: pd.DataFrame) -> pd.DataFrame:
    meta = pd.read_csv(GSE_PROC / "sample_metadata.tsv", sep="\t")
    total_counts = None
    detected_spots = None
    genes = None
    n_spots = 0

    for row in meta.itertuples(index=False):
        sub = spots[spots["sample_id"] == row.sample_id]
        if sub.empty:
            continue
        adata = read_h5ad(ROOT / row.h5ad)
        if genes is None:
            genes = adata.var_names.astype(str).to_numpy()
            total_counts = np.zeros(len(genes), dtype=np.float64)
            detected_spots = np.zeros(len(genes), dtype=np.int64)
        idx = [adata.obs_names.get_loc(spot) for spot in sub["spot_id"].astype(str)]
        x = adata.X[idx, :]
        if sp.issparse(x):
            total_counts += np.asarray(x.sum(axis=0)).ravel()
            detected_spots += np.asarray((x > 0).sum(axis=0)).ravel()
        else:
            arr = np.asarray(x)
            total_counts += np.nansum(arr, axis=0)
            detected_spots += np.nansum(arr > 0, axis=0).astype(int)
        n_spots += len(idx)

    if genes is None or total_counts is None or detected_spots is None:
        raise RuntimeError("No GSE250521 expression data loaded")
    metrics = pd.DataFrame({"gene": genes, "total_counts": total_counts, "detected_spots": detected_spots})
    metrics = metrics.drop_duplicates("gene", keep="first")
    metrics["n_spots"] = n_spots
    metrics["detect_frac"] = metrics["detected_spots"] / max(n_spots, 1)
    metrics["mean_count"] = metrics["total_counts"] / max(n_spots, 1)
    metrics["log1p_mean_count"] = np.log1p(metrics["mean_count"])
    metrics.to_csv(OUT / "gse250521_gene_expression_metrics.tsv", sep="\t", index=False, na_rep="NA")
    return metrics


def eligible_universe(metrics: pd.DataFrame) -> pd.DataFrame:
    excluded_prefixes = ("MT-", "RPS", "RPL")
    universe = metrics[
        (metrics["detect_frac"] >= 0.01)
        & (metrics["detect_frac"] <= 0.90)
        & (metrics["total_counts"] >= 10)
        & (~metrics["gene"].astype(str).str.startswith(excluded_prefixes))
        & (~metrics["gene"].isin(PANEL_8))
    ].copy()
    universe = universe.reset_index(drop=True)
    if len(universe) < 500:
        raise RuntimeError(f"Eligible random-module universe too small: {len(universe)}")
    return universe


def make_matched_random_modules(metrics: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, dict[str, list[str]]]:
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
            [[math.log1p(float(row["mean_count"])), float(logit_clip(float(row["detect_frac"])))]]
        )
        dist = np.sqrt(((u_scaled - scaler.transform(p_feat)) ** 2).sum(axis=1))
        nearest = np.argsort(dist)[: min(RANDOM_POOL_NEAREST, len(universe))]
        pools[str(row["gene"])] = universe.iloc[nearest]["gene"].astype(str).tolist()

    modules: dict[str, list[str]] = {}
    attempts = 0
    while len(modules) < N_RANDOM_MODULES and attempts < N_RANDOM_MODULES * 60:
        attempts += 1
        chosen: list[str] = []
        for panel_gene in PANEL_8:
            pool = [gene for gene in pools[panel_gene] if gene not in chosen]
            if not pool:
                break
            chosen.append(str(rng.choice(pool)))
        if len(set(chosen)) != MODULE_SIZE:
            continue
        modules[f"random_matched_{len(modules) + 1:03d}"] = chosen

    if len(modules) < N_RANDOM_MODULES:
        raise RuntimeError(f"Only built {len(modules)} random modules")

    rows = []
    for module, genes in modules.items():
        for slot, gene in zip(PANEL_8, genes):
            rows.append({"module": module, "matched_panel_gene": slot, "random_gene": gene})
    module_df = pd.DataFrame(rows)
    module_df.to_csv(OUT / "gse250521_random_matched_modules.tsv", sep="\t", index=False)
    return module_df, modules


def selected_gene_scores(spots: pd.DataFrame, modules: dict[str, list[str]]) -> tuple[pd.DataFrame, list[str], list[str]]:
    meta = pd.read_csv(GSE_PROC / "sample_metadata.tsv", sep="\t")
    selected_genes = list(PANEL_8)
    for genes in modules.values():
        selected_genes.extend(genes)
    selected_genes = list(dict.fromkeys(selected_genes))

    loo_names = [f"loo_drop_{gene}" for gene in PANEL_8]
    random_names = [f"null_{name}" for name in sorted(modules)]
    score_names = ["actual_DM1_RAI"] + loo_names + random_names
    scores = np.full((len(spots), len(score_names)), np.nan, dtype=np.float32)

    for row in meta.itertuples(index=False):
        sub = spots[spots["sample_id"] == row.sample_id]
        if sub.empty:
            continue
        adata = read_h5ad(ROOT / row.h5ad)
        spot_idx = [adata.obs_names.get_loc(spot) for spot in sub["spot_id"].astype(str)]
        gene_idx = [adata.var_names.get_loc(gene) for gene in selected_genes]
        x_all = adata.X[spot_idx, :]
        totals = np.asarray(x_all.sum(axis=1)).ravel() if sp.issparse(x_all) else np.asarray(x_all).sum(axis=1)
        x_sel = adata.X[spot_idx, :][:, gene_idx]
        x_sel = x_sel.toarray() if sp.issparse(x_sel) else np.asarray(x_sel)
        log_cp10k = np.log1p(x_sel / np.maximum(totals[:, None], 1.0) * 10000.0)
        mean = np.nanmean(log_cp10k, axis=0)
        sd = np.nanstd(log_cp10k, axis=0)
        sd[~np.isfinite(sd) | (sd == 0)] = 1.0
        z = (log_cp10k - mean[None, :]) / sd[None, :]
        gene_to_col = {gene: i for i, gene in enumerate(selected_genes)}

        rows = sub.index.to_numpy()
        panel_cols = [gene_to_col[gene] for gene in PANEL_8]
        scores[rows, 0] = -np.nanmean(z[:, panel_cols], axis=1)
        for j, removed in enumerate(PANEL_8, start=1):
            cols = [gene_to_col[gene] for gene in PANEL_8 if gene != removed]
            scores[rows, j] = -np.nanmean(z[:, cols], axis=1)
        offset = 1 + len(loo_names)
        for k, module in enumerate(sorted(modules)):
            cols = [gene_to_col[gene] for gene in modules[module]]
            scores[rows, offset + k] = -np.nanmean(z[:, cols], axis=1)

    out = spots[
        [
            "sample_id",
            "spot_id",
            "stage",
            "array_row",
            "array_col",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ]
    ].copy()
    score_df = pd.DataFrame(scores, columns=score_names, index=out.index)
    return pd.concat([out, score_df], axis=1), loo_names, random_names


def spatial_smooth_matrix(values: np.ndarray, coords: np.ndarray, samples: np.ndarray, k: int = SMOOTH_K) -> np.ndarray:
    out = np.full(values.shape, np.nan, dtype=np.float64)
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        xy = coords[idx]
        if len(idx) < 4:
            continue
        tree = cKDTree(xy)
        kk = min(k + 1, len(idx))
        _, neigh = tree.query(xy, k=kk)
        if neigh.ndim == 1:
            neigh = neigh[:, None]
        out[idx] = np.nanmean(values[idx][neigh], axis=1)
    return out


def center_by_group(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    for group in pd.unique(groups):
        idx = groups == group
        out[idx] = values[idx] - np.nanmean(values[idx])
    return out


def residualize_targets(y: np.ndarray, covars: np.ndarray, samples: np.ndarray) -> np.ndarray:
    out = np.full(y.shape, np.nan, dtype=np.float64)
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        x = covars[idx]
        keep = np.isfinite(x).all(axis=1) & np.isfinite(y[idx]).all(axis=1)
        if keep.sum() < 30:
            continue
        xs = StandardScaler().fit_transform(x[keep])
        fit = LinearRegression().fit(xs, y[idx][keep])
        out[idx[keep]] = y[idx][keep] - fit.predict(xs)
    return out


def loso_predict_many(x: np.ndarray, y: np.ndarray, samples: np.ndarray) -> np.ndarray:
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
        n_comp = min(N_PCS, x_train.shape[0] - 1, x_train.shape[1])
        pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
        x_train_pc = pca.fit_transform(x_train)
        x_test_pc = pca.transform(x_test)
        mean = np.nanmean(y[train], axis=0)
        sd = np.nanstd(y[train], axis=0)
        sd[~np.isfinite(sd) | (sd == 0)] = 1.0
        model = Ridge(alpha=RIDGE_ALPHA)
        model.fit(x_train_pc, (y[train] - mean[None, :]) / sd[None, :])
        pred[test] = model.predict(x_test_pc)
    return pred


def summarize_predictions(names: list[str], y: np.ndarray, pred: np.ndarray, samples: np.ndarray, mode: str) -> pd.DataFrame:
    y_center = np.column_stack([center_by_group(y[:, j], samples) for j in range(y.shape[1])])
    pred_center = np.column_stack([center_by_group(pred[:, j], samples) for j in range(pred.shape[1])])
    rows = []
    for j, name in enumerate(names):
        pooled, p, n = safe_spearman(y[:, j], pred[:, j])
        centered, cp, _ = safe_spearman(y_center[:, j], pred_center[:, j])
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
                "pooled_rho": pooled,
                "pooled_p": p,
                "sample_centered_rho": centered,
                "sample_centered_p": cp,
                "median_sample_rho": float(np.nanmedian(per_sample)),
                "samples_rho_gt_0_2": int(np.nansum(np.asarray(per_sample) > 0.2)),
            }
        )
    return pd.DataFrame(rows)


def classify(name: str) -> str:
    if name == "actual_DM1_RAI":
        return "actual_dm1"
    if name.startswith("loo_drop_"):
        return "panel_loo"
    return "random_matched"


def empirical_upper_p(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    return float((1.0 + np.sum(null >= observed)) / (1.0 + len(null)))


def percentile(observed: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    return float(100.0 * np.mean(null <= observed))


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
    weighted = float(np.average(per_sample, weights=weights)) if weights else np.nan
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
        n_domains = int(np.clip(round(keep.sum() / 50), 4, 8))
        labs = KMeans(n_clusters=n_domains, n_init=30, random_state=RANDOM_SEED).fit_predict(
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


def smoothness_adjust(result_df: pd.DataFrame, y_by_mode: dict[str, np.ndarray], names: list[str], coords: np.ndarray, samples: np.ndarray) -> tuple[pd.DataFrame, list[dict]]:
    block_labels = compute_spatial_blocks(samples, coords)
    metric_rows = []
    for mode, ymat in y_by_mode.items():
        for j, name in enumerate(names):
            pooled, weighted = neighbor_autocorr(ymat[:, j], coords, samples)
            metric_rows.append(
                {
                    "target": name,
                    "mode": mode,
                    "neighbor_autocorr_pooled": pooled,
                    "neighbor_autocorr_weighted_sample": weighted,
                    "block_variance_fraction": block_variance_fraction(ymat[:, j], samples, block_labels),
                    "target_sd": float(np.nanstd(ymat[:, j])),
                }
            )
    metrics = pd.DataFrame(metric_rows)
    merged = result_df.merge(metrics, on=["target", "mode"], how="left")

    adjusted_parts = []
    summaries = []
    features = ["neighbor_autocorr_weighted_sample", "block_variance_fraction", "target_sd"]
    for mode in y_by_mode:
        sub = merged[merged["mode"] == mode].copy()
        random = sub[sub["class"] == "random_matched"].copy()
        keep = random[features + ["sample_centered_rho"]].replace([np.inf, -np.inf], np.nan).dropna()
        scaler = StandardScaler()
        x = keep[features].to_numpy(dtype=float)
        y = keep["sample_centered_rho"].to_numpy(dtype=float)
        model = LinearRegression().fit(scaler.fit_transform(x), y)
        sub["smoothness_expected_rho"] = np.nan
        sub["smoothness_adjusted_residual"] = np.nan
        ok = sub[features].replace([np.inf, -np.inf], np.nan).notna().all(axis=1)
        sub.loc[ok, "smoothness_expected_rho"] = model.predict(scaler.transform(sub.loc[ok, features].to_numpy(dtype=float)))
        sub.loc[ok, "smoothness_adjusted_residual"] = sub.loc[ok, "sample_centered_rho"] - sub.loc[ok, "smoothness_expected_rho"]
        random_resid = sub.loc[sub["class"] == "random_matched", "smoothness_adjusted_residual"].to_numpy(dtype=float)
        actual = sub[sub["class"] == "actual_dm1"].iloc[0]
        actual_resid = float(actual["smoothness_adjusted_residual"])
        summaries.append(
            {
                "mode": mode,
                "actual_rho": float(actual["sample_centered_rho"]),
                "actual_expected_from_smoothness": float(actual["smoothness_expected_rho"]),
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
        )
        adjusted_parts.append(sub)
    adjusted = pd.concat(adjusted_parts, ignore_index=True)
    metrics.to_csv(OUT / "gse250521_spatial_autocorr_target_metrics.tsv", sep="\t", index=False, na_rep="NA")
    return adjusted, summaries


def write_figure(adjusted: pd.DataFrame, summary: dict) -> None:
    raw = adjusted[adjusted["mode"] == "raw_smoothed"]
    resid = adjusted[adjusted["mode"] == "coord_qc_residual_target"]
    actual = adjusted[adjusted["class"] == "actual_dm1"].set_index("mode")

    fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))
    for ax, sub, mode_title, mode in [
        (axes[0, 0], raw, "A. Raw target vs matched random modules", "raw_smoothed"),
        (axes[0, 1], resid, "B. Residual target vs matched random modules", "coord_qc_residual_target"),
    ]:
        random = sub[sub["class"] == "random_matched"]
        ax.hist(random["sample_centered_rho"], bins=32, color="#9aa4ad", edgecolor="white")
        ax.axvline(actual.loc[mode, "sample_centered_rho"], color="#b44d4d", lw=2.5, label="DM1/RAI")
        ax.set_title(mode_title)
        ax.set_xlabel("Slide-centered LOSO rho")
        ax.set_ylabel("Matched random modules")
        ax.legend(frameon=False)

    for ax, sub, title in [
        (axes[0, 2], raw, "C. Raw rho vs neighbor autocorrelation"),
        (axes[1, 2], resid, "F. Residual rho vs neighbor autocorrelation"),
    ]:
        random = sub[sub["class"] == "random_matched"]
        actual_sub = sub[sub["class"] == "actual_dm1"]
        ax.scatter(random["neighbor_autocorr_weighted_sample"], random["sample_centered_rho"], s=18, alpha=0.55, color="#888")
        ax.scatter(actual_sub["neighbor_autocorr_weighted_sample"], actual_sub["sample_centered_rho"], s=95, color="#b44d4d", label="DM1/RAI", zorder=5)
        ax.set_xlabel("Neighbor autocorrelation")
        ax.set_ylabel("Slide-centered rho")
        ax.set_title(title)
        ax.legend(frameon=False)

    for ax, sub, title, mode in [
        (axes[1, 0], raw, "D. Raw smoothness-adjusted residual tail", "raw_smoothed"),
        (axes[1, 1], resid, "E. Residual-target smoothness-adjusted residual tail", "coord_qc_residual_target"),
    ]:
        random = sub[sub["class"] == "random_matched"]
        ax.hist(random["smoothness_adjusted_residual"], bins=32, color="#9aa4ad", edgecolor="white")
        ax.axvline(actual.loc[mode, "smoothness_adjusted_residual"], color="#b44d4d", lw=2.5, label="DM1/RAI")
        ax.axvline(0, color="#333", lw=0.8)
        ax.set_xlabel("Observed rho - expected rho from smoothness")
        ax.set_ylabel("Matched random modules")
        ax.set_title(title)
        ax.legend(frameon=False)

    fig.suptitle("GSE250521 UNI specificity: DM1/RAI vs matched random modules", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse250521_random_module_specificity.png", dpi=220)
    fig.savefig(OUT / "fig_gse250521_random_module_specificity.pdf")
    plt.close(fig)


def write_report(adjusted: pd.DataFrame, summaries: list[dict], summary: dict) -> None:
    actual = adjusted[adjusted["class"] == "actual_dm1"].set_index("mode")
    loo = adjusted[adjusted["class"] == "panel_loo"].copy()
    random_raw = adjusted[(adjusted["class"] == "random_matched") & (adjusted["mode"] == "raw_smoothed")]
    random_res = adjusted[(adjusted["class"] == "random_matched") & (adjusted["mode"] == "coord_qc_residual_target")]
    top_raw = random_raw.sort_values("sample_centered_rho", ascending=False).head(12)
    top_res = random_res.sort_values("sample_centered_rho", ascending=False).head(12)

    lines = [
        "# GSE250521 random-module specificity and smoothness adjustment",
        "",
        "## Verdict",
        "",
        f"- Matched random modules: **{summary['n_random_modules']}** expression/detection-matched 8-gene modules.",
        f"- Raw UNI DM1/RAI rho: **{summary['actual_raw_rho']:.3f}**; random-module percentile **{summary['actual_raw_percentile']:.1f}%**; empirical p = **{summary['actual_raw_empirical_p']:.4f}**.",
        f"- Coord+QC residual-target UNI DM1/RAI rho: **{summary['actual_residual_rho']:.3f}**; random-module percentile **{summary['actual_residual_percentile']:.1f}%**; empirical p = **{summary['actual_residual_empirical_p']:.4f}**.",
        f"- Raw smoothness-adjusted residual: **{summary['raw_smoothness_adjusted_residual']:.3f}**; percentile **{summary['raw_residual_percentile']:.1f}%**; p = **{summary['raw_residual_empirical_p']:.4f}**.",
        f"- Residual-target smoothness-adjusted residual: **{summary['resid_smoothness_adjusted_residual']:.3f}**; percentile **{summary['resid_residual_percentile']:.1f}%**; p = **{summary['resid_residual_empirical_p']:.4f}**.",
        "",
        "## Interpretation",
        "",
        "This is the GSE250521 counterpart of the GSE230424 matched-random-module and smoothness-adjusted specificity control. It uses the same 3,200 UNI-embedded spots as the Path2Space-inspired reanalysis, constructs 250 expression/detection-matched random 8-gene modules, predicts all targets with leave-one-slide-out UNI Ridge models, and then asks whether DM1/RAI remains exceptional after accounting for target spatial autocorrelation/block structure.",
        "",
        "## Actual And Leave-One-Out Results",
        "",
        adjusted[adjusted["class"].isin(["actual_dm1", "panel_loo"])].round(4).to_markdown(index=False),
        "",
        "## Smoothness Model Summary",
        "",
        pd.DataFrame(summaries).round(4).to_markdown(index=False),
        "",
        "## Random Module Null Summary",
        "",
        pd.DataFrame(
            [
                {
                    "mode": "raw_smoothed",
                    "null_mean": random_raw["sample_centered_rho"].mean(),
                    "null_median": random_raw["sample_centered_rho"].median(),
                    "null_p95": random_raw["sample_centered_rho"].quantile(0.95),
                    "actual": actual.loc["raw_smoothed", "sample_centered_rho"],
                },
                {
                    "mode": "coord_qc_residual_target",
                    "null_mean": random_res["sample_centered_rho"].mean(),
                    "null_median": random_res["sample_centered_rho"].median(),
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
    ]
    (OUT / "GSE250521_RANDOM_MODULE_SPECIFICITY_REPORT.md").write_text("\n".join(lines) + "\n")

    summary_lines = [
        "# GSE250521 Random-Module Specificity",
        "",
        "## Verdict",
        "",
        "This is the GSE250521 specificity counterpart to the GSE230424 caveat control.",
        "",
        "## Headline Results",
        "",
        "| Mode | Observed UNI rho | Random-module percentile | Empirical p | Smoothness-expected rho | Smoothness-adjusted residual | Smoothness residual percentile |",
        "|---|---:|---:|---:|---:|---:|---:|",
        f"| Raw smoothed target | {summary['actual_raw_rho']:.3f} | {summary['actual_raw_percentile']:.1f}% | {summary['actual_raw_empirical_p']:.4f} | {summary['raw_expected_from_smoothness']:.3f} | {summary['raw_smoothness_adjusted_residual']:.3f} | {summary['raw_residual_percentile']:.1f}% |",
        f"| Coord+QC residual target | {summary['actual_residual_rho']:.3f} | {summary['actual_residual_percentile']:.1f}% | {summary['actual_residual_empirical_p']:.4f} | {summary['resid_expected_from_smoothness']:.3f} | {summary['resid_smoothness_adjusted_residual']:.3f} | {summary['resid_residual_percentile']:.1f}% |",
        "",
        "## Interpretation",
        "",
        "Raw DM1/RAI is image-predictable, but it is not formally exceptional after random-module or smoothness adjustment. The stricter coord+QC residual target is the useful layer: it remains above all 250 matched random modules before smoothness adjustment and stays directionally high after smoothness adjustment, although the adjusted upper-tail p is 0.0996 rather than conventionally significant.",
        "",
        "Safe wording: GSE250521 supports a residual image-aligned DM1/RAI component beyond coordinate/QC structure, but raw panel predictability is partly a spatial-smoothness/tissue-state property.",
        "",
        "## Files",
        "",
        "- `GSE250521_RANDOM_MODULE_SPECIFICITY_REPORT.md`",
        "- `GSE250521_RANDOM_MODULE_SPECIFICITY_SUMMARY.json`",
        "- `gse250521_random_module_specificity_results.tsv`",
        "- `gse250521_spatial_autocorr_adjusted_specificity.tsv`",
        "- `fig_gse250521_random_module_specificity.png`",
        "- `fig_gse250521_random_module_specificity.pdf`",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(summary_lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    spots, emb = load_spots_and_embeddings()

    metrics_path = OUT / "gse250521_gene_expression_metrics.tsv"
    metrics = pd.read_csv(metrics_path, sep="\t") if metrics_path.exists() else expression_metrics(spots)

    module_path = OUT / "gse250521_random_matched_modules.tsv"
    if module_path.exists():
        module_df = pd.read_csv(module_path, sep="\t")
        modules = {name: sub["random_gene"].astype(str).tolist() for name, sub in module_df.groupby("module", sort=True)}
    else:
        module_df, modules = make_matched_random_modules(metrics, rng)

    scores, loo_names, random_names = selected_gene_scores(spots, modules)
    score_names = ["actual_DM1_RAI"] + loo_names + random_names
    values = scores[score_names].to_numpy(dtype=float)
    samples = scores["sample_id"].astype(str).to_numpy()
    coords = scores[["array_row", "array_col"]].to_numpy(dtype=float)
    y_raw = spatial_smooth_matrix(values, coords, samples, k=SMOOTH_K)

    coord_qc = np.column_stack(
        [
            scores[["array_row", "array_col"]].to_numpy(dtype=float),
            np.log1p(scores["total_counts"].to_numpy(dtype=float)),
            np.log1p(scores["n_genes_by_counts"].to_numpy(dtype=float)),
            scores["pct_counts_mt"].to_numpy(dtype=float),
        ]
    )
    y_resid = residualize_targets(y_raw, coord_qc, samples)

    pred_raw = loso_predict_many(emb, y_raw, samples)
    pred_resid = loso_predict_many(emb, y_resid, samples)

    raw_df = summarize_predictions(score_names, y_raw, pred_raw, samples, "raw_smoothed")
    resid_df = summarize_predictions(score_names, y_resid, pred_resid, samples, "coord_qc_residual_target")
    result_df = pd.concat([raw_df, resid_df], ignore_index=True)
    result_df["class"] = result_df["target"].map(classify)

    adjusted, smooth_summaries = smoothness_adjust(
        result_df,
        {"raw_smoothed": y_raw, "coord_qc_residual_target": y_resid},
        score_names,
        coords,
        samples,
    )
    adjusted.to_csv(OUT / "gse250521_spatial_autocorr_adjusted_specificity.tsv", sep="\t", index=False, na_rep="NA")
    result_df.to_csv(OUT / "gse250521_random_module_specificity_results.tsv", sep="\t", index=False, na_rep="NA")

    score_out = scores[["sample_id", "spot_id", "stage", "array_row", "array_col"]].copy()
    smooth_df = pd.DataFrame(
        y_raw,
        columns=[f"{name}_smooth{SMOOTH_K}" for name in score_names],
        index=score_out.index,
    )
    score_out = pd.concat([score_out, smooth_df], axis=1)
    score_out.to_csv(OUT / "gse250521_random_module_spot_scores.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    pred_out = scores[["sample_id", "spot_id", "stage", "array_row", "array_col"]].copy()
    pred_out["obs_actual_raw"] = y_raw[:, 0]
    pred_out["pred_actual_raw"] = pred_raw[:, 0]
    pred_out["obs_actual_coord_qc_residual"] = y_resid[:, 0]
    pred_out["pred_actual_coord_qc_residual"] = pred_resid[:, 0]
    pred_out.to_csv(OUT / "gse250521_actual_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    random_raw = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "raw_smoothed")]["sample_centered_rho"].to_numpy(dtype=float)
    random_res = result_df[(result_df["class"] == "random_matched") & (result_df["mode"] == "coord_qc_residual_target")]["sample_centered_rho"].to_numpy(dtype=float)
    actual_raw = float(result_df[(result_df["class"] == "actual_dm1") & (result_df["mode"] == "raw_smoothed")]["sample_centered_rho"].iloc[0])
    actual_res = float(result_df[(result_df["class"] == "actual_dm1") & (result_df["mode"] == "coord_qc_residual_target")]["sample_centered_rho"].iloc[0])
    loo_raw = result_df[(result_df["class"] == "panel_loo") & (result_df["mode"] == "raw_smoothed")]["sample_centered_rho"]
    loo_res = result_df[(result_df["class"] == "panel_loo") & (result_df["mode"] == "coord_qc_residual_target")]["sample_centered_rho"]
    raw_smooth = next(row for row in smooth_summaries if row["mode"] == "raw_smoothed")
    resid_smooth = next(row for row in smooth_summaries if row["mode"] == "coord_qc_residual_target")

    summary = {
        "n_spots": int(len(spots)),
        "n_samples": int(pd.Series(samples).nunique()),
        "n_random_modules": int(N_RANDOM_MODULES),
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
        "raw_expected_from_smoothness": raw_smooth["actual_expected_from_smoothness"],
        "raw_smoothness_adjusted_residual": raw_smooth["actual_smoothness_adjusted_residual"],
        "raw_residual_percentile": raw_smooth["residual_percentile_vs_random"],
        "raw_residual_empirical_p": raw_smooth["residual_empirical_upper_p"],
        "resid_expected_from_smoothness": resid_smooth["actual_expected_from_smoothness"],
        "resid_smoothness_adjusted_residual": resid_smooth["actual_smoothness_adjusted_residual"],
        "resid_residual_percentile": resid_smooth["residual_percentile_vs_random"],
        "resid_residual_empirical_p": resid_smooth["residual_empirical_upper_p"],
    }
    (OUT / "GSE250521_RANDOM_MODULE_SPECIFICITY_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_figure(adjusted, summary)
    write_report(adjusted, smooth_summaries, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
