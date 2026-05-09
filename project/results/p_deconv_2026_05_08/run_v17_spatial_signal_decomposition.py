#!/usr/bin/env python3
"""v17: deeper decomposition of the GSE250521 spatial MAPK-panel failure.

v16 showed that the raw positive Visium MAPK x Panel-8 correlation attenuates
toward null after strict adjustment. v17 asks why: covariate family, rank-based
partial Spearman, detection/dropout, MAPK submodules, gene-pair consistency,
and a random-module specificity null.
"""
from __future__ import annotations

import html
import json
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from sklearn.linear_model import LinearRegression

warnings.filterwarnings("ignore")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROC = ROOT / "project" / "data" / "processed" / "GSE250521"
OUT = ROOT / "project" / "results" / "p_deconv_2026_05_08"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSETS = HUB / "assets" / "paper1"
LIVE_ASSETS = LIVE / "assets" / "paper1"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
MAPK_OUT = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]
MAPK_MODULES = {
    "MAPK_all9": MAPK_OUT,
    "feedback_DUSP_SPRY": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4"],
    "transcriptional_ETV_PHLDA1": ["ETV4", "ETV5", "PHLDA1"],
    "no_CCND1": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1"],
    "CCND1_only": ["CCND1"],
}
ADJUSTMENT_ORDER = [
    "raw",
    "epithelial",
    "qc",
    "microenv",
    "spatial",
    "detection",
    "epithelial_qc",
    "epithelial_qc_microenv",
    "epithelial_qc_microenv_spatial",
    "epithelial_qc_microenv_spatial_detection",
]
ADJUSTMENT_LABELS = {
    "raw": "Raw",
    "epithelial": "Epi",
    "qc": "QC",
    "microenv": "Microenv",
    "spatial": "Spatial",
    "detection": "Detect",
    "epithelial_qc": "Epi+QC",
    "epithelial_qc_microenv": "Full",
    "epithelial_qc_microenv_spatial": "Full+spatial",
    "epithelial_qc_microenv_spatial_detection": "Full+detect",
}
RANDOM_N = 50


def get_expr(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    col = adata.X[:, adata.var_names.get_loc(gene)]
    if sp.issparse(col):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return col.astype(float)


def module_z_and_detect(adata, genes: list[str]) -> tuple[np.ndarray, np.ndarray, int]:
    vals = []
    detected = np.zeros(adata.n_obs, dtype=float)
    n = 0
    for gene in genes:
        v = get_expr(adata, gene)
        if np.all(~np.isfinite(v)):
            continue
        n += 1
        detected += (v > 0).astype(float)
        sd = np.nanstd(v)
        vals.append((v - np.nanmean(v)) / (sd + 1e-9) if np.isfinite(sd) and sd > 0 else np.zeros_like(v))
    if not vals:
        return np.full(adata.n_obs, np.nan), np.full(adata.n_obs, np.nan), 0
    return np.nanmean(np.vstack(vals), axis=0), detected, n


def spearman_block(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 20:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def pearson_block(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 20:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.pearsonr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def standardize_columns(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    cols = []
    for j in range(x.shape[1]):
        col = x[:, j]
        sd = np.nanstd(col)
        if np.isfinite(sd) and sd > 1e-9:
            cols.append((col - np.nanmean(col)) / (sd + 1e-9))
    if not cols:
        return np.empty((x.shape[0], 0))
    return np.column_stack(cols)


def residualize(y: np.ndarray, covars: np.ndarray | None) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    out = np.full_like(y, np.nan, dtype=float)
    if covars is None:
        out[np.isfinite(y)] = y[np.isfinite(y)]
        return out
    covars = np.asarray(covars, dtype=float)
    keep = np.isfinite(y) & np.isfinite(covars).all(axis=1)
    if keep.sum() < 20:
        return out
    xs = standardize_columns(covars[keep])
    if xs.shape[1] == 0:
        out[keep] = y[keep] - np.nanmean(y[keep])
        return out
    out[keep] = y[keep] - LinearRegression().fit(xs, y[keep]).predict(xs)
    return out


def partial_spearman_residuals(y: np.ndarray, covars: np.ndarray | None) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    out = np.full_like(y, np.nan, dtype=float)
    if covars is None:
        keep = np.isfinite(y)
        out[keep] = stats.rankdata(y[keep])
        return out
    covars = np.asarray(covars, dtype=float)
    keep = np.isfinite(y) & np.isfinite(covars).all(axis=1)
    if keep.sum() < 20:
        return out
    y_rank = stats.rankdata(y[keep])
    cov_rank = np.column_stack([stats.rankdata(covars[keep, j]) for j in range(covars.shape[1])])
    xs = standardize_columns(cov_rank)
    if xs.shape[1] == 0:
        out[keep] = y_rank - np.nanmean(y_rank)
        return out
    out[keep] = y_rank - LinearRegression().fit(xs, y_rank).predict(xs)
    return out


def lin_r2(y: np.ndarray, covars: np.ndarray) -> tuple[float, int, int]:
    y = np.asarray(y, dtype=float)
    covars = np.asarray(covars, dtype=float)
    keep = np.isfinite(y) & np.isfinite(covars).all(axis=1)
    if keep.sum() < 20:
        return np.nan, int(keep.sum()), 0
    xs = standardize_columns(covars[keep])
    if xs.shape[1] == 0:
        return np.nan, int(keep.sum()), 0
    model = LinearRegression().fit(xs, y[keep])
    return float(model.score(xs, y[keep])), int(keep.sum()), int(xs.shape[1])


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.isfinite(p)
    if finite.sum() == 0:
        return pd.Series(out, index=pvals.index)
    order = np.argsort(p[finite])
    ranked = p[finite][order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.full(n, np.nan)
    tmp[order] = np.minimum(q, 1.0)
    out[np.where(finite)[0]] = tmp
    return pd.Series(out, index=pvals.index)


def covariate_sets(obs: pd.DataFrame, coords: np.ndarray, mapk_detect: np.ndarray, panel_detect: np.ndarray) -> dict[str, np.ndarray | None]:
    epithelial = obs["Epithelial_score"].astype(float).to_numpy()
    qc = np.column_stack(
        [
            np.log1p(obs["total_counts"].astype(float).to_numpy()),
            np.log1p(obs["n_genes_by_counts"].astype(float).to_numpy()),
            obs["pct_counts_mt"].astype(float).to_numpy(),
        ]
    )
    microenv = np.column_stack(
        [
            obs["CAF_ECM_score"].astype(float).to_numpy(),
            obs["EMT_score"].astype(float).to_numpy(),
            obs["Hypoxia_score"].astype(float).to_numpy(),
            obs["Proliferation_score"].astype(float).to_numpy(),
        ]
    )
    row = coords[:, 0]
    col = coords[:, 1]
    spatial = np.column_stack([row, col, row * row, col * col, row * col])
    detection = np.column_stack([mapk_detect, panel_detect])
    return {
        "raw": None,
        "epithelial": epithelial.reshape(-1, 1),
        "qc": qc,
        "microenv": microenv,
        "spatial": spatial,
        "detection": detection,
        "epithelial_qc": np.column_stack([epithelial, qc]),
        "epithelial_qc_microenv": np.column_stack([epithelial, qc, microenv]),
        "epithelial_qc_microenv_spatial": np.column_stack([epithelial, qc, microenv, spatial]),
        "epithelial_qc_microenv_spatial_detection": np.column_stack([epithelial, qc, microenv, spatial, detection]),
    }


def load_scored_h5ad(row: pd.Series):
    import scanpy as sc

    raw = ROOT / str(row["h5ad"])
    scored = raw.parent / (raw.stem.replace(".raw", "") + ".scored.h5ad")
    use_h5 = scored if scored.exists() else raw
    adata = sc.read_h5ad(use_h5)
    if "log1p" not in (adata.uns or {}):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    return adata


def eligible_random_genes(adata, forbidden: set[str]) -> list[str]:
    var = adata.var.copy()
    names = pd.Index(adata.var_names.astype(str))
    mt = var["mt"].astype(bool).to_numpy() if "mt" in var else names.str.upper().str.startswith("MT-").to_numpy()
    if "n_cells_by_counts" in var:
        det = pd.to_numeric(var["n_cells_by_counts"], errors="coerce").fillna(0).to_numpy()
        min_cells = max(20, int(0.01 * adata.n_obs))
        max_cells = int(0.95 * adata.n_obs)
        ok = (det >= min_cells) & (det <= max_cells)
    else:
        ok = np.ones(len(names), dtype=bool)
    special = names.isin(forbidden)
    genes = names[ok & ~mt & ~special].tolist()
    if len(genes) < 200:
        genes = names[~mt & ~special].tolist()
    return genes


def run_v17() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    meta = pd.read_csv(PROC / "sample_metadata.tsv", sep="\t")
    component_rows = []
    r2_rows = []
    qc_rows = []
    submodule_rows = []
    gene_rows = []
    random_rows = []
    pooled_chunks = []
    rng = np.random.default_rng(1701)

    for _, row in meta.iterrows():
        sid = row["sample_id"]
        stage = row["stage_inferred"]
        print(f"[v17] loading {sid} ({stage})")
        adata = load_scored_h5ad(row)
        obs = adata.obs.copy()
        coords = obs[["array_row", "array_col"]].astype(float).to_numpy()
        epi = obs["Epithelial_score"].astype(float).to_numpy()
        epi_mask = epi >= np.nanmedian(epi)
        mapk, mapk_detect, n_mapk = module_z_and_detect(adata, MAPK_OUT)
        panel, panel_detect, n_panel = module_z_and_detect(adata, PANEL_8)
        covars = covariate_sets(obs, coords, mapk_detect, panel_detect)
        subsets = {"all_spots": np.ones(adata.n_obs, dtype=bool), "epithelial_top50": epi_mask}

        for subset, mask in subsets.items():
            for adjustment in ADJUSTMENT_ORDER:
                cov = covars[adjustment]
                x_lin = residualize(mapk, cov)
                y_lin = residualize(panel, cov)
                rho_lin, p_lin, n_lin = spearman_block(x_lin[mask], y_lin[mask])
                x_rank = partial_spearman_residuals(mapk, cov)
                y_rank = partial_spearman_residuals(panel, cov)
                rho_rank, p_rank, n_rank = pearson_block(x_rank[mask], y_rank[mask])
                component_rows.append(
                    {
                        "sample_id": sid,
                        "stage": stage,
                        "subset": subset,
                        "adjustment": adjustment,
                        "rho_linear_resid_spearman": rho_lin,
                        "p_linear_resid_spearman": p_lin,
                        "n_linear": n_lin,
                        "rho_rank_partial_spearman": rho_rank,
                        "p_rank_partial_spearman": p_rank,
                        "n_rank": n_rank,
                        "n_MAPK_genes": n_mapk,
                        "n_panel_genes": n_panel,
                    }
                )

        for family in ["epithelial", "qc", "microenv", "spatial", "detection", "epithelial_qc_microenv", "epithelial_qc_microenv_spatial_detection"]:
            cov = covars[family]
            for target, y in {"MAPK_z": mapk, "Panel8_z": panel}.items():
                r2, n, k = lin_r2(y[epi_mask], cov[epi_mask])
                r2_rows.append({"sample_id": sid, "stage": stage, "covariate_family": family, "target": target, "r2": r2, "n": n, "n_covariates": k})

        qc_pairs = {
            "MAPK_z_vs_log1p_total": (mapk, np.log1p(obs["total_counts"].astype(float).to_numpy())),
            "Panel8_z_vs_log1p_total": (panel, np.log1p(obs["total_counts"].astype(float).to_numpy())),
            "MAPK_z_vs_MAPK_detect": (mapk, mapk_detect),
            "Panel8_z_vs_panel_detect": (panel, panel_detect),
            "MAPK_detect_vs_panel_detect": (mapk_detect, panel_detect),
            "MAPK_z_vs_Epithelial_score": (mapk, epi),
            "Panel8_z_vs_Epithelial_score": (panel, epi),
        }
        for metric, (x, y) in qc_pairs.items():
            rho, p, n = spearman_block(x[epi_mask], y[epi_mask])
            qc_rows.append({"sample_id": sid, "stage": stage, "metric": metric, "rho": rho, "p": p, "n": n})

        full_cov = covars["epithelial_qc_microenv"]
        for module_name, genes in MAPK_MODULES.items():
            mod, _, n_genes = module_z_and_detect(adata, genes)
            for adjustment, cov in {"raw": None, "full": full_cov}.items():
                x = residualize(mod, cov)
                y = residualize(panel, cov)
                rho, p, n = spearman_block(x[epi_mask], y[epi_mask])
                submodule_rows.append(
                    {
                        "sample_id": sid,
                        "stage": stage,
                        "module": module_name,
                        "adjustment": adjustment,
                        "rho_mapk_module_panel": rho,
                        "p": p,
                        "n": n,
                        "n_module_genes": n_genes,
                    }
                )

        if stage != "PT":
            full_cov_epi = full_cov
            mapk_resid = {gene: residualize(get_expr(adata, gene), full_cov_epi) for gene in MAPK_OUT}
            panel_resid = {gene: residualize(get_expr(adata, gene), full_cov_epi) for gene in PANEL_8}
            for mapk_gene, x in mapk_resid.items():
                for panel_gene, y in panel_resid.items():
                    rho, p, n = spearman_block(x[epi_mask], y[epi_mask])
                    gene_rows.append({"sample_id": sid, "stage": stage, "mapk_gene": mapk_gene, "panel_gene": panel_gene, "rho": rho, "p": p, "n": n})

            genes = eligible_random_genes(adata, set(PANEL_8 + MAPK_OUT))
            obs_raw, _, _ = spearman_block(mapk[epi_mask], panel[epi_mask])
            obs_full, _, _ = spearman_block(
                residualize(mapk, full_cov)[epi_mask],
                residualize(panel, full_cov)[epi_mask],
            )
            for i in range(RANDOM_N):
                if len(genes) < 17:
                    break
                picked = rng.choice(genes, size=17, replace=False)
                rand_mapk, _, _ = module_z_and_detect(adata, list(picked[:9]))
                rand_panel, _, _ = module_z_and_detect(adata, list(picked[9:]))
                rho_raw, p_raw, n_raw = spearman_block(rand_mapk[epi_mask], rand_panel[epi_mask])
                rho_full, p_full, n_full = spearman_block(
                    residualize(rand_mapk, full_cov)[epi_mask],
                    residualize(rand_panel, full_cov)[epi_mask],
                )
                random_rows.append(
                    {
                        "sample_id": sid,
                        "stage": stage,
                        "iteration": i,
                        "rho_random_raw": rho_raw,
                        "p_random_raw": p_raw,
                        "n_raw": n_raw,
                        "rho_random_full_resid": rho_full,
                        "p_random_full_resid": p_full,
                        "n_full": n_full,
                        "observed_raw_rho": obs_raw,
                        "observed_full_resid_rho": obs_full,
                    }
                )

        pooled_chunks.append(
            pd.DataFrame(
                {
                    "sample_id": sid,
                    "stage": stage,
                    "epithelial_top50": epi_mask.astype(int),
                    "MAPK_z": mapk,
                    "Panel8_z": panel,
                    "mapk_detect": mapk_detect,
                    "panel_detect": panel_detect,
                    "Epithelial_score": epi,
                    "log1p_total_counts": np.log1p(obs["total_counts"].astype(float).to_numpy()),
                    "log1p_n_genes": np.log1p(obs["n_genes_by_counts"].astype(float).to_numpy()),
                    "pct_counts_mt": obs["pct_counts_mt"].astype(float).to_numpy(),
                    "CAF_ECM_score": obs["CAF_ECM_score"].astype(float).to_numpy(),
                    "EMT_score": obs["EMT_score"].astype(float).to_numpy(),
                    "Hypoxia_score": obs["Hypoxia_score"].astype(float).to_numpy(),
                    "Proliferation_score": obs["Proliferation_score"].astype(float).to_numpy(),
                    "array_row": coords[:, 0],
                    "array_col": coords[:, 1],
                }
            )
        )

    components = pd.DataFrame(component_rows)
    r2_df = pd.DataFrame(r2_rows)
    qc_df = pd.DataFrame(qc_rows)
    submodules = pd.DataFrame(submodule_rows)
    gene_pairs = pd.DataFrame(gene_rows)
    random_null = pd.DataFrame(random_rows)
    pooled_spots = pd.concat(pooled_chunks, ignore_index=True)

    gene_summary = (
        gene_pairs.groupby(["mapk_gene", "panel_gene"], as_index=False)
        .agg(
            median_rho=("rho", "median"),
            mean_rho=("rho", "mean"),
            n_samples=("sample_id", "nunique"),
            n_negative=("rho", lambda s: int((s < 0).sum())),
            n_positive=("rho", lambda s: int((s > 0).sum())),
        )
        .sort_values(["median_rho", "n_negative"])
    )
    gene_summary["frac_negative"] = gene_summary["n_negative"] / gene_summary["n_samples"]
    gene_summary["binom_p_negative"] = gene_summary.apply(
        lambda r: stats.binomtest(int(r["n_negative"]), int(r["n_samples"]), p=0.5, alternative="greater").pvalue
        if r["n_samples"] > 0
        else np.nan,
        axis=1,
    )
    gene_summary["binom_fdr_negative"] = bh_fdr(gene_summary["binom_p_negative"])

    pooled_rows = []
    pooled_tumor_epi = pooled_spots[(pooled_spots["stage"] != "PT") & (pooled_spots["epithelial_top50"] == 1)].copy()
    pooled_atc_epi = pooled_spots[(pooled_spots["stage"] == "ATC") & (pooled_spots["epithelial_top50"] == 1)].copy()
    for subset, sub in {"tumor_epithelial_top50": pooled_tumor_epi, "ATC_epithelial_top50": pooled_atc_epi}.items():
        obs_proxy = pd.DataFrame(
            {
                "Epithelial_score": sub["Epithelial_score"],
                "total_counts": np.expm1(sub["log1p_total_counts"]),
                "n_genes_by_counts": np.expm1(sub["log1p_n_genes"]),
                "pct_counts_mt": sub["pct_counts_mt"],
                "CAF_ECM_score": sub["CAF_ECM_score"],
                "EMT_score": sub["EMT_score"],
                "Hypoxia_score": sub["Hypoxia_score"],
                "Proliferation_score": sub["Proliferation_score"],
            }
        )
        coords = sub[["array_row", "array_col"]].to_numpy()
        covs = covariate_sets(obs_proxy, coords, sub["mapk_detect"].to_numpy(), sub["panel_detect"].to_numpy())
        for adjustment in ADJUSTMENT_ORDER:
            cov = covs[adjustment]
            x_lin = residualize(sub["MAPK_z"].to_numpy(), cov)
            y_lin = residualize(sub["Panel8_z"].to_numpy(), cov)
            rho_lin, p_lin, n_lin = spearman_block(x_lin, y_lin)
            x_rank = partial_spearman_residuals(sub["MAPK_z"].to_numpy(), cov)
            y_rank = partial_spearman_residuals(sub["Panel8_z"].to_numpy(), cov)
            rho_rank, p_rank, n_rank = pearson_block(x_rank, y_rank)
            pooled_rows.append(
                {
                    "subset": subset,
                    "adjustment": adjustment,
                    "rho_linear_resid_spearman": rho_lin,
                    "p_linear_resid_spearman": p_lin,
                    "n_linear": n_lin,
                    "rho_rank_partial_spearman": rho_rank,
                    "p_rank_partial_spearman": p_rank,
                    "n_rank": n_rank,
                }
            )
    pooled_components = pd.DataFrame(pooled_rows)

    components.to_csv(OUT / "v17_spatial_component_ladder_per_sample.tsv", sep="\t", index=False)
    pooled_components.to_csv(OUT / "v17_spatial_component_ladder_pooled.tsv", sep="\t", index=False)
    r2_df.to_csv(OUT / "v17_spatial_covariate_r2.tsv", sep="\t", index=False)
    qc_df.to_csv(OUT / "v17_spatial_qc_detection_correlations.tsv", sep="\t", index=False)
    submodules.to_csv(OUT / "v17_spatial_mapk_submodule_correlations.tsv", sep="\t", index=False)
    gene_pairs.to_csv(OUT / "v17_spatial_gene_pair_sign_per_sample.tsv", sep="\t", index=False)
    gene_summary.to_csv(OUT / "v17_spatial_gene_pair_sign_summary.tsv", sep="\t", index=False)
    random_null.to_csv(OUT / "v17_spatial_random_module_null.tsv", sep="\t", index=False)

    focus = pooled_components.set_index(["subset", "adjustment"])
    full_key = ("tumor_epithelial_top50", "epithelial_qc_microenv")
    detect_key = ("tumor_epithelial_top50", "epithelial_qc_microenv_spatial_detection")
    random_raw = random_null["rho_random_raw"].dropna()
    random_full = random_null["rho_random_full_resid"].dropna()
    obs_raw_med = random_null.drop_duplicates("sample_id")["observed_raw_rho"].median()
    obs_full_med = random_null.drop_duplicates("sample_id")["observed_full_resid_rho"].median()
    r2_focus = r2_df[(r2_df["stage"] != "PT") & (r2_df["covariate_family"] == "epithelial_qc_microenv_spatial_detection")]
    qc_focus = qc_df[(qc_df["stage"] != "PT") & (qc_df["metric"] == "MAPK_detect_vs_panel_detect")]
    top_gene = gene_summary.sort_values(["median_rho", "binom_fdr_negative"]).iloc[0]
    metrics = {
        "n_samples": int(meta.shape[0]),
        "n_tumor_samples": int((meta["stage_inferred"] != "PT").sum()),
        "n_random_iterations_per_tumor_sample": RANDOM_N,
        "tumor_epi_raw_linear_rho": float(focus.loc[("tumor_epithelial_top50", "raw"), "rho_linear_resid_spearman"]),
        "tumor_epi_full_linear_rho": float(focus.loc[full_key, "rho_linear_resid_spearman"]),
        "tumor_epi_full_rank_partial_rho": float(focus.loc[full_key, "rho_rank_partial_spearman"]),
        "tumor_epi_full_spatial_detection_linear_rho": float(focus.loc[detect_key, "rho_linear_resid_spearman"]),
        "ATC_epi_full_rank_partial_rho": float(focus.loc[("ATC_epithelial_top50", "epithelial_qc_microenv"), "rho_rank_partial_spearman"]),
        "median_full_covariate_r2_MAPK_tumor_epi": float(r2_focus[r2_focus["target"] == "MAPK_z"]["r2"].median()),
        "median_full_covariate_r2_Panel8_tumor_epi": float(r2_focus[r2_focus["target"] == "Panel8_z"]["r2"].median()),
        "median_MAPK_detect_vs_panel_detect_rho_tumor_epi": float(qc_focus["rho"].median()),
        "observed_per_sample_median_raw_rho": float(obs_raw_med),
        "observed_per_sample_median_full_resid_rho": float(obs_full_med),
        "random_module_raw_median_rho": float(random_raw.median()),
        "random_module_full_resid_median_rho": float(random_full.median()),
        "observed_raw_percentile_vs_random": float((random_raw <= obs_raw_med).mean()),
        "observed_full_percentile_vs_random": float((random_full <= obs_full_med).mean()),
        "top_negative_gene_pair": f"{top_gene['mapk_gene']}_vs_{top_gene['panel_gene']}",
        "top_negative_gene_pair_median_rho": float(top_gene["median_rho"]),
        "top_negative_gene_pair_frac_negative": float(top_gene["frac_negative"]),
        "top_negative_gene_pair_fdr": float(top_gene["binom_fdr_negative"]),
    }
    (OUT / "v17_spatial_signal_decomposition_summary.json").write_text(json.dumps(metrics, indent=2))
    return components, pooled_components, r2_df, qc_df, submodules, gene_summary, random_null, metrics


def plot_v17(
    components: pd.DataFrame,
    pooled: pd.DataFrame,
    r2_df: pd.DataFrame,
    qc_df: pd.DataFrame,
    submodules: pd.DataFrame,
    gene_summary: pd.DataFrame,
    random_null: pd.DataFrame,
    metrics: dict,
) -> None:
    print("[v17] plotting")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    ax1, ax2, ax3, ax4, ax5, ax6 = axes.ravel()

    focus = pooled[pooled["subset"] == "tumor_epithelial_top50"].copy()
    focus["adjustment"] = pd.Categorical(focus["adjustment"], ADJUSTMENT_ORDER, ordered=True)
    focus = focus.sort_values("adjustment")
    x = np.arange(len(focus))
    ax1.plot(x, focus["rho_linear_resid_spearman"], marker="o", lw=2.5, color="#b2182b", label="linear residual + Spearman")
    ax1.plot(x, focus["rho_rank_partial_spearman"], marker="s", lw=2.5, color="#2166ac", label="rank-partial Spearman")
    ax1.axhline(0, color="#666", lw=0.8, ls=":")
    ax1.set_xticks(x)
    ax1.set_xticklabels([ADJUSTMENT_LABELS[a] for a in focus["adjustment"]], rotation=35, ha="right", fontsize=9)
    ax1.set_ylabel("Pooled tumor epithelial rho")
    ax1.set_title("A. Adjustment-family decomposition", loc="left", weight="bold")
    ax1.legend(frameon=False, fontsize=9)
    ax1.grid(axis="y", alpha=0.25)

    r2_focus = r2_df[r2_df["stage"] != "PT"].copy()
    r2_focus = r2_focus[r2_focus["covariate_family"].isin(["epithelial", "qc", "microenv", "spatial", "detection", "epithelial_qc_microenv_spatial_detection"])]
    families = ["epithelial", "qc", "microenv", "spatial", "detection", "epithelial_qc_microenv_spatial_detection"]
    width = 0.36
    for offset, target, color in [(-width / 2, "MAPK_z", "#ef8a62"), (width / 2, "Panel8_z", "#67a9cf")]:
        vals = [r2_focus[(r2_focus["target"] == target) & (r2_focus["covariate_family"] == fam)]["r2"].median() for fam in families]
        ax2.bar(np.arange(len(families)) + offset, vals, width=width, color=color, edgecolor="#111", label=target)
    ax2.set_xticks(np.arange(len(families)))
    ax2.set_xticklabels(["Epi", "QC", "Microenv", "Spatial", "Detect", "All"], rotation=35, ha="right")
    ax2.set_ylabel("Median per-sample R2")
    ax2.set_title("B. Covariates explain both modules", loc="left", weight="bold")
    ax2.legend(frameon=False, fontsize=9)
    ax2.grid(axis="y", alpha=0.25)

    qc_focus = qc_df[qc_df["stage"] != "PT"].copy()
    metrics_order = [
        "MAPK_z_vs_log1p_total",
        "Panel8_z_vs_log1p_total",
        "MAPK_z_vs_MAPK_detect",
        "Panel8_z_vs_panel_detect",
        "MAPK_detect_vs_panel_detect",
        "MAPK_z_vs_Epithelial_score",
        "Panel8_z_vs_Epithelial_score",
    ]
    labels = ["MAPK~depth", "Panel~depth", "MAPK~det", "Panel~det", "det~det", "MAPK~epi", "Panel~epi"]
    data = [qc_focus[qc_focus["metric"] == m]["rho"].dropna().to_numpy() for m in metrics_order]
    ax3.boxplot(data, labels=labels, showfliers=False, patch_artist=True, boxprops={"facecolor": "#f6c85f", "alpha": 0.75})
    ax3.axhline(0, color="#666", lw=0.8, ls=":")
    ax3.set_ylabel("Per-sample Spearman rho")
    ax3.set_title("C. Depth/detection structure drives co-localization", loc="left", weight="bold")
    ax3.tick_params(axis="x", rotation=35)
    ax3.grid(axis="y", alpha=0.25)

    sub_focus = submodules[(submodules["stage"] != "PT")].groupby(["module", "adjustment"], as_index=False)["rho_mapk_module_panel"].median()
    modules = list(MAPK_MODULES.keys())
    width = 0.36
    for offset, adj, color in [(-width / 2, "raw", "#b2182b"), (width / 2, "full", "#2166ac")]:
        vals = [sub_focus[(sub_focus["module"] == m) & (sub_focus["adjustment"] == adj)]["rho_mapk_module_panel"].iloc[0] for m in modules]
        ax4.bar(np.arange(len(modules)) + offset, vals, width=width, color=color, edgecolor="#111", label=adj)
    ax4.axhline(0, color="#666", lw=0.8, ls=":")
    ax4.set_xticks(np.arange(len(modules)))
    ax4.set_xticklabels(["All9", "DUSP/SPRY", "ETV/PHLDA1", "No CCND1", "CCND1"], rotation=35, ha="right")
    ax4.set_ylabel("Median per-sample rho vs Panel-8")
    ax4.set_title("D. MAPK submodules: no hidden module-level rescue", loc="left", weight="bold")
    ax4.legend(frameon=False, fontsize=9)
    ax4.grid(axis="y", alpha=0.25)

    gene_top = gene_summary.sort_values(["median_rho", "frac_negative"]).head(12).iloc[::-1]
    y = np.arange(len(gene_top))
    ax5.barh(y, gene_top["median_rho"], color="#2166ac", edgecolor="#111", alpha=0.9)
    ax5.axvline(0, color="#666", lw=0.8, ls=":")
    ax5.set_yticks(y)
    ax5.set_yticklabels([f"{r.mapk_gene} vs {r.panel_gene} ({int(r.n_negative)}/{int(r.n_samples)})" for _, r in gene_top.iterrows()], fontsize=8)
    ax5.set_xlabel("Median full-residual rho")
    ax5.set_title("E. Gene-pair sign consistency is weak", loc="left", weight="bold")
    ax5.grid(axis="x", alpha=0.25)

    ax6.hist(random_null["rho_random_raw"].dropna(), bins=40, alpha=0.62, color="#ef8a62", label="random raw")
    ax6.hist(random_null["rho_random_full_resid"].dropna(), bins=40, alpha=0.62, color="#67a9cf", label="random full-resid")
    ax6.axvline(metrics["observed_per_sample_median_raw_rho"], color="#7f0000", lw=2, label="observed raw median")
    ax6.axvline(metrics["observed_per_sample_median_full_resid_rho"], color="#08306b", lw=2, label="observed full median")
    ax6.axvline(0, color="#666", lw=0.8, ls=":")
    ax6.set_xlabel("Per-sample epithelial-top50 module-module rho")
    ax6.set_ylabel("Random draws")
    ax6.set_title("F. Random-module specificity null", loc="left", weight="bold")
    ax6.legend(frameon=False, fontsize=8)
    ax6.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "Figure SX_v17. Spatial MAPK-panel signal decomposition: positive Visium co-localization is broad detection/covariate structure, not a rescued anti-correlation.",
        fontsize=13,
        weight="bold",
        y=1.01,
    )
    fig.tight_layout()
    png = OUT / "Fig_SX_v17_spatial_signal_decomposition.png"
    pdf = OUT / "Fig_SX_v17_spatial_signal_decomposition.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    ASSETS.mkdir(parents=True, exist_ok=True)
    LIVE_ASSETS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(png, ASSETS / png.name)
    shutil.copy2(pdf, ASSETS / pdf.name)
    shutil.copy2(png, LIVE_ASSETS / png.name)
    shutil.copy2(pdf, LIVE_ASSETS / pdf.name)


def fmt(x: float, nd: int = 3) -> str:
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def table_html(df: pd.DataFrame, max_rows: int = 12) -> str:
    small = df.head(max_rows).copy()
    headers = "".join(f"<th>{html.escape(str(c))}</th>" for c in small.columns)
    rows = []
    for _, row in small.iterrows():
        cells = []
        for val in row:
            text = f"{val:.4g}" if isinstance(val, float) else str(val)
            cells.append(f"<td>{html.escape(text)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{headers}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def build_html(pooled: pd.DataFrame, r2_df: pd.DataFrame, qc_df: pd.DataFrame, submodules: pd.DataFrame, gene_summary: pd.DataFrame, random_null: pd.DataFrame, metrics: dict) -> None:
    print("[v17] building/deploying HTML")
    pooled_focus = pooled[pooled["subset"] == "tumor_epithelial_top50"].copy()
    pooled_focus["adjustment"] = pd.Categorical(pooled_focus["adjustment"], ADJUSTMENT_ORDER, ordered=True)
    pooled_focus = pooled_focus.sort_values("adjustment")[
        ["adjustment", "rho_linear_resid_spearman", "rho_rank_partial_spearman", "n_linear"]
    ]
    r2_summary = (
        r2_df[r2_df["stage"] != "PT"]
        .groupby(["covariate_family", "target"], as_index=False)["r2"]
        .median()
        .sort_values(["covariate_family", "target"])
    )
    qc_summary = (
        qc_df[qc_df["stage"] != "PT"]
        .groupby("metric", as_index=False)["rho"]
        .median()
        .sort_values("rho", ascending=False)
    )
    sub_summary = (
        submodules[submodules["stage"] != "PT"]
        .groupby(["module", "adjustment"], as_index=False)["rho_mapk_module_panel"]
        .median()
        .sort_values(["module", "adjustment"])
    )
    gene_top = gene_summary.sort_values(["median_rho", "binom_fdr_negative"]).head(12)
    random_summary = pd.DataFrame(
        [
            {
                "metric": "observed raw per-sample median",
                "rho": metrics["observed_per_sample_median_raw_rho"],
                "percentile_vs_random": metrics["observed_raw_percentile_vs_random"],
            },
            {
                "metric": "random raw median",
                "rho": metrics["random_module_raw_median_rho"],
                "percentile_vs_random": np.nan,
            },
            {
                "metric": "observed full-resid per-sample median",
                "rho": metrics["observed_per_sample_median_full_resid_rho"],
                "percentile_vs_random": metrics["observed_full_percentile_vs_random"],
            },
            {
                "metric": "random full-resid median",
                "rho": metrics["random_module_full_resid_median_rho"],
                "percentile_vs_random": np.nan,
            },
        ]
    )

    css = """
    :root{--bg:#0d1117;--panel:#151b23;--line:#30363d;--text:#e6edf3;--muted:#9ba7b4;--gold:#f6c85f;--blue:#67a9cf;--red:#ef8a62}
    *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--text);font-family:"JetBrains Mono",ui-monospace,Menlo,monospace;line-height:1.55}
    .hero{padding:56px 7vw 32px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#111a24,#0d1117)}
    .kicker{color:var(--gold);text-transform:uppercase;letter-spacing:.08em;font-size:12px;font-weight:700}
    h1{font-family:"Cormorant Garamond",Georgia,serif;font-size:44px;line-height:1.05;margin:10px 0 14px}
    .lead{max-width:1100px;color:#c8d1dc}.crumbs{margin-top:18px;color:var(--muted);font-size:12px}.crumbs a{color:var(--blue)}
    .stats{display:grid;grid-template-columns:repeat(6,minmax(120px,1fr));gap:10px;margin-top:24px}.stat{border:1px solid var(--line);background:#111821;padding:12px;border-radius:6px}
    .stat b{display:block;color:var(--gold);font-size:20px}.stat span{font-size:11px;color:var(--muted)}
    .wrap{display:grid;grid-template-columns:260px 1fr;gap:28px;padding:26px 7vw 60px}.toc{position:sticky;top:16px;align-self:start;border:1px solid var(--line);background:var(--panel);border-radius:6px;padding:14px}
    .toc a{display:block;color:#c8d1dc;text-decoration:none;font-size:12px;padding:6px 0;border-bottom:1px solid rgba(255,255,255,.04)}
    section{margin-bottom:28px;border-top:1px solid var(--line);padding-top:20px}.num{color:var(--gold);font-weight:700;margin-right:8px} h2{font-size:22px;margin:0 0 12px}
    .box{border:1px solid var(--line);background:var(--panel);border-radius:6px;padding:16px;margin:14px 0}.warn{border-color:#7c4b17;background:#21170c}.good{border-color:#345d40;background:#101d16}
    table{width:100%;border-collapse:collapse;font-size:12px;margin:12px 0} th,td{border:1px solid var(--line);padding:7px 8px;text-align:left} th{color:var(--gold);background:#111821}
    img{max-width:100%;border:1px solid var(--line);border-radius:6px;background:#fff}.cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}.card{border:1px solid var(--line);background:#111821;border-radius:6px;padding:13px;text-decoration:none;color:var(--text)}.card b{display:block;color:var(--gold);margin-bottom:4px}.card span{font-size:12px;color:var(--muted)}
    code{color:#ffd98a}.small{font-size:12px;color:var(--muted)}@media(max-width:900px){.wrap{grid-template-columns:1fr}.toc{position:relative}.stats,.cards{grid-template-columns:1fr}h1{font-size:34px}}
    """
    html_text = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Paper 1 v17 spatial signal decomposition</title><style>{css}</style></head>
<body>
<div class="hero">
  <div class="kicker">Paper 1 · Fig 8 reserve · v17</div>
  <h1>Spatial signal decomposition dossier</h1>
  <p class="lead">v17 decomposes the failed GSE250521 spatial anti-correlation test. The result is stronger as a caveat: the raw positive MAPK-panel co-localization is broad detection/covariate structure, and rank-partial plus random-module checks do not reveal a hidden module-level anti-correlation.</p>
  <div class="stats">
    <div class="stat"><b>{fmt(metrics['tumor_epi_raw_linear_rho'])}</b><span>raw tumor epi rho</span></div>
    <div class="stat"><b>{fmt(metrics['tumor_epi_full_rank_partial_rho'])}</b><span>full rank-partial rho</span></div>
    <div class="stat"><b>{fmt(metrics['median_MAPK_detect_vs_panel_detect_rho_tumor_epi'])}</b><span>detect-breadth rho</span></div>
    <div class="stat"><b>{fmt(metrics['random_module_raw_median_rho'])}</b><span>random raw median</span></div>
    <div class="stat"><b>{fmt(metrics['random_module_full_resid_median_rho'])}</b><span>random full median</span></div>
    <div class="stat"><b>{fmt(metrics['top_negative_gene_pair_median_rho'])}</b><span>{html.escape(metrics['top_negative_gene_pair'])}</span></div>
  </div>
  <div class="crumbs"><a href="index.html">hub</a> / <a href="paper1_spatial_failure_rescue_v16.html">v16 failure rescue</a> / v17 decomposition</div>
</div>
<div class="wrap">
  <nav class="toc">
    <a href="#tldr">TL;DR</a><a href="#figure">Figure</a><a href="#components">Components</a><a href="#qc">QC/detection</a><a href="#modules">Submodules</a><a href="#null">Random null</a><a href="#decision">Decision</a><a href="#sources">Sources</a>
  </nav>
  <main>
    <section id="tldr"><h2><span class="num">01</span>TL;DR</h2>
      <div class="box warn"><b>No hidden rescue.</b> Rank-partial Spearman and submodule checks do not recover the hypothesized within-spot MAPK-high / Panel-low relation. The honest result remains caveat/reserve.</div>
      <div class="box good"><b>Better defense.</b> The raw positive result behaves like generic detection/covariate co-localization: random module pairs also show positive raw correlations, then collapse toward null after the same adjustment.</div>
    </section>
    <section id="figure"><h2><span class="num">02</span>Composite Figure</h2>
      <img src="assets/paper1/Fig_SX_v17_spatial_signal_decomposition.png" alt="v17 spatial signal decomposition figure">
    </section>
    <section id="components"><h2><span class="num">03</span>Adjustment Components</h2>
      {table_html(pooled_focus, max_rows=20)}
      <h3>Covariate R2</h3>
      {table_html(r2_summary, max_rows=20)}
    </section>
    <section id="qc"><h2><span class="num">04</span>QC and Detection</h2>
      {table_html(qc_summary, max_rows=12)}
    </section>
    <section id="modules"><h2><span class="num">05</span>MAPK Submodules + Gene Pairs</h2>
      {table_html(sub_summary, max_rows=14)}
      {table_html(gene_top, max_rows=12)}
    </section>
    <section id="null"><h2><span class="num">06</span>Random-Module Null</h2>
      <p class="small">Random draws use 9-gene and 8-gene non-mitochondrial expressed modules per tumor sample, repeated {RANDOM_N} times per sample.</p>
      {table_html(random_summary, max_rows=8)}
    </section>
    <section id="decision"><h2><span class="num">07</span>Decision</h2>
      <table><thead><tr><th>Question</th><th>Answer</th><th>Disposition</th></tr></thead><tbody>
        <tr><td>Can GSE250521 be turned into a positive Fig 8 mechanism result?</td><td>No. Full rank-partial rho is {fmt(metrics['tumor_epi_full_rank_partial_rho'])}.</td><td>Keep out of positive mechanism chain.</td></tr>
        <tr><td>Can the failure be defended?</td><td>Yes. Detection breadth and covariates explain broad positive co-localization; random modules show the same behavior.</td><td>Use only as caveat/reserve.</td></tr>
        <tr><td>Is there any local foothold?</td><td>Weak gene-pair hints, top {html.escape(metrics['top_negative_gene_pair'])} rho {fmt(metrics['top_negative_gene_pair_median_rho'])}.</td><td>Do not overclaim.</td></tr>
      </tbody></table>
      <div class="cards">
        <a class="card" href="paper1_spatial_failure_rescue_v16.html"><b>v16 failure rescue</b><span>Adjustment ladder and near-null containment.</span></a>
        <a class="card" href="paper1_fig8_mechanism_dossier.html"><b>Fig 8 mechanism dossier</b><span>Primary positive mechanism layer.</span></a>
        <a class="card" href="paper1_methods_reproducibility_dossier.html"><b>Methods reproducibility</b><span>Gene lists, formulas, and cohort filters.</span></a>
      </div>
    </section>
    <section id="sources"><h2><span class="num">08</span>Sources + Paths</h2>
      <div class="box small">
        <p><code>project/data/processed/GSE250521/*/*.scored.h5ad</code></p>
        <p><code>project/results/p_deconv_2026_05_08/run_v17_spatial_signal_decomposition.py</code></p>
        <p><code>v17_spatial_component_ladder_pooled.tsv</code>, <code>v17_spatial_covariate_r2.tsv</code>, <code>v17_spatial_random_module_null.tsv</code>, <code>v17_spatial_signal_decomposition_summary.json</code></p>
      </div>
    </section>
  </main>
</div></body></html>
"""
    out_html = HUB / "paper1_spatial_signal_decomposition_v17.html"
    out_html.write_text(html_text)
    LIVE.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out_html, LIVE / out_html.name)


def update_index() -> None:
    index = HUB / "index.html"
    if not index.exists():
        return
    text = index.read_text()
    link = '<a href="paper1_spatial_signal_decomposition_v17.html" style="background:linear-gradient(180deg,#07101a 0%,#0d1422 100%);border-color:#67a9cf;color:#9dd8f2;font-weight:700">★ Paper 1 · Spatial signal decomposition (v17)</a>'
    if "paper1_spatial_signal_decomposition_v17.html" not in text:
        marker = '<a href="paper1_spatial_failure_rescue_v16.html" style="background:linear-gradient(180deg,#07101a 0%,#0d1422 100%);border-color:#67a9cf;color:#9dd8f2;font-weight:700">★ Paper 1 · Spatial failure rescue (v16)</a>'
        if marker in text:
            text = text.replace(marker, marker + "\n    " + link)
        else:
            text = text.replace("</div>\n</header>", f"    {link}\n  </div>\n</header>", 1)
        index.write_text(text)
        if LIVE.exists():
            shutil.copy2(index, LIVE / "index.html")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    components, pooled, r2_df, qc_df, submodules, gene_summary, random_null, metrics = run_v17()
    plot_v17(components, pooled, r2_df, qc_df, submodules, gene_summary, random_null, metrics)
    build_html(pooled, r2_df, qc_df, submodules, gene_summary, random_null, metrics)
    update_index()
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
