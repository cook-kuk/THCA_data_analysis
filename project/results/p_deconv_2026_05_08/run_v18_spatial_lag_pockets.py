#!/usr/bin/env python3
"""v18: spatial-lag and MAPK-high/Panel-low pocket analysis for GSE250521.

v17 showed that the same-spot MAPK x Panel-8 signal is mostly generic
co-detection/covariate structure. v18 tests the last plausible spatial rescue:
maybe MAPK-high spots do not anti-correlate in the same spot, but sit next to
Panel-low neighborhoods or form focal MAPK-high/Panel-low pockets.
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
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression

warnings.filterwarnings("ignore")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROC = ROOT / "project" / "data" / "processed" / "GSE250521"
OUT = ROOT / "project" / "results" / "p_deconv_2026_05_08"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
ASSETS = HUB / "assets" / "paper1"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
MAPK_OUT = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]
RINGS = {
    "same_spot": (0, 0),
    "knn_1_6": (1, 6),
    "knn_7_18": (7, 18),
    "knn_19_36": (19, 36),
}
POCKET_MODES = ["raw", "full_spatial_detection_resid"]


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


def standardize_columns(x: np.ndarray) -> np.ndarray:
    cols = []
    for j in range(x.shape[1]):
        col = x[:, j]
        sd = np.nanstd(col)
        if np.isfinite(sd) and sd > 1e-9:
            cols.append((col - np.nanmean(col)) / (sd + 1e-9))
    return np.column_stack(cols) if cols else np.empty((x.shape[0], 0))


def residualize(y: np.ndarray, covars: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    covars = np.asarray(covars, dtype=float)
    out = np.full_like(y, np.nan, dtype=float)
    keep = np.isfinite(y) & np.isfinite(covars).all(axis=1)
    if keep.sum() < 20:
        return out
    xs = standardize_columns(covars[keep])
    if xs.shape[1] == 0:
        out[keep] = y[keep] - np.nanmean(y[keep])
        return out
    out[keep] = y[keep] - LinearRegression().fit(xs, y[keep]).predict(xs)
    return out


def spearman_block(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 20:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def fisher_or(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    table = np.array([[a, b], [c, d]], dtype=float)
    if (table == 0).any():
        table = table + 0.5
    odds = (table[0, 0] * table[1, 1]) / (table[0, 1] * table[1, 0])
    _, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return float(odds), float(p)


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    s = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
    return float((np.mean(a) - np.mean(b)) / (s + 1e-12))


def covariates(obs: pd.DataFrame, coords: np.ndarray, mapk_detect: np.ndarray, panel_detect: np.ndarray) -> np.ndarray:
    row = coords[:, 0]
    col = coords[:, 1]
    return np.column_stack(
        [
            obs["Epithelial_score"].astype(float).to_numpy(),
            obs["CAF_ECM_score"].astype(float).to_numpy(),
            obs["EMT_score"].astype(float).to_numpy(),
            obs["Hypoxia_score"].astype(float).to_numpy(),
            obs["Proliferation_score"].astype(float).to_numpy(),
            np.log1p(obs["total_counts"].astype(float).to_numpy()),
            np.log1p(obs["n_genes_by_counts"].astype(float).to_numpy()),
            obs["pct_counts_mt"].astype(float).to_numpy(),
            mapk_detect,
            panel_detect,
            row,
            col,
            row * row,
            col * col,
            row * col,
        ]
    )


def neighbor_indices(coords: np.ndarray, max_k: int = 36) -> np.ndarray:
    k = min(max_k + 1, len(coords))
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k)
    if idx.ndim == 1:
        idx = idx[:, None]
    return idx


def ring_mean(values: np.ndarray, idx: np.ndarray, ring: tuple[int, int]) -> np.ndarray:
    lo, hi = ring
    if lo == 0 and hi == 0:
        return values
    if idx.shape[1] <= lo:
        return np.full(len(values), np.nan)
    hi = min(hi, idx.shape[1] - 1)
    use = idx[:, lo : hi + 1]
    return np.nanmean(values[use], axis=1)


def pocket_stats(x: np.ndarray, y: np.ndarray, idx: np.ndarray, mask: np.ndarray) -> dict:
    x2 = x[mask]
    y2 = y[mask]
    idx_masked_global = np.where(mask)[0]
    qx75 = np.nanquantile(x2, 0.75)
    qx25 = np.nanquantile(x2, 0.25)
    qy75 = np.nanquantile(y2, 0.75)
    qy25 = np.nanquantile(y2, 0.25)
    high_x = x >= qx75
    low_x = x <= qx25
    high_y = y >= qy75
    low_y = y <= qy25
    anti = high_x & low_y & mask
    concord = high_x & high_y & mask
    lowlow = low_x & low_y & mask
    high_count = int((high_x & mask).sum())
    low_panel_count = int((low_y & mask).sum())
    high_panel_count = int((high_y & mask).sum())
    n = int(mask.sum())
    anti_n = int(anti.sum())
    concord_n = int(concord.sum())
    lowlow_n = int(lowlow.sum())
    expected_anti = (high_count / n) * (low_panel_count / n) if n else np.nan
    expected_concord = (high_count / n) * (high_panel_count / n) if n else np.nan
    anti_frac = anti_n / n if n else np.nan
    concord_frac = concord_n / n if n else np.nan
    a = anti_n
    b = high_count - anti_n
    c = low_panel_count - anti_n
    d = n - a - b - c
    anti_or, anti_p = fisher_or(a, b, c, d)
    a2 = concord_n
    b2 = high_count - concord_n
    c2 = high_panel_count - concord_n
    d2 = n - a2 - b2 - c2
    concord_or, concord_p = fisher_or(a2, b2, c2, d2)

    # KNN clustering ratio for anti pockets inside the analyzed mask.
    ring = idx[:, 1 : min(idx.shape[1], 7)]
    anti_float = anti.astype(float)
    mask_float = mask.astype(float)
    neigh_mask_n = mask_float[ring].sum(axis=1)
    neigh_anti_n = anti_float[ring].sum(axis=1)
    neigh_frac = np.divide(neigh_anti_n, neigh_mask_n, out=np.full(len(anti_float), np.nan), where=neigh_mask_n > 0)
    pocket_neigh = np.nanmean(neigh_frac[anti]) if anti_n else np.nan
    baseline = anti_frac
    cluster_ratio = pocket_neigh / baseline if np.isfinite(pocket_neigh) and baseline > 0 else np.nan
    return {
        "n_spots": n,
        "anti_pocket_n": anti_n,
        "anti_pocket_frac": anti_frac,
        "anti_expected_frac": expected_anti,
        "anti_enrichment_vs_independence": anti_frac / expected_anti if expected_anti and expected_anti > 0 else np.nan,
        "anti_fisher_or": anti_or,
        "anti_fisher_p": anti_p,
        "anti_knn6_cluster_ratio": cluster_ratio,
        "concord_high_high_n": concord_n,
        "concord_high_high_frac": concord_frac,
        "concord_expected_frac": expected_concord,
        "concord_enrichment_vs_independence": concord_frac / expected_concord if expected_concord and expected_concord > 0 else np.nan,
        "concord_fisher_or": concord_or,
        "concord_fisher_p": concord_p,
        "low_low_n": lowlow_n,
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


def run_v18() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    meta = pd.read_csv(PROC / "sample_metadata.tsv", sep="\t")
    lag_rows = []
    pocket_rows = []
    covar_rows = []
    sample_rows = []

    for _, row in meta.iterrows():
        sid = row["sample_id"]
        stage = row["stage_inferred"]
        print(f"[v18] loading {sid} ({stage})")
        adata = load_scored_h5ad(row)
        obs = adata.obs.copy()
        coords = obs[["array_row", "array_col"]].astype(float).to_numpy()
        idx = neighbor_indices(coords, max_k=36)
        epi = obs["Epithelial_score"].astype(float).to_numpy()
        epi_mask = epi >= np.nanmedian(epi)
        mapk, mapk_detect, n_mapk = module_z_and_detect(adata, MAPK_OUT)
        panel, panel_detect, n_panel = module_z_and_detect(adata, PANEL_8)
        cov = covariates(obs, coords, mapk_detect, panel_detect)
        mapk_resid = residualize(mapk, cov)
        panel_resid = residualize(panel, cov)
        mode_values = {
            "raw": (mapk, panel),
            "full_spatial_detection_resid": (mapk_resid, panel_resid),
        }
        masks = {"all_spots": np.ones(adata.n_obs, dtype=bool), "epithelial_top50": epi_mask}

        for subset, mask in masks.items():
            for mode, (x, y) in mode_values.items():
                for ring_name, ring in RINGS.items():
                    lag_y = ring_mean(y, idx, ring)
                    lag_x = ring_mean(x, idx, ring)
                    rho_xy, p_xy, n_xy = spearman_block(x[mask], lag_y[mask])
                    rho_yx, p_yx, n_yx = spearman_block(y[mask], lag_x[mask])
                    lag_rows.append(
                        {
                            "sample_id": sid,
                            "stage": stage,
                            "subset": subset,
                            "mode": mode,
                            "ring": ring_name,
                            "rho_MAPK_spot_vs_panel_lag": rho_xy,
                            "p_MAPK_spot_vs_panel_lag": p_xy,
                            "n_MAPK_spot_vs_panel_lag": n_xy,
                            "rho_panel_spot_vs_MAPK_lag": rho_yx,
                            "p_panel_spot_vs_MAPK_lag": p_yx,
                            "n_panel_spot_vs_MAPK_lag": n_yx,
                            "n_MAPK_genes": n_mapk,
                            "n_panel_genes": n_panel,
                        }
                    )
                stats_row = pocket_stats(x, y, idx, mask)
                stats_row.update({"sample_id": sid, "stage": stage, "subset": subset, "mode": mode})
                pocket_rows.append(stats_row)

                # Characterize anti-pockets in tumor samples only.
                if stage != "PT" and subset == "epithelial_top50":
                    x2 = x[mask]
                    y2 = y[mask]
                    qx75 = np.nanquantile(x2, 0.75)
                    qy25 = np.nanquantile(y2, 0.25)
                    anti = (x >= qx75) & (y <= qy25) & mask
                    rest = mask & ~anti
                    covar_map = {
                        "Epithelial_score": epi,
                        "log1p_total_counts": np.log1p(obs["total_counts"].astype(float).to_numpy()),
                        "log1p_n_genes": np.log1p(obs["n_genes_by_counts"].astype(float).to_numpy()),
                        "pct_counts_mt": obs["pct_counts_mt"].astype(float).to_numpy(),
                        "CAF_ECM_score": obs["CAF_ECM_score"].astype(float).to_numpy(),
                        "EMT_score": obs["EMT_score"].astype(float).to_numpy(),
                        "Hypoxia_score": obs["Hypoxia_score"].astype(float).to_numpy(),
                        "Proliferation_score": obs["Proliferation_score"].astype(float).to_numpy(),
                        "MAPK_detect": mapk_detect,
                        "Panel_detect": panel_detect,
                    }
                    for covar_name, vals in covar_map.items():
                        covar_rows.append(
                            {
                                "sample_id": sid,
                                "stage": stage,
                                "mode": mode,
                                "covariate": covar_name,
                                "anti_pocket_n": int(anti.sum()),
                                "d_anti_vs_rest": cohen_d(vals[anti], vals[rest]),
                                "mean_anti": float(np.nanmean(vals[anti])) if anti.sum() else np.nan,
                                "mean_rest": float(np.nanmean(vals[rest])) if rest.sum() else np.nan,
                            }
                        )

        sample_rows.append(
            {
                "sample_id": sid,
                "stage": stage,
                "n_spots": int(adata.n_obs),
                "n_MAPK_genes": n_mapk,
                "n_panel_genes": n_panel,
                "mean_mapk_detect": float(np.nanmean(mapk_detect)),
                "mean_panel_detect": float(np.nanmean(panel_detect)),
            }
        )

    lag = pd.DataFrame(lag_rows)
    pockets = pd.DataFrame(pocket_rows)
    covars = pd.DataFrame(covar_rows)
    samples = pd.DataFrame(sample_rows)
    lag.to_csv(OUT / "v18_spatial_lag_correlations.tsv", sep="\t", index=False)
    pockets.to_csv(OUT / "v18_spatial_pocket_enrichment.tsv", sep="\t", index=False)
    covars.to_csv(OUT / "v18_spatial_pocket_covariate_contrasts.tsv", sep="\t", index=False)
    samples.to_csv(OUT / "v18_spatial_sample_qc_summary.tsv", sep="\t", index=False)

    tumor_lag = lag[(lag["stage"] != "PT") & (lag["subset"] == "epithelial_top50")]
    tumor_pockets = pockets[(pockets["stage"] != "PT") & (pockets["subset"] == "epithelial_top50")]
    lag_focus = tumor_lag.groupby(["mode", "ring"], as_index=False)["rho_MAPK_spot_vs_panel_lag"].median()
    pocket_focus = tumor_pockets.groupby(["mode"], as_index=False).agg(
        median_anti_enrichment=("anti_enrichment_vs_independence", "median"),
        median_anti_or=("anti_fisher_or", "median"),
        median_anti_frac=("anti_pocket_frac", "median"),
        median_concord_enrichment=("concord_enrichment_vs_independence", "median"),
        median_concord_or=("concord_fisher_or", "median"),
        median_concord_frac=("concord_high_high_frac", "median"),
        median_anti_cluster_ratio=("anti_knn6_cluster_ratio", "median"),
    )
    covar_focus = (
        covars[covars["mode"] == "raw"]
        .groupby("covariate", as_index=False)["d_anti_vs_rest"]
        .median()
        .sort_values("d_anti_vs_rest")
    )

    def get_lag(mode: str, ring: str) -> float:
        sub = lag_focus[(lag_focus["mode"] == mode) & (lag_focus["ring"] == ring)]
        return float(sub["rho_MAPK_spot_vs_panel_lag"].iloc[0]) if len(sub) else np.nan

    raw_p = pocket_focus[pocket_focus["mode"] == "raw"].iloc[0].to_dict()
    resid_p = pocket_focus[pocket_focus["mode"] == "full_spatial_detection_resid"].iloc[0].to_dict()
    metrics = {
        "n_samples": int(samples.shape[0]),
        "n_tumor_samples": int((samples["stage"] != "PT").sum()),
        "tumor_epi_raw_same_spot_median_rho": get_lag("raw", "same_spot"),
        "tumor_epi_raw_knn_1_6_median_rho": get_lag("raw", "knn_1_6"),
        "tumor_epi_raw_knn_7_18_median_rho": get_lag("raw", "knn_7_18"),
        "tumor_epi_resid_same_spot_median_rho": get_lag("full_spatial_detection_resid", "same_spot"),
        "tumor_epi_resid_knn_1_6_median_rho": get_lag("full_spatial_detection_resid", "knn_1_6"),
        "tumor_epi_resid_knn_7_18_median_rho": get_lag("full_spatial_detection_resid", "knn_7_18"),
        "raw_anti_pocket_median_enrichment": float(raw_p["median_anti_enrichment"]),
        "raw_anti_pocket_median_or": float(raw_p["median_anti_or"]),
        "raw_concord_high_high_median_enrichment": float(raw_p["median_concord_enrichment"]),
        "raw_concord_high_high_median_or": float(raw_p["median_concord_or"]),
        "resid_anti_pocket_median_enrichment": float(resid_p["median_anti_enrichment"]),
        "resid_anti_pocket_median_or": float(resid_p["median_anti_or"]),
        "resid_concord_high_high_median_enrichment": float(resid_p["median_concord_enrichment"]),
        "resid_concord_high_high_median_or": float(resid_p["median_concord_or"]),
        "raw_anti_pocket_cluster_ratio": float(raw_p["median_anti_cluster_ratio"]),
        "resid_anti_pocket_cluster_ratio": float(resid_p["median_anti_cluster_ratio"]),
        "strongest_raw_anti_pocket_covariate_depletion": str(covar_focus.iloc[0]["covariate"]) if len(covar_focus) else "",
        "strongest_raw_anti_pocket_covariate_d": float(covar_focus.iloc[0]["d_anti_vs_rest"]) if len(covar_focus) else np.nan,
    }
    (OUT / "v18_spatial_lag_pockets_summary.json").write_text(json.dumps(metrics, indent=2))
    return lag, pockets, covars, samples, metrics


def plot_v18(lag: pd.DataFrame, pockets: pd.DataFrame, covars: pd.DataFrame, metrics: dict) -> None:
    print("[v18] plotting")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    ax1, ax2, ax3, ax4, ax5, ax6 = axes.ravel()
    tumor_lag = lag[(lag["stage"] != "PT") & (lag["subset"] == "epithelial_top50")]
    rings = list(RINGS.keys())
    labels = ["same", "1-6", "7-18", "19-36"]
    width = 0.36
    for offset, mode, color, label in [
        (-width / 2, "raw", "#b2182b", "raw"),
        (width / 2, "full_spatial_detection_resid", "#2166ac", "full+spatial+detect residual"),
    ]:
        vals = [
            tumor_lag[(tumor_lag["mode"] == mode) & (tumor_lag["ring"] == ring)]["rho_MAPK_spot_vs_panel_lag"].median()
            for ring in rings
        ]
        ax1.bar(np.arange(len(rings)) + offset, vals, width=width, color=color, edgecolor="#111", alpha=0.9, label=label)
    ax1.axhline(0, color="#666", lw=0.8, ls=":")
    ax1.set_xticks(np.arange(len(rings)))
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("Median per-sample rho\nMAPK spot vs Panel neighborhood")
    ax1.set_title("A. Spatial-lag test: no negative neighborhood rescue", loc="left", weight="bold")
    ax1.legend(frameon=False, fontsize=9)
    ax1.grid(axis="y", alpha=0.25)

    tumor_pockets = pockets[(pockets["stage"] != "PT") & (pockets["subset"] == "epithelial_top50")]
    data = []
    tick = []
    colors = []
    for mode in POCKET_MODES:
        sub = tumor_pockets[tumor_pockets["mode"] == mode]
        data.extend([sub["anti_enrichment_vs_independence"].dropna().to_numpy(), sub["concord_enrichment_vs_independence"].dropna().to_numpy()])
        tick.extend([f"{mode}\nMAPKhi/PanelLo", f"{mode}\nMAPKhi/PanelHi"])
        colors.extend(["#67a9cf", "#ef8a62"])
    bp = ax2.boxplot(data, patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax2.axhline(1, color="#666", lw=0.8, ls=":")
    ax2.set_xticklabels(tick, rotation=25, ha="right", fontsize=8)
    ax2.set_ylabel("Enrichment vs independence")
    ax2.set_title("B. MAPK-high/Panel-low pockets are not enriched", loc="left", weight="bold")
    ax2.grid(axis="y", alpha=0.25)

    for stage, color in [("PTC", "#2c7fb8"), ("LPTC", "#e08214"), ("ATC", "#b2182b")]:
        sub = tumor_pockets[(tumor_pockets["stage"] == stage) & (tumor_pockets["mode"] == "raw")]
        ax3.scatter(
            sub["anti_pocket_frac"],
            sub["concord_high_high_frac"],
            s=85,
            color=color,
            edgecolor="#111",
            alpha=0.9,
            label=stage,
        )
    ax3.plot([0, 0.25], [0, 0.25], color="#666", lw=0.8, ls=":")
    ax3.set_xlabel("MAPK-high / Panel-low fraction")
    ax3.set_ylabel("MAPK-high / Panel-high fraction")
    ax3.set_title("C. Concordant high-high pockets dominate raw spots", loc="left", weight="bold")
    ax3.legend(frameon=False, fontsize=9)
    ax3.grid(alpha=0.25)

    cov_focus = (
        covars[covars["mode"] == "raw"]
        .groupby("covariate", as_index=False)["d_anti_vs_rest"]
        .median()
        .sort_values("d_anti_vs_rest")
    )
    top = pd.concat([cov_focus.head(5), cov_focus.tail(5)]).drop_duplicates("covariate")
    ax4.barh(np.arange(len(top)), top["d_anti_vs_rest"], color=np.where(top["d_anti_vs_rest"] < 0, "#2166ac", "#b2182b"), edgecolor="#111")
    ax4.axvline(0, color="#666", lw=0.8, ls=":")
    ax4.set_yticks(np.arange(len(top)))
    ax4.set_yticklabels(top["covariate"], fontsize=9)
    ax4.set_xlabel("Median Cohen's d\nanti-pocket vs rest")
    ax4.set_title("D. Raw anti-pockets are covariate-shaped", loc="left", weight="bold")
    ax4.grid(axis="x", alpha=0.25)

    resid_lag = tumor_lag[tumor_lag["mode"] == "full_spatial_detection_resid"]
    heat = resid_lag.pivot_table(index="sample_id", columns="ring", values="rho_MAPK_spot_vs_panel_lag").reindex(columns=rings)
    im = ax5.imshow(heat, cmap="RdBu_r", vmin=-0.25, vmax=0.25, aspect="auto")
    ax5.set_xticks(np.arange(len(rings)))
    ax5.set_xticklabels(labels)
    ax5.set_yticks(np.arange(len(heat.index)))
    ax5.set_yticklabels(heat.index, fontsize=7)
    ax5.set_title("E. Residual lag by tumor sample", loc="left", weight="bold")
    cbar = fig.colorbar(im, ax=ax5, fraction=0.046, pad=0.04)
    cbar.set_label("rho")

    stage_summary = (
        tumor_pockets[tumor_pockets["mode"] == "full_spatial_detection_resid"]
        .groupby("stage", as_index=False)
        .agg(
            anti_enrichment=("anti_enrichment_vs_independence", "median"),
            anti_cluster=("anti_knn6_cluster_ratio", "median"),
            concord_enrichment=("concord_enrichment_vs_independence", "median"),
        )
    )
    x = np.arange(len(stage_summary))
    w = 0.28
    ax6.bar(x - w, stage_summary["anti_enrichment"], width=w, color="#67a9cf", edgecolor="#111", label="anti enrichment")
    ax6.bar(x, stage_summary["anti_cluster"], width=w, color="#f6c85f", edgecolor="#111", label="anti cluster ratio")
    ax6.bar(x + w, stage_summary["concord_enrichment"], width=w, color="#ef8a62", edgecolor="#111", label="concord enrichment")
    ax6.axhline(1, color="#666", lw=0.8, ls=":")
    ax6.set_xticks(x)
    ax6.set_xticklabels(stage_summary["stage"])
    ax6.set_ylabel("Median ratio")
    ax6.set_title("F. Residual pockets by stage", loc="left", weight="bold")
    ax6.legend(frameon=False, fontsize=8)
    ax6.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "Figure SX_v18. Spatial lag and pocket analysis: no neighborhood-level rescue for the failed MAPK-panel anti-correlation.",
        fontsize=13,
        weight="bold",
        y=1.01,
    )
    fig.tight_layout()
    png = OUT / "Fig_SX_v18_spatial_lag_pockets.png"
    pdf = OUT / "Fig_SX_v18_spatial_lag_pockets.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    ASSETS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(png, ASSETS / png.name)
    shutil.copy2(pdf, ASSETS / pdf.name)


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


def build_html(lag: pd.DataFrame, pockets: pd.DataFrame, covars: pd.DataFrame, metrics: dict) -> None:
    print("[v18] building HTML")
    tumor_lag = lag[(lag["stage"] != "PT") & (lag["subset"] == "epithelial_top50")]
    lag_summary = (
        tumor_lag.groupby(["mode", "ring"], as_index=False)
        .agg(median_rho=("rho_MAPK_spot_vs_panel_lag", "median"), mean_rho=("rho_MAPK_spot_vs_panel_lag", "mean"))
        .sort_values(["mode", "ring"])
    )
    tumor_pockets = pockets[(pockets["stage"] != "PT") & (pockets["subset"] == "epithelial_top50")]
    pocket_summary = (
        tumor_pockets.groupby(["mode"], as_index=False)
        .agg(
            anti_enrichment=("anti_enrichment_vs_independence", "median"),
            anti_or=("anti_fisher_or", "median"),
            anti_frac=("anti_pocket_frac", "median"),
            anti_cluster=("anti_knn6_cluster_ratio", "median"),
            concord_enrichment=("concord_enrichment_vs_independence", "median"),
            concord_or=("concord_fisher_or", "median"),
            concord_frac=("concord_high_high_frac", "median"),
        )
    )
    covar_summary = (
        covars[covars["mode"] == "raw"]
        .groupby("covariate", as_index=False)["d_anti_vs_rest"]
        .median()
        .sort_values("d_anti_vs_rest")
    )

    css = """
    :root{--bg:#0d1117;--panel:#151b23;--line:#30363d;--text:#e6edf3;--muted:#9ba7b4;--gold:#f6c85f;--blue:#67a9cf}
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
<title>Paper 1 v18 spatial lag pockets</title><style>{css}</style></head>
<body>
<div class="hero">
  <div class="kicker">Paper 1 · Fig 8 reserve · v18</div>
  <h1>Spatial lag and pocket dossier</h1>
  <p class="lead">v18 tests the final spatial rescue: whether MAPK-high spots sit next to Panel-low neighborhoods or form focal MAPK-high/Panel-low pockets. The result remains caveat-only: no negative neighborhood-lag rescue and no enrichment of anti-pockets relative to independence.</p>
  <div class="stats">
    <div class="stat"><b>{fmt(metrics['tumor_epi_raw_same_spot_median_rho'])}</b><span>raw same-spot rho</span></div>
    <div class="stat"><b>{fmt(metrics['tumor_epi_raw_knn_1_6_median_rho'])}</b><span>raw KNN 1-6 lag rho</span></div>
    <div class="stat"><b>{fmt(metrics['tumor_epi_resid_knn_1_6_median_rho'])}</b><span>residual KNN 1-6 lag rho</span></div>
    <div class="stat"><b>{fmt(metrics['raw_anti_pocket_median_enrichment'])}</b><span>raw anti-pocket enrichment</span></div>
    <div class="stat"><b>{fmt(metrics['raw_concord_high_high_median_enrichment'])}</b><span>raw concord enrichment</span></div>
    <div class="stat"><b>{fmt(metrics['resid_anti_pocket_median_enrichment'])}</b><span>residual anti enrichment</span></div>
  </div>
  <div class="crumbs"><a href="index.html">hub</a> / <a href="paper1_spatial_signal_decomposition_v17.html">v17 decomposition</a> / v18 spatial lag</div>
</div>
<div class="wrap">
  <nav class="toc">
    <a href="#tldr">TL;DR</a><a href="#figure">Figure</a><a href="#lag">Lag</a><a href="#pockets">Pockets</a><a href="#covars">Covariates</a><a href="#decision">Decision</a><a href="#sources">Sources</a>
  </nav>
  <main>
    <section id="tldr"><h2><span class="num">01</span>TL;DR</h2>
      <div class="box warn"><b>No spatial rescue.</b> MAPK-high spots do not show a hidden negative relation to nearby Panel-low neighborhoods after residualization.</div>
      <div class="box good"><b>Reviewer value.</b> v18 closes the last spatial escape hatch: the failed GSE250521 test is a Visium limitation/caveat, not an alternative support pillar.</div>
    </section>
    <section id="figure"><h2><span class="num">02</span>Composite Figure</h2>
      <img src="assets/paper1/Fig_SX_v18_spatial_lag_pockets.png" alt="v18 spatial lag and pocket figure">
    </section>
    <section id="lag"><h2><span class="num">03</span>Spatial Lag</h2>
      {table_html(lag_summary, max_rows=20)}
    </section>
    <section id="pockets"><h2><span class="num">04</span>Pocket Enrichment</h2>
      {table_html(pocket_summary, max_rows=8)}
    </section>
    <section id="covars"><h2><span class="num">05</span>Anti-Pocket Covariates</h2>
      {table_html(covar_summary, max_rows=12)}
    </section>
    <section id="decision"><h2><span class="num">06</span>Decision</h2>
      <table><thead><tr><th>Question</th><th>Answer</th><th>Disposition</th></tr></thead><tbody>
        <tr><td>Does spatial lag reveal MAPK-high near Panel-low?</td><td>No. Raw KNN 1-6 median rho = {fmt(metrics['tumor_epi_raw_knn_1_6_median_rho'])}; residual KNN 1-6 median rho = {fmt(metrics['tumor_epi_resid_knn_1_6_median_rho'])}.</td><td>No mechanism support.</td></tr>
        <tr><td>Are MAPK-high/Panel-low pockets enriched?</td><td>Raw anti-pocket enrichment = {fmt(metrics['raw_anti_pocket_median_enrichment'])}; residual anti-pocket enrichment = {fmt(metrics['resid_anti_pocket_median_enrichment'])}.</td><td>No pocket rescue.</td></tr>
        <tr><td>What dominates raw spots?</td><td>Concordant MAPK-high/Panel-high enrichment = {fmt(metrics['raw_concord_high_high_median_enrichment'])}.</td><td>Matches co-detection structure.</td></tr>
      </tbody></table>
      <div class="cards">
        <a class="card" href="paper1_spatial_signal_decomposition_v17.html"><b>v17 decomposition</b><span>Random-module and detection-null autopsy.</span></a>
        <a class="card" href="paper1_spatial_failure_rescue_v16.html"><b>v16 failure rescue</b><span>Adjustment ladder and near-null containment.</span></a>
        <a class="card" href="paper1_fig8_mechanism_dossier.html"><b>Fig 8 mechanism dossier</b><span>Primary positive mechanism layer.</span></a>
      </div>
    </section>
    <section id="sources"><h2><span class="num">07</span>Sources + Paths</h2>
      <div class="box small">
        <p><code>project/data/processed/GSE250521/*/*.scored.h5ad</code></p>
        <p><code>project/results/p_deconv_2026_05_08/run_v18_spatial_lag_pockets.py</code></p>
        <p><code>v18_spatial_lag_correlations.tsv</code>, <code>v18_spatial_pocket_enrichment.tsv</code>, <code>v18_spatial_pocket_covariate_contrasts.tsv</code>, <code>v18_spatial_lag_pockets_summary.json</code></p>
      </div>
    </section>
  </main>
</div></body></html>
"""
    out_html = HUB / "paper1_spatial_lag_pockets_v18.html"
    out_html.write_text(html_text)


def update_index() -> None:
    index = HUB / "index.html"
    if not index.exists():
        return
    text = index.read_text()
    link = '<a href="paper1_spatial_lag_pockets_v18.html" style="background:linear-gradient(180deg,#07101a 0%,#0d1422 100%);border-color:#67a9cf;color:#9dd8f2;font-weight:700">★ Paper 1 · Spatial lag pockets (v18)</a>'
    if "paper1_spatial_lag_pockets_v18.html" not in text:
        marker = '<a href="paper1_spatial_signal_decomposition_v17.html" style="background:linear-gradient(180deg,#07101a 0%,#0d1422 100%);border-color:#67a9cf;color:#9dd8f2;font-weight:700">★ Paper 1 · Spatial signal decomposition (v17)</a>'
        if marker in text:
            text = text.replace(marker, marker + "\n    " + link)
        else:
            text = text.replace("</div>\n</header>", f"    {link}\n  </div>\n</header>", 1)
        index.write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    lag, pockets, covars, samples, metrics = run_v18()
    plot_v18(lag, pockets, covars, metrics)
    build_html(lag, pockets, covars, metrics)
    update_index()
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
