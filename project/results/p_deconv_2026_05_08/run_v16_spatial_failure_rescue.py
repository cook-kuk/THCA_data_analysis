#!/usr/bin/env python3
"""v16: rescue/autopsy of the GSE250521 spatial MAPK-panel failure.

The v15 spatial test expected MAPK output to anti-correlate with the
8-gene thyroid differentiation panel within Visium spots. It did not.
This script keeps the result honest and tests whether the positive
spot-level correlation is a broad compartment/QC/spatial-gradient effect
rather than direct support against the Paper 1 mechanism layer.
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
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSETS = HUB / "assets" / "paper1"
LIVE_ASSETS = LIVE / "assets" / "paper1"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
MAPK_OUT = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]
FULL_COVARS = [
    "Epithelial_score",
    "CAF_ECM_score",
    "EMT_score",
    "Hypoxia_score",
    "Proliferation_score",
    "log1p_total_counts",
    "pct_counts_mt",
]
KIND_ORDER = ["raw", "basic_resid", "full_resid", "spatial_poly_resid", "local_knn30_resid"]
KIND_LABELS = {
    "raw": "Raw",
    "basic_resid": "Epi/QC",
    "full_resid": "Full covar",
    "spatial_poly_resid": "Full+spatial",
    "local_knn30_resid": "Local KNN",
}


def get_expr(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    col = adata.X[:, adata.var_names.get_loc(gene)]
    if sp.issparse(col):
        col = col.toarray().ravel()
    else:
        col = np.asarray(col).ravel()
    return col.astype(float)


def module_z(adata, genes: list[str]) -> tuple[np.ndarray, int]:
    vals = []
    for gene in genes:
        v = get_expr(adata, gene)
        if np.all(~np.isfinite(v)):
            continue
        sd = np.nanstd(v)
        vals.append((v - np.nanmean(v)) / (sd + 1e-9) if np.isfinite(sd) and sd > 0 else np.zeros_like(v))
    if not vals:
        return np.full(adata.n_obs, np.nan), 0
    return np.nanmean(np.vstack(vals), axis=0), len(vals)


def spearman_block(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 20:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def residualize(y: np.ndarray, covars: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    covars = np.asarray(covars, dtype=float)
    keep = np.isfinite(y) & np.isfinite(covars).all(axis=1)
    out = np.full_like(y, np.nan, dtype=float)
    if keep.sum() < 20:
        return out
    cols = []
    x = covars[keep]
    for j in range(x.shape[1]):
        col = x[:, j]
        sd = np.nanstd(col)
        if np.isfinite(sd) and sd > 1e-9:
            cols.append((col - np.nanmean(col)) / (sd + 1e-9))
    if not cols:
        out[keep] = y[keep] - np.nanmean(y[keep])
        return out
    xs = np.column_stack(cols)
    out[keep] = y[keep] - LinearRegression().fit(xs, y[keep]).predict(xs)
    return out


def local_knn_residual(y: np.ndarray, coords: np.ndarray, k: int = 30) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    coords = np.asarray(coords, dtype=float)
    keep = np.isfinite(y) & np.isfinite(coords).all(axis=1)
    out = np.full_like(y, np.nan, dtype=float)
    if keep.sum() <= k + 5:
        return out
    tree = cKDTree(coords[keep])
    _, idx = tree.query(coords[keep], k=min(k + 1, keep.sum()))
    neigh = idx[:, 1:]
    out[keep] = y[keep] - np.nanmean(y[keep][neigh], axis=1)
    return out


def covariate_matrix(obs: pd.DataFrame, coords: np.ndarray, flavor: str) -> np.ndarray:
    epithelial = obs["Epithelial_score"].astype(float).to_numpy()
    log_total = np.log1p(obs["total_counts"].astype(float).to_numpy())
    pct_mt = obs["pct_counts_mt"].astype(float).to_numpy()
    if flavor == "basic":
        return np.column_stack([epithelial, log_total, pct_mt])
    full = np.column_stack(
        [
            epithelial,
            obs["CAF_ECM_score"].astype(float).to_numpy(),
            obs["EMT_score"].astype(float).to_numpy(),
            obs["Hypoxia_score"].astype(float).to_numpy(),
            obs["Proliferation_score"].astype(float).to_numpy(),
            log_total,
            pct_mt,
        ]
    )
    if flavor == "full":
        return full
    row = coords[:, 0]
    col = coords[:, 1]
    return np.column_stack([full, row, col, row * row, col * col, row * col])


def make_masks(obs: pd.DataFrame) -> dict[str, np.ndarray]:
    epi = obs["Epithelial_score"].astype(float).to_numpy()
    q25, q50, q75 = np.nanquantile(epi, [0.25, 0.50, 0.75])
    return {
        "all": np.ones(len(obs), dtype=bool),
        "epithelial_enriched_top50": epi >= q50,
        "epithelial_enriched_top75": epi >= q75,
        "epithelial_iqr_25_75": (epi >= q25) & (epi <= q75),
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
    return adata, use_h5


def run_spatial_rescue() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    meta = pd.read_csv(PROC / "sample_metadata.tsv", sep="\t")
    per_sample_rows = []
    pooled_chunks = []
    gene_rows = []

    for _, row in meta.iterrows():
        sid = row["sample_id"]
        stage = row["stage_inferred"]
        print(f"[v16] loading {sid} ({stage})")
        adata, _ = load_scored_h5ad(row)
        obs = adata.obs.copy()
        coords = obs[["array_row", "array_col"]].astype(float).to_numpy()
        masks = make_masks(obs)
        mapk, n_mapk = module_z(adata, MAPK_OUT)
        panel, n_panel = module_z(adata, PANEL_8)

        residual_sets = {
            "raw": (mapk, panel),
            "basic_resid": (residualize(mapk, covariate_matrix(obs, coords, "basic")), residualize(panel, covariate_matrix(obs, coords, "basic"))),
            "full_resid": (residualize(mapk, covariate_matrix(obs, coords, "full")), residualize(panel, covariate_matrix(obs, coords, "full"))),
            "spatial_poly_resid": (
                residualize(mapk, covariate_matrix(obs, coords, "spatial")),
                residualize(panel, covariate_matrix(obs, coords, "spatial")),
            ),
            "local_knn30_resid": (local_knn_residual(mapk, coords, 30), local_knn_residual(panel, coords, 30)),
        }

        for subset, mask in masks.items():
            for kind, (x, y) in residual_sets.items():
                rho, p, n = spearman_block(x[mask], y[mask])
                per_sample_rows.append(
                    {
                        "sample_id": sid,
                        "stage": stage,
                        "subset": subset,
                        "kind": kind,
                        "rho_mapk_panel": rho,
                        "rho_mapk_dm1_like": -rho if np.isfinite(rho) else np.nan,
                        "p": p,
                        "n": n,
                        "n_MAPK_genes": n_mapk,
                        "n_panel_genes": n_panel,
                    }
                )

        full_cov = covariate_matrix(obs, coords, "full")
        epi_mask = masks["epithelial_enriched_top50"]
        if stage != "PT":
            mapk_resid = {gene: residualize(get_expr(adata, gene), full_cov) for gene in MAPK_OUT}
            panel_resid = {gene: residualize(get_expr(adata, gene), full_cov) for gene in PANEL_8}
            for mapk_gene, x in mapk_resid.items():
                for panel_gene, y in panel_resid.items():
                    rho, p, n = spearman_block(x[epi_mask], y[epi_mask])
                    gene_rows.append(
                        {
                            "sample_id": sid,
                            "stage": stage,
                            "mapk_gene": mapk_gene,
                            "panel_gene": panel_gene,
                            "rho_full_resid_epi50": rho,
                            "p": p,
                            "n": n,
                        }
                    )

        pooled_chunks.append(
            pd.DataFrame(
                {
                    "sample_id": sid,
                    "stage": stage,
                    "spot_id": obs.index.astype(str).to_numpy(),
                    "epithelial_enriched_top50": masks["epithelial_enriched_top50"].astype(int),
                    "MAPK_z": mapk,
                    "panel8_z": panel,
                    "MAPK_basic_resid": residual_sets["basic_resid"][0],
                    "panel8_basic_resid": residual_sets["basic_resid"][1],
                    "MAPK_full_resid": residual_sets["full_resid"][0],
                    "panel8_full_resid": residual_sets["full_resid"][1],
                    "MAPK_spatial_poly_resid": residual_sets["spatial_poly_resid"][0],
                    "panel8_spatial_poly_resid": residual_sets["spatial_poly_resid"][1],
                    "MAPK_local_knn30_resid": residual_sets["local_knn30_resid"][0],
                    "panel8_local_knn30_resid": residual_sets["local_knn30_resid"][1],
                    "Epithelial_score": obs["Epithelial_score"].astype(float).to_numpy(),
                    "log1p_total_counts": np.log1p(obs["total_counts"].astype(float).to_numpy()),
                }
            )
        )

    per_sample = pd.DataFrame(per_sample_rows)
    spot_df = pd.concat(pooled_chunks, ignore_index=True)
    gene_pairs = pd.DataFrame(gene_rows)
    gene_summary = (
        gene_pairs.groupby(["mapk_gene", "panel_gene"])["rho_full_resid_epi50"]
        .agg(["median", "mean", "count"])
        .reset_index()
        .rename(columns={"median": "median_rho_full_resid_epi50", "mean": "mean_rho_full_resid_epi50", "count": "n_samples"})
        .sort_values("median_rho_full_resid_epi50")
    )

    pooled_rows = []
    pooled_defs = {
        "all_spots": spot_df,
        "tumor_all_spots": spot_df[spot_df["stage"] != "PT"],
        "tumor_epithelial_enriched_top50": spot_df[(spot_df["stage"] != "PT") & (spot_df["epithelial_enriched_top50"] == 1)],
        "ATC_all_spots": spot_df[spot_df["stage"] == "ATC"],
        "ATC_epithelial_enriched_top50": spot_df[(spot_df["stage"] == "ATC") & (spot_df["epithelial_enriched_top50"] == 1)],
        "PTC_LPTC_epithelial_enriched_top50": spot_df[
            (spot_df["stage"].isin(["PTC", "LPTC"])) & (spot_df["epithelial_enriched_top50"] == 1)
        ],
    }
    kind_cols = {
        "raw": ("MAPK_z", "panel8_z"),
        "basic_resid": ("MAPK_basic_resid", "panel8_basic_resid"),
        "full_resid": ("MAPK_full_resid", "panel8_full_resid"),
        "spatial_poly_resid": ("MAPK_spatial_poly_resid", "panel8_spatial_poly_resid"),
        "local_knn30_resid": ("MAPK_local_knn30_resid", "panel8_local_knn30_resid"),
    }
    for subset, sub in pooled_defs.items():
        for kind in KIND_ORDER:
            xcol, ycol = kind_cols[kind]
            rho, p, n = spearman_block(sub[xcol].to_numpy(), sub[ycol].to_numpy())
            pooled_rows.append(
                {
                    "subset": subset,
                    "kind": kind,
                    "rho_mapk_panel": rho,
                    "rho_mapk_dm1_like": -rho if np.isfinite(rho) else np.nan,
                    "p": p,
                    "n": n,
                }
            )
    pooled = pd.DataFrame(pooled_rows)

    stage_median = (
        per_sample.groupby(["stage", "subset", "kind"], as_index=False)
        .agg(
            median_rho_mapk_panel=("rho_mapk_panel", "median"),
            mean_rho_mapk_panel=("rho_mapk_panel", "mean"),
            n_samples=("sample_id", "nunique"),
        )
        .sort_values(["stage", "subset", "kind"])
    )

    per_sample.to_csv(OUT / "v16_spatial_adjustment_ladder_per_sample.tsv", sep="\t", index=False)
    pooled.to_csv(OUT / "v16_spatial_pooled_adjustment_summary.tsv", sep="\t", index=False)
    stage_median.to_csv(OUT / "v16_spatial_stage_median_ladder.tsv", sep="\t", index=False)
    gene_pairs.to_csv(OUT / "v16_spatial_gene_pair_residual_correlations.tsv", sep="\t", index=False)
    gene_summary.to_csv(OUT / "v16_spatial_gene_pair_residual_summary.tsv", sep="\t", index=False)

    tumor_epi = pooled[(pooled["subset"] == "tumor_epithelial_enriched_top50")].set_index("kind")
    atc_epi = pooled[(pooled["subset"] == "ATC_epithelial_enriched_top50")].set_index("kind")
    sample_epi = per_sample[(per_sample["stage"] != "PT") & (per_sample["subset"] == "epithelial_enriched_top50")]
    metrics = {
        "n_samples": int(meta.shape[0]),
        "n_tumor_samples": int((meta["stage_inferred"] != "PT").sum()),
        "n_spots_total": int(spot_df.shape[0]),
        "tumor_epi_raw_rho": float(tumor_epi.loc["raw", "rho_mapk_panel"]),
        "tumor_epi_full_resid_rho": float(tumor_epi.loc["full_resid", "rho_mapk_panel"]),
        "tumor_epi_spatial_poly_resid_rho": float(tumor_epi.loc["spatial_poly_resid", "rho_mapk_panel"]),
        "tumor_epi_raw_to_full_attenuation_pct": float(
            100
            * (1 - abs(tumor_epi.loc["full_resid", "rho_mapk_panel"]) / abs(tumor_epi.loc["raw", "rho_mapk_panel"]))
        ),
        "ATC_epi_raw_rho": float(atc_epi.loc["raw", "rho_mapk_panel"]),
        "ATC_epi_full_resid_rho": float(atc_epi.loc["full_resid", "rho_mapk_panel"]),
        "tumor_per_sample_median_raw_epi50": float(sample_epi[sample_epi["kind"] == "raw"]["rho_mapk_panel"].median()),
        "tumor_per_sample_median_full_resid_epi50": float(
            sample_epi[sample_epi["kind"] == "full_resid"]["rho_mapk_panel"].median()
        ),
        "tumor_samples_full_resid_le_zero_epi50": int(
            (sample_epi[sample_epi["kind"] == "full_resid"]["rho_mapk_panel"] <= 0).sum()
        ),
        "top_negative_gene_pair": str(
            gene_summary.iloc[0]["mapk_gene"] + "_vs_" + gene_summary.iloc[0]["panel_gene"]
        ),
        "top_negative_gene_pair_median_rho": float(gene_summary.iloc[0]["median_rho_full_resid_epi50"]),
    }
    (OUT / "v16_spatial_rescue_summary.json").write_text(json.dumps(metrics, indent=2))
    return per_sample, pooled, stage_median, gene_summary, metrics


def plot_v16(per_sample: pd.DataFrame, pooled: pd.DataFrame, stage_median: pd.DataFrame, gene_summary: pd.DataFrame, metrics: dict) -> None:
    print("[v16] plotting rescue/autopsy figure")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    ax1, ax2, ax3, ax4, ax5, ax6 = axes.ravel()

    colors = {
        "raw": "#b2182b",
        "basic_resid": "#ef8a62",
        "full_resid": "#f6c85f",
        "spatial_poly_resid": "#67a9cf",
        "local_knn30_resid": "#2166ac",
    }

    # A. Tumor epithelial per-sample ladder.
    sub = per_sample[(per_sample["stage"] != "PT") & (per_sample["subset"] == "epithelial_enriched_top50")]
    x = np.arange(len(KIND_ORDER))
    for sid, s2 in sub.groupby("sample_id"):
        vals = [s2.loc[s2["kind"] == kind, "rho_mapk_panel"].iloc[0] for kind in KIND_ORDER]
        ax1.plot(x, vals, color="#444", alpha=0.28, lw=1)
    med = [sub.loc[sub["kind"] == kind, "rho_mapk_panel"].median() for kind in KIND_ORDER]
    ax1.plot(x, med, color="#111", lw=3, marker="o", ms=7)
    ax1.axhline(0, color="#666", lw=0.8, ls=":")
    ax1.set_xticks(x)
    ax1.set_xticklabels([KIND_LABELS[k] for k in KIND_ORDER], rotation=25, ha="right")
    ax1.set_ylabel("Per-sample Spearman rho\nMAPK output vs Panel-8")
    ax1.set_title("A. Tumor epithelial spots: positive signal collapses after covariate adjustment", loc="left", weight="bold")
    ax1.grid(axis="y", alpha=0.25)

    # B. Pooled ladder.
    pooled_plot = pooled[pooled["subset"].isin(["tumor_epithelial_enriched_top50", "ATC_epithelial_enriched_top50"])].copy()
    width = 0.36
    for offset, subset, label in [
        (-width / 2, "tumor_epithelial_enriched_top50", "Tumor epi-enriched"),
        (width / 2, "ATC_epithelial_enriched_top50", "ATC epi-enriched"),
    ]:
        vals = [pooled_plot[(pooled_plot["subset"] == subset) & (pooled_plot["kind"] == k)]["rho_mapk_panel"].iloc[0] for k in KIND_ORDER]
        ax2.bar(x + offset, vals, width=width, label=label, color=[colors[k] for k in KIND_ORDER], edgecolor="#111", alpha=0.9)
    ax2.axhline(0, color="#666", lw=0.8, ls=":")
    ax2.set_xticks(x)
    ax2.set_xticklabels([KIND_LABELS[k] for k in KIND_ORDER], rotation=25, ha="right")
    ax2.set_ylabel("Pooled spot Spearman rho")
    ax2.set_title("B. Pooled rescue: ATC becomes null after full adjustment", loc="left", weight="bold")
    ax2.legend(frameon=False, fontsize=9)
    ax2.grid(axis="y", alpha=0.25)

    # C. Raw vs full residual per sample.
    raw = sub[sub["kind"] == "raw"].set_index("sample_id")
    full = sub[sub["kind"] == "full_resid"].set_index("sample_id")
    sample_stage = raw["stage"].to_dict()
    stage_color = {"PTC": "#2c7fb8", "LPTC": "#e08214", "ATC": "#b2182b"}
    common = raw.index.intersection(full.index)
    for sid in common:
        ax3.scatter(raw.loc[sid, "rho_mapk_panel"], full.loc[sid, "rho_mapk_panel"], s=85, color=stage_color[sample_stage[sid]], edgecolor="#111", alpha=0.9)
    lim = (-0.15, 0.55)
    ax3.plot(lim, lim, color="#777", lw=0.8, ls=":")
    ax3.axhline(0, color="#666", lw=0.8, ls=":")
    ax3.axvline(0, color="#666", lw=0.8, ls=":")
    ax3.set_xlim(lim)
    ax3.set_ylim(lim)
    ax3.set_xlabel("Raw rho")
    ax3.set_ylabel("Full-covariate residual rho")
    ax3.set_title("C. Sample-level neutralization after full adjustment", loc="left", weight="bold")
    for st, col in stage_color.items():
        ax3.scatter([], [], color=col, edgecolor="#111", label=st)
    ax3.legend(frameon=False, fontsize=9)
    ax3.grid(alpha=0.25)

    # D. Gene-pair heatmap.
    heat = gene_summary.pivot(index="mapk_gene", columns="panel_gene", values="median_rho_full_resid_epi50").reindex(MAPK_OUT)[PANEL_8]
    im = ax4.imshow(heat, cmap="RdBu_r", vmin=-0.12, vmax=0.12, aspect="auto")
    ax4.set_xticks(np.arange(len(PANEL_8)))
    ax4.set_xticklabels(PANEL_8, rotation=45, ha="right")
    ax4.set_yticks(np.arange(len(MAPK_OUT)))
    ax4.set_yticklabels(MAPK_OUT)
    ax4.set_title("D. Gene-pair footholds after full covariate adjustment", loc="left", weight="bold")
    cbar = fig.colorbar(im, ax=ax4, fraction=0.046, pad=0.04)
    cbar.set_label("Median rho across tumor samples")

    # E. Top negative gene-pair bars.
    top_neg = gene_summary.sort_values("median_rho_full_resid_epi50").head(12).iloc[::-1]
    labels = [f"{r.mapk_gene} vs {r.panel_gene}" for _, r in top_neg.iterrows()]
    ax5.barh(np.arange(len(top_neg)), top_neg["median_rho_full_resid_epi50"], color="#2166ac", edgecolor="#111", alpha=0.9)
    ax5.axvline(0, color="#666", lw=0.8, ls=":")
    ax5.set_yticks(np.arange(len(top_neg)))
    ax5.set_yticklabels(labels, fontsize=9)
    ax5.set_xlabel("Median residual Spearman rho")
    ax5.set_title("E. Weak gene-level anti-correlation footholds", loc="left", weight="bold")
    ax5.grid(axis="x", alpha=0.25)

    # F. Stage medians after full adjustment.
    st_sub = stage_median[
        (stage_median["subset"] == "epithelial_enriched_top50")
        & (stage_median["kind"].isin(["raw", "full_resid", "spatial_poly_resid"]))
    ].copy()
    stages = ["PTC", "LPTC", "ATC"]
    x2 = np.arange(len(stages))
    width2 = 0.25
    for i, kind in enumerate(["raw", "full_resid", "spatial_poly_resid"]):
        vals = [
            st_sub[(st_sub["stage"] == st) & (st_sub["kind"] == kind)]["median_rho_mapk_panel"].iloc[0]
            if len(st_sub[(st_sub["stage"] == st) & (st_sub["kind"] == kind)])
            else np.nan
            for st in stages
        ]
        ax6.bar(x2 + (i - 1) * width2, vals, width=width2, color=colors[kind], edgecolor="#111", label=KIND_LABELS[kind])
    ax6.axhline(0, color="#666", lw=0.8, ls=":")
    ax6.set_xticks(x2)
    ax6.set_xticklabels(stages)
    ax6.set_ylabel("Median per-sample rho")
    ax6.set_title("F. Stage view: strongest rescue in ATC/PTC, LPTC remains compartment-heavy", loc="left", weight="bold")
    ax6.legend(frameon=False, fontsize=9)
    ax6.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "Figure SX_v16. Spatial MAPK-panel failure autopsy: the Visium module signal is mostly compartment/QC/spatial structure, not a usable anti-correlation test.",
        fontsize=13,
        weight="bold",
        y=1.01,
    )
    fig.tight_layout()
    png = OUT / "Fig_SX_v16_spatial_failure_autopsy.png"
    pdf = OUT / "Fig_SX_v16_spatial_failure_autopsy.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    ASSETS.mkdir(parents=True, exist_ok=True)
    LIVE_ASSETS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(png, ASSETS / png.name)
    shutil.copy2(pdf, ASSETS / pdf.name)
    shutil.copy2(png, LIVE_ASSETS / png.name)
    shutil.copy2(pdf, LIVE_ASSETS / pdf.name)


def fmt_num(x: float, nd: int = 3) -> str:
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
            if isinstance(val, float):
                text = f"{val:.4g}"
            else:
                text = str(val)
            cells.append(f"<td>{html.escape(text)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table><thead><tr>{headers}</tr></thead><tbody>{''.join(rows)}</tbody></table>"


def build_html(per_sample: pd.DataFrame, pooled: pd.DataFrame, gene_summary: pd.DataFrame, metrics: dict) -> None:
    print("[v16] building/deploying HTML dossier")
    pooled_focus = pooled[pooled["subset"].isin(["tumor_epithelial_enriched_top50", "ATC_epithelial_enriched_top50"])].copy()
    pooled_focus["kind"] = pd.Categorical(pooled_focus["kind"], KIND_ORDER, ordered=True)
    pooled_focus = pooled_focus.sort_values(["subset", "kind"])[["subset", "kind", "rho_mapk_panel", "rho_mapk_dm1_like", "p", "n"]]
    top_gene = gene_summary.sort_values("median_rho_full_resid_epi50").head(12)[
        ["mapk_gene", "panel_gene", "median_rho_full_resid_epi50", "mean_rho_full_resid_epi50", "n_samples"]
    ]
    sample_focus = (
        per_sample[(per_sample["stage"] != "PT") & (per_sample["subset"] == "epithelial_enriched_top50")]
        .pivot_table(index=["sample_id", "stage"], columns="kind", values="rho_mapk_panel")
        .reset_index()
        [["sample_id", "stage", "raw", "basic_resid", "full_resid", "spatial_poly_resid", "local_knn30_resid"]]
        .sort_values(["stage", "sample_id"])
    )

    css = """
    :root{--bg:#0d1117;--panel:#151b23;--line:#30363d;--text:#e6edf3;--muted:#9ba7b4;--gold:#f6c85f;--red:#ef8a62;--blue:#67a9cf}
    *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--text);font-family:"JetBrains Mono",ui-monospace,Menlo,monospace;line-height:1.55}
    .hero{padding:56px 7vw 32px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#101824,#0d1117)}
    .kicker{color:var(--gold);text-transform:uppercase;letter-spacing:.08em;font-size:12px;font-weight:700}
    h1{font-family:"Cormorant Garamond",Georgia,serif;font-size:44px;line-height:1.05;margin:10px 0 14px;font-weight:700}
    .lead{max-width:1040px;color:#c8d1dc;font-size:16px}.crumbs{margin-top:18px;color:var(--muted);font-size:12px}.crumbs a{color:var(--blue)}
    .stats{display:grid;grid-template-columns:repeat(6,minmax(120px,1fr));gap:10px;margin-top:24px}.stat{border:1px solid var(--line);background:#111821;padding:12px;border-radius:6px}
    .stat b{display:block;color:var(--gold);font-size:20px}.stat span{font-size:11px;color:var(--muted)}
    .wrap{display:grid;grid-template-columns:260px 1fr;gap:28px;padding:26px 7vw 60px}.toc{position:sticky;top:16px;align-self:start;border:1px solid var(--line);background:var(--panel);border-radius:6px;padding:14px}
    .toc a{display:block;color:#c8d1dc;text-decoration:none;font-size:12px;padding:6px 0;border-bottom:1px solid rgba(255,255,255,.04)}
    section{margin-bottom:28px;border-top:1px solid var(--line);padding-top:20px}.num{color:var(--gold);font-weight:700;margin-right:8px} h2{font-size:22px;margin:0 0 12px}
    .box{border:1px solid var(--line);background:var(--panel);border-radius:6px;padding:16px;margin:14px 0}.warn{border-color:#7c4b17;background:#21170c}.good{border-color:#345d40;background:#101d16}
    .cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}.card{border:1px solid var(--line);background:#111821;border-radius:6px;padding:13px;text-decoration:none;color:var(--text)}.card b{display:block;color:var(--gold);margin-bottom:4px}.card span{font-size:12px;color:var(--muted)}
    table{width:100%;border-collapse:collapse;font-size:12px;margin:12px 0} th,td{border:1px solid var(--line);padding:7px 8px;text-align:left} th{color:var(--gold);background:#111821}
    img{max-width:100%;border:1px solid var(--line);border-radius:6px;background:#fff}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}.small{font-size:12px;color:var(--muted)}
    code{color:#ffd98a}@media(max-width:900px){.wrap{grid-template-columns:1fr}.toc{position:relative}.stats{grid-template-columns:repeat(2,1fr)}h1{font-size:34px}.grid,.cards{grid-template-columns:1fr}}
    """

    html_text = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Paper 1 v16 spatial failure rescue</title><style>{css}</style></head>
<body>
<div class="hero">
  <div class="kicker">Paper 1 · Fig 8 reserve · v16</div>
  <h1>Spatial MAPK-panel failure rescue dossier</h1>
  <p class="lead">The v15 GSE250521 Visium anti-correlation test failed in the literal module-level form. v16 turns that failure into a controlled reviewer-defense result: the positive spot-level signal largely collapses after epithelial, QC, microenvironment, and spatial-gradient adjustment, so the dataset should not be used as direct mechanism support but also should not be treated as a fatal contradiction.</p>
  <div class="stats">
    <div class="stat"><b>{fmt_num(metrics['tumor_epi_raw_rho'])}</b><span>tumor epithelial raw rho</span></div>
    <div class="stat"><b>{fmt_num(metrics['tumor_epi_full_resid_rho'])}</b><span>tumor epithelial full-resid rho</span></div>
    <div class="stat"><b>{fmt_num(metrics['tumor_epi_raw_to_full_attenuation_pct'],1)}%</b><span>raw-to-full attenuation</span></div>
    <div class="stat"><b>{fmt_num(metrics['ATC_epi_full_resid_rho'])}</b><span>ATC epithelial full-resid rho</span></div>
    <div class="stat"><b>{metrics['tumor_samples_full_resid_le_zero_epi50']}/{metrics['n_tumor_samples']}</b><span>tumor samples ≤0 after full residual</span></div>
    <div class="stat"><b>{fmt_num(metrics['top_negative_gene_pair_median_rho'])}</b><span>{html.escape(metrics['top_negative_gene_pair'])}</span></div>
  </div>
  <div class="crumbs"><a href="index.html">hub</a> / <a href="paper1_methods_reproducibility_dossier.html">methods dossier</a> / v16 spatial rescue</div>
</div>
<div class="wrap">
  <nav class="toc">
    <a href="#tldr">TL;DR</a><a href="#figure">Figure</a><a href="#ladder">Adjustment ladder</a><a href="#gene">Gene-level footholds</a><a href="#decision">Decision</a><a href="#siblings">Sibling pages</a><a href="#sources">Sources</a>
  </nav>
  <main>
    <section id="tldr"><h2><span class="num">01</span>TL;DR</h2>
      <div class="box warn"><b>Not a rescued positive mechanism figure.</b> The expected module-level spatial anti-correlation is still not observed. The usable win is different: v16 shows that the contrary positive result is mostly explained by compartment/QC/spatial structure and becomes near-null under stricter controls.</div>
      <div class="box good"><b>Paper use:</b> keep GSE250521 as a caveat/autopsy reserve. Do not cite it as Fig 8 mechanism support. Use v13/v14 for the mechanism layer and v16 only if a reviewer asks why spatial MAPK output did not anti-correlate with the 8-gene panel in Visium.</div>
    </section>
    <section id="figure"><h2><span class="num">02</span>Composite Figure</h2>
      <img src="assets/paper1/Fig_SX_v16_spatial_failure_autopsy.png" alt="v16 spatial failure autopsy figure">
    </section>
    <section id="ladder"><h2><span class="num">03</span>Adjustment Ladder</h2>
      <p class="small">Values are Spearman correlations between MAPK output module score and Panel-8 score. Negative would match the original anti-correlation hypothesis; near-zero after adjustment neutralizes the apparent contradiction.</p>
      {table_html(pooled_focus, max_rows=20)}
      <h3>Per-sample tumor epithelial-enriched ladder</h3>
      {table_html(sample_focus, max_rows=20)}
    </section>
    <section id="gene"><h2><span class="num">04</span>Gene-Level Footholds</h2>
      <p>The module-level effect is not rescued. A few residualized gene pairs show weak anti-correlations, led by CCND1/DUSP6 versus TPO/SLC5A5. This is a mechanistic hint only, not a replacement for the failed module-level test.</p>
      {table_html(top_gene, max_rows=12)}
    </section>
    <section id="decision"><h2><span class="num">05</span>Decision Matrix</h2>
      <table><thead><tr><th>Question</th><th>Result</th><th>Disposition</th></tr></thead><tbody>
        <tr><td>Does GSE250521 support direct spatial MAPK-output × Panel-8 anti-correlation?</td><td>No. Raw tumor epithelial-enriched rho = {fmt_num(metrics['tumor_epi_raw_rho'])}.</td><td>Do not use as support.</td></tr>
        <tr><td>Is the contrary positive result robust to compartment/QC/microenvironment adjustment?</td><td>Mostly no. Full-residual rho = {fmt_num(metrics['tumor_epi_full_resid_rho'])}; attenuation = {fmt_num(metrics['tumor_epi_raw_to_full_attenuation_pct'],1)}%.</td><td>Use as reviewer-defense autopsy.</td></tr>
        <tr><td>Does ATC retain a contradiction after strict adjustment?</td><td>No. ATC epithelial full-resid rho = {fmt_num(metrics['ATC_epi_full_resid_rho'])}.</td><td>Neutralized.</td></tr>
        <tr><td>Is there a small mechanistic foothold?</td><td>Weak gene-pair anti-correlations, top = {html.escape(metrics['top_negative_gene_pair'])} rho {fmt_num(metrics['top_negative_gene_pair_median_rho'])}.</td><td>Reserve only.</td></tr>
      </tbody></table>
    </section>
    <section id="siblings"><h2><span class="num">06</span>Sibling Pages</h2>
      <div class="cards">
        <a class="card" href="paper1_fig8_mechanism_dossier.html"><b>Fig 8 mechanism dossier</b><span>Primary positive mechanism layer: v6 to v15 forensic trail.</span></a>
        <a class="card" href="paper1_reviewer_defense_dashboard.html"><b>Reviewer defense dashboard</b><span>Q-numbered reviewer-facing evidence map.</span></a>
        <a class="card" href="paper1_methods_reproducibility_dossier.html"><b>Methods reproducibility</b><span>Gene lists, formulas, cohort filters, and v15 reserve outputs.</span></a>
      </div>
    </section>
    <section id="sources"><h2><span class="num">07</span>Sources + Paths</h2>
      <div class="box small">
        <p><code>project/data/processed/GSE250521/*/*.scored.h5ad</code></p>
        <p><code>project/results/p_deconv_2026_05_08/run_v16_spatial_failure_rescue.py</code></p>
        <p><code>v16_spatial_adjustment_ladder_per_sample.tsv</code>, <code>v16_spatial_pooled_adjustment_summary.tsv</code>, <code>v16_spatial_gene_pair_residual_summary.tsv</code>, <code>v16_spatial_rescue_summary.json</code></p>
      </div>
    </section>
  </main>
</div></body></html>
"""
    out_html = HUB / "paper1_spatial_failure_rescue_v16.html"
    out_html.write_text(html_text)
    LIVE.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out_html, LIVE / out_html.name)


def update_index() -> None:
    index = HUB / "index.html"
    if not index.exists():
        return
    text = index.read_text()
    link = '<a class="hero-action" href="paper1_spatial_failure_rescue_v16.html">★ Paper 1 · Spatial failure rescue (v16)</a>'
    if "paper1_spatial_failure_rescue_v16.html" not in text:
        marker = '<a class="hero-action" href="paper1_methods_reproducibility_dossier.html">★ Paper 1 · Methods reproducibility dossier (v15)</a>'
        if marker in text:
            text = text.replace(marker, marker + "\n          " + link)
        else:
            text = text.replace("</div>\n      </section>", f"{link}\n        </div>\n      </section>", 1)
        index.write_text(text)
        if LIVE.exists():
            shutil.copy2(index, LIVE / "index.html")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    per_sample, pooled, stage_median, gene_summary, metrics = run_spatial_rescue()
    plot_v16(per_sample, pooled, stage_median, gene_summary, metrics)
    build_html(per_sample, pooled, gene_summary, metrics)
    update_index()
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
