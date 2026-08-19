#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal, spearmanr
from sklearn.neighbors import NearestNeighbors

from pilot_utils import fdr_bh, pilot_root, setup_logging


CONDITION_ORDER = {
    "normal": 0,
    "PTC": 1,
    "locally_advanced_PTC": 2,
    "ATC": 3,
}

KEY_NICHE_ORDER = [
    "RAI-restorable niche",
    "RAI-low dedifferentiated niche",
    "HLA-visible inflamed niche",
    "HLA-low invisible tumor niche",
    "CD8-excluded / myeloid-CAF niche",
    "APC-rich niche",
    "Drug-delivery failure proxy niche",
    "Mixed/Other",
]


def safe_ratio(a: float, b: float) -> float:
    return float(a / b) if np.isfinite(a) and np.isfinite(b) and b > 0 else np.nan


def neighbor_indices(coords: np.ndarray, k: int) -> np.ndarray:
    k_eff = min(k + 1, len(coords))
    nn = NearestNeighbors(n_neighbors=k_eff, algorithm="auto")
    nn.fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)
    return idx[:, 1:] if idx.shape[1] > 1 else np.empty((len(coords), 0), dtype=int)


def adjacency_for_slide(sdf: pd.DataFrame, k: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    labels = sdf["niche_label"].fillna("Mixed/Other").astype(str).to_numpy()
    coords = sdf[["x", "y"]].to_numpy(dtype=float)
    neigh = neighbor_indices(coords, k)
    label_values = [x for x in KEY_NICHE_ORDER if x in set(labels)] + sorted(set(labels) - set(KEY_NICHE_ORDER))
    label_to_i = {v: i for i, v in enumerate(label_values)}
    codes = np.array([label_to_i[x] for x in labels], dtype=int)
    counts = np.zeros((len(label_values), len(label_values)), dtype=int)
    for source_idx in range(neigh.shape[0]):
        src = codes[source_idx]
        tgts = codes[neigh[source_idx]]
        counts[src] += np.bincount(tgts, minlength=len(label_values))

    target_freq = pd.Series(labels).value_counts(normalize=True).to_dict()
    rows = []
    for i, source in enumerate(label_values):
        n_edges = counts[i].sum()
        if n_edges == 0:
            continue
        for j, target in enumerate(label_values):
            obs = int(counts[i, j])
            obs_frac = obs / n_edges
            exp_frac = target_freq.get(target, 0.0)
            var = n_edges * exp_frac * (1 - exp_frac)
            z_proxy = (obs - n_edges * exp_frac) / np.sqrt(var) if var > 0 else np.nan
            rows.append(
                {
                    "source_niche": source,
                    "neighbor_niche": target,
                    "observed_edges": obs,
                    "source_edges": int(n_edges),
                    "observed_neighbor_fraction": obs_frac,
                    "expected_slide_label_fraction": exp_frac,
                    "log2_enrichment_proxy": np.log2((obs_frac + 1e-6) / (exp_frac + 1e-6)),
                    "binomial_z_proxy_not_independent": z_proxy,
                    "claim_boundary": "Proxy adjacency enrichment; spots are not independent biological replicates and p-values are not claimed.",
                }
            )
    adj = pd.DataFrame(rows)

    def qmask(col: str, q: float, direction: str = "high") -> np.ndarray:
        vals = pd.to_numeric(sdf[col], errors="coerce")
        thr = vals.quantile(q)
        return (vals >= thr).to_numpy() if direction == "high" else (vals <= thr).to_numpy()

    tumor_high = qmask("tumor_epithelial_score", 0.75, "high")
    cyt_high = qmask("cytotoxic_t_score", 0.75, "high")
    cyt_low = qmask("cytotoxic_t_score", 0.25, "low")
    hla_low = qmask("hla_i_apm_score", 0.25, "low")
    rai_low = qmask("rai_differentiation_score", 0.25, "low")
    mapk_high = qmask("mapk_score", 0.75, "high")
    barrier_high = qmask("spatial_myeloid_caf_barrier_score", 0.75, "high")
    delivery_high = qmask("spatial_drug_delivery_failure_proxy", 0.75, "high")
    aggressive_high = qmask("spatial_aggressive_dedifferentiation_score", 0.75, "high")

    def neigh_frac(source_mask: np.ndarray, target_mask: np.ndarray) -> float:
        source_idx = np.where(source_mask)[0]
        if len(source_idx) == 0 or neigh.shape[1] == 0:
            return np.nan
        return float(target_mask[neigh[source_idx]].mean())

    barrier_low = ~barrier_high
    metrics = pd.DataFrame(
        [
            {
                "sample_id": sdf["sample_id"].iloc[0],
                "condition": sdf["condition"].iloc[0],
                "n_spots": len(sdf),
                "n_tumor_high_spots": int(tumor_high.sum()),
                "slide_fraction_cd8_high": float(cyt_high.mean()),
                "tumor_high_neighbor_cd8_high_fraction": neigh_frac(tumor_high, cyt_high),
                "tumor_high_neighbor_cd8_access_ratio": safe_ratio(neigh_frac(tumor_high, cyt_high), float(cyt_high.mean())),
                "tumor_barrier_high_neighbor_cd8_high_fraction": neigh_frac(tumor_high & barrier_high, cyt_high),
                "tumor_barrier_low_neighbor_cd8_high_fraction": neigh_frac(tumor_high & barrier_low, cyt_high),
                "barrier_high_minus_low_cd8_neighbor_fraction": neigh_frac(tumor_high & barrier_high, cyt_high)
                - neigh_frac(tumor_high & barrier_low, cyt_high),
                "hla_low_tumor_neighbor_barrier_high_fraction": neigh_frac(tumor_high & hla_low, barrier_high),
                "hla_low_tumor_neighbor_delivery_high_fraction": neigh_frac(tumor_high & hla_low, delivery_high),
                "hla_low_tumor_neighbor_cd8_high_fraction": neigh_frac(tumor_high & hla_low, cyt_high),
                "rai_low_tumor_neighbor_mapk_high_fraction": neigh_frac(tumor_high & rai_low, mapk_high),
                "rai_low_tumor_neighbor_aggressive_high_fraction": neigh_frac(tumor_high & rai_low, aggressive_high),
                "rai_low_tumor_neighbor_barrier_high_fraction": neigh_frac(tumor_high & rai_low, barrier_high),
                "delivery_high_neighbor_cd8_low_fraction": neigh_frac(delivery_high, cyt_low),
                "delivery_high_neighbor_barrier_high_fraction": neigh_frac(delivery_high, barrier_high),
                "claim_boundary": "Slide-level summary of spot-neighborhood proxies; use for ROI hypothesis generation, not biological replicate p-values.",
            }
        ]
    )
    return adj, metrics


def condition_trends(slide: pd.DataFrame, metrics: pd.DataFrame) -> pd.DataFrame:
    merged = slide.merge(metrics.drop(columns=["condition", "n_spots", "claim_boundary"], errors="ignore"), on="sample_id", how="left")
    merged["condition_order"] = merged["condition"].map(CONDITION_ORDER)
    numeric_cols = [
        c
        for c in merged.columns
        if (
            c.startswith("fraction_")
            or c.startswith("mean_")
            or c.endswith("_fraction")
            or c.endswith("_ratio")
            or c.endswith("_z")
            or c == "barrier_high_minus_low_cd8_neighbor_fraction"
        )
    ]
    rows = []
    for col in numeric_cols:
        vals = pd.to_numeric(merged[col], errors="coerce")
        frame = pd.DataFrame({"condition": merged["condition"], "order": merged["condition_order"], "value": vals}).dropna()
        if frame["condition"].nunique() < 2:
            continue
        groups = [g["value"].to_numpy() for _, g in frame.groupby("condition") if len(g) > 0]
        try:
            kw_p = kruskal(*groups).pvalue if len(groups) >= 2 else np.nan
        except ValueError:
            kw_p = np.nan
        rho, sp_p = (np.nan, np.nan)
        if frame["order"].nunique() > 1 and frame["value"].nunique() > 1:
            rho, sp_p = spearmanr(frame["order"], frame["value"])
        rows.append(
            {
                "metric": col,
                "n_slides": len(frame),
                "n_conditions": frame["condition"].nunique(),
                "spearman_condition_order_rho": rho,
                "spearman_p": sp_p,
                "kruskal_p": kw_p,
                "normal_mean": frame.loc[frame["condition"].eq("normal"), "value"].mean(),
                "PTC_mean": frame.loc[frame["condition"].eq("PTC"), "value"].mean(),
                "locally_advanced_PTC_mean": frame.loc[frame["condition"].eq("locally_advanced_PTC"), "value"].mean(),
                "ATC_mean": frame.loc[frame["condition"].eq("ATC"), "value"].mean(),
                "claim_boundary": "Exploratory slide-level trend only; n=4 slides per condition and no clinical-response claim.",
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["spearman_fdr_q"] = fdr_bh(out["spearman_p"])
        out["kruskal_fdr_q"] = fdr_bh(out["kruskal_p"])
    return out.sort_values("spearman_p", na_position="last") if not out.empty else out


def plot_adjacency(adj: pd.DataFrame, fig_dir: Path) -> None:
    med = adj.groupby(["source_niche", "neighbor_niche"])["log2_enrichment_proxy"].median().reset_index()
    mat = med.pivot(index="source_niche", columns="neighbor_niche", values="log2_enrichment_proxy")
    rows = [x for x in KEY_NICHE_ORDER if x in mat.index]
    cols = [x for x in KEY_NICHE_ORDER if x in mat.columns]
    mat = mat.loc[rows, cols]
    fig, ax = plt.subplots(figsize=(11.5, 9.2))
    sns.heatmap(mat, cmap="vlag", center=0, vmin=-1.5, vmax=1.5, linewidths=0.35, linecolor="#222", ax=ax)
    ax.set_title("Spatial niche adjacency enrichment across GSE250521 slides", fontsize=14, weight="bold")
    ax.set_xlabel("Neighbor niche")
    ax.set_ylabel("Source niche")
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "spatial_niche_adjacency_enrichment_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "spatial_niche_adjacency_enrichment_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_condition(slide: pd.DataFrame, metrics: pd.DataFrame, fig_dir: Path) -> None:
    merged = slide.merge(metrics.drop(columns=["condition", "n_spots", "claim_boundary"], errors="ignore"), on="sample_id", how="left")
    condition_order = ["normal", "PTC", "locally_advanced_PTC", "ATC"]
    key_metrics = [
        "fraction_RAI-low dedifferentiated niche",
        "fraction_CD8-excluded / myeloid-CAF niche",
        "fraction_Drug-delivery failure proxy niche",
        "same_niche_z",
        "tumor_high_neighbor_cd8_access_ratio",
        "barrier_high_minus_low_cd8_neighbor_fraction",
        "hla_low_tumor_neighbor_barrier_high_fraction",
        "rai_low_tumor_neighbor_aggressive_high_fraction",
    ]
    plot_cols = [c for c in key_metrics if c in merged.columns]
    long = merged.melt(id_vars=["sample_id", "condition"], value_vars=plot_cols, var_name="metric", value_name="value")
    n = len(plot_cols)
    fig, axes = plt.subplots(int(np.ceil(n / 2)), 2, figsize=(15, max(7, 3.1 * np.ceil(n / 2))), squeeze=False)
    for ax, metric in zip(axes.ravel(), plot_cols):
        sub = long[long["metric"].eq(metric)]
        sns.boxplot(data=sub, x="condition", y="value", order=condition_order, color="#8dd3c7", fliersize=0, ax=ax)
        sns.stripplot(data=sub, x="condition", y="value", order=condition_order, color="#111111", size=4, ax=ax)
        ax.set_title(metric, fontsize=10, weight="bold")
        ax.tick_params(axis="x", rotation=30)
        ax.set_xlabel("")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle("Slide-level spatial niche and adjacency trends (exploratory)", fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(fig_dir / "spatial_condition_and_adjacency_trends.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "spatial_condition_and_adjacency_trends.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=pilot_root())
    parser.add_argument("--k", type=int, default=6)
    args = parser.parse_args()
    logger = setup_logging("16_spatial_adjacency_condition_trends")
    out = args.root / "results" / "extra_analyses"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    spot = pd.read_csv(args.root / "results" / "tables" / "spatial_spot_vulnerability_scores.tsv", sep="\t")
    slide = pd.read_csv(args.root / "results" / "tables" / "spatial_slide_niche_summary.tsv", sep="\t")
    adj_all = []
    metrics_all = []
    for sample_id, sdf in spot.groupby("sample_id", sort=False):
        logger.info("Spatial adjacency: %s (%d spots)", sample_id, len(sdf))
        adj, metrics = adjacency_for_slide(sdf.reset_index(drop=True), args.k)
        adj.insert(0, "sample_id", sample_id)
        adj.insert(1, "condition", sdf["condition"].iloc[0])
        adj_all.append(adj)
        metrics_all.append(metrics)
    adj_df = pd.concat(adj_all, ignore_index=True)
    metrics_df = pd.concat(metrics_all, ignore_index=True)
    trend_df = condition_trends(slide, metrics_df)

    adj_df.to_csv(tables / "spatial_niche_adjacency_enrichment.tsv", sep="\t", index=False)
    metrics_df.to_csv(tables / "spatial_functional_adjacency_metrics.tsv", sep="\t", index=False)
    trend_df.to_csv(tables / "spatial_condition_trend_summary.tsv", sep="\t", index=False)
    plot_adjacency(adj_df, figs)
    plot_condition(slide, metrics_df, figs)

    logger.info("Wrote %d adjacency rows, %d slide metrics, %d condition trends.", len(adj_df), len(metrics_df), len(trend_df))


if __name__ == "__main__":
    main()
