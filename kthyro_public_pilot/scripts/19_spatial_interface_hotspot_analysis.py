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
from scipy.stats import hypergeom
from sklearn.neighbors import NearestNeighbors

from pilot_utils import fdr_bh, pilot_root, setup_logging


PAIRS = [
    ("tumor_high", "cd8_high", "tumor_to_cd8_high"),
    ("hla_low_tumor", "cd8_high", "hla_low_tumor_to_cd8_high"),
    ("hla_low_tumor", "barrier_high", "hla_low_tumor_to_barrier_high"),
    ("rai_low_tumor", "delivery_high", "rai_low_tumor_to_delivery_failure"),
    ("rai_low_tumor", "aggressive_high", "rai_low_tumor_to_aggressive"),
    ("delivery_high", "vascular_high", "delivery_failure_to_vascular_high"),
    ("delivery_high", "cd8_low", "delivery_failure_to_cd8_low"),
    ("barrier_high", "cd8_high", "barrier_to_cd8_high"),
]

HOTSPOTS = [
    ("hla_low_tumor", "barrier_high", "HLA-low tumor + barrier-high"),
    ("hla_low_tumor", "cd8_low", "HLA-low tumor + CD8-low"),
    ("rai_low_tumor", "delivery_high", "RAI-low tumor + delivery-failure"),
    ("rai_low_tumor", "aggressive_high", "RAI-low tumor + aggressive-high"),
    ("delivery_high", "hypoxia_high", "delivery-failure + hypoxia-high"),
    ("barrier_high", "cd8_low", "barrier-high + CD8-low"),
    ("hla_high", "cd8_high", "HLA-high + CD8-high"),
    ("hla_high", "ifn_high", "HLA-high + IFN-high"),
]


def mask(df: pd.DataFrame, col: str, q: float, high: bool = True) -> np.ndarray:
    vals = pd.to_numeric(df[col], errors="coerce")
    thr = vals.quantile(q)
    return (vals >= thr).to_numpy() if high else (vals <= thr).to_numpy()


def masks(df: pd.DataFrame) -> dict[str, np.ndarray]:
    tumor = mask(df, "tumor_epithelial_score", 0.75, True)
    hla_low = mask(df, "hla_i_apm_score", 0.25, False)
    rai_low = mask(df, "rai_differentiation_score", 0.25, False)
    return {
        "tumor_high": tumor,
        "cd8_high": mask(df, "cytotoxic_t_score", 0.75, True),
        "cd8_low": mask(df, "cytotoxic_t_score", 0.25, False),
        "hla_high": mask(df, "hla_i_apm_score", 0.75, True),
        "hla_low_tumor": hla_low & tumor,
        "rai_low_tumor": rai_low & tumor,
        "barrier_high": mask(df, "spatial_myeloid_caf_barrier_score", 0.75, True),
        "delivery_high": mask(df, "spatial_drug_delivery_failure_proxy", 0.75, True),
        "aggressive_high": mask(df, "spatial_aggressive_dedifferentiation_score", 0.75, True),
        "hypoxia_high": mask(df, "hypoxia_score", 0.75, True),
        "vascular_high": mask(df, "vascular_delivery_proxy_score", 0.75, True),
        "ifn_high": mask(df, "ifn_activation_score", 0.75, True),
    }


def nearest_mean(coords: np.ndarray, source: np.ndarray, target: np.ndarray) -> float:
    src = np.where(source)[0]
    tgt = np.where(target)[0]
    if len(src) == 0 or len(tgt) == 0:
        return np.nan
    nn = NearestNeighbors(n_neighbors=1).fit(coords[tgt])
    return float(nn.kneighbors(coords[src], return_distance=True)[0][:, 0].mean())


def slide_interface(df: pd.DataFrame, n_perm: int, rng: np.random.Generator) -> pd.DataFrame:
    coords = df[["x", "y"]].to_numpy(float)
    mm = masks(df)
    all_idx = np.arange(len(df))
    rows = []
    for src_key, tgt_key, metric in PAIRS:
        src = mm[src_key]
        tgt = mm[tgt_key]
        obs = nearest_mean(coords, src, tgt)
        n_tgt = int(tgt.sum())
        null = []
        if src.sum() > 0 and n_tgt > 0:
            for _ in range(n_perm):
                rt = np.zeros(len(df), bool)
                rt[rng.choice(all_idx, n_tgt, replace=False)] = True
                null.append(nearest_mean(coords, src, rt))
        null = np.array(null, float)
        null = null[np.isfinite(null)]
        null_mean = float(null.mean()) if len(null) else np.nan
        null_sd = float(null.std(ddof=1)) if len(null) > 1 else np.nan
        z = (obs - null_mean) / null_sd if np.isfinite(obs) and np.isfinite(null_sd) and null_sd > 0 else np.nan
        p_close = (np.sum(null <= obs) + 1) / (len(null) + 1) if len(null) else np.nan
        rows.append(
            {
                "sample_id": df["sample_id"].iloc[0],
                "condition": df["condition"].iloc[0],
                "metric": metric,
                "source_mask": src_key,
                "target_mask": tgt_key,
                "n_source": int(src.sum()),
                "n_target": n_tgt,
                "observed_mean_nearest_distance_px": obs,
                "null_mean_distance_px": null_mean,
                "distance_delta_observed_minus_null": obs - null_mean if np.isfinite(obs) and np.isfinite(null_mean) else np.nan,
                "distance_z_observed_minus_null": z,
                "empirical_p_closer_than_random": p_close,
                "claim_boundary": "Nearest-distance proxy using spot coordinates; use as ROI/interface hypothesis, not mechanistic exclusion proof.",
            }
        )
    return pd.DataFrame(rows)


def slide_hotspots(df: pd.DataFrame) -> pd.DataFrame:
    mm = masks(df)
    n = len(df)
    rows = []
    for a_key, b_key, label in HOTSPOTS:
        a = mm[a_key]
        b = mm[b_key]
        ov = int((a & b).sum())
        na = int(a.sum())
        nb = int(b.sum())
        exp = na * nb / n if n else np.nan
        p = hypergeom.sf(ov - 1, n, na, nb) if na and nb else np.nan
        rows.append(
            {
                "sample_id": df["sample_id"].iloc[0],
                "condition": df["condition"].iloc[0],
                "hotspot": label,
                "mask_a": a_key,
                "mask_b": b_key,
                "n_spots": n,
                "n_a": na,
                "n_b": nb,
                "n_overlap": ov,
                "expected_overlap_independence": exp,
                "overlap_fraction_of_slide": ov / n if n else np.nan,
                "jaccard": ov / int((a | b).sum()) if int((a | b).sum()) else np.nan,
                "log2_overlap_enrichment": np.log2((ov + 0.5) / (exp + 0.5)) if np.isfinite(exp) else np.nan,
                "hypergeom_p_enrichment_spot_level_proxy": p,
                "claim_boundary": "Spot-level co-localization proxy; spatial autocorrelation means p-values are descriptive, not biological replicate inference.",
            }
        )
    return pd.DataFrame(rows)


def plot(dist: pd.DataFrame, hot: pd.DataFrame, fig_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 7))
    sns.barplot(data=dist, x="distance_z_observed_minus_null", y="metric", hue="condition", errorbar=None, ax=ax)
    ax.axvline(0, color="#333", lw=1)
    ax.set_xlabel("Z: observed distance minus random distance (negative = closer)")
    ax.set_ylabel("")
    ax.set_title("Spatial interface distance between therapeutic states", weight="bold")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "spatial_interface_distance_z_by_condition.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "spatial_interface_distance_z_by_condition.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

    med = hot.groupby(["hotspot", "condition"])["log2_overlap_enrichment"].median().reset_index()
    mat = med.pivot(index="hotspot", columns="condition", values="log2_overlap_enrichment")
    cols = [c for c in ["normal", "PTC", "locally_advanced_PTC", "ATC"] if c in mat.columns]
    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    sns.heatmap(mat[cols], cmap="vlag", center=0, linewidths=0.35, linecolor="#222", ax=ax)
    ax.set_title("Spatial hotspot co-localization enrichment", weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(fig_dir / "spatial_hotspot_colocalization_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(fig_dir / "spatial_hotspot_colocalization_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=pilot_root())
    parser.add_argument("--n-perm", type=int, default=75)
    parser.add_argument("--seed", type=int, default=20260509)
    args = parser.parse_args()
    logger = setup_logging("19_spatial_interface_hotspot_analysis")
    out = args.root / "results" / "extra_analyses"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    spot = pd.read_csv(args.root / "results" / "tables" / "spatial_spot_vulnerability_scores.tsv", sep="\t")
    rng = np.random.default_rng(args.seed)
    dist_all = []
    hot_all = []
    for sid, sdf in spot.groupby("sample_id", sort=False):
        logger.info("Interface/hotspot %s", sid)
        sdf = sdf.reset_index(drop=True)
        dist_all.append(slide_interface(sdf, args.n_perm, rng))
        hot_all.append(slide_hotspots(sdf))
    dist = pd.concat(dist_all, ignore_index=True)
    hot = pd.concat(hot_all, ignore_index=True)
    dist["empirical_fdr_q_closer"] = fdr_bh(dist["empirical_p_closer_than_random"])
    hot["hypergeom_fdr_q_spot_level_proxy"] = fdr_bh(hot["hypergeom_p_enrichment_spot_level_proxy"])
    dist.to_csv(tables / "spatial_interface_distance_tests.tsv", sep="\t", index=False)
    hot.to_csv(tables / "spatial_hotspot_colocalization_tests.tsv", sep="\t", index=False)
    plot(dist, hot, figs)
    logger.info("Wrote %d interface rows and %d hotspot rows.", len(dist), len(hot))


if __name__ == "__main__":
    main()
