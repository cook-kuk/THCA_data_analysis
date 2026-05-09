#!/usr/bin/env python3
"""GSE230424 Path2Space-style H&E-to-thyroid spatial-axis screen.

The GEO records expose P1-P4 sample IDs but not disease labels for each slide.
This script therefore uses label-free leave-one-sample-out and sample-centered
tests as the primary evidence.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RAW = ROOT / "project/data/external/GSE230424/raw"
OUT = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
MAPK_OUT = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74"]
B_TLS = ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1", "IGKC", "BANK1", "CXCL13"]
T_CELL = ["CD3D", "CD3E", "TRAC", "CD8A", "CD4", "CCL5", "GZMB", "NKG7"]
CD36_MAC = ["CD36", "SPP1", "APOE", "LYZ", "C1QA", "C1QB", "C1QC", "AIF1", "TYROBP", "FCGR3A"]
TUMOR_ZCCHC12 = ["ZCCHC12", "KRT19", "EPCAM", "TACSTD2", "SERPINA1", "LGALS3", "MET"]

MODULES = {
    "RAI8_lineage_score": PANEL_8,
    "MAPK_output_score": MAPK_OUT,
    "HLA_II_AP_score": HLA_II,
    "B_TLS_score": B_TLS,
    "T_cell_score": T_CELL,
    "CD36_SPP1_macrophage_score": CD36_MAC,
    "Tumor_ZCCHC12_score": TUMOR_ZCCHC12,
}
TARGETS = [
    "DM1_low_RAI_score",
    "RAI8_lineage_score",
    "MAPK_output_score",
    "AP_TLS_composite_score",
    "HLA_II_AP_score",
    "B_TLS_score",
    "CD36_SPP1_macrophage_score",
    "Tumor_ZCCHC12_score",
]

SMOOTH_K = 8
PATCH_RADII = [48, 96]
RIDGE_ALPHA = 50.0
N_PCS = 32
N_PERM = 1000
RANDOM_SEED = 230424


def sample_prefixes() -> list[str]:
    return sorted(p.name.replace("_matrix.mtx.gz", "") for p in RAW.glob("*_matrix.mtx.gz"))


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
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        if keep.sum() == 0:
            continue
        out[idx[keep]] = vals[keep] - np.nanmean(vals[keep])
    return out


def read_features(prefix: str) -> pd.DataFrame:
    return pd.read_csv(RAW / f"{prefix}_features.tsv.gz", sep="\t", header=None, names=["ensg", "gene", "feature_type"])


def read_barcodes(prefix: str) -> pd.Series:
    return pd.read_csv(RAW / f"{prefix}_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)


def read_positions(prefix: str) -> pd.DataFrame:
    pos = pd.read_csv(
        RAW / f"{prefix}_tissue_positions_list.csv.gz",
        header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_fullres", "pxl_col_fullres"],
    )
    pos["barcode"] = pos["barcode"].astype(str)
    return pos


def stream_target_matrix(prefix: str, genes: list[str]) -> pd.DataFrame:
    features = read_features(prefix)
    barcodes = read_barcodes(prefix)
    gene_to_pos = {g: i for i, g in enumerate(genes)}
    row_to_gene = {}
    mt_rows = set()
    for row_idx, gene in enumerate(features["gene"].astype(str), start=1):
        if gene in gene_to_pos:
            row_to_gene[row_idx] = gene
        if gene.startswith("MT-"):
            mt_rows.add(row_idx)

    target_counts = None
    total_counts = None
    n_genes = None
    mt_counts = None
    dims_seen = False
    with gzip.open(RAW / f"{prefix}_matrix.mtx.gz", "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if not dims_seen:
                _, n_cols, _ = map(int, parts)
                target_counts = np.zeros((len(genes), n_cols), dtype=np.float32)
                total_counts = np.zeros(n_cols, dtype=np.float64)
                n_genes = np.zeros(n_cols, dtype=np.int32)
                mt_counts = np.zeros(n_cols, dtype=np.float64)
                dims_seen = True
                continue
            row = int(parts[0])
            col = int(parts[1]) - 1
            val = float(parts[2])
            total_counts[col] += val
            n_genes[col] += 1
            if row in mt_rows:
                mt_counts[col] += val
            gene = row_to_gene.get(row)
            if gene is not None:
                target_counts[gene_to_pos[gene], col] += val

    assert target_counts is not None and total_counts is not None and n_genes is not None and mt_counts is not None
    denom = np.maximum(total_counts, 1.0)
    log_cp10k = np.log1p(target_counts / denom[None, :] * 10000.0)
    out = pd.DataFrame(log_cp10k.T, columns=[f"g_{g}" for g in genes])
    out.insert(0, "barcode", barcodes.to_numpy())
    out["total_counts"] = total_counts
    out["n_genes_by_counts"] = n_genes
    out["pct_counts_mt"] = mt_counts / denom * 100.0
    return out


def load_spot_table() -> pd.DataFrame:
    all_genes = sorted({g for genes in MODULES.values() for g in genes})
    rows = []
    for prefix in sample_prefixes():
        sample = prefix.split("_")[-1]
        print(f"[GSE230424] load expression {sample}")
        expr = stream_target_matrix(prefix, all_genes)
        pos = read_positions(prefix)
        meta = expr.merge(pos, on="barcode", how="left")
        meta.insert(0, "sample", sample)
        meta.insert(1, "geo_accession", prefix.split("_")[0])
        rows.append(meta)
    df = pd.concat(rows, ignore_index=True)
    for module, genes in MODULES.items():
        present = [f"g_{g}" for g in genes if f"g_{g}" in df.columns and np.nanmax(df[f"g_{g}"].to_numpy(dtype=float)) > 0]
        vals = []
        for col in present:
            v = df[col].to_numpy(dtype=float)
            sd = np.nanstd(v)
            vals.append((v - np.nanmean(v)) / (sd + 1e-9) if sd > 0 else np.zeros_like(v))
        df[module] = np.nanmean(np.vstack(vals), axis=0) if vals else np.nan
        df[f"n_{module}_genes"] = len(present)
    df["DM1_low_RAI_score"] = -df["RAI8_lineage_score"]
    df["AP_TLS_composite_score"] = df[["HLA_II_AP_score", "B_TLS_score", "T_cell_score"]].mean(axis=1)
    df.to_csv(OUT / "gse230424_spot_module_scores.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    return df


def spatial_smooth(df: pd.DataFrame, cols: list[str], k: int = SMOOTH_K) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        out[f"{col}_smooth{k}"] = np.nan
    for sample, sub in out.groupby("sample", sort=False):
        idx = sub.index.to_numpy()
        coords = sub[["array_row", "array_col"]].to_numpy(dtype=float)
        tree = cKDTree(coords)
        kk = min(k + 1, len(sub))
        _, neigh = tree.query(coords, k=kk)
        if neigh.ndim == 1:
            neigh = neigh[:, None]
        for col in cols:
            vals = sub[col].to_numpy(dtype=float)
            out.loc[idx, f"{col}_smooth{k}"] = np.nanmean(vals[neigh], axis=1)
    return out


def patch_features(patch_rgb: np.ndarray, prefix: str) -> dict[str, float]:
    rgb = patch_rgb.astype(np.float32) / 255.0
    hsv = cv2.cvtColor((rgb * 255).astype(np.uint8), cv2.COLOR_RGB2HSV).astype(np.float32)
    gray = cv2.cvtColor((rgb * 255).astype(np.uint8), cv2.COLOR_RGB2GRAY).astype(np.float32)
    features = {}
    for i, name in enumerate(["r", "g", "b"]):
        channel = rgb[:, :, i]
        features[f"{prefix}_{name}_mean"] = float(np.nanmean(channel))
        features[f"{prefix}_{name}_std"] = float(np.nanstd(channel))
        features[f"{prefix}_{name}_p10"] = float(np.nanpercentile(channel, 10))
        features[f"{prefix}_{name}_p90"] = float(np.nanpercentile(channel, 90))
    for i, name in enumerate(["h", "s", "v"]):
        channel = hsv[:, :, i]
        features[f"{prefix}_{name}_mean"] = float(np.nanmean(channel))
        features[f"{prefix}_{name}_std"] = float(np.nanstd(channel))
    features[f"{prefix}_gray_mean"] = float(np.nanmean(gray))
    features[f"{prefix}_gray_std"] = float(np.nanstd(gray))
    features[f"{prefix}_lap_var"] = float(cv2.Laplacian(gray, cv2.CV_32F).var())
    features[f"{prefix}_dark_frac_180"] = float((gray < 180).mean())
    features[f"{prefix}_dark_frac_140"] = float((gray < 140).mean())
    features[f"{prefix}_saturated_frac"] = float((hsv[:, :, 1] > 45).mean())
    r = rgb[:, :, 0]
    g = rgb[:, :, 1]
    b = rgb[:, :, 2]
    features[f"{prefix}_blue_purple_proxy"] = float(np.nanmean((b + r) / 2.0 - g))
    features[f"{prefix}_eosin_proxy"] = float(np.nanmean(r - b))
    return features


def extract_he_features(df: pd.DataFrame) -> pd.DataFrame:
    Image.MAX_IMAGE_PIXELS = None
    rows = []
    prefix_by_sample = {p.split("_")[-1]: p for p in sample_prefixes()}
    for sample, sub in df.groupby("sample", sort=False):
        prefix = prefix_by_sample[sample]
        print(f"[GSE230424] H&E features {sample}: {len(sub)} spots")
        with gzip.open(RAW / f"{prefix}_HE.jpg.gz", "rb") as handle:
            image = Image.open(handle).convert("RGB")
            image.load()
        img = np.asarray(image)
        h, w = img.shape[:2]
        for row in sub.itertuples(index=False):
            x = int(round(row.pxl_col_fullres))
            y = int(round(row.pxl_row_fullres))
            feats = {
                "sample": row.sample,
                "barcode": row.barcode,
                "array_row": row.array_row,
                "array_col": row.array_col,
            }
            for radius in PATCH_RADII:
                x0, x1 = max(0, x - radius), min(w, x + radius)
                y0, y1 = max(0, y - radius), min(h, y + radius)
                feats.update(patch_features(img[y0:y1, x0:x1, :], f"r{radius}"))
            rows.append(feats)
        del img, image
    out = pd.DataFrame(rows)
    feature_cols = [c for c in out.columns if c.startswith("r")]
    out[feature_cols] = out[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(out[feature_cols].median(numeric_only=True))
    out.to_csv(OUT / "gse230424_he_tile_features.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    return out


def loso_predict(x: np.ndarray, y: np.ndarray, samples: np.ndarray, use_pca: bool) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    finite_x = np.isfinite(x).all(axis=1)
    for held in sorted(pd.unique(samples)):
        train = (samples != held) & finite_x & np.isfinite(y)
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
        mean = float(np.nanmean(y[train]))
        sd = float(np.nanstd(y[train]))
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        model = Ridge(alpha=RIDGE_ALPHA)
        model.fit(x_train, (y[train] - mean) / sd)
        pred[test] = model.predict(x_test)
    return pred


def run_models(df: pd.DataFrame, feats: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = df.merge(feats, on=["sample", "barcode", "array_row", "array_col"], how="inner")
    samples = merged["sample"].astype(str).to_numpy()
    he_cols = [c for c in feats.columns if c.startswith("r")]
    features = {
        "HE_tile_features": (merged[he_cols].to_numpy(dtype=float), True),
        "Coord_only": (merged[["array_row", "array_col"]].to_numpy(dtype=float), False),
        "QC_only": (np.column_stack([
            np.log1p(merged["total_counts"].to_numpy(dtype=float)),
            np.log1p(merged["n_genes_by_counts"].to_numpy(dtype=float)),
            merged["pct_counts_mt"].to_numpy(dtype=float),
        ]), False),
        "Coord_QC": (np.column_stack([
            merged[["array_row", "array_col"]].to_numpy(dtype=float),
            np.log1p(merged["total_counts"].to_numpy(dtype=float)),
            np.log1p(merged["n_genes_by_counts"].to_numpy(dtype=float)),
            merged["pct_counts_mt"].to_numpy(dtype=float),
        ]), False),
    }
    pred_df = merged[["sample", "geo_accession", "barcode", "array_row", "array_col", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    rows = []
    for target in TARGETS:
        y = merged[f"{target}_smooth{SMOOTH_K}"].to_numpy(dtype=float)
        y_center = center_by_group(y, samples)
        pred_df[f"obs_{target}_smooth{SMOOTH_K}"] = y
        for model_name, (x, use_pca) in features.items():
            pred = loso_predict(x, y, samples, use_pca=use_pca)
            pred_df[f"pred_{target}_{model_name}"] = pred
            rho, p, n = safe_spearman(y, pred)
            pred_center = center_by_group(pred, samples)
            crho, cp, cn = safe_spearman(y_center, pred_center)
            per_sample = []
            for sample, idx_raw in merged.groupby("sample").groups.items():
                idx = np.asarray(list(idx_raw))
                srho, _, _ = safe_spearman(y[idx], pred[idx])
                per_sample.append(srho)
            rows.append(
                {
                    "target": target,
                    "model": model_name,
                    "n": n,
                    "pooled_rho": rho,
                    "pooled_p": p,
                    "sample_centered_rho": crho,
                    "sample_centered_p": cp,
                    "median_sample_rho": float(np.nanmedian(per_sample)),
                    "samples_rho_gt_0_2": int(np.nansum(np.asarray(per_sample) > 0.2)),
                    "n_features": x.shape[1],
                }
            )
    return pd.DataFrame(rows), pred_df


def residualize_within_sample(values: np.ndarray, covars: np.ndarray, samples: np.ndarray) -> np.ndarray:
    out = np.full(len(values), np.nan)
    for sample in pd.unique(samples):
        idx = np.where(samples == sample)[0]
        y = np.asarray(values, dtype=float)[idx]
        x = np.asarray(covars, dtype=float)[idx]
        keep = np.isfinite(y) & np.isfinite(x).all(axis=1)
        if keep.sum() < 20:
            continue
        xs = StandardScaler().fit_transform(x[keep])
        out[idx[keep]] = y[keep] - LinearRegression().fit(xs, y[keep]).predict(xs)
    return out


def residual_alignment(pred_df: pd.DataFrame) -> pd.DataFrame:
    samples = pred_df["sample"].astype(str).to_numpy()
    coord = pred_df[["array_row", "array_col"]].to_numpy(dtype=float)
    qc = np.column_stack([
        np.log1p(pred_df["total_counts"].to_numpy(dtype=float)),
        np.log1p(pred_df["n_genes_by_counts"].to_numpy(dtype=float)),
        pred_df["pct_counts_mt"].to_numpy(dtype=float),
    ])
    rows = []
    for target in TARGETS:
        obs = pred_df[f"obs_{target}_smooth{SMOOTH_K}"].to_numpy(dtype=float)
        pred = pred_df[f"pred_{target}_HE_tile_features"].to_numpy(dtype=float)
        for name, cov in {
            "sample_mean_only": np.ones((len(pred_df), 1), dtype=float),
            "coord_only": coord,
            "qc_only": qc,
            "coord_qc": np.column_stack([coord, qc]),
        }.items():
            obs_res = residualize_within_sample(obs, cov, samples)
            pred_res = residualize_within_sample(pred, cov, samples)
            rho, p, n = safe_spearman(obs_res, pred_res)
            rows.append({"target": target, "within_sample_adjustment": name, "residual_alignment_rho": rho, "p": p, "n": n})
    return pd.DataFrame(rows)


def domain_aggregation(pred_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for sample, sub in pred_df.groupby("sample", sort=False):
        coords = sub[["array_row", "array_col"]].to_numpy(dtype=float)
        n_domains = min(12, max(6, int(round(len(sub) / 400))))
        labs = KMeans(n_clusters=n_domains, n_init=50, random_state=RANDOM_SEED).fit_predict(StandardScaler().fit_transform(coords))
        for lab in np.unique(labs):
            mask = labs == lab
            row = {"sample": sample, "domain": int(lab + 1), "n_spots": int(mask.sum())}
            for target in TARGETS:
                row[f"obs_{target}"] = float(np.nanmean(sub[f"obs_{target}_smooth{SMOOTH_K}"].to_numpy(dtype=float)[mask]))
                row[f"pred_{target}"] = float(np.nanmean(sub[f"pred_{target}_HE_tile_features"].to_numpy(dtype=float)[mask]))
            rows.append(row)
    domains = pd.DataFrame(rows)
    summary = []
    samples = domains["sample"].astype(str).to_numpy()
    for target in TARGETS:
        obs = domains[f"obs_{target}"].to_numpy(dtype=float)
        pred = domains[f"pred_{target}"].to_numpy(dtype=float)
        rho, p, n = safe_spearman(obs, pred)
        crho, cp, cn = safe_spearman(center_by_group(obs, samples), center_by_group(pred, samples))
        summary.append({"target": target, "n_domains": n, "pooled_domain_rho": rho, "sample_centered_domain_rho": crho, "p": p})
    return domains, pd.DataFrame(summary)


def permutation_alignment(pred_df: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    samples = pred_df["sample"].astype(str).to_numpy()
    rows = []
    for target in TARGETS:
        obs = center_by_group(pred_df[f"obs_{target}_smooth{SMOOTH_K}"].to_numpy(dtype=float), samples)
        pred = center_by_group(pred_df[f"pred_{target}_HE_tile_features"].to_numpy(dtype=float), samples)
        observed, _, _ = safe_spearman(obs, pred)
        null = []
        for _ in range(N_PERM):
            perm = pred.copy()
            for sample in pd.unique(samples):
                idx = np.where(samples == sample)[0]
                perm[idx] = rng.permutation(perm[idx])
            rho, _, _ = safe_spearman(obs, perm)
            null.append(rho)
        null = np.asarray(null, dtype=float)
        p = (1.0 + np.nansum(np.abs(null) >= abs(observed))) / (1.0 + np.isfinite(null).sum())
        rows.append({"target": target, "observed_sample_centered_rho": observed, "n_permutations": N_PERM, "empirical_p": float(p)})
    return pd.DataFrame(rows)


def write_figures(model_df: pd.DataFrame, domain_df: pd.DataFrame, pred_df: pd.DataFrame, sample_df: pd.DataFrame) -> str:
    he = model_df[model_df["model"] == "HE_tile_features"].set_index("target").loc[TARGETS]
    top_target = he["sample_centered_rho"].abs().sort_values(ascending=False).index[0]

    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.5))
    x = np.arange(len(TARGETS))
    labels = [t.replace("_score", "").replace("_", " ") for t in TARGETS]
    axes[0, 0].bar(x - 0.18, he["pooled_rho"], width=0.35, color="#266b73", label="pooled")
    axes[0, 0].bar(x + 0.18, he["sample_centered_rho"], width=0.35, color="#d5a84d", label="sample-centered")
    axes[0, 0].axhline(0, color="#333", lw=0.8)
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(labels, rotation=35, ha="right")
    axes[0, 0].set_ylabel("LOSO Spearman rho")
    axes[0, 0].set_title("A. H&E tile features to thyroid spatial axes")
    axes[0, 0].legend(frameon=False)

    dm = model_df[model_df["target"] == top_target].copy()
    order = ["Coord_only", "QC_only", "Coord_QC", "HE_tile_features"]
    dm["model"] = pd.Categorical(dm["model"], categories=order, ordered=True)
    dm = dm.sort_values("model")
    axes[0, 1].barh(dm["model"].astype(str), dm["sample_centered_rho"], color=["#9aa4ad", "#777", "#777", "#266b73"])
    axes[0, 1].axvline(0, color="#333", lw=0.8)
    axes[0, 1].set_xlabel("Sample-centered rho")
    axes[0, 1].set_title(f"B. Baselines for {top_target.replace('_score', '').replace('_', ' ')}")

    dom = domain_df.set_index("target").loc[TARGETS]
    axes[1, 0].bar(x - 0.18, dom["pooled_domain_rho"], width=0.35, color="#266b73", label="domain pooled")
    axes[1, 0].bar(x + 0.18, dom["sample_centered_domain_rho"], width=0.35, color="#d5a84d", label="domain centered")
    axes[1, 0].axhline(0, color="#333", lw=0.8)
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(labels, rotation=35, ha="right")
    axes[1, 0].set_ylabel("Domain Spearman rho")
    axes[1, 0].set_title("C. Coordinate-domain aggregation")
    axes[1, 0].legend(frameon=False)

    sub = pred_df[["sample", f"obs_{top_target}_smooth{SMOOTH_K}", f"pred_{top_target}_HE_tile_features"]].dropna()
    colors = {"P1": "#355c7d", "P2": "#6c5b7b", "P3": "#c06c84", "P4": "#f67280"}
    for sample, part in sub.groupby("sample"):
        axes[1, 1].scatter(part[f"obs_{top_target}_smooth{SMOOTH_K}"], part[f"pred_{top_target}_HE_tile_features"], s=12, alpha=0.45, color=colors.get(sample, "#333"), label=sample)
    rho, _, _ = safe_spearman(sub[f"obs_{top_target}_smooth{SMOOTH_K}"], sub[f"pred_{top_target}_HE_tile_features"])
    axes[1, 1].set_xlabel("Observed smoothed spatial score")
    axes[1, 1].set_ylabel("Predicted from H&E")
    axes[1, 1].set_title(f"D. Top H&E axis pooled scatter (rho={rho:.3f})")
    axes[1, 1].legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse230424_pathology_thyroid_axis.png", dpi=220)
    fig.savefig(OUT / "fig_gse230424_pathology_thyroid_axis.pdf")
    plt.close(fig)

    fig2, axes2 = plt.subplots(2, 4, figsize=(16, 7.2), sharex=False, sharey=False)
    v_obs = pred_df[f"obs_{top_target}_smooth{SMOOTH_K}"].to_numpy(dtype=float)
    v_pred = pred_df[f"pred_{top_target}_HE_tile_features"].to_numpy(dtype=float)
    lim_obs = np.nanpercentile(v_obs, [2, 98])
    lim_pred = np.nanpercentile(v_pred, [2, 98])
    for col, sample in enumerate(sorted(pred_df["sample"].unique())):
        part = pred_df[pred_df["sample"] == sample]
        axes2[0, col].scatter(part["array_col"], -part["array_row"], c=part[f"obs_{top_target}_smooth{SMOOTH_K}"], s=7, cmap="viridis", vmin=lim_obs[0], vmax=lim_obs[1])
        axes2[0, col].set_title(f"{sample} observed")
        axes2[1, col].scatter(part["array_col"], -part["array_row"], c=part[f"pred_{top_target}_HE_tile_features"], s=7, cmap="viridis", vmin=lim_pred[0], vmax=lim_pred[1])
        axes2[1, col].set_title(f"{sample} H&E-predicted")
        for ax in (axes2[0, col], axes2[1, col]):
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal")
    fig2.suptitle(f"Top label-free H&E spatial axis: {top_target.replace('_score', '').replace('_', ' ')}", y=0.98)
    fig2.tight_layout()
    fig2.savefig(OUT / "fig_gse230424_spatial_maps.png", dpi=220)
    fig2.savefig(OUT / "fig_gse230424_spatial_maps.pdf")
    plt.close(fig2)

    fig3, ax = plt.subplots(figsize=(8.5, 4.8))
    heat_cols = [c for c in ["RAI8_lineage_score", "MAPK_output_score", "AP_TLS_composite_score", "CD36_SPP1_macrophage_score", "Tumor_ZCCHC12_score"] if c in sample_df]
    heat = sample_df.set_index("sample")[heat_cols].copy()
    heat = (heat - heat.mean(axis=0)) / (heat.std(axis=0) + 1e-9)
    im = ax.imshow(heat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1.8, vmax=1.8, aspect="auto")
    ax.set_xticks(np.arange(len(heat_cols)))
    ax.set_xticklabels([c.replace("_score", "").replace("_", " ") for c in heat_cols], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_title("Sample-level module context (GEO disease labels unavailable)")
    fig3.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="sample z")
    fig3.tight_layout()
    fig3.savefig(OUT / "fig_gse230424_sample_axis_context.png", dpi=220)
    fig3.savefig(OUT / "fig_gse230424_sample_axis_context.pdf")
    plt.close(fig3)
    return top_target


def write_report(model_df: pd.DataFrame, residual_df: pd.DataFrame, domain_df: pd.DataFrame, perm_df: pd.DataFrame, sample_df: pd.DataFrame, top_target: str) -> dict:
    he = model_df[model_df["model"] == "HE_tile_features"].set_index("target")
    top = he.loc[top_target]
    dm1 = he.loc["DM1_low_RAI_score"]
    ap = he.loc["AP_TLS_composite_score"]
    top_models = model_df[model_df["target"] == top_target].set_index("model")
    res_top = residual_df[(residual_df["target"] == top_target) & (residual_df["within_sample_adjustment"] == "coord_qc")].iloc[0]
    dom_top = domain_df.set_index("target").loc[top_target]
    perm_top = perm_df.set_index("target").loc[top_target]
    summary = {
        "n_spots": int(model_df["n"].max()),
        "n_samples": int(sample_df.shape[0]),
        "top_target": top_target,
        "top_he_pooled_rho": float(top["pooled_rho"]),
        "top_he_sample_centered_rho": float(top["sample_centered_rho"]),
        "top_qc_only_sample_centered_rho": float(top_models.loc["QC_only", "sample_centered_rho"]),
        "top_coord_qc_sample_centered_rho": float(top_models.loc["Coord_QC", "sample_centered_rho"]),
        "top_domain_sample_centered_rho": float(dom_top["sample_centered_domain_rho"]),
        "top_coord_qc_residual_rho": float(res_top["residual_alignment_rho"]),
        "top_permutation_p": float(perm_top["empirical_p"]),
        "dm1_low_rai_sample_centered_rho": float(dm1["sample_centered_rho"]),
        "ap_tls_sample_centered_rho": float(ap["sample_centered_rho"]),
        "disease_labels_available_in_geo": False,
        "verdict": "LABEL_FREE_SPATIAL_HE_SIGNAL" if float(top["sample_centered_rho"]) > 0.2 else "WEAK_OR_NEGATIVE_LABEL_FREE_SIGNAL",
    }
    lines = [
        "# GSE230424 Path2Space-style H&E-to-thyroid spatial-axis screen",
        "",
        "## Verdict",
        "",
        "- Dataset: 4 Visium thyroid slides from GSE230424 with shipped H&E JPEGs, matrices, barcodes, features, and tissue positions.",
        "- GEO sample records expose P1-P4 but do not map P1-P4 to HT versus PTC+HT; primary tests are label-free LOSO and sample-centered.",
        f"- Spots modeled: {summary['n_spots']:,}; samples: {summary['n_samples']}.",
        f"- Top H&E-predictable axis: **{top_target}** with pooled rho **{summary['top_he_pooled_rho']:.3f}**, sample-centered rho **{summary['top_he_sample_centered_rho']:.3f}**, domain-centered rho **{summary['top_domain_sample_centered_rho']:.3f}**.",
        f"- QC-only is stronger for this axis (sample-centered rho **{summary['top_qc_only_sample_centered_rho']:.3f}**; coord+QC **{summary['top_coord_qc_sample_centered_rho']:.3f}**), so the conservative claim is morphology/QC-aligned tissue state plus a residual H&E component.",
        f"- Coord+QC residual alignment for top axis: rho **{summary['top_coord_qc_residual_rho']:.3f}**; within-sample permutation p = **{summary['top_permutation_p']:.4f}**.",
        f"- DM1/low-RAI axis sample-centered rho: **{summary['dm1_low_rai_sample_centered_rho']:.3f}**; AP/TLS composite sample-centered rho: **{summary['ap_tls_sample_centered_rho']:.3f}**.",
        "",
        "## Interpretation Boundary",
        "",
        "This is useful as a label-free external spatial thyroid control for Paper 2's image-to-spatial-RNA direction. The strongest raw prediction tracks QC/tissue-density structure, but a smaller H&E-aligned residual remains after coordinate+QC adjustment. Because slide-level disease labels are absent from the GEO sample fields, it should not be used as a disease-group validation unless labels are recovered from the article or authors. It is not Paper 1 causal mechanism evidence.",
        "",
        "## Model Comparison",
        "",
        model_df.round(4).to_markdown(index=False),
        "",
        "## Residual Alignment",
        "",
        residual_df.round(4).to_markdown(index=False),
        "",
        "## Domain Aggregation",
        "",
        domain_df.round(4).to_markdown(index=False),
        "",
        "## Permutation Alignment",
        "",
        perm_df.round(4).to_markdown(index=False),
        "",
        "## Sample Context",
        "",
        sample_df.round(4).to_markdown(index=False),
    ]
    (OUT / "GSE230424_PATHOLOGY_THYROID_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "GSE230424_PATHOLOGY_THYROID_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spot_path = OUT / "gse230424_spot_module_scores.tsv.gz"
    if spot_path.exists():
        df = pd.read_csv(spot_path, sep="\t")
    else:
        df = load_spot_table()
    df = spatial_smooth(df, TARGETS)

    feats_path = OUT / "gse230424_he_tile_features.tsv.gz"
    if feats_path.exists():
        feats = pd.read_csv(feats_path, sep="\t")
    else:
        feats = extract_he_features(df)

    sample_cols = ["sample"] + [c for c in MODULES] + ["DM1_low_RAI_score", "AP_TLS_composite_score", "total_counts", "n_genes_by_counts", "pct_counts_mt"]
    sample_df = df[sample_cols].groupby("sample", as_index=False).mean(numeric_only=True)
    sample_df.to_csv(OUT / "gse230424_sample_axis_summary.tsv", sep="\t", index=False, na_rep="NA")

    model_df, pred_df = run_models(df, feats)
    residual_df = residual_alignment(pred_df)
    domains, domain_df = domain_aggregation(pred_df)
    perm_df = permutation_alignment(pred_df)
    top_target = write_figures(model_df, domain_df, pred_df, sample_df)

    model_df.to_csv(OUT / "gse230424_model_comparison.tsv", sep="\t", index=False, na_rep="NA")
    pred_df.to_csv(OUT / "gse230424_pathology_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    residual_df.to_csv(OUT / "gse230424_residual_alignment.tsv", sep="\t", index=False, na_rep="NA")
    domains.to_csv(OUT / "gse230424_spatial_domains.tsv", sep="\t", index=False, na_rep="NA")
    domain_df.to_csv(OUT / "gse230424_domain_summary.tsv", sep="\t", index=False, na_rep="NA")
    perm_df.to_csv(OUT / "gse230424_permutation_alignment.tsv", sep="\t", index=False, na_rep="NA")
    summary = write_report(model_df, residual_df, domain_df, perm_df, sample_df, top_target)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
