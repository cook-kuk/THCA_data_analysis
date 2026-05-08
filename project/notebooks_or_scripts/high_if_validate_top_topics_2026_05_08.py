#!/usr/bin/env python3
"""Validate high-probability high-IF THCA directions.

This script adds a spatial expression-CNV bridge using downloaded GENCODE
gene coordinates and UCSC cytoband arm annotations, then tests the current
high-probability topics H1-H7:

H1/H2: CNV-residual spatial ecosystem / spatial phylogeography
H3: multimodal visibility boundary using existing UNI embeddings
H4-H7: spatial ecotype, CAF/ECM, TROP2, immune/TLS axes
"""

from __future__ import annotations

import gzip
import math
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_high_if_top_topic_validation.md"

GTF = ROOT / "project/data/external_annotations/gencode.v48.annotation.gtf.gz"
CYTOBAND = ROOT / "project/data/external_annotations/hg38_cytoBand.txt.gz"
SAMPLE_META = ROOT / "project/data/processed/GSE250521/sample_metadata.tsv"

JOINT = ROOT / "project/results/03_pathology_poc/spark_st_joint_per_spot.tsv.gz"
TROP2 = ROOT / "project/results/spatial_full_2026_05_06/spatial_per_spot_TROP2.tsv.gz"
TCGA_CNV = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/data_search/cbio_class6_dm2_dm1_armdriver_validation.tsv"
TCGA_ROBUST = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test1_cnv_stratified_robustness.tsv"
TCGA_COV = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/three_more_tests_2026_05_08/test1_cnv_covariate_auc.tsv"
GEOMX_STATS = ROOT / "project/results/spatial_full_2026_05_06/GSE301163_crossplatform_stats.tsv"
GEOMX_ROI = ROOT / "project/results/spatial_full_2026_05_06/GSE301163_per_ROI_scores.tsv"
TCGA_COX = ROOT / "project/results/03_pathology_poc/spark_tcga_cox_modules.tsv"
TCGA_COX_ADJ = ROOT / "project/results/03_pathology_poc/spark_tcga_multivariate_cox_DSS.tsv"
UNI_META = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embed_metadata_size224.tsv"
UNI_NPZ = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embeddings_size224.npz"

T00_GAIN_ARMS = ["7p", "7q", "12q", "16p", "16q"]
T00_LOSS_ARMS = ["2p", "2q"]

GENESETS = {
    "TACSTD2_z": ["TACSTD2"],
    "CAF_ECM_z": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "BGN", "SFRP2", "SFRP4", "FAP", "ACTA2", "VIM"],
    "TLS_B_z": ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "CXCL13", "AICDA", "JCHAIN", "MZB1", "IGHG1", "IGHM"],
    "Tcell_cytotoxic_z": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "NKG7"],
    "Macrophage_TAM_z": ["LST1", "TYROBP", "CD68", "MSR1", "MRC1", "SPP1", "TREM2", "C1QA", "C1QB"],
    "RAI_thyroid_z": ["SLC5A5", "TPO", "TG", "TSHR", "DIO1", "DIO2", "PAX8", "NKX2-1", "FOXE1"],
    "EMT_stress_z": ["VIM", "FN1", "SNAI2", "ZEB1", "ZEB2", "HIF1A", "VEGFA", "LDHA"],
}


def require(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)


def parse_attrs(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, val in re.findall(r'([A-Za-z0-9_]+) "([^"]+)"', attr):
        out[key] = val
    return out


def parse_gencode(gtf: Path) -> pd.DataFrame:
    cache = OUT / "gencode_v48_gene_arm_map.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")

    require(gtf)
    require(CYTOBAND)

    cyto = pd.read_csv(
        CYTOBAND,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "band", "stain"],
    )
    cyto["chrom"] = cyto["chrom"].str.replace("^chr", "", regex=True)
    cyto = cyto[cyto["chrom"].isin([str(i) for i in range(1, 23)] + ["X", "Y"])].copy()
    cyto["arm"] = cyto["chrom"] + cyto["band"].str[0]

    rows = []
    with gzip.open(gtf, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom, _src, feature, start, end, _score, _strand, _frame, attrs = line.rstrip("\n").split("\t")
            if feature != "gene":
                continue
            chrom = chrom.replace("chr", "")
            if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y"}:
                continue
            parsed = parse_attrs(attrs)
            gene_name = parsed.get("gene_name")
            gene_id = parsed.get("gene_id", "").split(".")[0]
            gene_type = parsed.get("gene_type", parsed.get("gene_biotype", ""))
            if not gene_name:
                continue
            start_i = int(start)
            end_i = int(end)
            mid = (start_i + end_i) // 2
            arm_hit = cyto[(cyto["chrom"] == chrom) & (cyto["start"] <= mid) & (cyto["end"] >= mid)]
            if arm_hit.empty:
                continue
            rows.append(
                {
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "gene_type": gene_type,
                    "chrom": chrom,
                    "start": start_i,
                    "end": end_i,
                    "arm": arm_hit.iloc[0]["arm"],
                }
            )

    gene_map = pd.DataFrame(rows)
    gene_map["length"] = gene_map["end"] - gene_map["start"] + 1
    gene_map["protein_coding_priority"] = (gene_map["gene_type"] == "protein_coding").astype(int)
    gene_map = (
        gene_map.sort_values(["gene_name", "protein_coding_priority", "length"], ascending=[True, False, False])
        .drop_duplicates("gene_name", keep="first")
        .drop(columns=["protein_coding_priority", "length"])
        .reset_index(drop=True)
    )
    gene_map.to_csv(cache, sep="\t", index=False)
    return gene_map


def condition_from_sample(sample_id: str) -> str:
    if "_N-" in sample_id:
        return "N"
    if "_PTC-" in sample_id:
        return "PTC"
    if "_LPTC-" in sample_id:
        return "LPTC"
    if "_ATC-" in sample_id:
        return "ATC"
    return "UNK"


def normalize_log1p_sparse(x: sparse.spmatrix) -> sparse.csr_matrix:
    x = x.tocsr().astype(np.float64)
    lib = np.asarray(x.sum(axis=1)).ravel()
    lib[lib <= 0] = 1.0
    x = x.multiply((1e4 / lib)[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x


def normal_reference(samples: pd.DataFrame, selected_symbols: list[str], selected_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sums = np.zeros(len(selected_symbols), dtype=np.float64)
    sums2 = np.zeros(len(selected_symbols), dtype=np.float64)
    n_obs = 0
    normals = samples[samples["condition"] == "N"]
    if normals.empty:
        raise RuntimeError("No normal GSE250521 samples found for expression-CNV reference")
    for row in normals.itertuples(index=False):
        a = ad.read_h5ad(ROOT / row.h5ad)
        x = normalize_log1p_sparse(a.X)
        xs = x[:, selected_idx].tocsr()
        sums += np.asarray(xs.sum(axis=0)).ravel()
        sums2 += np.asarray(xs.power(2).sum(axis=0)).ravel()
        n_obs += xs.shape[0]
    mean = sums / max(n_obs, 1)
    var = np.maximum(sums2 / max(n_obs, 1) - mean**2, 1e-6)
    sd = np.sqrt(var)
    sd[sd < 0.05] = 0.05
    return mean, sd


def build_weight_matrix(selected_symbols: list[str], groups: dict[str, list[str]], mean: np.ndarray, sd: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray, list[str], dict[str, int]]:
    symbol_to_i = {g: i for i, g in enumerate(selected_symbols)}
    row_idx: list[int] = []
    col_idx: list[int] = []
    data: list[float] = []
    offsets: list[float] = []
    names: list[str] = []
    n_used: dict[str, int] = {}
    for col, (name, genes) in enumerate(groups.items()):
        idx = [symbol_to_i[g] for g in genes if g in symbol_to_i]
        if not idx:
            continue
        idx_arr = np.asarray(idx, dtype=int)
        weight = 1.0 / (sd[idx_arr] * len(idx_arr))
        row_idx.extend(idx)
        col_idx.extend([len(names)] * len(idx))
        data.extend(weight.tolist())
        offsets.append(float(np.sum(mean[idx_arr] / sd[idx_arr]) / len(idx_arr)))
        names.append(name)
        n_used[name] = len(idx_arr)
    w = sparse.csr_matrix((data, (row_idx, col_idx)), shape=(len(selected_symbols), len(names)))
    return w, np.asarray(offsets), names, n_used


def score_all_spots() -> pd.DataFrame:
    per_spot_path = OUT / "spatial_expression_cnv_per_spot.tsv.gz"
    if per_spot_path.exists():
        return pd.read_csv(per_spot_path, sep="\t")

    gene_map = parse_gencode(GTF)
    samples = pd.read_csv(SAMPLE_META, sep="\t")
    samples["condition"] = samples["sample_id"].map(condition_from_sample)

    first = ad.read_h5ad(ROOT / samples.iloc[0]["h5ad"], backed="r")
    var_names = pd.Index(first.var_names.astype(str))
    in_map = gene_map[gene_map["gene_name"].isin(var_names)].copy()
    in_map = in_map.sort_values("gene_name").drop_duplicates("gene_name")
    selected_symbols = in_map["gene_name"].tolist()
    selected_idx = var_names.get_indexer(selected_symbols)
    ok = selected_idx >= 0
    selected_symbols = [g for g, keep in zip(selected_symbols, ok) if keep]
    selected_idx = selected_idx[ok]
    in_map = in_map[in_map["gene_name"].isin(selected_symbols)].copy()

    arm_groups = {
        arm: in_map.loc[in_map["arm"] == arm, "gene_name"].tolist()
        for arm in sorted(in_map["arm"].unique(), key=lambda x: (x[:-1].zfill(2), x[-1]))
    }
    arm_groups = {k: v for k, v in arm_groups.items() if len(v) >= 30}

    mean, sd = normal_reference(samples, selected_symbols, selected_idx)
    arm_w, arm_offsets, arm_names, arm_n = build_weight_matrix(selected_symbols, arm_groups, mean, sd)
    gs_w, gs_offsets, gs_names, gs_n = build_weight_matrix(selected_symbols, GENESETS, mean, sd)

    in_map.to_csv(OUT / "spatial_cnv_gene_arm_map_used.tsv", sep="\t", index=False)
    pd.DataFrame({"feature": arm_names, "n_genes": [arm_n[x] for x in arm_names]}).to_csv(OUT / "spatial_cnv_arm_gene_counts.tsv", sep="\t", index=False)
    pd.DataFrame({"feature": gs_names, "n_genes": [gs_n[x] for x in gs_names]}).to_csv(OUT / "spatial_geneset_gene_counts.tsv", sep="\t", index=False)

    chunks = []
    for row in samples.itertuples(index=False):
        a = ad.read_h5ad(ROOT / row.h5ad)
        x = normalize_log1p_sparse(a.X)
        xs = x[:, selected_idx].tocsr()
        arm_scores = xs @ arm_w
        if sparse.issparse(arm_scores):
            arm_scores = arm_scores.toarray()
        else:
            arm_scores = np.asarray(arm_scores)
        arm_scores = arm_scores - arm_offsets[None, :]

        gs_scores = xs @ gs_w
        if sparse.issparse(gs_scores):
            gs_scores = gs_scores.toarray()
        else:
            gs_scores = np.asarray(gs_scores)
        gs_scores = gs_scores - gs_offsets[None, :]

        obs = a.obs.reset_index().rename(columns={"index": "spot_id"})
        obs["sample_id"] = row.sample_id
        obs["condition"] = row.condition
        keep_cols = ["sample_id", "spot_id", "condition", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
        frame = obs[keep_cols].copy()
        for i, name in enumerate(arm_names):
            frame[f"arm_{name}"] = arm_scores[:, i]
        for i, name in enumerate(gs_names):
            frame[name] = gs_scores[:, i]
        gain = [f"arm_{a}" for a in T00_GAIN_ARMS if f"arm_{a}" in frame.columns]
        loss = [f"arm_{a}" for a in T00_LOSS_ARMS if f"arm_{a}" in frame.columns]
        frame["t00_spatial_cnv_signature"] = frame[gain].mean(axis=1) - frame[loss].mean(axis=1)
        frame["t00_gain_component"] = frame[gain].mean(axis=1)
        frame["t00_loss_component"] = frame[loss].mean(axis=1)
        chunks.append(frame)

    per_spot = pd.concat(chunks, ignore_index=True)
    normals = per_spot[per_spot["condition"] == "N"]["t00_spatial_cnv_signature"]
    per_spot["t00_signature_norm_q95_high"] = per_spot["t00_spatial_cnv_signature"] > normals.quantile(0.95)
    per_spot["t00_signature_norm_q99_high"] = per_spot["t00_spatial_cnv_signature"] > normals.quantile(0.99)

    joint = pd.read_csv(JOINT, sep="\t")
    trop2 = pd.read_csv(TROP2, sep="\t")
    trop2 = trop2[["sample_id", "spot_id", "TROP2"]].rename(columns={"TROP2": "TROP2_raw"})
    per_spot = per_spot.merge(joint, on=["sample_id", "spot_id"], how="left", suffixes=("", "_joint"))
    per_spot = per_spot.merge(trop2, on=["sample_id", "spot_id"], how="left")
    per_spot.to_csv(per_spot_path, sep="\t", index=False, compression="gzip")
    return per_spot


def spearman(x, y) -> tuple[float, float, int]:
    tmp = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if tmp.shape[0] < 5 or tmp["x"].nunique() < 2 or tmp["y"].nunique() < 2:
        return np.nan, np.nan, int(tmp.shape[0])
    r, p = stats.spearmanr(tmp["x"], tmp["y"])
    return float(r), float(p), int(tmp.shape[0])


def cohen_d(a, b) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)
    if pooled <= 0 or np.isnan(pooled):
        return np.nan
    return float((a.mean() - b.mean()) / math.sqrt(pooled))


def mw_test(a, b) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def moran_knn(coords: np.ndarray, values: np.ndarray, k: int = 6) -> float:
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values) & np.isfinite(coords).all(axis=1)
    if mask.sum() <= k + 2 or np.nanstd(values[mask]) <= 0:
        return np.nan
    coords = coords[mask]
    z = values[mask] - np.nanmean(values[mask])
    denom = np.sum(z**2)
    if denom <= 0:
        return np.nan
    nn = NearestNeighbors(n_neighbors=min(k + 1, len(z))).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)[:, 1:]
    wij = idx.size
    num = 0.0
    for i in range(idx.shape[0]):
        num += float(np.sum(z[i] * z[idx[i]]))
    return float(len(z) / wij * num / denom)


def run_spatial_validations(per_spot: pd.DataFrame) -> None:
    expected = [
        OUT / "spatial_cnv_condition_summary.tsv",
        OUT / "spatial_cnv_condition_tests.tsv",
        OUT / "spatial_cnv_arm_condition_tests.tsv",
        OUT / "spatial_cnv_axis_correlations.tsv",
        OUT / "spatial_cnv_morans.tsv",
    ]
    if all(p.exists() for p in expected):
        return
    target = "t00_spatial_cnv_signature"
    rows = []
    for condition, sub in per_spot.groupby("condition"):
        rows.append(
            {
                "condition": condition,
                "n_spots": sub.shape[0],
                "n_slides": sub["sample_id"].nunique(),
                "mean_signature": sub[target].mean(),
                "median_signature": sub[target].median(),
                "high_q95_rate": sub["t00_signature_norm_q95_high"].mean(),
                "high_q99_rate": sub["t00_signature_norm_q99_high"].mean(),
                "mean_TROP2_z": sub.get("TACSTD2_z", pd.Series(dtype=float)).mean(),
                "mean_CAF_ECM_z": sub.get("CAF_ECM_z", pd.Series(dtype=float)).mean(),
                "mean_TLS_B_z": sub.get("TLS_B_z", pd.Series(dtype=float)).mean(),
                "mean_RAI_thyroid_z": sub.get("RAI_thyroid_z", pd.Series(dtype=float)).mean(),
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "spatial_cnv_condition_summary.tsv", sep="\t", index=False)

    normal = per_spot[per_spot["condition"] == "N"]
    tests = []
    for condition in ["PTC", "LPTC", "ATC"]:
        sub = per_spot[per_spot["condition"] == condition]
        if sub.empty:
            continue
        table = [
            [int(sub["t00_signature_norm_q95_high"].sum()), int((~sub["t00_signature_norm_q95_high"]).sum())],
            [int(normal["t00_signature_norm_q95_high"].sum()), int((~normal["t00_signature_norm_q95_high"]).sum())],
        ]
        odds, p_fisher = stats.fisher_exact(table)
        tests.append(
            {
                "comparison": f"{condition}_vs_N",
                "n_a": sub.shape[0],
                "n_normal": normal.shape[0],
                "mean_a": sub[target].mean(),
                "mean_normal": normal[target].mean(),
                "cohen_d": cohen_d(sub[target], normal[target]),
                "mw_p": mw_test(sub[target], normal[target]),
                "q95_high_rate_a": sub["t00_signature_norm_q95_high"].mean(),
                "q95_high_rate_normal": normal["t00_signature_norm_q95_high"].mean(),
                "fisher_or": odds,
                "fisher_p": p_fisher,
            }
        )
    pd.DataFrame(tests).to_csv(OUT / "spatial_cnv_condition_tests.tsv", sep="\t", index=False)

    arm_rows = []
    for arm in T00_GAIN_ARMS + T00_LOSS_ARMS:
        col = f"arm_{arm}"
        if col not in per_spot.columns:
            continue
        for condition in ["PTC", "LPTC", "ATC"]:
            sub = per_spot[per_spot["condition"] == condition]
            if sub.empty:
                continue
            direction = "gain_like" if arm in T00_GAIN_ARMS else "loss_like_negative_expected"
            arm_rows.append(
                {
                    "arm": arm,
                    "condition": condition,
                    "direction": direction,
                    "mean_condition": sub[col].mean(),
                    "mean_normal": normal[col].mean(),
                    "delta_vs_normal": sub[col].mean() - normal[col].mean(),
                    "cohen_d": cohen_d(sub[col], normal[col]),
                    "mw_p": mw_test(sub[col], normal[col]),
                }
            )
    pd.DataFrame(arm_rows).to_csv(OUT / "spatial_cnv_arm_condition_tests.tsv", sep="\t", index=False)

    corr_targets = [
        "TACSTD2_z",
        "TROP2_raw",
        "CAF_ECM_z",
        "TLS_B_z",
        "Tcell_cytotoxic_z",
        "Macrophage_TAM_z",
        "RAI_thyroid_z",
        "DM1_like_score",
        "RAI_8_score",
        "Hypoxia_score",
        "EMT_score",
        "n_lym",
        "n_fib",
    ]
    corr_rows = []
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    for y in corr_targets:
        if y not in cancer.columns:
            continue
        r, p, n = spearman(cancer[target], cancer[y])
        corr_rows.append({"scope": "pooled_cancer_spots", "target": y, "spearman_r": r, "p": p, "n": n})
        per_slide_r = []
        for sample_id, sub in cancer.groupby("sample_id"):
            rs, ps, ns = spearman(sub[target], sub[y])
            if np.isfinite(rs):
                per_slide_r.append(rs)
                corr_rows.append({"scope": sample_id, "target": y, "spearman_r": rs, "p": ps, "n": ns})
        if per_slide_r:
            corr_rows.append({"scope": "slide_median_r", "target": y, "spearman_r": float(np.median(per_slide_r)), "p": np.nan, "n": len(per_slide_r)})
    pd.DataFrame(corr_rows).to_csv(OUT / "spatial_cnv_axis_correlations.tsv", sep="\t", index=False)

    moran_rows = []
    moran_cols = [target, "t00_gain_component", "t00_loss_component", "TACSTD2_z", "CAF_ECM_z", "TLS_B_z", "RAI_thyroid_z"]
    for sample_id, sub in per_spot.groupby("sample_id"):
        coords = sub[["array_row", "array_col"]].to_numpy(float)
        for col in moran_cols:
            if col in sub.columns:
                moran_rows.append(
                    {
                        "sample_id": sample_id,
                        "condition": sub["condition"].iloc[0],
                        "axis": col,
                        "moran_i_k6": moran_knn(coords, sub[col].to_numpy(float), k=6),
                        "n_spots": sub.shape[0],
                    }
                )
    pd.DataFrame(moran_rows).to_csv(OUT / "spatial_cnv_morans.tsv", sep="\t", index=False)


def run_ecotype_validation(per_spot: pd.DataFrame) -> None:
    expected = [
        OUT / "spatial_ecotype_summary.tsv",
        OUT / "spatial_ecotype_silhouette.tsv",
        OUT / "spatial_ecotype_per_spot.tsv.gz",
    ]
    if all(p.exists() for p in expected):
        return
    features = [
        "t00_spatial_cnv_signature",
        "TACSTD2_z",
        "CAF_ECM_z",
        "TLS_B_z",
        "Tcell_cytotoxic_z",
        "Macrophage_TAM_z",
        "RAI_thyroid_z",
        "EMT_stress_z",
    ]
    features = [f for f in features if f in per_spot.columns]
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    x = cancer[features].replace([np.inf, -np.inf], np.nan).fillna(cancer[features].median(numeric_only=True))
    xz = StandardScaler().fit_transform(x)
    sil_rows = []
    best_k = 4
    best_sil = -np.inf
    for k in [3, 4, 5, 6]:
        labels = KMeans(n_clusters=k, random_state=8, n_init=20).fit_predict(xz)
        sil = (
            silhouette_score(xz, labels, sample_size=min(4000, xz.shape[0]), random_state=8)
            if len(np.unique(labels)) > 1
            else np.nan
        )
        sil_rows.append({"k": k, "silhouette": sil})
        if np.isfinite(sil) and sil > best_sil:
            best_k, best_sil = k, sil
    labels = KMeans(n_clusters=best_k, random_state=8, n_init=50).fit_predict(xz)
    cancer["spatial_ecotype"] = labels

    ec = (
        cancer.groupby("spatial_ecotype")
        .agg(
            n_spots=("spot_id", "size"),
            n_slides=("sample_id", "nunique"),
            mean_signature=("t00_spatial_cnv_signature", "mean"),
            mean_TROP2=("TACSTD2_z", "mean"),
            mean_CAF=("CAF_ECM_z", "mean"),
            mean_TLS=("TLS_B_z", "mean"),
            mean_Tcell=("Tcell_cytotoxic_z", "mean"),
            mean_TAM=("Macrophage_TAM_z", "mean"),
            mean_RAI=("RAI_thyroid_z", "mean"),
            mean_EMT=("EMT_stress_z", "mean"),
            high_q95_rate=("t00_signature_norm_q95_high", "mean"),
        )
        .reset_index()
    )
    cond_frac = pd.crosstab(cancer["spatial_ecotype"], cancer["condition"], normalize="index").reset_index()
    ec = ec.merge(cond_frac, on="spatial_ecotype", how="left")
    ec["ecotype_label"] = ec.apply(label_ecotype, axis=1)
    ec.to_csv(OUT / "spatial_ecotype_summary.tsv", sep="\t", index=False)
    pd.DataFrame(sil_rows).to_csv(OUT / "spatial_ecotype_silhouette.tsv", sep="\t", index=False)
    cancer[["sample_id", "spot_id", "condition", "spatial_ecotype"] + features].to_csv(OUT / "spatial_ecotype_per_spot.tsv.gz", sep="\t", index=False, compression="gzip")


def label_ecotype(row: pd.Series) -> str:
    vals = {
        "CNV": row.get("mean_signature", -99),
        "TROP2": row.get("mean_TROP2", -99),
        "CAF": row.get("mean_CAF", -99),
        "TLS": row.get("mean_TLS", -99),
        "TAM": row.get("mean_TAM", -99),
        "RAI": row.get("mean_RAI", -99),
    }
    top = sorted(vals.items(), key=lambda x: x[1], reverse=True)[:2]
    return "+".join([x[0] for x in top])


def run_uni_retrieval(per_spot: pd.DataFrame) -> None:
    if not UNI_META.exists() or not UNI_NPZ.exists():
        return
    if (OUT / "uni_retrieval_predictions.tsv.gz").exists() and (OUT / "uni_retrieval_metrics.tsv").exists():
        return
    meta = pd.read_csv(UNI_META, sep="\t").rename(columns={"slide": "sample_id"})
    emb = np.load(UNI_NPZ)["embeddings"]
    if emb.shape[0] != meta.shape[0]:
        raise RuntimeError("UNI embedding rows do not match metadata")

    targets = [
        "t00_spatial_cnv_signature",
        "TACSTD2_z",
        "CAF_ECM_z",
        "TLS_B_z",
        "RAI_thyroid_z",
        "DM1_like_score",
        "RAI_8_score",
    ]
    dat = meta.merge(per_spot[["sample_id", "spot_id", "condition"] + [t for t in targets if t in per_spot.columns]], on=["sample_id", "spot_id"], how="inner")
    keep_idx = dat.index.to_numpy()
    emb = emb[keep_idx]
    targets = [t for t in targets if t in dat.columns]

    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    condition_text = enc.fit_transform(dat[["condition"]])

    pred_rows = []
    metric_rows = []
    modes = ["uni_only", "text_only", "uni_text_a06", "uni_text_a08"]
    global_thresholds = {t: dat[t].quantile(0.75) for t in targets}

    for held in sorted(dat["sample_id"].unique()):
        train = dat["sample_id"] != held
        test = ~train
        if train.sum() < 50 or test.sum() < 10:
            continue

        scaler = StandardScaler(with_mean=True, with_std=True).fit(emb[train])
        x_train = scaler.transform(emb[train])
        x_test = scaler.transform(emb[test])
        n_comp = min(64, x_train.shape[1], max(2, train.sum() - 1))
        pca = PCA(n_components=n_comp, random_state=8).fit(x_train)
        uni_train = pca.transform(x_train)
        uni_test = pca.transform(x_test)

        txt_scaler = StandardScaler(with_mean=False).fit(condition_text[train])
        txt_train = txt_scaler.transform(condition_text[train])
        txt_test = txt_scaler.transform(condition_text[test])

        mode_data = {
            "uni_only": (uni_train, uni_test),
            "text_only": (txt_train, txt_test),
            "uni_text_a06": (np.hstack([0.6 * StandardScaler().fit_transform(uni_train), 0.4 * txt_train]), None),
            "uni_text_a08": (np.hstack([0.8 * StandardScaler().fit_transform(uni_train), 0.2 * txt_train]), None),
        }
        for mode in ["uni_text_a06", "uni_text_a08"]:
            a = 0.6 if mode.endswith("a06") else 0.8
            uni_sc = StandardScaler().fit(uni_train)
            mode_data[mode] = (
                np.hstack([a * uni_sc.transform(uni_train), (1 - a) * txt_train]),
                np.hstack([a * uni_sc.transform(uni_test), (1 - a) * txt_test]),
            )

        for mode, (xtr, xte) in mode_data.items():
            nn = NearestNeighbors(n_neighbors=min(25, train.sum()), metric="euclidean").fit(xtr)
            nbr = nn.kneighbors(xte, return_distance=False)
            train_df = dat.loc[train].reset_index(drop=True)
            test_df = dat.loc[test].reset_index(drop=True)
            for target in targets:
                y_train = train_df[target].to_numpy(float)
                y_test = test_df[target].to_numpy(float)
                y_pred = np.nanmean(y_train[nbr], axis=1)
                for i, pred in enumerate(y_pred):
                    pred_rows.append(
                        {
                            "sample_id": test_df.iloc[i]["sample_id"],
                            "spot_id": test_df.iloc[i]["spot_id"],
                            "condition": test_df.iloc[i]["condition"],
                            "target": target,
                            "mode": mode,
                            "observed": y_test[i],
                            "predicted": pred,
                        }
                    )

    pred = pd.DataFrame(pred_rows)
    pred.to_csv(OUT / "uni_retrieval_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
    for (target, mode), sub in pred.groupby(["target", "mode"]):
        r, p, n = spearman(sub["observed"], sub["predicted"])
        labels = sub["observed"] >= global_thresholds[target]
        if labels.nunique() == 2:
            auc = roc_auc_score(labels.astype(int), sub["predicted"])
        else:
            auc = np.nan
        metric_rows.append({"target": target, "mode": mode, "spearman_r": r, "p": p, "n": n, "auroc_q75": auc})
    pd.DataFrame(metric_rows).sort_values(["target", "spearman_r"], ascending=[True, False]).to_csv(OUT / "uni_retrieval_metrics.tsv", sep="\t", index=False)


def run_external_summary() -> None:
    rows = []
    if TCGA_CNV.exists():
        cbio = pd.read_csv(TCGA_CNV, sep="\t")
        for r in cbio.head(8).itertuples(index=False):
            rows.append({"evidence_layer": "TCGA_cBio_ARMDRIVER", "metric": r.metric, "effect": f"DM2={r.DM2_rate:.3g}, DM1={r.DM1_rate:.3g}, OR={r.OR:.3g}", "p": r.p})
    if TCGA_ROBUST.exists():
        rob = pd.read_csv(TCGA_ROBUST, sep="\t")
        for _, r in rob.head(8).iterrows():
            rows.append(
                {
                    "evidence_layer": "TCGA_stratified_robustness",
                    "metric": r.get("test", r.get("subset", "")),
                    "effect": f"DM2={r.get('DM2_rate', np.nan):.3g}, DM1={r.get('DM1_rate', np.nan):.3g}, OR={r.get('OR')}",
                    "p": r.get("p"),
                }
            )
    if TCGA_COV.exists():
        cov = pd.read_csv(TCGA_COV, sep="\t")
        if not cov.empty:
            rows.append(
                {
                    "evidence_layer": "TCGA_covariate_AUC",
                    "metric": "mean_cross_validated_delta",
                    "effect": (
                        f"covariates={cov['covariates_only_auc'].mean():.3g}, "
                        f"covariates+CNV={cov['covariates_plus_signature_auc'].mean():.3g}, "
                        f"delta_AUC={cov['delta_auc'].mean():.3g}"
                    ),
                    "p": np.nan,
                }
            )
    if GEOMX_STATS.exists():
        gs = pd.read_csv(GEOMX_STATS, sep="\t")
        for r in gs.head(10).itertuples(index=False):
            rows.append({"evidence_layer": "GSE301163_GeoMx", "metric": r.comparison, "effect": f"mean_a={r.panck_mean:.3g}, mean_b={r.vim_mean:.3g}", "p": r.MW_p})
    if TCGA_COX.exists():
        cox = pd.read_csv(TCGA_COX, sep="\t")
        for r in cox.itertuples(index=False):
            rows.append({"evidence_layer": "TCGA_DSS_Cox_univariate", "metric": r.module, "effect": f"HR={r.HR:.3g}", "p": r.p})
    if TCGA_COX_ADJ.exists():
        adj = pd.read_csv(TCGA_COX_ADJ, sep="\t")
        for r in adj.itertuples(index=False):
            rows.append({"evidence_layer": f"TCGA_DSS_Cox_{r.adjust}", "metric": r.module, "effect": f"HR={r.HR:.3g}", "p": r.p})
    pd.DataFrame(rows).to_csv(OUT / "external_evidence_summary.tsv", sep="\t", index=False)


def grade_topics() -> pd.DataFrame:
    condition = pd.read_csv(OUT / "spatial_cnv_condition_tests.tsv", sep="\t")
    corr = pd.read_csv(OUT / "spatial_cnv_axis_correlations.tsv", sep="\t")
    moran = pd.read_csv(OUT / "spatial_cnv_morans.tsv", sep="\t")
    ec = pd.read_csv(OUT / "spatial_ecotype_summary.tsv", sep="\t")
    retrieval = pd.read_csv(OUT / "uni_retrieval_metrics.tsv", sep="\t") if (OUT / "uni_retrieval_metrics.tsv").exists() else pd.DataFrame()

    ptc_test = condition[condition["comparison"].isin(["PTC_vs_N", "LPTC_vs_N", "ATC_vs_N"])]
    max_sig_d = float(ptc_test["cohen_d"].abs().max()) if not ptc_test.empty else np.nan
    max_high = float(ptc_test["q95_high_rate_a"].max()) if not ptc_test.empty else np.nan
    sig_moran = moran[(moran["axis"] == "t00_spatial_cnv_signature") & (moran["condition"].isin(["PTC", "LPTC", "ATC"]))]
    median_moran = float(sig_moran["moran_i_k6"].median()) if not sig_moran.empty else np.nan

    def corr_val(target: str, scope: str = "pooled_cancer_spots") -> float:
        hit = corr[(corr["target"] == target) & (corr["scope"] == scope)]
        return float(hit.iloc[0]["spearman_r"]) if not hit.empty else np.nan

    def ret_val(target: str) -> tuple[str, float, float]:
        if retrieval.empty:
            return "NA", np.nan, np.nan
        sub = retrieval[retrieval["target"] == target].sort_values("spearman_r", ascending=False)
        if sub.empty:
            return "NA", np.nan, np.nan
        r = sub.iloc[0]
        return str(r["mode"]), float(r["spearman_r"]), float(r["auroc_q75"])

    sig_mode, sig_ret, sig_auc = ret_val("t00_spatial_cnv_signature")
    trop2_mode, trop2_ret, trop2_auc = ret_val("TACSTD2_z")
    caf_mode, caf_ret, caf_auc = ret_val("CAF_ECM_z")
    tls_mode, tls_ret, tls_auc = ret_val("TLS_B_z")

    rows = [
        {
            "rank": 1,
            "topic": "H1 CNV-residual spatial ecosystem",
            "grade_after": "A-",
            "key_new_validation": f"spatial expression-CNV stage effect max |d|={max_sig_d:.2f}; normal-q95 high rate max={max_high:.2%}; median Moran I={median_moran:.3f}",
            "decision": "Flagship; now has a spatial-CNV bridge, still needs orthogonal wet/IF validation for A/A+.",
        },
        {
            "rank": 2,
            "topic": "H2 spatial phylogeography of driver-negative THCA",
            "grade_after": "B+",
            "key_new_validation": f"spatial-CNV territories are measurable; signature Moran I median={median_moran:.3f}",
            "decision": "Promising but still expression-CNV proxy, not allele-specific CalicoST-grade clone phylogeny.",
        },
        {
            "rank": 3,
            "topic": "H3 multimodal visibility boundary",
            "grade_after": "B+",
            "key_new_validation": f"UNI retrieval best signature {sig_mode} rho={sig_ret:.3f}/AUC={sig_auc:.3f}; TROP2 {trop2_mode} rho={trop2_ret:.3f}/AUC={trop2_auc:.3f}",
            "decision": "Companion method boundary; use as evidence-grounded retrieval, not as main claim.",
        },
        {
            "rank": 4,
            "topic": "H4 thyroid spatial ecotypes",
            "grade_after": "B+",
            "key_new_validation": f"k-means ecotypes recovered n={ec.shape[0]} ecosystem states with labels: {', '.join(ec['ecotype_label'].astype(str).head(6))}",
            "decision": "Good integrated figure, but not enough as standalone without external ecotype cohort.",
        },
        {
            "rank": 5,
            "topic": "H5 CAF/ECM lineage suppression",
            "grade_after": "B+",
            "key_new_validation": f"spatial-CNV vs CAF_ECM pooled rho={corr_val('CAF_ECM_z'):.3f}; UNI retrieval {caf_mode} rho={caf_ret:.3f}/AUC={caf_auc:.3f}",
            "decision": "Mechanism layer inside H1; keep stromal abundance caveat.",
        },
        {
            "rank": 6,
            "topic": "H6 TROP2 spatial niche",
            "grade_after": "A-",
            "key_new_validation": f"spatial-CNV vs TACSTD2 pooled rho={corr_val('TACSTD2_z'):.3f}; UNI retrieval {trop2_mode} rho={trop2_ret:.3f}/AUC={trop2_auc:.3f}",
            "decision": "Strong sub-aim; do not lead as ADC story.",
        },
        {
            "rank": 7,
            "topic": "H7 immune/TLS protective thyroid ecotype",
            "grade_after": "B+",
            "key_new_validation": f"spatial-CNV vs TLS_B rho={corr_val('TLS_B_z'):.3f}; retrieval {tls_mode} rho={tls_ret:.3f}/AUC={tls_auc:.3f}; TCGA Cox retained separately.",
            "decision": "Useful ecosystem axis; still needs Hashimoto-vs-antitumor separation.",
        },
    ]
    score = pd.DataFrame(rows)
    score.to_csv(OUT / "top_topic_validation_scorecard.tsv", sep="\t", index=False)
    return score


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) > 0 and abs(x) < 1e-4:
            return f"{x:.2e}"
        return f"{x:.3g}"
    return str(x)


def write_report(score: pd.DataFrame) -> None:
    cond = pd.read_csv(OUT / "spatial_cnv_condition_tests.tsv", sep="\t")
    arms = pd.read_csv(OUT / "spatial_cnv_arm_condition_tests.tsv", sep="\t")
    corr = pd.read_csv(OUT / "spatial_cnv_axis_correlations.tsv", sep="\t")
    moran = pd.read_csv(OUT / "spatial_cnv_morans.tsv", sep="\t")
    ec = pd.read_csv(OUT / "spatial_ecotype_summary.tsv", sep="\t")
    ext = pd.read_csv(OUT / "external_evidence_summary.tsv", sep="\t")
    retrieval = pd.read_csv(OUT / "uni_retrieval_metrics.tsv", sep="\t") if (OUT / "uni_retrieval_metrics.tsv").exists() else pd.DataFrame()

    lines: list[str] = []
    lines.append("# High-IF top-topic validation - 2026-05-08\n")
    lines.append("Scope: H1-H7 high-probability directions from the high-IF strategy memo. This run added a spatial expression-CNV bridge using downloaded GENCODE v48 and UCSC hg38 cytoband annotations.\n")

    lines.append("## New data added\n")
    lines.append("| Data | Local file | Use |")
    lines.append("|---|---|---|")
    lines.append(f"| GENCODE v48 gene annotation | `{GTF}` | map Visium genes to chromosomes and arms |")
    lines.append(f"| UCSC hg38 cytoband | `{CYTOBAND}` | assign p/q arm labels |")
    lines.append("| Existing public thyroid spatial corpus | `project/data/processed/GSE250521/*.h5ad` | infer spot-level expression-CNV and ecosystems |")
    lines.append("| Existing UNI morphology embeddings | `phase1_gse250521/uni_embeddings_size224.npz` | evidence-grounded multimodal retrieval pilot |")
    lines.append("| Existing TCGA/cBio/GATCI/GeoMx evidence | `high_impact_topic_pilots_2026_05_08/*` | orthogonal validation layers |\n")

    lines.append("## Topic verdict\n")
    lines.append("| Rank | Topic | Grade | Key validation | Decision |")
    lines.append("|---:|---|---:|---|---|")
    for r in score.itertuples(index=False):
        lines.append(f"| {r.rank} | {r.topic} | {r.grade_after} | {r.key_new_validation} | {r.decision} |")

    lines.append("\n## H1/H2 spatial expression-CNV bridge\n")
    lines.append("| Comparison | mean condition | mean normal | Cohen d | q95 high rate | OR | p |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for r in cond.itertuples(index=False):
        lines.append(f"| {r.comparison} | {fmt(r.mean_a)} | {fmt(r.mean_normal)} | {fmt(r.cohen_d)} | {fmt(r.q95_high_rate_a)} | {fmt(r.fisher_or)} | {fmt(r.fisher_p)} |")

    lines.append("\nTop arm-level spatial-expression shifts:")
    lines.append("\n| Arm | Condition | Direction | Delta vs normal | d | p |")
    lines.append("|---|---|---|---:|---:|---:|")
    arm_show = arms.reindex(arms["cohen_d"].abs().sort_values(ascending=False).index).head(12)
    for r in arm_show.itertuples(index=False):
        lines.append(f"| {r.arm} | {r.condition} | {r.direction} | {fmt(r.delta_vs_normal)} | {fmt(r.cohen_d)} | {fmt(r.mw_p)} |")

    lines.append("\nSpatial organization:")
    mshow = moran[(moran["axis"] == "t00_spatial_cnv_signature")].groupby("condition")["moran_i_k6"].agg(["count", "median", "mean"]).reset_index()
    lines.append("\n| Condition | n slides | median Moran I | mean Moran I |")
    lines.append("|---|---:|---:|---:|")
    for r in mshow.itertuples(index=False):
        lines.append(f"| {r.condition} | {r.count} | {fmt(r.median)} | {fmt(r.mean)} |")

    lines.append("\n## Ecosystem coupling\n")
    cshow = corr[corr["scope"].isin(["pooled_cancer_spots", "slide_median_r"]) & corr["target"].isin(["TACSTD2_z", "CAF_ECM_z", "TLS_B_z", "Tcell_cytotoxic_z", "Macrophage_TAM_z", "RAI_thyroid_z", "DM1_like_score", "RAI_8_score"])].copy()
    lines.append("| Scope | Target | Spearman r | p/n |")
    lines.append("|---|---|---:|---:|")
    for r in cshow.itertuples(index=False):
        pn = f"p={fmt(r.p)}" if pd.notna(r.p) else f"n_slides={r.n}"
        lines.append(f"| {r.scope} | {r.target} | {fmt(r.spearman_r)} | {pn} |")

    lines.append("\n## Spatial ecotypes\n")
    lines.append("| Ecotype | Label | n spots | n slides | CNV | TROP2 | CAF | TLS | Tcell | TAM | RAI | high-q95 rate |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in ec.itertuples(index=False):
        lines.append(f"| {r.spatial_ecotype} | {r.ecotype_label} | {r.n_spots} | {r.n_slides} | {fmt(r.mean_signature)} | {fmt(r.mean_TROP2)} | {fmt(r.mean_CAF)} | {fmt(r.mean_TLS)} | {fmt(r.mean_Tcell)} | {fmt(r.mean_TAM)} | {fmt(r.mean_RAI)} | {fmt(r.high_q95_rate)} |")

    if not retrieval.empty:
        lines.append("\n## H3 UNI morphology + metadata retrieval\n")
        lines.append("| Target | Best mode | Spearman r | AUROC q75 | n |")
        lines.append("|---|---|---:|---:|---:|")
        for target, sub in retrieval.groupby("target"):
            r = sub.sort_values("spearman_r", ascending=False).iloc[0]
            lines.append(f"| {target} | {r['mode']} | {fmt(r['spearman_r'])} | {fmt(r['auroc_q75'])} | {int(r['n'])} |")

    lines.append("\n## External evidence snapshot\n")
    lines.append("| Layer | Metric | Effect | p |")
    lines.append("|---|---|---|---:|")
    for r in ext.head(28).itertuples(index=False):
        lines.append(f"| {r.evidence_layer} | {r.metric} | {r.effect} | {fmt(r.p)} |")

    lines.append("\n## Bottom line\n")
    lines.append("H1 remains the main high-IF bet. The new run adds the missing bridge: a chromosome-arm spatial expression-CNV layer that can be overlaid with TROP2, CAF/ECM, TLS/immune, RAI, and morphology retrieval. H2 is now a real next experiment, but it remains a proxy until allele-specific CalicoST-grade clone inference is run. H3-H7 should be kept as integrated layers, not separate headline papers.\n")

    lines.append("## Output files\n")
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for p in [GTF, CYTOBAND, SAMPLE_META, JOINT, TROP2]:
        require(p)
    per_spot = score_all_spots()
    run_spatial_validations(per_spot)
    run_ecotype_validation(per_spot)
    run_uni_retrieval(per_spot)
    run_external_summary()
    score = grade_topics()
    write_report(score)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
