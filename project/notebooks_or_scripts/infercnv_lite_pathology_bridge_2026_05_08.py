#!/usr/bin/env python3
"""Spatial transcriptome inferCNV-lite plus pathology bridge.

This is a practical next gate after the high-IF topic screen:

1. Use GENCODE gene order + normal thyroid Visium spots to create a
   gene-order/bin-level inferCNV-like profile for each GSE250521 spot.
2. Cluster spots within each slide into CNV-like territories.
3. Rank territories by the T00 arm signature and by overall aneuploidy.
4. Test whether pathology image embeddings can recover the ST-derived
   CNV territories in leave-one-slide-out evaluation.
"""

from __future__ import annotations

import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_pathology_bridge_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_infercnv_lite_pathology_bridge.md"

SAMPLE_META = ROOT / "project/data/processed/GSE250521/sample_metadata.tsv"
GENE_ARM = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_gene_arm_map_used.tsv"
ST_PER_SPOT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_expression_cnv_per_spot.tsv.gz"
UNI_META = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embed_metadata_size224.tsv"
UNI_NPZ = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embeddings_size224.npz"

GAIN_ARMS = ["7p", "7q", "12q", "16p", "16q"]
LOSS_ARMS = ["2p", "2q"]
ARM_COLS = [f"arm_{a}" for a in GAIN_ARMS + LOSS_ARMS]
ECOSYSTEM_AXES = [
    "TACSTD2_z",
    "TROP2_raw",
    "CAF_ECM_z",
    "TLS_B_z",
    "Tcell_cytotoxic_z",
    "Macrophage_TAM_z",
    "RAI_thyroid_z",
    "EMT_stress_z",
    "DM1_like_score",
    "RAI_8_score",
    "n_fib",
    "n_lym",
]
MORPH_AXES = [
    "B1_stromal_encased_thyrocyte_idx",
    "B3_nuclear_eccentricity_var",
    "B4_tumor_stroma_ratio",
    "C2_spatial_celltype_entropy",
    "B5_thyrocyte_cluster_med",
    "n_thy",
    "n_lym",
    "n_fib",
]


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


def cohen_d(a, b) -> float:
    a = pd.Series(a).dropna().astype(float)
    b = pd.Series(b).dropna().astype(float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)
    if pooled <= 0 or np.isnan(pooled):
        return np.nan
    return float((a.mean() - b.mean()) / math.sqrt(pooled))


def spearman(x, y) -> tuple[float, float, int]:
    tmp = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if tmp.shape[0] < 5 or tmp["x"].nunique() < 2 or tmp["y"].nunique() < 2:
        return np.nan, np.nan, int(tmp.shape[0])
    r, p = stats.spearmanr(tmp["x"], tmp["y"])
    return float(r), float(p), int(tmp.shape[0])


def build_gene_reference(samples: pd.DataFrame, gene_names: list[str], gene_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    cache = OUT / "normal_reference_gene_mean_sd.tsv"
    if cache.exists():
        ref = pd.read_csv(cache, sep="\t")
        ref = ref.set_index("gene").reindex(gene_names)
        return ref["normal_mean"].to_numpy(float), ref["normal_sd"].to_numpy(float)

    sums = np.zeros(len(gene_names), dtype=np.float64)
    sums2 = np.zeros(len(gene_names), dtype=np.float64)
    n_obs = 0
    for row in samples[samples["condition"] == "N"].itertuples(index=False):
        a = ad.read_h5ad(ROOT / row.h5ad)
        x = normalize_log1p_sparse(a.X)[:, gene_idx].tocsr()
        sums += np.asarray(x.sum(axis=0)).ravel()
        sums2 += np.asarray(x.power(2).sum(axis=0)).ravel()
        n_obs += x.shape[0]
    mean = sums / max(n_obs, 1)
    var = np.maximum(sums2 / max(n_obs, 1) - mean**2, 1e-6)
    sd = np.sqrt(var)
    sd[sd < 0.05] = 0.05
    pd.DataFrame({"gene": gene_names, "normal_mean": mean, "normal_sd": sd}).to_csv(cache, sep="\t", index=False)
    return mean, sd


def make_bins(gene_map: pd.DataFrame, gene_names: list[str], genes_per_bin: int = 100) -> pd.DataFrame:
    keep = gene_map[gene_map["gene_name"].isin(gene_names)].copy()
    keep = keep[keep["chrom"].astype(str).isin([str(i) for i in range(1, 23)])].copy()
    keep["chrom_i"] = keep["chrom"].astype(int)
    keep = keep.sort_values(["chrom_i", "start", "gene_name"]).reset_index(drop=True)
    rows = []
    for chrom, sub in keep.groupby("chrom_i", sort=True):
        sub = sub.reset_index(drop=True)
        for i, start in enumerate(range(0, sub.shape[0], genes_per_bin)):
            part = sub.iloc[start : start + genes_per_bin]
            if part.shape[0] < 30:
                continue
            rows.append(
                pd.DataFrame(
                    {
                        "gene_name": part["gene_name"].to_numpy(),
                        "chrom": chrom,
                        "bin": f"chr{chrom:02d}_bin{i:03d}",
                        "bin_index": len(rows),
                        "bin_start": int(part["start"].min()),
                        "bin_end": int(part["end"].max()),
                    }
                )
            )
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(OUT / "infercnv_lite_gene_bins.tsv", sep="\t", index=False)
    return out


def build_weight_matrix(gene_names: list[str], groups: dict[str, list[str]], mean: np.ndarray, sd: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray, list[str]]:
    gene_to_i = {g: i for i, g in enumerate(gene_names)}
    row_idx = []
    col_idx = []
    data = []
    offsets = []
    names = []
    for name, genes in groups.items():
        idx = [gene_to_i[g] for g in genes if g in gene_to_i]
        if not idx:
            continue
        idx_arr = np.asarray(idx, dtype=int)
        weight = 1.0 / (sd[idx_arr] * len(idx_arr))
        row_idx.extend(idx)
        col_idx.extend([len(names)] * len(idx))
        data.extend(weight.tolist())
        offsets.append(float(np.sum(mean[idx_arr] / sd[idx_arr]) / len(idx_arr)))
        names.append(name)
    w = sparse.csr_matrix((data, (row_idx, col_idx)), shape=(len(gene_names), len(names)))
    return w, np.asarray(offsets), names


def smooth_bins(mat: np.ndarray, bin_meta: pd.DataFrame, bin_names: list[str]) -> np.ndarray:
    out = mat.copy()
    meta = bin_meta.drop_duplicates("bin").set_index("bin").loc[bin_names].reset_index()
    for chrom, idx in meta.groupby("chrom").groups.items():
        cols = np.asarray(list(idx), dtype=int)
        if len(cols) < 3:
            continue
        block = mat[:, cols]
        sm = block.copy()
        sm[:, 1:-1] = 0.25 * block[:, :-2] + 0.5 * block[:, 1:-1] + 0.25 * block[:, 2:]
        sm[:, 0] = 0.67 * block[:, 0] + 0.33 * block[:, 1]
        sm[:, -1] = 0.67 * block[:, -1] + 0.33 * block[:, -2]
        out[:, cols] = sm
    return out


def knn_edges(coords: np.ndarray, k: int = 6) -> np.ndarray:
    n = coords.shape[0]
    nn = NearestNeighbors(n_neighbors=min(k + 1, n)).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)[:, 1:]
    src = np.repeat(np.arange(n), idx.shape[1])
    dst = idx.reshape(-1)
    keep = src < dst
    return np.column_stack([src[keep], dst[keep]])


def same_label_fraction(edges: np.ndarray, labels: np.ndarray) -> float:
    if edges.size == 0:
        return np.nan
    return float(np.mean(labels[edges[:, 0]] == labels[edges[:, 1]]))


def assign_slide_clones(frame: pd.DataFrame, bin_scores: np.ndarray, bin_names: list[str]) -> pd.DataFrame:
    sub = frame.copy().reset_index(drop=True)
    x = StandardScaler().fit_transform(np.nan_to_num(bin_scores, nan=0.0))
    labels = MiniBatchKMeans(n_clusters=3, random_state=8, n_init=12, batch_size=2048, max_iter=200).fit_predict(x)
    sub["infercnv_clone_raw"] = labels
    sub["infercnv_aneuploidy_score"] = np.mean(np.abs(bin_scores), axis=1)

    gain = [f"arm_{a}" for a in GAIN_ARMS if f"arm_{a}" in sub.columns]
    loss = [f"arm_{a}" for a in LOSS_ARMS if f"arm_{a}" in sub.columns]
    sub["infercnv_t00_signature"] = sub[gain].mean(axis=1) - sub[loss].mean(axis=1)
    sub["infercnv_bin_sd"] = np.std(bin_scores, axis=1)

    clone = (
        sub.groupby("infercnv_clone_raw")
        .agg(mean_aneuploidy=("infercnv_aneuploidy_score", "mean"), mean_t00=("infercnv_t00_signature", "mean"))
        .reset_index()
    )
    aneu_rank = {raw: f"A{rank}_{name}" for rank, (raw, name) in enumerate(zip(clone.sort_values("mean_aneuploidy")["infercnv_clone_raw"], ["low", "mid", "high"]))}
    t00_rank = {raw: f"T{rank}_{name}" for rank, (raw, name) in enumerate(zip(clone.sort_values("mean_t00")["infercnv_clone_raw"], ["low", "mid", "high"]))}
    sub["aneuploidy_territory"] = sub["infercnv_clone_raw"].map(aneu_rank)
    sub["t00_territory"] = sub["infercnv_clone_raw"].map(t00_rank)
    sub["aneuploidy_clone_high"] = sub["aneuploidy_territory"].str.contains("high")
    sub["t00_clone_high"] = sub["t00_territory"].str.contains("high")
    return sub


def compute_infercnv_lite() -> pd.DataFrame:
    per_spot_path = OUT / "infercnv_lite_per_spot.tsv.gz"
    if per_spot_path.exists():
        return pd.read_csv(per_spot_path, sep="\t")

    samples = pd.read_csv(SAMPLE_META, sep="\t")
    samples["condition"] = samples["sample_id"].map(condition_from_sample)
    gene_map = pd.read_csv(GENE_ARM, sep="\t")
    gene_map = gene_map[gene_map["chrom"].astype(str).isin([str(i) for i in range(1, 23)])].copy()

    first = ad.read_h5ad(ROOT / samples.iloc[0]["h5ad"], backed="r")
    var_names = pd.Index(first.var_names.astype(str))
    use_genes = gene_map[gene_map["gene_name"].isin(var_names)].copy()
    use_genes = use_genes.sort_values(["chrom", "start", "gene_name"]).drop_duplicates("gene_name")
    gene_names = use_genes["gene_name"].tolist()
    gene_idx = var_names.get_indexer(gene_names)
    ok = gene_idx >= 0
    gene_names = [g for g, keep in zip(gene_names, ok) if keep]
    gene_idx = gene_idx[ok]
    use_genes = use_genes[use_genes["gene_name"].isin(gene_names)].copy()

    mean, sd = build_gene_reference(samples, gene_names, gene_idx)
    bin_meta = make_bins(use_genes, gene_names, genes_per_bin=100)
    bin_groups = bin_meta.groupby("bin")["gene_name"].apply(list).to_dict()
    arm_groups = use_genes.groupby("arm")["gene_name"].apply(list).to_dict()
    arm_groups = {k: v for k, v in arm_groups.items() if k in GAIN_ARMS + LOSS_ARMS}
    bin_w, bin_offsets, bin_names = build_weight_matrix(gene_names, bin_groups, mean, sd)
    arm_w, arm_offsets, arm_names = build_weight_matrix(gene_names, arm_groups, mean, sd)

    chunks = []
    clone_rows = []
    coherence_rows = []
    stability_rows = []
    rng = np.random.default_rng(8)

    for row in samples.itertuples(index=False):
        a = ad.read_h5ad(ROOT / row.h5ad)
        x = normalize_log1p_sparse(a.X)[:, gene_idx].tocsr()
        bin_scores = x @ bin_w
        bin_scores = bin_scores.toarray() if sparse.issparse(bin_scores) else np.asarray(bin_scores)
        bin_scores = smooth_bins(bin_scores - bin_offsets[None, :], bin_meta, bin_names)

        arm_scores = x @ arm_w
        arm_scores = arm_scores.toarray() if sparse.issparse(arm_scores) else np.asarray(arm_scores)
        arm_scores = arm_scores - arm_offsets[None, :]

        obs = a.obs.reset_index().rename(columns={"index": "spot_id"})
        frame = obs[["spot_id", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]].copy()
        frame["sample_id"] = row.sample_id
        frame["condition"] = row.condition
        for i, arm in enumerate(arm_names):
            frame[f"arm_{arm}"] = arm_scores[:, i]

        frame = assign_slide_clones(frame, bin_scores, bin_names)

        clone_summary = (
            frame.groupby(["sample_id", "condition", "infercnv_clone_raw", "aneuploidy_territory", "t00_territory"])
            .agg(
                n_spots=("spot_id", "size"),
                mean_aneuploidy=("infercnv_aneuploidy_score", "mean"),
                mean_t00_signature=("infercnv_t00_signature", "mean"),
                mean_bin_sd=("infercnv_bin_sd", "mean"),
            )
            .reset_index()
        )
        for arm in ARM_COLS:
            if arm in frame.columns:
                clone_summary[arm] = frame.groupby("infercnv_clone_raw")[arm].mean().reindex(clone_summary["infercnv_clone_raw"]).to_numpy()
        clone_rows.append(clone_summary)

        coords = frame[["array_row", "array_col"]].to_numpy(float)
        edges = knn_edges(coords, k=6)
        for label_col in ["t00_territory", "aneuploidy_territory"]:
            labels = frame[label_col].to_numpy()
            obs_same = same_label_fraction(edges, labels)
            null = np.array([same_label_fraction(edges, rng.permutation(labels)) for _ in range(300)])
            coherence_rows.append(
                {
                    "sample_id": row.sample_id,
                    "condition": row.condition,
                    "label": label_col,
                    "same_knn_obs": obs_same,
                    "same_knn_null_mean": float(null.mean()),
                    "same_knn_z": float((obs_same - null.mean()) / null.std(ddof=1)) if null.std(ddof=1) > 0 else np.nan,
                    "perm_p": float((1 + np.sum(null >= obs_same)) / (len(null) + 1)),
                }
            )

        base = frame["t00_territory"].str.extract(r"T([0-9])_")[0].astype(int).to_numpy()
        aris = []
        xz = StandardScaler().fit_transform(np.nan_to_num(bin_scores, nan=0.0))
        for seed in range(8):
            cols = rng.choice(np.arange(xz.shape[1]), size=max(50, int(xz.shape[1] * 0.75)), replace=False)
            lab = MiniBatchKMeans(n_clusters=3, random_state=200 + seed, n_init=5, batch_size=2048, max_iter=120).fit_predict(xz[:, cols])
            tmp = pd.DataFrame({"lab": lab, "sig": frame["infercnv_t00_signature"]})
            means = tmp.groupby("lab")["sig"].mean().sort_values()
            rank = {raw: i for i, raw in enumerate(means.index)}
            ordered = np.array([rank[x] for x in lab])
            aris.append(adjusted_rand_score(base, ordered))
        stability_rows.append(
            {
                "sample_id": row.sample_id,
                "condition": row.condition,
                "median_ARI_bin_bootstrap": float(np.median(aris)),
                "min_ARI_bin_bootstrap": float(np.min(aris)),
                "max_ARI_bin_bootstrap": float(np.max(aris)),
            }
        )
        chunks.append(frame)

    per_spot = pd.concat(chunks, ignore_index=True)
    st = pd.read_csv(ST_PER_SPOT, sep="\t")
    keep = ["sample_id", "spot_id"] + [c for c in ECOSYSTEM_AXES + MORPH_AXES if c in st.columns]
    per_spot = per_spot.merge(st[keep].drop_duplicates(["sample_id", "spot_id"]), on=["sample_id", "spot_id"], how="left")

    per_spot.to_csv(per_spot_path, sep="\t", index=False, compression="gzip")
    pd.concat(clone_rows, ignore_index=True).to_csv(OUT / "infercnv_lite_clone_summary.tsv", sep="\t", index=False)
    pd.DataFrame(coherence_rows).to_csv(OUT / "infercnv_lite_clone_coherence.tsv", sep="\t", index=False)
    pd.DataFrame(stability_rows).to_csv(OUT / "infercnv_lite_clone_stability.tsv", sep="\t", index=False)
    pd.DataFrame({"bin": bin_names}).to_csv(OUT / "infercnv_lite_bin_names.tsv", sep="\t", index=False)
    return per_spot


def ecosystem_enrichment(per_spot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    for scope, sub in [("pooled_cancer_spots", cancer)] + list(cancer.groupby("sample_id")):
        high = sub[sub["t00_clone_high"]]
        low = sub[sub["t00_territory"].str.contains("low")]
        for axis in [a for a in ECOSYSTEM_AXES if a in sub.columns]:
            if high[axis].dropna().shape[0] < 2 or low[axis].dropna().shape[0] < 2:
                continue
            rows.append(
                {
                    "scope": scope,
                    "condition": sub["condition"].iloc[0] if scope != "pooled_cancer_spots" else "PTC_LPTC_ATC",
                    "axis": axis,
                    "n_high": high[axis].dropna().shape[0],
                    "n_low": low[axis].dropna().shape[0],
                    "mean_high": high[axis].mean(),
                    "mean_low": low[axis].mean(),
                    "delta_high_low": high[axis].mean() - low[axis].mean(),
                    "cohen_d": cohen_d(high[axis], low[axis]),
                    "mw_p": float(stats.mannwhitneyu(high[axis].dropna(), low[axis].dropna()).pvalue),
                }
            )

    df = pd.DataFrame(rows)
    med = []
    for axis, sub in df[~df["scope"].eq("pooled_cancer_spots")].groupby("axis"):
        med.append(
            {
                "scope": "slide_median_effect",
                "condition": "PTC_LPTC_ATC",
                "axis": axis,
                "n_high": sub.shape[0],
                "n_low": sub.shape[0],
                "mean_high": np.nan,
                "mean_low": np.nan,
                "delta_high_low": float(sub["delta_high_low"].median()),
                "cohen_d": float(sub["cohen_d"].median()),
                "mw_p": np.nan,
            }
        )
    df = pd.concat([df, pd.DataFrame(med)], ignore_index=True)
    df.to_csv(OUT / "infercnv_lite_ecosystem_enrichment.tsv", sep="\t", index=False)
    return df


def pathology_bridge(per_spot: pd.DataFrame) -> pd.DataFrame:
    if not UNI_META.exists() or not UNI_NPZ.exists():
        return pd.DataFrame()

    meta = pd.read_csv(UNI_META, sep="\t").rename(columns={"slide": "sample_id"})
    emb = np.load(UNI_NPZ)["embeddings"]
    if emb.shape[0] != meta.shape[0]:
        raise RuntimeError("UNI embedding row count does not match metadata")

    label_cols = [
        "condition",
        "infercnv_t00_signature",
        "infercnv_aneuploidy_score",
        "infercnv_bin_sd",
        "t00_clone_high",
        "aneuploidy_clone_high",
    ] + [c for c in MORPH_AXES if c in per_spot.columns]
    dat = meta.merge(per_spot[["sample_id", "spot_id"] + label_cols], on=["sample_id", "spot_id"], how="inner")
    emb = emb[dat.index.to_numpy()]

    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    text_all = enc.fit_transform(dat[["condition"]])
    morph_cols = [c for c in MORPH_AXES if c in dat.columns]
    morph_all = dat[morph_cols].replace([np.inf, -np.inf], np.nan).fillna(dat[morph_cols].median(numeric_only=True)).to_numpy(float)

    pred_rows = []
    metric_rows = []
    continuous_targets = ["infercnv_t00_signature", "infercnv_aneuploidy_score", "infercnv_bin_sd"]
    binary_targets = ["t00_clone_high", "aneuploidy_clone_high"]

    for held in sorted(dat["sample_id"].unique()):
        train = dat["sample_id"] != held
        test = ~train
        if train.sum() < 100 or test.sum() < 20:
            continue

        uni_scaler = StandardScaler().fit(emb[train])
        uni_train = uni_scaler.transform(emb[train])
        uni_test = uni_scaler.transform(emb[test])
        n_comp = min(64, uni_train.shape[1], train.sum() - 1)
        pca = PCA(n_components=n_comp, random_state=8).fit(uni_train)
        uni_train = pca.transform(uni_train)
        uni_test = pca.transform(uni_test)

        morph_scaler = StandardScaler().fit(morph_all[train])
        morph_train = morph_scaler.transform(morph_all[train])
        morph_test = morph_scaler.transform(morph_all[test])
        text_train = text_all[train]
        text_test = text_all[test]

        modes = {
            "UNI_image_only": (uni_train, uni_test),
            "cell_morph_only": (morph_train, morph_test),
            "stage_text_only": (text_train, text_test),
            "UNI_plus_cell_morph": (np.hstack([uni_train, morph_train]), np.hstack([uni_test, morph_test])),
            "UNI_plus_stage_text": (np.hstack([uni_train, text_train]), np.hstack([uni_test, text_test])),
            "UNI_cell_stage": (np.hstack([uni_train, morph_train, text_train]), np.hstack([uni_test, morph_test, text_test])),
        }

        for mode, (xtr, xte) in modes.items():
            for target in continuous_targets:
                ytr = dat.loc[train, target].to_numpy(float)
                yte = dat.loc[test, target].to_numpy(float)
                model = Ridge(alpha=10.0)
                model.fit(xtr, ytr)
                pred = model.predict(xte)
                for i, p in enumerate(pred):
                    pred_rows.append(
                        {
                            "sample_id": dat.loc[test].iloc[i]["sample_id"],
                            "spot_id": dat.loc[test].iloc[i]["spot_id"],
                            "condition": dat.loc[test].iloc[i]["condition"],
                            "mode": mode,
                            "target": target,
                            "observed": yte[i],
                            "predicted": p,
                        }
                    )
            for target in binary_targets:
                ytr = dat.loc[train, target].astype(int).to_numpy()
                yte = dat.loc[test, target].astype(int).to_numpy()
                if len(np.unique(ytr)) < 2:
                    continue
                model = LogisticRegression(max_iter=2000, class_weight="balanced", C=0.3)
                model.fit(xtr, ytr)
                prob = model.predict_proba(xte)[:, 1]
                for i, p in enumerate(prob):
                    pred_rows.append(
                        {
                            "sample_id": dat.loc[test].iloc[i]["sample_id"],
                            "spot_id": dat.loc[test].iloc[i]["spot_id"],
                            "condition": dat.loc[test].iloc[i]["condition"],
                            "mode": mode,
                            "target": target,
                            "observed": yte[i],
                            "predicted": p,
                        }
                    )

    pred = pd.DataFrame(pred_rows)
    pred.to_csv(OUT / "pathology_to_infercnv_predictions.tsv.gz", sep="\t", index=False, compression="gzip")

    for (mode, target), sub in pred.groupby(["mode", "target"]):
        if target in continuous_targets:
            r, p, n = spearman(sub["observed"], sub["predicted"])
            metric_rows.append({"mode": mode, "target": target, "metric": "spearman", "value": r, "p": p, "n": n})
        else:
            labels = sub["observed"].astype(int)
            auc = roc_auc_score(labels, sub["predicted"]) if labels.nunique() == 2 else np.nan
            metric_rows.append({"mode": mode, "target": target, "metric": "AUROC", "value": auc, "p": np.nan, "n": sub.shape[0]})

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(OUT / "pathology_to_infercnv_metrics.tsv", sep="\t", index=False)
    return metrics


def arm_direction_audit(per_spot: pd.DataFrame) -> pd.DataFrame:
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    high = cancer[cancer["t00_clone_high"]]
    low = cancer[cancer["t00_territory"].str.contains("low")]
    rows = []
    for arm in GAIN_ARMS + LOSS_ARMS:
        col = f"arm_{arm}"
        if col not in cancer.columns:
            continue
        expected = "high>low" if arm in GAIN_ARMS else "high<low"
        delta = high[col].mean() - low[col].mean()
        rows.append(
            {
                "arm": arm,
                "expected": expected,
                "mean_high": high[col].mean(),
                "mean_low": low[col].mean(),
                "delta_high_low": delta,
                "direction_matches": (delta > 0 and expected == "high>low") or (delta < 0 and expected == "high<low"),
                "cohen_d": cohen_d(high[col], low[col]),
                "mw_p": float(stats.mannwhitneyu(high[col].dropna(), low[col].dropna()).pvalue),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "infercnv_lite_t00_arm_direction_audit.tsv", sep="\t", index=False)
    return out


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) > 0 and abs(x) < 1e-4:
            return f"{x:.2e}"
        return f"{x:.3g}"
    return str(x)


def write_report(per_spot: pd.DataFrame, eco: pd.DataFrame, path_metrics: pd.DataFrame, arms: pd.DataFrame) -> None:
    clone = pd.read_csv(OUT / "infercnv_lite_clone_summary.tsv", sep="\t")
    coh = pd.read_csv(OUT / "infercnv_lite_clone_coherence.tsv", sep="\t")
    stab = pd.read_csv(OUT / "infercnv_lite_clone_stability.tsv", sep="\t")
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    cancer_coh = coh[(coh["condition"].isin(["PTC", "LPTC", "ATC"])) & (coh["label"] == "t00_territory")]

    lines = []
    lines.append("# inferCNV-lite and pathology bridge - 2026-05-08\n")
    lines.append("Goal: create a spatial transcriptome-based CNV-like clone/territory layer, then test whether pathology image embeddings can predict it.\n")

    lines.append("## Bottom line\n")
    lines.append(f"- ST inferCNV-lite territories are spatially coherent: cancer-slide median same-label KNN z **{fmt(cancer_coh['same_knn_z'].median())}**.")
    lines.append(f"- ST territory stability is usable: cancer-slide median bin-bootstrap ARI **{fmt(stab[stab['condition'].isin(['PTC','LPTC','ATC'])]['median_ARI_bin_bootstrap'].median())}**.")
    t00_auc = path_metrics[(path_metrics["mode"] == "UNI_image_only") & (path_metrics["target"] == "t00_clone_high") & (path_metrics["metric"] == "AUROC")]
    t00_r = path_metrics[(path_metrics["mode"] == "UNI_image_only") & (path_metrics["target"] == "infercnv_t00_signature") & (path_metrics["metric"] == "spearman")]
    lines.append(f"- Pathology-only bridge to T00 territory: UNI AUROC **{fmt(t00_auc.iloc[0]['value'] if not t00_auc.empty else np.nan)}**, signature rho **{fmt(t00_r.iloc[0]['value'] if not t00_r.empty else np.nan)}**.")
    lines.append("- Interpretation: ST-CNV territories are worth a figure. Pathology image prediction is only useful if UNI_image_only is clearly above stage/text baselines; otherwise it belongs as a negative boundary.\n")

    lines.append("## ST inferCNV-lite clone/territory evidence\n")
    lines.append("| Metric | Value |")
    lines.append("|---|---:|")
    lines.append(f"| spots analyzed | {per_spot.shape[0]} |")
    lines.append(f"| cancer spots analyzed | {cancer.shape[0]} |")
    lines.append(f"| cancer slides | {cancer['sample_id'].nunique()} |")
    lines.append(f"| median t00 same-label KNN z | {fmt(cancer_coh['same_knn_z'].median())} |")
    lines.append(f"| median t00 permutation p | {fmt(cancer_coh['perm_p'].median())} |")
    lines.append(f"| median bin-bootstrap ARI | {fmt(stab[stab['condition'].isin(['PTC','LPTC','ATC'])]['median_ARI_bin_bootstrap'].median())} |")

    lines.append("\n## T00-high territory ecosystem enrichment\n")
    show = eco[eco["scope"].isin(["pooled_cancer_spots", "slide_median_effect"])].copy()
    axes = ["TROP2_raw", "TACSTD2_z", "Macrophage_TAM_z", "TLS_B_z", "CAF_ECM_z", "RAI_thyroid_z", "DM1_like_score", "RAI_8_score"]
    lines.append("| Scope | Axis | delta high-low | d | p |")
    lines.append("|---|---|---:|---:|---:|")
    for axis in axes:
        for scope in ["pooled_cancer_spots", "slide_median_effect"]:
            sub = show[(show["axis"] == axis) & (show["scope"] == scope)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            lines.append(f"| {scope} | {axis} | {fmt(r['delta_high_low'])} | {fmt(r['cohen_d'])} | {fmt(r['mw_p'])} |")

    lines.append("\n## T00 arm direction audit\n")
    lines.append("| Arm | Expected | delta high-low | matches? | d | p |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for r in arms.itertuples(index=False):
        lines.append(f"| {r.arm} | {r.expected} | {fmt(r.delta_high_low)} | {r.direction_matches} | {fmt(r.cohen_d)} | {fmt(r.mw_p)} |")

    lines.append("\n## Pathology-to-ST-CNV bridge\n")
    lines.append("| Mode | Target | Metric | Value | n |")
    lines.append("|---|---|---|---:|---:|")
    order = ["UNI_image_only", "cell_morph_only", "stage_text_only", "UNI_plus_cell_morph", "UNI_plus_stage_text", "UNI_cell_stage"]
    targets = ["t00_clone_high", "infercnv_t00_signature", "aneuploidy_clone_high", "infercnv_aneuploidy_score"]
    for target in targets:
        for mode in order:
            sub = path_metrics[(path_metrics["mode"] == mode) & (path_metrics["target"] == target)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            lines.append(f"| {mode} | {target} | {r['metric']} | {fmt(r['value'])} | {int(r['n'])} |")

    lines.append("\n## Paper interpretation\n")
    lines.append("The paper-worthy version is: bulk TCGA/cBio ARMDRIVER defines the true CNV residual class; spatial transcriptomics projects a gene-order inferCNV-like territory layer; T00-high territories are enriched for TROP2/TAM/TLS/RAI programs. The pathology bridge is a separate gate: if image-only AUROC is modest, present it as a boundary showing that ST-derived CNV territories are not reliably visible from H&E alone.\n")

    lines.append("## Output files\n")
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    per_spot = compute_infercnv_lite()
    eco = ecosystem_enrichment(per_spot)
    arms = arm_direction_audit(per_spot)
    path_metrics = pathology_bridge(per_spot)
    write_report(per_spot, eco, path_metrics, arms)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
