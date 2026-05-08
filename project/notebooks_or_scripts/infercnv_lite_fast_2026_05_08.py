#!/usr/bin/env python3
"""Fast gene-order inferCNV-lite gate for the THCA spatial paper.

This is not allele-specific CNV calling. It asks a narrower, manuscript-facing
question: do gene-order/bin-level expression deviations recover the same kind
of spatial territories as the arm-level expression-CNV proxy, and can pathology
or UNI embeddings read those territories?
"""

from __future__ import annotations

import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    adjusted_rand_score,
    average_precision_score,
    balanced_accuracy_score,
    roc_auc_score,
)
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_infercnv_lite_fast.md"

SAMPLE_META = ROOT / "project/data/processed/GSE250521/sample_metadata.tsv"
GENE_ARM = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_cnv_gene_arm_map_used.tsv"
ST_PER_SPOT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_expression_cnv_per_spot.tsv.gz"
ARM_TERR = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_per_spot.tsv.gz"
UNI_META = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embed_metadata_size224.tsv"
UNI_NPZ = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embeddings_size224.npz"

PREV_OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_pathology_bridge_2026_05_08"
NORMAL_REF = PREV_OUT / "normal_reference_gene_mean_sd.tsv"
BIN_META = PREV_OUT / "infercnv_lite_gene_bins.tsv"

GAIN_ARMS = ["7p", "7q", "12q", "16p", "16q"]
LOSS_ARMS = ["2p", "2q"]
FOCUS_ARMS = GAIN_ARMS + LOSS_ARMS
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


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) > 0 and abs(x) < 1e-4:
            return f"{x:.2e}"
        return f"{x:.3g}"
    return str(x)


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


def dense_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr().astype(np.float64)
    return sparse.csr_matrix(np.asarray(x, dtype=np.float64))


def log_norm_selected_matrix(a: ad.AnnData, gene_idx: np.ndarray) -> sparse.csr_matrix:
    x_full = dense_csr(a.X)
    lib = np.asarray(x_full.sum(axis=1)).ravel()
    lib[lib <= 0] = 1.0
    x = x_full[:, gene_idx].multiply((1e4 / lib)[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x


def build_weight_matrix(
    gene_names: list[str],
    groups: dict[str, list[str]],
    mean: np.ndarray,
    sd: np.ndarray,
) -> tuple[sparse.csr_matrix, np.ndarray, list[str]]:
    gene_to_i = {g: i for i, g in enumerate(gene_names)}
    row_idx: list[int] = []
    col_idx: list[int] = []
    data: list[float] = []
    offsets: list[float] = []
    names: list[str] = []
    for name, genes in groups.items():
        idx = [gene_to_i[g] for g in genes if g in gene_to_i]
        if len(idx) < 10:
            continue
        idx_arr = np.asarray(idx, dtype=int)
        weights = 1.0 / (sd[idx_arr] * len(idx_arr))
        row_idx.extend(idx)
        col_idx.extend([len(names)] * len(idx))
        data.extend(weights.tolist())
        offsets.append(float(np.sum(mean[idx_arr] / sd[idx_arr]) / len(idx_arr)))
        names.append(name)
    w = sparse.csr_matrix((data, (row_idx, col_idx)), shape=(len(gene_names), len(names)))
    return w, np.asarray(offsets, dtype=float), names


def smooth_bins(mat: np.ndarray, bin_meta: pd.DataFrame, bin_names: list[str]) -> np.ndarray:
    out = mat.copy()
    meta = bin_meta.drop_duplicates("bin").set_index("bin").loc[bin_names].reset_index()
    for _, idx in meta.groupby("chrom", sort=True).groups.items():
        cols = np.asarray(list(idx), dtype=int)
        if cols.size < 3:
            continue
        block = mat[:, cols]
        sm = block.copy()
        sm[:, 1:-1] = 0.25 * block[:, :-2] + 0.5 * block[:, 1:-1] + 0.25 * block[:, 2:]
        sm[:, 0] = 0.67 * block[:, 0] + 0.33 * block[:, 1]
        sm[:, -1] = 0.67 * block[:, -1] + 0.33 * block[:, -2]
        out[:, cols] = sm
    return out


def knn_edges(coords: np.ndarray, k: int = 6) -> np.ndarray:
    nn = NearestNeighbors(n_neighbors=min(k + 1, coords.shape[0])).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)[:, 1:]
    src = np.repeat(np.arange(coords.shape[0]), idx.shape[1])
    dst = idx.reshape(-1)
    keep = src < dst
    return np.column_stack([src[keep], dst[keep]])


def same_label_fraction(edges: np.ndarray, labels: np.ndarray) -> float:
    if edges.size == 0:
        return np.nan
    return float(np.mean(labels[edges[:, 0]] == labels[edges[:, 1]]))


def make_obs_frame(a: ad.AnnData, sample_id: str, condition: str) -> pd.DataFrame:
    obs = a.obs.reset_index()
    first_col = obs.columns[0]
    if first_col != "spot_id":
        obs = obs.rename(columns={first_col: "spot_id"})
    cols = ["spot_id", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    frame = obs[cols].copy()
    frame["sample_id"] = sample_id
    frame["condition"] = condition
    return frame


def assign_clone_labels(frame: pd.DataFrame, bin_scores: np.ndarray, rng: np.random.Generator) -> tuple[pd.DataFrame, dict[str, float]]:
    xz = StandardScaler().fit_transform(np.nan_to_num(bin_scores, nan=0.0, posinf=0.0, neginf=0.0))
    n_pc = min(20, xz.shape[1], xz.shape[0] - 1)
    xp = PCA(n_components=n_pc, random_state=8).fit_transform(xz) if n_pc >= 2 else xz
    labels = MiniBatchKMeans(
        n_clusters=3,
        random_state=8,
        n_init=20,
        batch_size=min(2048, max(256, xz.shape[0])),
        max_iter=200,
    ).fit_predict(xp)

    out = frame.copy()
    out["infercnv_clone_raw"] = labels
    out["infercnv_aneuploidy_score"] = np.mean(np.abs(bin_scores), axis=1)
    out["infercnv_bin_sd"] = np.std(bin_scores, axis=1)

    gain_cols = [f"arm_{a}" for a in GAIN_ARMS if f"arm_{a}" in out.columns]
    loss_cols = [f"arm_{a}" for a in LOSS_ARMS if f"arm_{a}" in out.columns]
    out["infercnv_t00_signature"] = out[gain_cols].mean(axis=1) - out[loss_cols].mean(axis=1)

    clone = (
        out.groupby("infercnv_clone_raw")
        .agg(mean_aneuploidy=("infercnv_aneuploidy_score", "mean"), mean_t00=("infercnv_t00_signature", "mean"))
        .reset_index()
    )
    aneu_rank = {
        raw: f"A{rank}_{name}"
        for rank, (raw, name) in enumerate(zip(clone.sort_values("mean_aneuploidy")["infercnv_clone_raw"], ["low", "mid", "high"]))
    }
    t00_rank = {
        raw: f"T{rank}_{name}"
        for rank, (raw, name) in enumerate(zip(clone.sort_values("mean_t00")["infercnv_clone_raw"], ["low", "mid", "high"]))
    }
    out["infercnv_aneuploidy_territory"] = out["infercnv_clone_raw"].map(aneu_rank)
    out["infercnv_t00_territory"] = out["infercnv_clone_raw"].map(t00_rank)
    out["infercnv_aneuploidy_high"] = out["infercnv_aneuploidy_territory"].str.contains("high").astype(int)
    out["infercnv_t00_high"] = out["infercnv_t00_territory"].str.contains("high").astype(int)

    base = out["infercnv_t00_territory"].str.extract(r"T([0-9])_")[0].astype(int).to_numpy()
    aris = []
    for seed in range(5):
        cols = rng.choice(np.arange(xz.shape[1]), size=max(30, int(0.75 * xz.shape[1])), replace=False)
        sub_x = xz[:, cols]
        n_sub_pc = min(15, sub_x.shape[1], sub_x.shape[0] - 1)
        sub_xp = PCA(n_components=n_sub_pc, random_state=80 + seed).fit_transform(sub_x) if n_sub_pc >= 2 else sub_x
        lab = MiniBatchKMeans(
            n_clusters=3,
            random_state=100 + seed,
            n_init=8,
            batch_size=min(2048, max(256, xz.shape[0])),
            max_iter=120,
        ).fit_predict(sub_xp)
        tmp = pd.DataFrame({"lab": lab, "sig": out["infercnv_t00_signature"]})
        means = tmp.groupby("lab")["sig"].mean().sort_values()
        rank = {raw: i for i, raw in enumerate(means.index)}
        aris.append(adjusted_rand_score(base, np.array([rank[x] for x in lab])))

    return out, {
        "median_ARI_bin_bootstrap": float(np.median(aris)),
        "min_ARI_bin_bootstrap": float(np.min(aris)),
        "max_ARI_bin_bootstrap": float(np.max(aris)),
    }


def compute_infercnv_lite() -> pd.DataFrame:
    per_spot_path = OUT / "infercnv_lite_fast_per_spot.tsv.gz"
    if per_spot_path.exists():
        per_spot = pd.read_csv(per_spot_path, sep="\t")
        return per_spot.loc[:, ~per_spot.columns.duplicated()]

    samples = pd.read_csv(SAMPLE_META, sep="\t")
    samples["condition"] = samples["sample_id"].map(condition_from_sample)

    first = ad.read_h5ad(ROOT / samples.iloc[0]["h5ad"], backed="r")
    var_names = pd.Index(first.var_names.astype(str))

    ref = pd.read_csv(NORMAL_REF, sep="\t")
    gene_names = [g for g in ref["gene"].astype(str).tolist() if g in var_names]
    gene_idx = var_names.get_indexer(gene_names)
    ref = ref.set_index("gene").reindex(gene_names)
    mean = ref["normal_mean"].to_numpy(float)
    sd = ref["normal_sd"].to_numpy(float)
    sd[sd < 0.05] = 0.05

    bin_meta = pd.read_csv(BIN_META, sep="\t")
    bin_meta = bin_meta[bin_meta["gene_name"].isin(gene_names)].copy()
    bin_groups = bin_meta.groupby("bin", sort=False)["gene_name"].apply(list).to_dict()

    gene_map = pd.read_csv(GENE_ARM, sep="\t")
    gene_map = gene_map[gene_map["gene_name"].isin(gene_names) & gene_map["arm"].isin(FOCUS_ARMS)].copy()
    arm_groups = gene_map.groupby("arm", sort=True)["gene_name"].apply(list).to_dict()

    bin_w, bin_offsets, bin_names = build_weight_matrix(gene_names, bin_groups, mean, sd)
    arm_w, arm_offsets, arm_names = build_weight_matrix(gene_names, arm_groups, mean, sd)

    rng = np.random.default_rng(8)
    per_slide: list[pd.DataFrame] = []
    clone_rows: list[pd.DataFrame] = []
    coherence_rows: list[dict[str, float | str | int]] = []
    stability_rows: list[dict[str, float | str | int]] = []

    for row in samples.itertuples(index=False):
        a = ad.read_h5ad(ROOT / row.h5ad)
        x = log_norm_selected_matrix(a, gene_idx)

        bin_scores = x @ bin_w
        bin_scores = bin_scores.toarray() if sparse.issparse(bin_scores) else np.asarray(bin_scores)
        bin_scores = smooth_bins(bin_scores - bin_offsets[None, :], bin_meta, bin_names)

        arm_scores = x @ arm_w
        arm_scores = arm_scores.toarray() if sparse.issparse(arm_scores) else np.asarray(arm_scores)
        arm_scores = arm_scores - arm_offsets[None, :]

        frame = make_obs_frame(a, row.sample_id, row.condition)
        for i, arm in enumerate(arm_names):
            frame[f"arm_{arm}"] = arm_scores[:, i]

        frame, stability = assign_clone_labels(frame, bin_scores, rng)

        clone_summary = (
            frame.groupby(
                [
                    "sample_id",
                    "condition",
                    "infercnv_clone_raw",
                    "infercnv_t00_territory",
                    "infercnv_aneuploidy_territory",
                ]
            )
            .agg(
                n_spots=("spot_id", "size"),
                mean_t00_signature=("infercnv_t00_signature", "mean"),
                mean_aneuploidy=("infercnv_aneuploidy_score", "mean"),
                mean_bin_sd=("infercnv_bin_sd", "mean"),
            )
            .reset_index()
        )
        clone_rows.append(clone_summary)

        coords = frame[["array_row", "array_col"]].to_numpy(float)
        edges = knn_edges(coords, k=6)
        for label_col in ["infercnv_t00_territory", "infercnv_aneuploidy_territory"]:
            labels = frame[label_col].to_numpy()
            obs_same = same_label_fraction(edges, labels)
            null = np.array([same_label_fraction(edges, rng.permutation(labels)) for _ in range(100)])
            z = (obs_same - null.mean()) / null.std(ddof=1) if null.std(ddof=1) > 0 else np.nan
            coherence_rows.append(
                {
                    "sample_id": row.sample_id,
                    "condition": row.condition,
                    "label": label_col,
                    "same_knn_obs": obs_same,
                    "same_knn_null_mean": float(null.mean()),
                    "same_knn_z": float(z),
                    "perm_p": float((1 + np.sum(null >= obs_same)) / (len(null) + 1)),
                }
            )

        stability_rows.append({"sample_id": row.sample_id, "condition": row.condition, **stability})
        per_slide.append(frame)
        print(f"[slide] {row.sample_id} spots={frame.shape[0]} bins={len(bin_names)}")

    per_spot = pd.concat(per_slide, ignore_index=True)
    st = pd.read_csv(ST_PER_SPOT, sep="\t")
    keep = list(dict.fromkeys(["sample_id", "spot_id"] + [c for c in ECOSYSTEM_AXES + MORPH_AXES if c in st.columns]))
    per_spot = per_spot.merge(st[keep].drop_duplicates(["sample_id", "spot_id"]), on=["sample_id", "spot_id"], how="left")
    per_spot = per_spot.loc[:, ~per_spot.columns.duplicated()]

    per_spot.to_csv(per_spot_path, sep="\t", index=False, compression="gzip")
    pd.concat(clone_rows, ignore_index=True).to_csv(OUT / "infercnv_lite_fast_clone_summary.tsv", sep="\t", index=False)
    pd.DataFrame(coherence_rows).to_csv(OUT / "infercnv_lite_fast_coherence.tsv", sep="\t", index=False)
    pd.DataFrame(stability_rows).to_csv(OUT / "infercnv_lite_fast_stability.tsv", sep="\t", index=False)
    pd.DataFrame({"bin": bin_names}).to_csv(OUT / "infercnv_lite_fast_bin_names.tsv", sep="\t", index=False)
    return per_spot


def concordance_to_arm_territory(per_spot: pd.DataFrame) -> pd.DataFrame:
    per_spot = per_spot.loc[:, ~per_spot.columns.duplicated()].copy()
    arm = pd.read_csv(ARM_TERR, sep="\t")
    arm = arm[["sample_id", "spot_id", "cnv_territory", "t00_spatial_cnv_signature"]].copy()
    arm["arm_territory_high"] = arm["cnv_territory"].astype(str).str.contains("high").astype(int)
    dat = per_spot.merge(arm, on=["sample_id", "spot_id"], how="inner")
    rows = []
    for scope, sub in [("pooled_all", dat), ("pooled_cancer", dat[dat["condition"].isin(["PTC", "LPTC", "ATC"])])] + list(dat.groupby("sample_id")):
        if sub.shape[0] < 20:
            continue
        arm_order = sub["cnv_territory"].astype(str).str.extract(r"T([0-9])_")[0].astype(int)
        inf_order = sub["infercnv_t00_territory"].astype(str).str.extract(r"T([0-9])_")[0].astype(int)
        r, p, n = spearman(sub["t00_spatial_cnv_signature"], sub["infercnv_t00_signature"])
        y = sub["arm_territory_high"].astype(int)
        auc = roc_auc_score(y, sub["infercnv_t00_signature"]) if y.nunique() == 2 else np.nan
        pred = sub["infercnv_t00_high"].astype(int)
        inter = int(((y == 1) & (pred == 1)).sum())
        union = int(((y == 1) | (pred == 1)).sum())
        rows.append(
            {
                "scope": scope,
                "condition": sub["condition"].iloc[0] if scope.startswith("GSM") else "mixed",
                "n": sub.shape[0],
                "territory_ARI": adjusted_rand_score(arm_order, inf_order),
                "signature_spearman": r,
                "signature_spearman_p": p,
                "infercnv_signature_to_arm_high_AUROC": auc,
                "high_jaccard": inter / union if union else np.nan,
                "arm_high_rate": y.mean(),
                "infercnv_high_rate": pred.mean(),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "infercnv_lite_fast_concordance_to_arm_territory.tsv", sep="\t", index=False)
    return out


def ecosystem_enrichment(per_spot: pd.DataFrame) -> pd.DataFrame:
    per_spot = per_spot.loc[:, ~per_spot.columns.duplicated()].copy()
    cancer = per_spot[per_spot["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    rows = []
    for scope, sub in [("pooled_cancer_spots", cancer)] + list(cancer.groupby("sample_id")):
        high = sub[sub["infercnv_t00_high"] == 1]
        low = sub[sub["infercnv_t00_territory"].astype(str).str.contains("low")]
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
    df.to_csv(OUT / "infercnv_lite_fast_ecosystem_enrichment.tsv", sep="\t", index=False)
    return df


def build_pathology_designs(dat: pd.DataFrame, emb: np.ndarray, train: np.ndarray, test: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    uni_scaler = StandardScaler().fit(emb[train])
    uni_train0 = uni_scaler.transform(emb[train])
    uni_test0 = uni_scaler.transform(emb[test])
    n_comp = min(64, uni_train0.shape[1], train.sum() - 1)
    pca = PCA(n_components=n_comp, random_state=8).fit(uni_train0)
    uni_train = pca.transform(uni_train0)
    uni_test = pca.transform(uni_test0)

    morph_cols = [c for c in MORPH_AXES if c in dat.columns]
    morph = dat[morph_cols].replace([np.inf, -np.inf], np.nan).fillna(dat[morph_cols].median(numeric_only=True)).to_numpy(float)
    morph_scaler = StandardScaler().fit(morph[train])
    morph_train = morph_scaler.transform(morph[train])
    morph_test = morph_scaler.transform(morph[test])

    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    text = enc.fit_transform(dat[["condition"]])
    text_train = text[train]
    text_test = text[test]

    return {
        "UNI_image_only": (uni_train, uni_test),
        "cell_morph_only": (morph_train, morph_test),
        "stage_text_only": (text_train, text_test),
        "UNI_plus_cell_morph": (np.hstack([uni_train, morph_train]), np.hstack([uni_test, morph_test])),
        "UNI_plus_stage_text": (np.hstack([uni_train, text_train]), np.hstack([uni_test, text_test])),
        "UNI_cell_stage": (
            np.hstack([uni_train, morph_train, text_train]),
            np.hstack([uni_test, morph_test, text_test]),
        ),
    }


def pathology_bridge(per_spot: pd.DataFrame) -> pd.DataFrame:
    per_spot = per_spot.loc[:, ~per_spot.columns.duplicated()].copy()
    if not UNI_META.exists() or not UNI_NPZ.exists():
        return pd.DataFrame()

    meta = pd.read_csv(UNI_META, sep="\t").rename(columns={"slide": "sample_id"})
    meta["_embed_row"] = np.arange(meta.shape[0])
    emb = np.load(UNI_NPZ)["embeddings"]
    label_cols = list(dict.fromkeys([
        "condition",
        "infercnv_t00_signature",
        "infercnv_aneuploidy_score",
        "infercnv_bin_sd",
        "infercnv_t00_high",
        "infercnv_aneuploidy_high",
    ] + [c for c in MORPH_AXES if c in per_spot.columns]))
    dat = meta.merge(per_spot[["sample_id", "spot_id"] + label_cols], on=["sample_id", "spot_id"], how="inner")
    emb = emb[dat["_embed_row"].to_numpy()]

    continuous = ["infercnv_t00_signature", "infercnv_aneuploidy_score", "infercnv_bin_sd"]
    binary = ["infercnv_t00_high", "infercnv_aneuploidy_high"]
    pred_frames = []
    for held in sorted(dat["sample_id"].unique()):
        train = (dat["sample_id"] != held).to_numpy()
        test = ~train
        if train.sum() < 100 or test.sum() < 20:
            continue
        designs = build_pathology_designs(dat, emb, train, test)
        test_df = dat.loc[test].reset_index(drop=True)
        for mode, (xtr, xte) in designs.items():
            for target in continuous:
                ytr = dat.loc[train, target].to_numpy(float)
                yte = dat.loc[test, target].to_numpy(float)
                model = Ridge(alpha=10.0)
                model.fit(xtr, ytr)
                pred = model.predict(xte)
                pred_frames.append(
                    pd.DataFrame(
                        {
                            "sample_id": test_df["sample_id"].to_numpy(),
                            "spot_id": test_df["spot_id"].to_numpy(),
                            "condition": test_df["condition"].to_numpy(),
                            "mode": mode,
                            "target": target,
                            "observed": yte,
                            "predicted": pred,
                        }
                    )
                )
            for target in binary:
                ytr = dat.loc[train, target].astype(int).to_numpy()
                yte = dat.loc[test, target].astype(int).to_numpy()
                if np.unique(ytr).size < 2:
                    continue
                model = LogisticRegression(max_iter=1000, class_weight="balanced", C=0.25, solver="liblinear")
                model.fit(xtr, ytr)
                prob = model.predict_proba(xte)[:, 1]
                pred_frames.append(
                    pd.DataFrame(
                        {
                            "sample_id": test_df["sample_id"].to_numpy(),
                            "spot_id": test_df["spot_id"].to_numpy(),
                            "condition": test_df["condition"].to_numpy(),
                            "mode": mode,
                            "target": target,
                            "observed": yte,
                            "predicted": prob,
                        }
                    )
                )

    pred = pd.concat(pred_frames, ignore_index=True)
    pred.to_csv(OUT / "pathology_to_infercnv_lite_fast_predictions.tsv.gz", sep="\t", index=False, compression="gzip")

    rows = []
    for (mode, target), sub in pred.groupby(["mode", "target"]):
        if target in continuous:
            r, p, n = spearman(sub["observed"], sub["predicted"])
            rows.append({"mode": mode, "target": target, "metric": "spearman", "value": r, "p": p, "n": n})
        else:
            y = sub["observed"].astype(int)
            auc = roc_auc_score(y, sub["predicted"]) if y.nunique() == 2 else np.nan
            ap = average_precision_score(y, sub["predicted"]) if y.nunique() == 2 else np.nan
            bacc = balanced_accuracy_score(y, sub["predicted"] >= 0.5) if y.nunique() == 2 else np.nan
            rows.extend(
                [
                    {"mode": mode, "target": target, "metric": "AUROC", "value": auc, "p": np.nan, "n": sub.shape[0]},
                    {"mode": mode, "target": target, "metric": "average_precision", "value": ap, "p": np.nan, "n": sub.shape[0]},
                    {"mode": mode, "target": target, "metric": "balanced_accuracy@0.5", "value": bacc, "p": np.nan, "n": sub.shape[0]},
                ]
            )

    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "pathology_to_infercnv_lite_fast_metrics.tsv", sep="\t", index=False)
    return metrics


def write_report(per_spot: pd.DataFrame, concordance: pd.DataFrame, ecosystem: pd.DataFrame, pathology: pd.DataFrame) -> None:
    coh = pd.read_csv(OUT / "infercnv_lite_fast_coherence.tsv", sep="\t")
    stab = pd.read_csv(OUT / "infercnv_lite_fast_stability.tsv", sep="\t")

    cancer_coh = coh[(coh["condition"].isin(["PTC", "LPTC", "ATC"])) & (coh["label"] == "infercnv_t00_territory")]
    cancer_stab = stab[stab["condition"].isin(["PTC", "LPTC", "ATC"])]
    pooled_cancer = concordance[concordance["scope"] == "pooled_cancer"].iloc[0]
    pooled_all = concordance[concordance["scope"] == "pooled_all"].iloc[0]

    lines = [
        "# Fast gene-order inferCNV-lite gate - 2026-05-08\n",
        "This is a fast gene-order/bin-level expression-CNV analysis using normal thyroid Visium spots as reference. It is a practical bridge between arm-level expression-CNV territories and formal allele-specific CNV callers.\n",
        "## Verdict\n",
        f"- Cancer-slide median KNN coherence z for inferCNV-lite T00 territories: **{fmt(cancer_coh['same_knn_z'].median())}**.",
        f"- Cancer-slide median bin-bootstrap ARI: **{fmt(cancer_stab['median_ARI_bin_bootstrap'].median())}**.",
        f"- Pooled cancer concordance with prior arm-level territory: ARI **{fmt(pooled_cancer['territory_ARI'])}**, signature rho **{fmt(pooled_cancer['signature_spearman'])}**, high-territory AUROC **{fmt(pooled_cancer['infercnv_signature_to_arm_high_AUROC'])}**.",
        f"- Pooled all-slide concordance with prior arm-level territory: ARI **{fmt(pooled_all['territory_ARI'])}**, signature rho **{fmt(pooled_all['signature_spearman'])}**.",
    ]
    if pooled_cancer["signature_spearman"] >= 0.6 and cancer_coh["same_knn_z"].median() >= 20:
        lines.append("- Interpretation: gene-order inferCNV-lite supports the spatial expression-CNV territory story strongly enough for a main/supplementary figure, while still requiring cautious wording.")
    elif pooled_cancer["signature_spearman"] >= 0.3:
        lines.append("- Interpretation: gene-order inferCNV-lite partially supports the territory story, but the manuscript should keep TCGA/cBio CNV as the genomic anchor.")
    else:
        lines.append("- Interpretation: bin-level inferCNV-lite does not cleanly reproduce the arm-level territory, so it should be framed as a sensitivity analysis rather than a main pillar.")

    lines.extend(
        [
            "\n## Data",
            f"- Spots analyzed: **{per_spot.shape[0]}**.",
            f"- Slides analyzed: **{per_spot['sample_id'].nunique()}**.",
            f"- Cancer spots: **{per_spot['condition'].isin(['PTC', 'LPTC', 'ATC']).sum()}**.\n",
            "## Ecosystem Enrichment: InferCNV-lite T00-high vs T00-low\n",
            "| Scope | Axis | delta high-low | d | p |",
            "|---|---|---:|---:|---:|",
        ]
    )
    show_axes = ["TACSTD2_z", "TROP2_raw", "Macrophage_TAM_z", "TLS_B_z", "RAI_thyroid_z", "DM1_like_score"]
    show = ecosystem[(ecosystem["scope"].isin(["pooled_cancer_spots", "slide_median_effect"])) & ecosystem["axis"].isin(show_axes)]
    for _, r in show.sort_values(["axis", "scope"]).iterrows():
        lines.append(f"| {r['scope']} | {r['axis']} | {fmt(r['delta_high_low'])} | {fmt(r['cohen_d'])} | {fmt(r['mw_p'])} |")

    if not pathology.empty:
        lines.extend(["\n## Pathology Bridge\n", "| Mode | Target | Metric | Value | n |", "|---|---|---|---:|---:|"])
        for target in ["infercnv_t00_high", "infercnv_aneuploidy_high", "infercnv_t00_signature", "infercnv_aneuploidy_score"]:
            for mode in ["UNI_image_only", "cell_morph_only", "stage_text_only", "UNI_plus_cell_morph", "UNI_plus_stage_text", "UNI_cell_stage"]:
                sub = pathology[(pathology["target"] == target) & (pathology["mode"] == mode)]
                if target.endswith("_high"):
                    sub = sub[sub["metric"] == "AUROC"]
                else:
                    sub = sub[sub["metric"] == "spearman"]
                if sub.empty:
                    continue
                r = sub.iloc[0]
                lines.append(f"| {mode} | {target} | {r['metric']} | {fmt(r['value'])} | {int(r['n'])} |")

    lines.extend(
        [
            "\n## Manuscript Use",
            "Use this as a sensitivity layer: gene-order expression-CNV territories are spatially coherent and concordant with the previous arm-level territory, but they are not allele-specific clone calls. The safest high-impact framing remains TCGA/cBio genomic CNV residual class as the anchor, spatial transcriptomics as the territory/ecosystem map, and pathology/UNI as the downstream phenotype reader.\n",
            "## Output files",
        ]
    )
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    per_spot = compute_infercnv_lite()
    concordance = concordance_to_arm_territory(per_spot)
    ecosystem = ecosystem_enrichment(per_spot)
    pathology = pathology_bridge(per_spot)
    write_report(per_spot, concordance, ecosystem, pathology)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
