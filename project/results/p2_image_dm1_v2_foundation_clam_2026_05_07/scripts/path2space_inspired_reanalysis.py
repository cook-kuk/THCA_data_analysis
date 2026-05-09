#!/usr/bin/env python3
"""Path2Space-inspired reanalysis of GSE250521 spatial pathology axes.

This is a measured-ST + H&E-embedding reanalysis, not a full Path2Space
replication. It ports the reusable ideas from Shulman et al. Cell 2026:
spot-neighborhood smoothing, leave-slide-out image-to-expression prediction,
SPAND-like spatial heterogeneity, and patient/slide-level cluster composition.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler


REPO = Path(os.environ.get("THCA_REPO_ROOT", "/home/seungho/personal/THCA_data_analysis"))
P2_ROOT = REPO / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
PHASE1 = P2_ROOT / "phase1_gse250521"
SCORES = REPO / "project/results/01_spatial_score/all_spots_scored.tsv.gz"
GSE_PROC = REPO / "project/data/processed/GSE250521"
OUT = P2_ROOT / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
MAPK_OUT = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]
BASE_TARGETS = [
    "DM1_like_score",
    "TDS_like_score",
    "Epithelial_score",
    "CAF_ECM_score",
    "Hypoxia_score",
    "Proliferation_score",
    "EMT_score",
    "MAPK_output_score",
]
CLUSTER_FEATURES = [
    "pred_DM1_like_score_smooth8",
    "pred_TDS_like_score_smooth8",
    "pred_Epithelial_score_smooth8",
    "pred_CAF_ECM_score_smooth8",
    "pred_Hypoxia_score_smooth8",
    "pred_Proliferation_score_smooth8",
    "pred_MAPK_output_score_smooth8",
]
STAGE_ORDER = ["PT", "PTC", "LPTC", "ATC"]
N_PCS = 64
RIDGE_ALPHA = 100.0
SMOOTH_K = 8
RANDOM_SEED = 9001


def get_expr(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    col = adata.X[:, adata.var_names.get_loc(gene)]
    if sp.issparse(col):
        return col.toarray().ravel().astype(float)
    return np.asarray(col).ravel().astype(float)


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


def add_mapk_score(scores: pd.DataFrame) -> pd.DataFrame:
    meta = pd.read_csv(GSE_PROC / "sample_metadata.tsv", sep="\t")
    chunks = []
    for row in meta.itertuples(index=False):
        import scanpy as sc

        raw = REPO / row.h5ad
        scored = raw.parent / (raw.stem.replace(".raw", "") + ".scored.h5ad")
        use_h5 = scored if scored.exists() else raw
        adata = sc.read_h5ad(use_h5)
        if "log1p" not in (adata.uns or {}):
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        mapk, n_mapk = module_z(adata, MAPK_OUT)
        panel, n_panel = module_z(adata, PANEL_8)
        chunks.append(
            pd.DataFrame(
                {
                    "sample_id": row.sample_id,
                    "spot_id": adata.obs.index.astype(str).to_numpy(),
                    "MAPK_output_score": mapk,
                    "Panel8_module_z": panel,
                    "n_MAPK_genes": n_mapk,
                    "n_panel_genes": n_panel,
                }
            )
        )
    mapk_df = pd.concat(chunks, ignore_index=True)
    out = scores.merge(mapk_df, on=["sample_id", "spot_id"], how="left")
    return out


def spatial_smooth_by_slide(df: pd.DataFrame, cols: list[str], k: int = SMOOTH_K) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        out[f"{col}_smooth{k}"] = np.nan
    for slide, sub in out.groupby("sample_id", sort=False):
        idx = sub.index.to_numpy()
        coords = sub[["array_row", "array_col"]].astype(float).to_numpy()
        if len(idx) < 4:
            continue
        tree = cKDTree(coords)
        kk = min(k + 1, len(idx))
        _, neigh = tree.query(coords, k=kk)
        if neigh.ndim == 1:
            neigh = neigh[:, None]
        for col in cols:
            vals = pd.to_numeric(sub[col], errors="coerce").to_numpy(dtype=float)
            out.loc[idx, f"{col}_smooth{k}"] = np.nanmean(vals[neigh], axis=1)
    return out


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def predict_loso(x: np.ndarray, y: np.ndarray, slides: np.ndarray) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    for held in sorted(pd.unique(slides)):
        train = slides != held
        test = ~train
        if train.sum() < 20 or test.sum() == 0:
            continue
        yt = y[train]
        mean = float(np.nanmean(yt))
        sd = float(np.nanstd(yt))
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        scaler = StandardScaler()
        x_train = scaler.fit_transform(x[train])
        x_test = scaler.transform(x[test])
        n_comp = min(N_PCS, x_train.shape[0] - 1, x_train.shape[1])
        pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
        x_train_pc = pca.fit_transform(x_train)
        x_test_pc = pca.transform(x_test)
        model = Ridge(alpha=RIDGE_ALPHA)
        model.fit(x_train_pc, (yt - mean) / sd)
        pred[test] = model.predict(x_test_pc)
    return pred


def run_prediction(df: pd.DataFrame, emb: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pred_base = df[["sample_id", "spot_id", "stage", "array_row", "array_col"]].copy()
    slides = df["sample_id"].astype(str).to_numpy()
    target_defs = []
    for target in BASE_TARGETS:
        for mode, y_col in [("raw", target), (f"smooth{SMOOTH_K}", f"{target}_smooth{SMOOTH_K}")]:
            y_raw = pd.to_numeric(df[y_col], errors="coerce").to_numpy(dtype=float)
            keep = np.isfinite(y_raw) & np.isfinite(emb).all(axis=1)
            y = np.full(len(y_raw), np.nan)
            y[keep] = (y_raw[keep] - np.nanmean(y_raw[keep])) / (np.nanstd(y_raw[keep]) + 1e-9)
            target_defs.append({"target": target, "mode": mode, "y": y, "keep": keep})
            pred_base[f"obs_{target}_{mode}"] = y

    # PCA depends only on the held-out slide, not on the target. Cache each
    # fold transform and fit one Ridge head per target on top of it.
    x_keep = np.isfinite(emb).all(axis=1)
    pc_cache = {}
    for held in sorted(pd.unique(slides)):
        train = (slides != held) & x_keep
        test = (slides == held) & x_keep
        scaler = StandardScaler()
        x_train = scaler.fit_transform(emb[train])
        x_test = scaler.transform(emb[test])
        n_comp = min(N_PCS, x_train.shape[0] - 1, x_train.shape[1])
        pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
        pc_cache[held] = {
            "train_idx": np.where(train)[0],
            "test_idx": np.where(test)[0],
            "x_train_pc": pca.fit_transform(x_train),
            "x_test_pc": pca.transform(x_test),
        }

    for item in target_defs:
        y = item["y"]
        full_pred = np.full(len(y), np.nan)
        for held, cache in pc_cache.items():
            train_idx = cache["train_idx"]
            test_idx = cache["test_idx"]
            train_ok = np.isfinite(y[train_idx])
            if train_ok.sum() < 20:
                continue
            yt = y[train_idx][train_ok]
            mean = float(np.nanmean(yt))
            sd = float(np.nanstd(yt))
            if not np.isfinite(sd) or sd == 0:
                sd = 1.0
            model = Ridge(alpha=RIDGE_ALPHA)
            model.fit(cache["x_train_pc"][train_ok], (yt - mean) / sd)
            full_pred[test_idx] = model.predict(cache["x_test_pc"])
        pred_base[f"pred_{item['target']}_{item['mode']}"] = full_pred

    summary_rows = []
    slide_rows = []
    for item in target_defs:
        target = item["target"]
        mode = item["mode"]
        y = pred_base[f"obs_{target}_{mode}"].to_numpy(dtype=float)
        full_pred = pred_base[f"pred_{target}_{mode}"].to_numpy(dtype=float)
        rho, p, n = safe_spearman(y, full_pred)
        per_slide_rhos = []
        for slide, idx_raw in df.groupby("sample_id").groups.items():
            idx = np.asarray(list(idx_raw))
            srho, sp, sn = safe_spearman(y[idx], full_pred[idx])
            per_slide_rhos.append(srho)
            slide_rows.append(
                {
                    "target": target,
                    "mode": mode,
                    "sample_id": slide,
                    "stage": str(df.loc[idx[0], "stage"]),
                    "n_spots": sn,
                    "rho": srho,
                    "p": sp,
                }
            )
        summary_rows.append(
            {
                "target": target,
                "mode": mode,
                "n_spots": n,
                "overall_spearman_rho": rho,
                "overall_spearman_p": p,
                "median_slide_rho": float(np.nanmedian(per_slide_rhos)),
                "slides_rho_gt_0_2": int(np.nansum(np.asarray(per_slide_rhos) > 0.2)),
                "slides_rho_gt_0_3": int(np.nansum(np.asarray(per_slide_rhos) > 0.3)),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values(["mode", "overall_spearman_rho"], ascending=[True, False])
    slide_summary = pd.DataFrame(slide_rows)
    return summary, slide_summary, pred_base


def moran_i(values: np.ndarray, coords: np.ndarray, k: int = 4) -> float:
    values = np.asarray(values, dtype=float)
    coords = np.asarray(coords, dtype=float)
    keep = np.isfinite(values) & np.isfinite(coords).all(axis=1)
    if keep.sum() < k + 3:
        return np.nan
    vals = values[keep]
    xy = coords[keep]
    z = vals - np.nanmean(vals)
    denom = float(np.sum(z * z))
    if denom <= 1e-12:
        return np.nan
    tree = cKDTree(xy)
    kk = min(k + 1, len(vals))
    _, neigh = tree.query(xy, k=kk)
    if neigh.ndim == 1:
        neigh = neigh[:, None]
    neigh = neigh[:, 1:]
    local = np.nanmean(z[neigh], axis=1)
    return float(np.sum(z * local) / denom)


def rank01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.full(len(x), np.nan)
    keep = np.isfinite(x)
    if keep.sum() < 3:
        return out
    out[keep] = (stats.rankdata(x[keep], method="average") - 0.5) / keep.sum()
    return out


def spand_table(df: pd.DataFrame, preds: pd.DataFrame) -> pd.DataFrame:
    merged = df[["sample_id", "spot_id", "stage", "array_row", "array_col"]].merge(
        preds, on=["sample_id", "spot_id", "stage", "array_row", "array_col"], how="left"
    )
    rows = []
    axes = ["DM1_like_score", "TDS_like_score", "Epithelial_score", "CAF_ECM_score", "MAPK_output_score"]
    for slide, sub in merged.groupby("sample_id", sort=False):
        coords = sub[["array_row", "array_col"]].astype(float).to_numpy()
        for axis in axes:
            for source, col in [
                ("measured_smooth8", f"obs_{axis}_smooth{SMOOTH_K}"),
                ("predicted_from_HE", f"pred_{axis}_smooth{SMOOTH_K}"),
            ]:
                vals = rank01(pd.to_numeric(sub[col], errors="coerce").to_numpy(dtype=float))
                mi = moran_i(vals, coords, k=4)
                mean_val = float(np.nanmean(vals))
                rows.append(
                    {
                        "sample_id": slide,
                        "stage": str(sub["stage"].iloc[0]),
                        "axis": axis,
                        "source": source,
                        "n_spots": int(np.isfinite(vals).sum()),
                        "moran_i": mi,
                        "spand_neg_moran_over_mean": float(-mi / (mean_val + 1e-9)) if np.isfinite(mi) else np.nan,
                        "one_minus_moran": float(1.0 - mi) if np.isfinite(mi) else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def cluster_architecture(preds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    use = preds.copy()
    for feature in CLUSTER_FEATURES:
        old = feature.replace("pred_", "").replace("_smooth8", "_smooth8")
        src = "pred_" + old
        if src in use.columns and feature not in use.columns:
            use[feature] = use[src]
    # The prediction columns are named pred_<target>_smooth8.
    actual_features = [f"pred_{t}_smooth{SMOOTH_K}" for t in [
        "DM1_like_score",
        "TDS_like_score",
        "Epithelial_score",
        "CAF_ECM_score",
        "Hypoxia_score",
        "Proliferation_score",
        "MAPK_output_score",
    ]]
    mat = use[actual_features].apply(pd.to_numeric, errors="coerce")
    keep = np.isfinite(mat.to_numpy()).all(axis=1)
    x = StandardScaler().fit_transform(mat.loc[keep].to_numpy(dtype=float))
    k_values = range(3, 9)
    sil_rows = []
    best_k = 5
    best_sil = -np.inf
    for k in k_values:
        labs = KMeans(n_clusters=k, n_init=50, random_state=RANDOM_SEED).fit_predict(x)
        sil = silhouette_score(x, labs) if len(np.unique(labs)) > 1 else np.nan
        sil_rows.append({"k": k, "silhouette": sil})
        if np.isfinite(sil) and sil > best_sil:
            best_sil = sil
            best_k = k
    labels = KMeans(n_clusters=best_k, n_init=100, random_state=RANDOM_SEED).fit_predict(x)
    per_spot = use.loc[keep, ["sample_id", "spot_id", "stage", "array_row", "array_col"]].copy()
    per_spot["p2s_cluster"] = labels + 1
    for col in actual_features:
        per_spot[col] = mat.loc[keep, col].to_numpy(dtype=float)

    props = (
        per_spot.assign(value=1)
        .pivot_table(index=["sample_id", "stage"], columns="p2s_cluster", values="value", aggfunc="sum", fill_value=0)
        .reset_index()
    )
    cluster_cols = [c for c in props.columns if isinstance(c, (int, np.integer))]
    denom = props[cluster_cols].sum(axis=1).replace(0, np.nan)
    for c in cluster_cols:
        props[f"cluster_{c}_prop"] = props[c] / denom
    props = props.drop(columns=cluster_cols)
    prop_cols = [c for c in props.columns if str(c).endswith("_prop")]

    char = per_spot.groupby("p2s_cluster")[actual_features].mean().reset_index()
    char["n_spots"] = per_spot.groupby("p2s_cluster").size().reindex(char["p2s_cluster"]).to_numpy()
    for stage in STAGE_ORDER:
        char[f"frac_{stage}"] = (
            per_spot.assign(is_stage=(per_spot["stage"] == stage).astype(float))
            .groupby("p2s_cluster")["is_stage"]
            .mean()
            .reindex(char["p2s_cluster"])
            .to_numpy()
        )

    stage_numeric = props["stage"].map({s: i for i, s in enumerate(STAGE_ORDER)}).to_numpy()
    dominant = props[prop_cols].idxmax(axis=1).astype(str).str.extract(r"cluster_(\d+)_prop")[0].astype(int).to_numpy()
    ari = adjusted_rand_score(stage_numeric, dominant)
    stage_eta_rows = []
    for col in prop_cols:
        groups = [props.loc[props["stage"] == stage, col].to_numpy(dtype=float) for stage in STAGE_ORDER]
        groups = [g[np.isfinite(g)] for g in groups if len(g) > 0]
        if len(groups) >= 2:
            f, p = stats.f_oneway(*groups)
        else:
            f, p = np.nan, np.nan
        stage_eta_rows.append({"cluster_prop": col, "stage_anova_f": f, "stage_anova_p": p})
    stage_tests = pd.DataFrame(stage_eta_rows)
    meta = {
        "best_k": int(best_k),
        "best_silhouette": float(best_sil),
        "dominant_cluster_stage_ari": float(ari),
    }
    return per_spot, props, char, {"silhouette": sil_rows, "meta": meta, "stage_tests": stage_eta_rows}


def plot_outputs(summary: pd.DataFrame, preds: pd.DataFrame, spand: pd.DataFrame, props: pd.DataFrame, char: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.5))

    pivot = summary.pivot(index="target", columns="mode", values="overall_spearman_rho").reindex(BASE_TARGETS)
    y = np.arange(len(pivot))
    axes[0, 0].barh(y - 0.18, pivot["raw"], height=0.34, color="#9aa4ad", label="raw target")
    axes[0, 0].barh(y + 0.18, pivot[f"smooth{SMOOTH_K}"], height=0.34, color="#246b7a", label="spatial-smoothed target")
    axes[0, 0].axvline(0, color="#444", lw=0.8)
    axes[0, 0].set_yticks(y)
    axes[0, 0].set_yticklabels(pivot.index)
    axes[0, 0].set_xlabel("Leave-slide-out Spearman rho")
    axes[0, 0].set_title("A. H&E UNI to measured spatial programs")
    axes[0, 0].legend(frameon=False, loc="lower right")

    top = "DM1_like_score"
    sub = preds[["stage", f"obs_{top}_smooth{SMOOTH_K}", f"pred_{top}_smooth{SMOOTH_K}"]].dropna()
    colors = {"PT": "#718096", "PTC": "#2b6cb0", "LPTC": "#b7791f", "ATC": "#9b2c2c"}
    for stage, stage_sub in sub.groupby("stage"):
        axes[0, 1].scatter(
            stage_sub[f"obs_{top}_smooth{SMOOTH_K}"],
            stage_sub[f"pred_{top}_smooth{SMOOTH_K}"],
            s=18,
            alpha=0.65,
            color=colors.get(stage, "#333333"),
            label=stage,
        )
    rho, _, _ = safe_spearman(sub[f"obs_{top}_smooth{SMOOTH_K}"], sub[f"pred_{top}_smooth{SMOOTH_K}"])
    axes[0, 1].set_xlabel("Observed smoothed DM1/RAI axis z")
    axes[0, 1].set_ylabel("Predicted from H&E UNI")
    axes[0, 1].set_title(f"B. Top axis prediction (rho={rho:.3f})")
    axes[0, 1].legend(frameon=False, markerscale=1.5)

    sp = spand[(spand["source"] == "predicted_from_HE") & (spand["axis"].isin(["DM1_like_score", "TDS_like_score", "MAPK_output_score"]))]
    sp_pivot = sp.pivot_table(index=["sample_id", "stage"], columns="axis", values="one_minus_moran").reset_index()
    sp_pivot["stage_rank"] = sp_pivot["stage"].map({s: i for i, s in enumerate(STAGE_ORDER)})
    sp_pivot = sp_pivot.sort_values(["stage_rank", "sample_id"])
    im = axes[1, 0].imshow(sp_pivot[["DM1_like_score", "TDS_like_score", "MAPK_output_score"]].to_numpy(dtype=float), aspect="auto", cmap="viridis")
    axes[1, 0].set_yticks(np.arange(len(sp_pivot)))
    axes[1, 0].set_yticklabels(sp_pivot["sample_id"], fontsize=7)
    axes[1, 0].set_xticks(np.arange(3))
    axes[1, 0].set_xticklabels(["DM1/RAI", "TDS", "MAPK"], rotation=25, ha="right")
    axes[1, 0].set_title("C. Predicted spatial heterogeneity (1 - Moran's I)")
    fig.colorbar(im, ax=axes[1, 0], fraction=0.046, pad=0.04)

    prop_cols = [c for c in props.columns if str(c).endswith("_prop")]
    prop_df = props.copy()
    prop_df["stage_rank"] = prop_df["stage"].map({s: i for i, s in enumerate(STAGE_ORDER)})
    prop_df = prop_df.sort_values(["stage_rank", "sample_id"])
    if len(prop_df) > 2 and len(prop_cols) > 1:
        try:
            order = leaves_list(linkage(prop_df[prop_cols].to_numpy(dtype=float), method="ward"))
            prop_df = prop_df.iloc[order]
        except Exception:
            pass
    im2 = axes[1, 1].imshow(prop_df[prop_cols].to_numpy(dtype=float), aspect="auto", cmap="magma", vmin=0)
    axes[1, 1].set_yticks(np.arange(len(prop_df)))
    axes[1, 1].set_yticklabels(prop_df["sample_id"] + " (" + prop_df["stage"] + ")", fontsize=7)
    axes[1, 1].set_xticks(np.arange(len(prop_cols)))
    axes[1, 1].set_xticklabels([c.replace("_prop", "").replace("cluster_", "C") for c in prop_cols], rotation=0)
    axes[1, 1].set_title("D. Predicted spatial-cluster composition")
    fig.colorbar(im2, ax=axes[1, 1], fraction=0.046, pad=0.04)

    fig.suptitle("Path2Space-inspired GSE250521 reanalysis: smoothing, SPAND, and spatial architecture", y=0.995, fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_inspired_reanalysis.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_inspired_reanalysis.pdf")
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(11, 5.8))
    features = [c for c in char.columns if c.startswith("pred_")]
    heat = char.set_index("p2s_cluster")[features]
    heat = (heat - heat.mean(axis=0)) / (heat.std(axis=0) + 1e-9)
    im = ax.imshow(heat.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-2, vmax=2)
    ax.set_yticks(np.arange(len(heat)))
    ax.set_yticklabels([f"C{int(c)}" for c in heat.index])
    ax.set_xticks(np.arange(len(features)))
    ax.set_xticklabels([f.replace("pred_", "").replace(f"_smooth{SMOOTH_K}", "") for f in features], rotation=35, ha="right")
    ax.set_title("Predicted ST cluster expression programs")
    fig2.colorbar(im, ax=ax, fraction=0.03, pad=0.03)
    fig2.tight_layout()
    fig2.savefig(OUT / "fig_path2space_cluster_programs.png", dpi=220)
    fig2.savefig(OUT / "fig_path2space_cluster_programs.pdf")
    plt.close(fig2)


def write_report(summary: pd.DataFrame, slide_summary: pd.DataFrame, spand: pd.DataFrame, props: pd.DataFrame, char: pd.DataFrame, cluster_info: dict) -> dict:
    raw = summary[summary["mode"] == "raw"].set_index("target")
    smooth = summary[summary["mode"] == f"smooth{SMOOTH_K}"].set_index("target")
    top = smooth.sort_values("overall_spearman_rho", ascending=False).iloc[0]
    dm1_rho = float(smooth.loc["DM1_like_score", "overall_spearman_rho"])
    mapk_rho = float(smooth.loc["MAPK_output_score", "overall_spearman_rho"])
    dm1_gain = dm1_rho - float(raw.loc["DM1_like_score", "overall_spearman_rho"])
    cluster_meta = cluster_info["meta"]

    sp_stage = (
        spand[spand["source"] == "predicted_from_HE"]
        .groupby(["axis", "stage"])["one_minus_moran"]
        .median()
        .reset_index()
        .pivot(index="axis", columns="stage", values="one_minus_moran")
    )

    lines = [
        "# Path2Space-inspired spatial pathology reanalysis",
        "",
        "## Strategy Ported From The Cell 2026 Paper",
        "",
        "- Keep the useful pattern, not the breast-specific claim: H&E embeddings are used as a low-cost surrogate for local spatial RNA programs.",
        f"- Replace raw spot labels with KNN-{SMOOTH_K} spatial-smoothed labels before model evaluation, matching the paper's noise-reduction logic.",
        "- Evaluate by leave-one-slide-out validation; no same-slide leakage.",
        "- Move beyond direct MAPK x Panel spot correlation by adding SPAND-like heterogeneity and spatial-cluster composition.",
        "- Treat this as GSE250521/Paper-2 support or caveat, not as causal Paper-1 mechanism prose.",
        "",
        "## Headline Results",
        "",
        f"- Top smoothed target: **{top.name}**, LOSO Spearman rho = **{top.overall_spearman_rho:.3f}**.",
        f"- DM1/RAI axis: raw rho = {raw.loc['DM1_like_score', 'overall_spearman_rho']:.3f}; smoothed rho = **{dm1_rho:.3f}**; gain = {dm1_gain:+.3f}.",
        f"- MAPK output axis from H&E: smoothed rho = **{mapk_rho:.3f}**.",
        f"- Predicted spatial clusters: best k = **{cluster_meta['best_k']}**, silhouette = {cluster_meta['best_silhouette']:.3f}, dominant-cluster/stage ARI = {cluster_meta['dominant_cluster_stage_ari']:.3f}.",
        "",
        "## Ranked Smoothed Targets",
        "",
        "| rank | target | rho | p | median slide rho | slides rho>0.2 |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    ranked = summary[summary["mode"] == f"smooth{SMOOTH_K}"].sort_values("overall_spearman_rho", ascending=False)
    for i, row in enumerate(ranked.itertuples(index=False), 1):
        lines.append(
            f"| {i} | {row.target} | {row.overall_spearman_rho:.3f} | {row.overall_spearman_p:.2e} | "
            f"{row.median_slide_rho:.3f} | {int(row.slides_rho_gt_0_2)}/16 |"
        )
    lines += [
        "",
        "## Predicted SPAND-Like Heterogeneity",
        "",
        "Values are `1 - Moran's I` on rank-normalized predicted scores, so higher values mean more local intermixing / less spatial autocorrelation.",
        "",
        sp_stage.reindex([a for a in ["DM1_like_score", "TDS_like_score", "MAPK_output_score", "Epithelial_score", "CAF_ECM_score"] if a in sp_stage.index])
        .reindex(columns=[s for s in STAGE_ORDER if s in sp_stage.columns])
        .round(3)
        .to_markdown(),
        "",
        "## Boundary",
        "",
        "This reanalysis improves the pathology-to-spatial-RNA framing but does not rescue the direct same-spot MAPK x Panel anti-correlation closed in v15-v18. The safe use is: spatial H&E carries a weak-to-moderate DM1/RAI signal and can support Paper 2 image-to-DM1 strategy; GSE250521 remains caveat-only for Paper 1 MAPK-silencing mechanism.",
    ]
    (OUT / "PATH2SPACE_REANALYSIS_REPORT.md").write_text("\n".join(lines) + "\n")
    summary_json = {
        "top_smoothed_target": str(top.name),
        "top_smoothed_rho": float(top.overall_spearman_rho),
        "dm1_raw_rho": float(raw.loc["DM1_like_score", "overall_spearman_rho"]),
        "dm1_smoothed_rho": dm1_rho,
        "dm1_smoothing_gain": float(dm1_gain),
        "mapk_smoothed_rho": mapk_rho,
        "best_spatial_cluster_k": int(cluster_meta["best_k"]),
        "best_spatial_cluster_silhouette": float(cluster_meta["best_silhouette"]),
        "dominant_cluster_stage_ari": float(cluster_meta["dominant_cluster_stage_ari"]),
        "n_spots_with_embeddings": int(summary["n_spots"].max()),
        "n_slides": int(slide_summary["sample_id"].nunique()),
    }
    (OUT / "PATH2SPACE_REANALYSIS_SUMMARY.json").write_text(json.dumps(summary_json, indent=2) + "\n")
    return summary_json


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    emb = np.load(PHASE1 / "uni_embeddings_size224.npz")["embeddings"].astype(np.float32)
    meta = pd.read_csv(PHASE1 / "uni_embed_metadata_size224.tsv", sep="\t").reset_index().rename(columns={"index": "embed_idx"})
    scores = pd.read_csv(SCORES, sep="\t")
    scores = add_mapk_score(scores)
    targets = [c for c in BASE_TARGETS if c in scores.columns]
    scores = spatial_smooth_by_slide(scores, targets, k=SMOOTH_K)
    meta["key"] = meta["slide"].astype(str) + "_" + meta["spot_id"].astype(str)
    scores["key"] = scores["sample_id"].astype(str) + "_" + scores["spot_id"].astype(str)
    cols = [
        "key",
        "sample_id",
        "spot_id",
        "stage",
        "array_row",
        "array_col",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ] + targets + [f"{t}_smooth{SMOOTH_K}" for t in targets]
    df = meta.merge(scores[cols], on="key", how="inner")
    df = df.rename(columns={"spot_id_x": "spot_id"}).drop(columns=[c for c in ["spot_id_y"] if c in df.columns])
    df = df.sort_values("embed_idx").reset_index(drop=True)
    emb = emb[df["embed_idx"].to_numpy()]
    df.to_csv(OUT / "path2space_input_spots.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    summary, slide_summary, preds = run_prediction(df, emb)
    summary.to_csv(OUT / "path2space_loso_target_summary.tsv", sep="\t", index=False, na_rep="NA")
    slide_summary.to_csv(OUT / "path2space_loso_slide_summary.tsv", sep="\t", index=False, na_rep="NA")
    preds.to_csv(OUT / "path2space_spot_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    spand = spand_table(df, preds)
    spand.to_csv(OUT / "path2space_spand_by_slide.tsv", sep="\t", index=False, na_rep="NA")

    cluster_spots, cluster_props, cluster_char, cluster_info = cluster_architecture(preds)
    cluster_spots.to_csv(OUT / "path2space_spatial_cluster_spots.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    cluster_props.to_csv(OUT / "path2space_spatial_cluster_proportions.tsv", sep="\t", index=False, na_rep="NA")
    cluster_char.to_csv(OUT / "path2space_spatial_cluster_characterization.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(cluster_info["silhouette"]).to_csv(OUT / "path2space_cluster_k_silhouette.tsv", sep="\t", index=False)
    pd.DataFrame(cluster_info["stage_tests"]).to_csv(OUT / "path2space_cluster_stage_tests.tsv", sep="\t", index=False, na_rep="NA")

    plot_outputs(summary, preds, spand, cluster_props, cluster_char)
    summary_json = write_report(summary, slide_summary, spand, cluster_props, cluster_char, cluster_info)
    print(json.dumps(summary_json, indent=2))


if __name__ == "__main__":
    main()
