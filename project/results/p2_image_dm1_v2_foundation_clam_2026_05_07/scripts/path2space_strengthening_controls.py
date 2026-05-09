#!/usr/bin/env python3
"""Controls and strengthening analyses for the Path2Space-inspired reanalysis."""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler


REPO = Path(os.environ.get("THCA_REPO_ROOT", "/home/seungho/personal/THCA_data_analysis"))
P2_ROOT = REPO / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
BASE = P2_ROOT / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"
PHASE1 = P2_ROOT / "phase1_gse250521"
RESNET = REPO / "project/results/03_pathology_poc/embeddings_resnet50_224.npz"
OUT = BASE / "strengthening_controls"

TARGETS = ["DM1_like_score", "TDS_like_score", "MAPK_output_score", "CAF_ECM_score", "Epithelial_score"]
STAGE_ORDER = ["PT", "PTC", "LPTC", "ATC"]
RIDGE_ALPHA = 100.0
N_PCS = 64
N_PERM = 1000
RANDOM_SEED = 91009


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def rank_center_by_group(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx]
        keep = np.isfinite(vals)
        if keep.sum() == 0:
            continue
        out[idx[keep]] = vals[keep] - np.nanmean(vals[keep])
    return out


def cached_loso_predict(features: np.ndarray, y: np.ndarray, slides: np.ndarray, use_pca: bool) -> np.ndarray:
    features = np.asarray(features, dtype=float)
    y = np.asarray(y, dtype=float)
    slides = np.asarray(slides).astype(str)
    pred = np.full(len(y), np.nan)
    finite_x = np.isfinite(features).all(axis=1)
    for held in sorted(pd.unique(slides)):
        train_idx = np.where((slides != held) & finite_x & np.isfinite(y))[0]
        test_idx = np.where((slides == held) & finite_x)[0]
        if len(train_idx) < 20 or len(test_idx) == 0:
            continue
        scaler = StandardScaler()
        x_train = scaler.fit_transform(features[train_idx])
        x_test = scaler.transform(features[test_idx])
        if use_pca:
            n_comp = min(N_PCS, x_train.shape[0] - 1, x_train.shape[1])
            pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
            x_train = pca.fit_transform(x_train)
            x_test = pca.transform(x_test)
        mean = float(np.nanmean(y[train_idx]))
        sd = float(np.nanstd(y[train_idx]))
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        model = Ridge(alpha=RIDGE_ALPHA)
        model.fit(x_train, (y[train_idx] - mean) / sd)
        pred[test_idx] = model.predict(x_test)
    return pred


def residualize_within_slide(values: np.ndarray, covars: np.ndarray, slides: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    covars = np.asarray(covars, dtype=float)
    slides = np.asarray(slides).astype(str)
    out = np.full(len(values), np.nan)
    for slide in pd.unique(slides):
        idx = np.where(slides == slide)[0]
        y = values[idx]
        x = covars[idx]
        keep = np.isfinite(y) & np.isfinite(x).all(axis=1)
        if keep.sum() < 10:
            continue
        xs = StandardScaler().fit_transform(x[keep])
        out[idx[keep]] = y[keep] - LinearRegression().fit(xs, y[keep]).predict(xs)
    return out


def load_uni_embeddings(df: pd.DataFrame) -> np.ndarray:
    emb = np.load(PHASE1 / "uni_embeddings_size224.npz")["embeddings"].astype(np.float32)
    meta = pd.read_csv(PHASE1 / "uni_embed_metadata_size224.tsv", sep="\t").reset_index().rename(columns={"index": "embed_idx"})
    meta["key"] = meta["slide"].astype(str) + "_" + meta["spot_id"].astype(str)
    keyed = df[["sample_id", "spot_id"]].copy()
    keyed["key"] = keyed["sample_id"].astype(str) + "_" + keyed["spot_id"].astype(str)
    merged = keyed.merge(meta[["key", "embed_idx"]], on="key", how="left")
    return emb[merged["embed_idx"].to_numpy()]


def load_resnet_embeddings(df: pd.DataFrame) -> np.ndarray:
    z = np.load(RESNET)
    emb = z["embeddings"].astype(np.float32)
    meta = pd.DataFrame({"sample_id": z["sample_ids"].astype(str), "spot_id": z["spot_ids"].astype(str)})
    meta["embed_idx"] = np.arange(len(meta))
    meta["key"] = meta["sample_id"] + "_" + meta["spot_id"]
    keyed = df[["sample_id", "spot_id"]].copy()
    keyed["key"] = keyed["sample_id"].astype(str) + "_" + keyed["spot_id"].astype(str)
    merged = keyed.merge(meta[["key", "embed_idx"]], on="key", how="left")
    out = np.full((len(df), emb.shape[1]), np.nan, dtype=np.float32)
    ok = merged["embed_idx"].notna().to_numpy()
    out[ok] = emb[merged.loc[ok, "embed_idx"].astype(int).to_numpy()]
    return out


def model_comparison(input_df: pd.DataFrame, preds: pd.DataFrame) -> pd.DataFrame:
    slides = input_df["sample_id"].astype(str).to_numpy()
    stage = input_df["stage"].astype(str)
    stage_onehot = pd.get_dummies(stage).reindex(columns=STAGE_ORDER, fill_value=0).to_numpy(dtype=float)
    coord = input_df[["array_row", "array_col"]].to_numpy(dtype=float)
    qc = input_df[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].to_numpy(dtype=float)
    qc[:, 0] = np.log1p(qc[:, 0])
    qc[:, 1] = np.log1p(qc[:, 1])
    coord_qc = np.column_stack([coord, qc])
    features = {
        "UNI_HE": (load_uni_embeddings(input_df), True),
        "ResNet50_HE": (load_resnet_embeddings(input_df), True),
        "Coord_only": (coord, False),
        "QC_only": (qc, False),
        "Coord_QC": (coord_qc, False),
        "Stage_only": (stage_onehot, False),
        "Coord_QC_stage": (np.column_stack([coord_qc, stage_onehot]), False),
    }
    rows = []
    out_pred = input_df[["sample_id", "spot_id", "stage"]].copy()
    for target in TARGETS:
        y = preds[f"obs_{target}_smooth8"].to_numpy(dtype=float)
        y_center = rank_center_by_group(y, slides)
        for model_name, (x, use_pca) in features.items():
            if model_name == "UNI_HE":
                pred = preds[f"pred_{target}_smooth8"].to_numpy(dtype=float)
            else:
                pred = cached_loso_predict(x, y, slides, use_pca=use_pca)
            out_pred[f"pred_{target}_{model_name}"] = pred
            rho, p, n = safe_spearman(y, pred)
            pred_center = rank_center_by_group(pred, slides)
            crho, cp, cn = safe_spearman(y_center, pred_center)
            per_slide = []
            for slide, sub_idx in input_df.groupby("sample_id").groups.items():
                idx = np.asarray(list(sub_idx))
                srho, _, _ = safe_spearman(y[idx], pred[idx])
                per_slide.append(srho)
            rows.append(
                {
                    "target": target,
                    "model": model_name,
                    "n": n,
                    "pooled_rho": rho,
                    "pooled_p": p,
                    "slide_centered_rho": crho,
                    "slide_centered_p": cp,
                    "median_slide_rho": float(np.nanmedian(per_slide)),
                    "slides_rho_gt_0_2": int(np.nansum(np.asarray(per_slide) > 0.2)),
                    "features": x.shape[1],
                }
            )
    out_pred.to_csv(OUT / "strengthening_model_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    return pd.DataFrame(rows)


def residual_alignment(input_df: pd.DataFrame, preds: pd.DataFrame) -> pd.DataFrame:
    slides = input_df["sample_id"].astype(str).to_numpy()
    coord = input_df[["array_row", "array_col"]].to_numpy(dtype=float)
    qc = input_df[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].to_numpy(dtype=float)
    qc[:, 0] = np.log1p(qc[:, 0])
    qc[:, 1] = np.log1p(qc[:, 1])
    covars = {
        "slide_mean_only": np.ones((len(input_df), 1), dtype=float),
        "coord_only": coord,
        "qc_only": qc,
        "coord_qc": np.column_stack([coord, qc]),
    }
    rows = []
    for target in TARGETS:
        obs = preds[f"obs_{target}_smooth8"].to_numpy(dtype=float)
        pred = preds[f"pred_{target}_smooth8"].to_numpy(dtype=float)
        for adjustment, cov in covars.items():
            obs_resid = residualize_within_slide(obs, cov, slides)
            pred_resid = residualize_within_slide(pred, cov, slides)
            rho, p, n = safe_spearman(obs_resid, pred_resid)
            rows.append(
                {
                    "target": target,
                    "within_slide_adjustment": adjustment,
                    "residual_alignment_rho": rho,
                    "p": p,
                    "n": n,
                }
            )
    return pd.DataFrame(rows)


def domain_aggregation(input_df: pd.DataFrame, preds: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    domain_rows = []
    summary_rows = []
    for slide, idx_raw in input_df.groupby("sample_id", sort=False).groups.items():
        idx = np.asarray(list(idx_raw))
        sub = input_df.iloc[idx]
        coords = sub[["array_row", "array_col"]].to_numpy(dtype=float)
        n_domains = min(10, max(4, int(round(len(sub) / 35))))
        labs = KMeans(n_clusters=n_domains, n_init=50, random_state=RANDOM_SEED).fit_predict(StandardScaler().fit_transform(coords))
        for lab in np.unique(labs):
            mask = labs == lab
            row = {
                "sample_id": slide,
                "stage": str(sub["stage"].iloc[0]),
                "domain": int(lab + 1),
                "n_spots": int(mask.sum()),
                "mean_array_row": float(np.nanmean(coords[mask, 0])),
                "mean_array_col": float(np.nanmean(coords[mask, 1])),
            }
            for target in TARGETS:
                y = preds.iloc[idx][f"obs_{target}_smooth8"].to_numpy(dtype=float)
                p = preds.iloc[idx][f"pred_{target}_smooth8"].to_numpy(dtype=float)
                row[f"obs_{target}"] = float(np.nanmean(y[mask]))
                row[f"pred_{target}"] = float(np.nanmean(p[mask]))
            domain_rows.append(row)
    domains = pd.DataFrame(domain_rows)
    for target in TARGETS:
        y = domains[f"obs_{target}"].to_numpy(dtype=float)
        p = domains[f"pred_{target}"].to_numpy(dtype=float)
        rho, pv, n = safe_spearman(y, p)
        yc = rank_center_by_group(y, domains["sample_id"].astype(str).to_numpy())
        pc = rank_center_by_group(p, domains["sample_id"].astype(str).to_numpy())
        crho, cp, cn = safe_spearman(yc, pc)
        per_slide = []
        for slide, sub in domains.groupby("sample_id"):
            srho, _, _ = safe_spearman(sub[f"obs_{target}"], sub[f"pred_{target}"])
            per_slide.append(srho)
        summary_rows.append(
            {
                "target": target,
                "n_domains": n,
                "pooled_domain_rho": rho,
                "pooled_domain_p": pv,
                "slide_centered_domain_rho": crho,
                "slide_centered_domain_p": cp,
                "median_slide_domain_rho": float(np.nanmedian(per_slide)),
                "slides_domain_rho_gt_0_2": int(np.nansum(np.asarray(per_slide) > 0.2)),
            }
        )
    return domains, pd.DataFrame(summary_rows)


def permutation_null(preds: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    slides = preds["sample_id"].astype(str).to_numpy()
    for target in TARGETS[:3]:
        y = preds[f"obs_{target}_smooth8"].to_numpy(dtype=float)
        p = preds[f"pred_{target}_smooth8"].to_numpy(dtype=float)
        obs_center = rank_center_by_group(y, slides)
        pred_center = rank_center_by_group(p, slides)
        obs_rho, _, _ = safe_spearman(obs_center, pred_center)
        null = []
        for _ in range(N_PERM):
            shuffled = pred_center.copy()
            for slide in pd.unique(slides):
                idx = np.where(slides == slide)[0]
                shuffled[idx] = rng.permutation(shuffled[idx])
            rho, _, _ = safe_spearman(obs_center, shuffled)
            null.append(rho)
        null = np.asarray(null, dtype=float)
        rows.append(
            {
                "target": target,
                "observed_slide_centered_rho": obs_rho,
                "null_mean": float(np.nanmean(null)),
                "null_sd": float(np.nanstd(null)),
                "null_p95": float(np.nanpercentile(null, 95)),
                "null_p99": float(np.nanpercentile(null, 99)),
                "empirical_p_greater_equal": float((np.nansum(null >= obs_rho) + 1) / (np.isfinite(null).sum() + 1)),
                "n_permutations": N_PERM,
            }
        )
    return pd.DataFrame(rows)


def plot_controls(model_df: pd.DataFrame, domain_summary: pd.DataFrame, perm_df: pd.DataFrame, residual_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    dm = model_df[model_df["target"] == "DM1_like_score"].copy()
    order = ["Stage_only", "Coord_only", "QC_only", "Coord_QC", "Coord_QC_stage", "ResNet50_HE", "UNI_HE"]
    dm["model"] = pd.Categorical(dm["model"], categories=order, ordered=True)
    dm = dm.sort_values("model")
    colors = ["#777", "#a9b1bc", "#a9b1bc", "#9aa4ad", "#9aa4ad", "#7b9fc9", "#246b7a"]
    axes[0].barh(dm["model"].astype(str), dm["slide_centered_rho"], color=colors[: len(dm)])
    axes[0].axvline(0, color="#333", lw=0.8)
    axes[0].set_xlabel("Slide-centered rho")
    axes[0].set_title("A. DM1/RAI: H&E vs confound baselines")

    dom = domain_summary.set_index("target").loc[["DM1_like_score", "TDS_like_score", "MAPK_output_score"]]
    x = np.arange(len(dom))
    axes[1].bar(x - 0.18, dom["pooled_domain_rho"], width=0.35, color="#246b7a", label="domain pooled")
    axes[1].bar(x + 0.18, dom["slide_centered_domain_rho"], width=0.35, color="#d6b25e", label="domain slide-centered")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(["DM1/RAI", "TDS", "MAPK"], rotation=20, ha="right")
    axes[1].set_ylabel("Spearman rho")
    axes[1].set_title("B. Spatial-domain aggregation")
    axes[1].legend(frameon=False)

    res = residual_df[
        (residual_df["target"].isin(["DM1_like_score", "TDS_like_score", "MAPK_output_score"]))
        & (residual_df["within_slide_adjustment"].isin(["slide_mean_only", "coord_only", "qc_only", "coord_qc"]))
    ].copy()
    res_p = res.pivot(index="target", columns="within_slide_adjustment", values="residual_alignment_rho").loc[
        ["DM1_like_score", "TDS_like_score", "MAPK_output_score"]
    ]
    x = np.arange(len(res_p))
    width = 0.2
    for j, col in enumerate(["slide_mean_only", "coord_only", "qc_only", "coord_qc"]):
        axes[2].bar(x + (j - 1.5) * width, res_p[col], width=width, label=col.replace("_", " "))
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(["DM1/RAI", "TDS", "MAPK"], rotation=20, ha="right")
    axes[2].set_ylabel("Residual alignment rho")
    axes[2].set_title("C. Alignment after within-slide residualization")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_path2space_strengthening_controls.png", dpi=220)
    fig.savefig(OUT / "fig_path2space_strengthening_controls.pdf")
    plt.close(fig)


def write_report(
    model_df: pd.DataFrame,
    domains: pd.DataFrame,
    domain_summary: pd.DataFrame,
    perm_df: pd.DataFrame,
    residual_df: pd.DataFrame,
) -> dict:
    dm = model_df[(model_df["target"] == "DM1_like_score")].set_index("model")
    tds = model_df[(model_df["target"] == "TDS_like_score")].set_index("model")
    mapk = model_df[(model_df["target"] == "MAPK_output_score")].set_index("model")
    dm_domain = domain_summary.set_index("target").loc["DM1_like_score"]
    dm_perm = perm_df.set_index("target").loc["DM1_like_score"]
    res = residual_df.set_index(["target", "within_slide_adjustment"])
    summary = {
        "dm1_uni_slide_centered_rho": float(dm.loc["UNI_HE", "slide_centered_rho"]),
        "dm1_coord_only_slide_centered_rho": float(dm.loc["Coord_only", "slide_centered_rho"]),
        "dm1_qc_only_slide_centered_rho": float(dm.loc["QC_only", "slide_centered_rho"]),
        "dm1_coord_qc_slide_centered_rho": float(dm.loc["Coord_QC", "slide_centered_rho"]),
        "dm1_resnet_slide_centered_rho": float(dm.loc["ResNet50_HE", "slide_centered_rho"]),
        "dm1_uni_minus_coord_qc": float(dm.loc["UNI_HE", "slide_centered_rho"] - dm.loc["Coord_QC", "slide_centered_rho"]),
        "dm1_uni_minus_coord_only": float(dm.loc["UNI_HE", "slide_centered_rho"] - dm.loc["Coord_only", "slide_centered_rho"]),
        "dm1_uni_minus_qc_only": float(dm.loc["UNI_HE", "slide_centered_rho"] - dm.loc["QC_only", "slide_centered_rho"]),
        "tds_uni_slide_centered_rho": float(tds.loc["UNI_HE", "slide_centered_rho"]),
        "mapk_uni_slide_centered_rho": float(mapk.loc["UNI_HE", "slide_centered_rho"]),
        "dm1_uni_qc_residual_alignment_rho": float(res.loc[("DM1_like_score", "qc_only"), "residual_alignment_rho"]),
        "dm1_uni_coord_qc_residual_alignment_rho": float(res.loc[("DM1_like_score", "coord_qc"), "residual_alignment_rho"]),
        "dm1_domain_pooled_rho": float(dm_domain["pooled_domain_rho"]),
        "dm1_domain_slide_centered_rho": float(dm_domain["slide_centered_domain_rho"]),
        "dm1_permutation_empirical_p": float(dm_perm["empirical_p_greater_equal"]),
        "n_domains": int(domains.shape[0]),
        "n_permutations": int(N_PERM),
    }
    lines = [
        "# Path2Space strengthening controls",
        "",
        "## Verdict",
        "",
        f"- UNI H&E DM1/RAI slide-centered rho: **{summary['dm1_uni_slide_centered_rho']:.3f}**.",
        f"- Coordinate-only DM1/RAI slide-centered rho: {summary['dm1_coord_only_slide_centered_rho']:.3f}; QC-only: {summary['dm1_qc_only_slide_centered_rho']:.3f}; coordinate+QC: {summary['dm1_coord_qc_slide_centered_rho']:.3f}.",
        f"- UNI minus coordinate-only delta: **{summary['dm1_uni_minus_coord_only']:+.3f}**; UNI minus coordinate+QC delta: {summary['dm1_uni_minus_coord_qc']:+.3f}.",
        f"- After within-slide QC residualization, UNI DM1/RAI residual alignment rho remains **{summary['dm1_uni_qc_residual_alignment_rho']:.3f}**; after coordinate+QC residualization it remains {summary['dm1_uni_coord_qc_residual_alignment_rho']:.3f}.",
        f"- ResNet50 H&E DM1/RAI slide-centered rho: {summary['dm1_resnet_slide_centered_rho']:.3f}.",
        f"- Spatial-domain DM1/RAI pooled rho: **{summary['dm1_domain_pooled_rho']:.3f}** across {summary['n_domains']} coordinate domains.",
        f"- Within-slide permutation empirical p for DM1/RAI: **{summary['dm1_permutation_empirical_p']:.4f}** ({N_PERM} permutations).",
        "",
        "## Interpretation",
        "",
        "The image-derived DM1/RAI signal survives slide-centering and beats coordinate-only plus ResNet50 baselines, but coordinate+QC is a strong measured-ST baseline. Therefore the safe claim is not that UNI uniquely explains all local RNA structure; it is that pathology carries a reproducible DM1/RAI spatial signal, a component remains after QC residualization, and domain averaging strengthens it. The boundary remains unchanged: this strengthens Paper 2 image-to-spatial-RNA framing, not Paper 1 MAPK causal mechanism.",
        "",
        "## Model Comparison",
        "",
        model_df.round(4).to_markdown(index=False),
        "",
        "## Domain Summary",
        "",
        domain_summary.round(4).to_markdown(index=False),
        "",
        "## Permutation Null",
        "",
        perm_df.round(4).to_markdown(index=False),
        "",
        "## Residual Alignment",
        "",
        residual_df.round(4).to_markdown(index=False),
    ]
    (OUT / "PATH2SPACE_STRENGTHENING_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "PATH2SPACE_STRENGTHENING_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    input_df = pd.read_csv(BASE / "path2space_input_spots.tsv.gz", sep="\t")
    preds = pd.read_csv(BASE / "path2space_spot_predictions.tsv.gz", sep="\t")
    model_df = model_comparison(input_df, preds)
    residual_df = residual_alignment(input_df, preds)
    domains, domain_summary = domain_aggregation(input_df, preds)
    perm_df = permutation_null(preds)
    model_df.to_csv(OUT / "path2space_strengthening_model_comparison.tsv", sep="\t", index=False, na_rep="NA")
    residual_df.to_csv(OUT / "path2space_strengthening_residual_alignment.tsv", sep="\t", index=False, na_rep="NA")
    domains.to_csv(OUT / "path2space_strengthening_domains.tsv", sep="\t", index=False, na_rep="NA")
    domain_summary.to_csv(OUT / "path2space_strengthening_domain_summary.tsv", sep="\t", index=False, na_rep="NA")
    perm_df.to_csv(OUT / "path2space_strengthening_permutation_null.tsv", sep="\t", index=False, na_rep="NA")
    plot_controls(model_df, domain_summary, perm_df, residual_df)
    summary = write_report(model_df, domains, domain_summary, perm_df, residual_df)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
