#!/usr/bin/env python3
"""GSE248205 H&E-to-AITD spatial-axis scout.

Fast Path2Space-style screen using handcrafted H&E tile features. This is not
a foundation-model result; it asks whether local thyroid histology carries the
AITD antigen-presentation/TLS axes already measured in Visium spots.
"""
from __future__ import annotations

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
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC = ROOT / "project/data/external/GSE248205/Processed files"
MODULES = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09/gse248205_aitd_spatial_validation/tables/T02_gse248205_spot_module_scores.tsv.gz"
OUT = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/gse248205_pathology_aitd_axis_2026_05_09"

TARGETS = ["AP_TLS_composite", "HLA_II_AP", "B_TLS", "CD74_MIF_axis", "Thyrocyte"]
GROUP_ORDER = ["Control", "GD", "HT"]
RANDOM_SEED = 248205
SMOOTH_K = 8
N_PCS = 48
RIDGE_ALPHA = 50.0
PATCH_RADII = [32, 72]


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def spatial_smooth(df: pd.DataFrame, cols: list[str], k: int = SMOOTH_K) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        out[f"{col}_smooth{k}"] = np.nan
    for sample, sub in out.groupby("sample", sort=False):
        idx = sub.index.to_numpy()
        coords = sub[["array_row", "array_col"]].to_numpy(dtype=float)
        tree = cKDTree(coords)
        kk = min(k + 1, len(idx))
        _, neigh = tree.query(coords, k=kk)
        if neigh.ndim == 1:
            neigh = neigh[:, None]
        for col in cols:
            vals = sub[col].to_numpy(dtype=float)
            out.loc[idx, f"{col}_smooth{k}"] = np.nanmean(vals[neigh], axis=1)
    return out


def rank_center_by_group(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    out = np.full(len(values), np.nan)
    for group in pd.unique(groups):
        idx = np.where(groups == group)[0]
        vals = values[idx].astype(float)
        keep = np.isfinite(vals)
        if keep.sum() == 0:
            continue
        out[idx[keep]] = vals[keep] - np.nanmean(vals[keep])
    return out


def patch_features(patch_rgb: np.ndarray, prefix: str) -> dict[str, float]:
    if patch_rgb.size == 0:
        return {}
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
    # Simple hematoxylin/eosin proxies: blue-purple and pink balance.
    r = rgb[:, :, 0]
    g = rgb[:, :, 1]
    b = rgb[:, :, 2]
    features[f"{prefix}_blue_purple_proxy"] = float(np.nanmean((b + r) / 2.0 - g))
    features[f"{prefix}_eosin_proxy"] = float(np.nanmean(r - b))
    return features


def extract_he_features(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, sub in df.groupby("sample", sort=False):
        img = np.asarray(Image.open(SRC / sample / f"{sample}_tissue_hires_image.png").convert("RGB"))
        h, w = img.shape[:2]
        print(f"[GSE248205] H&E features {sample}: {len(sub)} spots")
        for row in sub.itertuples(index=False):
            x = int(round(row.pxl_col_hires))
            y = int(round(row.pxl_row_hires))
            feats = {"sample": row.sample, "barcode": row.barcode, "group": row.group, "array_row": row.array_row, "array_col": row.array_col}
            for radius in PATCH_RADII:
                x0, x1 = max(0, x - radius), min(w, x + radius)
                y0, y1 = max(0, y - radius), min(h, y + radius)
                feats.update(patch_features(img[y0:y1, x0:x1, :], f"r{radius}"))
            rows.append(feats)
    out = pd.DataFrame(rows)
    feature_cols = [c for c in out.columns if c.startswith("r")]
    out[feature_cols] = out[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(out[feature_cols].median(numeric_only=True))
    return out


def load_data() -> pd.DataFrame:
    long = pd.read_csv(MODULES, sep="\t")
    wide = long.pivot_table(
        index=["sample", "group", "barcode", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"],
        columns="module",
        values="spot_log1p_cp10k",
        aggfunc="mean",
    ).reset_index()
    # The coordinates in the prior table are full-resolution; scale them to the
    # shipped high-resolution images.
    scales = {}
    for sample in sorted(wide["sample"].unique()):
        import json as _json

        sf = _json.load(open(SRC / sample / f"{sample}_scalefactors_json.json"))
        scales[sample] = float(sf["tissue_hires_scalef"])
    wide["pxl_col_hires"] = [wide.loc[i, "pxl_col_in_fullres"] * scales[wide.loc[i, "sample"]] for i in wide.index]
    wide["pxl_row_hires"] = [wide.loc[i, "pxl_row_in_fullres"] * scales[wide.loc[i, "sample"]] for i in wide.index]
    wide = spatial_smooth(wide, TARGETS, k=SMOOTH_K)
    return wide


def loso_predict(features: np.ndarray, y: np.ndarray, samples: np.ndarray, use_pca: bool) -> np.ndarray:
    pred = np.full(len(y), np.nan)
    finite_x = np.isfinite(features).all(axis=1)
    for held in sorted(pd.unique(samples)):
        train = (samples != held) & finite_x & np.isfinite(y)
        test = (samples == held) & finite_x
        if train.sum() < 20 or test.sum() == 0:
            continue
        scaler = StandardScaler()
        x_train = scaler.fit_transform(features[train])
        x_test = scaler.transform(features[test])
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
    merged = df.merge(feats, on=["sample", "barcode", "group", "array_row", "array_col"], how="inner")
    samples = merged["sample"].astype(str).to_numpy()
    groups = merged["group"].astype(str)
    group_onehot = pd.get_dummies(groups).reindex(columns=GROUP_ORDER, fill_value=0).to_numpy(dtype=float)
    coord = merged[["array_row", "array_col"]].to_numpy(dtype=float)
    he_cols = [c for c in feats.columns if c.startswith("r")]
    he = merged[he_cols].to_numpy(dtype=float)
    feature_sets = {
        "HE_tile_features": (he, True),
        "Coord_only": (coord, False),
        "Group_only": (group_onehot, False),
        "Coord_group": (np.column_stack([coord, group_onehot]), False),
    }
    pred_df = merged[["sample", "group", "barcode", "array_row", "array_col"]].copy()
    rows = []
    for target in TARGETS:
        y = merged[f"{target}_smooth{SMOOTH_K}"].to_numpy(dtype=float)
        pred_df[f"obs_{target}_smooth{SMOOTH_K}"] = y
        y_center = rank_center_by_group(y, samples)
        for model_name, (x, use_pca) in feature_sets.items():
            pred = loso_predict(x, y, samples, use_pca=use_pca)
            pred_df[f"pred_{target}_{model_name}"] = pred
            rho, p, n = safe_spearman(y, pred)
            pred_center = rank_center_by_group(pred, samples)
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


def domain_aggregation(pred_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, sub in pred_df.groupby("sample", sort=False):
        coords = sub[["array_row", "array_col"]].to_numpy(dtype=float)
        n_domains = min(12, max(5, int(round(len(sub) / 220))))
        labs = KMeans(n_clusters=n_domains, n_init=50, random_state=RANDOM_SEED).fit_predict(StandardScaler().fit_transform(coords))
        for lab in np.unique(labs):
            mask = labs == lab
            row = {"sample": sample, "group": sub["group"].iloc[0], "domain": int(lab + 1), "n_spots": int(mask.sum())}
            for target in TARGETS:
                row[f"obs_{target}"] = float(np.nanmean(sub[f"obs_{target}_smooth{SMOOTH_K}"].to_numpy(dtype=float)[mask]))
                row[f"pred_{target}"] = float(np.nanmean(sub[f"pred_{target}_HE_tile_features"].to_numpy(dtype=float)[mask]))
            rows.append(row)
    domains = pd.DataFrame(rows)
    summary_rows = []
    for target in TARGETS:
        rho, p, n = safe_spearman(domains[f"obs_{target}"], domains[f"pred_{target}"])
        obs_c = rank_center_by_group(domains[f"obs_{target}"].to_numpy(dtype=float), domains["sample"].astype(str).to_numpy())
        pred_c = rank_center_by_group(domains[f"pred_{target}"].to_numpy(dtype=float), domains["sample"].astype(str).to_numpy())
        crho, cp, cn = safe_spearman(obs_c, pred_c)
        summary_rows.append({"target": target, "n_domains": n, "pooled_domain_rho": rho, "sample_centered_domain_rho": crho, "p": p})
    domains.to_csv(OUT / "gse248205_pathology_domains.tsv", sep="\t", index=False, na_rep="NA")
    return pd.DataFrame(summary_rows)


def plot_summary(model_df: pd.DataFrame, domain_df: pd.DataFrame, pred_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    he = model_df[model_df["model"] == "HE_tile_features"].set_index("target").loc[TARGETS]
    x = np.arange(len(TARGETS))
    axes[0].bar(x - 0.18, he["pooled_rho"], width=0.35, color="#246b7a", label="pooled")
    axes[0].bar(x + 0.18, he["sample_centered_rho"], width=0.35, color="#d6b25e", label="sample-centered")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(TARGETS, rotation=30, ha="right")
    axes[0].set_ylabel("LOSO Spearman rho")
    axes[0].set_title("A. H&E tile features to AITD spatial axes")
    axes[0].legend(frameon=False)

    dm = model_df[model_df["target"] == "AP_TLS_composite"].copy()
    order = ["Group_only", "Coord_only", "Coord_group", "HE_tile_features"]
    dm["model"] = pd.Categorical(dm["model"], categories=order, ordered=True)
    dm = dm.sort_values("model")
    axes[1].barh(dm["model"].astype(str), dm["sample_centered_rho"], color=["#777", "#9aa4ad", "#9aa4ad", "#246b7a"])
    axes[1].set_xlabel("Sample-centered rho")
    axes[1].set_title("B. AP/TLS axis against baselines")

    dom = domain_df.set_index("target").loc[TARGETS]
    axes[2].bar(x - 0.18, dom["pooled_domain_rho"], width=0.35, color="#246b7a", label="domain pooled")
    axes[2].bar(x + 0.18, dom["sample_centered_domain_rho"], width=0.35, color="#d6b25e", label="domain centered")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(TARGETS, rotation=30, ha="right")
    axes[2].set_ylabel("Domain Spearman rho")
    axes[2].set_title("C. Domain aggregation")
    axes[2].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse248205_pathology_aitd_axis.png", dpi=220)
    fig.savefig(OUT / "fig_gse248205_pathology_aitd_axis.pdf")
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(6.4, 5.8))
    sub = pred_df[["group", "obs_AP_TLS_composite_smooth8", "pred_AP_TLS_composite_HE_tile_features"]].dropna()
    colors = {"Control": "#718096", "GD": "#2b6cb0", "HT": "#b7791f"}
    for group, part in sub.groupby("group"):
        ax.scatter(part["obs_AP_TLS_composite_smooth8"], part["pred_AP_TLS_composite_HE_tile_features"], s=13, alpha=0.5, color=colors.get(group, "#333"), label=group)
    rho, _, _ = safe_spearman(sub["obs_AP_TLS_composite_smooth8"], sub["pred_AP_TLS_composite_HE_tile_features"])
    ax.set_xlabel("Observed smoothed AP/TLS composite")
    ax.set_ylabel("Predicted from H&E tile features")
    ax.set_title(f"AP/TLS H&E prediction (rho={rho:.3f})")
    ax.legend(frameon=False)
    fig2.tight_layout()
    fig2.savefig(OUT / "fig_gse248205_ap_tls_scatter.png", dpi=220)
    fig2.savefig(OUT / "fig_gse248205_ap_tls_scatter.pdf")
    plt.close(fig2)


def write_report(model_df: pd.DataFrame, domain_df: pd.DataFrame) -> dict:
    he = model_df[model_df["model"] == "HE_tile_features"].set_index("target")
    ap = he.loc["AP_TLS_composite"]
    hla = he.loc["HLA_II_AP"]
    b = he.loc["B_TLS"]
    thy = he.loc["Thyrocyte"]
    ap_base = model_df[(model_df["target"] == "AP_TLS_composite")].set_index("model")
    ap_dom = domain_df.set_index("target").loc["AP_TLS_composite"]
    summary = {
        "n_spots": int(ap["n"]),
        "n_samples": 8,
        "ap_tls_he_pooled_rho": float(ap["pooled_rho"]),
        "ap_tls_he_sample_centered_rho": float(ap["sample_centered_rho"]),
        "hla2_he_sample_centered_rho": float(hla["sample_centered_rho"]),
        "b_tls_he_sample_centered_rho": float(b["sample_centered_rho"]),
        "thyrocyte_he_sample_centered_rho": float(thy["sample_centered_rho"]),
        "ap_tls_coord_only_centered_rho": float(ap_base.loc["Coord_only", "sample_centered_rho"]),
        "ap_tls_group_only_centered_rho": float(ap_base.loc["Group_only", "sample_centered_rho"]),
        "ap_tls_domain_sample_centered_rho": float(ap_dom["sample_centered_domain_rho"]),
        "verdict": "NO_GO_LOCAL_HE_CAVEAT_ONLY",
    }
    lines = [
        "# GSE248205 H&E-to-AITD spatial-axis scout",
        "",
        "## Verdict",
        "",
        f"- Spots: {summary['n_spots']:,} across 8 thyroid Visium samples (Control n=2, GD n=3, HT n=3).",
        f"- H&E tile features -> AP/TLS composite: pooled rho **{summary['ap_tls_he_pooled_rho']:.3f}**, sample-centered rho **{summary['ap_tls_he_sample_centered_rho']:.3f}**.",
        f"- AP/TLS coordinate-only centered rho {summary['ap_tls_coord_only_centered_rho']:.3f}; group-only centered rho {summary['ap_tls_group_only_centered_rho']:.3f}.",
        f"- Domain-level AP/TLS sample-centered rho **{summary['ap_tls_domain_sample_centered_rho']:.3f}**.",
        f"- HLA-II/AP centered rho {summary['hla2_he_sample_centered_rho']:.3f}; B/TLS {summary['b_tls_he_sample_centered_rho']:.3f}; Thyrocyte {summary['thyrocyte_he_sample_centered_rho']:.3f}.",
        "",
        "## Interpretation",
        "",
        "This is a fast handcrafted-feature screen, not a UNI/foundation-model result. The pooled association is driven by sample/group structure and does not survive the sample-centered test. Treat GSE248205 as a caveat-only negative control for local H&E-to-AITD spatial prediction, not as positive Paper 2 evidence and not as Paper 1 mechanism support.",
        "",
        "## Model comparison",
        "",
        model_df.round(4).to_markdown(index=False),
        "",
        "## Domain aggregation",
        "",
        domain_df.round(4).to_markdown(index=False),
    ]
    (OUT / "GSE248205_PATHOLOGY_AITD_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "GSE248205_PATHOLOGY_AITD_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_data()
    feats_path = OUT / "gse248205_he_tile_features.tsv.gz"
    if feats_path.exists():
        feats = pd.read_csv(feats_path, sep="\t")
    else:
        feats = extract_he_features(df)
        feats.to_csv(feats_path, sep="\t", index=False, compression="gzip", na_rep="NA")
    model_df, pred_df = run_models(df, feats)
    domain_df = domain_aggregation(pred_df)
    model_df.to_csv(OUT / "gse248205_pathology_model_summary.tsv", sep="\t", index=False, na_rep="NA")
    pred_df.to_csv(OUT / "gse248205_pathology_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    domain_df.to_csv(OUT / "gse248205_pathology_domain_summary.tsv", sep="\t", index=False, na_rep="NA")
    plot_summary(model_df, domain_df, pred_df)
    summary = write_report(model_df, domain_df)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
