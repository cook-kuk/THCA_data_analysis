#!/usr/bin/env python3
"""Feature-family ablations for GSE230424 H&E-to-DM1/RAI signal."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
BASE_SCRIPT = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/scripts/gse230424_pathology_thyroid_axis.py"
BASE_OUT = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp/gse230424_pathology_thyroid_axis_2026_05_09"
OUT = BASE_OUT / "feature_ablation_controls"

spec = importlib.util.spec_from_file_location("gse230424_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)


def load_data() -> pd.DataFrame:
    df = pd.read_csv(BASE_OUT / "gse230424_spot_module_scores.tsv.gz", sep="\t")
    if "disease_group" not in df.columns:
        df.insert(1, "disease_group", df["sample"].map(base.DISEASE_BY_SAMPLE))
    df = base.spatial_smooth(df, ["DM1_low_RAI_score", "MAPK_output_score", "HLA_II_AP_score", "Tumor_ZCCHC12_score"])
    feats = pd.read_csv(BASE_OUT / "gse230424_he_tile_features.tsv.gz", sep="\t")
    return df.merge(feats, on=["sample", "barcode", "array_row", "array_col"], how="inner")


def feature_sets(df: pd.DataFrame) -> dict[str, list[str]]:
    he_cols = [c for c in df.columns if c.startswith("r") and any(token in c for token in ["_mean", "_std", "_p10", "_p90", "_frac", "_var", "_proxy"])]
    radius_48 = [c for c in he_cols if c.startswith("r48_")]
    radius_96 = [c for c in he_cols if c.startswith("r96_")]
    rgb_distribution = [c for c in he_cols if any(tok in c for tok in ["_r_", "_g_", "_b_"])]
    hsv_gray = [c for c in he_cols if any(tok in c for tok in ["_h_", "_s_", "_v_", "_gray_"])]
    density_texture = [c for c in he_cols if any(tok in c for tok in ["dark_frac", "saturated_frac", "lap_var", "gray_std"])]
    stain_proxy = [c for c in he_cols if any(tok in c for tok in ["blue_purple_proxy", "eosin_proxy"])]
    color_means_only = [c for c in he_cols if c.endswith("_mean") and any(tok in c for tok in ["_r_", "_g_", "_b_", "_h_", "_s_", "_v_", "_gray_"])]
    return {
        "all_HE_features": he_cols,
        "radius48_only": radius_48,
        "radius96_only": radius_96,
        "rgb_distribution": rgb_distribution,
        "hsv_gray_summary": hsv_gray,
        "density_texture": density_texture,
        "stain_proxy_only": stain_proxy,
        "color_means_only": color_means_only,
    }


def covariates(df: pd.DataFrame) -> np.ndarray:
    coord = df[["array_row", "array_col"]].to_numpy(dtype=float)
    qc = np.column_stack([
        np.log1p(df["total_counts"].to_numpy(dtype=float)),
        np.log1p(df["n_genes_by_counts"].to_numpy(dtype=float)),
        df["pct_counts_mt"].to_numpy(dtype=float),
    ])
    return np.column_stack([coord, qc])


def evaluate(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = df["sample"].astype(str).to_numpy()
    sets = feature_sets(df)
    rows = []
    pred_cols = df[["sample", "disease_group", "barcode", "array_row", "array_col"]].copy()
    for target in ["DM1_low_RAI_score", "MAPK_output_score", "HLA_II_AP_score", "Tumor_ZCCHC12_score"]:
        y_raw = df[f"{target}_smooth{base.SMOOTH_K}"].to_numpy(dtype=float)
        y_resid = base.residualize_within_sample(y_raw, covariates(df), samples)
        pred_cols[f"obs_{target}"] = y_raw
        pred_cols[f"obs_{target}_coord_qc_resid"] = y_resid
        for target_mode, y in {"raw_smoothed": y_raw, "coord_qc_residual_target": y_resid}.items():
            y_center = base.center_by_group(y, samples)
            for set_name, cols in sets.items():
                if not cols:
                    continue
                x = df[cols].to_numpy(dtype=float)
                pred = base.loso_predict(x, y, samples, use_pca=True)
                if target == "DM1_low_RAI_score":
                    pred_cols[f"pred_{target_mode}_{set_name}"] = pred
                rho, p, n = base.safe_spearman(y, pred)
                pred_center = base.center_by_group(pred, samples)
                crho, cp, cn = base.safe_spearman(y_center, pred_center)
                per_sample = []
                for sample, idx_raw in df.groupby("sample").groups.items():
                    idx = np.asarray(list(idx_raw))
                    srho, _, _ = base.safe_spearman(y[idx], pred[idx])
                    per_sample.append(srho)
                rows.append({
                    "target": target,
                    "target_mode": target_mode,
                    "feature_set": set_name,
                    "n": n,
                    "n_features": len(cols),
                    "pooled_rho": rho,
                    "sample_centered_rho": crho,
                    "median_sample_rho": float(np.nanmedian(per_sample)),
                    "samples_rho_gt_0_1": int(np.nansum(np.asarray(per_sample) > 0.1)),
                })
    return pd.DataFrame(rows), pred_cols


def write_outputs(result: pd.DataFrame, preds: pd.DataFrame) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT / "gse230424_feature_ablation_model_comparison.tsv", sep="\t", index=False, na_rep="NA")
    preds.to_csv(OUT / "gse230424_feature_ablation_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    dm1_raw = result[(result["target"] == "DM1_low_RAI_score") & (result["target_mode"] == "raw_smoothed")].set_index("feature_set")
    dm1_resid = result[(result["target"] == "DM1_low_RAI_score") & (result["target_mode"] == "coord_qc_residual_target")].set_index("feature_set")
    top_raw = dm1_raw["sample_centered_rho"].idxmax()
    top_resid = dm1_resid["sample_centered_rho"].idxmax()
    summary = {
        "n_spots": int(result["n"].max()),
        "n_samples": int(preds["sample"].nunique()),
        "dm1_raw_top_feature_set": top_raw,
        "dm1_raw_top_rho": float(dm1_raw.loc[top_raw, "sample_centered_rho"]),
        "dm1_all_he_raw_rho": float(dm1_raw.loc["all_HE_features", "sample_centered_rho"]),
        "dm1_density_texture_raw_rho": float(dm1_raw.loc["density_texture", "sample_centered_rho"]),
        "dm1_rgb_raw_rho": float(dm1_raw.loc["rgb_distribution", "sample_centered_rho"]),
        "dm1_residual_top_feature_set": top_resid,
        "dm1_residual_top_rho": float(dm1_resid.loc[top_resid, "sample_centered_rho"]),
        "dm1_all_he_residual_rho": float(dm1_resid.loc["all_HE_features", "sample_centered_rho"]),
        "dm1_density_texture_residual_rho": float(dm1_resid.loc["density_texture", "sample_centered_rho"]),
        "dm1_rgb_residual_rho": float(dm1_resid.loc["rgb_distribution", "sample_centered_rho"]),
    }

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.2))
    order = ["all_HE_features", "radius48_only", "radius96_only", "rgb_distribution", "hsv_gray_summary", "density_texture", "stain_proxy_only", "color_means_only"]
    colors = ["#266b73", "#7aa6a1", "#7aa6a1", "#d5a84d", "#d5a84d", "#b75f69", "#c8913f", "#9aa4ad"]
    axes[0].bar(order, dm1_raw.loc[order, "sample_centered_rho"], color=colors)
    axes[0].axhline(0, color="#333", lw=0.8)
    axes[0].set_ylabel("Sample-centered LOSO rho")
    axes[0].set_title("A. DM1/low-RAI raw spatial target")
    axes[0].tick_params(axis="x", rotation=35)
    axes[1].bar(order, dm1_resid.loc[order, "sample_centered_rho"], color=colors)
    axes[1].axhline(0, color="#333", lw=0.8)
    axes[1].set_ylabel("Sample-centered LOSO rho")
    axes[1].set_title("B. DM1/low-RAI coord+QC residual target")
    axes[1].tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse230424_feature_ablation_controls.png", dpi=220)
    fig.savefig(OUT / "fig_gse230424_feature_ablation_controls.pdf")
    plt.close(fig)

    lines = [
        "# GSE230424 feature-family ablation controls",
        "",
        "## Verdict",
        "",
        f"- DM1/low-RAI raw target all-feature rho: **{summary['dm1_all_he_raw_rho']:.3f}**.",
        f"- Best raw feature family: **{summary['dm1_raw_top_feature_set']}** rho **{summary['dm1_raw_top_rho']:.3f}**.",
        f"- DM1/low-RAI coord+QC residual target all-feature rho: **{summary['dm1_all_he_residual_rho']:.3f}**.",
        f"- Best residual feature family: **{summary['dm1_residual_top_feature_set']}** rho **{summary['dm1_residual_top_rho']:.3f}**.",
        f"- Density/texture residual rho: {summary['dm1_density_texture_residual_rho']:.3f}; RGB-distribution residual rho: {summary['dm1_rgb_residual_rho']:.3f}.",
        "",
        "## Interpretation",
        "",
        "The DM1/RAI signal is not a single-feature artifact, but density/texture and stain/color summaries carry much of the raw signal. After coordinate+QC residualization, the retained signal is smaller and feature-family dependent; this supports a cautious Paper 2 image-to-spatial-RNA claim with explicit QC/stain caveats.",
        "",
        "## Full Results",
        "",
        result.round(4).to_markdown(index=False),
    ]
    (OUT / "GSE230424_FEATURE_ABLATION_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "GSE230424_FEATURE_ABLATION_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    df = load_data()
    result, preds = evaluate(df)
    summary = write_outputs(result, preds)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
