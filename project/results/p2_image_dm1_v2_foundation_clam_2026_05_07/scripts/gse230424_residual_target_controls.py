#!/usr/bin/env python3
"""Residual-target controls for GSE230424 H&E-to-DM1/RAI signal.

This asks a stricter question than residual alignment: after removing
coordinate/QC structure from the observed spatial RNA target within each slide,
can H&E tile features directly predict the remaining residual target?
"""
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
OUT = BASE_OUT / "residual_target_controls"


spec = importlib.util.spec_from_file_location("gse230424_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)


STRICT_TARGETS = [
    "DM1_low_RAI_score",
    "RAI8_lineage_score",
    "MAPK_output_score",
    "HLA_II_AP_score",
    "CD36_SPP1_macrophage_score",
    "Tumor_ZCCHC12_score",
]


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(BASE_OUT / "gse230424_spot_module_scores.tsv.gz", sep="\t")
    df = base.spatial_smooth(df, base.TARGETS)
    feats = pd.read_csv(BASE_OUT / "gse230424_he_tile_features.tsv.gz", sep="\t")
    return df, feats


def build_covariates(df: pd.DataFrame) -> dict[str, np.ndarray]:
    coord = df[["array_row", "array_col"]].to_numpy(dtype=float)
    qc = np.column_stack(
        [
            np.log1p(df["total_counts"].to_numpy(dtype=float)),
            np.log1p(df["n_genes_by_counts"].to_numpy(dtype=float)),
            df["pct_counts_mt"].to_numpy(dtype=float),
        ]
    )
    return {
        "sample_mean_only": np.ones((len(df), 1), dtype=float),
        "coord_only": coord,
        "qc_only": qc,
        "coord_qc": np.column_stack([coord, qc]),
    }


def model_residual_targets(df: pd.DataFrame, feats: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = df.merge(feats, on=["sample", "barcode", "array_row", "array_col"], how="inner")
    samples = merged["sample"].astype(str).to_numpy()
    covars = build_covariates(merged)
    he_cols = [c for c in feats.columns if c.startswith("r")]
    feature_sets = {
        "HE_tile_features": (merged[he_cols].to_numpy(dtype=float), True),
        "Coord_only": (covars["coord_only"], False),
        "QC_only": (covars["qc_only"], False),
        "Coord_QC": (covars["coord_qc"], False),
    }
    rows = []
    pred_df = merged[["sample", "geo_accession", "barcode", "array_row", "array_col"]].copy()
    for target in STRICT_TARGETS:
        y_raw = merged[f"{target}_smooth{base.SMOOTH_K}"].to_numpy(dtype=float)
        for adjustment, cov in covars.items():
            y_resid = base.residualize_within_sample(y_raw, cov, samples)
            pred_df[f"obs_resid_{target}_{adjustment}"] = y_resid
            y_resid_center = base.center_by_group(y_resid, samples)
            for model_name, (x, use_pca) in feature_sets.items():
                pred = base.loso_predict(x, y_resid, samples, use_pca=use_pca)
                if adjustment == "coord_qc":
                    pred_df[f"pred_resid_{target}_{model_name}"] = pred
                rho, p, n = base.safe_spearman(y_resid, pred)
                pred_center = base.center_by_group(pred, samples)
                crho, cp, cn = base.safe_spearman(y_resid_center, pred_center)
                per_sample = []
                for sample, idx_raw in merged.groupby("sample").groups.items():
                    idx = np.asarray(list(idx_raw))
                    srho, _, _ = base.safe_spearman(y_resid[idx], pred[idx])
                    per_sample.append(srho)
                rows.append(
                    {
                        "target": target,
                        "residual_target_adjustment": adjustment,
                        "model": model_name,
                        "n": n,
                        "pooled_rho": rho,
                        "pooled_p": p,
                        "sample_centered_rho": crho,
                        "sample_centered_p": cp,
                        "median_sample_rho": float(np.nanmedian(per_sample)),
                        "samples_rho_gt_0_1": int(np.nansum(np.asarray(per_sample) > 0.1)),
                    }
                )
    return pd.DataFrame(rows), pred_df


def permutation_for_coord_qc(pred_df: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(base.RANDOM_SEED + 17)
    samples = pred_df["sample"].astype(str).to_numpy()
    rows = []
    for target in STRICT_TARGETS:
        obs = pred_df[f"obs_resid_{target}_coord_qc"].to_numpy(dtype=float)
        pred = pred_df[f"pred_resid_{target}_HE_tile_features"].to_numpy(dtype=float)
        observed, _, _ = base.safe_spearman(obs, pred)
        null = []
        for _ in range(base.N_PERM):
            perm = pred.copy()
            for sample in pd.unique(samples):
                idx = np.where(samples == sample)[0]
                perm[idx] = rng.permutation(perm[idx])
            rho, _, _ = base.safe_spearman(obs, perm)
            null.append(rho)
        null = np.asarray(null, dtype=float)
        p = (1.0 + np.nansum(np.abs(null) >= abs(observed))) / (1.0 + np.isfinite(null).sum())
        rows.append({"target": target, "coord_qc_residual_target_he_rho": observed, "n_permutations": base.N_PERM, "empirical_p": float(p)})
    return pd.DataFrame(rows)


def write_outputs(model_df: pd.DataFrame, pred_df: pd.DataFrame, perm_df: pd.DataFrame) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    model_df.to_csv(OUT / "gse230424_residual_target_model_comparison.tsv", sep="\t", index=False, na_rep="NA")
    pred_df.to_csv(OUT / "gse230424_residual_target_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    perm_df.to_csv(OUT / "gse230424_residual_target_permutation.tsv", sep="\t", index=False, na_rep="NA")

    strict = model_df[(model_df["residual_target_adjustment"] == "coord_qc") & (model_df["model"] == "HE_tile_features")].set_index("target")
    dm1 = strict.loc["DM1_low_RAI_score"]
    mapk = strict.loc["MAPK_output_score"]
    tumor = strict.loc["Tumor_ZCCHC12_score"]
    hla = strict.loc["HLA_II_AP_score"]
    summary = {
        "n_spots": int(model_df["n"].max()),
        "n_samples": int(pred_df["sample"].nunique()),
        "dm1_coord_qc_residual_target_he_rho": float(dm1["sample_centered_rho"]),
        "dm1_coord_qc_residual_target_median_sample_rho": float(dm1["median_sample_rho"]),
        "dm1_coord_qc_residual_target_permutation_p": float(perm_df.set_index("target").loc["DM1_low_RAI_score", "empirical_p"]),
        "mapk_coord_qc_residual_target_he_rho": float(mapk["sample_centered_rho"]),
        "hla_coord_qc_residual_target_he_rho": float(hla["sample_centered_rho"]),
        "tumor_zcchc12_coord_qc_residual_target_he_rho": float(tumor["sample_centered_rho"]),
        "verdict": "STRICT_RESIDUAL_DM1_SIGNAL_PRESENT" if float(dm1["sample_centered_rho"]) > 0.15 else "STRICT_RESIDUAL_DM1_SIGNAL_WEAK",
    }

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.2))
    plot_df = model_df[model_df["model"] == "HE_tile_features"].copy()
    piv = plot_df.pivot(index="target", columns="residual_target_adjustment", values="sample_centered_rho").loc[STRICT_TARGETS]
    x = np.arange(len(piv.index))
    width = 0.2
    colors = {"sample_mean_only": "#266b73", "coord_only": "#7aa6a1", "qc_only": "#d5a84d", "coord_qc": "#b75f69"}
    for i, adjustment in enumerate(["sample_mean_only", "coord_only", "qc_only", "coord_qc"]):
        axes[0].bar(x + (i - 1.5) * width, piv[adjustment], width=width, color=colors[adjustment], label=adjustment)
    axes[0].axhline(0, color="#333", lw=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([t.replace("_score", "").replace("_", " ") for t in piv.index], rotation=35, ha="right")
    axes[0].set_ylabel("LOSO Spearman rho")
    axes[0].set_title("A. H&E prediction of residualized spatial targets")
    axes[0].legend(frameon=False, fontsize=9)

    sub = pred_df[["sample", "obs_resid_DM1_low_RAI_score_coord_qc", "pred_resid_DM1_low_RAI_score_HE_tile_features"]].dropna()
    for sample, part in sub.groupby("sample"):
        axes[1].scatter(part["obs_resid_DM1_low_RAI_score_coord_qc"], part["pred_resid_DM1_low_RAI_score_HE_tile_features"], s=12, alpha=0.45, label=sample)
    rho, _, _ = base.safe_spearman(sub["obs_resid_DM1_low_RAI_score_coord_qc"], sub["pred_resid_DM1_low_RAI_score_HE_tile_features"])
    axes[1].set_xlabel("Observed DM1/low-RAI residual after coord+QC")
    axes[1].set_ylabel("Predicted residual from H&E")
    axes[1].set_title(f"B. Strict residual DM1/RAI prediction (rho={rho:.3f})")
    axes[1].legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse230424_residual_target_controls.png", dpi=220)
    fig.savefig(OUT / "fig_gse230424_residual_target_controls.pdf")
    plt.close(fig)

    lines = [
        "# GSE230424 residual-target controls",
        "",
        "## Verdict",
        "",
        f"- Spots: {summary['n_spots']:,}; samples: {summary['n_samples']}.",
        "- Strict test: residualize observed spatial targets within each sample, then train leave-one-sample-out H&E models to predict the residual targets.",
        f"- DM1/low-RAI residual after coordinate+QC adjustment: H&E sample-centered rho **{summary['dm1_coord_qc_residual_target_he_rho']:.3f}**, median sample rho **{summary['dm1_coord_qc_residual_target_median_sample_rho']:.3f}**, permutation p **{summary['dm1_coord_qc_residual_target_permutation_p']:.4f}**.",
        f"- MAPK residual target rho {summary['mapk_coord_qc_residual_target_he_rho']:.3f}; HLA-II/AP {summary['hla_coord_qc_residual_target_he_rho']:.3f}; Tumor/ZCCHC12 {summary['tumor_zcchc12_coord_qc_residual_target_he_rho']:.3f}.",
        "",
        "## Interpretation",
        "",
        "This strengthens the conservative GSE230424 claim: the major DM1/RAI signal is QC/tissue-density aligned, but H&E still predicts a smaller residual component after coordinate+QC removal. This remains Paper 2 image-to-spatial-RNA support, not Paper 1 causal mechanism evidence.",
        "",
        "## Model Comparison",
        "",
        model_df.round(4).to_markdown(index=False),
        "",
        "## Permutation",
        "",
        perm_df.round(4).to_markdown(index=False),
    ]
    (OUT / "GSE230424_RESIDUAL_TARGET_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "GSE230424_RESIDUAL_TARGET_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    df, feats = load_inputs()
    model_df, pred_df = model_residual_targets(df, feats)
    perm_df = permutation_for_coord_qc(pred_df)
    summary = write_outputs(model_df, pred_df, perm_df)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
