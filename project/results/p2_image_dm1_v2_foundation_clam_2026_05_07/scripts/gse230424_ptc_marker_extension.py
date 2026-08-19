#!/usr/bin/env python3
"""GSE230424 extension using article Supplementary Table S6 PTC markers."""
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
OUT = BASE_OUT / "ptc_marker_extension"

spec = importlib.util.spec_from_file_location("gse230424_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(base)

PTC_MARKERS = [
    "ZCCHC12", "CHI3L1", "LAMB3", "IGSF1", "LRRK2", "CITED1", "EMILIN2", "LRP4", "PCSK2", "HMGA2", "MRC2", "CAMK2N1",
    "PROS1", "LGALS3", "TUSC3", "PDZK1IP1", "TGFA", "DPP4", "EPHA4", "PDE5A", "QPCT", "ITGA6", "GALE", "ENTPD1",
    "PLXNC1", "PON2", "FN1", "EPS8", "RAB27A", "NT5E", "ARHGAP36", "ATP11A", "ETV4", "STK32A", "FAM20A", "NRP2",
    "DCBLD2", "DOCK9", "IGF2R", "TGFBR1", "GGCT", "IGF2BP2", "POSTN", "CYP1B1", "MALL", "NOD1", "CDH6", "ABCC3",
    "PSD3", "CRABP2", "LIPH", "CLDN16", "ADORA1", "CHST2", "PRICKLE1", "ETV5", "SLC17A5",
]


def module_z(df: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [g for g in genes if f"g_{g}" in df.columns and np.nanmax(df[f"g_{g}"].to_numpy(dtype=float)) > 0]
    vals = []
    for gene in present:
        v = df[f"g_{gene}"].to_numpy(dtype=float)
        sd = np.nanstd(v)
        vals.append((v - np.nanmean(v)) / (sd + 1e-9) if sd > 0 else np.zeros_like(v))
    if not vals:
        return np.full(len(df), np.nan), present
    return np.nanmean(np.vstack(vals), axis=0), present


def load_ptc_marker_table() -> pd.DataFrame:
    rows = []
    for prefix in base.sample_prefixes():
        sample = prefix.split("_")[-1]
        print(f"[GSE230424 PTC markers] load {sample}")
        expr = base.stream_target_matrix(prefix, PTC_MARKERS)
        pos = base.read_positions(prefix)
        meta = expr.merge(pos, on="barcode", how="left")
        meta.insert(0, "sample", sample)
        meta.insert(1, "disease_group", base.DISEASE_BY_SAMPLE.get(sample, "unknown"))
        meta.insert(2, "geo_accession", prefix.split("_")[0])
        rows.append(meta)
    df = pd.concat(rows, ignore_index=True)
    score, present = module_z(df, PTC_MARKERS)
    df["PTC_specific_marker_score"] = score
    df["n_PTC_marker_genes"] = len(present)
    df.to_csv(OUT / "gse230424_ptc_marker_spot_scores.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    (OUT / "gse230424_ptc_marker_genes_used.txt").write_text("\n".join(present) + "\n")
    return df


def run_models(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feats = pd.read_csv(BASE_OUT / "gse230424_he_tile_features.tsv.gz", sep="\t")
    df = base.spatial_smooth(df, ["PTC_specific_marker_score"])
    merged = df.merge(feats, on=["sample", "barcode", "array_row", "array_col"], how="inner")
    samples = merged["sample"].astype(str).to_numpy()
    disease = merged["disease_group"].astype(str).to_numpy()
    he_cols = [c for c in feats.columns if c.startswith("r")]
    coord = merged[["array_row", "array_col"]].to_numpy(dtype=float)
    qc = np.column_stack([
        np.log1p(merged["total_counts"].to_numpy(dtype=float)),
        np.log1p(merged["n_genes_by_counts"].to_numpy(dtype=float)),
        merged["pct_counts_mt"].to_numpy(dtype=float),
    ])
    feature_sets = {
        "HE_tile_features": (merged[he_cols].to_numpy(dtype=float), True),
        "Coord_only": (coord, False),
        "QC_only": (qc, False),
        "Coord_QC": (np.column_stack([coord, qc]), False),
    }
    y = merged[f"PTC_specific_marker_score_smooth{base.SMOOTH_K}"].to_numpy(dtype=float)
    pred_df = merged[["sample", "disease_group", "geo_accession", "barcode", "array_row", "array_col", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    pred_df["obs_PTC_specific_marker_score_smooth8"] = y
    rows = []
    for model_name, (x, use_pca) in feature_sets.items():
        pred = base.loso_predict(x, y, samples, use_pca=use_pca)
        pred_df[f"pred_PTC_specific_marker_score_{model_name}"] = pred
        rho, p, n = base.safe_spearman(y, pred)
        crho, cp, cn = base.safe_spearman(base.center_by_group(y, samples), base.center_by_group(pred, samples))
        per_sample = []
        for sample, idx_raw in merged.groupby("sample").groups.items():
            idx = np.asarray(list(idx_raw))
            srho, _, _ = base.safe_spearman(y[idx], pred[idx])
            per_sample.append(srho)
        rows.append({
            "target": "PTC_specific_marker_score",
            "model": model_name,
            "n": n,
            "pooled_rho": rho,
            "pooled_p": p,
            "sample_centered_rho": crho,
            "sample_centered_p": cp,
            "median_sample_rho": float(np.nanmedian(per_sample)),
            "samples_rho_gt_0_2": int(np.nansum(np.asarray(per_sample) > 0.2)),
        })

    # Strict residual target after coordinate+QC removal.
    y_resid = base.residualize_within_sample(y, np.column_stack([coord, qc]), samples)
    pred_resid = base.loso_predict(merged[he_cols].to_numpy(dtype=float), y_resid, samples, use_pca=True)
    rrho, rp, rn = base.safe_spearman(y_resid, pred_resid)
    rcrho, rcp, rcn = base.safe_spearman(base.center_by_group(y_resid, samples), base.center_by_group(pred_resid, samples))
    residual_df = pd.DataFrame([{
        "target": "PTC_specific_marker_score",
        "residual_target_adjustment": "coord_qc",
        "model": "HE_tile_features",
        "n": rn,
        "pooled_rho": rrho,
        "sample_centered_rho": rcrho,
        "p": rp,
    }])
    pred_df["obs_resid_PTC_specific_marker_score_coord_qc"] = y_resid
    pred_df["pred_resid_PTC_specific_marker_score_HE_tile_features"] = pred_resid

    sample_summary = merged[["sample", "disease_group", "PTC_specific_marker_score"]].groupby(["sample", "disease_group"], as_index=False).mean()
    group_summary = sample_summary.groupby("disease_group", as_index=False).mean(numeric_only=True)
    group_summary["n_samples"] = sample_summary.groupby("disease_group")["sample"].count().to_numpy()
    sample_summary.to_csv(OUT / "gse230424_ptc_marker_sample_summary.tsv", sep="\t", index=False, na_rep="NA")
    group_summary.to_csv(OUT / "gse230424_ptc_marker_group_summary.tsv", sep="\t", index=False, na_rep="NA")
    return pd.DataFrame(rows), pred_df, residual_df


def write_outputs(model_df: pd.DataFrame, pred_df: pd.DataFrame, residual_df: pd.DataFrame) -> dict:
    model_df.to_csv(OUT / "gse230424_ptc_marker_model_comparison.tsv", sep="\t", index=False, na_rep="NA")
    pred_df.to_csv(OUT / "gse230424_ptc_marker_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    residual_df.to_csv(OUT / "gse230424_ptc_marker_residual_target.tsv", sep="\t", index=False, na_rep="NA")
    he = model_df.set_index("model").loc["HE_tile_features"]
    qc = model_df.set_index("model").loc["QC_only"]
    coord_qc = model_df.set_index("model").loc["Coord_QC"]
    resid = residual_df.iloc[0]
    group_summary = pd.read_csv(OUT / "gse230424_ptc_marker_group_summary.tsv", sep="\t")
    summary = {
        "n_spots": int(he["n"]),
        "n_samples": int(pred_df["sample"].nunique()),
        "n_ptc_marker_genes": len([g for g in (OUT / "gse230424_ptc_marker_genes_used.txt").read_text().splitlines() if g.strip()]),
        "he_sample_centered_rho": float(he["sample_centered_rho"]),
        "qc_only_sample_centered_rho": float(qc["sample_centered_rho"]),
        "coord_qc_sample_centered_rho": float(coord_qc["sample_centered_rho"]),
        "coord_qc_residual_target_he_rho": float(resid["sample_centered_rho"]),
        "ptc_ht_mean": float(group_summary.loc[group_summary["disease_group"] == "PTC+HT", "PTC_specific_marker_score"].iloc[0]),
        "ht_mean": float(group_summary.loc[group_summary["disease_group"] == "HT", "PTC_specific_marker_score"].iloc[0]),
    }
    summary["ptc_ht_minus_ht_mean"] = summary["ptc_ht_mean"] - summary["ht_mean"]

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))
    axes[0].bar(model_df["model"], model_df["sample_centered_rho"], color=["#266b73", "#9aa4ad", "#777", "#777"])
    axes[0].axhline(0, color="#333", lw=0.8)
    axes[0].set_ylabel("Sample-centered LOSO rho")
    axes[0].set_title("A. PTC marker H&E prediction")
    axes[0].tick_params(axis="x", rotation=25)

    sample_summary = pd.read_csv(OUT / "gse230424_ptc_marker_sample_summary.tsv", sep="\t")
    colors = {"PTC+HT": "#b75f69", "HT": "#266b73"}
    axes[1].bar(sample_summary["sample"], sample_summary["PTC_specific_marker_score"], color=[colors[g] for g in sample_summary["disease_group"]])
    axes[1].axhline(0, color="#333", lw=0.8)
    axes[1].set_title("B. Sample mean PTC marker score")
    axes[1].set_ylabel("Mean module z")

    sub = pred_df[["sample", "disease_group", "obs_PTC_specific_marker_score_smooth8", "pred_PTC_specific_marker_score_HE_tile_features"]].dropna()
    for group, part in sub.groupby("disease_group"):
        axes[2].scatter(part["obs_PTC_specific_marker_score_smooth8"], part["pred_PTC_specific_marker_score_HE_tile_features"], s=12, alpha=0.45, color=colors[group], label=group)
    rho, _, _ = base.safe_spearman(sub["obs_PTC_specific_marker_score_smooth8"], sub["pred_PTC_specific_marker_score_HE_tile_features"])
    axes[2].set_xlabel("Observed PTC marker spatial score")
    axes[2].set_ylabel("Predicted from H&E")
    axes[2].set_title(f"C. Pooled scatter (rho={rho:.3f})")
    axes[2].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse230424_ptc_marker_extension.png", dpi=220)
    fig.savefig(OUT / "fig_gse230424_ptc_marker_extension.pdf")
    plt.close(fig)

    lines = [
        "# GSE230424 PTC-marker extension",
        "",
        "## Verdict",
        "",
        f"- PTC marker genes used: {summary['n_ptc_marker_genes']} from article Supplementary Table S6.",
        f"- H&E sample-centered rho for the PTC-marker spatial score: **{summary['he_sample_centered_rho']:.3f}**.",
        f"- QC-only sample-centered rho: **{summary['qc_only_sample_centered_rho']:.3f}**; coord+QC rho: **{summary['coord_qc_sample_centered_rho']:.3f}**.",
        f"- Strict coordinate+QC residual-target H&E rho: **{summary['coord_qc_residual_target_he_rho']:.3f}**.",
        f"- PTC+HT minus HT sample-mean PTC marker score: **{summary['ptc_ht_minus_ht_mean']:.3f}** (n=2/group, descriptive only).",
        "",
        "## Interpretation",
        "",
        "The article's own PTC-specific spatial marker program is H&E-predictable, but QC/tissue-density baselines remain strong. Treat this as a descriptive extension of Paper 2 image-to-spatial-RNA support, not a robust disease-group validation.",
        "",
        "## Model Comparison",
        "",
        model_df.round(4).to_markdown(index=False),
        "",
        "## Residual Target",
        "",
        residual_df.round(4).to_markdown(index=False),
    ]
    (OUT / "GSE230424_PTC_MARKER_EXTENSION_REPORT.md").write_text("\n".join(lines) + "\n")
    (OUT / "GSE230424_PTC_MARKER_EXTENSION_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spot_path = OUT / "gse230424_ptc_marker_spot_scores.tsv.gz"
    if spot_path.exists():
        df = pd.read_csv(spot_path, sep="\t")
    else:
        df = load_ptc_marker_table()
    model_df, pred_df, residual_df = run_models(df)
    summary = write_outputs(model_df, pred_df, residual_df)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
