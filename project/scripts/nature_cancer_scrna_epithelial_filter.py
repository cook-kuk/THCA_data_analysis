#!/usr/bin/env python3
"""Marker-filtered epithelial pseudo-bulk scoring for GSE232237 scRNA counts.

This is a conservative no-annotation pass. It does not claim CNV-defined
malignancy; it enriches for epithelial/tumor-compartment cells using marker
expression and then recomputes the Paper 1 lineage-state scores.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "nature_cancer_feasibility"
RAW = OUT / "raw" / "GSE232237"
FIG = OUT / "figures"

RAI_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
THYROID_NONOVERLAP = ["SLC26A4", "IYD", "DUOX1", "DUOX2", "TFF3", "HHEX", "GLIS3", "DIO2"]
TDS_TF = ["FOXE1", "NKX2-1", "PAX8", "HHEX"]
MECHANISM = ["DNMT1", "DNMT3B", "STAT3", "FOSL1", "JUNB", "TACSTD2"]

EPITHELIAL_MARKERS = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "TACSTD2", "MUC1", "CDH1"]
THYROID_MARKERS = ["TG", "TPO", "TSHR", "SLC5A5", "PAX8", "NKX2-1", "FOXE1", "DIO1", "IYD", "DUOX2"]
IMMUNE_MARKERS = ["PTPRC", "CD3D", "CD3E", "CD79A", "MS4A1", "NKG7", "LYZ", "LST1"]
STROMAL_MARKERS = ["COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "TAGLN", "PECAM1", "VWF", "RGS5"]

TARGET_GENES = sorted(set(
    RAI_8 + THYROID_NONOVERLAP + TDS_TF + MECHANISM +
    EPITHELIAL_MARKERS + THYROID_MARKERS + IMMUNE_MARKERS + STROMAL_MARKERS
))


def fnum(x):
    if x is None or pd.isna(x):
        return None
    return float(x)


def cohen_d(a, b):
    a = pd.Series(a).dropna().to_numpy(dtype=float)
    b = pd.Series(b).dropna().to_numpy(dtype=float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return np.nan if sp == 0 else float((a.mean() - b.mean()) / sp)


def mwu_p(a, b):
    a = pd.Series(a).dropna()
    b = pd.Series(b).dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def zscore_by_gene(expr: pd.DataFrame) -> pd.DataFrame:
    arr = expr.to_numpy(dtype=float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, keepdims=True, ddof=0)
    sd[sd == 0] = np.nan
    return pd.DataFrame((arr - mu) / sd, index=expr.index, columns=expr.columns)


def panel_score(z: pd.DataFrame, genes: list[str]) -> pd.Series:
    found = [g for g in genes if g in z.index]
    if not found:
        return pd.Series(np.nan, index=z.columns)
    return z.loc[found].mean(axis=0)


def mean_marker(log_expr: pd.DataFrame, genes: list[str]) -> np.ndarray:
    found = [g for g in genes if g in log_expr.index]
    if not found:
        return np.full(log_expr.shape[1], np.nan)
    return log_expr.loc[found].mean(axis=0).to_numpy(dtype=float)


def classify_cells(log_expr: pd.DataFrame, lib: np.ndarray) -> pd.DataFrame:
    epi = mean_marker(log_expr, EPITHELIAL_MARKERS)
    thyroid = mean_marker(log_expr, THYROID_MARKERS)
    immune = mean_marker(log_expr, IMMUNE_MARKERS)
    stromal = mean_marker(log_expr, STROMAL_MARKERS)
    non_epi = np.nanmax(np.vstack([immune, stromal]), axis=0)
    epithelial_enriched = (
        ((epi >= 0.12) | (thyroid >= 0.12)) &
        (np.nanmax(np.vstack([epi, thyroid]), axis=0) > non_epi)
    )
    immune_like = (immune >= stromal) & (immune > np.nanmax(np.vstack([epi, thyroid]), axis=0))
    stromal_like = (stromal > immune) & (stromal > np.nanmax(np.vstack([epi, thyroid]), axis=0))
    label = np.full(log_expr.shape[1], "other_low_marker", dtype=object)
    label[immune_like] = "immune_like"
    label[stromal_like] = "stromal_like"
    label[epithelial_enriched] = "epithelial_enriched"
    return pd.DataFrame({
        "cell_barcode": log_expr.columns,
        "library_size": lib,
        "epithelial_marker_score": epi,
        "thyroid_marker_score": thyroid,
        "immune_marker_score": immune,
        "stromal_marker_score": stromal,
        "marker_label": label,
    })


def read_target_matrix(path: Path) -> tuple[str, str, pd.DataFrame, np.ndarray]:
    name = path.name.replace(".count.tsv.gz", "")
    sample = name.split("_", 1)[1]
    histology = "ATC" if sample.startswith("AT") else "PTC"
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cells = header[1:]
        lib = np.zeros(len(cells), dtype=float)
        rows: dict[str, np.ndarray] = {}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0]
            vals = np.fromiter((float(x) for x in parts[1:]), dtype=float, count=len(cells))
            lib += vals
            if gene in TARGET_GENES:
                rows[gene] = vals
    lib_safe = lib.copy()
    lib_safe[lib_safe <= 0] = np.nan
    data = {}
    for gene in TARGET_GENES:
        counts = rows.get(gene, np.zeros(len(cells), dtype=float))
        data[gene] = np.log1p((counts / lib_safe) * 10000.0)
    log_expr = pd.DataFrame(data, index=cells).T
    log_expr.index.name = "gene_symbol"
    return sample, histology, log_expr, lib


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    files = sorted(RAW.glob("GSM*.count.tsv.gz"))

    sample_expr_rows = {}
    sample_meta = []
    cell_qc_parts = []
    cell_score_parts = []

    for path in files:
        sample, histology, log_expr, lib = read_target_matrix(path)
        cell_qc = classify_cells(log_expr, lib)
        cell_qc.insert(0, "sample_id", sample)
        cell_qc.insert(1, "histology", histology)

        epi_mask = cell_qc["marker_label"].eq("epithelial_enriched").to_numpy()
        sample_meta.append({
            "sample_id": sample,
            "histology": histology,
            "source_file": path.name,
            "n_cells": int(len(cell_qc)),
            "n_epithelial_enriched": int(epi_mask.sum()),
            "pct_epithelial_enriched": fnum(epi_mask.mean()),
            "n_immune_like": int(cell_qc["marker_label"].eq("immune_like").sum()),
            "n_stromal_like": int(cell_qc["marker_label"].eq("stromal_like").sum()),
            "median_umi_all": fnum(np.nanmedian(lib)),
            "median_umi_epithelial": fnum(np.nanmedian(lib[epi_mask])) if epi_mask.any() else None,
        })
        cell_qc_parts.append(cell_qc)

        if epi_mask.sum() >= 20:
            sample_expr_rows[sample] = log_expr.loc[sorted(set(RAI_8 + THYROID_NONOVERLAP + TDS_TF + MECHANISM)), epi_mask].mean(axis=1)
        else:
            sample_expr_rows[sample] = pd.Series(np.nan, index=sorted(set(RAI_8 + THYROID_NONOVERLAP + TDS_TF + MECHANISM)))

        # Per-cell scores inside each sample for distribution plots; gene z is within sample.
        z_cell = zscore_by_gene(log_expr.loc[sorted(set(RAI_8 + THYROID_NONOVERLAP + TDS_TF + MECHANISM))])
        cell_scores = pd.DataFrame({
            "sample_id": sample,
            "histology": histology,
            "cell_barcode": log_expr.columns,
            "marker_label": cell_qc["marker_label"].values,
            "RAI_8_cell_score": panel_score(z_cell, RAI_8).values,
            "DM1_like_cell_score": (-panel_score(z_cell, RAI_8)).values,
            "THYROID_NONOVERLAP_cell_score": panel_score(z_cell, THYROID_NONOVERLAP).values,
            "TF_collapse_cell_score": panel_score(z_cell, TDS_TF).values,
            "STAT3_AP1_DNMT_cell_score": panel_score(z_cell, MECHANISM).values,
        })
        cell_score_parts.append(cell_scores)

    meta = pd.DataFrame(sample_meta).set_index("sample_id")
    expr = pd.DataFrame(sample_expr_rows)
    expr.index.name = "gene_symbol"
    expr = expr.loc[:, meta.index]
    z_sample = zscore_by_gene(expr)
    scores = pd.DataFrame({
        "sample_id": z_sample.columns,
        "histology": meta.loc[z_sample.columns, "histology"].values,
        "n_cells": meta.loc[z_sample.columns, "n_cells"].values,
        "n_epithelial_enriched": meta.loc[z_sample.columns, "n_epithelial_enriched"].values,
        "pct_epithelial_enriched": meta.loc[z_sample.columns, "pct_epithelial_enriched"].values,
        "RAI_8_score": panel_score(z_sample, RAI_8).values,
        "DM1_like_score": (-panel_score(z_sample, RAI_8)).values,
        "THYROID_NONOVERLAP_score": panel_score(z_sample, THYROID_NONOVERLAP).values,
        "TF_collapse_score": panel_score(z_sample, TDS_TF).values,
        "STAT3_AP1_DNMT_score": panel_score(z_sample, MECHANISM).values,
    })

    expected = {
        "RAI_8_score": -1,
        "DM1_like_score": 1,
        "THYROID_NONOVERLAP_score": -1,
        "TF_collapse_score": -1,
        "STAT3_AP1_DNMT_score": 1,
    }
    contrast_rows = []
    for score, sign in expected.items():
        atc = scores.loc[scores["histology"] == "ATC", score]
        ptc = scores.loc[scores["histology"] == "PTC", score]
        d = cohen_d(atc, ptc)
        contrast_rows.append({
            "dataset": "GSE232237",
            "contrast": "ATC_vs_PTC_epithelial_enriched_pseudobulk",
            "score": score,
            "n_atc": int(atc.notna().sum()),
            "n_ptc": int(ptc.notna().sum()),
            "mean_atc": fnum(atc.mean()),
            "mean_ptc": fnum(ptc.mean()),
            "cohen_d_atc_minus_ptc": fnum(d),
            "p_mwu": fnum(mwu_p(atc, ptc)),
            "expected_sign": sign,
            "direction_match": bool(np.sign(d) == sign) if not pd.isna(d) else False,
            "note": "Marker-filtered epithelial-enriched cells; not CNV-defined malignancy.",
        })
    contrasts = pd.DataFrame(contrast_rows)

    cell_qc_df = pd.concat(cell_qc_parts, ignore_index=True)
    cell_scores_df = pd.concat(cell_score_parts, ignore_index=True)

    meta.to_csv(OUT / "gse232237_epithelial_sample_metadata.tsv", sep="\t")
    expr.to_csv(OUT / "gse232237_epithelial_target_gene_log1p_cpm10k.tsv", sep="\t")
    scores.to_csv(OUT / "gse232237_epithelial_sample_scores.tsv", sep="\t", index=False)
    contrasts.to_csv(OUT / "gse232237_epithelial_score_contrasts.tsv", sep="\t", index=False)
    cell_qc_df.to_csv(OUT / "gse232237_epithelial_cell_marker_labels.tsv.gz", sep="\t", index=False, compression="gzip")
    cell_scores_df.to_csv(OUT / "gse232237_epithelial_cell_scores.tsv.gz", sep="\t", index=False, compression="gzip")

    # Figures
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    counts = meta[["n_epithelial_enriched", "n_immune_like", "n_stromal_like"]].copy()
    counts.index = [f"{idx}\n{meta.loc[idx, 'histology']}" for idx in counts.index]
    counts.plot(kind="bar", stacked=True, ax=axes[0], color=["#506f59", "#2e5b87", "#a97922"], width=0.85)
    axes[0].set_title("GSE232237 marker-filtered cell compartments")
    axes[0].set_ylabel("Cells")
    axes[0].tick_params(axis="x", labelsize=7, rotation=55)
    axes[0].legend(fontsize=7)

    score_order = ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TF_collapse_score", "STAT3_AP1_DNMT_score"]
    x = np.arange(len(score_order))
    ptc_means = [scores.loc[scores["histology"] == "PTC", s].mean() for s in score_order]
    atc_means = [scores.loc[scores["histology"] == "ATC", s].mean() for s in score_order]
    axes[1].bar(x - 0.18, ptc_means, width=0.36, label="PTC", color="#2e5b87")
    axes[1].bar(x + 0.18, atc_means, width=0.36, label="ATC", color="#8b2635")
    axes[1].axhline(0, color="#999", linewidth=0.7)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([s.replace("_score", "").replace("_", "\n") for s in score_order], fontsize=7)
    axes[1].set_title("Epithelial-enriched pseudo-bulk score means")
    axes[1].legend(fontsize=8)

    epi_cells = cell_scores_df[cell_scores_df["marker_label"] == "epithelial_enriched"]
    vals = [
        epi_cells.loc[epi_cells["histology"] == "PTC", "DM1_like_cell_score"].dropna().sample(
            min(5000, (epi_cells["histology"] == "PTC").sum()), random_state=1
        ),
        epi_cells.loc[epi_cells["histology"] == "ATC", "DM1_like_cell_score"].dropna().sample(
            min(5000, (epi_cells["histology"] == "ATC").sum()), random_state=1
        ),
    ]
    axes[2].boxplot(vals, labels=["PTC epi", "ATC epi"], patch_artist=True,
                    boxprops={"facecolor": "#f7f5ef"}, medianprops={"color": "#8b2635", "linewidth": 2})
    axes[2].set_title("Cell-level DM1_like within epithelial-enriched cells")
    axes[2].set_ylabel("Within-sample cell score")
    fig.tight_layout()
    fig.savefig(FIG / "gse232237_epithelial_filter_validation.png")
    plt.close(fig)

    summary = {
        "dataset": "GSE232237",
        "analysis": "marker-filtered epithelial-enriched pseudo-bulk",
        "n_samples": int(len(scores)),
        "n_ptc": int((scores["histology"] == "PTC").sum()),
        "n_atc": int((scores["histology"] == "ATC").sum()),
        "n_cells_total": int(meta["n_cells"].sum()),
        "n_epithelial_enriched_total": int(meta["n_epithelial_enriched"].sum()),
        "pct_epithelial_enriched_total": fnum(meta["n_epithelial_enriched"].sum() / meta["n_cells"].sum()),
        "direction_matches": int(contrasts["direction_match"].sum()),
        "direction_cells": int(len(contrasts)),
        "dm1_cohen_d_atc_minus_ptc": fnum(contrasts.loc[contrasts["score"] == "DM1_like_score", "cohen_d_atc_minus_ptc"].iloc[0]),
        "rai8_cohen_d_atc_minus_ptc": fnum(contrasts.loc[contrasts["score"] == "RAI_8_score", "cohen_d_atc_minus_ptc"].iloc[0]),
        "nonoverlap_cohen_d_atc_minus_ptc": fnum(contrasts.loc[contrasts["score"] == "THYROID_NONOVERLAP_score", "cohen_d_atc_minus_ptc"].iloc[0]),
        "caveat": "Marker-enriched epithelial compartment, not CNV-inferred malignant cells.",
    }
    (OUT / "gse232237_epithelial_filter_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
