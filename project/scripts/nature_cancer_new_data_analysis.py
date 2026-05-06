#!/usr/bin/env python3
"""Additional public-data validation for Paper 1 Nature Cancer feasibility board.

Inputs:
  - GSE146003 EPIC methylation signal-intensity matrix (10 ATC, 4 normal thyroid)
  - GSE60542 GPL570 series matrix (PTC, nodal metastasis, matched controls)

Outputs are intentionally kept separate from the original validation bundle.
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "nature_cancer_feasibility"
RAW = OUT / "raw"
FIG = OUT / "figures"
EXT_SCRIPT = ROOT / "results" / "p_external_expression_validation" / "scripts"
sys.path.insert(0, str(EXT_SCRIPT))

from run_external_validation import (  # noqa: E402
    MECHANISM,
    RAI_8,
    TDS_TF,
    THYROID_NONOVERLAP,
    cohen_d,
    load_gpl570_probe2gene,
    panel_score,
    parse_series_matrix,
    probe_to_gene,
    zscore_within_dataset,
)

plt.rcParams.update({
    "figure.dpi": 140,
    "savefig.dpi": 160,
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def fnum(x: float | int | None) -> float | None:
    if x is None:
        return None
    if pd.isna(x):
        return None
    return float(x)


def p_mwu(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    try:
        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def read_gse146003() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = RAW / "GSE146003" / "GSE146003_Matrix_signal_intensities.txt.gz"
    df = pd.read_csv(path, sep="\t")
    sample_names = []
    for col in df.columns:
        if col == "ID_REF":
            continue
        m = re.match(r"(.+?) (Unmethylated Signal|Methylated signal|Detection Pval)$", col)
        if m:
            sample_names.append(m.group(1))
    sample_names = sorted(set(sample_names), key=sample_names.index)

    unm = pd.DataFrame(index=df["ID_REF"])
    met = pd.DataFrame(index=df["ID_REF"])
    pval = pd.DataFrame(index=df["ID_REF"])
    for sample in sample_names:
        unm[sample] = pd.to_numeric(df[f"{sample} Unmethylated Signal"], errors="coerce").to_numpy()
        met[sample] = pd.to_numeric(df[f"{sample} Methylated signal"], errors="coerce").to_numpy()
        pval[sample] = pd.to_numeric(df[f"{sample} Detection Pval"], errors="coerce").to_numpy()
    beta = met / (met + unm + 100.0)
    meta = pd.DataFrame({
        "sample_id": sample_names,
        "group": ["normal_thyroid" if s.startswith("NT") else "ATC" for s in sample_names],
    }).set_index("sample_id")
    return meta, beta, pval, df[["ID_REF"]]


def analyze_gse146003() -> dict:
    meta, beta, pval, _ = read_gse146003()
    normal = meta.index[meta["group"] == "normal_thyroid"].tolist()
    atc = meta.index[meta["group"] == "ATC"].tolist()

    sample_qc = pd.DataFrame({
        "sample_id": beta.columns,
        "group": meta.loc[beta.columns, "group"].values,
        "n_probes": beta.notna().sum(axis=0).values,
        "mean_beta": beta.mean(axis=0).values,
        "median_beta": beta.median(axis=0).values,
        "detection_fail_rate_p001": (pval > 0.01).mean(axis=0).values,
        "detection_fail_rate_p005": (pval > 0.05).mean(axis=0).values,
    })
    sample_qc.to_csv(OUT / "gse146003_sample_qc.tsv", sep="\t", index=False)

    probe_fail_rate = (pval > 0.01).mean(axis=1)
    qc_mask = probe_fail_rate <= 0.2
    beta_qc = beta.loc[qc_mask]

    mean_normal = beta_qc[normal].mean(axis=1)
    mean_atc = beta_qc[atc].mean(axis=1)
    delta = mean_atc - mean_normal
    sd_normal = beta_qc[normal].std(axis=1, ddof=1)
    sd_atc = beta_qc[atc].std(axis=1, ddof=1)
    pooled = np.sqrt(((len(normal) - 1) * sd_normal.pow(2) + (len(atc) - 1) * sd_atc.pow(2)) / (len(normal) + len(atc) - 2))
    d = delta / pooled.replace(0, np.nan)

    diff = pd.DataFrame({
        "probe_id": beta_qc.index,
        "mean_beta_normal": mean_normal.values,
        "mean_beta_atc": mean_atc.values,
        "delta_beta_atc_minus_normal": delta.values,
        "cohen_d_atc_minus_normal": d.values,
        "detection_fail_rate_p001": probe_fail_rate.loc[beta_qc.index].values,
    }).sort_values("delta_beta_atc_minus_normal", key=lambda s: s.abs(), ascending=False)
    diff.to_csv(OUT / "gse146003_probe_delta_beta.tsv.gz", sep="\t", index=False, compression="gzip")

    top = pd.read_csv(ROOT / "results" / "ml" / "methylation_top_probes.tsv", sep="\t")
    overlap = top.merge(diff, on="probe_id", how="inner")
    overlap["direction_match_tcga"] = np.sign(overlap["delta_beta"]) == np.sign(overlap["delta_beta_atc_minus_normal"])
    overlap["abs_delta_ratio_gse146003_to_tcga"] = overlap["delta_beta_atc_minus_normal"].abs() / overlap["delta_beta"].abs()
    overlap["p_mwu_gse146003"] = [
        p_mwu(beta_qc.loc[row.probe_id, atc], beta_qc.loc[row.probe_id, normal])
        for row in overlap.itertuples(index=False)
    ]
    overlap.to_csv(OUT / "gse146003_tcga_top_probe_overlap.tsv", sep="\t", index=False)

    rho = float("nan")
    rho_p = float("nan")
    if len(overlap) >= 3:
        rho, rho_p = spearmanr(overlap["delta_beta"], overlap["delta_beta_atc_minus_normal"], nan_policy="omit")

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    groups = ["normal_thyroid", "ATC"]
    vals = [sample_qc.loc[sample_qc["group"] == g, "mean_beta"].values for g in groups]
    axes[0].boxplot(vals, labels=["Normal", "ATC"], patch_artist=True,
                    boxprops={"facecolor": "#e6f2ec"}, medianprops={"color": "#8b2635", "linewidth": 2})
    for i, arr in enumerate(vals, start=1):
        axes[0].scatter(np.full(len(arr), i) + np.linspace(-0.04, 0.04, len(arr)), arr, s=24, color="#2e5b87", zorder=3)
    axes[0].set_title("GSE146003 global mean beta")
    axes[0].set_ylabel("Mean beta")

    axes[1].axhline(0, color="#999", linewidth=0.8)
    axes[1].axvline(0, color="#999", linewidth=0.8)
    axes[1].scatter(overlap["delta_beta"], overlap["delta_beta_atc_minus_normal"],
                    s=42, color=np.where(overlap["direction_match_tcga"], "#506f59", "#8b2635"))
    for _, row in overlap.head(8).iterrows():
        axes[1].text(row["delta_beta"], row["delta_beta_atc_minus_normal"], row["probe_id"], fontsize=6)
    axes[1].set_title(f"TCGA top CpGs vs GSE146003 (rho={rho:.2f})")
    axes[1].set_xlabel("TCGA THCA tumor-normal delta beta")
    axes[1].set_ylabel("GSE146003 ATC-normal delta beta")
    fig.tight_layout()
    fig.savefig(FIG / "gse146003_methylation_validation.png")
    plt.close(fig)

    summary = {
        "dataset": "GSE146003",
        "n_samples": int(len(meta)),
        "n_normal": int(len(normal)),
        "n_atc": int(len(atc)),
        "n_probes_raw": int(len(beta)),
        "n_probes_qc": int(len(beta_qc)),
        "median_sample_detection_fail_rate_p001": fnum(sample_qc["detection_fail_rate_p001"].median()),
        "global_mean_beta_normal": fnum(sample_qc.loc[sample_qc["group"] == "normal_thyroid", "mean_beta"].mean()),
        "global_mean_beta_atc": fnum(sample_qc.loc[sample_qc["group"] == "ATC", "mean_beta"].mean()),
        "global_mean_beta_p": fnum(p_mwu(
            sample_qc.loc[sample_qc["group"] == "ATC", "mean_beta"],
            sample_qc.loc[sample_qc["group"] == "normal_thyroid", "mean_beta"],
        )),
        "n_tcga_top_probes": int(len(top)),
        "n_tcga_top_overlap_qc": int(len(overlap)),
        "tcga_top_direction_match": int(overlap["direction_match_tcga"].sum()),
        "tcga_top_direction_match_rate": fnum(overlap["direction_match_tcga"].mean()),
        "tcga_delta_spearman_rho": fnum(rho),
        "tcga_delta_spearman_p": fnum(rho_p),
        "median_abs_delta_top_overlap": fnum(overlap["delta_beta_atc_minus_normal"].abs().median()),
    }
    return summary


def label_gse60542(meta: pd.DataFrame) -> pd.DataFrame:
    out = meta.copy()
    title = out["Sample_title"].fillna("")
    groups = []
    nodal = []
    patient = []
    for t in title:
        tl = t.lower()
        if "normal lymph node" in tl:
            groups.append("normal_lymph_node")
        elif "normal thyroid" in tl:
            groups.append("normal_thyroid")
        elif "papillary thyroid carcinoma" in tl:
            groups.append("PTC_primary")
        elif "lymph node metastasis" in tl:
            groups.append("LNM")
        elif "pleural metastasis" in tl:
            groups.append("PLM")
        elif "recurrence" in tl:
            groups.append("recurrence")
        else:
            groups.append("other")
        nodal.append("N1" if ",N1" in t else ("N0" if ",N0" in t else "unknown"))
        patient.append(t.split(",", 1)[0].split("-", 1)[0])
    out["histology_clean"] = groups
    out["nodal_status"] = nodal
    out["patient_token"] = patient
    return out


def analyze_gse60542() -> dict:
    path = RAW / "GSE60542" / "GSE60542_series_matrix.txt.gz"
    meta, expr_probe = parse_series_matrix(path)
    meta = label_gse60542(meta)
    probe2gene = load_gpl570_probe2gene(ROOT / "results" / "p_external_expression_validation" / "raw" / "GPL570" / "GPL570.annot.gz")
    gene = probe_to_gene(expr_probe, probe2gene)
    z = zscore_within_dataset(gene)

    rai, rai_found, rai_missing = panel_score(z, RAI_8)
    non, non_found, non_missing = panel_score(z, THYROID_NONOVERLAP)
    tf, tf_found, tf_missing = panel_score(z, TDS_TF)
    mech, mech_found, mech_missing = panel_score(z, MECHANISM)
    scores = pd.DataFrame({
        "sample_id": z.columns,
        "title": meta.loc[z.columns, "Sample_title"].values,
        "histology_clean": meta.loc[z.columns, "histology_clean"].values,
        "nodal_status": meta.loc[z.columns, "nodal_status"].values,
        "patient_token": meta.loc[z.columns, "patient_token"].values,
        "RAI_8_score": rai.loc[z.columns].values,
        "DM1_like_score": -rai.loc[z.columns].values,
        "THYROID_NONOVERLAP_score": non.loc[z.columns].values,
        "TF_collapse_score": tf.loc[z.columns].values,
        "STAT3_AP1_DNMT_score": mech.loc[z.columns].values,
    })
    scores.to_csv(OUT / "gse60542_sample_scores.tsv", sep="\t", index=False)

    contrasts = [
        ("PTC_primary_vs_normal_thyroid", "PTC_primary", "normal_thyroid", "PTC primary vs normal thyroid"),
        ("LNM_vs_PTC_primary", "LNM", "PTC_primary", "Lymph-node metastasis vs primary PTC; lymphoid confounded"),
        ("PTC_N1_vs_PTC_N0", "PTC_primary:N1", "PTC_primary:N0", "N1 primary PTC vs N0 primary PTC"),
        ("LNM_vs_normal_lymph_node", "LNM", "normal_lymph_node", "LNM vs normal lymph node; tissue-composition control"),
        ("metastatic_or_recurrent_vs_PTC_primary", "MET", "PTC_primary", "LNM/PLM/recurrence vs primary PTC; tissue-composition confounded"),
    ]
    expected = {
        "RAI_8_score": -1,
        "DM1_like_score": 1,
        "THYROID_NONOVERLAP_score": -1,
        "TF_collapse_score": -1,
        "STAT3_AP1_DNMT_score": 1,
    }

    def mask(spec: str) -> pd.Series:
        if spec == "MET":
            return scores["histology_clean"].isin(["LNM", "PLM", "recurrence"])
        if ":" in spec:
            group, nod = spec.split(":", 1)
            return (scores["histology_clean"] == group) & (scores["nodal_status"] == nod)
        return scores["histology_clean"] == spec

    rows = []
    for contrast_id, a_spec, b_spec, note in contrasts:
        ma = mask(a_spec)
        mb = mask(b_spec)
        for score in expected:
            a = scores.loc[ma, score]
            b = scores.loc[mb, score]
            d = cohen_d(a.to_numpy(dtype=float), b.to_numpy(dtype=float))
            rows.append({
                "dataset": "GSE60542",
                "contrast": contrast_id,
                "score": score,
                "group_a": a_spec,
                "group_b": b_spec,
                "n_a": int(a.notna().sum()),
                "n_b": int(b.notna().sum()),
                "mean_a": fnum(a.mean()),
                "mean_b": fnum(b.mean()),
                "cohen_d_a_minus_b": fnum(d),
                "p_mwu": fnum(p_mwu(a, b)),
                "expected_sign_for_a_minus_b": expected[score],
                "direction_match": bool(np.sign(d) == expected[score]) if not pd.isna(d) else False,
                "note": note,
            })
    tests = pd.DataFrame(rows)
    tests.to_csv(OUT / "gse60542_score_contrasts.tsv", sep="\t", index=False)

    keep = scores[["DM1_like_score", "THYROID_NONOVERLAP_score", "TF_collapse_score"]].dropna()
    rho_dm1_non, p_dm1_non = spearmanr(keep["DM1_like_score"], keep["THYROID_NONOVERLAP_score"])
    rho_dm1_tf, p_dm1_tf = spearmanr(keep["DM1_like_score"], keep["TF_collapse_score"])
    pd.DataFrame([
        {"dataset": "GSE60542", "x": "DM1_like_score", "y": "THYROID_NONOVERLAP_score", "rho": rho_dm1_non, "p": p_dm1_non, "n": len(keep)},
        {"dataset": "GSE60542", "x": "DM1_like_score", "y": "TF_collapse_score", "rho": rho_dm1_tf, "p": p_dm1_tf, "n": len(keep)},
    ]).to_csv(OUT / "gse60542_score_spearman.tsv", sep="\t", index=False)

    group_order = ["normal_thyroid", "PTC_primary", "LNM", "normal_lymph_node", "recurrence", "PLM"]
    score_order = ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TF_collapse_score"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), sharex=True)
    colors = ["#506f59", "#2e5b87", "#8b2635", "#a97922", "#7a546a", "#333333"]
    for ax, score in zip(axes.ravel(), score_order):
        vals = [scores.loc[scores["histology_clean"] == g, score].dropna().values for g in group_order]
        ax.boxplot(vals, labels=[g.replace("_", "\n") for g in group_order], patch_artist=True,
                   boxprops={"facecolor": "#f7f5ef"}, medianprops={"color": "#8b2635", "linewidth": 2})
        for i, arr in enumerate(vals, start=1):
            if len(arr):
                jitter = np.linspace(-0.08, 0.08, len(arr))
                ax.scatter(np.full(len(arr), i) + jitter, arr, s=12, color=colors[i - 1], alpha=0.75)
        ax.axhline(0, color="#999", linewidth=0.7)
        ax.set_title(score)
        ax.tick_params(axis="x", labelrotation=0, labelsize=7)
    fig.suptitle("GSE60542 lineage-state scores; interpret nodal metastasis contrasts with lymphoid confounding")
    fig.tight_layout()
    fig.savefig(FIG / "gse60542_lineage_scores.png")
    plt.close(fig)

    counts = scores["histology_clean"].value_counts().to_dict()
    headline = tests[tests["contrast"].isin(["PTC_primary_vs_normal_thyroid", "LNM_vs_PTC_primary", "PTC_N1_vs_PTC_N0"])]
    summary = {
        "dataset": "GSE60542",
        "n_samples": int(len(scores)),
        "group_counts": {str(k): int(v) for k, v in counts.items()},
        "gene_rows_after_collapse": int(gene.shape[0]),
        "rai_found": len(rai_found),
        "rai_missing": rai_missing,
        "thyroid_nonoverlap_found": len(non_found),
        "thyroid_nonoverlap_missing": non_missing,
        "tf_found": len(tf_found),
        "tf_missing": tf_missing,
        "mechanism_found": len(mech_found),
        "mechanism_missing": mech_missing,
        "dm1_vs_nonoverlap_rho": fnum(rho_dm1_non),
        "dm1_vs_nonoverlap_p": fnum(p_dm1_non),
        "dm1_vs_tf_rho": fnum(rho_dm1_tf),
        "dm1_vs_tf_p": fnum(p_dm1_tf),
        "headline_direction_matches": int(headline["direction_match"].sum()),
        "headline_direction_cells": int(len(headline)),
    }
    return summary


def analyze_gse137697() -> dict:
    raw_dir = RAW / "GSE137697" / "extracted"
    files = sorted(raw_dir.glob("GSM*.txt.gz"))
    matrices = []
    meta_rows = []
    for path in files:
        df = pd.read_csv(path, sep="\t")
        value_col = [c for c in df.columns if c != "Symbol"][0]
        sample = path.name.replace(".txt.gz", "")
        gsm = sample.split("_", 1)[0]
        title = value_col
        treatment = "DMSO" if "DMSO" in title else ("THZ531_100nM" if title.startswith("100") else "THZ531_400nM")
        rep = "rep1" if "rep1" in title else "rep2"
        one = df[["Symbol", value_col]].copy()
        one["Symbol"] = one["Symbol"].astype(str).str.strip()
        one[value_col] = pd.to_numeric(one[value_col], errors="coerce")
        one = one.dropna(subset=["Symbol"]).groupby("Symbol", as_index=True)[value_col].mean()
        one.name = gsm
        matrices.append(one)
        meta_rows.append({"sample_id": gsm, "source_file": path.name, "title": title, "treatment": treatment, "replicate": rep})

    expr = pd.concat(matrices, axis=1).fillna(0.0)
    meta = pd.DataFrame(meta_rows).set_index("sample_id").loc[expr.columns]
    log_expr = np.log2(expr + 1.0)
    log_expr.index.name = "gene_symbol"
    log_expr.to_csv(OUT / "gse137697_log2_rpkm_matrix.tsv.gz", sep="\t", compression="gzip")
    meta.to_csv(OUT / "gse137697_sample_metadata.tsv", sep="\t")

    z = zscore_within_dataset(log_expr)
    rai, rai_found, rai_missing = panel_score(z, RAI_8)
    non, non_found, non_missing = panel_score(z, THYROID_NONOVERLAP)
    tf, tf_found, tf_missing = panel_score(z, TDS_TF)
    mech, mech_found, mech_missing = panel_score(z, MECHANISM)
    scores = pd.DataFrame({
        "sample_id": z.columns,
        "treatment": meta.loc[z.columns, "treatment"].values,
        "replicate": meta.loc[z.columns, "replicate"].values,
        "RAI_8_score": rai.loc[z.columns].values,
        "DM1_like_score": -rai.loc[z.columns].values,
        "THYROID_NONOVERLAP_score": non.loc[z.columns].values,
        "TF_collapse_score": tf.loc[z.columns].values,
        "STAT3_AP1_DNMT_score": mech.loc[z.columns].values,
    })
    scores.to_csv(OUT / "gse137697_sample_scores.tsv", sep="\t", index=False)

    score_cols = ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TF_collapse_score", "STAT3_AP1_DNMT_score"]
    contrast_rows = []
    for treatment in ["THZ531_100nM", "THZ531_400nM"]:
        for score in score_cols:
            a = scores.loc[scores["treatment"] == treatment, score]
            b = scores.loc[scores["treatment"] == "DMSO", score]
            contrast_rows.append({
                "dataset": "GSE137697",
                "contrast": f"{treatment}_vs_DMSO",
                "score": score,
                "n_treatment": int(a.notna().sum()),
                "n_dmso": int(b.notna().sum()),
                "mean_treatment": fnum(a.mean()),
                "mean_dmso": fnum(b.mean()),
                "delta_treatment_minus_dmso": fnum(a.mean() - b.mean()),
                "cohen_d_treatment_minus_dmso": fnum(cohen_d(a.to_numpy(dtype=float), b.to_numpy(dtype=float))),
                "p_mwu": fnum(p_mwu(a, b)),
                "note": "n=2 per arm; use as perturbation-direction sidecar only.",
            })
    score_contrasts = pd.DataFrame(contrast_rows)
    score_contrasts.to_csv(OUT / "gse137697_score_contrasts.tsv", sep="\t", index=False)

    genes = sorted(set(RAI_8 + THYROID_NONOVERLAP + TDS_TF + MECHANISM))
    gene_rows = []
    for gene in genes:
        if gene not in log_expr.index:
            gene_rows.append({"gene_symbol": gene, "present": False})
            continue
        dmso = log_expr.loc[gene, meta.index[meta["treatment"] == "DMSO"]]
        row = {"gene_symbol": gene, "present": True, "mean_log2_rpkm_dmso": fnum(dmso.mean())}
        for treatment in ["THZ531_100nM", "THZ531_400nM"]:
            vals = log_expr.loc[gene, meta.index[meta["treatment"] == treatment]]
            row[f"mean_log2_rpkm_{treatment}"] = fnum(vals.mean())
            row[f"delta_{treatment}_minus_dmso"] = fnum(vals.mean() - dmso.mean())
        gene_rows.append(row)
    gene_delta = pd.DataFrame(gene_rows)
    gene_delta.to_csv(OUT / "gse137697_panel_gene_delta.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    order = ["DMSO", "THZ531_100nM", "THZ531_400nM"]
    for score in ["RAI_8_score", "THYROID_NONOVERLAP_score", "STAT3_AP1_DNMT_score"]:
        axes[0].plot(order, [scores.loc[scores["treatment"] == t, score].mean() for t in order],
                     marker="o", label=score)
    axes[0].axhline(0, color="#999", linewidth=0.7)
    axes[0].set_title("GSE137697 panel scores by THZ531 dose")
    axes[0].tick_params(axis="x", labelrotation=18)
    axes[0].legend(fontsize=7)

    present_delta = gene_delta[gene_delta["present"] == True].copy()  # noqa: E712
    plot_genes = [g for g in RAI_8 + THYROID_NONOVERLAP if g in set(present_delta["gene_symbol"])]
    y = np.arange(len(plot_genes))
    d400 = present_delta.set_index("gene_symbol").loc[plot_genes, "delta_THZ531_400nM_minus_dmso"]
    axes[1].barh(y, d400, color=np.where(d400 >= 0, "#506f59", "#8b2635"))
    axes[1].axvline(0, color="#999", linewidth=0.7)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(plot_genes, fontsize=7)
    axes[1].set_title("400 nM THZ531: lineage gene log2(RPKM+1) delta")
    axes[1].set_xlabel("delta vs DMSO")
    fig.tight_layout()
    fig.savefig(FIG / "gse137697_thz531_sidecar.png")
    plt.close(fig)

    dmso = scores[scores["treatment"] == "DMSO"]
    high = scores[scores["treatment"] == "THZ531_400nM"]
    summary = {
        "dataset": "GSE137697",
        "n_samples": int(len(scores)),
        "cell_line": "CAL-62",
        "treatments": order,
        "genes_in_matrix": int(log_expr.shape[0]),
        "rai_found": len(rai_found),
        "rai_missing": rai_missing,
        "thyroid_nonoverlap_found": len(non_found),
        "thyroid_nonoverlap_missing": non_missing,
        "tf_found": len(tf_found),
        "tf_missing": tf_missing,
        "mechanism_found": len(mech_found),
        "mechanism_missing": mech_missing,
        "delta_400nM_rai8_vs_dmso": fnum(high["RAI_8_score"].mean() - dmso["RAI_8_score"].mean()),
        "delta_400nM_dm1_vs_dmso": fnum(high["DM1_like_score"].mean() - dmso["DM1_like_score"].mean()),
        "delta_400nM_mechanism_vs_dmso": fnum(high["STAT3_AP1_DNMT_score"].mean() - dmso["STAT3_AP1_DNMT_score"].mean()),
        "interpretation": "THZ531 perturbation is a small ATC cell-line sidecar, not independent clinical validation.",
    }
    return summary


def analyze_gse232237() -> dict:
    raw_dir = RAW / "GSE232237"
    files = sorted(raw_dir.glob("GSM*.count.tsv.gz"))
    target_genes = sorted(set(RAI_8 + THYROID_NONOVERLAP + TDS_TF + MECHANISM))
    sample_expr = {}
    qc_rows = []
    meta_rows = []

    for path in files:
        name = path.name.replace(".count.tsv.gz", "")
        sample_label = name.split("_", 1)[1]
        histology = "ATC" if sample_label.startswith("AT") else "PTC"
        with gzip.open(path, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            n_cells = len(header) - 1
            lib = np.zeros(n_cells, dtype=float)
            target_counts = {}
            n_genes_streamed = 0
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if not parts:
                    continue
                gene = parts[0]
                vals = np.fromiter((float(x) for x in parts[1:]), dtype=float, count=n_cells)
                lib += vals
                n_genes_streamed += 1
                if gene in target_genes:
                    target_counts[gene] = vals
        lib_safe = lib.copy()
        lib_safe[lib_safe <= 0] = np.nan
        expr_vals = {}
        for gene in target_genes:
            counts = target_counts.get(gene, np.zeros(n_cells, dtype=float))
            norm = np.log1p((counts / lib_safe) * 10000.0)
            expr_vals[gene] = float(np.nanmean(norm))
            qc_rows.append({
                "dataset": "GSE232237",
                "sample_id": sample_label,
                "histology": histology,
                "gene_symbol": gene,
                "present_in_file": gene in target_counts,
                "n_cells": int(n_cells),
                "mean_raw_count": fnum(np.nanmean(counts)),
                "pct_cells_detected": fnum(np.mean(counts > 0)),
                "mean_log1p_cpm10k": fnum(expr_vals[gene]),
            })
        sample_expr[sample_label] = expr_vals
        meta_rows.append({
            "sample_id": sample_label,
            "histology": histology,
            "source_file": path.name,
            "n_cells": int(n_cells),
            "n_genes_streamed": int(n_genes_streamed),
            "median_umi_per_cell": fnum(np.nanmedian(lib)),
            "mean_umi_per_cell": fnum(np.nanmean(lib)),
        })

    expr = pd.DataFrame(sample_expr).loc[target_genes]
    expr.index.name = "gene_symbol"
    meta = pd.DataFrame(meta_rows).set_index("sample_id").loc[expr.columns]
    pd.DataFrame(qc_rows).to_csv(OUT / "gse232237_target_gene_qc.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "gse232237_sample_metadata.tsv", sep="\t")
    expr.to_csv(OUT / "gse232237_target_gene_log1p_cpm10k.tsv", sep="\t")

    z = zscore_within_dataset(expr)
    rai, rai_found, rai_missing = panel_score(z, RAI_8)
    non, non_found, non_missing = panel_score(z, THYROID_NONOVERLAP)
    tf, tf_found, tf_missing = panel_score(z, TDS_TF)
    mech, mech_found, mech_missing = panel_score(z, MECHANISM)
    scores = pd.DataFrame({
        "sample_id": z.columns,
        "histology": meta.loc[z.columns, "histology"].values,
        "n_cells": meta.loc[z.columns, "n_cells"].values,
        "median_umi_per_cell": meta.loc[z.columns, "median_umi_per_cell"].values,
        "RAI_8_score": rai.loc[z.columns].values,
        "DM1_like_score": -rai.loc[z.columns].values,
        "THYROID_NONOVERLAP_score": non.loc[z.columns].values,
        "TF_collapse_score": tf.loc[z.columns].values,
        "STAT3_AP1_DNMT_score": mech.loc[z.columns].values,
    })
    scores.to_csv(OUT / "gse232237_sample_scores.tsv", sep="\t", index=False)

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
        d = cohen_d(atc.to_numpy(dtype=float), ptc.to_numpy(dtype=float))
        contrast_rows.append({
            "dataset": "GSE232237",
            "contrast": "ATC_vs_PTC_pseudobulk_all_cells",
            "score": score,
            "n_atc": int(atc.notna().sum()),
            "n_ptc": int(ptc.notna().sum()),
            "mean_atc": fnum(atc.mean()),
            "mean_ptc": fnum(ptc.mean()),
            "cohen_d_atc_minus_ptc": fnum(d),
            "p_mwu": fnum(p_mwu(atc, ptc)),
            "expected_sign": sign,
            "direction_match": bool(np.sign(d) == sign) if not pd.isna(d) else False,
            "note": "Pseudo-bulk across all cells; no malignant epithelial annotation used.",
        })
    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(OUT / "gse232237_score_contrasts.tsv", sep="\t", index=False)

    gene_delta_rows = []
    for gene in target_genes:
        atc = expr.loc[gene, meta.index[meta["histology"] == "ATC"]]
        ptc = expr.loc[gene, meta.index[meta["histology"] == "PTC"]]
        gene_delta_rows.append({
            "gene_symbol": gene,
            "mean_log1p_cpm10k_atc": fnum(atc.mean()),
            "mean_log1p_cpm10k_ptc": fnum(ptc.mean()),
            "delta_atc_minus_ptc": fnum(atc.mean() - ptc.mean()),
            "p_mwu": fnum(p_mwu(atc, ptc)),
        })
    gene_delta = pd.DataFrame(gene_delta_rows)
    gene_delta.to_csv(OUT / "gse232237_panel_gene_delta.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    score_order = ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TF_collapse_score", "STAT3_AP1_DNMT_score"]
    x = np.arange(len(score_order))
    atc_means = [scores.loc[scores["histology"] == "ATC", s].mean() for s in score_order]
    ptc_means = [scores.loc[scores["histology"] == "PTC", s].mean() for s in score_order]
    axes[0].bar(x - 0.18, ptc_means, width=0.36, label="PTC", color="#2e5b87")
    axes[0].bar(x + 0.18, atc_means, width=0.36, label="ATC", color="#8b2635")
    axes[0].axhline(0, color="#999", linewidth=0.7)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([s.replace("_score", "").replace("_", "\n") for s in score_order], fontsize=7)
    axes[0].set_title("GSE232237 pseudo-bulk score means")
    axes[0].legend(fontsize=8)

    plot_genes = [g for g in RAI_8 + THYROID_NONOVERLAP if g in set(gene_delta["gene_symbol"])]
    gd = gene_delta.set_index("gene_symbol").loc[plot_genes, "delta_atc_minus_ptc"]
    y = np.arange(len(plot_genes))
    axes[1].barh(y, gd, color=np.where(gd >= 0, "#506f59", "#8b2635"))
    axes[1].axvline(0, color="#999", linewidth=0.7)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(plot_genes, fontsize=7)
    axes[1].set_title("ATC - PTC target-gene pseudo-bulk delta")
    axes[1].set_xlabel("delta mean log1p(CPM10k)")
    fig.tight_layout()
    fig.savefig(FIG / "gse232237_scrna_pseudobulk.png")
    plt.close(fig)

    summary = {
        "dataset": "GSE232237",
        "n_samples": int(len(scores)),
        "n_ptc": int((scores["histology"] == "PTC").sum()),
        "n_atc": int((scores["histology"] == "ATC").sum()),
        "n_cells_total": int(meta["n_cells"].sum()),
        "analysis_level": "sample pseudo-bulk across all cells",
        "rai_found": len(rai_found),
        "rai_missing": rai_missing,
        "thyroid_nonoverlap_found": len(non_found),
        "thyroid_nonoverlap_missing": non_missing,
        "tf_found": len(tf_found),
        "tf_missing": tf_missing,
        "mechanism_found": len(mech_found),
        "mechanism_missing": mech_missing,
        "direction_matches": int(contrasts["direction_match"].sum()),
        "direction_cells": int(len(contrasts)),
        "dm1_cohen_d_atc_minus_ptc": fnum(contrasts.loc[contrasts["score"] == "DM1_like_score", "cohen_d_atc_minus_ptc"].iloc[0]),
        "rai8_cohen_d_atc_minus_ptc": fnum(contrasts.loc[contrasts["score"] == "RAI_8_score", "cohen_d_atc_minus_ptc"].iloc[0]),
        "note": "Useful Korean scRNA-derived sample-level check, but not malignant-cell-resolved without cell annotation.",
    }
    return summary


def main() -> None:
    ensure_dirs()
    summaries = {
        "gse146003": analyze_gse146003(),
        "gse60542": analyze_gse60542(),
        "gse137697": analyze_gse137697(),
        "gse232237": analyze_gse232237(),
        "paperclip_status": {
            "installed": True,
            "auth": False,
            "note": "paperclip binary is installed, but search requires login or PAPERCLIP_API_KEY on this server.",
        },
    }
    (OUT / "new_data_pull_summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
