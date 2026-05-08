#!/usr/bin/env python3
"""Compare QC sensitivity of arm-level and gene-bin expression-CNV territories."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
ARM = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_per_spot.tsv.gz"
INFER = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_per_spot.tsv.gz"
QC = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_qc_residual_gate_2026_05_08/spot_qc.tsv.gz"
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_qc_compare_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_spatial_cnv_territory_qc_compare.md"

QC_AXES = ["log_total_counts", "n_genes_by_counts", "pct_mito"]
ECOSYSTEM_AXES = [
    "TACSTD2_z",
    "TROP2_raw",
    "Macrophage_TAM_z",
    "TLS_B_z",
    "Tcell_cytotoxic_z",
    "RAI_thyroid_z",
    "EMT_stress_z",
    "DM1_like_score",
    "RAI_8_score",
    "n_fib",
]


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


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    model = LinearRegression()
    model.fit(x, y)
    return y - model.predict(x)


def load_label_sets() -> pd.DataFrame:
    arm = pd.read_csv(ARM, sep="\t")
    arm = arm.loc[:, ~arm.columns.duplicated()].copy()
    arm["territory_label_set"] = "arm_level_expression_cnv"
    arm["territory"] = arm["cnv_territory"].astype(str)
    arm["territory_high"] = arm["territory"].str.contains("high").astype(int)
    arm["territory_signature"] = arm["t00_spatial_cnv_signature"]

    infer = pd.read_csv(INFER, sep="\t")
    infer = infer.loc[:, ~infer.columns.duplicated()].copy()
    infer["territory_label_set"] = "gene_bin_infercnv_lite"
    infer["territory"] = infer["infercnv_t00_territory"].astype(str)
    infer["territory_high"] = infer["infercnv_t00_high"].astype(int)
    infer["territory_signature"] = infer["infercnv_t00_signature"]

    keep = ["sample_id", "spot_id", "condition", "territory_label_set", "territory", "territory_high", "territory_signature"]
    axis_keep = [c for c in ECOSYSTEM_AXES if c in arm.columns or c in infer.columns]
    frames = []
    for df in [arm, infer]:
        cols = keep + [c for c in axis_keep if c in df.columns]
        frames.append(df[cols])
    all_labels = pd.concat(frames, ignore_index=True)

    qc = pd.read_csv(QC, sep="\t")
    return all_labels.merge(qc, on=["sample_id", "spot_id", "condition"], how="inner")


def qc_highlow(dat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cancer = dat[dat["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    for label_set, d0 in cancer.groupby("territory_label_set"):
        for scope, sub in [("pooled_cancer_spots", d0)] + list(d0.groupby("sample_id")):
            high = sub[sub["territory_high"] == 1]
            low = sub[sub["territory"].str.contains("low")]
            for axis in QC_AXES:
                if high[axis].dropna().shape[0] < 2 or low[axis].dropna().shape[0] < 2:
                    continue
                rows.append(
                    {
                        "territory_label_set": label_set,
                        "scope": scope,
                        "condition": sub["condition"].iloc[0] if scope != "pooled_cancer_spots" else "PTC_LPTC_ATC",
                        "axis": axis,
                        "n_high": high[axis].dropna().shape[0],
                        "n_low": low[axis].dropna().shape[0],
                        "delta_high_low": high[axis].mean() - low[axis].mean(),
                        "cohen_d": cohen_d(high[axis], low[axis]),
                        "mw_p": float(stats.mannwhitneyu(high[axis].dropna(), low[axis].dropna()).pvalue),
                    }
                )
    df = pd.DataFrame(rows)
    med = []
    for (label_set, axis), sub in df[~df["scope"].eq("pooled_cancer_spots")].groupby(["territory_label_set", "axis"]):
        med.append(
            {
                "territory_label_set": label_set,
                "scope": "slide_median_effect",
                "condition": "PTC_LPTC_ATC",
                "axis": axis,
                "n_high": sub.shape[0],
                "n_low": sub.shape[0],
                "delta_high_low": float(sub["delta_high_low"].median()),
                "cohen_d": float(sub["cohen_d"].median()),
                "mw_p": np.nan,
            }
        )
    df = pd.concat([df, pd.DataFrame(med)], ignore_index=True)
    df.to_csv(OUT / "territory_qc_highlow_compare.tsv", sep="\t", index=False)
    return df


def residual_ecosystem(dat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cancer = dat[dat["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    for label_set, d0 in cancer.groupby("territory_label_set"):
        enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        sample_x = enc.fit_transform(d0[["sample_id"]])
        qc_x = d0[QC_AXES].replace([np.inf, -np.inf], np.nan).fillna(d0[QC_AXES].median(numeric_only=True)).to_numpy(float)
        qc_x = StandardScaler().fit_transform(qc_x)
        cov = np.hstack([sample_x, qc_x])
        high = d0["territory_high"].astype(float).to_numpy()
        for axis in [a for a in ECOSYSTEM_AXES if a in d0.columns]:
            mask = d0[axis].notna().to_numpy()
            if mask.sum() < 50:
                continue
            y = d0.loc[mask, axis].astype(float).to_numpy()
            high_sub = high[mask]
            cov_sub = cov[mask]
            high_res = residualize(high_sub, cov_sub)
            y_res = residualize(y, cov_sub)
            if np.std(high_res) <= 0 or np.std(y_res) <= 0:
                r = p = beta = np.nan
            else:
                r, p = stats.pearsonr(high_res, y_res)
                beta = float(np.dot(high_res, y_res) / np.dot(high_res, high_res))
            rows.append(
                {
                    "territory_label_set": label_set,
                    "axis": axis,
                    "n": int(mask.sum()),
                    "raw_delta_high_low": float(y[high_sub == 1].mean() - y[high_sub == 0].mean()),
                    "raw_cohen_d": cohen_d(y[high_sub == 1], y[high_sub == 0]),
                    "residual_beta_high": beta,
                    "residual_pearson_r": float(r) if not pd.isna(r) else np.nan,
                    "residual_p": float(p) if not pd.isna(p) else np.nan,
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "territory_ecosystem_residualized_compare.tsv", sep="\t", index=False)
    return df


def write_report(dat: pd.DataFrame, qc: pd.DataFrame, resid: pd.DataFrame) -> None:
    lines = [
        "# Spatial CNV territory QC comparison - 2026-05-08\n",
        "Purpose: compare whether the original arm-level expression-CNV territory and the new gene-bin inferCNV-lite territory survive QC/slide confound checks.\n",
        "## Verdict\n",
    ]
    for label_set in ["arm_level_expression_cnv", "gene_bin_infercnv_lite"]:
        q = qc[(qc["territory_label_set"] == label_set) & (qc["scope"] == "slide_median_effect")].set_index("axis")
        r = resid[resid["territory_label_set"] == label_set].set_index("axis")
        lines.append(f"### {label_set}")
        lines.append(f"- Slide-median log-total-count delta: **{fmt(q.loc['log_total_counts', 'delta_high_low'] if 'log_total_counts' in q.index else np.nan)}**.")
        lines.append(f"- Slide-median n-gene delta: **{fmt(q.loc['n_genes_by_counts', 'delta_high_low'] if 'n_genes_by_counts' in q.index else np.nan)}**.")
        for axis in ["TROP2_raw", "TACSTD2_z", "Macrophage_TAM_z", "TLS_B_z", "Tcell_cytotoxic_z"]:
            if axis in r.index:
                lines.append(
                    f"- Residualized {axis}: beta **{fmt(r.loc[axis, 'residual_beta_high'])}**, r **{fmt(r.loc[axis, 'residual_pearson_r'])}**, p **{fmt(r.loc[axis, 'residual_p'])}**."
                )
    lines.append("\nInterpretation: high-CNV expression territories are strongly coupled to library depth/n-gene counts. TROP2/TACSTD2 enrichment is therefore not robust as a direct CNV-territory consequence after QC residualization; TAM/T cell/TLS immune ecosystem signals are more defensible.\n")

    lines.extend(["## QC High-Low Summary\n", "| Label set | Scope | Axis | delta high-low | d | p |", "|---|---|---|---:|---:|---:|"])
    show_qc = qc[qc["scope"].isin(["pooled_cancer_spots", "slide_median_effect"])]
    for _, r in show_qc.iterrows():
        lines.append(f"| {r['territory_label_set']} | {r['scope']} | {r['axis']} | {fmt(r['delta_high_low'])} | {fmt(r['cohen_d'])} | {fmt(r['mw_p'])} |")

    lines.extend(["\n## Residualized Ecosystem\n", "| Label set | Axis | n | raw delta | raw d | residual beta | residual r | p |", "|---|---|---:|---:|---:|---:|---:|---:|"])
    for _, r in resid.iterrows():
        lines.append(
            f"| {r['territory_label_set']} | {r['axis']} | {int(r['n'])} | {fmt(r['raw_delta_high_low'])} | {fmt(r['raw_cohen_d'])} | {fmt(r['residual_beta_high'])} | {fmt(r['residual_pearson_r'])} | {fmt(r['residual_p'])} |"
        )

    lines.extend(["\n## Output files"])
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dat = load_label_sets()
    qc = qc_highlow(dat)
    resid = residual_ecosystem(dat)
    write_report(dat, qc, resid)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
