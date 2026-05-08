#!/usr/bin/env python3
"""QC and residual confound gate for inferCNV-lite territories."""

from __future__ import annotations

import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
INFER = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_per_spot.tsv.gz"
SAMPLE_META = ROOT / "project/data/processed/GSE250521/sample_metadata.tsv"
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_qc_residual_gate_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_infercnv_lite_qc_residual_gate.md"

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
QC_AXES = ["log_total_counts", "n_genes_by_counts", "pct_mito"]


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


def dense_csr(x) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(np.asarray(x))


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


def compute_spot_qc() -> pd.DataFrame:
    out_path = OUT / "spot_qc.tsv.gz"
    if out_path.exists():
        return pd.read_csv(out_path, sep="\t")

    samples = pd.read_csv(SAMPLE_META, sep="\t")
    rows = []
    for row in samples.itertuples(index=False):
        a = ad.read_h5ad(ROOT / row.h5ad)
        x = dense_csr(a.X)
        obs = a.obs.reset_index()
        first_col = obs.columns[0]
        if first_col != "spot_id":
            obs = obs.rename(columns={first_col: "spot_id"})
        total = np.asarray(x.sum(axis=1)).ravel()
        n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
        mito_mask = np.asarray(pd.Index(a.var_names.astype(str)).str.upper().str.startswith("MT-"))
        if mito_mask.any():
            mito = np.asarray(x[:, mito_mask].sum(axis=1)).ravel()
            pct_mito = mito / np.maximum(total, 1) * 100.0
        else:
            pct_mito = np.full(x.shape[0], np.nan)
        rows.append(
            pd.DataFrame(
                {
                    "sample_id": row.sample_id,
                    "spot_id": obs["spot_id"].to_numpy(),
                    "condition": condition_from_sample(row.sample_id),
                    "total_counts": total,
                    "log_total_counts": np.log1p(total),
                    "n_genes_by_counts": n_genes,
                    "pct_mito": pct_mito,
                }
            )
        )
        print(f"[qc] {row.sample_id} spots={x.shape[0]}")
    qc = pd.concat(rows, ignore_index=True)
    qc.to_csv(out_path, sep="\t", index=False, compression="gzip")
    return qc


def high_low_qc(dat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    cancer = dat[dat["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    for scope, sub in [("pooled_cancer_spots", cancer)] + list(cancer.groupby("sample_id")):
        high = sub[sub["infercnv_t00_high"] == 1]
        low = sub[sub["infercnv_t00_territory"].astype(str).str.contains("low")]
        for axis in QC_AXES:
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
    df.to_csv(OUT / "infercnv_lite_qc_highlow.tsv", sep="\t", index=False)
    return df


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    model = LinearRegression()
    model.fit(x, y)
    return y - model.predict(x)


def ecosystem_residual_gate(dat: pd.DataFrame) -> pd.DataFrame:
    cancer = dat[dat["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    sample_x = enc.fit_transform(cancer[["sample_id"]])
    qc_x = cancer[QC_AXES].replace([np.inf, -np.inf], np.nan).fillna(cancer[QC_AXES].median(numeric_only=True)).to_numpy(float)
    qc_x = StandardScaler().fit_transform(qc_x)
    cov = np.hstack([sample_x, qc_x])

    rows = []
    high = cancer["infercnv_t00_high"].astype(float).to_numpy()
    high_res = residualize(high, cov)
    for axis in [a for a in ECOSYSTEM_AXES if a in cancer.columns]:
        submask = cancer[axis].notna().to_numpy()
        if submask.sum() < 50:
            continue
        y = cancer.loc[submask, axis].astype(float).to_numpy()
        cov_sub = cov[submask]
        high_sub = high[submask]
        high_res_sub = residualize(high_sub.astype(float), cov_sub)
        y_res = residualize(y, cov_sub)
        if np.std(high_res_sub) <= 0 or np.std(y_res) <= 0:
            r = p = beta = np.nan
        else:
            r, p = stats.pearsonr(high_res_sub, y_res)
            beta = float(np.dot(high_res_sub, y_res) / np.dot(high_res_sub, high_res_sub))
        raw_high = y[high_sub == 1]
        raw_low = y[high_sub == 0]
        rows.append(
            {
                "axis": axis,
                "n": int(submask.sum()),
                "raw_delta_high_low": float(raw_high.mean() - raw_low.mean()) if raw_high.size and raw_low.size else np.nan,
                "raw_cohen_d": cohen_d(raw_high, raw_low),
                "residual_beta_high": beta,
                "residual_pearson_r": float(r) if not pd.isna(r) else np.nan,
                "residual_p": float(p) if not pd.isna(p) else np.nan,
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "infercnv_lite_ecosystem_residualized_for_slide_qc.tsv", sep="\t", index=False)
    return df


def write_report(dat: pd.DataFrame, qc_highlow: pd.DataFrame, resid: pd.DataFrame) -> None:
    qc_med = qc_highlow[qc_highlow["scope"] == "slide_median_effect"].set_index("axis")
    trop2 = resid[resid["axis"] == "TROP2_raw"].iloc[0] if (resid["axis"] == "TROP2_raw").any() else None
    tac = resid[resid["axis"] == "TACSTD2_z"].iloc[0] if (resid["axis"] == "TACSTD2_z").any() else None
    tam = resid[resid["axis"] == "Macrophage_TAM_z"].iloc[0] if (resid["axis"] == "Macrophage_TAM_z").any() else None

    lines = [
        "# inferCNV-lite QC/residual gate - 2026-05-08\n",
        "Purpose: test whether inferCNV-lite T00-high territories are mostly library-size/quality artifacts, and whether ecosystem enrichment survives slide and QC residualization.\n",
        "## Verdict\n",
    ]
    lines.append(f"- Cancer spots evaluated: **{dat['condition'].isin(['PTC', 'LPTC', 'ATC']).sum()}**.")
    lines.append(f"- Slide-median high-low log-total-count delta: **{fmt(qc_med.loc['log_total_counts', 'delta_high_low'] if 'log_total_counts' in qc_med.index else np.nan)}**.")
    lines.append(f"- Slide-median high-low n-gene delta: **{fmt(qc_med.loc['n_genes_by_counts', 'delta_high_low'] if 'n_genes_by_counts' in qc_med.index else np.nan)}**.")
    lines.append(f"- Slide-median high-low pct-mito delta: **{fmt(qc_med.loc['pct_mito', 'delta_high_low'] if 'pct_mito' in qc_med.index else np.nan)}**.")
    if trop2 is not None:
        lines.append(f"- Residualized TROP2_raw association with T00-high after slide + QC: beta **{fmt(trop2['residual_beta_high'])}**, r **{fmt(trop2['residual_pearson_r'])}**, p **{fmt(trop2['residual_p'])}**.")
    if tac is not None:
        lines.append(f"- Residualized TACSTD2_z association: beta **{fmt(tac['residual_beta_high'])}**, r **{fmt(tac['residual_pearson_r'])}**, p **{fmt(tac['residual_p'])}**.")
    if tam is not None:
        lines.append(f"- Residualized Macrophage_TAM_z association: beta **{fmt(tam['residual_beta_high'])}**, r **{fmt(tam['residual_pearson_r'])}**, p **{fmt(tam['residual_p'])}**.")
    lines.append("- Interpretation: if residualized ecosystem betas remain positive, the territory ecosystem claim is not just a slide/QC artifact. QC deltas should still be disclosed because expression-CNV methods are inherently sensitive to transcriptome quality and cell composition.\n")

    lines.extend(["## QC High-Low Checks\n", "| Scope | Axis | delta high-low | d | p |", "|---|---|---:|---:|---:|"])
    for _, r in qc_highlow[qc_highlow["scope"].isin(["pooled_cancer_spots", "slide_median_effect"])].iterrows():
        lines.append(f"| {r['scope']} | {r['axis']} | {fmt(r['delta_high_low'])} | {fmt(r['cohen_d'])} | {fmt(r['mw_p'])} |")

    lines.extend(["\n## Residualized Ecosystem Associations\n", "| Axis | n | raw delta | raw d | residual beta | residual r | p |", "|---|---:|---:|---:|---:|---:|---:|"])
    for _, r in resid.iterrows():
        lines.append(
            f"| {r['axis']} | {int(r['n'])} | {fmt(r['raw_delta_high_low'])} | {fmt(r['raw_cohen_d'])} | {fmt(r['residual_beta_high'])} | {fmt(r['residual_pearson_r'])} | {fmt(r['residual_p'])} |"
        )

    lines.extend(["\n## Output files"])
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    per = pd.read_csv(INFER, sep="\t")
    per = per.loc[:, ~per.columns.duplicated()].copy()
    qc = compute_spot_qc()
    dat = per.merge(qc, on=["sample_id", "spot_id", "condition"], how="inner")
    qc_highlow = high_low_qc(dat)
    resid = ecosystem_residual_gate(dat)
    write_report(dat, qc_highlow, resid)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
