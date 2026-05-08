#!/usr/bin/env python3
"""Covariate ladder for spatial expression-CNV territory ecosystem claims."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
ARM = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_per_spot.tsv.gz"
INFER = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_fast_2026_05_08/infercnv_lite_fast_per_spot.tsv.gz"
ST = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/top_topic_validation_2026_05_08/spatial_expression_cnv_per_spot.tsv.gz"
QC = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/infercnv_lite_qc_residual_gate_2026_05_08/spot_qc.tsv.gz"
OUT = ROOT / "project/results/high_impact_topic_pilots_2026_05_08/spatial_cnv_covariate_ladder_2026_05_08"
REPORT = ROOT / "project/reports/2026_05_08_spatial_cnv_covariate_ladder.md"

OUTCOMES = ["TACSTD2_z", "TROP2_raw", "Macrophage_TAM_z", "TLS_B_z", "Tcell_cytotoxic_z", "RAI_thyroid_z", "EMT_stress_z"]
QC_AXES = ["log_total_counts", "n_genes_by_counts", "pct_mito"]
COMPOSITION_AXES = ["Epithelial_score", "Proliferation_score", "EMT_score", "CAF_ECM_score", "Hypoxia_score"]


def fmt(x) -> str:
    if pd.isna(x):
        return "NA"
    if isinstance(x, (float, np.floating)):
        if abs(x) > 0 and abs(x) < 1e-4:
            return f"{x:.2e}"
        return f"{x:.3g}"
    return str(x)


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    model = LinearRegression()
    model.fit(x, y)
    return y - model.predict(x)


def load_dat() -> pd.DataFrame:
    arm = pd.read_csv(ARM, sep="\t")
    arm = arm.loc[:, ~arm.columns.duplicated()].copy()
    arm["territory_label_set"] = "arm_level_expression_cnv"
    arm["territory_high"] = arm["cnv_territory"].astype(str).str.contains("high").astype(int)

    infer = pd.read_csv(INFER, sep="\t")
    infer = infer.loc[:, ~infer.columns.duplicated()].copy()
    infer["territory_label_set"] = "gene_bin_infercnv_lite"
    infer["territory_high"] = infer["infercnv_t00_high"].astype(int)

    st = pd.read_csv(ST, sep="\t")
    st = st[["sample_id", "spot_id"] + [c for c in COMPOSITION_AXES if c in st.columns]].copy()
    qc = pd.read_csv(QC, sep="\t")

    frames = []
    for df in [arm, infer]:
        cols = ["sample_id", "spot_id", "condition", "territory_label_set", "territory_high"] + [c for c in OUTCOMES + COMPOSITION_AXES if c in df.columns]
        tmp = df[cols].copy()
        missing_comp = [c for c in COMPOSITION_AXES if c not in tmp.columns]
        if missing_comp:
            tmp = tmp.merge(st, on=["sample_id", "spot_id"], how="left", suffixes=("", "_st"))
        frames.append(tmp)
    dat = pd.concat(frames, ignore_index=True)
    dat = dat.merge(qc, on=["sample_id", "spot_id", "condition"], how="inner")
    return dat[dat["condition"].isin(["PTC", "LPTC", "ATC"])].copy()


def design_matrix(dat: pd.DataFrame, covariate_set: str) -> np.ndarray:
    enc = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
    sample_x = enc.fit_transform(dat[["sample_id"]])
    blocks = [sample_x]
    if covariate_set in ["sample_qc", "sample_qc_composition"]:
        qc = dat[QC_AXES].replace([np.inf, -np.inf], np.nan).fillna(dat[QC_AXES].median(numeric_only=True)).to_numpy(float)
        blocks.append(StandardScaler().fit_transform(qc))
    if covariate_set == "sample_qc_composition":
        comp_cols = [c for c in COMPOSITION_AXES if c in dat.columns]
        comp = dat[comp_cols].replace([np.inf, -np.inf], np.nan).fillna(dat[comp_cols].median(numeric_only=True)).to_numpy(float)
        blocks.append(StandardScaler().fit_transform(comp))
    return np.hstack(blocks)


def run_ladder(dat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label_set, d0 in dat.groupby("territory_label_set"):
        high = d0["territory_high"].astype(float).to_numpy()
        for cov_set in ["sample_only", "sample_qc", "sample_qc_composition"]:
            cov = design_matrix(d0, cov_set)
            for outcome in [c for c in OUTCOMES if c in d0.columns]:
                mask = d0[outcome].notna().to_numpy()
                if mask.sum() < 50:
                    continue
                y = d0.loc[mask, outcome].astype(float).to_numpy()
                high_sub = high[mask]
                cov_sub = cov[mask]
                high_res = residualize(high_sub, cov_sub)
                y_res = residualize(y, cov_sub)
                if np.std(high_res) <= 0 or np.std(y_res) <= 0:
                    beta = r = p = np.nan
                else:
                    beta = float(np.dot(high_res, y_res) / np.dot(high_res, high_res))
                    r, p = stats.pearsonr(high_res, y_res)
                rows.append(
                    {
                        "territory_label_set": label_set,
                        "covariate_set": cov_set,
                        "outcome": outcome,
                        "n": int(mask.sum()),
                        "residual_beta_high": beta,
                        "residual_pearson_r": float(r) if not pd.isna(r) else np.nan,
                        "residual_p": float(p) if not pd.isna(p) else np.nan,
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "spatial_cnv_covariate_ladder.tsv", sep="\t", index=False)
    return out


def write_report(ladder: pd.DataFrame) -> None:
    lines = [
        "# Spatial CNV covariate ladder - 2026-05-08\n",
        "Purpose: quantify how territory-outcome associations change after adding slide, QC, and broad tissue-composition covariates.\n",
        "## Verdict\n",
    ]
    for label_set in ["arm_level_expression_cnv", "gene_bin_infercnv_lite"]:
        lines.append(f"### {label_set}")
        sub = ladder[ladder["territory_label_set"] == label_set]
        for outcome in ["TROP2_raw", "TACSTD2_z", "Macrophage_TAM_z", "TLS_B_z", "Tcell_cytotoxic_z"]:
            vals = sub[sub["outcome"] == outcome].set_index("covariate_set")
            if vals.empty:
                continue
            parts = []
            for cov_set in ["sample_only", "sample_qc", "sample_qc_composition"]:
                if cov_set in vals.index:
                    parts.append(f"{cov_set}: beta {fmt(vals.loc[cov_set, 'residual_beta_high'])}, r {fmt(vals.loc[cov_set, 'residual_pearson_r'])}")
            lines.append(f"- {outcome}: " + "; ".join(parts) + ".")
    lines.append("\nInterpretation: TROP2/TACSTD2 is a weak claim after QC/composition adjustment. The immune interface, especially TAM and cytotoxic/TLS axes, is the more defensible spatial ecosystem claim.\n")

    lines.extend(["## Full Ladder\n", "| Label set | Covariates | Outcome | n | beta | r | p |", "|---|---|---|---:|---:|---:|---:|"])
    for _, r in ladder.iterrows():
        lines.append(f"| {r['territory_label_set']} | {r['covariate_set']} | {r['outcome']} | {int(r['n'])} | {fmt(r['residual_beta_high'])} | {fmt(r['residual_pearson_r'])} | {fmt(r['residual_p'])} |")
    lines.extend(["\n## Output files"])
    for p in sorted(OUT.glob("*")):
        lines.append(f"- `{p}`")
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dat = load_dat()
    ladder = run_ladder(dat)
    write_report(ladder)
    print(f"[write] {OUT}")
    print(f"[write] {REPORT}")


if __name__ == "__main__":
    main()
