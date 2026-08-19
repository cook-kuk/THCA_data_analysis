#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal, spearmanr

from pilot_utils import fdr_bh, pilot_root, setup_logging


AXES = [
    "rai_differentiation_score",
    "hla_i_apm_score",
    "immune_visibility_score",
    "cd8_exclusion_proxy",
    "myeloid_caf_barrier_score",
    "drug_delivery_failure_proxy",
    "aggressive_dedifferentiation_score",
    "proliferation_score",
]

LABELS = {
    "rai_differentiation_score": "RAI differentiation",
    "hla_i_apm_score": "HLA-I/APM",
    "immune_visibility_score": "Immune visibility",
    "cd8_exclusion_proxy": "CD8 exclusion",
    "myeloid_caf_barrier_score": "Myeloid/CAF barrier",
    "drug_delivery_failure_proxy": "Drug-delivery failure",
    "aggressive_dedifferentiation_score": "Aggressive dediff.",
    "proliferation_score": "Proliferation",
}


def stage_order(x: object) -> float:
    s = str(x).upper()
    if "STAGE IV" in s:
        return 4
    if "STAGE III" in s:
        return 3
    if "STAGE II" in s:
        return 2
    if "STAGE I" in s:
        return 1
    m = re.search(r"\b(IV|III|II|I)\b", s)
    return {"I": 1, "II": 2, "III": 3, "IV": 4}.get(m.group(1), np.nan) if m else np.nan


def design(df: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    pieces = []
    for term in terms:
        if term == "driver":
            pieces.append(pd.get_dummies(df["driver_group"], prefix="driver", drop_first=True, dtype=float))
        elif term == "stage":
            pieces.append(pd.get_dummies(df["stage_group"], prefix="stage", drop_first=True, dtype=float))
        elif term == "sex":
            pieces.append(pd.get_dummies(df["gender"].fillna("unknown"), prefix="sex", drop_first=True, dtype=float))
        elif term in df.columns:
            pieces.append(pd.to_numeric(df[term], errors="coerce").rename(term).to_frame())
    if pieces:
        x = pd.concat(pieces, axis=1)
    else:
        x = pd.DataFrame(index=df.index)
    x = x.replace([np.inf, -np.inf], np.nan)
    x.insert(0, "intercept", 1.0)
    return x


def fit_r2(y: pd.Series, x: pd.DataFrame) -> tuple[float, pd.Series]:
    frame = pd.concat([y.rename("y"), x], axis=1).dropna()
    if len(frame) < 20:
        return np.nan, pd.Series(np.nan, index=y.index)
    yy = frame["y"].to_numpy(float)
    xx = frame.drop(columns="y").to_numpy(float)
    beta = np.linalg.lstsq(xx, yy, rcond=None)[0]
    pred = xx @ beta
    ss_res = float(((yy - pred) ** 2).sum())
    ss_tot = float(((yy - yy.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    resid = pd.Series(np.nan, index=y.index)
    resid.loc[frame.index] = yy - pred
    return r2, resid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=pilot_root())
    args = parser.parse_args()
    logger = setup_logging("20_tcga_multivariable_confounding_checks")
    out = args.root / "results" / "extra_analyses"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.root / "results" / "tables" / "tcga_thca_patient_vulnerability_scores.tsv", sep="\t")
    df["driver_group"] = df["driver_anchor"].fillna("No_BRAF_RAS_call")
    stage_source = df["stage"] if "stage" in df.columns else df["ajcc_pathologic_tumor_stage"]
    df["stage_order"] = stage_source.map(stage_order)
    df["stage_group"] = np.where(df["stage_order"].ge(3), "III_IV", np.where(df["stage_order"].isin([1, 2]), "I_II", "unknown"))
    if "age_at_diagnosis" not in df.columns:
        df["age_at_diagnosis"] = df.get("age_at_initial_pathologic_diagnosis", np.nan)

    models = {
        "driver_only": ["driver"],
        "driver_stage": ["driver", "stage"],
        "driver_stage_epithelial": ["driver", "stage", "tumor_epithelial_score"],
        "driver_stage_epithelial_age_sex": ["driver", "stage", "tumor_epithelial_score", "age_at_diagnosis", "sex"],
    }
    r2_rows = []
    residual_tests = []
    corr_rows = []
    for axis in AXES:
        y = pd.to_numeric(df[axis], errors="coerce")
        residuals = {}
        for model, terms in models.items():
            r2, resid = fit_r2(y, design(df, terms))
            residuals[model] = resid
            r2_rows.append(
                {
                    "axis": axis,
                    "axis_label": LABELS[axis],
                    "model": model,
                    "r_squared": r2,
                    "unexplained_fraction": 1 - r2 if np.isfinite(r2) else np.nan,
                    "claim_boundary": "Descriptive variance check; omitted fusions, CNV, methylation, purity, and treatment variables may contribute.",
                }
            )
        resid = residuals["driver_stage_epithelial_age_sex"]
        tmp = pd.DataFrame({"label": df["primary_vulnerability_label"], "resid": resid}).dropna()
        groups = [g["resid"].to_numpy() for _, g in tmp.groupby("label") if len(g) >= 2]
        p = kruskal(*groups).pvalue if len(groups) >= 2 else np.nan
        residual_tests.append(
            {
                "axis": axis,
                "axis_label": LABELS[axis],
                "residual_model": "driver_stage_epithelial_age_sex",
                "kruskal_p_residual_by_vulnerability_label": p,
                "n_samples": len(tmp),
                "claim_boundary": "Residual label separation suggests vulnerability labels are not just driver/stage/epithelial score, but this is descriptive.",
            }
        )
        for cov in ["tumor_epithelial_score", "age_at_diagnosis", "stage_order"]:
            frame = pd.DataFrame({"axis": y, "cov": pd.to_numeric(df[cov], errors="coerce")}).dropna()
            rho, p = (np.nan, np.nan)
            if frame["axis"].nunique() > 1 and frame["cov"].nunique() > 1:
                rho, p = spearmanr(frame["axis"], frame["cov"])
            corr_rows.append({"axis": axis, "axis_label": LABELS[axis], "covariate": cov, "spearman_rho": rho, "p_value": p, "n": len(frame)})

    r2 = pd.DataFrame(r2_rows)
    res = pd.DataFrame(residual_tests)
    res["fdr_q"] = fdr_bh(res["kruskal_p_residual_by_vulnerability_label"])
    corr = pd.DataFrame(corr_rows)
    corr["fdr_q"] = fdr_bh(corr["p_value"])
    r2.to_csv(tables / "tcga_axis_variance_models.tsv", sep="\t", index=False)
    res.to_csv(tables / "tcga_residual_vulnerability_label_tests.tsv", sep="\t", index=False)
    corr.to_csv(tables / "tcga_axis_covariate_correlations.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(11, 6))
    sns.barplot(data=r2, x="r_squared", y="axis_label", hue="model", ax=ax)
    ax.set_title("How much do driver/stage/epithelial covariates explain K-Thyro axes?", weight="bold")
    ax.set_xlabel("R-squared")
    ax.set_ylabel("")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(figs / "tcga_axis_variance_models.png", dpi=300, bbox_inches="tight")
    fig.savefig(figs / "tcga_axis_variance_models.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

    pivot = corr.pivot(index="axis_label", columns="covariate", values="spearman_rho")
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    sns.heatmap(pivot, cmap="vlag", center=0, annot=True, fmt=".2f", linewidths=0.4, linecolor="#222", ax=ax)
    ax.set_title("Axis correlations with potential confounders", weight="bold")
    fig.tight_layout()
    fig.savefig(figs / "tcga_axis_covariate_correlation_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(figs / "tcga_axis_covariate_correlation_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Wrote TCGA multivariable/confounding checks.")


if __name__ == "__main__":
    main()
