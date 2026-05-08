#!/usr/bin/env python3
"""TCGA-THCA DM1/DM2 survival forest — Paper 1 paper-blocking.

Joins:
  - DM call + clinical (BRAF/RAS/TERT) from d4p2/tcga_with_clinical_mutations.tsv
  - OS/DSS/PFI from /data/thca/repo_data/raw/TCGA_pancan/survival.tsv

Cox regression (lifelines) for univariate + multivariate models.
Outputs:
  - project/results/dm1_robustness_v2026_05_08/survival_forest.{png,pdf}
  - project/results/dm1_robustness_v2026_05_08/survival_cox_hr.tsv
  - project/results/dm1_robustness_v2026_05_08/survival_summary.json

NO HLA × cancer-outcome joining (boundary). DM call is RNA-derived, not HLA-allele.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08"
OUT.mkdir(parents=True, exist_ok=True)

CLIN = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_with_clinical_mutations.tsv"
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")


def short_id(s: str) -> str:
    """TCGA-XX-XXXX-01A -> TCGA-XX-XXXX-01"""
    return "-".join(s.split("-")[:4])[:15] if pd.notna(s) else s


def main() -> None:
    clin = pd.read_csv(CLIN, sep="\t", index_col=0)
    surv = pd.read_csv(SURV, sep="\t")

    clin["sample_short"] = clin.index.to_series().apply(short_id)
    surv["sample_short"] = surv["sample"].apply(short_id)

    merged = clin.merge(surv[["sample_short", "OS", "OS.time", "DSS", "DSS.time", "PFI", "PFI.time"]],
                        on="sample_short", how="left")
    merged = merged.dropna(subset=["DM", "OS.time"])

    # Filter to TCGA-THCA only (samples beginning with TCGA - they all should be)
    print(f"merged sample count: {len(merged)} (DM1: {(merged['DM']=='DM1').sum()}, DM2: {(merged['DM']=='DM2').sum()})")

    # Univariate Cox per outcome
    cox_rows = []
    for outcome in ("OS", "DSS", "PFI"):
        time_col = f"{outcome}.time"
        evt = merged[outcome].astype(float)
        T = merged[time_col].astype(float)
        df = pd.DataFrame({
            "DM1": (merged["DM"] == "DM1").astype(int),
            "age": pd.to_numeric(merged.get("age_at_diagnosis", np.nan), errors="coerce"),
            "braf": pd.to_numeric(merged.get("has_braf_v600e", np.nan), errors="coerce"),
            "ras": pd.to_numeric(merged.get("has_ras_mut", np.nan), errors="coerce"),
            "hashi_otsu": pd.to_numeric(merged.get("hashi_otsu", np.nan), errors="coerce"),
            "T": T, "E": evt,
        }).dropna(subset=["T", "E"])

        # Univariate
        for cov in ("DM1", "age", "braf", "ras", "hashi_otsu"):
            sub = df[[cov, "T", "E"]].dropna()
            if sub[cov].nunique() < 2 or len(sub) < 20:
                continue
            try:
                cph = CoxPHFitter().fit(sub, duration_col="T", event_col="E")
                row = cph.summary.iloc[0]
                cox_rows.append({
                    "outcome": outcome,
                    "covariate": cov,
                    "model": "univariate",
                    "n": int(len(sub)),
                    "n_events": int(sub["E"].sum()),
                    "HR": float(row["exp(coef)"]),
                    "HR_lo95": float(row["exp(coef) lower 95%"]),
                    "HR_hi95": float(row["exp(coef) upper 95%"]),
                    "p": float(row["p"]),
                })
            except Exception as e:
                cox_rows.append({"outcome": outcome, "covariate": cov, "model": "univariate", "n": int(len(sub)), "HR": np.nan, "p": np.nan, "note": str(e)[:60]})

        # Multivariate: DM1 + age + BRAF + RAS
        sub = df[["DM1", "age", "braf", "ras", "T", "E"]].dropna()
        if len(sub) > 30 and sub["DM1"].nunique() == 2:
            try:
                cph = CoxPHFitter().fit(sub, duration_col="T", event_col="E")
                for cov in ("DM1", "age", "braf", "ras"):
                    if cov in cph.summary.index:
                        row = cph.summary.loc[cov]
                        cox_rows.append({
                            "outcome": outcome,
                            "covariate": cov,
                            "model": "multivariate (DM1+age+BRAF+RAS)",
                            "n": int(len(sub)),
                            "n_events": int(sub["E"].sum()),
                            "HR": float(row["exp(coef)"]),
                            "HR_lo95": float(row["exp(coef) lower 95%"]),
                            "HR_hi95": float(row["exp(coef) upper 95%"]),
                            "p": float(row["p"]),
                        })
            except Exception as e:
                cox_rows.append({"outcome": outcome, "covariate": "all", "model": "multivariate", "n": int(len(sub)), "HR": np.nan, "p": np.nan, "note": str(e)[:60]})

    cox = pd.DataFrame(cox_rows)
    cox.to_csv(OUT / "survival_cox_hr.tsv", sep="\t", index=False)

    # ----- Forest plot -----
    plot_df = cox.dropna(subset=["HR"]).copy()
    plot_df = plot_df[(plot_df["HR_lo95"] > 0) & (plot_df["HR_hi95"] < 100)]  # filter wild estimates
    fig, ax = plt.subplots(figsize=(11, max(4, 0.4 * len(plot_df) + 1.5)))
    y = np.arange(len(plot_df))
    colors = []
    for _, r in plot_df.iterrows():
        if r["p"] < 0.05 and r["HR"] > 1:
            colors.append("#d62728")
        elif r["p"] < 0.05 and r["HR"] < 1:
            colors.append("#1f77b4")
        else:
            colors.append("#7f7f7f")
    for i, (_, r) in enumerate(plot_df.iterrows()):
        ax.errorbar(r["HR"], i,
                    xerr=[[max(r["HR"] - r["HR_lo95"], 0)], [max(r["HR_hi95"] - r["HR"], 0)]],
                    fmt="s", color=colors[i], ecolor=colors[i], capsize=4, markersize=6)
    ax.axvline(1.0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yticks(y)
    labels = []
    for _, r in plot_df.iterrows():
        p = r["p"]
        p_str = f"p<{10**(-int(-np.log10(p))):.0e}" if (p < 1e-3) else f"p={p:.3g}"
        labels.append(f"{r['outcome']} | {r['covariate']} ({r['model'][:13]})  n={r['n']}, ev={r.get('n_events', '?')}, {p_str}")
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Hazard ratio (log scale)\nleft = lower risk in covariate=1, right = higher risk", fontsize=9)
    ax.set_title(
        "TCGA-THCA Cox HR forest — DM1 vs DM2 + age + BRAF + RAS + HT-overlap (2026-05-08)\n"
        "Red = HR>1, p<0.05; Blue = HR<1, p<0.05; Grey = ns", fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "survival_forest.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "survival_forest.pdf", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "generated_at": "2026-05-08",
        "n_merged_with_DM_and_survival": int(len(merged)),
        "n_DM1": int((merged["DM"] == "DM1").sum()),
        "n_DM2": int((merged["DM"] == "DM2").sum()),
        "n_OS_events": int(merged["OS"].sum()),
        "n_DSS_events": int(merged["DSS"].sum()),
        "n_PFI_events": int(merged["PFI"].sum()),
        "n_cox_rows": int(len(cox)),
    }
    (OUT / "survival_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"wrote survival_forest.{{png,pdf}}, survival_cox_hr.tsv, survival_summary.json")
    print(json.dumps(summary, indent=2))
    print()
    print("Significant covariates (p<0.05):")
    sig = cox[cox["p"] < 0.05].sort_values("p")
    if len(sig):
        print(sig[["outcome", "covariate", "model", "n", "HR", "HR_lo95", "HR_hi95", "p"]].to_string(index=False))


if __name__ == "__main__":
    main()
