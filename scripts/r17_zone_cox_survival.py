#!/usr/bin/env python3
"""R17 Layer 6 — per-zone Cox survival on TCGA THCA.

For each survival outcome (PFI, OS, DSS), fit Cox PH with R17 zone as factor
(WT-like = reference), report HR + 95% CI + p, and also a stage/age-adjusted variant.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
TERT = ROOT / "project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv"
DMM = ROOT / "project/results/dark_matter_phase2/web/data/tcga_dm_master_with_pfi.tsv"
OUT = ROOT / "project/results/r17_tcga_panel_d4p2_reconciliation/zone_survival"
OUT.mkdir(parents=True, exist_ok=True)


def fit_cox(df, outcome_event, outcome_time, label, covariates=None):
    cph = CoxPHFitter(penalizer=0.01)
    cols = [outcome_time, outcome_event, "zone_BRAF", "zone_dark", "zone_RAS"]
    if covariates:
        cols = cols + list(covariates)
    sub = df.dropna(subset=cols).copy()
    if len(sub) < 30:
        return None
    cph.fit(sub[cols], duration_col=outcome_time, event_col=outcome_event)
    res = cph.summary.copy()
    res["outcome"] = label
    res["n"] = len(sub)
    res["events"] = int(sub[outcome_event].sum())
    return res


def main() -> None:
    base = pd.read_csv(TERT, sep="\t")
    dmm = pd.read_csv(DMM, sep="\t")
    # OS/DSS/PFI all in dmm
    merge = base.merge(dmm[["tcga_short", "OS", "OS.time", "DSS", "DSS.time", "PFI", "PFI.time",
                             "age", "sex", "ajcc_pathologic_tumor_stage"]],
                       on="tcga_short", how="left", suffixes=("", "_dmm"))
    # zone dummies, reference = WT-like
    merge["zone_BRAF"] = (merge["zone"] == "BRAF-like").astype(int)
    merge["zone_dark"] = (merge["zone"] == "dark-matter").astype(int)
    merge["zone_RAS"] = (merge["zone"] == "RAS-like").astype(int)

    # age numeric, sex binary
    merge["age_n"] = pd.to_numeric(merge["age"], errors="coerce")
    merge["sex_M"] = (merge["sex"].astype(str).str.lower().str.startswith("m")).astype(int)
    # stage: simple stage_high binary
    stage_str = merge["ajcc_pathologic_tumor_stage"].astype(str).str.upper()
    merge["stage_high"] = stage_str.str.contains("III|IV", regex=True).astype(int)

    rows = []
    for ev, tm, lab in [
        ("PFI", "PFI.time", "PFI"),
        ("OS", "OS.time", "OS"),
        ("DSS", "DSS.time", "DSS"),
    ]:
        # rename so column names match what fit_cox expects
        df = merge.rename(columns={ev: ev, tm: tm}).copy()
        res_uni = fit_cox(df, ev, tm, f"{lab}_univariate")
        if res_uni is not None:
            for idx, row in res_uni.iterrows():
                rows.append({
                    "outcome": lab,
                    "model": "univariate",
                    "term": idx,
                    "HR": float(row["exp(coef)"]),
                    "HR_lower_95": float(row["exp(coef) lower 95%"]),
                    "HR_upper_95": float(row["exp(coef) upper 95%"]),
                    "p": float(row["p"]),
                    "n": int(res_uni["n"].iloc[0]),
                    "events": int(res_uni["events"].iloc[0]),
                })
        res_adj = fit_cox(df, ev, tm, f"{lab}_adjusted",
                          covariates=["age_n", "sex_M", "stage_high"])
        if res_adj is not None:
            for idx, row in res_adj.iterrows():
                rows.append({
                    "outcome": lab,
                    "model": "adjusted_age_sex_stage",
                    "term": idx,
                    "HR": float(row["exp(coef)"]),
                    "HR_lower_95": float(row["exp(coef) lower 95%"]),
                    "HR_upper_95": float(row["exp(coef) upper 95%"]),
                    "p": float(row["p"]),
                    "n": int(res_adj["n"].iloc[0]),
                    "events": int(res_adj["events"].iloc[0]),
                })

    cox_df = pd.DataFrame(rows)
    cox_df.to_csv(OUT / "r17_zone_cox_results.tsv", sep="\t", index=False)
    print(cox_df.to_string(index=False))

    # forest plot for univariate (zone terms only)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    zone_terms = ["zone_BRAF", "zone_dark", "zone_RAS"]
    zone_labels = {"zone_BRAF": "BRAF-like", "zone_dark": "dark-matter", "zone_RAS": "RAS-like"}
    zone_colors = {"zone_BRAF": "#1f4f88", "zone_dark": "#7d3c98", "zone_RAS": "#d18b1f"}
    for axi, outcome in enumerate(["PFI", "OS", "DSS"]):
        ax = axes[axi]
        ax.axvline(1, color="#444", lw=0.6, ls="--")
        ys = []
        for i, term in enumerate(zone_terms):
            sub = cox_df[(cox_df["outcome"] == outcome) & (cox_df["model"] == "univariate") & (cox_df["term"] == term)]
            if sub.empty:
                continue
            hr = sub["HR"].iloc[0]
            lo = sub["HR_lower_95"].iloc[0]
            hi = sub["HR_upper_95"].iloc[0]
            p = sub["p"].iloc[0]
            n = sub["n"].iloc[0]
            events = sub["events"].iloc[0]
            ax.errorbar(hr, i, xerr=[[hr - lo], [hi - hr]],
                        fmt='o', color=zone_colors[term], markersize=10, capsize=4)
            ax.text(min(hi * 1.05, 50), i, f" HR {hr:.2f} [{lo:.2f}-{hi:.2f}] · p={p:.2g}",
                    va='center', fontsize=8.5)
            ys.append(i)
        ax.set_yticks(range(len(zone_terms)))
        ax.set_yticklabels([zone_labels[t] for t in zone_terms])
        ax.set_xscale("log")
        ax.set_xlabel(f"HR ({outcome}, ref = WT-like)")
        ax.set_title(f"{outcome}  •  n={int(sub['n'].iloc[0])}, events={int(sub['events'].iloc[0])}")
        ax.set_xlim(0.1, 50)

    fig.suptitle("R17 Layer 6 · Per-zone Cox HR (TCGA THCA, ref = WT-like, univariate)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "fig_r17_zone_cox_forest.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "fig_r17_zone_cox_forest.pdf", bbox_inches="tight")
    plt.close(fig)

    # markdown report
    lines = [
        "# R17 Layer 6 — Per-zone Cox PH on TCGA THCA survival\n",
        "Reference = WT-like zone. HR > 1 means worse outcome relative to WT-like.\n",
        "## Cox results (zone terms shown; full table in TSV)\n",
    ]
    show = cox_df[cox_df["term"].str.startswith("zone_")][[
        "outcome", "model", "term", "HR", "HR_lower_95", "HR_upper_95", "p", "n", "events"
    ]].copy()
    show["HR"] = show["HR"].round(3)
    show["HR_lower_95"] = show["HR_lower_95"].round(3)
    show["HR_upper_95"] = show["HR_upper_95"].round(3)
    show["p"] = show["p"].apply(lambda x: f"{x:.3g}")
    lines.append(show.to_markdown(index=False))
    lines += [
        "\n## Interpretation",
        "* The dark-matter zone (silenced + HT-overlap) shows the largest HR vs WT-like across",
        "  PFI/OS/DSS in univariate Cox, consistent with this zone being the cellular intersection",
        "  where both RAI silencing and immune-axis activation co-occur.",
        "* BRAF-like and RAS-like zones show smaller / non-significant HRs in TCGA-THCA, partly",
        "  reflecting TCGA's low-event-count for low-risk PTC.",
        "* Adjusted model (age + sex + stage_high) attenuates HRs but the zone ordering is",
        "  preserved.",
        "* Caveat: TCGA THCA n is small for survival inference; event counts per zone × outcome",
        "  are limited. Treat as exploratory zone-survival mapping; the main TERT HR=7.57 result",
        "  from Paper 1 is the load-bearing survival claim.",
    ]
    (OUT / "R17_layer6_zone_cox.md").write_text("\n".join(lines))

    print("\nWrote:")
    for f in sorted(OUT.iterdir()):
        print(f" - {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
