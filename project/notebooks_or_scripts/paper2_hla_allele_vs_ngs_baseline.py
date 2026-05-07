#!/usr/bin/env python3
"""Paper 2 HLA allele-vs-allele validation against Korean NGS controls.

This script fixes the earlier carrier-vs-allele metric mismatch by counting
PTC allele copies directly from two-column HLA calls per locus, then comparing
those allele counts with Baek 2021 Korean healthy NGS allele counts.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PTC_FILE = ROOT / "project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv"
NGS_SOURCE_FILE = ROOT / "project/results/p2_pillar1_forest_v2/paper2_korean_ngs_control_source_bundle.tsv"
OUT_DIR = ROOT / "project/results/p2_pillar1_forest_v2"
ASSET_DIR = ROOT / "project/papers_hub_2026_05_04/assets/paper2_hla"
REPORT_FILE = ROOT / "project/reports/2026_05_06_paper2_allele_vs_ngs_baseline_validation.md"


TARGETS = [
    ("A*02:07", "A"),
    ("B*46:01", "B"),
    ("C*01:02", "C"),
    ("DPB1*05:01", "DPB1"),
    ("DQB1*02:01", "DQB1"),
    ("DRB1*07:01", "DRB1"),
]


def normalize_4digit(value: object) -> str | None:
    if pd.isna(value):
        return None
    s = str(value).strip()
    if not s or s.upper() in {"NA", "NAN", "NONE"}:
        return None
    m = re.match(r"^([A-Z0-9]+)\*(\d+):(\d+)", s)
    if not m:
        return s
    return f"{m.group(1)}*{m.group(2)}:{m.group(3)}"


def haldane_or_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float, bool]:
    used_haldane = min(a, b, c, d) == 0
    aa, bb, cc, dd = (a, b, c, d)
    if used_haldane:
        aa, bb, cc, dd = (a + 0.5, b + 0.5, c + 0.5, d + 0.5)
    odds_ratio = (aa * dd) / (bb * cc)
    se = math.sqrt((1 / aa) + (1 / bb) + (1 / cc) + (1 / dd))
    lo = math.exp(math.log(odds_ratio) - 1.96 * se)
    hi = math.exp(math.log(odds_ratio) + 1.96 * se)
    return odds_ratio, lo, hi, used_haldane


def fmt_p(p: float) -> str:
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)

    ptc = pd.read_csv(PTC_FILE, sep="\t")
    ngs = pd.read_csv(NGS_SOURCE_FILE, sep="\t")
    baek = ngs[ngs["source_id"] == "BAEK2021_PLOS_NGS"].copy()

    dosage_rows = []
    result_rows = []
    for allele, locus in TARGETS:
        cols = [f"{locus}_a1_4d", f"{locus}_a2_4d"]
        normalized = ptc[cols].apply(lambda col: col.map(normalize_4digit))
        nonmissing = normalized.notna()
        ptc_den = int(nonmissing.to_numpy().sum())
        ptc_count = int((normalized == allele).to_numpy().sum())
        ptc_non = ptc_den - ptc_count
        callable_both = int(nonmissing.all(axis=1).sum())
        callable_any = int(nonmissing.any(axis=1).sum())

        ctrl = baek[baek["allele_4digit"] == allele]
        if ctrl.empty:
            raise RuntimeError(f"Missing Baek 2021 NGS row for {allele}")
        ctrl_row = ctrl.iloc[0]
        ctrl_count = int(ctrl_row["reported_count"])
        ctrl_den = int(ctrl_row["reported_denominator"])
        ctrl_non = ctrl_den - ctrl_count

        fisher_or, fisher_p = fisher_exact([[ptc_count, ptc_non], [ctrl_count, ctrl_non]], alternative="two-sided")
        or_h, lo, hi, used_h = haldane_or_ci(ptc_count, ptc_non, ctrl_count, ctrl_non)

        ptc_freq = ptc_count / ptc_den if ptc_den else float("nan")
        ctrl_freq = ctrl_count / ctrl_den if ctrl_den else float("nan")
        direction = "enriched_in_ptc" if or_h > 1 else "depleted_in_ptc"
        if lo <= 1 <= hi:
            direction = "no_clear_difference"

        dosage_rows.append(
            {
                "allele_4digit": allele,
                "locus": locus,
                "ptc_total_rows": len(ptc),
                "ptc_callable_individuals_both_alleles": callable_both,
                "ptc_callable_individuals_at_least_one_allele": callable_any,
                "ptc_allele_count": ptc_count,
                "ptc_allele_denominator": ptc_den,
                "ptc_allele_frequency": ptc_freq,
                "ptc_non_target_allele_count": ptc_non,
                "source_file": str(PTC_FILE.relative_to(ROOT)),
            }
        )
        result_rows.append(
            {
                "allele_4digit": allele,
                "locus": locus,
                "ptc_metric_type": "allele_frequency",
                "ptc_allele_count": ptc_count,
                "ptc_allele_denominator": ptc_den,
                "ptc_allele_frequency": ptc_freq,
                "control_source_id": "BAEK2021_PLOS_NGS",
                "control_citation": "Baek IC et al. PLoS ONE 2021;16:e0253619",
                "control_n_individuals": 173,
                "control_metric_type": "allele_frequency",
                "control_allele_count": ctrl_count,
                "control_allele_denominator": ctrl_den,
                "control_allele_frequency": ctrl_freq,
                "metric_compatibility": "allele_vs_allele",
                "fisher_exact_or_raw": fisher_or,
                "haldane_or_for_ci": or_h,
                "ci95_low": lo,
                "ci95_high": hi,
                "fisher_p": fisher_p,
                "haldane_correction_used": used_h,
                "direction": direction,
                "interpretation_status": "metric_matched_allele_frequency_validation",
                "caveat": "PTC HLA calls and Baek 2021 NGS healthy-control calls are metric-matched as allele counts, but cohort/platform/control matching remains imperfect.",
            }
        )

    dosage = pd.DataFrame(dosage_rows)
    results = pd.DataFrame(result_rows)
    dosage.to_csv(OUT_DIR / "paper2_ptc_allele_dosage_counts.tsv", sep="\t", index=False)
    results.to_csv(OUT_DIR / "paper2_allele_vs_baek2021_ngs_validation.tsv", sep="\t", index=False)

    plot_forest(results)
    plot_inputs(results)
    plot_flow()
    write_report(results, dosage)


def plot_forest(results: pd.DataFrame) -> None:
    df = results.sort_values("haldane_or_for_ci", ascending=True).reset_index(drop=True)
    colors = ["#2d6f9f" if x < 1 else "#a83c32" for x in df["haldane_or_for_ci"]]
    y = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ax.axvline(1, color="#17212f", lw=1.2, ls="--")
    ax.errorbar(
        df["haldane_or_for_ci"],
        y,
        xerr=[
            df["haldane_or_for_ci"] - df["ci95_low"],
            df["ci95_high"] - df["haldane_or_for_ci"],
        ],
        fmt="none",
        ecolor="#384454",
        elinewidth=1.6,
        capsize=4,
        zorder=1,
    )
    ax.scatter(df["haldane_or_for_ci"], y, s=95, c=colors, edgecolors="#111827", zorder=2)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(df["allele_4digit"], fontsize=11)
    ax.set_xlabel("Allele-vs-allele odds ratio, PTC vs Baek 2021 Korean NGS control (log scale)")
    ax.set_title("F19. Metric-fixed validation: PTC allele dosage vs Baek 2021 NGS allele frequency", loc="left", fontsize=15, weight="bold")
    ax.text(
        0.01,
        1.02,
        "Both sides are allele counts over called alleles. This fixes the previous carrier-vs-allele mismatch, but remains cohort/platform sensitivity evidence.",
        transform=ax.transAxes,
        fontsize=10,
        color="#52606f",
    )
    xmin = min(0.03, df["ci95_low"].min() * 0.7)
    xmax = max(4, df["ci95_high"].max() * 1.25)
    ax.set_xlim(xmin, xmax)
    for i, row in df.iterrows():
        label = f"OR {row['haldane_or_for_ci']:.2g}, p={fmt_p(row['fisher_p'])}"
        ax.text(row["ci95_high"] * 1.08, i, label, va="center", fontsize=9)
    ax.grid(axis="x", color="#e5dfd4", lw=0.8)
    fig.tight_layout()
    fig.savefig(ASSET_DIR / "F19_allele_vs_baek2021_ngs_forest.png", dpi=220)
    plt.close(fig)


def plot_inputs(results: pd.DataFrame) -> None:
    df = results.copy()
    labels = df["allele_4digit"].tolist()
    x = np.arange(len(df))
    w = 0.38
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.bar(x - w / 2, df["ptc_allele_frequency"] * 100, width=w, label="Korean PTC allele frequency", color="#a83c32")
    ax.bar(x + w / 2, df["control_allele_frequency"] * 100, width=w, label="Baek 2021 Korean NGS control allele frequency", color="#426b50")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("Allele frequency (%)")
    ax.set_title("F20. Input values after metric fix: allele frequency vs allele frequency", loc="left", fontsize=15, weight="bold")
    ax.text(
        0.01,
        1.02,
        "Unlike F15, both red and green bars use allele-count denominators. No carrier-frequency conversion was used.",
        transform=ax.transAxes,
        fontsize=10,
        color="#52606f",
    )
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", color="#e5dfd4", lw=0.8)
    for idx, row in df.iterrows():
        ax.text(idx - w / 2, row["ptc_allele_frequency"] * 100 + 0.8, f"{row['ptc_allele_count']}/{row['ptc_allele_denominator']}", ha="center", fontsize=8, rotation=90)
        ax.text(idx + w / 2, row["control_allele_frequency"] * 100 + 0.8, f"{row['control_allele_count']}/{row['control_allele_denominator']}", ha="center", fontsize=8, rotation=90)
    fig.tight_layout()
    fig.savefig(ASSET_DIR / "F20_allele_frequency_inputs_baek2021.png", dpi=220)
    plt.close(fig)


def plot_flow() -> None:
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.axis("off")
    steps = [
        ("PTC HLA calls", "Count allele copies from a1/a2 columns\nper callable locus"),
        ("PTC allele denominator", "Use called alleles, not n carriers\nA/B=1748; C=1746; DRB1=1722;\nDQB1=1262; DPB1=1558"),
        ("Baek 2021 NGS control", "Healthy South Koreans n=173\npublished allele counts / 346"),
        ("Metric-fixed table", "allele count vs allele count\nno carrier/allele conversion"),
        ("Validation forest", "Fisher + Haldane CI\nsensitivity evidence, not final claim"),
    ]
    xs = np.linspace(0.08, 0.92, len(steps))
    for i, ((title, body), x) in enumerate(zip(steps, xs), start=1):
        rect = plt.Rectangle((x - 0.085, 0.36), 0.17, 0.36, facecolor="#fffaf2", edgecolor="#17212f", lw=1.3)
        ax.add_patch(rect)
        ax.text(x - 0.07, 0.67, f"{i}", ha="center", va="center", fontsize=10, color="white", bbox=dict(boxstyle="circle", facecolor="#8f2d25", edgecolor="none"))
        ax.text(x, 0.60, title, ha="center", va="center", fontsize=11, weight="bold", color="#102033")
        ax.text(x, 0.48, body, ha="center", va="center", fontsize=9.2, color="#52606f")
        if i < len(steps):
            ax.annotate("", xy=(xs[i] - 0.09, 0.54), xytext=(x + 0.09, 0.54), arrowprops=dict(arrowstyle="->", lw=1.6, color="#244e73"))
    ax.text(0.5, 0.88, "F21. Statistical fix flow: make the PTC side allele-based before testing", ha="center", fontsize=16, weight="bold", color="#17212f")
    ax.text(0.5, 0.18, "Green meaning: the main metric mismatch is fixed. Red caveat: this still needs independent matched Korean controls and platform/QC review before final association language.", ha="center", fontsize=10.5, color="#8f2d25")
    fig.tight_layout()
    fig.savefig(ASSET_DIR / "F21_metric_fix_flow_allele_vs_allele.png", dpi=220)
    plt.close(fig)


def write_report(results: pd.DataFrame, dosage: pd.DataFrame) -> None:
    rows = []
    for _, r in results.iterrows():
        rows.append(
            f"| {r['allele_4digit']} | {r['ptc_allele_count']}/{r['ptc_allele_denominator']} "
            f"({r['ptc_allele_frequency']*100:.2f}%) | {r['control_allele_count']}/{r['control_allele_denominator']} "
            f"({r['control_allele_frequency']*100:.2f}%) | {r['haldane_or_for_ci']:.3g} | "
            f"{r['ci95_low']:.3g}-{r['ci95_high']:.3g} | {fmt_p(r['fisher_p'])} | {r['direction']} |"
        )
    report = f"""# Paper 2 HLA Allele-vs-Baek 2021 NGS Baseline Validation

Date: 2026-05-06

## Executive Verdict

The carrier-vs-allele mismatch was addressed by recounting Korean PTC HLA calls as allele copies from the two allele columns per locus, then comparing those counts with Baek 2021 healthy Korean NGS allele counts. This creates a metric-matched allele-frequency validation table.

This does not make the result final clinical or population-genetic proof. It is still a sensitivity/validation layer because the PTC calls and Baek 2021 controls differ by cohort design, platform, callable-locus completeness, and control matching.

## Primary Results

| Allele | PTC allele count | Baek 2021 NGS control count | OR | 95% CI | Fisher p | Direction |
|---|---:|---:|---:|---:|---:|---|
{chr(10).join(rows)}

## Files Created

- `project/results/p2_pillar1_forest_v2/paper2_ptc_allele_dosage_counts.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_allele_vs_baek2021_ngs_validation.tsv`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F19_allele_vs_baek2021_ngs_forest.png`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F20_allele_frequency_inputs_baek2021.png`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F21_metric_fix_flow_allele_vs_allele.png`

## Method Boundary

- PTC denominator is the number of non-missing allele calls per locus.
- Baek 2021 denominator is 346 alleles from 173 healthy South Koreans.
- No carrier-to-allele conversion was performed.
- No allele-to-carrier conversion was performed.
- Fisher exact tests were run on allele-count tables.
- Haldane correction was used only for CI/OR display when a zero cell occurred.

## Remaining Caveats

- PTC HLA calls are derived from the Korean PTC pool, not from a matched case-control genotyping study.
- Baek 2021 is an open Korean healthy NGS baseline, but it is not matched by age, sex, region, platform, or technical pipeline to the PTC calls.
- DQB1 and DPB1 have lower PTC callable denominators than A/B/C.
- Final claim still needs matched Korean control carrier or allele-count validation and HLA QC review.
"""
    REPORT_FILE.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
