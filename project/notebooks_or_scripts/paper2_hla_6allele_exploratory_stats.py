#!/usr/bin/env python3
"""Paper 2 HLA 6-allele exploratory OR/Fisher/forest.

This script is intentionally scoped to Paper 2 only.

Analysis rule approved by the latest user instruction:
- Korean PTC values are carrier counts / individuals.
- Korean baseline values are published allele frequencies over 2n alleles.
- We reconstruct baseline allele counts from published AF and 2n.
- We do not convert carrier frequency to allele frequency or vice versa.
- Result is labelled metric-mismatched exploratory statistics.

No Paper 4 GD meta-analysis is executed here.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[2]
P2 = ROOT / "project" / "results" / "p2_pillar1_forest_v2"
V1 = ROOT / "project" / "results" / "p2_pillar1_forest"
WEB_ASSETS = ROOT / "project" / "papers_hub_2026_05_04" / "assets" / "paper2_hla"
REPORTS = ROOT / "project" / "reports"

SOURCE_TABLE = P2 / "paper2_6allele_forest_source_table.tsv"
PTC_TABLE = V1 / "korean_PTC_pool_per_subcohort.tsv"

PRIMARY_OUT = P2 / "paper2_6allele_or_fisher_primary.tsv"
JSON_OUT = P2 / "paper2_6allele_or_fisher_primary.json"
SUMMARY_OUT = REPORTS / "2026_05_06_paper2_6allele_exploratory_stats_report.md"
FOREST_PNG = P2 / "paper2_6allele_exploratory_or_fisher_forest.png"
FOREST_PDF = P2 / "paper2_6allele_exploratory_or_fisher_forest.pdf"


def parse_frequency(value: str) -> float:
    s = str(value).strip()
    if s.endswith("%"):
        return float(s[:-1]) / 100.0
    return float(s)


def odds_ratio_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """Haldane-corrected OR and Wald CI for stable zero-cell reporting.

    Table layout:
      PTC      a carriers, b non-carriers
      Baseline c allele-present, d allele-absent
    """
    aa, bb, cc, dd = map(float, (a, b, c, d))
    if min(aa, bb, cc, dd) == 0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
    log_or = math.log((aa * dd) / (bb * cc))
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return math.exp(log_or), math.exp(log_or - 1.96 * se), math.exp(log_or + 1.96 * se)


def build_primary_stats() -> pd.DataFrame:
    src = pd.read_csv(SOURCE_TABLE, sep="\t")
    ptc = pd.read_csv(PTC_TABLE, sep="\t")
    ptc = ptc[ptc["cohort"].eq("Combined")].copy()
    ptc["carriers"] = pd.to_numeric(ptc["carriers"], errors="coerce")
    ptc["n"] = pd.to_numeric(ptc["n"], errors="coerce")

    rows: list[dict] = []
    for _, r in src.iterrows():
        allele = r["allele_4digit"]
        p = ptc[ptc["allele"].eq(allele)].iloc[0]
        ptc_n = int(p["n"])
        ptc_carriers = int(round(float(p["carriers"])))
        ptc_noncarriers = ptc_n - ptc_carriers

        baseline_n_individuals = int(r["baseline_primary_n"])
        baseline_af = parse_frequency(r["baseline_primary_value"])
        baseline_denominator = 2 * baseline_n_individuals
        baseline_allele_count_float = baseline_af * baseline_denominator
        baseline_allele_count = int(round(baseline_allele_count_float))
        baseline_nonallele_count = baseline_denominator - baseline_allele_count

        table = [[ptc_carriers, ptc_noncarriers],
                 [baseline_allele_count, baseline_nonallele_count]]
        fisher_or, fisher_p = fisher_exact(table, alternative="two-sided")
        or_h, ci_lo, ci_hi = odds_ratio_ci(
            ptc_carriers, ptc_noncarriers,
            baseline_allele_count, baseline_nonallele_count,
        )
        rows.append({
            "allele_4digit": allele,
            "locus": r["locus"],
            "ptc_metric_type": "carrier_frequency",
            "ptc_n_individuals": ptc_n,
            "ptc_carriers": ptc_carriers,
            "ptc_noncarriers": ptc_noncarriers,
            "ptc_frequency": round(ptc_carriers / ptc_n, 6),
            "baseline_source_id": r["baseline_primary_source_id"],
            "baseline_citation": r["baseline_primary_citation"],
            "baseline_metric_type": "allele_frequency",
            "baseline_n_individuals": baseline_n_individuals,
            "baseline_denominator_2n": baseline_denominator,
            "baseline_allele_frequency": round(baseline_af, 6),
            "baseline_allele_count_float": round(baseline_allele_count_float, 3),
            "baseline_allele_count_rounded": baseline_allele_count,
            "baseline_nonallele_count_rounded": baseline_nonallele_count,
            "fisher_exact_or_raw": fisher_or,
            "fisher_exact_p": fisher_p,
            "or_haldane_for_ci": or_h,
            "ci95_low": ci_lo,
            "ci95_high": ci_hi,
            "metric_rule": "PTC carrier count vs baseline published allele count over 2n; metric-mismatched exploratory",
            "interpretation_status": "exploratory_metric_mismatch",
        })
    return pd.DataFrame(rows)


def plot_forest(df: pd.DataFrame) -> None:
    data = df.copy().sort_values("or_haldane_for_ci")
    y = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    colors = ["#8f2d25" if x >= 1 else "#244e73" for x in data["or_haldane_for_ci"]]
    for i, (_, r) in enumerate(data.iterrows()):
        ax.plot([r["ci95_low"], r["ci95_high"]], [i, i], color=colors[i], lw=2.5)
        ax.scatter([r["or_haldane_for_ci"]], [i], s=90, color=colors[i], edgecolor="#17212f", zorder=3)
    ax.axvline(1, color="#52606f", lw=1.2, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(data["allele_4digit"])
    ax.set_xlabel("Exploratory OR (PTC carrier count vs baseline allele count)")
    ax.set_title("Paper 2 HLA 6-allele exploratory Fisher/OR forest\nMetric mismatch explicitly retained; not carrier-to-allele converted")
    ax.grid(axis="x", color="#eadfcc", lw=0.8)
    for i, (_, r) in enumerate(data.iterrows()):
        ax.text(r["ci95_high"] * 1.08, i, f"p={r['fisher_exact_p']:.2g}", va="center", fontsize=9)
    fig.tight_layout()
    FOREST_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FOREST_PNG, dpi=220)
    fig.savefig(FOREST_PDF)
    WEB_ASSETS.mkdir(parents=True, exist_ok=True)
    fig.savefig(WEB_ASSETS / "F13_exploratory_or_fisher_forest.png", dpi=220)
    plt.close(fig)


def plot_stats_table(df: pd.DataFrame) -> None:
    small = df[[
        "allele_4digit", "ptc_carriers", "ptc_n_individuals",
        "baseline_allele_count_rounded", "baseline_denominator_2n",
        "or_haldane_for_ci", "ci95_low", "ci95_high", "fisher_exact_p",
    ]].copy()
    for col in ["or_haldane_for_ci", "ci95_low", "ci95_high"]:
        small[col] = small[col].map(lambda x: f"{x:.3g}")
    small["fisher_exact_p"] = small["fisher_exact_p"].map(lambda x: f"{x:.2g}")
    small["PTC"] = small["ptc_carriers"].astype(str) + "/" + small["ptc_n_individuals"].astype(str)
    small["Baseline"] = small["baseline_allele_count_rounded"].astype(str) + "/" + small["baseline_denominator_2n"].astype(str)
    small = small[["allele_4digit", "PTC", "Baseline", "or_haldane_for_ci", "ci95_low", "ci95_high", "fisher_exact_p"]]
    small.columns = ["Allele", "PTC carriers", "Baseline alleles", "OR", "CI low", "CI high", "Fisher p"]

    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.axis("off")
    tbl = ax.table(cellText=small.values, colLabels=small.columns, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.45)
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#d4c6b3")
        if row == 0:
            cell.set_facecolor("#17212f")
            cell.set_text_props(color="white", weight="bold")
        else:
            cell.set_facecolor("#fffaf1" if row % 2 else "#f6f0e6")
    ax.set_title("Paper 2 6-allele exploratory OR/Fisher table\nPTC carrier count vs baseline allele count; no metric conversion", pad=18)
    fig.tight_layout()
    fig.savefig(WEB_ASSETS / "F14_exploratory_stats_table.png", dpi=220)
    plt.close(fig)


def plot_metric_count_diagram(df: pd.DataFrame) -> None:
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ptc_pct = df["ptc_frequency"] * 100
    base_pct = df["baseline_allele_frequency"] * 100
    width = 0.38
    ax.bar(x - width / 2, ptc_pct, width, label="Korean PTC carrier frequency", color="#8f2d25")
    ax.bar(x + width / 2, base_pct, width, label="Korean baseline allele frequency", color="#426b50")
    ax.set_xticks(x)
    ax.set_xticklabels(df["allele_4digit"], rotation=25, ha="right")
    ax.set_ylabel("Percent")
    ax.set_title("Exploratory values used for OR/Fisher\nDisplayed side-by-side, not metric-converted")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#eadfcc")
    fig.tight_layout()
    fig.savefig(WEB_ASSETS / "F15_exploratory_metric_inputs.png", dpi=220)
    plt.close(fig)


def write_summary(df: pd.DataFrame) -> None:
    rows = []
    for _, r in df.iterrows():
        rows.append(
            f"| {r['allele_4digit']} | {r['ptc_carriers']}/{r['ptc_n_individuals']} | "
            f"{r['baseline_allele_count_rounded']}/{r['baseline_denominator_2n']} | "
            f"{r['or_haldane_for_ci']:.3g} ({r['ci95_low']:.3g}-{r['ci95_high']:.3g}) | "
            f"{r['fisher_exact_p']:.3g} | {r['baseline_source_id']} |"
        )
    text = "\n".join([
        "---",
        'title: "Paper 2 HLA 6-allele exploratory statistics"',
        "date: 2026-05-06",
        'status: "COMPLETE - Paper 2 exploratory OR/Fisher/forest only"',
        "---",
        "",
        "# Executive verdict",
        "",
        "The 6-allele Paper 2 HLA exploratory OR/Fisher/forest has been executed after approval.",
        "The analysis keeps the metric mismatch explicit: Korean PTC uses carrier counts over individuals, while Korean baseline uses published allele frequencies reconstructed as allele counts over 2n.",
        "No carrier-to-allele or allele-to-carrier conversion was performed.",
        "No Paper 4 GD meta-analysis was executed.",
        "",
        "# Primary exploratory table",
        "",
        "| Allele | PTC carriers/n | Baseline alleles/2n | OR, Haldane CI | Fisher exact p | Baseline source |",
        "|---|---:|---:|---:|---:|---|",
        *rows,
        "",
        "# Files created",
        "",
        f"- `{PRIMARY_OUT.relative_to(ROOT)}`",
        f"- `{JSON_OUT.relative_to(ROOT)}`",
        f"- `{FOREST_PNG.relative_to(ROOT)}`",
        f"- `{FOREST_PDF.relative_to(ROOT)}`",
        f"- `{(WEB_ASSETS / 'F13_exploratory_or_fisher_forest.png').relative_to(ROOT)}`",
        f"- `{(WEB_ASSETS / 'F14_exploratory_stats_table.png').relative_to(ROOT)}`",
        f"- `{(WEB_ASSETS / 'F15_exploratory_metric_inputs.png').relative_to(ROOT)}`",
        "",
        "# Caveat",
        "",
        "These are exploratory metric-mismatched statistics. They are useful as a decision/visualization layer, not as a clean population-genetic case-control estimate unless the metric rule is accepted in writing.",
    ])
    SUMMARY_OUT.write_text(text + "\n")


def main() -> int:
    P2.mkdir(parents=True, exist_ok=True)
    WEB_ASSETS.mkdir(parents=True, exist_ok=True)
    df = build_primary_stats()
    df.to_csv(PRIMARY_OUT, sep="\t", index=False)
    df.to_json(JSON_OUT, orient="records", indent=2)
    plot_forest(df)
    plot_stats_table(df)
    plot_metric_count_diagram(df)
    write_summary(df)
    print(df[["allele_4digit", "ptc_carriers", "ptc_n_individuals", "baseline_allele_count_rounded", "baseline_denominator_2n", "or_haldane_for_ci", "ci95_low", "ci95_high", "fisher_exact_p"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
