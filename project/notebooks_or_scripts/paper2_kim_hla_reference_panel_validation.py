#!/usr/bin/env python3
"""Paper 2 HLA validation using Kim 2014 Korean HLA reference panel.

This is a metric-matched carrier-vs-carrier validation candidate:
- Cases: Korean PTC carrier counts from the existing Paper 2 table.
- Controls: Kim et al. Korean HLA reference panel, 413 unrelated Koreans,
  parsed from phased SNP2HLA BGL markers.

The Kim panel is a reference-panel control-like baseline, not a clinical
matched healthy cohort. Outputs must keep that caveat.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "project" / "external_refs" / "kim_korean_hla_ref" / "KOR_REF"
P2_OUT = ROOT / "project" / "results" / "p2_pillar1_forest_v2"
ASSET = ROOT / "project" / "papers_hub_2026_05_04" / "assets" / "paper2_hla"
P2_STATS = P2_OUT / "paper2_6allele_or_fisher_primary.tsv"

TARGETS = {
    "A*02:07": ("A", "HLA_A_0207"),
    "B*46:01": ("B", "HLA_B_4601"),
    "C*01:02": ("C", "HLA_C_0102"),
    "DPB1*05:01": ("DPB1", "HLA_DPB1_0501"),
    "DQB1*02:01": ("DQB1", "HLA_DQB1_0201"),
    "DRB1*07:01": ("DRB1", "HLA_DRB1_0701"),
}


def parse_bgl_carriers() -> pd.DataFrame:
    fam = pd.read_csv(REF / "Kim_KOR_HLA.fam", sep=r"\s+", header=None)
    n = len(fam)
    target_markers = {marker: allele for allele, (_, marker) in TARGETS.items()}
    rows = []
    with (REF / "Kim_KOR_HLA.bgl.phased").open() as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 3 or parts[0] != "M" or parts[1] not in target_markers:
                continue
            marker = parts[1]
            values = parts[2:]
            if len(values) != 2 * n:
                raise ValueError(f"{marker} has {len(values)} haplotypes; expected {2*n}")
            allele_copies = sum(v == "P" for v in values)
            carriers = 0
            homozygous = 0
            for i in range(0, len(values), 2):
                dose = int(values[i] == "P") + int(values[i + 1] == "P")
                if dose > 0:
                    carriers += 1
                if dose == 2:
                    homozygous += 1
            allele = target_markers[marker]
            locus, _ = TARGETS[allele]
            rows.append(
                {
                    "allele_4digit": allele,
                    "locus": locus,
                    "kim_marker": marker,
                    "kim_n_individuals": n,
                    "kim_carriers": carriers,
                    "kim_noncarriers": n - carriers,
                    "kim_carrier_frequency": carriers / n,
                    "kim_allele_copies": allele_copies,
                    "kim_denominator_2n": 2 * n,
                    "kim_allele_frequency_from_bgl": allele_copies / (2 * n),
                    "kim_homozygous_count": homozygous,
                }
            )
    return pd.DataFrame(rows)


def add_frq_check(kim: pd.DataFrame) -> pd.DataFrame:
    frq = pd.read_csv(REF / "Kim_KOR_HLA.FRQ.frq", sep=r"\s+")
    frq = frq.rename(columns={"SNP": "kim_marker", "MAF": "kim_frq_maf", "NCHROBS": "kim_frq_nchrobs"})
    return kim.merge(frq[["kim_marker", "kim_frq_maf", "kim_frq_nchrobs"]], on="kim_marker", how="left")


def ci_haldane(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    aa, bb, cc, dd = a, b, c, d
    if min(aa, bb, cc, dd) == 0:
        aa += 0.5
        bb += 0.5
        cc += 0.5
        dd += 0.5
    log_or = np.log((aa * dd) / (bb * cc))
    se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return float(np.exp(log_or)), float(np.exp(log_or - 1.959963984540054 * se)), float(np.exp(log_or + 1.959963984540054 * se))


def validate_against_ptc(kim: pd.DataFrame) -> pd.DataFrame:
    ptc = pd.read_csv(P2_STATS, sep="\t")
    ptc = ptc[["allele_4digit", "ptc_n_individuals", "ptc_carriers", "ptc_noncarriers", "ptc_frequency"]]
    out = kim.merge(ptc, on="allele_4digit", how="left")
    rows = []
    for _, r in out.iterrows():
        a = int(r["ptc_carriers"])
        b = int(r["ptc_noncarriers"])
        c = int(r["kim_carriers"])
        d = int(r["kim_noncarriers"])
        fisher_or, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        or_h, lo, hi = ci_haldane(a, b, c, d)
        rr = r.to_dict()
        rr.update(
            {
                "control_source_id": "KIM2014_KOREAN_HLA_REFERENCE_PANEL_V1",
                "control_citation": "Kim K et al. PLoS ONE 2014; Korean HLA Reference Panel v1.0",
                "control_metric_type": "carrier_frequency_from_phased_reference_panel",
                "ptc_metric_type": "carrier_frequency",
                "metric_compatibility": "carrier_vs_carrier",
                "fisher_exact_or": fisher_or,
                "fisher_exact_p": p,
                "or_haldane_for_ci": or_h,
                "ci95_low": lo,
                "ci95_high": hi,
                "interpretation_status": "metric_matched_reference_panel_validation_candidate",
                "caveat": "Kim panel contains 413 unrelated Koreans and is an imputation reference panel/control-like baseline; phenotype-matched healthy-control status must be stated conservatively.",
            }
        )
        rows.append(rr)
    return pd.DataFrame(rows)


def plot_forest(df: pd.DataFrame) -> None:
    data = df.sort_values("or_haldane_for_ci").reset_index(drop=True)
    y = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    colors = ["#244e73" if x < 1 else "#8f2d25" for x in data["or_haldane_for_ci"]]
    for i, r in data.iterrows():
        ax.plot([r["ci95_low"], r["ci95_high"]], [i, i], color=colors[i], lw=2.8)
        ax.scatter([r["or_haldane_for_ci"]], [i], s=110, color=colors[i], edgecolor="#17212f", zorder=3)
        ax.text(r["ci95_high"] * 1.08, i, f"p={r['fisher_exact_p']:.1e}", va="center", fontsize=9)
    ax.axvline(1, color="#52606f", ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(data["allele_4digit"])
    ax.set_xlabel("PTC carrier vs Kim 2014 Korean reference-panel carrier OR")
    ax.set_title("Paper 2 HLA validation candidate: metric-matched carrier-vs-carrier forest\nControl-like Kim 2014 Korean HLA Reference Panel v1.0, n=413")
    ax.grid(axis="x", color="#eadfcc", lw=0.8)
    fig.tight_layout()
    fig.savefig(ASSET / "F17_kim2014_reference_panel_carrier_validation_forest.png", dpi=220)
    fig.savefig(P2_OUT / "paper2_kim2014_reference_panel_carrier_validation_forest.png", dpi=220)
    plt.close(fig)


def plot_carrier_inputs(df: pd.DataFrame) -> None:
    data = df.set_index("allele_4digit").loc[list(TARGETS)].reset_index()
    x = np.arange(len(data))
    width = 0.38
    fig, ax = plt.subplots(figsize=(11, 5.4))
    ax.bar(x - width / 2, data["ptc_frequency"], width, label="Korean PTC carrier frequency", color="#8f2d25")
    ax.bar(x + width / 2, data["kim_carrier_frequency"], width, label="Kim reference-panel carrier frequency", color="#426b50")
    ax.set_xticks(x)
    ax.set_xticklabels(data["allele_4digit"], rotation=25, ha="right")
    ax.set_ylabel("Carrier frequency")
    ax.set_title("Paper 2 metric-matched inputs\nBoth sides are carrier frequencies over individuals")
    ax.legend(frameon=False)
    for i, r in data.iterrows():
        ax.text(i - width / 2, r["ptc_frequency"] + 0.012, f"{int(r['ptc_carriers'])}/{int(r['ptc_n_individuals'])}", ha="center", fontsize=8)
        ax.text(i + width / 2, r["kim_carrier_frequency"] + 0.012, f"{int(r['kim_carriers'])}/{int(r['kim_n_individuals'])}", ha="center", fontsize=8)
    ax.set_ylim(0, max(data["ptc_frequency"].max(), data["kim_carrier_frequency"].max()) * 1.24)
    fig.tight_layout()
    fig.savefig(ASSET / "F18_kim2014_metric_matched_carrier_inputs.png", dpi=220)
    plt.close(fig)


def write_report(kim: pd.DataFrame, validation: pd.DataFrame) -> None:
    report = ROOT / "project" / "reports" / "2026_05_06_paper2_kim2014_hla_reference_panel_validation.md"
    result_cols = [
        "allele_4digit",
        "ptc_carriers",
        "ptc_n_individuals",
        "kim_carriers",
        "kim_n_individuals",
        "kim_carrier_frequency",
        "kim_allele_frequency_from_bgl",
        "or_haldane_for_ci",
        "ci95_low",
        "ci95_high",
        "fisher_exact_p",
        "interpretation_status",
    ]
    lines = ["| " + " | ".join(result_cols) + " |", "| " + " | ".join(["---"] * len(result_cols)) + " |"]
    for _, r in validation[result_cols].iterrows():
        vals = []
        for c in result_cols:
            v = r[c]
            if isinstance(v, float):
                vals.append(f"{v:.4g}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    report.write_text(
        "# Paper 2 Kim 2014 Korean HLA reference-panel validation\n\n"
        "## Executive verdict\n"
        "- Kim et al. Korean HLA Reference Panel v1.0 is usable as a Korean control-like carrier-frequency baseline candidate for Paper 2.\n"
        "- The panel contains 413 unrelated Korean subjects and covers all six target loci: HLA-A, B, C, DRB1, DPB1, DQB1.\n"
        "- Unlike the previous In/AFND baseline comparison, this validation is metric-matched: PTC carrier frequency vs Kim reference-panel carrier frequency.\n"
        "- Caveat: the panel is an imputation reference panel, not a phenotype-matched healthy-control cohort. Use as validation/sensitivity unless Yu approves primary use.\n\n"
        "## Result table\n\n"
        + "\n".join(lines)
        + "\n\n## Files\n"
        "- `project/results/p2_pillar1_forest_v2/paper2_kim2014_reference_panel_carrier_counts.tsv`\n"
        "- `project/results/p2_pillar1_forest_v2/paper2_kim2014_reference_panel_validation.tsv`\n"
        "- `project/papers_hub_2026_05_04/assets/paper2_hla/F17_kim2014_reference_panel_carrier_validation_forest.png`\n"
        "- `project/papers_hub_2026_05_04/assets/paper2_hla/F18_kim2014_metric_matched_carrier_inputs.png`\n",
        encoding="utf-8",
    )


def main() -> int:
    ASSET.mkdir(parents=True, exist_ok=True)
    P2_OUT.mkdir(parents=True, exist_ok=True)
    kim = add_frq_check(parse_bgl_carriers()).sort_values("allele_4digit")
    validation = validate_against_ptc(kim).sort_values("allele_4digit")
    kim.to_csv(P2_OUT / "paper2_kim2014_reference_panel_carrier_counts.tsv", sep="\t", index=False)
    validation.to_csv(P2_OUT / "paper2_kim2014_reference_panel_validation.tsv", sep="\t", index=False)
    plot_forest(validation)
    plot_carrier_inputs(validation)
    write_report(kim, validation)
    print(validation[["allele_4digit", "ptc_carriers", "ptc_n_individuals", "kim_carriers", "kim_n_individuals", "or_haldane_for_ci", "ci95_low", "ci95_high", "fisher_exact_p"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
