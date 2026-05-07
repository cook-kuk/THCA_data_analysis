#!/usr/bin/env python3
"""Update Paper 2/Paper 4 HLA professor pages with validation/meta assets.

This script intentionally separates:
- Paper 2 PTC HLA: exploratory carrier-vs-allele results plus validation gate.
- Paper 4 GD HLA: extractable-source screening meta-analysis for Graves disease.

It does not touch Paper 3 or manuscript prose.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[2]
P2_OUT = ROOT / "project" / "results" / "p2_pillar1_forest_v2"
P4_OUT = ROOT / "project" / "results" / "paper4_gd_hla"
P2_ASSET = ROOT / "project" / "papers_hub_2026_05_04" / "assets" / "paper2_hla"
P4_ASSET = ROOT / "project" / "papers_hub_2026_05_04" / "assets" / "paper4_hla"
P2_STATS = P2_OUT / "paper2_6allele_or_fisher_primary.tsv"


TARGETS = ["A*02:07", "B*46:01", "C*01:02", "DPB1*05:01", "DQB1*02:01", "DRB1*07:01"]


@dataclass
class MetaRow:
    source: str
    country: str
    ancestry: str
    allele: str
    or_value: float
    p_value: float | None
    ci_low: float | None
    ci_high: float | None
    n_case: int | None
    n_control: int | None
    se_source: str
    caveat: str


def ensure_dirs() -> None:
    P2_OUT.mkdir(parents=True, exist_ok=True)
    P4_OUT.mkdir(parents=True, exist_ok=True)
    P2_ASSET.mkdir(parents=True, exist_ok=True)
    P4_ASSET.mkdir(parents=True, exist_ok=True)


def p_to_se(or_value: float, p_value: float) -> float:
    """Approximate log(OR) SE from a two-sided p-value."""
    z = norm.isf(p_value / 2.0)
    return abs(np.log(or_value)) / z


def ci_to_se(ci_low: float, ci_high: float) -> float:
    return (np.log(ci_high) - np.log(ci_low)) / (2.0 * 1.959963984540054)


def fixed_and_random_effects(rows: pd.DataFrame) -> dict[str, float]:
    yi = rows["log_or"].to_numpy(float)
    vi = rows["se"].to_numpy(float) ** 2
    wi = 1.0 / vi
    fixed = float(np.sum(wi * yi) / np.sum(wi))
    q = float(np.sum(wi * (yi - fixed) ** 2))
    df = len(rows) - 1
    c = float(np.sum(wi) - np.sum(wi**2) / np.sum(wi))
    tau2 = max(0.0, (q - df) / c) if df > 0 and c > 0 else 0.0
    wr = 1.0 / (vi + tau2)
    pooled = float(np.sum(wr * yi) / np.sum(wr))
    se = float(np.sqrt(1.0 / np.sum(wr)))
    lo = pooled - 1.959963984540054 * se
    hi = pooled + 1.959963984540054 * se
    p = float(2.0 * norm.sf(abs(pooled / se)))
    i2 = max(0.0, ((q - df) / q) * 100.0) if q > 0 and df > 0 else 0.0
    return {
        "k": int(len(rows)),
        "pooled_log_or": pooled,
        "pooled_or": float(np.exp(pooled)),
        "ci_low": float(np.exp(lo)),
        "ci_high": float(np.exp(hi)),
        "p_random": p,
        "tau2": tau2,
        "q": q,
        "i2_pct": i2,
    }


def build_paper2_validation_gate() -> None:
    stats = pd.read_csv(P2_STATS, sep="\t")
    summary = stats[
        [
            "allele_4digit",
            "ptc_carriers",
            "ptc_n_individuals",
            "ptc_metric_type",
            "baseline_allele_count_rounded",
            "baseline_denominator_2n",
            "baseline_metric_type",
            "fisher_exact_p",
            "or_haldane_for_ci",
        ]
    ].copy()
    summary["current_status"] = "exploratory_metric_mismatch"
    summary["validation_needed"] = "matched_control_carrier_frequency_or_individual_level_control_HLA"
    summary["local_independent_control_genotypes_found"] = "no"
    summary["adjacent_normal_hla_calls_note"] = (
        "GSE213647 adjacent-normal arcasHLA calls exist locally, but they are cancer-patient paired normal tissues; "
        "they are not independent healthy population controls."
    )
    summary.to_csv(P2_OUT / "paper2_matched_control_validation_gate.tsv", sep="\t", index=False)

    labels = [
        "Current PTC carrier signal",
        "Published Korean baseline",
        "Metric mismatch",
        "Independent matched carrier controls",
        "Final case-control association",
    ]
    colors = ["#426b50", "#426b50", "#b58534", "#8f2d25", "#8f2d25"]
    notes = [
        "available\nn=874",
        "available\nallele frequency",
        "blocks final claim",
        "not found locally",
        "blocked",
    ]
    fig, ax = plt.subplots(figsize=(11, 4.7))
    ax.axis("off")
    xs = np.linspace(0.08, 0.92, len(labels))
    for i, (x, lab, color, note) in enumerate(zip(xs, labels, colors, notes)):
        ax.scatter([x], [0.58], s=1800, color=color, edgecolor="#17212f", zorder=3)
        ax.text(x, 0.58, str(i + 1), color="white", ha="center", va="center", fontsize=18, weight="bold")
        ax.text(x, 0.26, lab, ha="center", va="center", fontsize=11, weight="bold", wrap=True)
        ax.text(x, 0.10, note, ha="center", va="center", fontsize=10, color="#52606f")
        if i < len(xs) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.055, 0.58), xytext=(x + 0.055, 0.58),
                        arrowprops=dict(arrowstyle="->", lw=2, color="#52606f"))
    ax.text(
        0.5,
        0.88,
        "Paper 2 validation gate: what must happen before the exploratory HLA forest becomes a final association claim",
        ha="center",
        fontsize=15,
        weight="bold",
    )
    ax.text(
        0.5,
        0.78,
        "Adjacent-normal HLA calls can be a technical sanity check, but independent Korean healthy/control carrier data are required for population association.",
        ha="center",
        fontsize=11,
        color="#8f2d25",
    )
    fig.tight_layout()
    fig.savefig(P2_ASSET / "F16_paper2_matched_control_validation_gate.png", dpi=220)
    plt.close(fig)


def screening_rows_from_sources() -> pd.DataFrame:
    chu = pd.read_csv(ROOT / "project" / "results" / "p2_pillar1_forest" / "chu2018_allele_summary.tsv", sep="\t")
    rows: list[MetaRow] = []
    for _, r in chu[chu["allele"].isin(TARGETS)].iterrows():
        rows.append(
            MetaRow(
                source="Chu 2018",
                country="China",
                ancestry="Han Chinese",
                allele=r["allele"],
                or_value=float(r["OR"]),
                p_value=float(r["p"]),
                ci_low=float(r["ci_lo"]),
                ci_high=float(r["ci_hi"]),
                n_case=int(r["gd_n"]),
                n_control=int(r["ctrl_n"]),
                se_source="reported_CI",
                caveat="count-backed anchor from local Chu 2018 summary table",
            )
        )

    rows.extend(
        [
            MetaRow("Shin 2019", "Korea", "Korean children", "B*46:01", 3.96, 0.008, None, None, 71, 142, "corrected_p_approx", "Pc used to approximate SE; small pediatric GD cohort"),
            MetaRow("Shin 2019", "Korea", "Korean children", "C*01:02", 2.51, 0.04, None, None, 71, 142, "corrected_p_approx", "Pc used to approximate SE; small pediatric GD cohort"),
            MetaRow("Shin 2019", "Korea", "Korean children", "DPB1*05:01", 4.60, 0.003, None, None, 71, 142, "corrected_p_approx", "Pc used to approximate SE; small pediatric GD cohort"),
            MetaRow("Chen 2011", "Taiwan", "Taiwan ethnic Chinese", "B*46:01", 1.33, 1.17e-2, None, None, None, None, "corrected_p_approx", "Bonferroni-corrected p from registry; exact counts still need table extraction"),
            MetaRow("Chen 2011", "Taiwan", "Taiwan ethnic Chinese", "DPB1*05:01", 2.34, 2.58e-10, None, None, None, None, "corrected_p_approx", "Bonferroni-corrected p from registry; exact counts still need table extraction"),
        ]
    )
    data = pd.DataFrame([r.__dict__ for r in rows])
    data["log_or"] = np.log(data["or_value"])
    ses = []
    for _, r in data.iterrows():
        if pd.notna(r["ci_low"]) and pd.notna(r["ci_high"]):
            ses.append(ci_to_se(float(r["ci_low"]), float(r["ci_high"])))
        elif pd.notna(r["p_value"]):
            ses.append(p_to_se(float(r["or_value"]), float(r["p_value"])))
        else:
            ses.append(np.nan)
    data["se"] = ses
    data["ci_low_derived"] = np.exp(data["log_or"] - 1.959963984540054 * data["se"])
    data["ci_high_derived"] = np.exp(data["log_or"] + 1.959963984540054 * data["se"])
    return data


def build_screening_meta() -> tuple[pd.DataFrame, pd.DataFrame]:
    source = screening_rows_from_sources()
    source.to_csv(P4_OUT / "paper4_screening_panasian_meta_source_rows.tsv", sep="\t", index=False)
    meta_rows = []
    for allele, g in source.groupby("allele", sort=False):
        usable = g.dropna(subset=["log_or", "se"])
        if len(usable) >= 2:
            out = fixed_and_random_effects(usable)
            out["allele"] = allele
            out["sources"] = "; ".join(usable["source"].tolist())
            out["meta_status"] = "screening_random_effects"
        else:
            r = usable.iloc[0]
            out = {
                "allele": allele,
                "k": 1,
                "pooled_log_or": float(r["log_or"]),
                "pooled_or": float(r["or_value"]),
                "ci_low": float(r["ci_low_derived"]),
                "ci_high": float(r["ci_high_derived"]),
                "p_random": float(r["p_value"]) if pd.notna(r["p_value"]) else np.nan,
                "tau2": np.nan,
                "q": np.nan,
                "i2_pct": np.nan,
                "sources": str(r["source"]),
                "meta_status": "single_source_anchor_only",
            }
        meta_rows.append(out)
    meta = pd.DataFrame(meta_rows)
    meta = meta[["allele", "k", "pooled_or", "ci_low", "ci_high", "p_random", "tau2", "q", "i2_pct", "sources", "meta_status"]]
    meta.to_csv(P4_OUT / "paper4_screening_panasian_meta.tsv", sep="\t", index=False)
    return source, meta


def plot_screening_meta(source: pd.DataFrame, meta: pd.DataFrame) -> None:
    order = ["DRB1*07:01", "DQB1*02:01", "C*01:02", "A*02:07", "B*46:01", "DPB1*05:01"]
    data = meta.set_index("allele").loc[[a for a in order if a in set(meta["allele"])]].reset_index()
    y = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    colors = ["#244e73" if x < 1 else "#8f2d25" for x in data["pooled_or"]]
    for i, r in data.iterrows():
        ax.plot([r["ci_low"], r["ci_high"]], [i, i], color=colors[i], lw=2.8)
        ax.scatter([r["pooled_or"]], [i], s=110, color=colors[i], edgecolor="#17212f", zorder=3)
        label = "screening meta" if r["k"] >= 2 else "single source"
        ax.text(r["ci_high"] * 1.07, i, f"k={int(r['k'])}; {label}", va="center", fontsize=9)
    ax.axvline(1, color="#52606f", ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(data["allele"])
    ax.set_xlabel("GD/control OR, random-effects where k>=2")
    ax.set_title("Paper 4 screening Pan-Asian GD HLA meta-analysis\nExtractable Chu/Shin/Chen rows only; not the final full-source meta-analysis")
    ax.grid(axis="x", color="#eadfcc", lw=0.8)
    fig.tight_layout()
    fig.savefig(P4_ASSET / "P4_F15_screening_panasian_meta_forest.png", dpi=220)
    fig.savefig(P4_OUT / "paper4_screening_panasian_meta_forest.png", dpi=220)
    plt.close(fig)

    readiness = pd.DataFrame(
        {
            "allele": TARGETS,
            "count_backed_chu": [1, 1, 1, 1, 1, 1],
            "korean_shin_extractable": [0, 1, 1, 1, 0, 0],
            "taiwan_chen_extractable": [0, 1, 0, 1, 0, 0],
            "random_effects_now": [0, 1, 1, 1, 0, 0],
            "full_meta_ready": [0, 0, 0, 0, 0, 0],
        }
    )
    readiness.to_csv(P4_OUT / "paper4_meta_readiness_by_allele.tsv", sep="\t", index=False)
    mat = readiness.drop(columns=["allele"]).to_numpy(int)
    fig, ax = plt.subplots(figsize=(9.5, 4.7))
    cmap = plt.matplotlib.colors.ListedColormap(["#fff0ed", "#edf6ed"])
    ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(np.arange(len(readiness)))
    ax.set_yticklabels(readiness["allele"])
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(["Chu count", "Shin row", "Chen row", "screening meta", "full meta"], rotation=25, ha="right")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, "yes" if mat[i, j] else "gap", ha="center", va="center", fontsize=8)
    ax.set_title("Paper 4 meta-analysis readiness by allele\nScreening meta can run for B*46:01, C*01:02, DPB1*05:01; full meta still blocked")
    fig.tight_layout()
    fig.savefig(P4_ASSET / "P4_F16_meta_readiness_by_allele.png", dpi=220)
    plt.close(fig)


def write_report(source: pd.DataFrame, meta: pd.DataFrame) -> None:
    def md_table(df: pd.DataFrame) -> str:
        cols = df.columns.tolist()
        lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
        for _, row in df.iterrows():
            vals = []
            for col in cols:
                val = row[col]
                if isinstance(val, float):
                    vals.append(f"{val:.4g}")
                else:
                    vals.append(str(val))
            lines.append("| " + " | ".join(vals) + " |")
        return "\n".join(lines)

    report = ROOT / "project" / "reports" / "2026_05_06_hla_validation_and_screening_meta_update.md"
    report.write_text(
        "# HLA validation and screening meta update\n\n"
        "## Executive verdict\n"
        "- Paper 2 matched-control validation is the correct next step, but a final population case-control claim still requires independent Korean healthy/control carrier-frequency or individual-level HLA genotype data.\n"
        "- Local GSE213647 adjacent-normal arcasHLA calls can support a technical sanity check only; they are cancer-patient normal tissues, not independent healthy controls.\n"
        "- Paper 4 screening Pan-Asian GD meta-analysis was executed only for extractable Chu/Shin/Chen source rows. It is not the final full-source meta-analysis.\n\n"
        "## Files created\n"
        "- `project/results/p2_pillar1_forest_v2/paper2_matched_control_validation_gate.tsv`\n"
        "- `project/papers_hub_2026_05_04/assets/paper2_hla/F16_paper2_matched_control_validation_gate.png`\n"
        "- `project/results/paper4_gd_hla/paper4_screening_panasian_meta_source_rows.tsv`\n"
        "- `project/results/paper4_gd_hla/paper4_screening_panasian_meta.tsv`\n"
        "- `project/results/paper4_gd_hla/paper4_meta_readiness_by_allele.tsv`\n"
        "- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F15_screening_panasian_meta_forest.png`\n"
        "- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F16_meta_readiness_by_allele.png`\n\n"
        "## Paper 4 screening meta results\n\n"
        + md_table(meta)
        + "\n\n## Caveats\n"
        "- Shin 2019 and Chen 2011 rows use corrected p-values to approximate SE where exact count/CI rows are not locally extracted.\n"
        "- A*02:07, DQB1*02:01, and DRB1*07:01 remain Chu-only anchors in this update.\n"
        "- Liao 2022, Park 2005, Ueda 2014, and older Japanese/Hong Kong/Taiwan/Korean studies still need source table extraction before final Paper 4 meta-analysis.\n",
        encoding="utf-8",
    )


def main() -> int:
    ensure_dirs()
    build_paper2_validation_gate()
    source, meta = build_screening_meta()
    plot_screening_meta(source, meta)
    write_report(source, meta)
    print("Paper 2 validation gate written")
    print(meta.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
