#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RES = ROOT / "project/results/p2_pillar1_forest_v2"
SUB = ROOT / "project/results/p2_pillar1_forest/korean_PTC_pool_per_subcohort.tsv"
ASSET = ROOT / "project/papers_hub_2026_05_04/assets/paper2_hla"
ASSET.mkdir(parents=True, exist_ok=True)


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        running = min(running, p[idx] * n / (i + 1))
        q[idx] = min(running, 1.0)
    return q.tolist()


def or_ci(a: float, b: float, c: float, d: float, correction: bool = False) -> tuple[float, float, float]:
    if correction:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    if a == 0 or b == 0 or c == 0 or d == 0:
        return np.nan, np.nan, np.nan
    odds = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo = math.exp(math.log(odds) - 1.96 * se)
    hi = math.exp(math.log(odds) + 1.96 * se)
    return float(odds), float(lo), float(hi)


def fisher_row(a: int, b: int, c: int, d: int, correction_for_ci: bool = True) -> dict:
    raw_or, p = fisher_exact([[a, b], [c, d]])
    corr = correction_for_ci and (a == 0 or b == 0 or c == 0 or d == 0)
    ci_or, lo, hi = or_ci(a, b, c, d, correction=corr)
    if not corr:
        ci_or = raw_or
    return {
        "or_raw": float(raw_or) if np.isfinite(raw_or) else np.nan,
        "or_for_ci": float(ci_or) if np.isfinite(ci_or) else np.nan,
        "ci95_low": lo,
        "ci95_high": hi,
        "fisher_p": float(p),
        "haldane_ci_used": bool(corr),
    }


def md_table(df: pd.DataFrame) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
        else:
            d[col] = d[col].map(lambda x: "" if pd.isna(x) else str(x))
    return "\n".join(
        [
            "| " + " | ".join(d.columns) + " |",
            "| " + " | ".join(["---"] * len(d.columns)) + " |",
            *["| " + " | ".join(map(str, row)) + " |" for row in d.to_numpy()],
        ]
    )


def main() -> None:
    primary = pd.read_csv(RES / "paper2_6allele_or_fisher_primary.tsv", sep="\t")
    source = pd.read_csv(RES / "paper2_6allele_forest_source_table.tsv", sep="\t")
    sub = pd.read_csv(SUB, sep="\t")
    alleles = primary["allele_4digit"].tolist()

    sensitivity = []
    for _, r in primary.iterrows():
        allele = r["allele_4digit"]
        ptc_n = int(r["ptc_n_individuals"])
        ptc_car = int(r["ptc_carriers"])
        base_n = int(r["baseline_n_individuals"])
        base_2n = int(r["baseline_denominator_2n"])
        af = float(r["baseline_allele_frequency"])
        ptc_freq = ptc_car / ptc_n

        scenarios = []
        scenarios.append(
            (
                "A_current_carrier_vs_baseline_allele",
                ptc_car,
                ptc_n - ptc_car,
                int(r["baseline_allele_count_rounded"]),
                int(r["baseline_nonallele_count_rounded"]),
                "PTC carrier count over individuals vs baseline allele count over 2n",
            )
        )
        exp_carrier = 1 - (1 - af) ** 2
        exp_base_car = int(round(exp_carrier * base_n))
        scenarios.append(
            (
                "B_baseline_AF_to_HWE_expected_carrier",
                ptc_car,
                ptc_n - ptc_car,
                exp_base_car,
                base_n - exp_base_car,
                "Baseline allele frequency converted to expected carrier frequency under HWE",
            )
        )
        implied_af = 1 - math.sqrt(max(0.0, 1 - ptc_freq))
        ptc_alleles = int(round(implied_af * 2 * ptc_n))
        base_alleles = int(round(af * base_2n))
        scenarios.append(
            (
                "C_PTC_carrier_to_HWE_implied_AF",
                ptc_alleles,
                2 * ptc_n - ptc_alleles,
                base_alleles,
                base_2n - base_alleles,
                "PTC carrier frequency converted to implied allele frequency under HWE",
            )
        )

        for scenario, a, b, c, d, note in scenarios:
            vals = fisher_row(a, b, c, d, correction_for_ci=True)
            sensitivity.append(
                {
                    "allele_4digit": allele,
                    "scenario": scenario,
                    "ptc_event_count": a,
                    "ptc_nonevent_count": b,
                    "baseline_event_count": c,
                    "baseline_nonevent_count": d,
                    "or": vals["or_for_ci"],
                    "ci95_low": vals["ci95_low"],
                    "ci95_high": vals["ci95_high"],
                    "fisher_p": vals["fisher_p"],
                    "haldane_ci_used": vals["haldane_ci_used"],
                    "direction": "enrichment" if vals["or_for_ci"] > 1 else "depletion",
                    "direction_vs_current": "",
                    "note": note,
                }
            )

    sens = pd.DataFrame(sensitivity)
    current_dir = sens[sens["scenario"].eq("A_current_carrier_vs_baseline_allele")].set_index("allele_4digit")["direction"].to_dict()
    sens["direction_vs_current"] = sens.apply(
        lambda r: "same" if r["direction"] == current_dir.get(r["allele_4digit"]) else "flips", axis=1
    )
    sens.to_csv(RES / "paper2_hla_metric_sensitivity.tsv", sep="\t", index=False)

    exclude_d = primary[primary["allele_4digit"] != "DPB1*05:01"].copy()
    d_report = {
        "scenario": "D_exclude_DPB1_source_sensitivity",
        "reason": "DPB1*05:01 baseline source is AFND South Korea, while the other five primary baselines are IN2015.",
        "n_alleles_retained": int(len(exclude_d)),
        "retained_alleles": exclude_d["allele_4digit"].tolist(),
        "excluded_allele": "DPB1*05:01",
    }
    zero = primary[primary["allele_4digit"].eq("DQB1*02:01")].iloc[0]
    zero_report = {
        "scenario": "E_zero_cell_stress_DQB1_02_01",
        "ptc_carriers": int(zero["ptc_carriers"]),
        "ptc_n": int(zero["ptc_n_individuals"]),
        "baseline_allele_count": int(zero["baseline_allele_count_rounded"]),
        "baseline_denominator_2n": int(zero["baseline_denominator_2n"]),
        "raw_fisher_or": float(zero["fisher_exact_or_raw"]),
        "exact_fisher_p": float(zero["fisher_exact_p"]),
        "haldane_or_for_ci": float(zero["or_haldane_for_ci"]),
        "ci95_low": float(zero["ci95_low"]),
        "ci95_high": float(zero["ci95_high"]),
        "warning": "Zero PTC carriers make the depletion visually strong but technically fragile; recheck HLA typing/provenance and validate with matched controls.",
    }

    metric_json = {
        "source_provenance": source.to_dict("records"),
        "scenario_D": d_report,
        "scenario_E": zero_report,
        "metric_sensitivity": sens.to_dict("records"),
    }
    (RES / "paper2_hla_metric_sensitivity.json").write_text(json.dumps(metric_json, indent=2))

    stable = sens.pivot_table(index="allele_4digit", columns="scenario", values="direction_vs_current", aggfunc="first")
    report_lines = [
        "# Paper 2 HLA metric sensitivity report",
        "",
        "Primary result is exploratory because Korean PTC values are carrier counts over individuals, while baseline values are published allele frequencies over 2n.",
        "",
        "## Source/provenance audit",
        md_table(
            source[
                [
                    "allele_4digit",
                    "korean_ptc_metric_type",
                    "korean_ptc_value",
                    "baseline_primary_source_id",
                    "baseline_primary_metric_type",
                    "baseline_primary_value",
                    "caveat",
                ]
            ]
        ),
        "",
        "## Metric sensitivity",
        md_table(
            sens[
                [
                    "allele_4digit",
                    "scenario",
                    "or",
                    "ci95_low",
                    "ci95_high",
                    "fisher_p",
                    "haldane_ci_used",
                    "direction",
                    "direction_vs_current",
                ]
            ]
        ),
        "",
        "## DPB1 source sensitivity",
        f"DPB1*05:01 is flagged because its baseline source is {primary.loc[primary['allele_4digit'].eq('DPB1*05:01'), 'baseline_source_id'].iloc[0]}, not IN2015. Excluding DPB1 leaves five IN2015-backed alleles for sensitivity discussion.",
        "",
        "## DQB1 zero-cell stress",
        f"DQB1*02:01 has 0/874 PTC carriers versus {zero_report['baseline_allele_count']}/{zero_report['baseline_denominator_2n']} baseline alleles. Exact Fisher p={zero_report['exact_fisher_p']:.3g}; Haldane OR for CI={zero_report['haldane_or_for_ci']:.3g}. This is a depletion signal with zero-cell caution.",
    ]
    (RES / "paper2_hla_metric_sensitivity_report.md").write_text("\n".join(report_lines))

    # Multiple testing.
    mt = primary[["allele_4digit", "fisher_exact_p", "fisher_exact_or_raw", "or_haldane_for_ci"]].copy()
    mt["bh_q_six_candidates"] = bh(mt["fisher_exact_p"].tolist())
    mt["significant_at_0_05_within_six"] = mt["bh_q_six_candidates"] < 0.05
    mt["broader_fdr_available"] = False
    mt["broader_fdr_note"] = "Full searched HLA allele universe with p-values is not available locally; only six-candidate FDR was computed."
    mt.to_csv(RES / "paper2_hla_multiple_testing.tsv", sep="\t", index=False)
    (RES / "paper2_hla_multiple_testing_report.md").write_text(
        "\n".join(
            [
                "# Paper 2 HLA multiple-testing report",
                "",
                "BH/FDR was applied to the six candidate alleles only.",
                "",
                md_table(mt),
                "",
                "Full searched allele-universe FDR was not computed because a full local universe of all tested HLA alleles with p-values was not available.",
            ]
        )
    )

    # Subcohort heterogeneity.
    sub6 = sub[sub["allele"].isin(alleles) & ~sub["cohort"].eq("Combined")].copy()
    sub6["carriers"] = pd.to_numeric(sub6["carriers"], errors="coerce")
    sub6["n"] = pd.to_numeric(sub6["n"], errors="coerce")
    sub6["freq"] = pd.to_numeric(sub6["freq"], errors="coerce")
    het_rows = []
    for allele, g in sub6.groupby("allele"):
        gg = g.dropna(subset=["carriers", "n"]).copy()
        total_car = gg["carriers"].sum()
        total_n = gg["n"].sum()
        pooled = total_car / total_n if total_n else np.nan
        max_row = gg.loc[gg["freq"].idxmax()] if len(gg) else None
        min_row = gg.loc[gg["freq"].idxmin()] if len(gg) else None
        tab = gg[["carriers", "n"]].copy()
        if len(gg) >= 2:
            contingency = np.array([[int(r.carriers), int(r.n - r.carriers)] for _, r in gg.iterrows()])
            chi2, p, dof, exp = (np.nan, np.nan, np.nan, None)
            try:
                chi2, p, dof, exp = __import__("scipy.stats").stats.chi2_contingency(contingency)
            except Exception:
                pass
        else:
            p = np.nan
        dominant_share = float(gg["carriers"].max() / total_car) if total_car else np.nan
        dominant_cohort = str(gg.loc[gg["carriers"].idxmax(), "cohort"]) if total_car else ""
        flag = []
        if len(gg) and gg["n"].min() < 20:
            flag.append("small_subcohort_present")
        if dominant_share > 0.75 and total_car > 0:
            flag.append("carrier_count_dominated_by_one_subcohort")
        if pd.notna(p) and p < 0.05:
            flag.append("heterogeneity_p_lt_0_05")
        het_rows.append(
            {
                "allele_4digit": allele,
                "n_subcohorts": int(len(gg)),
                "pooled_carriers": int(total_car),
                "pooled_n": int(total_n),
                "pooled_frequency": pooled,
                "min_frequency": float(min_row["freq"]) if min_row is not None else np.nan,
                "min_frequency_cohort": str(min_row["cohort"]) if min_row is not None else "",
                "max_frequency": float(max_row["freq"]) if max_row is not None else np.nan,
                "max_frequency_cohort": str(max_row["cohort"]) if max_row is not None else "",
                "dominant_carrier_share": dominant_share,
                "dominant_cohort": dominant_cohort,
                "chi_square_p": float(p) if pd.notna(p) else np.nan,
                "flags": ";".join(flag) if flag else "no_major_flag",
            }
        )
    het = pd.DataFrame(het_rows)
    het.to_csv(RES / "paper2_hla_subcohort_heterogeneity.tsv", sep="\t", index=False)
    (RES / "paper2_hla_subcohort_heterogeneity_report.md").write_text(
        "\n".join(
            [
                "# Paper 2 HLA subcohort heterogeneity report",
                "",
                "Subcohorts available locally: K2, Lee2024, and GSE286332_PTC. GSE286332_PTC has n=9, so heterogeneity flags should not be over-interpreted.",
                "",
                md_table(het),
                "",
                "The six-allele PTC pool is not dominated by one large subcohort for most carriers, but the n=9 GSE286332_PTC row produces unstable high frequencies for several alleles. Use subcohort plots as transparency, not validation.",
            ]
        )
    )

    # Claim grades.
    grade_rows = []
    qmap = mt.set_index("allele_4digit")["bh_q_six_candidates"].to_dict()
    for _, r in primary.iterrows():
        allele = r["allele_4digit"]
        orv = float(r["or_haldane_for_ci"])
        grade = "needs_matched_control_validation"
        if allele == "DQB1*02:01":
            grade = "depletion_signal_with_zero_cell_caution"
        elif allele == "DPB1*05:01":
            grade = "source_sensitive"
        elif orv >= 2 and qmap[allele] < 0.05:
            grade = "strong_exploratory_signal"
        elif orv > 1 and qmap[allele] < 0.05:
            grade = "moderate_exploratory_signal"
        risks = []
        if allele == "DPB1*05:01":
            risks.append("baseline source differs from IN2015")
        if allele == "DQB1*02:01":
            risks.append("zero PTC carriers")
        risks.append("PTC carrier vs baseline allele metric mismatch")
        grade_rows.append(
            {
                "allele_4digit": allele,
                "claim_grade": grade,
                "six_candidate_bh_q": qmap[allele],
                "metric_mismatch_risk": "high",
                "source_robustness": "source_sensitive" if allele == "DPB1*05:01" else "IN2015_primary",
                "zero_cell_risk": "high" if allele == "DQB1*02:01" else "low",
                "next_validation_priority": "high" if allele in ["A*02:07", "B*46:01", "DPB1*05:01", "DQB1*02:01"] else "medium",
                "claim_boundary": "; ".join(risks),
            }
        )
    grades = pd.DataFrame(grade_rows)
    grades.to_csv(RES / "paper2_hla_claim_grade.tsv", sep="\t", index=False)
    (RES / "paper2_hla_claim_boundary.md").write_text(
        "\n".join(
            [
                "# Paper 2 HLA claim boundary",
                "",
                "Allowed: exploratory prioritization, professor discussion, hypothesis generation, and validation-roadmap design.",
                "",
                "Forbidden: validated association, causal susceptibility, clinical risk prediction, patient selection, or treating carrier frequency and allele frequency as interchangeable.",
                "",
                md_table(grades),
            ]
        )
    )

    make_figures(primary, sens, sub6, grades, het)


def make_figures(primary: pd.DataFrame, sens: pd.DataFrame, sub6: pd.DataFrame, grades: pd.DataFrame, het: pd.DataFrame) -> None:
    # F16 current vs HWE sensitivity forest.
    fig, ax = plt.subplots(figsize=(11.5, 6.5), dpi=180)
    alleles = primary["allele_4digit"].tolist()
    y = np.arange(len(alleles))
    cur = sens[sens["scenario"].eq("A_current_carrier_vs_baseline_allele")].set_index("allele_4digit").loc[alleles]
    hwe = sens[sens["scenario"].eq("B_baseline_AF_to_HWE_expected_carrier")].set_index("allele_4digit").loc[alleles]
    ax.axvline(1, color="#333", lw=1, ls="--")
    for i, allele in enumerate(alleles):
        ax.plot([cur.loc[allele, "ci95_low"], cur.loc[allele, "ci95_high"]], [i + 0.12, i + 0.12], color="#8f2d25", lw=2)
        ax.scatter(cur.loc[allele, "or"], i + 0.12, color="#8f2d25", s=45, label="Current carrier vs allele" if i == 0 else None)
        ax.plot([hwe.loc[allele, "ci95_low"], hwe.loc[allele, "ci95_high"]], [i - 0.12, i - 0.12], color="#244e73", lw=2)
        ax.scatter(hwe.loc[allele, "or"], i - 0.12, color="#244e73", s=45, label="Baseline AF to HWE carrier" if i == 0 else None)
        stable = "stable" if cur.loc[allele, "direction"] == hwe.loc[allele, "direction"] else "flips"
        ax.text(8.5, i, stable, va="center", fontsize=9, color="#426b50" if stable == "stable" else "#8f2d25")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(alleles)
    ax.invert_yaxis()
    ax.set_xlabel("Exploratory odds ratio (log scale)")
    ax.set_title("F16. Metric sensitivity: current OR vs HWE carrier approximation")
    ax.legend(loc="lower right", frameon=False)
    ax.set_xlim(0.005, 12)
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(ASSET / "F16_metric_sensitivity_forest.png")
    plt.close(fig)

    # F17 subcohort frequency heatmap.
    if len(sub6):
        mat = sub6.pivot(index="allele", columns="cohort", values="freq").loc[alleles]
        fig, ax = plt.subplots(figsize=(8.5, 5.2), dpi=180)
        im = ax.imshow(mat.values, cmap="YlOrRd", vmin=0, vmax=max(0.6, np.nanmax(mat.values)))
        ax.set_xticks(np.arange(len(mat.columns)))
        ax.set_xticklabels(mat.columns, rotation=25, ha="right")
        ax.set_yticks(np.arange(len(mat.index)))
        ax.set_yticklabels(mat.index)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                val = mat.iloc[i, j]
                txt = "" if pd.isna(val) else f"{val:.3f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=8, color="#111")
        ax.set_title("F17. Korean PTC carrier frequency by subcohort")
        fig.colorbar(im, ax=ax, label="Carrier frequency")
        fig.tight_layout()
        fig.savefig(ASSET / "F17_subcohort_frequency_heatmap.png")
        plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(8, 4), dpi=180)
        ax.axis("off")
        ax.text(0.5, 0.5, "Subcohort data not available / not analyzable", ha="center", va="center", fontsize=18)
        fig.savefig(ASSET / "F17_subcohort_frequency_heatmap.png")
        plt.close(fig)

    # F18 claim boundary matrix.
    cols = ["signal_strength", "source_robustness", "zero_cell_risk", "metric_mismatch_risk", "next_validation_priority"]
    matrix = []
    labels = []
    for _, g in grades.iterrows():
        labels.append(g["allele_4digit"])
        strength = {
            "strong_exploratory_signal": 3,
            "moderate_exploratory_signal": 2,
            "source_sensitive": 2,
            "depletion_signal_with_zero_cell_caution": 2,
            "needs_matched_control_validation": 1,
        }[g["claim_grade"]]
        source = 1 if g["source_robustness"] == "source_sensitive" else 3
        zero = 1 if g["zero_cell_risk"] == "high" else 3
        metric = 1
        priority = 3 if g["next_validation_priority"] == "high" else 2
        matrix.append([strength, source, zero, metric, priority])
    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=180)
    im = ax.imshow(matrix, cmap="RdYlGn", vmin=1, vmax=3)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(["Signal", "Source", "Zero-cell", "Metric", "Priority"], rotation=0)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    for i in range(len(labels)):
        for j in range(len(cols)):
            ax.text(j, i, str(matrix[i][j]), ha="center", va="center", fontsize=10, weight="bold")
    ax.set_title("F18. Claim boundary matrix (3=strong/low-risk/high-priority; 1=limited/high-risk)")
    fig.tight_layout()
    fig.savefig(ASSET / "F18_claim_boundary_matrix.png")
    plt.close(fig)

    # F19 roadmap.
    items = [
        ("Matched control\ncarrier data", "highest"),
        ("Individual-level\nHLA genotypes", "highest"),
        ("Carrier-vs-carrier\ncomparison", "highest"),
        ("External Korean\nPTC validation", "high"),
        ("AITD/Hashimoto\nstratification", "high"),
        ("LD / haplotype\nanalysis", "medium"),
        ("HLA typing / imputation\nQC details", "high"),
        ("DPB1 source\nharmonization", "high"),
    ]
    colors = {"highest": "#8f2d25", "high": "#b58534", "medium": "#244e73"}
    fig, ax = plt.subplots(figsize=(12, 6), dpi=180)
    ax.axis("off")
    ax.set_title("F19. Missing evidence roadmap before any validated association claim", fontsize=16, weight="bold", pad=18)
    for idx, (text, level) in enumerate(items):
        x = (idx % 4) * 0.245 + 0.02
        y = 0.66 - (idx // 4) * 0.36
        rect = plt.Rectangle((x, y), 0.21, 0.22, facecolor="#fffaf2", edgecolor=colors[level], lw=2)
        ax.add_patch(rect)
        ax.text(x + 0.105, y + 0.13, text, ha="center", va="center", fontsize=11, weight="bold")
        ax.text(x + 0.105, y + 0.04, level, ha="center", va="center", fontsize=9, color=colors[level])
    fig.tight_layout()
    fig.savefig(ASSET / "F19_missing_evidence_roadmap.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
