#!/usr/bin/env python3
"""Build preregistered impact v5 for CROSS-Neo.

v5 converts the v4 value-of-information plate into a preregistered assay
statistics package: endpoints, binomial power/false-unlock risk, interpretation
templates, reviewer-risk safeguards, and a live dashboard.
"""

from __future__ import annotations

import json
import math
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10"
V4_DIR = BASE / "translational_impact_v4_2026_05_10"
OUT_DIR = BASE / "preregistered_impact_v5_2026_05_10"
FIG_DIR = OUT_DIR / "figures"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_preregistered_v5"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_preregistered_v5"

PLATE_V4 = V4_DIR / "wetlab_plate_v4_value_of_information.tsv"
WELLMAP_V4 = V4_DIR / "assay_96well_map_v4.tsv"
LADDER_V4 = V4_DIR / "claim_unlock_ladder_v4.tsv"
RULES_V4 = V4_DIR / "assay_decision_rules_v4.tsv"
SUMMARY_V4 = V4_DIR / "translational_impact_v4_summary.json"


ARM_CONFIG = {
    "A_clean_discovery": {
        "endpoint": "mutant_pMHC_binding_with_WT_decoy_specificity",
        "threshold_direction": "at_least",
        "threshold_successes": 2,
        "confirmatory_threshold_successes": 3,
        "null_rate": 0.10,
        "target_rate": 0.45,
        "claim_layer": "wetlab-supported clean antigen discovery lane",
        "primary_risk": "false discovery from overfit clean score",
    },
    "B_mechanism_TCR_MD": {
        "endpoint": "pMHC_binding_plus_TCR_readout",
        "threshold_direction": "at_least",
        "threshold_successes": 1,
        "confirmatory_threshold_successes": 2,
        "null_rate": 0.05,
        "target_rate": 0.55,
        "claim_layer": "orthogonal TCR/MD mechanistic support",
        "primary_risk": "TCR evidence is source/overlap diagnostic only",
    },
    "C_label_rescue": {
        "endpoint": "provenance_resolved_positive_binding_rescue",
        "threshold_direction": "at_least",
        "threshold_successes": 1,
        "confirmatory_threshold_successes": 2,
        "null_rate": 0.05,
        "target_rate": 0.35,
        "claim_layer": "label-noise or assay-mismatch evidence",
        "primary_risk": "rescue case misread as relabeling without provenance",
    },
    "D_specificity_moat": {
        "endpoint": "hard_negative_remains_weak",
        "threshold_direction": "at_least",
        "threshold_successes": 3,
        "confirmatory_threshold_successes": 4,
        "null_rate": 0.50,
        "target_rate": 0.85,
        "claim_layer": "specificity-aware control moat",
        "primary_risk": "hard negatives bind as often as positives",
    },
    "E_positive_QC_control": {
        "endpoint": "positive_control_assay_dynamic_range",
        "threshold_direction": "at_least",
        "threshold_successes": 3,
        "confirmatory_threshold_successes": 3,
        "null_rate": 0.20,
        "target_rate": 0.90,
        "claim_layer": "assay QC only",
        "primary_risk": "positive controls fail or nonspecific binding dominates",
    },
    "F_model_boundary": {
        "endpoint": "TCR_structure_disagreement_resolved",
        "threshold_direction": "at_least",
        "threshold_successes": 2,
        "confirmatory_threshold_successes": 3,
        "null_rate": 0.20,
        "target_rate": 0.55,
        "claim_layer": "model-boundary clarification",
        "primary_risk": "disagreement remains unresolved",
    },
}


def read_tsv(path: Path, required: bool = False) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def write_tsv(df: pd.DataFrame, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    df.to_csv(path, sep="\t", index=False)
    return path


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="float64")
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def binom_p_ge(n: int, k: int, p: float) -> float:
    if n <= 0:
        return float("nan")
    return float(sum(math.comb(n, i) * (p**i) * ((1 - p) ** (n - i)) for i in range(k, n + 1)))


def binom_p_le(n: int, k: int, p: float) -> float:
    if n <= 0:
        return float("nan")
    return float(sum(math.comb(n, i) * (p**i) * ((1 - p) ** (n - i)) for i in range(0, k + 1)))


def unlock_probability(n: int, threshold: int, p: float, direction: str = "at_least") -> float:
    if direction == "at_most":
        return binom_p_le(n, threshold, p)
    return binom_p_ge(n, threshold, p)


def build_endpoint_plan(plate: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for arm, adf in plate.groupby("plate_v4_arm", sort=True):
        config = ARM_CONFIG.get(arm, {})
        n = int(len(adf))
        threshold = int(config.get("threshold_successes", max(1, round(n * 0.5))))
        confirmatory_threshold = int(config.get("confirmatory_threshold_successes", threshold))
        direction = config.get("threshold_direction", "at_least")
        null_rate = float(config.get("null_rate", 0.10))
        target_rate = float(config.get("target_rate", float(adf["v4_lane_specific_prior"].mean())))
        mean_prior = float(adf["v4_lane_specific_prior"].mean())
        mean_unlock = float(adf["v4_claim_unlock_score"].mean())
        mean_eig = float(adf["v4_expected_information_gain"].mean())
        false_unlock = unlock_probability(n, threshold, null_rate, direction)
        power_target = unlock_probability(n, threshold, target_rate, direction)
        power_prior = unlock_probability(n, threshold, mean_prior, direction)
        confirmatory_false_unlock = unlock_probability(n, confirmatory_threshold, null_rate, direction)
        confirmatory_power_target = unlock_probability(n, confirmatory_threshold, target_rate, direction)
        confirmatory_power_prior = unlock_probability(n, confirmatory_threshold, mean_prior, direction)
        rows.append(
            {
                "plate_v4_arm": arm,
                "impact_lane": "|".join(sorted(adf["impact_lane"].unique())),
                "n_candidates": n,
                "endpoint": config.get("endpoint", "candidate_review"),
                "threshold_direction": direction,
                "threshold_successes": threshold,
                "confirmatory_threshold_successes": confirmatory_threshold,
                "null_rate": null_rate,
                "target_rate": target_rate,
                "mean_lane_prior": mean_prior,
                "mean_claim_unlock_score": mean_unlock,
                "mean_expected_information_gain": mean_eig,
                "false_unlock_risk_under_null": false_unlock,
                "power_at_target_rate": power_target,
                "power_at_mean_prior": power_prior,
                "confirmatory_false_unlock_risk_under_null": confirmatory_false_unlock,
                "confirmatory_power_at_target_rate": confirmatory_power_target,
                "confirmatory_power_at_mean_prior": confirmatory_power_prior,
                "claim_layer_after_unlock": config.get("claim_layer", "supporting claim"),
                "primary_risk": config.get("primary_risk", "unspecified"),
                "preregistered_unlock_rule": f"{direction.replace('_', ' ')} {threshold}/{n} candidates meet endpoint",
                "confirmatory_unlock_rule": f"{direction.replace('_', ' ')} {confirmatory_threshold}/{n} candidates meet endpoint",
                "main_text_claim_tier": "confirmatory_ready"
                if confirmatory_false_unlock <= 0.05
                else "confirmatory_with_caveat"
                if confirmatory_false_unlock <= 0.10
                else "exploratory_only",
                "interpretation_if_unlocked": config.get("claim_layer", "supporting claim"),
                "interpretation_if_not_unlocked": "retain computational/diagnostic boundary; do not promote claim",
            }
        )
    out = pd.DataFrame(rows)
    out["endpoint_priority_score"] = (
        0.42 * out["power_at_mean_prior"]
        + 0.28 * out["mean_expected_information_gain"]
        + 0.20 * out["mean_claim_unlock_score"]
        - 0.10 * out["false_unlock_risk_under_null"]
    )
    return out.sort_values("endpoint_priority_score", ascending=False)


def build_power_curves(endpoint: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rates = np.linspace(0.0, 1.0, 101)
    for _, row in endpoint.iterrows():
        for rate in rates:
            rows.append(
                {
                    "plate_v4_arm": row["plate_v4_arm"],
                    "candidate_success_rate": float(rate),
                    "unlock_probability": unlock_probability(
                        int(row["n_candidates"]),
                        int(row["threshold_successes"]),
                        float(rate),
                        row["threshold_direction"],
                    ),
                    "threshold_successes": int(row["threshold_successes"]),
                    "n_candidates": int(row["n_candidates"]),
                    "endpoint": row["endpoint"],
                }
            )
    return pd.DataFrame(rows)


def build_endpoint_wellmap(wellmap: pd.DataFrame, endpoint: pd.DataFrame) -> pd.DataFrame:
    merged = wellmap.merge(
        endpoint[
            [
                "plate_v4_arm",
                "endpoint",
                "preregistered_unlock_rule",
                "claim_layer_after_unlock",
                "false_unlock_risk_under_null",
                "power_at_mean_prior",
            ]
        ],
        on="plate_v4_arm",
        how="left",
    )
    merged["result_entry_placeholder"] = ""
    merged["predefined_result_code"] = np.select(
        [
            merged["assay_condition"].eq("MUT_pMHC"),
            merged["assay_condition"].str.contains("WT|DECOY", regex=True),
            merged["assay_condition"].str.contains("TCR|PROVENANCE", regex=True),
        ],
        [
            "POS/NEG binding call",
            "specificity control call",
            "manual/TCR/MD resolution call",
        ],
        default="QC call",
    )
    return merged


def build_reviewer_safeguards(endpoint: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "reviewer_risk": "public/overlap leakage",
            "v5_safeguard": "separate clean discovery from overlap-blocked positive controls",
            "evidence_file": "wetlab_plate_v4_value_of_information.tsv",
            "residual_risk": "positive-control arm cannot support benchmark novelty",
            "severity": "high",
        },
        {
            "reviewer_risk": "overfit fine-tuning",
            "v5_safeguard": "score booster selected only after source/HLA clean OOF; failed method-matrix stacker retained as negative diagnostic",
            "evidence_file": "score_booster_cv_metrics.tsv + v4 summary",
            "residual_risk": "small clean positive set n=35",
            "severity": "medium",
        },
        {
            "reviewer_risk": "no wetlab interpretation rule",
            "v5_safeguard": "binomial endpoint plan plus predeclared go/no-go thresholds",
            "evidence_file": "preregistered_endpoint_plan_v5.tsv",
            "residual_risk": "assay-specific thresholds still need lab calibration",
            "severity": "medium",
        },
        {
            "reviewer_risk": "TCR/MD evidence overclaimed",
            "v5_safeguard": "TCR/MD arm unlocks mechanistic support only, not clean benchmark novelty or clinical claim",
            "evidence_file": "claim_unlock_ladder_v5.tsv",
            "residual_risk": "requires TCR readout or repeat MD to promote mechanism",
            "severity": "medium",
        },
        {
            "reviewer_risk": "label-noise rescue relabels negatives too casually",
            "v5_safeguard": "rescue requires provenance audit plus pMHC binding",
            "evidence_file": "assay_96well_map_preregistered_v5.tsv",
            "residual_risk": "provenance may fail and rescue claim must be dropped",
            "severity": "medium",
        },
        {
            "reviewer_risk": "clinical utility overclaim",
            "v5_safeguard": "clinical/patient layer remains locked in claim ladder",
            "evidence_file": "claim_unlock_ladder_v5.tsv",
            "residual_risk": "patient metadata and external validation still missing",
            "severity": "high",
        },
    ]
    out = pd.DataFrame(rows)
    false_sum = endpoint["false_unlock_risk_under_null"].sum()
    confirmatory_false_sum = endpoint["confirmatory_false_unlock_risk_under_null"].sum()
    out.loc[len(out)] = {
        "reviewer_risk": "multiple endpoint fishing",
        "v5_safeguard": f"exploratory thresholds fixed; stricter confirmatory thresholds added; null-risk sums exploratory={false_sum:.3f}, confirmatory={confirmatory_false_sum:.3f}",
        "evidence_file": "endpoint_power_curves_v5.tsv",
        "residual_risk": "descriptive exploratory package; not a pivotal clinical trial",
        "severity": "medium",
    }
    return out


def build_claim_ladder_v5(ladder_v4: pd.DataFrame, endpoint: pd.DataFrame) -> pd.DataFrame:
    ladder = ladder_v4.copy()
    power_lookup = endpoint.set_index("plate_v4_arm")["power_at_mean_prior"].to_dict()
    risk_lookup = endpoint.set_index("plate_v4_arm")["false_unlock_risk_under_null"].to_dict()
    conf_power_lookup = endpoint.set_index("plate_v4_arm")["confirmatory_power_at_mean_prior"].to_dict()
    conf_risk_lookup = endpoint.set_index("plate_v4_arm")["confirmatory_false_unlock_risk_under_null"].to_dict()
    ladder["v5_statistical_status"] = ladder["claim_layer"].map(
        {
            "computational score recovery": "already verified",
            "clean antigen discovery": f"exploratory power {power_lookup.get('A_clean_discovery', np.nan):.3f}, risk {risk_lookup.get('A_clean_discovery', np.nan):.3f}; confirmatory power {conf_power_lookup.get('A_clean_discovery', np.nan):.3f}, risk {conf_risk_lookup.get('A_clean_discovery', np.nan):.3f}",
            "TCR/MD mechanism": f"exploratory power {power_lookup.get('B_mechanism_TCR_MD', np.nan):.3f}, risk {risk_lookup.get('B_mechanism_TCR_MD', np.nan):.3f}; confirmatory power {conf_power_lookup.get('B_mechanism_TCR_MD', np.nan):.3f}, risk {conf_risk_lookup.get('B_mechanism_TCR_MD', np.nan):.3f}",
            "label-noise rescue": f"exploratory power {power_lookup.get('C_label_rescue', np.nan):.3f}, risk {risk_lookup.get('C_label_rescue', np.nan):.3f}; confirmatory power {conf_power_lookup.get('C_label_rescue', np.nan):.3f}, risk {conf_risk_lookup.get('C_label_rescue', np.nan):.3f}",
            "specificity moat": f"exploratory power {power_lookup.get('D_specificity_moat', np.nan):.3f}, risk {risk_lookup.get('D_specificity_moat', np.nan):.3f}; confirmatory power {conf_power_lookup.get('D_specificity_moat', np.nan):.3f}, risk {conf_risk_lookup.get('D_specificity_moat', np.nan):.3f}",
            "clinical/patient relevance": "locked; no assay endpoint unlocks this alone",
        }
    ).fillna("supporting")
    return ladder


def build_result_template(plate: pd.DataFrame, endpoint: pd.DataFrame) -> pd.DataFrame:
    merged = plate.merge(
        endpoint[
            [
                "plate_v4_arm",
                "endpoint",
                "preregistered_unlock_rule",
                "confirmatory_unlock_rule",
                "main_text_claim_tier",
                "interpretation_if_unlocked",
                "interpretation_if_not_unlocked",
            ]
        ],
        on="plate_v4_arm",
        how="left",
    )
    merged["observed_candidate_success"] = ""
    merged["observed_specificity_pass"] = ""
    merged["observed_tcr_or_provenance_pass"] = ""
    merged["predeclared_interpretation"] = np.where(
        merged["impact_lane"].eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
        "positive assay QC only; do not count as novelty evidence",
        merged["interpretation_if_unlocked"],
    )
    cols = [
        "plate_v4_slot",
        "plate_v4_arm",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "impact_lane",
        "endpoint",
        "preregistered_unlock_rule",
        "observed_candidate_success",
        "observed_specificity_pass",
        "observed_tcr_or_provenance_pass",
        "main_text_claim_tier",
        "confirmatory_unlock_rule",
        "predeclared_interpretation",
        "interpretation_if_not_unlocked",
        "claim_boundary_v3",
    ]
    return merged[safe_cols(merged, cols)]


def make_figures(endpoint: pd.DataFrame, curves: pd.DataFrame, wellmap: pd.DataFrame, safeguards: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures = []
    colors = {
        "A_clean_discovery": "#15803d",
        "B_mechanism_TCR_MD": "#d97706",
        "C_label_rescue": "#b45309",
        "D_specificity_moat": "#7c3aed",
        "E_positive_QC_control": "#64748b",
        "F_model_boundary": "#dc2626",
    }

    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    for arm, adf in curves.groupby("plate_v4_arm"):
        ax.plot(adf["candidate_success_rate"], adf["unlock_probability"], label=arm.replace("_", " "), color=colors.get(arm), linewidth=2)
    ax.set_xlabel("Candidate-level success rate")
    ax.set_ylabel("Probability endpoint unlocks")
    ax.set_title("Pre-registered arm-level power curves")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig1_preregistered_power_curves.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    e = endpoint.sort_values("endpoint_priority_score")
    ax.barh(e["plate_v4_arm"].str.replace("_", " ", regex=False), e["power_at_mean_prior"], color="#2563eb", label="power at mean prior")
    ax.scatter(e["false_unlock_risk_under_null"], e["plate_v4_arm"].str.replace("_", " ", regex=False), color="#dc2626", label="null unlock risk", zorder=3)
    ax.set_xlabel("Probability")
    ax.set_title("Power vs false-unlock risk by arm")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig2_power_vs_false_unlock_risk.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    counts = pd.crosstab(wellmap["plate_v4_arm"], wellmap["assay_condition"])
    counts = counts.loc[[arm for arm in colors if arm in counts.index]]
    bottom = np.zeros(len(counts))
    condition_colors = {
        "MUT_pMHC": "#22c55e",
        "WT_CONTROL": "#60a5fa",
        "DECOY_CONTROL": "#818cf8",
        "TECH_REPLICATE": "#facc15",
        "TCR_READOUT": "#f97316",
        "QC_REPLICATE": "#facc15",
        "NEG_REPLICATE": "#a78bfa",
        "TCR_REVIEW": "#ef4444",
        "PROVENANCE_REVIEW": "#fb923c",
    }
    for cond in counts.columns:
        ax.bar(
            counts.index.str.replace("_", " ", regex=False),
            counts[cond],
            bottom=bottom,
            label=cond.replace("_", " "),
            color=condition_colors.get(cond, "#94a3b8"),
        )
        bottom += counts[cond].to_numpy()
    ax.set_ylabel("Well count")
    ax.set_title("96-well endpoint allocation")
    ax.tick_params(axis="x", rotation=28)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig3_endpoint_well_allocation.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(9.8, 5.6))
    severity_order = {"high": 3, "medium": 2, "low": 1}
    safeguards["severity_score"] = safeguards["severity"].map(severity_order).fillna(1)
    y = np.arange(len(safeguards))
    ax.scatter(safeguards["severity_score"], y, s=220, color="#f59e0b", edgecolor="white")
    ax.set_yticks(y)
    ax.set_yticklabels(safeguards["reviewer_risk"], fontsize=8)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["low", "medium", "high"])
    ax.set_xlabel("Residual risk severity")
    ax.set_title("Reviewer-risk safeguards after preregistration")
    ax.grid(axis="x", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig4_reviewer_risk_safeguards.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)
    return figures


def html_table(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    cols = safe_cols(df, cols)
    view = df[cols].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_report(summary: dict, endpoint: pd.DataFrame, ladder: pd.DataFrame, rules: pd.DataFrame, safeguards: pd.DataFrame) -> str:
    return f"""# CROSS-Neo preregistered impact v5

Generated: {summary["generated_at"]}

## Impact move

v5 turns the 96-well plan into a preregistered statistics package. Each experimental arm now has a fixed endpoint, unlock threshold, estimated power, and false-unlock risk before any assay result is observed.

## Headline metrics

- Endpoint arms: {summary["endpoint_arms"]}
- Candidates: {summary["plate_candidates"]}
- Wells: {summary["well_count"]}
- Mean arm power at lane prior: {summary["mean_power_at_prior"]:.3f}
- Mean false-unlock risk under null: {summary["mean_false_unlock_risk"]:.3f}
- Family sum of arm-level null unlock risks: {summary["sum_false_unlock_risk"]:.3f}
- Confirmatory family null-risk sum: {summary["confirmatory_sum_false_unlock_risk"]:.3f}

## Endpoint plan

{endpoint.to_markdown(index=False, floatfmt=".3f")}

## Claim ladder v5

{ladder.to_markdown(index=False)}

## Assay rules

{rules.to_markdown(index=False)}

## Reviewer safeguards

{safeguards.drop(columns=[c for c in ["severity_score"] if c in safeguards.columns]).to_markdown(index=False)}

## Decision

1. This is now a pre-specified experimental-statistical package, not a post-hoc candidate list.
2. It raises impact because a single 96-well run can unlock clean discovery, TCR/MD mechanism, label-noise rescue, specificity, and assay QC claims under fixed rules.
3. Clinical vaccine-selection remains explicitly locked.
"""


def make_html(summary: dict, endpoint: pd.DataFrame, curves: pd.DataFrame, wellmap: pd.DataFrame, ladder: pd.DataFrame, rules: pd.DataFrame, safeguards: pd.DataFrame, result_template: pd.DataFrame, figures: list[Path]) -> str:
    assets = {fig.name: f"assets/cross_neo_preregistered_v5/{fig.name}" for fig in figures}
    stat_cards = [
        ("Endpoint arms", summary["endpoint_arms"], "predefined"),
        ("Candidates", summary["plate_candidates"], "plate v5"),
        ("Wells", summary["well_count"], "96-well map"),
        ("Mean power", summary["mean_power_at_prior"], "at lane prior"),
        ("Mean null risk", summary["mean_false_unlock_risk"], "false unlock"),
        ("Conf. null sum", summary["confirmatory_sum_false_unlock_risk"], "main-text thresholds"),
        ("Clinical claim", "locked", "patient metadata needed"),
    ]
    cards = "\n".join(
        f"<div class='stat'><b>{v:.3f}</b><span>{k}</span><small>{n}</small></div>"
        if isinstance(v, float)
        else f"<div class='stat'><b>{v}</b><span>{k}</span><small>{n}</small></div>"
        for k, v, n in stat_cards
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo Preregistered Impact v5</title>
<style>
:root {{ --bg:#07111f; --panel:#111827; --line:#263244; --text:#e5e7eb; --muted:#94a3b8; --gold:#f5c542; }}
body {{ margin:0; background:var(--bg); color:var(--text); font-family:Inter, system-ui, -apple-system, Segoe UI, sans-serif; }}
header {{ padding:42px 5vw 30px; background:#0f172a; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); font-size:12px; font-weight:800; letter-spacing:.08em; text-transform:uppercase; }}
h1 {{ font-size:clamp(34px,5vw,66px); line-height:1.02; margin:10px 0 14px; }}
.lead {{ color:#cbd5e1; max-width:990px; font-size:18px; line-height:1.5; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:10px; margin-top:24px; }}
.stat {{ background:#111827; border:1px solid var(--line); border-radius:8px; padding:14px; }}
.stat b {{ display:block; font-size:24px; }}
.stat span {{ display:block; font-size:13px; margin-top:4px; }}
.stat small {{ display:block; color:var(--muted); font-size:11px; margin-top:5px; }}
main {{ display:grid; grid-template-columns:260px minmax(0,1fr); gap:28px; padding:30px 5vw 64px; }}
nav {{ position:sticky; top:18px; align-self:start; border-right:1px solid var(--line); padding-right:16px; }}
nav a {{ display:block; color:#cbd5e1; text-decoration:none; padding:7px 0; font-size:14px; }}
section {{ margin-bottom:30px; }}
.panel {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:18px; }}
h2 {{ margin:0 0 14px; }}
.grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:18px; }}
img {{ max-width:100%; border-radius:8px; background:white; }}
.table-wrap {{ overflow:auto; max-height:560px; border:1px solid var(--line); border-radius:8px; }}
.data-table {{ width:100%; border-collapse:collapse; font-size:12px; }}
.data-table th,.data-table td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
.data-table th {{ color:var(--gold); background:#0f172a; }}
.callout {{ border-left:4px solid var(--gold); padding:12px 14px; background:#0b1220; color:#dbeafe; }}
@media(max-width:900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:static; border-right:0; border-bottom:1px solid var(--line); }} .stats,.grid2 {{ grid-template-columns:1fr 1fr; }} }}
@media(max-width:600px) {{ .stats,.grid2 {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
<div class="kicker">CROSS-Neo · Preregistered impact v5 · 2026-05-10</div>
<h1>Assay statistics before the assay</h1>
<p class="lead">v5 fixes endpoints, thresholds, power, and false-unlock risk before wetlab execution. This converts the portfolio into a reviewer-defensible preregistered experimental package.</p>
<div class="stats">{cards}</div>
</header>
<main>
<nav>
<a href="#decision">Decision</a>
<a href="#figures">Figures</a>
<a href="#endpoints">Endpoints</a>
<a href="#result-template">Result Template</a>
<a href="#ladder">Claim Ladder</a>
<a href="#safeguards">Safeguards</a>
</nav>
<div>
<section id="decision" class="panel">
<h2>Decision</h2>
<div class="callout">This is now pre-specified: if the plate succeeds, claim unlocks are already defined; if it fails, the downgrade path is also defined. Clinical selection remains locked.</div>
</section>
<section id="figures">
<h2>Figures</h2>
<div class="grid2">
<div class="panel"><img src="{assets['fig1_preregistered_power_curves.png']}" alt="power curves"></div>
<div class="panel"><img src="{assets['fig2_power_vs_false_unlock_risk.png']}" alt="power risk"></div>
<div class="panel"><img src="{assets['fig3_endpoint_well_allocation.png']}" alt="well allocation"></div>
<div class="panel"><img src="{assets['fig4_reviewer_risk_safeguards.png']}" alt="reviewer safeguards"></div>
</div>
</section>
<section id="endpoints" class="panel">
<h2>Endpoint Plan</h2>
<div class="table-wrap">{html_table(endpoint, ["plate_v4_arm","n_candidates","endpoint","preregistered_unlock_rule","confirmatory_unlock_rule","main_text_claim_tier","false_unlock_risk_under_null","confirmatory_false_unlock_risk_under_null","power_at_mean_prior","confirmatory_power_at_mean_prior","claim_layer_after_unlock","primary_risk"], 12)}</div>
</section>
<section id="result-template" class="panel">
<h2>Result Template</h2>
<div class="table-wrap">{html_table(result_template, ["plate_v4_slot","plate_v4_arm","candidate_id","peptide","hla_allele_4digit","endpoint","preregistered_unlock_rule","confirmatory_unlock_rule","main_text_claim_tier","observed_candidate_success","predeclared_interpretation","claim_boundary_v3"], 24)}</div>
</section>
<section id="ladder" class="panel">
<h2>Claim Ladder v5</h2>
<div class="table-wrap">{html_table(ladder, ["ladder_step","claim_layer","current_status","v5_statistical_status","unlock_rule","claim_after_unlock","claim_boundary"], 12)}</div>
</section>
<section id="safeguards" class="panel">
<h2>Reviewer Safeguards</h2>
<div class="table-wrap">{html_table(safeguards, ["reviewer_risk","v5_safeguard","residual_risk","severity"], 12)}</div>
</section>
</div>
</main>
</body>
</html>
"""


def copy_to_hub(figures: list[Path], html: str) -> Path:
    HUB_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for fig in figures:
        shutil.copy2(fig, HUB_ASSET_DIR / fig.name)
    html_path = HUB_DIR / "cross_neo_preregistered_impact_v5.html"
    html_path.write_text(html)
    if LIVE_HUB_DIR.exists():
        LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
        for fig in figures:
            shutil.copy2(fig, LIVE_ASSET_DIR / fig.name)
        shutil.copy2(html_path, LIVE_HUB_DIR / html_path.name)
    return html_path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plate = read_tsv(PLATE_V4, required=True)
    wellmap = read_tsv(WELLMAP_V4, required=True)
    ladder_v4 = read_tsv(LADDER_V4, required=True)
    rules_v4 = read_tsv(RULES_V4, required=True)
    summary_v4 = json.loads(SUMMARY_V4.read_text()) if SUMMARY_V4.exists() else {}

    endpoint = build_endpoint_plan(plate)
    curves = build_power_curves(endpoint)
    endpoint_wellmap = build_endpoint_wellmap(wellmap, endpoint)
    ladder = build_claim_ladder_v5(ladder_v4, endpoint)
    safeguards = build_reviewer_safeguards(endpoint)
    result_template = build_result_template(plate, endpoint)
    figures = make_figures(endpoint, curves, endpoint_wellmap, safeguards)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(OUT_DIR.relative_to(ROOT)),
        "endpoint_arms": int(len(endpoint)),
        "plate_candidates": int(len(plate)),
        "well_count": int(len(endpoint_wellmap)),
        "mean_power_at_prior": float(endpoint["power_at_mean_prior"].mean()),
        "mean_false_unlock_risk": float(endpoint["false_unlock_risk_under_null"].mean()),
        "sum_false_unlock_risk": float(endpoint["false_unlock_risk_under_null"].sum()),
        "confirmatory_mean_power_at_prior": float(endpoint["confirmatory_power_at_mean_prior"].mean()),
        "confirmatory_mean_false_unlock_risk": float(endpoint["confirmatory_false_unlock_risk_under_null"].mean()),
        "confirmatory_sum_false_unlock_risk": float(endpoint["confirmatory_false_unlock_risk_under_null"].sum()),
        "max_endpoint_priority_arm": endpoint.sort_values("endpoint_priority_score", ascending=False).iloc[0]["plate_v4_arm"],
        "score_recovery_vs_method_matrix": float(summary_v4.get("score_recovery_vs_method_matrix", 0.0)),
        "top_endpoint_plan": endpoint.head(6)[
            ["plate_v4_arm", "endpoint_priority_score", "power_at_mean_prior", "false_unlock_risk_under_null"]
        ].to_dict("records"),
    }

    write_tsv(endpoint, "preregistered_endpoint_plan_v5.tsv")
    write_tsv(curves, "endpoint_power_curves_v5.tsv")
    write_tsv(endpoint_wellmap, "assay_96well_map_preregistered_v5.tsv")
    write_tsv(ladder, "claim_unlock_ladder_v5.tsv")
    write_tsv(rules_v4, "assay_decision_rules_v5.tsv")
    write_tsv(safeguards, "reviewer_safeguards_v5.tsv")
    write_tsv(result_template, "assay_result_template_v5.tsv")
    (OUT_DIR / "PREREGISTERED_IMPACT_V5_REPORT_KR.md").write_text(make_report(summary, endpoint, ladder, rules_v4, safeguards))

    html = make_html(summary, endpoint, curves, endpoint_wellmap, ladder, rules_v4, safeguards, result_template, figures)
    html_path = copy_to_hub(figures, html)
    summary["html_path"] = str(html_path.relative_to(ROOT))
    (OUT_DIR / "preregistered_impact_v5_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
