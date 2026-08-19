#!/usr/bin/env python3
"""Build a high-impact CROSS-Neo candidate portfolio and plate v3 package.

This turns the scoring/fine-tune outputs into a reviewer-facing experimental
portfolio: clean BAR-Neo leads, TCR/MD structural leads, label-noise rescue
cases, and assay controls. It keeps claim boundaries explicit.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10"
BOOST_DIR = BASE / "score_booster_finetune_2026_05_10"
P0_DIR = BASE / "merged_p0_experiments"
OUT_DIR = BASE / "impact_portfolio_v3_2026_05_10"
FIG_DIR = OUT_DIR / "figures"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_impact_v3"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_impact_v3"

SCORES = BOOST_DIR / "score_booster_candidate_scores.tsv"
BOOST_CV = BOOST_DIR / "score_booster_cv_metrics.tsv"
MATRIX_CV = BASE / "score_finetune_2026_05_10/clean_cv_finetune_metrics.tsv"
WETLAB_V2 = P0_DIR / "cross_neo_wetlab_plate_v2.tsv"
TCR_DISAGREE = P0_DIR / "tcr_structure_disagreement_queue.tsv"
LABEL_NOISE = P0_DIR / "cross_neo_label_noise_audit_queue.tsv"


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


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="float64")
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def text(df: pd.DataFrame, col: str, default: str = "") -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="object")
    return df[col].fillna(default).astype(str)


def truthy(df: pd.DataFrame, col: str, default: bool = False) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype="bool")
    raw = df[col]
    if raw.dtype == bool:
        return raw.fillna(default)
    filled = raw.astype("object").where(raw.notna(), default)
    return filled.astype(str).str.lower().str.strip().isin({"1", "true", "t", "yes", "y"})


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def clip01(values):
    return np.clip(values, 0.0, 1.0)


def prepare_wetlab(wetlab: pd.DataFrame) -> pd.DataFrame:
    if wetlab.empty:
        return pd.DataFrame()
    w = wetlab.copy()
    if "row_id" in w.columns:
        if "candidate_id" in w.columns:
            w = w.drop(columns=["candidate_id"])
        w = w.rename(columns={"row_id": "candidate_id"})
    keep = [
        "candidate_id",
        "plate_v2_order_score",
        "plate_v2_tier",
        "plate_v2_assay_bundle",
        "tcr_augmented_score_mean",
        "paired_tcr_evidence_count",
        "tcr_evidence_count",
        "md_label",
        "md_score",
        "control_readiness_score",
        "wt_or_decoy_ready",
        "main_dl_score",
        "baker_structural_score",
        "prediction_outcome",
        "next_action",
        "wildtype_peptide",
        "wt_status",
        "anchor_preserved_decoy",
        "wt_status_control",
        "anchor_preserved_decoy_control",
    ]
    return w[safe_cols(w, keep)].drop_duplicates("candidate_id")


def prepare_tcr_disagreement(tcr: pd.DataFrame) -> pd.DataFrame:
    if tcr.empty:
        return pd.DataFrame()
    t = tcr.copy()
    if "row_id" in t.columns:
        if "candidate_id" in t.columns:
            t = t.drop(columns=["candidate_id"])
        t = t.rename(columns={"row_id": "candidate_id"})
    keep = [
        "candidate_id",
        "tcr_structure_disagreement_score",
        "tcr_structure_action",
        "pmhc_low_tcr_high",
        "pmhc_high_tcr_low",
        "tcr_high_md_low_or_missing",
        "model_high_md_low",
        "md_high_model_low",
    ]
    return t[safe_cols(t, keep)].drop_duplicates("candidate_id")


def prepare_label_noise(noise: pd.DataFrame) -> pd.DataFrame:
    if noise.empty:
        return pd.DataFrame()
    keep = [
        "candidate_id",
        "label_noise_priority_score",
        "label_noise_bucket",
        "label_noise_next_action",
    ]
    return noise[safe_cols(noise, keep)].drop_duplicates("candidate_id")


def build_portfolio() -> pd.DataFrame:
    scores = read_tsv(SCORES, required=True)
    wetlab = prepare_wetlab(read_tsv(WETLAB_V2))
    tcr = prepare_tcr_disagreement(read_tsv(TCR_DISAGREE))
    noise = prepare_label_noise(read_tsv(LABEL_NOISE))

    df = scores.copy()
    for add in [wetlab, tcr, noise]:
        if not add.empty:
            df = df.merge(add, on="candidate_id", how="left")

    paired_tcr = numeric(df, "paired_tcr_evidence_count", 0.0) > 0
    any_tcr = numeric(df, "tcr_evidence_count", 0.0) > 0
    md_score = numeric(df, "md_score", 0.0)
    tcr_score = numeric(df, "tcr_augmented_score_mean", 0.0)
    control = numeric(df, "control_readiness_score", 0.0)
    dl = numeric(df, "main_dl_score", 0.0)
    baker = numeric(df, "baker_structural_score", 0.0)
    tcr_md = (
        0.34 * tcr_score
        + 0.24 * md_score
        + 0.14 * control
        + 0.12 * dl
        + 0.08 * baker
        + 0.05 * paired_tcr.astype(float)
        + 0.03 * any_tcr.astype(float)
    )
    df["tcr_md_integrated_score"] = clip01(tcr_md)

    df["source_hla_novelty_score"] = clip01(
        0.25 * truthy(df, "is_underrepresented_hla").astype(float)
        + 0.18 * truthy(df, "is_low_prevalence_context").astype(float)
        + 0.18 * truthy(df, "is_korean_hla_focus").astype(float)
        + 0.20 * (1.0 - numeric(df, "exact_peptide_train_overlap", 0.0).clip(0, 1))
        + 0.19 * numeric(df, "source_hla_generalization_cap", 0.0)
    )

    df["impact_portfolio_score"] = clip01(
        0.34 * numeric(df, "finetuned_experiment_priority_score", 0.0)
        + 0.18 * numeric(df, "bma_v2_discovery_score", 0.0)
        + 0.22 * df["tcr_md_integrated_score"]
        + 0.12 * numeric(df, "label_noise_priority_score", 0.0)
        + 0.08 * df["source_hla_novelty_score"]
        + 0.06 * control
    )
    df["claim_safe_impact_score"] = np.minimum(
        df["impact_portfolio_score"],
        numeric(df, "finetuned_claim_capped_score", 0.0),
    )

    leakage = text(df, "leakage_risk_level").str.lower()
    plate_tier = text(df, "plate_v2_tier")
    label_noise_bucket = text(df, "label_noise_bucket")
    high_clean_score = numeric(df, "finetuned_experiment_priority_score", 0.0) >= 0.78
    low_or_medium_leakage = leakage.isin(["low", "medium"])

    df["impact_lane"] = np.select(
        [
            plate_tier.eq("TIER_1_TCR_MD_CONTROL_READY"),
            low_or_medium_leakage & high_clean_score & numeric(df, "label", 0).eq(1),
            label_noise_bucket.eq("possible_false_negative_or_assay_miss")
            & (numeric(df, "finetuned_score_booster_prob", 0.0) >= 0.70),
            leakage.eq("high") & numeric(df, "label", 0).eq(1) & (numeric(df, "finetuned_score_booster_prob", 0.0) >= 0.85),
            low_or_medium_leakage & high_clean_score & numeric(df, "label", 0).eq(0),
            numeric(df, "tcr_structure_disagreement_score", 0.0) >= 0.55,
        ],
        [
            "TCR_MD_STRUCTURAL_LEAD",
            "CLEAN_BARNEO_LEAD",
            "LABEL_NOISE_RESCUE",
            "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED",
            "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT",
            "TCR_STRUCTURE_DISAGREEMENT_AUDIT",
        ],
        default="SUPPORTING_RESERVE",
    )

    df["claim_boundary_v3"] = np.select(
        [
            df["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD"),
            df["impact_lane"].eq("CLEAN_BARNEO_LEAD"),
            df["impact_lane"].eq("LABEL_NOISE_RESCUE"),
            df["impact_lane"].eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
            df["impact_lane"].eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
            df["impact_lane"].eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        ],
        [
            "TCR/MD diagnostic wetlab priority; not clean benchmark or clinical claim",
            "clean BAR-Neo experiment priority; patient metadata still blocks clinical claim",
            "rescue/audit case; could expose false-negative label or assay mismatch",
            "positive assay control only; overlap/leakage blocks benchmark claim",
            "hard-negative or false-negative audit; use to sharpen specificity",
            "TCR/structure disagreement audit; no positive claim until resolved",
        ],
        default="supporting reserve only",
    )
    df["recommended_next_step_v3"] = np.select(
        [
            df["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD"),
            df["impact_lane"].eq("CLEAN_BARNEO_LEAD"),
            df["impact_lane"].eq("LABEL_NOISE_RESCUE"),
            df["impact_lane"].eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
            df["impact_lane"].eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
            df["impact_lane"].eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        ],
        [
            "mutant-vs-WT/decoy pMHC binding + TCR tetramer + activation",
            "mutant-vs-WT/decoy pMHC binding first, then TCR expansion if positive",
            "manual provenance audit + low-cost pMHC binding rescue",
            "use as assay positive control, not model validation evidence",
            "paired negative-control assay or manual label audit",
            "repeat/extend MD or TCR specificity review before plate claim",
        ],
        default="hold as reserve",
    )

    front = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "leakage_risk_level",
        "impact_lane",
        "impact_portfolio_score",
        "claim_safe_impact_score",
        "finetuned_score_booster_prob",
        "finetuned_experiment_priority_score",
        "finetuned_claim_capped_score",
        "bma_v2_discovery_score",
        "bma_v2_claim_safe_score",
        "tcr_md_integrated_score",
        "tcr_augmented_score_mean",
        "paired_tcr_evidence_count",
        "md_label",
        "md_score",
        "control_readiness_score",
        "label_noise_priority_score",
        "label_noise_bucket",
        "source_hla_novelty_score",
        "primary_claim_blocker",
        "claim_boundary_v3",
        "recommended_next_step_v3",
        "plate_v2_tier",
        "plate_v2_order_score",
    ]
    front = safe_cols(df, front)
    return df.sort_values(["impact_portfolio_score", "finetuned_experiment_priority_score"], ascending=False)[
        front + [col for col in df.columns if col not in front]
    ]


def pick_rows(
    portfolio: pd.DataFrame,
    lane: str,
    n: int,
    exclude: set[str],
    require_label: int | None = None,
    sort_col: str = "impact_portfolio_score",
) -> pd.DataFrame:
    sub = portfolio[portfolio["impact_lane"].eq(lane)].copy()
    if require_label is not None and "label" in sub.columns:
        sub = sub[numeric(sub, "label", -1).eq(require_label)]
    sub = sub[~sub["candidate_id"].isin(exclude)]
    sub = sub.sort_values(sort_col, ascending=False).head(n)
    exclude.update(sub["candidate_id"].tolist())
    return sub


def build_plate_v3(portfolio: pd.DataFrame) -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    exclude: set[str] = set()
    plan = [
        ("CLEAN_BARNEO_LEAD", 6, 1, "finetuned_experiment_priority_score"),
        ("TCR_MD_STRUCTURAL_LEAD", 2, None, "plate_v2_order_score"),
        ("LABEL_NOISE_RESCUE", 4, 0, "label_noise_priority_score"),
        ("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT", 4, 0, "impact_portfolio_score"),
        ("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED", 4, 1, "finetuned_experiment_priority_score"),
        ("TCR_STRUCTURE_DISAGREEMENT_AUDIT", 4, None, "tcr_structure_disagreement_score"),
    ]
    for lane, n, label, sort_col in plan:
        picked = pick_rows(portfolio, lane, n, exclude, label, sort_col)
        selected.append(picked)

    plate = pd.concat([df for df in selected if not df.empty], ignore_index=True)
    if len(plate) < 24:
        reserve = portfolio[~portfolio["candidate_id"].isin(exclude)].copy()
        reserve = reserve.sort_values("impact_portfolio_score", ascending=False).head(24 - len(plate))
        reserve["impact_lane"] = reserve["impact_lane"].where(
            reserve["impact_lane"].ne("SUPPORTING_RESERVE"), "SUPPORTING_RESERVE_FILLER"
        )
        plate = pd.concat([plate, reserve], ignore_index=True)

    plate = plate.head(24).copy()
    plate["plate_v3_slot"] = np.arange(1, len(plate) + 1)
    plate["plate_v3_block"] = np.select(
        [
            plate["impact_lane"].eq("CLEAN_BARNEO_LEAD"),
            plate["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD"),
            plate["impact_lane"].eq("LABEL_NOISE_RESCUE"),
            plate["impact_lane"].eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
            plate["impact_lane"].eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
            plate["impact_lane"].eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        ],
        [
            "A_clean_generalizable_antigen",
            "B_TCR_MD_mechanistic",
            "C_false_negative_rescue",
            "D_hard_negative_specificity",
            "E_assay_positive_control",
            "F_TCR_structure_audit",
        ],
        default="G_reserve",
    )
    plate["plate_v3_assay_bundle"] = np.select(
        [
            plate["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD"),
            plate["impact_lane"].eq("CLEAN_BARNEO_LEAD"),
            plate["impact_lane"].eq("LABEL_NOISE_RESCUE"),
            plate["impact_lane"].eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
            plate["impact_lane"].eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
            plate["impact_lane"].eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        ],
        [
            "mutant_vs_WT_decoy_pMHC_binding + TCR_tetramer + activation",
            "mutant_vs_WT_decoy_pMHC_binding; TCR expansion if binder",
            "provenance_audit + low_cost_pMHC_binding_rescue",
            "matched_negative_specificity_control",
            "assay_QC_positive_control_only",
            "TCR_specificity_review + repeat_MD_or_pMHC_binding_triage",
        ],
        default="reserve_review",
    )
    front = [
        "plate_v3_slot",
        "plate_v3_block",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "leakage_risk_level",
        "impact_lane",
        "impact_portfolio_score",
        "finetuned_experiment_priority_score",
        "tcr_md_integrated_score",
        "label_noise_priority_score",
        "plate_v3_assay_bundle",
        "claim_boundary_v3",
        "recommended_next_step_v3",
    ]
    front = safe_cols(plate, front)
    return plate[front + [col for col in plate.columns if col not in front]]


def get_score_summary(portfolio: pd.DataFrame) -> dict:
    boost_cv = read_tsv(BOOST_CV)
    matrix_cv = read_tsv(MATRIX_CV)
    boost_oof = boost_cv[boost_cv["protocol"].isin(["source_heldout_oof", "hla_heldout_oof"])].copy()
    selected = boost_oof[
        boost_oof["feature_set"].eq("stress_only") & boost_oof["model_name"].eq("logistic_C1.0")
    ]
    matrix_best = None
    if not matrix_cv.empty:
        mm = matrix_cv[matrix_cv["protocol"].isin(["source_heldout_oof", "hla_heldout_oof"])]
        if not mm.empty:
            matrix_best = float(mm["auprc"].max())
    summary = {
        "score_booster_source_auprc": float(
            selected.loc[selected["protocol"].eq("source_heldout_oof"), "auprc"].iloc[0]
        )
        if not selected.empty and selected["protocol"].eq("source_heldout_oof").any()
        else None,
        "score_booster_hla_auprc": float(
            selected.loc[selected["protocol"].eq("hla_heldout_oof"), "auprc"].iloc[0]
        )
        if not selected.empty and selected["protocol"].eq("hla_heldout_oof").any()
        else None,
        "score_booster_mean_auprc": float(selected["auprc"].mean()) if not selected.empty else None,
        "score_booster_mean_auroc": float(selected["auroc"].mean()) if not selected.empty else None,
        "method_matrix_best_oof_auprc": matrix_best,
        "portfolio_rows": int(len(portfolio)),
        "lane_counts": portfolio["impact_lane"].value_counts().to_dict(),
    }
    if matrix_best and summary["score_booster_mean_auprc"]:
        summary["auprc_recovery_vs_method_matrix"] = summary["score_booster_mean_auprc"] / matrix_best
    return summary


def make_figures(portfolio: pd.DataFrame, plate: pd.DataFrame, summary: dict) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")
    figures: list[Path] = []

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    labels = ["Method-matrix\nbest OOF", "Score booster\nsource OOF", "Score booster\nHLA OOF"]
    values = [
        summary.get("method_matrix_best_oof_auprc") or 0,
        summary.get("score_booster_source_auprc") or 0,
        summary.get("score_booster_hla_auprc") or 0,
    ]
    colors = ["#7a869a", "#1f9d8a", "#1f78b4"]
    bars = ax.bar(labels, values, color=colors)
    ax.set_ylim(0, max(values) * 1.25 if values else 1)
    ax.set_ylabel("AUPRC")
    ax.set_title("Clean-CV score recovery")
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.015, f"{val:.3f}", ha="center", va="bottom", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig1_clean_cv_score_recovery.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    lanes = [
        "CLEAN_BARNEO_LEAD",
        "TCR_MD_STRUCTURAL_LEAD",
        "LABEL_NOISE_RESCUE",
        "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT",
        "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED",
        "TCR_STRUCTURE_DISAGREEMENT_AUDIT",
        "SUPPORTING_RESERVE",
    ]
    palette = {
        "CLEAN_BARNEO_LEAD": "#1f9d55",
        "TCR_MD_STRUCTURAL_LEAD": "#d97706",
        "LABEL_NOISE_RESCUE": "#b45309",
        "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT": "#7c3aed",
        "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED": "#64748b",
        "TCR_STRUCTURE_DISAGREEMENT_AUDIT": "#dc2626",
        "SUPPORTING_RESERVE": "#94a3b8",
    }
    for lane in lanes:
        sub = portfolio[portfolio["impact_lane"].eq(lane)]
        if sub.empty:
            continue
        plot_sub = sub.sample(min(len(sub), 450), random_state=13) if len(sub) > 450 else sub
        ax.scatter(
            numeric(plot_sub, "finetuned_experiment_priority_score", 0),
            numeric(plot_sub, "tcr_md_integrated_score", 0),
            s=28 + 90 * numeric(plot_sub, "label_noise_priority_score", 0),
            alpha=0.68,
            label=lane.replace("_", " "),
            color=palette[lane],
            edgecolor="white",
            linewidth=0.3,
        )
    ax.set_xlabel("Fine-tuned experiment-priority score")
    ax.set_ylabel("TCR/MD integrated score")
    ax.set_title("Two-axis vaccine impact portfolio")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig2_two_axis_impact_portfolio.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    counts = plate["plate_v3_block"].value_counts().sort_index()
    ax.barh(counts.index.str.replace("_", " ", regex=False), counts.values, color="#2563eb")
    ax.set_xlabel("Plate v3 candidate count")
    ax.set_title("Plate v3 balances discovery, mechanism, rescue, and controls")
    for y, val in enumerate(counts.values):
        ax.text(val + 0.08, y, str(val), va="center")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig3_plate_v3_composition.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    heat = pd.crosstab(portfolio["impact_lane"], portfolio["finetuned_score_action"])
    heat = heat.loc[[lane for lane in lanes if lane in heat.index]]
    im = ax.imshow(np.log1p(heat.values), cmap="viridis")
    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_xticklabels([c.replace("_", " ") for c in heat.columns], rotation=35, ha="right", fontsize=8)
    ax.set_yticklabels([i.replace("_", " ") for i in heat.index], fontsize=8)
    ax.set_title("Claim boundary audit: lanes vs score actions")
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            ax.text(j, i, str(int(heat.iloc[i, j])), ha="center", va="center", color="white", fontsize=7)
    fig.colorbar(im, ax=ax, label="log1p(count)")
    fig.tight_layout()
    path = FIG_DIR / "fig4_claim_boundary_lane_heatmap.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    return figures


def html_table(df: pd.DataFrame, cols: list[str], n: int = 12) -> str:
    cols = safe_cols(df, cols)
    if not cols:
        return "<p>No rows.</p>"
    view = df[cols].head(n).copy()
    for col in view.columns:
        if pd.api.types.is_numeric_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    return view.to_html(index=False, classes="data-table", escape=False)


def make_report(portfolio: pd.DataFrame, plate: pd.DataFrame, summary: dict, figures: list[Path]) -> str:
    top_clean = portfolio[portfolio["impact_lane"].eq("CLEAN_BARNEO_LEAD")].sort_values(
        "finetuned_experiment_priority_score", ascending=False
    )
    top_struct = portfolio[portfolio["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD")].sort_values(
        "tcr_md_integrated_score", ascending=False
    )
    top_rescue = portfolio[portfolio["impact_lane"].eq("LABEL_NOISE_RESCUE")].sort_values(
        "label_noise_priority_score", ascending=False
    )

    def md_table(df: pd.DataFrame, cols: list[str], n: int = 10) -> str:
        cols = safe_cols(df, cols)
        if not cols:
            return "_no rows_"
        return df[cols].head(n).to_markdown(index=False, floatfmt=".3f")

    report = f"""# CROSS-Neo impact portfolio v3

Generated: {summary["generated_at"]}

## Impact move

Single-score story에서 빠져나와 **two-lane portfolio**로 올렸다: clean BAR-Neo lead는 generalizable antigen discovery 축, TCR/MD lead는 mechanistic wetlab 축이다. 두 축은 같은 claim으로 섞지 않고 plate v3에서 같이 검증한다.

## Score recovery

- method-matrix stacker best OOF AUPRC: {summary.get("method_matrix_best_oof_auprc", 0):.3f}
- score booster source OOF AUPRC: {summary.get("score_booster_source_auprc", 0):.3f}
- score booster HLA OOF AUPRC: {summary.get("score_booster_hla_auprc", 0):.3f}
- score booster mean AUROC: {summary.get("score_booster_mean_auroc", 0):.3f}

## Lane counts

```json
{json.dumps(summary["lane_counts"], ensure_ascii=False, indent=2)}
```

## Top clean BAR-Neo leads

{md_table(top_clean, ["candidate_id", "peptide", "hla_allele_4digit", "label", "source_name", "leakage_risk_level", "finetuned_score_booster_prob", "finetuned_experiment_priority_score", "impact_portfolio_score", "claim_boundary_v3"], 12)}

## TCR/MD structural leads

{md_table(top_struct, ["candidate_id", "peptide", "hla_allele_4digit", "label", "source_name", "tcr_md_integrated_score", "tcr_augmented_score_mean", "md_label", "md_score", "control_readiness_score", "claim_boundary_v3"], 8)}

## Label-noise rescue leads

{md_table(top_rescue, ["candidate_id", "peptide", "hla_allele_4digit", "label", "source_name", "label_noise_priority_score", "finetuned_score_booster_prob", "label_noise_bucket", "recommended_next_step_v3"], 8)}

## Plate v3

{md_table(plate, ["plate_v3_slot", "plate_v3_block", "candidate_id", "peptide", "hla_allele_4digit", "label", "impact_lane", "impact_portfolio_score", "plate_v3_assay_bundle", "claim_boundary_v3"], 24)}

## Decision

1. Impact framing: not a single immunogenicity predictor; this is a controlled vaccine-candidate portfolio with orthogonal evidence lanes.
2. Wetlab plate v3 now has discovery leads, structural/TCR leads, rescue cases, hard negatives, and positive assay controls.
3. Claim boundary remains intact: no clinical vaccine-selection claim without patient metadata and external/wetlab validation.

## Files

- `impact_candidate_portfolio_v3.tsv`
- `wetlab_plate_v3_impact_design.tsv`
- `impact_portfolio_v3_summary.json`
- `figures/fig1_clean_cv_score_recovery.png`
- `figures/fig2_two_axis_impact_portfolio.png`
- `figures/fig3_plate_v3_composition.png`
- `figures/fig4_claim_boundary_lane_heatmap.png`
"""
    return report


def make_html(portfolio: pd.DataFrame, plate: pd.DataFrame, summary: dict, figures: list[Path]) -> str:
    top_clean = portfolio[portfolio["impact_lane"].eq("CLEAN_BARNEO_LEAD")].sort_values(
        "finetuned_experiment_priority_score", ascending=False
    )
    top_struct = portfolio[portfolio["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD")].sort_values(
        "tcr_md_integrated_score", ascending=False
    )
    top_rescue = portfolio[portfolio["impact_lane"].eq("LABEL_NOISE_RESCUE")].sort_values(
        "label_noise_priority_score", ascending=False
    )
    assets = {fig.name: f"assets/cross_neo_impact_v3/{fig.name}" for fig in figures}
    stat_cards = [
        ("Score booster mean AUPRC", summary.get("score_booster_mean_auprc"), "clean source/HLA OOF"),
        ("Mean AUROC", summary.get("score_booster_mean_auroc"), "clean source/HLA OOF"),
        ("Plate v3 rows", len(plate), "balanced experimental design"),
        ("Clean leads", summary["lane_counts"].get("CLEAN_BARNEO_LEAD", 0), "generalizable antigen lane"),
        ("TCR/MD leads", summary["lane_counts"].get("TCR_MD_STRUCTURAL_LEAD", 0), "mechanistic lane"),
        ("Rescue cases", summary["lane_counts"].get("LABEL_NOISE_RESCUE", 0), "label-noise audit lane"),
    ]
    cards = "\n".join(
        f"<div class='stat'><b>{val:.3f}</b><span>{label}</span><small>{note}</small></div>"
        if isinstance(val, float)
        else f"<div class='stat'><b>{val}</b><span>{label}</span><small>{note}</small></div>"
        for label, val, note in stat_cards
    )
    lane_json = json.dumps(summary["lane_counts"], ensure_ascii=False, indent=2)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo Impact Portfolio v3</title>
<style>
:root {{ --bg:#0b1020; --panel:#111827; --muted:#94a3b8; --text:#e5e7eb; --gold:#f5c542; --line:#263244; --green:#34d399; --blue:#60a5fa; }}
body {{ margin:0; background:var(--bg); color:var(--text); font-family:Inter, system-ui, -apple-system, Segoe UI, sans-serif; line-height:1.5; }}
header {{ padding:42px 5vw 28px; border-bottom:1px solid var(--line); background:linear-gradient(180deg,#111827,#0b1020); }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; font-weight:700; }}
h1 {{ font-size:clamp(34px,5vw,68px); line-height:1.02; margin:10px 0 14px; max-width:1180px; }}
.lead {{ color:#cbd5e1; max-width:980px; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(6,minmax(120px,1fr)); gap:10px; margin-top:26px; }}
.stat {{ background:#0f172a; border:1px solid var(--line); border-radius:8px; padding:14px; }}
.stat b {{ display:block; font-size:24px; color:#fff; }}
.stat span {{ display:block; font-size:13px; margin-top:4px; }}
.stat small {{ color:var(--muted); display:block; margin-top:6px; font-size:11px; }}
main {{ display:grid; grid-template-columns:260px minmax(0,1fr); gap:28px; padding:28px 5vw 60px; }}
nav {{ position:sticky; top:18px; align-self:start; border-right:1px solid var(--line); padding-right:16px; }}
nav a {{ display:block; color:#cbd5e1; text-decoration:none; padding:7px 0; font-size:14px; }}
section {{ margin-bottom:38px; }}
h2 {{ font-size:24px; margin:0 0 14px; }}
.panel {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:18px; margin-bottom:18px; }}
.grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:18px; }}
img {{ max-width:100%; border-radius:8px; background:#fff; }}
pre {{ white-space:pre-wrap; background:#0f172a; border:1px solid var(--line); border-radius:8px; padding:14px; color:#cbd5e1; }}
.data-table {{ width:100%; border-collapse:collapse; font-size:12px; }}
.data-table th,.data-table td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; }}
.data-table th {{ color:var(--gold); background:#0f172a; position:sticky; top:0; }}
.table-wrap {{ overflow:auto; max-height:520px; border:1px solid var(--line); border-radius:8px; }}
.callout {{ border-left:4px solid var(--gold); padding:12px 14px; background:#101827; color:#dbeafe; }}
@media (max-width:900px) {{ main {{ grid-template-columns:1fr; }} nav {{ position:static; border-right:0; border-bottom:1px solid var(--line); }} .stats,.grid2 {{ grid-template-columns:1fr 1fr; }} }}
@media (max-width:600px) {{ .stats,.grid2 {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
<div class="kicker">CROSS-Neo · Kaggle-style impact upgrade · 2026-05-10</div>
<h1>Two-axis cancer vaccine candidate portfolio</h1>
<p class="lead">Clean BAR-Neo scoring, TCR/MD structural evidence, label-noise rescue, and assay controls are separated into claim-safe lanes. The impact move is portfolio design, not an unsafe single immunogenicity headline.</p>
<div class="stats">{cards}</div>
</header>
<main>
<nav>
<a href="#decision">Decision</a>
<a href="#figures">Figures</a>
<a href="#clean">Clean BAR-Neo Leads</a>
<a href="#structural">TCR/MD Leads</a>
<a href="#rescue">Label-Noise Rescue</a>
<a href="#plate">Plate v3</a>
<a href="#counts">Lane Counts</a>
</nav>
<div>
<section id="decision" class="panel">
<h2>Decision</h2>
<div class="callout">Plate v3 now tests two orthogonal claims: clean generalizable antigen discovery and TCR/MD mechanistic support. Patient metadata and overlap gates remain active, so clinical vaccine-selection claims stay blocked.</div>
</section>
<section id="figures">
<h2>Figures</h2>
<div class="grid2">
<div class="panel"><img src="{assets['fig1_clean_cv_score_recovery.png']}" alt="score recovery"></div>
<div class="panel"><img src="{assets['fig2_two_axis_impact_portfolio.png']}" alt="portfolio scatter"></div>
<div class="panel"><img src="{assets['fig3_plate_v3_composition.png']}" alt="plate composition"></div>
<div class="panel"><img src="{assets['fig4_claim_boundary_lane_heatmap.png']}" alt="claim boundary heatmap"></div>
</div>
</section>
<section id="clean" class="panel">
<h2>Clean BAR-Neo Leads</h2>
<div class="table-wrap">{html_table(top_clean, ["candidate_id","peptide","hla_allele_4digit","label","source_name","leakage_risk_level","finetuned_score_booster_prob","finetuned_experiment_priority_score","impact_portfolio_score","claim_boundary_v3"], 16)}</div>
</section>
<section id="structural" class="panel">
<h2>TCR/MD Structural Leads</h2>
<div class="table-wrap">{html_table(top_struct, ["candidate_id","peptide","hla_allele_4digit","label","source_name","tcr_md_integrated_score","tcr_augmented_score_mean","md_label","md_score","control_readiness_score","claim_boundary_v3"], 10)}</div>
</section>
<section id="rescue" class="panel">
<h2>Label-Noise Rescue</h2>
<div class="table-wrap">{html_table(top_rescue, ["candidate_id","peptide","hla_allele_4digit","label","source_name","label_noise_priority_score","finetuned_score_booster_prob","label_noise_bucket","recommended_next_step_v3"], 12)}</div>
</section>
<section id="plate" class="panel">
<h2>Wetlab Plate v3</h2>
<div class="table-wrap">{html_table(plate, ["plate_v3_slot","plate_v3_block","candidate_id","peptide","hla_allele_4digit","label","impact_lane","impact_portfolio_score","plate_v3_assay_bundle","claim_boundary_v3"], 24)}</div>
</section>
<section id="counts" class="panel">
<h2>Lane Counts</h2>
<pre>{lane_json}</pre>
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
    html_path = HUB_DIR / "cross_neo_impact_portfolio_v3.html"
    html_path.write_text(html)

    if LIVE_HUB_DIR.exists():
        LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
        for fig in figures:
            shutil.copy2(fig, LIVE_ASSET_DIR / fig.name)
        shutil.copy2(html_path, LIVE_HUB_DIR / html_path.name)
    return html_path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    portfolio = build_portfolio()
    plate = build_plate_v3(portfolio)
    summary = get_score_summary(portfolio)
    summary.update(
        {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "output_dir": str(OUT_DIR.relative_to(ROOT)),
            "plate_v3_rows": int(len(plate)),
            "top_clean_leads": portfolio[portfolio["impact_lane"].eq("CLEAN_BARNEO_LEAD")]
            .sort_values("finetuned_experiment_priority_score", ascending=False)
            .head(8)[["candidate_id", "peptide", "hla_allele_4digit", "finetuned_experiment_priority_score"]]
            .to_dict("records"),
            "top_tcr_md_leads": portfolio[portfolio["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD")]
            .sort_values("tcr_md_integrated_score", ascending=False)
            .head(4)[["candidate_id", "peptide", "hla_allele_4digit", "tcr_md_integrated_score"]]
            .to_dict("records"),
        }
    )

    write_tsv(portfolio, "impact_candidate_portfolio_v3.tsv")
    write_tsv(plate, "wetlab_plate_v3_impact_design.tsv")
    figures = make_figures(portfolio, plate, summary)
    report = make_report(portfolio, plate, summary, figures)
    (OUT_DIR / "IMPACT_PORTFOLIO_V3_REPORT_KR.md").write_text(report)
    (OUT_DIR / "impact_portfolio_v3_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")

    html = make_html(portfolio, plate, summary, figures)
    html_path = copy_to_hub(figures, html)
    summary["html_path"] = str(html_path.relative_to(ROOT))
    (OUT_DIR / "impact_portfolio_v3_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
