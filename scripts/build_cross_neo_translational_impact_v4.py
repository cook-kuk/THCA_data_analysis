#!/usr/bin/env python3
"""Build CROSS-Neo translational impact v4.

v4 adds value-of-information scoring, claim unlock ladders, and a 96-well
assay map on top of the v3 portfolio. This is still claim-boundary scaffolding:
it plans experiments and decision rules, but does not assert clinical utility.
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
V3_DIR = BASE / "impact_portfolio_v3_2026_05_10"
OUT_DIR = BASE / "translational_impact_v4_2026_05_10"
FIG_DIR = OUT_DIR / "figures"

HUB_DIR = ROOT / "project/papers_hub_2026_05_04"
HUB_ASSET_DIR = HUB_DIR / "assets/cross_neo_translational_v4"
LIVE_HUB_DIR = Path("/var/www/papers/papers_hub_2026_05_04")
LIVE_ASSET_DIR = LIVE_HUB_DIR / "assets/cross_neo_translational_v4"

PORTFOLIO = V3_DIR / "impact_candidate_portfolio_v3.tsv"
PLATE_V3 = V3_DIR / "wetlab_plate_v3_impact_design.tsv"
V3_SUMMARY = V3_DIR / "impact_portfolio_v3_summary.json"


LANE_WEIGHTS = {
    "CLEAN_BARNEO_LEAD": 1.00,
    "TCR_MD_STRUCTURAL_LEAD": 1.25,
    "LABEL_NOISE_RESCUE": 0.92,
    "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT": 0.84,
    "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED": 0.70,
    "TCR_STRUCTURE_DISAGREEMENT_AUDIT": 0.78,
    "SUPPORTING_RESERVE": 0.30,
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


def clip01(values):
    return np.clip(values, 0.0, 1.0)


def safe_cols(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [col for col in cols if col in df.columns]


def entropy01(prob: pd.Series | np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(prob, dtype=float), 1e-6, 1 - 1e-6)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))


def lane_prior(df: pd.DataFrame) -> pd.Series:
    lane = text(df, "impact_lane")
    booster = numeric(df, "finetuned_score_booster_prob", 0.0)
    exp_score = numeric(df, "finetuned_experiment_priority_score", 0.0)
    tcr_md = numeric(df, "tcr_md_integrated_score", 0.0)
    md = numeric(df, "md_score", 0.0)
    bma = numeric(df, "bma_v2_discovery_score", 0.0)
    noise = numeric(df, "label_noise_priority_score", 0.0)
    prior = pd.Series(0.35, index=df.index, dtype="float64")
    prior = prior.where(~lane.eq("CLEAN_BARNEO_LEAD"), 0.68 * booster + 0.22 * bma + 0.10 * exp_score)
    prior = prior.where(
        ~lane.eq("TCR_MD_STRUCTURAL_LEAD"),
        0.42 * tcr_md + 0.22 * bma + 0.18 * md + 0.18 * booster,
    )
    prior = prior.where(~lane.eq("LABEL_NOISE_RESCUE"), 0.52 * booster + 0.25 * noise + 0.23 * bma)
    prior = prior.where(
        ~lane.eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
        0.50 * booster + 0.25 * exp_score + 0.25 * noise,
    )
    prior = prior.where(~lane.eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"), 0.70 * booster + 0.30 * bma)
    prior = prior.where(
        ~lane.eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        0.45 * tcr_md + 0.25 * bma + 0.20 * booster + 0.10 * numeric(df, "tcr_structure_disagreement_score", 0.0),
    )
    return pd.Series(clip01(prior), index=df.index)


def claim_unlock_text(lane: str) -> tuple[str, str, str]:
    mapping = {
        "CLEAN_BARNEO_LEAD": (
            "generalizable antigen discovery",
            "mutant pMHC binding plus WT/decoy specificity",
            "elevates clean BAR-Neo lane from computational priority to wetlab-supported antigen discovery",
        ),
        "TCR_MD_STRUCTURAL_LEAD": (
            "mechanistic TCR/MD support",
            "pMHC binding plus tetramer/activation or strong TCR-specific readout",
            "elevates mechanistic lane without making clean benchmark or clinical claim",
        ),
        "LABEL_NOISE_RESCUE": (
            "false-negative rescue / assay mismatch",
            "manual provenance audit plus positive pMHC binding",
            "turns a negative label conflict into reviewer-facing label-noise evidence",
        ),
        "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT": (
            "specificity moat",
            "negative or weak binding under matched assay conditions",
            "shows high-score model can be challenged by hard negatives and controls",
        ),
        "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED": (
            "assay QC positive control",
            "strong positive binding in known/overlap-blocked positive",
            "validates assay dynamic range but cannot support benchmark novelty claim",
        ),
        "TCR_STRUCTURE_DISAGREEMENT_AUDIT": (
            "model-boundary clarification",
            "resolved TCR/MD disagreement by repeat MD, binding, or TCR review",
            "defines where TCR recognition and pMHC presentation disagree",
        ),
    }
    return mapping.get(lane, ("reserve support", "secondary review", "supporting evidence only"))


def build_value_portfolio() -> pd.DataFrame:
    portfolio = read_tsv(PORTFOLIO, required=True)
    df = portfolio.copy()
    lane = text(df, "impact_lane")
    prior = lane_prior(df)
    df["v4_lane_specific_prior"] = prior
    df["v4_prior_entropy"] = entropy01(prior)
    df["v4_lane_weight"] = lane.map(LANE_WEIGHTS).fillna(0.30)
    feasibility = clip01(
        0.34
        + 0.22 * numeric(df, "control_readiness_score", 0.0)
        + 0.18 * numeric(df, "source_hla_novelty_score", 0.0)
        + 0.16 * numeric(df, "validity_dag_cap", 0.0)
        + 0.10 * numeric(df, "stress_guarded_confidence", 0.0)
    )
    feasibility = np.where(lane.eq("TCR_MD_STRUCTURAL_LEAD"), np.maximum(feasibility, 0.70), feasibility)
    feasibility = np.where(lane.eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"), np.maximum(feasibility, 0.62), feasibility)
    df["v4_assay_feasibility_score"] = feasibility
    df["v4_expected_information_gain"] = clip01(
        (0.58 * df["v4_prior_entropy"] + 0.42 * prior)
        * df["v4_lane_weight"]
        * (0.62 + 0.38 * feasibility)
        * (0.78 + 0.22 * numeric(df, "source_hla_novelty_score", 0.0))
    )
    df["v4_claim_unlock_score"] = clip01(
        0.38 * df["impact_portfolio_score"]
        + 0.32 * df["v4_expected_information_gain"]
        + 0.18 * feasibility
        + 0.12 * df["v4_lane_weight"]
    )
    unlocks = [claim_unlock_text(x) for x in lane]
    df["v4_claim_unlock"] = [x[0] for x in unlocks]
    df["v4_unlock_condition"] = [x[1] for x in unlocks]
    df["v4_unlock_impact"] = [x[2] for x in unlocks]
    df["v4_translational_action"] = np.select(
        [
            lane.eq("CLEAN_BARNEO_LEAD"),
            lane.eq("TCR_MD_STRUCTURAL_LEAD"),
            lane.eq("LABEL_NOISE_RESCUE"),
            lane.eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
            lane.eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
            lane.eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        ],
        [
            "plate_v4_clean_discovery_arm",
            "plate_v4_mechanism_arm",
            "plate_v4_label_rescue_arm",
            "plate_v4_specificity_control_arm",
            "plate_v4_qc_positive_control_arm",
            "plate_v4_tcr_boundary_arm",
        ],
        default="reserve",
    )
    front = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "leakage_risk_level",
        "impact_lane",
        "v4_claim_unlock_score",
        "v4_expected_information_gain",
        "v4_lane_specific_prior",
        "v4_prior_entropy",
        "v4_assay_feasibility_score",
        "v4_claim_unlock",
        "v4_unlock_condition",
        "v4_unlock_impact",
        "v4_translational_action",
        "impact_portfolio_score",
        "finetuned_experiment_priority_score",
        "tcr_md_integrated_score",
        "label_noise_priority_score",
        "claim_boundary_v3",
    ]
    front = safe_cols(df, front)
    return df.sort_values(["v4_claim_unlock_score", "impact_portfolio_score"], ascending=False)[
        front + [col for col in df.columns if col not in front]
    ]


def pick_lane(df: pd.DataFrame, lane: str, n: int, exclude: set[str], sort_col: str = "v4_claim_unlock_score") -> pd.DataFrame:
    sub = df[df["impact_lane"].eq(lane) & ~df["candidate_id"].isin(exclude)].copy()
    sub = sub.sort_values(sort_col, ascending=False).head(n)
    exclude.update(sub["candidate_id"].tolist())
    return sub


def build_plate_v4(value_df: pd.DataFrame) -> pd.DataFrame:
    exclude: set[str] = set()
    pieces = [
        pick_lane(value_df, "CLEAN_BARNEO_LEAD", 6, exclude),
        pick_lane(value_df, "TCR_MD_STRUCTURAL_LEAD", 2, exclude),
        pick_lane(value_df, "LABEL_NOISE_RESCUE", 4, exclude),
        pick_lane(value_df, "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT", 4, exclude),
        pick_lane(value_df, "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED", 4, exclude),
        pick_lane(value_df, "TCR_STRUCTURE_DISAGREEMENT_AUDIT", 4, exclude),
    ]
    plate = pd.concat([x for x in pieces if not x.empty], ignore_index=True)
    if len(plate) < 24:
        reserve = value_df[~value_df["candidate_id"].isin(exclude)].head(24 - len(plate)).copy()
        plate = pd.concat([plate, reserve], ignore_index=True)
    plate = plate.head(24).copy()
    plate["plate_v4_slot"] = np.arange(1, len(plate) + 1)
    plate["plate_v4_arm"] = np.select(
        [
            plate["impact_lane"].eq("CLEAN_BARNEO_LEAD"),
            plate["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD"),
            plate["impact_lane"].eq("LABEL_NOISE_RESCUE"),
            plate["impact_lane"].eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
            plate["impact_lane"].eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
            plate["impact_lane"].eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        ],
        [
            "A_clean_discovery",
            "B_mechanism_TCR_MD",
            "C_label_rescue",
            "D_specificity_moat",
            "E_positive_QC_control",
            "F_model_boundary",
        ],
        default="G_reserve",
    )
    plate["plate_v4_primary_readout"] = np.select(
        [
            plate["impact_lane"].eq("TCR_MD_STRUCTURAL_LEAD"),
            plate["impact_lane"].eq("CLEAN_BARNEO_LEAD"),
            plate["impact_lane"].eq("LABEL_NOISE_RESCUE"),
            plate["impact_lane"].eq("HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT"),
            plate["impact_lane"].eq("ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED"),
            plate["impact_lane"].eq("TCR_STRUCTURE_DISAGREEMENT_AUDIT"),
        ],
        [
            "pMHC binding + TCR tetramer/activation",
            "mutant-vs-WT/decoy pMHC binding",
            "pMHC binding rescue after provenance audit",
            "matched negative binding specificity",
            "assay dynamic range/QC",
            "repeat MD/TCR review plus pMHC triage",
        ],
        default="reserve",
    )
    front = [
        "plate_v4_slot",
        "plate_v4_arm",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "label",
        "source_name",
        "leakage_risk_level",
        "impact_lane",
        "v4_claim_unlock_score",
        "v4_expected_information_gain",
        "v4_lane_specific_prior",
        "v4_claim_unlock",
        "plate_v4_primary_readout",
        "v4_unlock_condition",
        "v4_unlock_impact",
        "claim_boundary_v3",
    ]
    front = safe_cols(plate, front)
    return plate[front + [col for col in plate.columns if col not in front]]


def build_claim_ladder(plate: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "ladder_step": 1,
            "claim_layer": "computational score recovery",
            "current_status": "achieved",
            "evidence_now": "score booster source/HLA clean OOF mean AUPRC 0.623, AUROC 0.909",
            "unlock_rule": "already generated and verified",
            "claim_after_unlock": "validated computational prioritization layer",
            "claim_boundary": "benchmark/prioritization only",
        },
        {
            "ladder_step": 2,
            "claim_layer": "clean antigen discovery",
            "current_status": "plate-ready",
            "evidence_now": f"{int((plate['impact_lane'] == 'CLEAN_BARNEO_LEAD').sum())} clean BAR-Neo leads on plate v4",
            "unlock_rule": ">=2 clean leads show mutant pMHC binding with WT/decoy specificity",
            "claim_after_unlock": "wetlab-supported generalizable antigen discovery lane",
            "claim_boundary": "no clinical vaccine-selection claim",
        },
        {
            "ladder_step": 3,
            "claim_layer": "TCR/MD mechanism",
            "current_status": "plate-ready",
            "evidence_now": f"{int((plate['impact_lane'] == 'TCR_MD_STRUCTURAL_LEAD').sum())} TCR/MD structural leads on plate v4",
            "unlock_rule": ">=1 structural lead passes pMHC binding plus TCR tetramer/activation",
            "claim_after_unlock": "orthogonal mechanistic support for candidate prioritization",
            "claim_boundary": "not clean benchmark novelty because source overlap remains",
        },
        {
            "ladder_step": 4,
            "claim_layer": "label-noise rescue",
            "current_status": "audit-ready",
            "evidence_now": f"{int((plate['impact_lane'] == 'LABEL_NOISE_RESCUE').sum())} high-score negative-label rescue cases",
            "unlock_rule": ">=1 rescue case passes provenance audit and pMHC binding",
            "claim_after_unlock": "label-noise/assay-mismatch evidence explaining benchmark ceiling",
            "claim_boundary": "do not relabel without provenance review",
        },
        {
            "ladder_step": 5,
            "claim_layer": "specificity moat",
            "current_status": "control-ready",
            "evidence_now": f"{int((plate['impact_lane'] == 'HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT').sum())} hard negative/specificity cases plus positive controls",
            "unlock_rule": "hard negatives stay weak while positive controls pass",
            "claim_after_unlock": "specificity-aware assay and model-boundary defense",
            "claim_boundary": "assay control evidence, not patient outcome evidence",
        },
        {
            "ladder_step": 6,
            "claim_layer": "clinical/patient relevance",
            "current_status": "locked",
            "evidence_now": "patient metadata gate remains limiting",
            "unlock_rule": "complete patient context plus external/wetlab validation",
            "claim_after_unlock": "patient-level prioritization hypothesis",
            "claim_boundary": "clinical vaccine-selection claim remains unavailable now",
        },
    ]
    return pd.DataFrame(rows)


def condition_set(lane: str) -> list[tuple[str, str, str]]:
    if lane == "TCR_MD_STRUCTURAL_LEAD":
        return [
            ("MUT_pMHC", "binding", "positive supports presentation"),
            ("WT_CONTROL", "binding", "weak WT supports specificity"),
            ("DECOY_CONTROL", "binding", "weak decoy supports specificity"),
            ("TCR_READOUT", "tetramer_or_activation", "positive supports mechanism"),
        ]
    if lane == "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED":
        return [
            ("MUT_pMHC", "binding", "positive confirms assay dynamic range"),
            ("WT_CONTROL", "binding", "specificity context"),
            ("DECOY_CONTROL", "binding", "specificity context"),
            ("QC_REPLICATE", "binding", "replicate QC"),
        ]
    if lane == "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT":
        return [
            ("MUT_pMHC", "binding", "weak result supports specificity"),
            ("WT_CONTROL", "binding", "matched baseline"),
            ("DECOY_CONTROL", "binding", "matched baseline"),
            ("NEG_REPLICATE", "binding", "specificity replicate"),
        ]
    if lane == "TCR_STRUCTURE_DISAGREEMENT_AUDIT":
        return [
            ("MUT_pMHC", "binding", "triages presentation"),
            ("WT_CONTROL", "binding", "specificity context"),
            ("DECOY_CONTROL", "binding", "specificity context"),
            ("TCR_REVIEW", "manual_or_repeat_MD", "resolves disagreement"),
        ]
    if lane == "LABEL_NOISE_RESCUE":
        return [
            ("MUT_pMHC", "binding", "positive supports rescue"),
            ("WT_CONTROL", "binding", "specificity context"),
            ("DECOY_CONTROL", "binding", "specificity context"),
            ("PROVENANCE_REVIEW", "manual_audit", "confirms label/source issue"),
        ]
    return [
        ("MUT_pMHC", "binding", "positive supports discovery"),
        ("WT_CONTROL", "binding", "weak WT supports specificity"),
        ("DECOY_CONTROL", "binding", "weak decoy supports specificity"),
        ("TECH_REPLICATE", "binding", "replicate QC"),
    ]


def build_96well_map(plate: pd.DataFrame) -> pd.DataFrame:
    rows = list("ABCDEFGH")
    wells = [f"{row}{col}" for col in range(1, 13) for row in rows]
    records = []
    well_idx = 0
    for _, row in plate.sort_values("plate_v4_slot").iterrows():
        for condition, readout, expected in condition_set(str(row["impact_lane"])):
            well = wells[well_idx]
            well_idx += 1
            records.append(
                {
                    "well": well,
                    "plate_row": well[0],
                    "plate_col": int(well[1:]),
                    "plate_v4_slot": int(row["plate_v4_slot"]),
                    "plate_v4_arm": row["plate_v4_arm"],
                    "candidate_id": row["candidate_id"],
                    "peptide": row["peptide"],
                    "hla_allele_4digit": row["hla_allele_4digit"],
                    "impact_lane": row["impact_lane"],
                    "assay_condition": condition,
                    "readout": readout,
                    "expected_interpretation": expected,
                    "v4_claim_unlock": row["v4_claim_unlock"],
                    "v4_unlock_condition": row["v4_unlock_condition"],
                    "claim_boundary_v3": row["claim_boundary_v3"],
                }
            )
    return pd.DataFrame(records)


def build_decision_rules() -> pd.DataFrame:
    rows = [
        {
            "decision_node": "assay_qc",
            "go_rule": ">=3/4 positive-control wells positive and blank/decoy controls weak",
            "no_go_rule": "positive controls fail or widespread nonspecific binding",
            "action_if_go": "interpret candidate arms",
            "action_if_no_go": "repeat assay before candidate interpretation",
        },
        {
            "decision_node": "clean_BARNeo_leads",
            "go_rule": ">=2/6 clean leads pass mutant binding and WT/decoy specificity",
            "no_go_rule": "0/6 pass or all fail specificity",
            "action_if_go": "promote clean antigen discovery lane to wetlab-supported",
            "action_if_no_go": "downgrade clean booster score to computational-only",
        },
        {
            "decision_node": "TCR_MD_mechanism",
            "go_rule": ">=1/2 structural leads passes binding plus TCR readout",
            "no_go_rule": "binding or TCR readout fails in both structural leads",
            "action_if_go": "promote orthogonal mechanism support",
            "action_if_no_go": "keep TCR/MD lane diagnostic-only",
        },
        {
            "decision_node": "label_noise_rescue",
            "go_rule": ">=1/4 rescue cases passes provenance audit and binding",
            "no_go_rule": "no rescue case passes after provenance audit",
            "action_if_go": "use as label-noise explanation for benchmark ceiling",
            "action_if_no_go": "remove rescue claim and keep as error-analysis reserve",
        },
        {
            "decision_node": "specificity_moat",
            "go_rule": "hard negatives remain weak while positives pass",
            "no_go_rule": "hard negatives bind strongly at similar rate as positives",
            "action_if_go": "strengthen specificity and control narrative",
            "action_if_no_go": "flag model/assay specificity failure",
        },
    ]
    return pd.DataFrame(rows)


def make_figures(value_df: pd.DataFrame, plate: pd.DataFrame, wellmap: pd.DataFrame, ladder: pd.DataFrame) -> list[Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []
    palette = {
        "CLEAN_BARNEO_LEAD": "#15803d",
        "TCR_MD_STRUCTURAL_LEAD": "#d97706",
        "LABEL_NOISE_RESCUE": "#b45309",
        "HARD_NEGATIVE_OR_FALSE_NEGATIVE_AUDIT": "#7c3aed",
        "ASSAY_POSITIVE_CONTROL_OVERLAP_BLOCKED": "#64748b",
        "TCR_STRUCTURE_DISAGREEMENT_AUDIT": "#dc2626",
        "SUPPORTING_RESERVE": "#94a3b8",
    }

    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    lanes = plate.groupby("impact_lane")["v4_expected_information_gain"].mean().sort_values()
    ax.barh([x.replace("_", " ") for x in lanes.index], lanes.values, color=[palette.get(x, "#94a3b8") for x in lanes.index])
    ax.set_xlabel("Mean expected information gain")
    ax.set_title("Value-of-information by experimental lane")
    for y, val in enumerate(lanes.values):
        ax.text(val + 0.01, y, f"{val:.2f}", va="center")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig1_value_of_information_by_lane.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    y = np.arange(len(ladder))
    status_color = ["#16a34a" if s == "achieved" else "#f59e0b" if s != "locked" else "#64748b" for s in ladder["current_status"]]
    ax.scatter(ladder["ladder_step"], y, s=360, color=status_color, edgecolor="white", linewidth=1.2)
    for _, row in ladder.iterrows():
        ax.text(row["ladder_step"] + 0.08, row["ladder_step"] - 1, row["claim_layer"], va="center", fontsize=9)
    ax.set_yticks([])
    ax.set_xticks(ladder["ladder_step"])
    ax.set_xlabel("Claim ladder step")
    ax.set_title("Claim unlock ladder: computation to wetlab-supported portfolio")
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig2_claim_unlock_ladder.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(10.5, 6.2))
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
    grid = np.zeros((8, 12, 3))
    for _, row in wellmap.iterrows():
        r = ord(row["plate_row"]) - ord("A")
        c = int(row["plate_col"]) - 1
        hex_color = condition_colors.get(row["assay_condition"], "#94a3b8").lstrip("#")
        grid[r, c, :] = [int(hex_color[i : i + 2], 16) / 255 for i in (0, 2, 4)]
    ax.imshow(grid, aspect="equal")
    ax.set_xticks(np.arange(12))
    ax.set_xticklabels(np.arange(1, 13))
    ax.set_yticks(np.arange(8))
    ax.set_yticklabels(list("ABCDEFGH"))
    ax.set_xlabel("Plate column")
    ax.set_ylabel("Plate row")
    ax.set_title("96-well assay map: 24 candidates x 4 conditions")
    for _, row in wellmap.iterrows():
        r = ord(row["plate_row"]) - ord("A")
        c = int(row["plate_col"]) - 1
        ax.text(c, r, str(int(row["plate_v4_slot"])), ha="center", va="center", fontsize=6, color="black")
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=color, label=label.replace("_", " "))
        for label, color in condition_colors.items()
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8)
    fig.tight_layout()
    path = FIG_DIR / "fig3_96well_assay_map.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    figures.append(path)

    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    ax.scatter(
        plate["v4_lane_specific_prior"],
        plate["v4_expected_information_gain"],
        s=160 * plate["v4_claim_unlock_score"],
        c=[palette.get(x, "#94a3b8") for x in plate["impact_lane"]],
        edgecolor="white",
        linewidth=0.7,
        alpha=0.88,
    )
    for _, row in plate.head(10).iterrows():
        ax.text(row["v4_lane_specific_prior"] + 0.008, row["v4_expected_information_gain"] + 0.004, row["candidate_id"], fontsize=7)
    ax.set_xlabel("Lane-specific prior")
    ax.set_ylabel("Expected information gain")
    ax.set_title("Plate v4 candidate value: probability vs information")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = FIG_DIR / "fig4_prior_vs_information_gain.png"
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


def make_report(value_df: pd.DataFrame, plate: pd.DataFrame, wellmap: pd.DataFrame, ladder: pd.DataFrame, rules: pd.DataFrame, summary: dict) -> str:
    return f"""# CROSS-Neo translational impact v4

Generated: {summary["generated_at"]}

## Impact move

v4 converts the portfolio into a **value-of-information wetlab plan**. The key question is no longer just which peptide scores highest; it is which assay result unlocks the strongest claim while preserving overlap, patient-context, and clinical-use boundaries.

## Headline metrics

- Plate v4 candidates: {summary["plate_v4_rows"]}
- 96-well assay map rows: {summary["well_map_rows"]}
- Mean plate claim-unlock score: {summary["mean_plate_claim_unlock_score"]:.3f}
- Mean expected information gain: {summary["mean_plate_expected_information_gain"]:.3f}
- Score recovery vs method-matrix stacker: {summary["score_recovery_vs_method_matrix"]:.2f}x

## Plate v4

{plate[safe_cols(plate, ["plate_v4_slot", "plate_v4_arm", "candidate_id", "peptide", "hla_allele_4digit", "label", "impact_lane", "v4_claim_unlock_score", "v4_expected_information_gain", "v4_claim_unlock", "plate_v4_primary_readout", "claim_boundary_v3"])].to_markdown(index=False, floatfmt=".3f")}

## Claim unlock ladder

{ladder.to_markdown(index=False)}

## Assay decision rules

{rules.to_markdown(index=False)}

## Decision

1. The next impact layer is a 96-well experiment design, not another unrestricted model sprint.
2. The manuscript-safe frame is a controlled translational prioritization portfolio.
3. Clinical vaccine-selection remains locked until patient metadata and external/wetlab evidence are complete.

## Files

- `translational_value_portfolio_v4.tsv`
- `wetlab_plate_v4_value_of_information.tsv`
- `assay_96well_map_v4.tsv`
- `claim_unlock_ladder_v4.tsv`
- `assay_decision_rules_v4.tsv`
- `TRANSLATIONAL_IMPACT_V4_REPORT_KR.md`
"""


def make_html(value_df: pd.DataFrame, plate: pd.DataFrame, wellmap: pd.DataFrame, ladder: pd.DataFrame, rules: pd.DataFrame, summary: dict, figures: list[Path]) -> str:
    asset_map = {fig.name: f"assets/cross_neo_translational_v4/{fig.name}" for fig in figures}
    cards = [
        ("Plate candidates", summary["plate_v4_rows"], "value-of-information selected"),
        ("96-well rows", summary["well_map_rows"], "24 candidates x 4 conditions"),
        ("Mean unlock", summary["mean_plate_claim_unlock_score"], "claim-unlock score"),
        ("Mean EIG", summary["mean_plate_expected_information_gain"], "expected information gain"),
        ("Recovery", summary["score_recovery_vs_method_matrix"], "AUPRC ratio vs failed stacker"),
        ("Locked claim", "clinical", "patient gate remains active"),
    ]
    stat_html = "\n".join(
        f"<div class='stat'><b>{v:.3f}</b><span>{k}</span><small>{n}</small></div>"
        if isinstance(v, float)
        else f"<div class='stat'><b>{v}</b><span>{k}</span><small>{n}</small></div>"
        for k, v, n in cards
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo Translational Impact v4</title>
<style>
:root {{ --bg:#08111f; --panel:#111827; --line:#263244; --text:#e5e7eb; --muted:#94a3b8; --gold:#f5c542; }}
body {{ margin:0; background:var(--bg); color:var(--text); font-family:Inter, system-ui, -apple-system, Segoe UI, sans-serif; }}
header {{ padding:42px 5vw 30px; background:#0f172a; border-bottom:1px solid var(--line); }}
.kicker {{ color:var(--gold); font-size:12px; font-weight:800; letter-spacing:.08em; text-transform:uppercase; }}
h1 {{ font-size:clamp(34px,5vw,66px); line-height:1.02; margin:10px 0 14px; }}
.lead {{ color:#cbd5e1; max-width:980px; font-size:18px; line-height:1.5; }}
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
<div class="kicker">CROSS-Neo · Translational impact v4 · 2026-05-10</div>
<h1>Value-of-information wetlab plan</h1>
<p class="lead">v4 turns the candidate portfolio into a 96-well assay map with explicit claim-unlock rules. The plan raises translational impact while keeping clinical and patient-level claims locked.</p>
<div class="stats">{stat_html}</div>
</header>
<main>
<nav>
<a href="#decision">Decision</a>
<a href="#figures">Figures</a>
<a href="#plate">Plate v4</a>
<a href="#wellmap">96-Well Map</a>
<a href="#ladder">Claim Ladder</a>
<a href="#rules">Decision Rules</a>
</nav>
<div>
<section id="decision" class="panel">
<h2>Decision</h2>
<div class="callout">The impact upgrade is not another score. It is a claim-unlocking experimental design: clean discovery, mechanism, rescue, specificity, assay QC, and model-boundary arms are all on the same 96-well map.</div>
</section>
<section id="figures">
<h2>Figures</h2>
<div class="grid2">
<div class="panel"><img src="{asset_map['fig1_value_of_information_by_lane.png']}" alt="value by lane"></div>
<div class="panel"><img src="{asset_map['fig2_claim_unlock_ladder.png']}" alt="claim ladder"></div>
<div class="panel"><img src="{asset_map['fig3_96well_assay_map.png']}" alt="96 well map"></div>
<div class="panel"><img src="{asset_map['fig4_prior_vs_information_gain.png']}" alt="prior vs information gain"></div>
</div>
</section>
<section id="plate" class="panel">
<h2>Plate v4</h2>
<div class="table-wrap">{html_table(plate, ["plate_v4_slot","plate_v4_arm","candidate_id","peptide","hla_allele_4digit","label","impact_lane","v4_claim_unlock_score","v4_expected_information_gain","v4_claim_unlock","plate_v4_primary_readout","claim_boundary_v3"], 24)}</div>
</section>
<section id="wellmap" class="panel">
<h2>96-Well Map</h2>
<div class="table-wrap">{html_table(wellmap, ["well","plate_v4_slot","candidate_id","peptide","hla_allele_4digit","assay_condition","readout","expected_interpretation","v4_claim_unlock"], 96)}</div>
</section>
<section id="ladder" class="panel">
<h2>Claim Ladder</h2>
<div class="table-wrap">{html_table(ladder, ["ladder_step","claim_layer","current_status","unlock_rule","claim_after_unlock","claim_boundary"], 12)}</div>
</section>
<section id="rules" class="panel">
<h2>Decision Rules</h2>
<div class="table-wrap">{html_table(rules, ["decision_node","go_rule","no_go_rule","action_if_go","action_if_no_go"], 12)}</div>
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
    html_path = HUB_DIR / "cross_neo_translational_impact_v4.html"
    html_path.write_text(html)
    if LIVE_HUB_DIR.exists():
        LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
        for fig in figures:
            shutil.copy2(fig, LIVE_ASSET_DIR / fig.name)
        shutil.copy2(html_path, LIVE_HUB_DIR / html_path.name)
    return html_path


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    value_df = build_value_portfolio()
    plate = build_plate_v4(value_df)
    ladder = build_claim_ladder(plate)
    wellmap = build_96well_map(plate)
    rules = build_decision_rules()
    figures = make_figures(value_df, plate, wellmap, ladder)

    v3_summary = json.loads(V3_SUMMARY.read_text()) if V3_SUMMARY.exists() else {}
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(OUT_DIR.relative_to(ROOT)),
        "plate_v4_rows": int(len(plate)),
        "well_map_rows": int(len(wellmap)),
        "mean_plate_claim_unlock_score": float(plate["v4_claim_unlock_score"].mean()),
        "mean_plate_expected_information_gain": float(plate["v4_expected_information_gain"].mean()),
        "score_recovery_vs_method_matrix": float(v3_summary.get("auprc_recovery_vs_method_matrix", 0.0)),
        "plate_arm_counts": plate["plate_v4_arm"].value_counts().to_dict(),
        "top_claim_unlocks": plate.sort_values("v4_claim_unlock_score", ascending=False).head(8)[
            ["candidate_id", "peptide", "hla_allele_4digit", "v4_claim_unlock_score", "v4_claim_unlock"]
        ].to_dict("records"),
    }

    write_tsv(value_df, "translational_value_portfolio_v4.tsv")
    write_tsv(plate, "wetlab_plate_v4_value_of_information.tsv")
    write_tsv(wellmap, "assay_96well_map_v4.tsv")
    write_tsv(ladder, "claim_unlock_ladder_v4.tsv")
    write_tsv(rules, "assay_decision_rules_v4.tsv")
    (OUT_DIR / "TRANSLATIONAL_IMPACT_V4_REPORT_KR.md").write_text(make_report(value_df, plate, wellmap, ladder, rules, summary))

    html = make_html(value_df, plate, wellmap, ladder, rules, summary, figures)
    html_path = copy_to_hub(figures, html)
    summary["html_path"] = str(html_path.relative_to(ROOT))
    (OUT_DIR / "translational_impact_v4_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
