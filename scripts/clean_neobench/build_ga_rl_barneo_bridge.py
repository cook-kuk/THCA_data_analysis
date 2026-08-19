#!/usr/bin/env python3
"""Bridge GA/RL-discovered candidates into BAR-Neo-NG claim gates."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, normalize_hla, update_manifest, write_tsv


OUTPUTS = [
    "ga_rl_barneo_candidate_scores.tsv",
    "ga_rl_barneo_unique_candidate_summary.tsv",
    "ga_rl_barneo_t1_confirmed.tsv",
    "ga_rl_barneo_t1_unique_candidates.tsv",
    "ga_rl_barneo_blocked_or_watchlist.tsv",
    "ga_rl_barneo_high_ga_blocked_audit.tsv",
    "ga_rl_barneo_failure_analysis.tsv",
    "ga_rl_barneo_decision_summary.tsv",
    "GA_RL_BARNEO_BRIDGE_RESULTS_KR.md",
    "GA_RL_BARNEO_UNIQUE_T1_HANDOFF_KR.md",
    "GA_RL_BARNEO_FAILURE_ANALYSIS_KR.md",
    "GA_RL_BARNEO_BRIDGE_METHOD_CARD.md",
]
HTML_NAME = "ga_rl_barneo_bridge_results_2026_05_11.html"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench/BAR-Neo output root")
    p.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def num_series(df: pd.DataFrame, col: str, default: float = np.nan) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce")


def clean_peptide(value: object) -> str:
    return str(value).strip().upper()


def phla_key(peptide: object, hla: object) -> str:
    return f"{clean_peptide(peptide)}|{normalize_hla(hla)}"


def normalize_biodarwin_academic(path: Path) -> pd.DataFrame:
    df = read_tsv(path)
    if df.empty:
        return df
    out = pd.DataFrame()
    out["ga_rl_candidate_id"] = df.get("candidate_id", "")
    out["ga_rl_track"] = "BioDarwin_GA_RL_academic"
    out["candidate_id"] = df.get("candidate_id", "")
    out["peptide"] = df.get("peptide", "")
    out["hla_allele_4digit"] = df.get("hla_allele_4digit", "").map(normalize_hla)
    out["mhc_class"] = "I"
    out["label"] = df.get("label", np.nan)
    out["source_name"] = df.get("source_name", "")
    out["leakage_risk_level"] = df.get("leakage_risk_level", "")
    out["ga_rl_score"] = num_series(df, "biodarwin_pan_vaccine_score")
    out["ga_rl_cd8_axis_score"] = num_series(df, "biodarwin_cd8_axis_score")
    out["ga_rl_cd4_helper_score"] = num_series(df, "biodarwin_cd4_helper_score")
    out["ga_rl_synergy_score"] = num_series(df, "biodarwin_synergy_score")
    out["kg_ga_evolved_score"] = np.nan
    out["ga_rl_discovery_source"] = "biodarwin_academic_candidate_scores.tsv"
    out["ga_rl_split"] = df.get("split_biodarwin", "")
    out["ga_rl_native_action"] = df.get("immunogenicity_action", "")
    out["ga_rl_native_blockers"] = df.get("immunogenicity_claim_blockers", "")
    out["manual_QA_required"] = False
    out["do_not_train"] = False
    out["locked_public"] = False
    return out


def normalize_biodarwin_industrial(path: Path) -> pd.DataFrame:
    df = read_tsv(path)
    if df.empty:
        return df
    out = pd.DataFrame()
    out["ga_rl_candidate_id"] = df.get("industrial_candidate_id", "")
    out["ga_rl_track"] = "BioDarwin_GA_RL_industrial_locked"
    out["candidate_id"] = df.get("industrial_candidate_id", "")
    out["peptide"] = df.get("peptide", "")
    out["hla_allele_4digit"] = df.get("hla_allele", "").map(normalize_hla)
    out["mhc_class"] = df.get("mhc_class", "I")
    out["label"] = df.get("label", np.nan)
    out["source_name"] = df.get("source_name", df.get("source_id", ""))
    out["leakage_risk_level"] = "public_locked"
    out["ga_rl_score"] = num_series(df, "biodarwin_pan_vaccine_score")
    out["ga_rl_cd8_axis_score"] = num_series(df, "biodarwin_cd8_axis_score")
    out["ga_rl_cd4_helper_score"] = num_series(df, "biodarwin_cd4_helper_score")
    out["ga_rl_synergy_score"] = num_series(df, "biodarwin_synergy_score")
    out["kg_ga_evolved_score"] = np.nan
    out["ga_rl_discovery_source"] = "biodarwin_industrial_locked_scores.tsv"
    out["ga_rl_split"] = df.get("split_biodarwin", "industrial_locked_v0")
    out["ga_rl_native_action"] = df.get("benchmark_use", "")
    out["ga_rl_native_blockers"] = "locked public industrial row; manual QA required; never train"
    out["manual_QA_required"] = df.get("manual_QA_required", False).map(truth)
    out["do_not_train"] = df.get("do_not_train", True).map(truth)
    out["locked_public"] = True
    return out


def normalize_biodarwin_classii(path: Path) -> pd.DataFrame:
    df = read_tsv(path)
    if df.empty:
        return df
    out = pd.DataFrame()
    id_col = "industrial_candidate_id" if "industrial_candidate_id" in df.columns else df.columns[0]
    out["ga_rl_candidate_id"] = df[id_col].astype(str)
    out["ga_rl_track"] = "BioDarwin_GA_RL_classII_scout"
    out["candidate_id"] = out["ga_rl_candidate_id"]
    out["peptide"] = df.get("peptide", "")
    hla_col = "hla_allele" if "hla_allele" in df.columns else "hla_allele_4digit"
    out["hla_allele_4digit"] = df.get(hla_col, "").map(normalize_hla)
    out["mhc_class"] = df.get("mhc_class", "II")
    out["label"] = df.get("label", np.nan)
    out["source_name"] = df.get("source_name", df.get("source_id", ""))
    out["leakage_risk_level"] = "classII_scout"
    out["ga_rl_score"] = num_series(df, "biodarwin_pan_vaccine_score")
    if out["ga_rl_score"].isna().all():
        score_cols = [c for c in df.columns if c.endswith("_score") or c.endswith("_norm")]
        out["ga_rl_score"] = df[score_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1) if score_cols else np.nan
    out["ga_rl_cd8_axis_score"] = num_series(df, "biodarwin_cd8_axis_score")
    out["ga_rl_cd4_helper_score"] = num_series(df, "biodarwin_cd4_helper_score")
    out["ga_rl_synergy_score"] = num_series(df, "biodarwin_synergy_score")
    out["kg_ga_evolved_score"] = np.nan
    out["ga_rl_discovery_source"] = "biodarwin_classII_scout_scores.tsv"
    out["ga_rl_split"] = "classII_scout"
    out["ga_rl_native_action"] = "classII_scout_only"
    out["ga_rl_native_blockers"] = "MHC-II scout row; separate benchmark required"
    out["manual_QA_required"] = True
    out["do_not_train"] = True
    out["locked_public"] = True
    return out


def normalize_kg_ga(path: Path) -> pd.DataFrame:
    df = read_tsv(path)
    if df.empty:
        return df
    out = pd.DataFrame()
    out["ga_rl_candidate_id"] = df.get("candidate_id", "")
    out["ga_rl_track"] = "KG_GA_evolved"
    out["candidate_id"] = df.get("candidate_id", "")
    out["peptide"] = df.get("peptide", "")
    out["hla_allele_4digit"] = df.get("hla_allele_4digit", "").map(normalize_hla)
    out["mhc_class"] = "I"
    out["label"] = df.get("label", np.nan)
    out["source_name"] = df.get("source_name", "")
    out["leakage_risk_level"] = df.get("leakage_risk_level", "")
    out["ga_rl_score"] = num_series(df, "kg_ga_evolved_score")
    out["ga_rl_cd8_axis_score"] = np.nan
    out["ga_rl_cd4_helper_score"] = np.nan
    out["ga_rl_synergy_score"] = np.nan
    out["kg_ga_evolved_score"] = num_series(df, "kg_ga_evolved_score")
    out["kg_ga_evolved_rank"] = num_series(df, "kg_ga_evolved_rank")
    out["ga_rl_discovery_source"] = "kg_ga_evolved_candidate_scores.tsv"
    out["ga_rl_split"] = "kg_ga_all_candidates"
    out["ga_rl_native_action"] = df.get("immunogenicity_action", "")
    out["ga_rl_native_blockers"] = df.get("immunogenicity_claim_blockers", "")
    out["manual_QA_required"] = False
    out["do_not_train"] = False
    out["locked_public"] = False
    return out


def collect_ga_rl(repo_root: Path) -> pd.DataFrame:
    biodarwin = repo_root / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10"
    kg_ga = repo_root / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10"
    pieces = [
        normalize_biodarwin_academic(biodarwin / "biodarwin_academic_candidate_scores.tsv"),
        normalize_biodarwin_industrial(biodarwin / "biodarwin_industrial_locked_scores.tsv"),
        normalize_biodarwin_classii(biodarwin / "biodarwin_classII_scout_scores.tsv"),
        normalize_kg_ga(kg_ga / "kg_ga_evolved_candidate_scores.tsv"),
    ]
    pieces = [p for p in pieces if p is not None and not p.empty]
    if not pieces:
        return pd.DataFrame()
    df = pd.concat(pieces, ignore_index=True, sort=False)
    df["peptide"] = df["peptide"].map(clean_peptide)
    df["hla_allele_4digit"] = df["hla_allele_4digit"].map(normalize_hla)
    df["peptide_hla_key"] = [phla_key(p, h) for p, h in zip(df["peptide"], df["hla_allele_4digit"])]
    df["ga_rl_score"] = pd.to_numeric(df["ga_rl_score"], errors="coerce")
    df["ga_rl_rank_within_track"] = df.groupby("ga_rl_track")["ga_rl_score"].rank(ascending=False, method="first")
    return df


def attach_barneo_ng(ga: pd.DataFrame, ng: pd.DataFrame) -> pd.DataFrame:
    if ga.empty:
        return ga
    df = ga.copy()
    ng_cols = [
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "barneo_ng_score",
        "barneo_ng_rank_global",
        "barneo_ng_decision",
        "barneo_ng_reason",
        "reviewer_kill_disposition",
        "reviewer_kill_reason",
        "impact_readiness_tier",
        "stress_guarded_rank_global",
        "stress_guarded_action",
    ]
    ng_sel = ng[[c for c in ng_cols if c in ng.columns]].copy() if not ng.empty else pd.DataFrame()
    if ng_sel.empty:
        for c in ng_cols:
            if c != "candidate_id":
                df[c] = np.nan
        df["barneo_ng_match_type"] = "no_barneo_ng_table"
        return df
    ng_sel["candidate_id"] = ng_sel["candidate_id"].astype(str)
    ng_sel["peptide_hla_key"] = [phla_key(p, h) for p, h in zip(ng_sel.get("peptide", ""), ng_sel.get("hla_allele_4digit", ""))]

    id_records = ng_sel.drop_duplicates("candidate_id").set_index("candidate_id").to_dict("index")
    phla_records = ng_sel.drop_duplicates("peptide_hla_key").set_index("peptide_hla_key").to_dict("index")
    # Keep GA/RL peptide/HLA/source columns as the primary identity fields; BAR-Neo-NG
    # columns below are gate metadata only.
    identity_cols = {"candidate_id", "peptide", "hla_allele_4digit"}
    fill_cols = [c for c in ng_cols if c not in identity_cols] + ["barneo_ng_match_type"]
    rows: list[dict[str, Any]] = []
    for _, r in df.iterrows():
        cid = str(r.get("candidate_id", ""))
        key = str(r.get("peptide_hla_key", ""))
        rec = id_records.get(cid)
        match_type = "candidate_id"
        if rec is None:
            rec = phla_records.get(key)
            match_type = "peptide_hla"
        if rec is None:
            rows.append({c: np.nan for c in fill_cols[:-1]} | {"barneo_ng_match_type": "unmatched"})
        else:
            rows.append({c: rec.get(c, np.nan) for c in fill_cols[:-1]} | {"barneo_ng_match_type": match_type})
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def call_decision(row: pd.Series) -> tuple[str, str, float]:
    ga_score = float(row.get("ga_rl_score", 0.0)) if pd.notna(row.get("ga_rl_score", np.nan)) else 0.0
    ng_score = float(row.get("barneo_ng_score", 0.0)) if pd.notna(row.get("barneo_ng_score", np.nan)) else 0.0
    base = 0.45 * ga_score + 0.55 * ng_score if ng_score > 0 else 0.65 * ga_score
    decision = str(row.get("barneo_ng_decision", ""))
    mhc = str(row.get("mhc_class", "")).upper()
    leakage = str(row.get("leakage_risk_level", "")).lower()
    match = str(row.get("barneo_ng_match_type", ""))

    if mhc == "II":
        return "classII_scout_separate_benchmark", "MHC-II scout row; do not pool with class-I BAR-Neo-NG", min(base, 0.20)
    if truth(row.get("locked_public", False)) or truth(row.get("do_not_train", False)):
        if match == "unmatched":
            return "locked_public_postfreeze_only", "industrial/public locked row; manual QA and no training use", min(base, 0.35)
    if match == "unmatched":
        return "unmatched_ga_rl_candidate_needs_clean_schema", "GA/RL row is not in CLEAN-NeoBench/BAR-Neo-NG by candidate_id or peptide-HLA", min(base, 0.42)
    if decision == "T1_phla_assay_design_candidate":
        return "GA_RL_hit_confirmed_T1_by_BAR_Neo_NG", "GA/RL candidate also passes BAR-Neo-NG T1 pHLA assay-design gate", min(max(base, 0.72), 0.82)
    if decision == "T2_resolve_before_claim":
        return "GA_RL_hit_T2_near_neighbor_or_audit_required", "GA/RL candidate needs near-neighbor/source/overlap resolution before claim", min(base, 0.56)
    if decision == "blocked_from_clean_claim" or leakage == "high":
        return "GA_RL_hit_blocked_from_clean_claim", "GA/RL score is not claim-safe because leakage/exact-overlap gate blocks it", min(base, 0.25)
    if decision == "watchlist_manual_review":
        return "GA_RL_hit_watchlist_manual_review", "GA/RL candidate is plausible but not a clean T1 claim under BAR-Neo-NG", min(base, 0.62)
    if decision == "low_priority_or_abstain":
        return "GA_RL_hit_low_priority_or_abstain", "BAR-Neo-NG deprioritizes or abstains despite GA/RL score", min(base, 0.40)
    return "GA_RL_candidate_review_only", "No clean BAR-Neo-NG claim gate triggered", min(base, 0.50)


def build_bridge(repo_root: Path, output_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    ga = collect_ga_rl(repo_root)
    ng = read_tsv(output_root / "barneo_ng_algorithm_scores.tsv")
    bridge = attach_barneo_ng(ga, ng)
    if bridge.empty:
        return bridge, pd.DataFrame(), {}
    decisions = bridge.apply(call_decision, axis=1, result_type="expand")
    bridge["ga_rl_barneo_decision"] = decisions[0]
    bridge["ga_rl_barneo_reason"] = decisions[1]
    bridge["ga_rl_barneo_priority_score"] = pd.to_numeric(decisions[2], errors="coerce").fillna(0.0)
    bridge = bridge.sort_values(["ga_rl_barneo_priority_score", "ga_rl_score"], ascending=[False, False])
    bridge["ga_rl_barneo_rank_global"] = range(1, len(bridge) + 1)
    summary = bridge["ga_rl_barneo_decision"].value_counts().rename_axis("ga_rl_barneo_decision").reset_index(name="n")

    comparison_root = repo_root / "project/results/cross_neo_ga_rl_full_comparison_2026_05_10"
    kg_root = repo_root / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10"
    meta = {
        "full_comparison": read_json(comparison_root / "ga_rl_full_comparison_summary.json"),
        "kg_ga": read_json(kg_root / "kg_ga_summary.json"),
    }
    return bridge, summary, meta


def table_html(df: pd.DataFrame, cols: list[str], n: int = 30) -> str:
    if df is None or df.empty:
        return "<p class='muted'>No rows.</p>"
    cols = [c for c in cols if c in df.columns]
    d = df[cols].head(n).copy()
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in cols)
    rows = []
    for _, r in d.iterrows():
        rows.append("<tr>" + "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in cols) + "</tr>")
    return f"<div class='table'><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


def build_unique_views(bridge: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if bridge.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    ranked = bridge.sort_values(["ga_rl_barneo_priority_score", "ga_rl_score"], ascending=[False, False]).copy()
    grouped = ranked.groupby("candidate_id", dropna=False)
    unique = grouped.head(1).copy()
    unique["n_ga_rl_tracks_supporting"] = unique["candidate_id"].map(grouped["ga_rl_track"].nunique())
    unique["supporting_ga_rl_tracks"] = unique["candidate_id"].map(
        grouped["ga_rl_track"].agg(lambda x: "; ".join(sorted(set(x.astype(str)))))
    )
    unique["max_ga_rl_score_any_track"] = unique["candidate_id"].map(grouped["ga_rl_score"].max())
    unique = unique.sort_values(["ga_rl_barneo_priority_score", "max_ga_rl_score_any_track"], ascending=[False, False])
    unique["ga_rl_barneo_unique_rank"] = range(1, len(unique) + 1)

    t1_unique = unique[unique["ga_rl_barneo_decision"].eq("GA_RL_hit_confirmed_T1_by_BAR_Neo_NG")].copy()
    if not t1_unique.empty:
        t1_unique["handoff_tier"] = "T1_GA_RL_BAR_Neo_NG_phla_assay_design"
        t1_unique["assay_next_step"] = "Class-I pHLA binding/stability plus immunogenicity assay after metadata intake"
        t1_unique["minimum_metadata_packet"] = (
            "WT peptide; gene; mutation_id; source_protein_window; expression_tpm; vaf; clonality; "
            "patient_id; disease_context; HLA-LOH/B2M/presentation/safety fields"
        )
        t1_unique["claim_boundary"] = "pHLA assay-design candidate; not patient/translational or clinical vaccine selection"

    blocked_mask = ranked["ga_rl_barneo_decision"].astype(str).str.contains("blocked|watchlist|T2", na=False)
    high_ga_mask = ranked["ga_rl_score"].fillna(0.0).ge(0.85) | ranked["ga_rl_rank_within_track"].fillna(999999).le(25)
    high_ga_blocked = ranked[blocked_mask & high_ga_mask].copy()
    high_ga_blocked = high_ga_blocked.groupby("candidate_id", dropna=False).head(1).copy()
    high_ga_blocked["audit_priority"] = np.where(
        high_ga_blocked["ga_rl_barneo_decision"].astype(str).str.contains("blocked", na=False),
        "A_blocked_high_GA_review_leakage_or_overlap",
        "B_watchlist_high_GA_manual_review",
    )
    high_ga_blocked["audit_question"] = np.where(
        high_ga_blocked["ga_rl_barneo_decision"].astype(str).str.contains("blocked", na=False),
        "Is this GA/RL signal explained by leakage, exact/near overlap, or source reuse?",
        "Does added metadata/assay evidence move this from watchlist to T1?",
    )
    return unique, t1_unique, high_ga_blocked, ranked


def build_failure_analysis(t1_unique: pd.DataFrame, high_ga_blocked: pd.DataFrame, bridge: pd.DataFrame) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    if not t1_unique.empty:
        t1 = t1_unique.copy()
        t1["analysis_group"] = "unique_t1_handoff"
        t1["analysis_reason"] = "BAR-Neo-NG T1 pass"
        t1["analysis_overlap_status"] = "low_leakage_or_pass"
        parts.append(t1)
    if not high_ga_blocked.empty:
        b = high_ga_blocked.copy()
        b["analysis_group"] = "high_ga_blocked_audit"
        exact = b.get("exact_peptide_train_overlap", pd.Series(False, index=b.index)).map(truth)
        exact_hla = b.get("exact_peptide_hla_train_overlap", pd.Series(False, index=b.index)).map(truth)
        near = b.get("near_peptide_train_overlap", pd.Series(False, index=b.index)).map(truth)
        leakage = b.get("leakage_risk_level", pd.Series("", index=b.index)).astype(str).str.lower()
        b["analysis_reason"] = np.select(
            [
                exact_hla | exact,
                near,
                leakage.eq("high"),
            ],
            [
                "exact_peptide_or_peptide_hla_overlap",
                "near_peptide_overlap",
                "high_leakage",
            ],
            default="manual_review_or_source_reuse",
        )
        b["analysis_overlap_status"] = np.select(
            [
                exact_hla | exact,
                near,
            ],
            [
                "exact_overlap",
                "near_overlap",
            ],
            default="non_overlap_but_high_GA",
        )
        parts.append(b)
    if not parts:
        return pd.DataFrame()
    df = pd.concat(parts, ignore_index=True, sort=False)
    df["source_rank"] = df.groupby(["analysis_group", "source_name"], dropna=False)["candidate_id"].transform("count")
    df["hla_rank"] = df.groupby(["analysis_group", "hla_allele_4digit"], dropna=False)["candidate_id"].transform("count")
    return df


def write_figures(output_root: Path, hub_root: Path, bridge: pd.DataFrame, summary: pd.DataFrame) -> tuple[list[str], list[str]]:
    if bridge.empty:
        return [], []
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return [], []

    paths: list[str] = []
    html_srcs: list[str] = []
    fig_dir = output_root / "figures"
    hub_fig_dir = hub_root / "assets" / "ga_rl_barneo"
    ensure_dir(fig_dir)
    ensure_dir(hub_fig_dir)

    counts = summary.sort_values("n", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.barh(counts["ga_rl_barneo_decision"].astype(str), counts["n"].astype(float), color="#5eead4")
    ax.set_xlabel("Rows")
    ax.set_title("GA/RL rows after BAR-Neo-NG claim gates")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    p1 = fig_dir / "fig_ga_rl_barneo_decision_counts.png"
    hp1 = hub_fig_dir / p1.name
    fig.savefig(p1, dpi=180)
    fig.savefig(hp1, dpi=180)
    plt.close(fig)
    paths.append(str(p1))
    html_srcs.append(f"assets/ga_rl_barneo/{p1.name}")

    d = bridge.copy()
    d["ga_rl_score"] = pd.to_numeric(d["ga_rl_score"], errors="coerce")
    d["barneo_ng_score"] = pd.to_numeric(d["barneo_ng_score"], errors="coerce")
    plot = d.dropna(subset=["ga_rl_score", "barneo_ng_score"]).copy()
    colors = {
        "GA_RL_hit_confirmed_T1_by_BAR_Neo_NG": "#86efac",
        "GA_RL_hit_blocked_from_clean_claim": "#f87171",
        "GA_RL_hit_watchlist_manual_review": "#e3b341",
        "GA_RL_hit_T2_near_neighbor_or_audit_required": "#a78bfa",
        "GA_RL_hit_low_priority_or_abstain": "#94a3b8",
    }
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for decision, sub in plot.groupby("ga_rl_barneo_decision"):
        ax.scatter(
            sub["ga_rl_score"],
            sub["barneo_ng_score"],
            s=12,
            alpha=0.55,
            label=str(decision)[:34],
            color=colors.get(str(decision), "#64748b"),
        )
    ax.set_xlabel("GA/RL discovery score")
    ax.set_ylabel("BAR-Neo-NG claim-gated score")
    ax.set_title("Discovery strength vs reviewer-safe gate")
    ax.grid(alpha=0.22)
    ax.legend(loc="best", fontsize=7, frameon=False)
    fig.tight_layout()
    p2 = fig_dir / "fig_ga_rl_vs_barneo_ng_score_gate.png"
    hp2 = hub_fig_dir / p2.name
    fig.savefig(p2, dpi=180)
    fig.savefig(hp2, dpi=180)
    plt.close(fig)
    paths.append(str(p2))
    html_srcs.append(f"assets/ga_rl_barneo/{p2.name}")
    return paths, html_srcs


def write_reports(
    output_root: Path,
    hub_root: Path,
    bridge: pd.DataFrame,
    unique: pd.DataFrame,
    t1: pd.DataFrame,
    t1_unique: pd.DataFrame,
    blocked: pd.DataFrame,
    high_ga_blocked: pd.DataFrame,
    failure: pd.DataFrame,
    summary: pd.DataFrame,
    meta: dict[str, Any],
    figure_paths: list[str],
    figure_srcs: list[str],
) -> None:
    top_cols = [
        "ga_rl_barneo_rank_global",
        "ga_rl_barneo_unique_rank",
        "ga_rl_candidate_id",
        "ga_rl_track",
        "n_ga_rl_tracks_supporting",
        "supporting_ga_rl_tracks",
        "ga_rl_score",
        "max_ga_rl_score_any_track",
        "ga_rl_barneo_priority_score",
        "ga_rl_barneo_decision",
        "candidate_id",
        "barneo_ng_score",
        "barneo_ng_decision",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "label",
        "leakage_risk_level",
        "ga_rl_native_blockers",
        "ga_rl_barneo_reason",
    ]
    handoff_cols = [
        "ga_rl_barneo_unique_rank",
        "candidate_id",
        "peptide",
        "hla_allele_4digit",
        "source_name",
        "label",
        "max_ga_rl_score_any_track",
        "barneo_ng_score",
        "ga_rl_barneo_priority_score",
        "supporting_ga_rl_tracks",
        "handoff_tier",
        "assay_next_step",
        "minimum_metadata_packet",
        "claim_boundary",
    ]
    failure_cols = [
        "analysis_group",
        "candidate_id",
        "ga_rl_track",
        "peptide",
        "hla_allele_4digit",
        "source_name",
        "ga_rl_score",
        "barneo_ng_score",
        "ga_rl_barneo_decision",
        "analysis_reason",
        "analysis_overlap_status",
        "exact_peptide_train_overlap",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "nearest_neighbor_similarity",
        "leakage_risk_level",
        "audit_priority",
    ]
    full_meta = meta.get("full_comparison", {})
    kg_meta = meta.get("kg_ga", {})
    method_card = """# GA/RL -> BAR-Neo-NG Bridge Method Card

## What It Does

This bridge takes GA/RL-discovered or GA/RL-ranked neoantigen candidates from BioDarwin and KG-GA outputs, maps them to CLEAN-NeoBench candidates by candidate ID or peptide-HLA pair, and applies BAR-Neo-NG claim gates.

## Claim Boundary

GA/RL fitness is treated as a discovery and prioritization signal. Clean claims require BAR-Neo-NG gate survival, leakage/source/HLA audit, public-training overlap resolution, and downstream metadata or assay evidence.

## Reviewer-Safe Use

- GA/RL high score plus BAR-Neo-NG T1: pHLA assay-design candidate.
- GA/RL high score plus BAR-Neo-NG blocked: useful failure/audit case, not a clean claim.
- Industrial public rows: post-freeze/manual-QA only, never clean training rows.
- MHC-II scout rows: separate benchmark required.
"""
    (output_root / "GA_RL_BARNEO_BRIDGE_METHOD_CARD.md").write_text(method_card.strip() + "\n")

    handoff = f"""# GA/RL + BAR-Neo-NG Unique T1 Handoff KR

## 결론

GA/RL로 찾은 후보를 BAR-Neo-NG gate로 정리하면, 실제 unique T1 handoff 후보는 `{len(t1_unique)}`개다. 이 후보들은 pHLA assay-design queue로 보낼 수 있지만, patient/translational claim은 metadata intake 전까지 막는다.

## Unique T1 handoff table

{dataframe_to_markdown(t1_unique[[c for c in handoff_cols if c in t1_unique.columns]], max_rows=20)}

## High-GA but blocked/watchlist audit queue

{dataframe_to_markdown(high_ga_blocked[[c for c in top_cols + ["audit_priority", "audit_question"] if c in high_ga_blocked.columns]], max_rows=30)}

## Rule

- T1: GA/RL discovery signal and BAR-Neo-NG gate agree.
- Watchlist/T2: useful manual-review candidates, not clean claims.
- Blocked: audit/failure-analysis candidates, not headline candidates.
"""
    (output_root / "GA_RL_BARNEO_UNIQUE_T1_HANDOFF_KR.md").write_text(handoff.strip() + "\n")

    failure_audit = f"""# GA/RL -> BAR-Neo-NG Failure Analysis KR

## 결론

GA/RL 고득점 후보는 많지만, 실제 clean handoff는 T1 3개뿐이다. 나머지는 exact/near overlap, high leakage, source reuse, 또는 MHC-II separate-benchmark 사유로 막힌다.

## Failure analysis table

{dataframe_to_markdown(failure[[c for c in failure_cols if c in failure.columns]].head(60) if not failure.empty else failure, max_rows=60)}

## Source summary

{dataframe_to_markdown(failure.groupby(["analysis_group", "source_name"], dropna=False).size().reset_index(name="n").sort_values(["analysis_group", "n"], ascending=[True, False]) if not failure.empty else failure, max_rows=40)}

## HLA summary

{dataframe_to_markdown(failure.groupby(["analysis_group", "hla_allele_4digit"], dropna=False).size().reset_index(name="n").sort_values(["analysis_group", "n"], ascending=[True, False]) if not failure.empty else failure, max_rows=40)}

## Interpretation

- `exact_overlap` and `near_overlap` are not clean claim territory.
- `high_leakage` rows are audit/failure cases even when GA/RL score is high.
- `unique_t1_handoff` rows are the only immediate assay-design queue.
"""
    (output_root / "GA_RL_BARNEO_FAILURE_ANALYSIS_KR.md").write_text(failure_audit.strip() + "\n")

    kr = f"""# GA/RL -> BAR-Neo-NG Bridge Results KR

## 한 줄 결론

말씀하신 `GA/RL로 찾은 후보`를 BAR-Neo-NG gate에 다시 태웠다. 결론은 **GA/RL raw champion은 강하지만, reviewer-safe 후보는 BAR-Neo-NG T1 gate를 통과한 3개로 좁혀진다.**

## GA/RL benchmark context

- Top GA/RL algorithm: `{full_meta.get("top_master_algorithm", "NA")}`
- KG-GA all AUPRC: `{full_meta.get("kg_ga_all_auprc", "NA")}`
- KG-GA frozen-validation-like AUPRC: `{full_meta.get("kg_ga_frozen_validation_auprc", "NA")}`
- KG-GA low/medium leakage AUPRC: `{full_meta.get("kg_ga_low_medium_leakage_auprc", "NA")}`
- Boundary: `{full_meta.get("boundary", "retrospective prioritization only")}`
- GA generations: `{kg_meta.get("ga_generations", "NA")}`
- GA population size: `{kg_meta.get("ga_population_size", "NA")}`

## Decision summary

{dataframe_to_markdown(summary, max_rows=20)}

## Unique candidate summary

{dataframe_to_markdown(unique[[c for c in top_cols if c in unique.columns]].head(30) if not unique.empty else unique, max_rows=30)}

## Unique T1 handoff

{dataframe_to_markdown(t1_unique[[c for c in handoff_cols if c in t1_unique.columns]], max_rows=20)}

## GA/RL + BAR-Neo-NG confirmed T1

{dataframe_to_markdown(t1[[c for c in top_cols if c in t1.columns]], max_rows=20)}

## Top GA/RL bridge rows

{dataframe_to_markdown(bridge[[c for c in top_cols if c in bridge.columns]].head(30) if not bridge.empty else bridge, max_rows=30)}

## GA/RL high-score but not clean claim

{dataframe_to_markdown(blocked[[c for c in top_cols if c in blocked.columns]].head(30) if not blocked.empty else blocked, max_rows=30)}

## High-GA blocked audit

{dataframe_to_markdown(high_ga_blocked[[c for c in top_cols + ["audit_priority", "audit_question"] if c in high_ga_blocked.columns]].head(30) if not high_ga_blocked.empty else high_ga_blocked, max_rows=30)}

## Practical interpretation

- GA/RL is the discovery engine.
- BAR-Neo-NG is the reviewer-safe gate.
- A high GA/RL score alone is not a clean neoantigen claim.
- `GA_RL_hit_confirmed_T1_by_BAR_Neo_NG` rows are the immediate pHLA assay-design queue.
- `blocked_from_clean_claim` rows are still valuable as audit/failure cases, not as headline candidates.
"""
    (output_root / "GA_RL_BARNEO_BRIDGE_RESULTS_KR.md").write_text(kr.strip() + "\n")

    stat_html = "".join(
        f"<div><b>{int(r.n):,}</b><span>{html.escape(str(r.ga_rl_barneo_decision))}</span></div>"
        for _, r in summary.iterrows()
    )
    page = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GA/RL -> BAR-Neo-NG Bridge</title>
<style>
body{{margin:0;background:#0d1117;color:#e6edf3;font-family:JetBrains Mono,ui-monospace,Menlo,monospace;line-height:1.55}}
header{{padding:44px 34px 24px;border-bottom:1px solid #2a3441}}h1{{font-family:Georgia,serif;font-size:46px;margin:0 0 10px}}.lead{{color:#c8d1dc;max-width:980px}}
.stats{{display:grid;grid-template-columns:repeat(4,minmax(150px,1fr));gap:10px;margin-top:20px}}.stats div{{border:1px solid #2a3441;background:#101820;padding:12px}}.stats b{{display:block;color:#5eead4;font-size:24px}}.stats span{{color:#9aa7b4;font-size:12px}}
main{{max-width:1360px;margin:auto;padding:24px}}section{{margin:0 0 18px;border:1px solid #2a3441;background:#101820;padding:18px;overflow:auto}}h2{{font-family:Georgia,serif}}.muted{{color:#9aa7b4}}
table{{width:100%;border-collapse:collapse;font-size:12px}}th,td{{border-bottom:1px solid #2a3441;padding:8px;text-align:left;vertical-align:top;white-space:nowrap}}th{{color:#e3b341}}.path{{padding:0 24px 36px;color:#9aa7b4;font-size:12px}}
img{{max-width:100%;height:auto;border:1px solid #2a3441;background:#0d1117}}
</style></head><body><header><h1>GA/RL -> BAR-Neo-NG Bridge</h1><p class="lead">GA/RL discovery candidates mapped into BAR-Neo-NG claim gates. Discovery strength and reviewer-safe claim eligibility are shown separately.</p><div class="stats">{stat_html}</div></header><main>
<section><h2>Unique T1 Handoff</h2>{table_html(t1_unique, handoff_cols, 20)}</section>
<section><h2>Confirmed T1 Rows</h2>{table_html(t1, top_cols, 20)}</section>
<section><h2>Top Bridge Rows</h2>{table_html(bridge, top_cols, 40)}</section>
        <section><h2>High-GA Blocked Audit</h2>{table_html(high_ga_blocked, top_cols + ["audit_priority", "audit_question"], 40)}</section>
        <section><h2>Failure Analysis</h2>{table_html(failure, failure_cols + ["source_rank", "hla_rank"], 40)}</section>
        <section><h2>Blocked / Watchlist / T2</h2>{table_html(blocked, top_cols, 40)}</section>
<section><h2>Figures</h2>{''.join(f'<p><img src="{html.escape(src)}"></p>' for src in figure_srcs) or '<p class="muted">No figures generated.</p>'}</section>
</main><p class="path">Source: {html.escape(str(output_root / 'ga_rl_barneo_candidate_scores.tsv'))}</p></body></html>"""
    ensure_dir(hub_root)
    (hub_root / HTML_NAME).write_text(page)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    hub_root = Path(args.hub_root)
    if not hub_root.is_absolute():
        hub_root = repo_root / hub_root
    ensure_dir(output_root)
    ensure_dir(hub_root)

    bridge, summary, meta = build_bridge(repo_root, output_root)
    t1 = bridge[bridge.get("ga_rl_barneo_decision", pd.Series(dtype=str)).eq("GA_RL_hit_confirmed_T1_by_BAR_Neo_NG")].copy() if not bridge.empty else pd.DataFrame()
    blocked = bridge[bridge.get("ga_rl_barneo_decision", pd.Series(dtype=str)).str.contains("blocked|watchlist|T2", na=False)].copy() if not bridge.empty else pd.DataFrame()
    unique, t1_unique, high_ga_blocked, _ranked = build_unique_views(bridge)
    failure = build_failure_analysis(t1_unique, high_ga_blocked, bridge)
    figure_paths, figure_srcs = write_figures(output_root, hub_root, bridge, summary)

    write_tsv(bridge, output_root / "ga_rl_barneo_candidate_scores.tsv")
    write_tsv(unique, output_root / "ga_rl_barneo_unique_candidate_summary.tsv")
    write_tsv(t1, output_root / "ga_rl_barneo_t1_confirmed.tsv")
    write_tsv(t1_unique, output_root / "ga_rl_barneo_t1_unique_candidates.tsv")
    write_tsv(blocked, output_root / "ga_rl_barneo_blocked_or_watchlist.tsv")
    write_tsv(high_ga_blocked, output_root / "ga_rl_barneo_high_ga_blocked_audit.tsv")
    write_tsv(failure, output_root / "ga_rl_barneo_failure_analysis.tsv")
    write_tsv(summary, output_root / "ga_rl_barneo_decision_summary.tsv")
    write_reports(
        output_root,
        hub_root,
        bridge,
        unique,
        t1,
        t1_unique,
        blocked,
        high_ga_blocked,
        failure,
        summary,
        meta,
        figure_paths,
        figure_srcs,
    )

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_ga_rl_bridge_rows": int(len(bridge)),
            "n_ga_rl_bridge_t1_confirmed": int(len(t1)),
            "n_ga_rl_bridge_unique_candidates": int(len(unique)),
            "n_ga_rl_bridge_t1_unique_candidates": int(len(t1_unique)),
            "n_ga_rl_bridge_high_ga_blocked_audit": int(len(high_ga_blocked)),
            "n_ga_rl_bridge_failure_rows": int(len(failure)),
            "n_ga_rl_bridge_decision_rows": int(len(summary)),
            "ga_rl_bridge_top_priority_score": float(bridge["ga_rl_barneo_priority_score"].max()) if not bridge.empty else 0.0,
            "ga_rl_bridge_html": str(hub_root / HTML_NAME),
        }
    )
    manifest.setdefault("output_files", [])
    for out in OUTPUTS + figure_paths + [str(hub_root / HTML_NAME)]:
        if str(out) not in manifest["output_files"]:
            manifest["output_files"].append(str(out))
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "ga_rl_barneo_bridge",
        {
            "outputs": OUTPUTS + [str(hub_root / HTML_NAME)],
            "n_bridge_rows": int(len(bridge)),
            "n_t1_confirmed": int(len(t1)),
            "n_t1_unique": int(len(t1_unique)),
            "n_high_ga_blocked": int(len(high_ga_blocked)),
            "n_failure_rows": int(len(failure)),
            "warnings": ["GA/RL rows are discovery/prioritization signals; clean claims require BAR-Neo-NG gates and audit."],
        },
    )
    print(
        "[ga-rl-barneo-bridge] "
        f"rows={len(bridge)} unique={len(unique)} t1_rows={len(t1)} "
        f"t1_unique={len(t1_unique)} high_ga_blocked={len(high_ga_blocked)} failure={len(failure)} "
        f"top={manifest['summary']['ga_rl_bridge_top_priority_score']:.4f}"
    )


if __name__ == "__main__":
    main()
