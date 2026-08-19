#!/usr/bin/env python3
"""Write a reviewer-safe CLEAN-NeoBench + BAR-Neo-BMA HTML dossier."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd

from common import update_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    parser.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub directory")
    parser.add_argument("--page-name", default="clean_neobench_barneo_dossier_2026_05_09.html")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def fmt_num(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        x = float(value)
        if abs(x) >= 1000:
            return f"{x:,.0f}"
        return f"{x:.{digits}f}"
    except Exception:
        return html.escape(str(value))


def table_html(df: pd.DataFrame, columns: list[str], limit: int = 20) -> str:
    if df.empty:
        return "<p class=\"muted\">No rows available.</p>"
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return "<p class=\"muted\">Requested columns unavailable.</p>"
    view = df.loc[:, cols].head(limit).copy()
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    rows = []
    for _, row in view.iterrows():
        cells = []
        for c in cols:
            value = row[c]
            if isinstance(value, float):
                value = fmt_num(value)
            cells.append(f"<td>{html.escape(str(value))}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f"<div class=\"table-wrap\"><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


def stat(label: str, value: object, note: str = "") -> str:
    return (
        "<div class=\"stat\">"
        f"<div class=\"stat-value\">{html.escape(str(value))}</div>"
        f"<div class=\"stat-label\">{html.escape(label)}</div>"
        f"<div class=\"stat-note\">{html.escape(note)}</div>"
        "</div>"
    )


def build_context_summary(selector: pd.DataFrame) -> pd.DataFrame:
    if selector.empty:
        return pd.DataFrame()
    rows = []
    for ctx, group in selector.groupby("selector_context", dropna=False):
        ordered = group.sort_values("selector_weight", ascending=False)
        rows.append(
            {
                "selector_context": ctx,
                "n_selected": len(ordered),
                "top_methods": ", ".join(ordered["method_name"].astype(str).head(7)),
                "top_roles": ", ".join(ordered["method_role"].astype(str).head(7)),
            }
        )
    return pd.DataFrame(rows)


def build_abstention_summary(bma: pd.DataFrame) -> pd.DataFrame:
    if bma.empty or "bma_abstention_reason_primary" not in bma.columns:
        return pd.DataFrame()
    counts = (
        bma["bma_abstention_reason_primary"]
        .fillna("")
        .replace("", "No abstention")
        .value_counts()
    )
    return counts.rename_axis("reason").reset_index(name="n")


def build_review_tier_summary(failure_aware: pd.DataFrame) -> pd.DataFrame:
    if failure_aware.empty or "failure_aware_review_tier" not in failure_aware.columns:
        return pd.DataFrame()
    counts = (
        failure_aware["failure_aware_review_tier"]
        .fillna("")
        .replace("", "no_review_tier")
        .value_counts()
    )
    return counts.rename_axis("review_tier").reset_index(name="n")


def build_split_summary(split_metrics: pd.DataFrame) -> pd.DataFrame:
    if split_metrics.empty or "split_contract" not in split_metrics.columns:
        return pd.DataFrame()
    agg = (
        split_metrics.groupby("split_contract")
        .agg(
            n_metric_rows=("method_name", "size"),
            n_methods=("method_name", "nunique"),
            median_AUPRC=("AUPRC", "median"),
            median_top10_precision=("top10_precision", "median"),
        )
        .reset_index()
        .sort_values("n_metric_rows", ascending=False)
    )
    return agg


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    hub_root = Path(args.hub_root)
    if not hub_root.is_absolute():
        hub_root = repo_root / hub_root
    hub_root.mkdir(parents=True, exist_ok=True)

    leaderboard = read_tsv(output_root / "clean_neobench_leaderboard.tsv")
    split_metrics = read_tsv(output_root / "clean_neobench_split_metrics.tsv")
    winloss_summary = read_tsv(output_root / "clean_neobench_winloss_method_summary.tsv")
    split_winners = read_tsv(output_root / "clean_neobench_split_winner_board.tsv")
    stress_methods = read_tsv(output_root / "clean_neobench_source_hla_stress_method_summary.tsv")
    stress_slices = read_tsv(output_root / "clean_neobench_source_hla_stress_slice_board.tsv")
    korean_hla_board = read_tsv(output_root / "clean_neobench_korean_hla_method_board.tsv")
    bma = read_tsv(output_root / "barneo_bma_candidate_scores.tsv")
    weights = read_tsv(output_root / "barneo_bma_method_weights.tsv")
    selector = read_tsv(output_root / "barneo_bma_selector_audit.tsv")
    stress_guarded = read_tsv(output_root / "barneo_stress_guarded_candidate_scores.tsv")
    stress_guarded_weights = read_tsv(output_root / "barneo_stress_guarded_method_weights.tsv")
    high_impact_leads = read_tsv(output_root / "barneo_high_impact_lead_candidates.tsv")
    high_impact_topk = read_tsv(output_root / "barneo_high_impact_topk_audit.tsv")
    high_impact_ablation = read_tsv(output_root / "barneo_high_impact_ablation_audit.tsv")
    reviewer_kill = read_tsv(output_root / "barneo_high_impact_reviewer_kill_audit.tsv")
    t1_readiness = read_tsv(output_root / "barneo_t1_translational_readiness.tsv")
    t1_impact_assay = read_tsv(output_root / "barneo_t1_assay_design_matrix.tsv")
    ng_scores = read_tsv(output_root / "barneo_ng_algorithm_scores.tsv")
    ng_t1 = read_tsv(output_root / "barneo_ng_t1_candidates.tsv")
    ng_summary = read_tsv(output_root / "barneo_ng_decision_summary.tsv")
    ga_rl_bridge = read_tsv(output_root / "ga_rl_barneo_candidate_scores.tsv")
    ga_rl_t1 = read_tsv(output_root / "ga_rl_barneo_t1_confirmed.tsv")
    ga_rl_t1_unique = read_tsv(output_root / "ga_rl_barneo_t1_unique_candidates.tsv")
    ga_rl_high_blocked = read_tsv(output_root / "ga_rl_barneo_high_ga_blocked_audit.tsv")
    ga_rl_failure = read_tsv(output_root / "ga_rl_barneo_failure_analysis.tsv")
    ga_rl_patient_gate = read_tsv(output_root / "ga_rl_patient_gate_priority.tsv")
    ga_rl_patient_gate_summary = read_tsv(output_root / "ga_rl_patient_gate_priority_summary.tsv")
    ga_rl_heatmap = read_tsv(output_root / "ga_rl_t1_failure_heatmap.tsv")
    ga_rl_neighbors = read_tsv(output_root / "ga_rl_t1_nearest_neighbors.tsv")
    ga_rl_distribution = read_tsv(output_root / "ga_rl_t1_distribution_summary.tsv")
    ga_rl_atlas = read_tsv(output_root / "ga_rl_failure_atlas.tsv")
    ga_rl_atlas_summary = read_tsv(output_root / "ga_rl_failure_atlas_summary.tsv")
    ga_rl_atlas_reason = read_tsv(output_root / "ga_rl_failure_atlas_reason_summary.tsv")
    ga_rl_rescue = read_tsv(output_root / "ga_rl_rescue_quadrant.tsv")
    ga_rl_rescue_summary = read_tsv(output_root / "ga_rl_rescue_quadrant_summary.tsv")
    ga_rl_rescue_reason = read_tsv(output_root / "ga_rl_rescue_quadrant_reason_summary.tsv")
    ga_rl_rescue_evidence = read_tsv(output_root / "ga_rl_rescue_evidence.tsv")
    ga_rl_rescue_evidence_summary = read_tsv(output_root / "ga_rl_rescue_evidence_summary.tsv")
    ga_rl_summary = read_tsv(output_root / "ga_rl_barneo_decision_summary.tsv")
    nature_claims = read_tsv(output_root / "clean_neobench_nature_claim_ladder.tsv")
    nature_gaps = read_tsv(output_root / "clean_neobench_nature_validation_gap_table.tsv")
    public_audit = read_tsv(output_root / "clean_neobench_public_tool_overlap_audit.tsv")
    public_intake = read_tsv(output_root / "clean_neobench_public_overlap_audit_intake.tsv")
    public_training_row_summary = read_tsv(output_root / "clean_neobench_public_training_row_overlap_summary.tsv")
    public_training_requirements = read_tsv(output_root / "clean_neobench_public_training_corpus_requirements.tsv")
    failure_aware = read_tsv(output_root / "barneo_failure_aware_candidate_scores.tsv")
    manual_queue = read_tsv(output_root / "barneo_manual_review_queue.tsv")
    challenge_queue = read_tsv(output_root / "barneo_distribution_challenge_queue.tsv")
    distribution_vulnerability = read_tsv(output_root / "clean_neobench_method_distribution_vulnerability.tsv")
    shift_error_link = read_tsv(output_root / "clean_neobench_train_distribution_error_link.tsv")
    remediation_plan = read_tsv(output_root / "clean_neobench_distribution_remediation_plan.tsv")
    contextual_scores = read_tsv(output_root / "barneo_contextual_bma_candidate_scores.tsv")
    contextual_weights = read_tsv(output_root / "barneo_contextual_bma_method_weights.tsv")
    contextual_summary = read_tsv(output_root / "barneo_contextual_bma_context_summary.tsv")
    challenge_pack = read_tsv(output_root / "clean_neobench_challenge_pack.tsv")
    challenge_summary = read_tsv(output_root / "clean_neobench_challenge_axis_summary.tsv")
    experiment_plan = read_tsv(output_root / "clean_neobench_next_experiment_plan.tsv")
    patient_requirements = read_tsv(output_root / "patient_gated_clean_neo_metadata_requirements.tsv")
    patient_scenarios = read_tsv(output_root / "patient_gated_clean_neo_demo_scenarios.tsv")
    patient_queue = read_tsv(output_root / "patient_gated_clean_neo_candidate_queue.tsv")
    master = read_tsv(output_root / "clean_neobench_master.tsv")
    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}

    n_candidates = len(master) if not master.empty else len(bma)
    n_methods = int(leaderboard["method_name"].nunique()) if "method_name" in leaderboard.columns else 0
    n_metric_rows = len(split_metrics)
    n_bma_abstain = int(bma["bma_abstain"].sum()) if "bma_abstain" in bma.columns else 0
    n_bma_nonabstain = len(bma) - n_bma_abstain if not bma.empty else 0
    n_stress_claim_manual = int((stress_guarded.get("stress_guarded_action", pd.Series(dtype=str)) == "claim_safe_candidate_after_manual_audit").sum()) if not stress_guarded.empty else 0
    n_reviewer_kill_pass = int((reviewer_kill.get("reviewer_kill_disposition", pd.Series(dtype=str)) == "passes_current_reviewer_kill_audit").sum()) if not reviewer_kill.empty else 0
    n_t1_assay_ready = int(t1_readiness.get("peptide_hla_assay_design_ready", pd.Series(dtype=bool)).astype(bool).sum()) if not t1_readiness.empty else 0
    n_ng_t1 = int((ng_scores.get("barneo_ng_decision", pd.Series(dtype=str)) == "T1_phla_assay_design_candidate").sum()) if not ng_scores.empty else 0
    n_ga_rl_t1_unique = len(ga_rl_t1_unique) if not ga_rl_t1_unique.empty else int(ga_rl_t1.get("candidate_id", pd.Series(dtype=str)).nunique()) if not ga_rl_t1.empty else 0
    n_ga_rl_failure = len(ga_rl_failure) if not ga_rl_failure.empty else 0
    n_ga_rl_patient_gate = len(ga_rl_patient_gate) if not ga_rl_patient_gate.empty else 0
    n_ga_rl_patient_gate_candidates = len(ga_rl_patient_gate_summary) if not ga_rl_patient_gate_summary.empty else 0
    high_conf = int((bma.get("bma_confidence_bin", pd.Series(dtype=str)) == "high").sum()) if not bma.empty else 0
    selected_methods = bma.get("selected_methods", pd.Series("", index=bma.index)).fillna("").astype(str) if not bma.empty else pd.Series(dtype=str)
    fallback_mask = (
        selected_methods.str.startswith("BAR_Neo_fallback")
        | bma.get("selected_method_count", pd.Series(1, index=bma.index)).fillna(1).eq(0)
        if not bma.empty
        else pd.Series(dtype=bool)
    )
    fallback_rows = int(fallback_mask.sum()) if not bma.empty else 0
    failure_nonabstain = int((~failure_aware["failure_aware_abstain"].astype(bool)).sum()) if not failure_aware.empty and "failure_aware_abstain" in failure_aware.columns else 0

    context_summary = build_context_summary(selector)
    abstention_summary = build_abstention_summary(bma)
    review_tier_summary = build_review_tier_summary(failure_aware)
    split_summary = build_split_summary(split_metrics)
    if not patient_requirements.empty and "metadata_status" in patient_requirements.columns:
        patient_blocked = patient_requirements[
            patient_requirements["metadata_status"].isin(["absent_from_master", "present_but_empty"])
        ].copy()
    else:
        patient_blocked = pd.DataFrame()
    top_bma = bma.sort_values("bma_rank_global") if "bma_rank_global" in bma.columns else bma
    fallback_min_rank = int(top_bma.loc[fallback_mask, "bma_rank_global"].min()) if fallback_rows and "bma_rank_global" in top_bma.columns else "NA"
    clean_review = bma[(~bma.get("bma_abstain", pd.Series([True] * len(bma))))].copy() if not bma.empty else pd.DataFrame()
    if not clean_review.empty and "patient_gated_bma_score" in clean_review.columns:
        clean_review = clean_review.sort_values("patient_gated_bma_score", ascending=False)

    missing_metadata = manifest.get("missing_metadata", {})
    missing_text = ", ".join(f"{k}: {v}" for k, v in sorted(missing_metadata.items())) if isinstance(missing_metadata, dict) else "Recorded in manifest."

    css = """
    :root { --bg:#0d1117; --panel:#151b23; --panel2:#101820; --ink:#e6edf3; --muted:#9aa7b4; --line:#2a3441; --gold:#e3b341; --cyan:#5eead4; --red:#f87171; --green:#86efac; }
    * { box-sizing:border-box; }
    body { margin:0; background:var(--bg); color:var(--ink); font-family:"JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, monospace; line-height:1.55; }
    a { color:var(--cyan); text-decoration:none; }
    .hero { min-height:72vh; display:flex; align-items:flex-end; padding:56px 48px 38px; background:linear-gradient(160deg, rgba(94,234,212,.10), transparent 38%), radial-gradient(circle at 82% 12%, rgba(227,179,65,.13), transparent 30%), #0d1117; border-bottom:1px solid var(--line); }
    .hero-inner { max-width:1180px; width:100%; }
    .kicker { color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:13px; font-weight:700; }
    h1 { font-family:Georgia, "Times New Roman", serif; font-size:clamp(42px, 6vw, 78px); line-height:.98; margin:14px 0 18px; letter-spacing:0; max-width:1000px; }
    .lead { max-width:980px; font-size:18px; color:#c8d1dc; }
    .hero-links { display:flex; flex-wrap:wrap; gap:10px; margin-top:18px; }
    .hero-links a { border:1px solid var(--line); padding:8px 11px; background:#101820; color:#d7fff6; }
    .stats { display:grid; grid-template-columns:repeat(6, minmax(130px,1fr)); gap:12px; margin-top:28px; }
    .stat { border:1px solid var(--line); background:rgba(21,27,35,.72); padding:14px; min-height:104px; }
    .stat-value { font-size:24px; color:white; font-weight:800; }
    .stat-label { color:var(--cyan); font-size:12px; margin-top:4px; }
    .stat-note { color:var(--muted); font-size:11px; margin-top:8px; }
    .layout { display:grid; grid-template-columns:260px minmax(0,1fr); gap:30px; max-width:1400px; margin:0 auto; padding:34px 28px 80px; }
    nav { position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:18px 0; }
    nav a { display:block; padding:8px 10px; color:#b9c4d0; border-left:2px solid transparent; font-size:13px; }
    nav a:hover { border-left-color:var(--cyan); color:white; }
    section { padding:22px 0 34px; border-bottom:1px solid var(--line); }
    h2 { font-family:Georgia, "Times New Roman", serif; font-size:32px; margin:0 0 14px; letter-spacing:0; }
    h3 { font-size:18px; margin:24px 0 10px; color:var(--gold); }
    .num { color:var(--gold); margin-right:10px; font-size:18px; }
    .grid { display:grid; grid-template-columns:repeat(2, minmax(0,1fr)); gap:14px; }
    .box { border:1px solid var(--line); background:var(--panel); padding:16px; }
    .box strong { color:white; }
    .muted { color:var(--muted); }
    .warn { color:var(--red); font-weight:700; }
    .ok { color:var(--green); font-weight:700; }
    .gold { color:var(--gold); font-weight:700; }
    .table-wrap { overflow:auto; border:1px solid var(--line); background:var(--panel2); margin:14px 0; }
    table { width:100%; border-collapse:collapse; font-size:12px; }
    th, td { padding:9px 10px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; white-space:nowrap; }
    th { color:var(--gold); background:#10151c; position:sticky; top:0; }
    td { color:#d6dee7; }
    .cards { display:grid; grid-template-columns:repeat(3, minmax(0,1fr)); gap:14px; }
    .card { border:1px solid var(--line); background:var(--panel); padding:16px; }
    .path { color:#b9c4d0; font-size:12px; overflow-wrap:anywhere; }
    @media (max-width: 920px) { .layout { grid-template-columns:1fr; padding:24px 16px 64px; } nav { position:relative; max-height:none; } .stats, .cards, .grid { grid-template-columns:1fr; } .hero { padding:42px 20px 28px; } }
    """

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CLEAN-NeoBench + BAR-Neo-BMA Dossier</title>
  <style>{css}</style>
</head>
<body>
  <header class="hero">
    <div class="hero-inner">
      <div class="kicker">Leakage-aware neoantigen benchmarking · 2026-05-09</div>
      <h1>CLEAN-NeoBench + BAR-Neo-BMA</h1>
      <p class="lead">Reviewer-safe benchmark contract, benchmark-adaptive reliability ranking, and practical agentic ensemble weighting for neoantigen candidate triage. This page deliberately avoids SOTA, clinical vaccine-selection, public-clean-baseline, and QK headline claims.</p>
      <div class="hero-links">
        <a href="barneo_high_impact_candidate_pack_2026_05_10.html">BAR-Neo high-impact candidate pack</a>
        <a href="barneo_top_pass_reviewer_evidence_2026_05_10.html">Top-pass reviewer evidence</a>
        <a href="barneo_t1_translational_readiness_2026_05_10.html">T1 translational readiness</a>
        <a href="barneo_t1_impact_upgrade_plan_2026_05_10.html">T1 impact upgrade plan</a>
        <a href="barneo_ng_algorithm_results_2026_05_10.html">BAR-Neo-NG algorithm results</a>
        <a href="ga_rl_barneo_bridge_results_2026_05_11.html">GA/RL bridge results</a>
        <a href="barneo_nature_grade_evolution_2026_05_10.html">Nature-grade evolution</a>
        <a href="cancer_vaccine_barneo_x_interpretability_2026_05_09.html">BAR-Neo-X claim-safe XAI layer</a>
        <a href="cancer_vaccine_full_dossier.html">Cancer vaccine full dossier</a>
        <a href="index.html">Hub index</a>
      </div>
      <div class="stats">
        {stat("Candidates", f"{n_candidates:,}", "Unified benchmark rows")}
        {stat("Methods", f"{n_methods:,}", "Internal + caveated comparator outputs")}
        {stat("Split metrics", f"{n_metric_rows:,}", "Exact, near, source, HLA, low prevalence")}
        {stat("BMA abstain", f"{n_bma_abstain:,}", "Reviewer-safe conservative mode")}
        {stat("BMA review rows", f"{n_bma_nonabstain:,}", "Non-abstain candidate prompts")}
        {stat("Stress claim-audit", f"{n_stress_claim_manual:,}", "Stress-guarded manual-audit rows")}
        {stat("Kill-audit pass", f"{n_reviewer_kill_pass:,}", "Automated high-impact reviewer checks")}
        {stat("T1 assay-ready", f"{n_t1_assay_ready:,}", "pHLA design-ready, metadata-blocked")}
        {stat("BAR-Neo-NG T1", f"{n_ng_t1:,}", "Claim-gated assay-design candidates")}
        {stat("GA/RL T1 unique", f"{n_ga_rl_t1_unique:,}", "GA/RL hits surviving BAR-Neo-NG")}
        {stat("GA/RL failure rows", f"{n_ga_rl_failure:,}", "Blocked / overlap / leakage audit")}
        {stat("GA/RL patient gate", f"{n_ga_rl_patient_gate:,}", "Patient-gated research triage rows")}
        {stat("GA/RL patient candidates", f"{n_ga_rl_patient_gate_candidates:,}", "Unique patient-gated candidates")}
      </div>
    </div>
  </header>
  <div class="layout">
    <nav>
      <a href="#tldr">01 TL;DR</a>
      <a href="#position">02 Algorithm Position</a>
      <a href="#leaderboard">03 Leaderboard</a>
      <a href="#bma">04 BMA Weights</a>
      <a href="#selector">05 Selector Audit</a>
      <a href="#candidates">06 Candidate Scores</a>
      <a href="#splits">07 Split Robustness</a>
      <a href="#failure-aware">08 Failure-Aware</a>
      <a href="#distribution-error">09 Distribution Error</a>
      <a href="#contextual-bma">10 Contextual BMA</a>
      <a href="#challenge-pack">11 Challenge Pack</a>
      <a href="#manual-queue">12 Manual Queue</a>
      <a href="#public-audit">13 Public Audit</a>
      <a href="#patient-gated">14 Patient Gates</a>
      <a href="#caveats">15 Caveats</a>
      <a href="#visuals">16 Visuals</a>
      <a href="#paths">17 Paths</a>
    </nav>
    <main>
      <section id="tldr">
        <h2><span class="num">01</span>TL;DR</h2>
        <div class="grid">
          <div class="box"><strong>What is built:</strong> CLEAN-NeoBench normalizes candidate/method/split tables; BAR-Neo adds reliability score and abstention; BAR-Neo-BMA turns benchmark behavior into posterior expert weights.</div>
          <div class="box"><strong>What is not claimed:</strong> no clinical vaccine selection, no new SOTA predictor, no external validation proven, no QK headline claim, no clean public baseline without overlap audit.</div>
          <div class="box"><strong>Current practical use:</strong> use BMA scores as research review prompts. High score + abstention means inspect manually, not claim validity.</div>
          <div class="box"><strong>BAR-Neo-X upgrade:</strong> the separate claim-safe XAI layer explains each row and downranks leakage-heavy top rows for reviewer-facing candidate triage.</div>
          <div class="box"><strong>Stress-guarded upgrade:</strong> source/HLA stress behavior is now converted into candidate-specific method weights, so fragile methods are downweighted in the contexts where they fail.</div>
          <div class="box"><strong>Fallback control:</strong> fallback-only rows are capped and start at rank {fallback_min_rank}, so missing expert support cannot dominate.</div>
        </div>
      </section>

      <section id="position">
        <h2><span class="num">02</span>Algorithm Position</h2>
        <p><span class="gold">Our safe position:</span> benchmark-adaptive reliability controller, not a standalone black-box SOTA predictor.</p>
        {table_html(pd.DataFrame([
            {"Layer": "Structure_LR", "Role": "local anchor", "Allowed claim": "honest internal anchor"},
            {"Layer": "RF / PU / source-balanced models", "Role": "internal candidates", "Allowed claim": "benchmarked internal comparators"},
            {"Layer": "W7A / W7B / stacked models", "Role": "fusion candidates", "Allowed claim": "internal candidate ensemble"},
            {"Layer": "ESM2 Bayesian", "Role": "uncertainty branch", "Allowed claim": "uncertainty/OOD support"},
            {"Layer": "QK branch", "Role": "bounded fallback", "Allowed claim": "bounded fallback/fusion component"},
            {"Layer": "BAR-Neo", "Role": "reliability ranker", "Allowed claim": "calibrated prioritization with abstention"},
            {"Layer": "BAR-Neo-BMA", "Role": "agentic ensemble controller", "Allowed claim": "Bayesian-style benchmark-derived expert weighting"},
            {"Layer": "BAR-Neo-X", "Role": "post-hoc XAI + claim-safe reranker", "Allowed claim": "leakage-aware reviewer-facing triage, not SOTA"},
            {"Layer": "BAR-Neo-NG", "Role": "claim-gated algorithm layer", "Allowed claim": "research triage / pHLA assay-design ranking with abstention"},
        ]), ["Layer", "Role", "Allowed claim"], 20)}
      </section>

      <section id="leaderboard">
        <h2><span class="num">03</span>CLEAN-NeoBench Leaderboard</h2>
        <p class="muted">Overall table is reviewer-safe: public pretrained tools remain caveated unless overlap audit exists.</p>
        {table_html(leaderboard, ["method_name", "method_role", "method_family", "mean_AUPRC", "mean_top10_precision", "mean_ECE", "reviewer_safe_score"], 15)}
        <h3>Win/Loss vs Structure_LR</h3>
        <p class="muted">This table answers the practical question: which methods beat the honest local anchor on matched split slices, and where source/HLA stress still breaks them.</p>
        {table_html(winloss_summary, ["method_name", "method_role", "method_disposition", "n_matched_anchor_splits", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"], 20)}
        <h3>Split Winner Board</h3>
        {table_html(split_winners, ["split_contract", "split_group", "overall_winner", "overall_winner_role", "overall_winner_AUPRC", "clean_internal_winner", "clean_internal_winner_AUPRC", "caveated_public_winner", "public_winner_clean_claim_allowed"], 35)}
      </section>

      <section id="bma">
        <h2><span class="num">04</span>BAR-Neo-BMA Posterior Weights</h2>
        <p>BMA weights combine AUPRC, top-k behavior, calibration, reviewer-safe score, source-collapse penalty, method role priors, and public-overlap/QK/uncertainty penalties.</p>
        {table_html(weights, ["method_name", "method_role", "method_family", "mean_AUPRC", "mean_top10_precision", "mean_ECE", "source_collapse_rate", "utility_score", "posterior_weight_global"], 20)}
      </section>

      <section id="selector">
        <h2><span class="num">05</span>QUBO-Style Selector Audit</h2>
        <p>The selector is QUBO-style in the practical engineering sense: maximize expert utility and family diversity while penalizing redundancy and unsafe over-reliance. It is not a QK headline claim.</p>
        {table_html(context_summary, ["selector_context", "n_selected", "top_methods", "top_roles"], 10)}
        {table_html(selector.sort_values(["selector_context", "selector_weight"], ascending=[True, False]) if not selector.empty else selector, ["selector_context", "method_name", "method_role", "method_family", "selector_weight", "selector_objective", "redundancy_penalty"], 35)}
      </section>

      <section id="candidates">
        <h2><span class="num">06</span>Candidate Scores</h2>
        <h3>Stress-guarded BAR-Neo review queue</h3>
        <p class="muted">This is the practical upgrade layer: source/HLA/leakage stress audits control candidate-specific method weights before ranking.</p>
        {table_html(stress_guarded, ["candidate_id", "stress_guarded_final_review_score", "stress_guarded_claim_safe_score", "stress_guarded_discovery_score", "stress_guarded_confidence", "stress_guarded_action", "stress_guarded_reason_primary", "source_name", "hla_allele_4digit", "peptide", "label", "stress_guarded_clean_method_count", "public_weight_fraction", "fallback_weight_fraction"], 35)}
        <h3>Stress-guarded method weights</h3>
        {table_html(stress_guarded_weights, ["method_name", "method_role", "method_family", "mean_AUPRC", "reviewer_safe_score", "source_heldout_median_delta_AUPRC", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "stress_guarded_base_weight", "recommended_use"], 25)}
        <h3>BAR-Neo-NG algorithm layer</h3>
        <p class="muted">BAR-Neo-NG is the new claim-gated algorithm layer. It combines stress-guarded discovery, clean-claim score, confidence, method agreement, HLA support, assay readiness, translational metadata, and explicit penalties for leakage, public dependency, fallback dependency, near neighbors, and MHC-II class mismatch.</p>
        {table_html(ng_summary, ["barneo_ng_decision", "n"], 20)}
        {table_html(ng_scores, ["candidate_id", "barneo_ng_rank_global", "barneo_ng_score", "barneo_ng_decision", "source_name", "hla_allele_4digit", "peptide", "label", "leakage_risk_level", "reviewer_kill_disposition", "impact_readiness_tier", "barneo_ng_reason"], 25)}
        <h3>BAR-Neo-NG T1 candidates</h3>
        {table_html(ng_t1, ["candidate_id", "barneo_ng_rank_global", "barneo_ng_score", "barneo_ng_decision", "source_name", "hla_allele_4digit", "peptide", "label", "barneo_ng_reason"], 12)}
        <h3>GA/RL -> BAR-Neo-NG bridge</h3>
        <p class="muted">This is the correction for the GA/RL-discovered candidate pool: GA/RL score is treated as discovery strength, then BAR-Neo-NG decides whether the row is T1 assay-design, watchlist, blocked, class-II separate benchmark, or locked public post-freeze only.</p>
        {table_html(ga_rl_summary, ["ga_rl_barneo_decision", "n"], 20)}
        {table_html(ga_rl_t1_unique, ["ga_rl_barneo_unique_rank", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "max_ga_rl_score_any_track", "barneo_ng_score", "ga_rl_barneo_priority_score", "supporting_ga_rl_tracks", "handoff_tier", "assay_next_step"], 12)}
        {table_html(ga_rl_t1, ["ga_rl_barneo_rank_global", "candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_priority_score", "ga_rl_barneo_decision", "ga_rl_barneo_reason"], 12)}
        <h3>GA/RL high-score blocked audit</h3>
        {table_html(ga_rl_high_blocked, ["candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_decision", "leakage_risk_level", "audit_priority", "audit_question"], 20)}
        <h3>GA/RL failure analysis</h3>
        {table_html(ga_rl_failure, ["analysis_group", "candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "source_name", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_decision", "analysis_reason", "analysis_overlap_status"], 20)}
        <h3>GA/RL patient gate</h3>
        {table_html(ga_rl_patient_gate_summary, ["patient_gate_unique_rank", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "patient_gate_best_score", "patient_gate_best_weighted_priority", "patient_gate_best_scenario_name", "patient_gate_best_disease", "patient_gate_clinical_use"], 10)}
        {table_html(ga_rl_patient_gate, ["patient_gate_priority_rank", "candidate_id", "scenario_id", "scenario_name", "disease", "research_priority", "demo_patient_gated_score", "patient_gate_weighted_priority", "patient_gated_rank_within_scenario", "patient_gated_confidence_bin", "hard_gate_status", "abstain", "abstention_reason_primary"], 20)}
        <h3>GA/RL failure heatmap</h3>
        <p class="muted">This pack separates exact / near / public overlap, held-out source and HLA contracts, and nearest-neighbor distribution context for the three unique GA/RL T1 candidates.</p>
        {table_html(ga_rl_distribution, ["candidate_id", "peptide", "hla_allele_4digit", "source_name", "label", "ga_rl_barneo_priority_score", "barneo_ng_score", "leakage_risk_level", "split_source_heldout", "split_hla_heldout", "source_positive_prevalence", "source_support_count", "hla_positive_prevalence", "hla_support_count", "top_nn_candidate_id", "top_nn_peptide", "top_nn_similarity", "top_nn_label"], 10)}
        {table_html(ga_rl_heatmap, ["candidate_id", "axis", "status", "value", "detail"], 30)}
        {table_html(ga_rl_neighbors, ["query_candidate_id", "candidate_id", "peptide", "source_name", "hla_allele_4digit", "label", "similarity", "same_source", "same_hla", "same_supertype"], 18)}
        <h3>GA/RL failure atlas</h3>
        <p class="muted">This is the higher-impact board: 3 clean T1 survivors plus 72 high-GA blocked rows on the same audit surface.</p>
        {table_html(ga_rl_atlas_summary, ["atlas_group", "n_candidates", "n_pos", "mean_priority_score", "mean_barneo_ng_score", "n_high_leakage", "mean_source_prev", "mean_hla_prev", "mean_top_nn_similarity"], 10)}
        {table_html(ga_rl_atlas_reason, ["atlas_group", "atlas_decision", "audit_priority", "n_candidates", "mean_priority_score", "mean_barneo_ng_score", "mean_source_prev", "mean_hla_prev"], 20)}
        {table_html(ga_rl_atlas, ["atlas_order", "atlas_group", "candidate_id", "candidate_tier", "atlas_decision", "peptide", "hla_allele_4digit", "source_name", "ga_rl_barneo_priority_score", "barneo_ng_score", "leakage_risk_level", "audit_priority", "audit_question", "top_nn_candidate_id", "top_nn_similarity"], 25)}
        <h3>GA/RL rescue quadrant</h3>
        <p class="muted">This is the action queue: 3 T1 survivors, 3 recoverable watchlist candidates, and the rest as hard blocks.</p>
        {table_html(ga_rl_rescue_summary, ["rescue_category", "n_candidates", "n_pos", "mean_priority_score", "mean_rescue_priority_score", "mean_top_nn_similarity"], 10)}
        {table_html(ga_rl_rescue_reason, ["rescue_category", "rescue_reason_primary", "n_candidates", "mean_priority_score", "mean_top_nn_similarity"], 20)}
        {table_html(ga_rl_rescue, ["rescue_order", "rescue_category", "candidate_id", "peptide", "hla_allele_4digit", "source_name", "label", "ga_rl_barneo_priority_score", "ga_rl_barneo_decision", "rescue_action", "rescue_reason_primary", "top_nn_candidate_id", "top_nn_similarity"], 25)}
        <h3>GA/RL rescue evidence</h3>
        <p class="muted">The recoverable watchlist rows are still research triage only. This table is the missing-metadata intake packet that would be needed to consider promoting them.</p>
        {table_html(ga_rl_rescue_evidence_summary, ["evidence_grade", "n_candidates", "mean_missing", "mean_priority", "mean_source_prev", "mean_hla_prev"], 10)}
        {table_html(ga_rl_rescue_evidence, ["candidate_id", "peptide", "hla_allele_4digit", "source_name", "missing_metadata_count", "present_metadata_count", "evidence_grade", "missing_metadata", "present_metadata", "evidence_action"], 10)}
        {table_html(ga_rl_bridge, ["ga_rl_barneo_rank_global", "candidate_id", "ga_rl_track", "peptide", "hla_allele_4digit", "ga_rl_score", "barneo_ng_score", "ga_rl_barneo_priority_score", "ga_rl_barneo_decision", "leakage_risk_level", "ga_rl_native_blockers"], 25)}
        <h3>High-impact lead evidence pack</h3>
        <p class="muted">Lead candidates are checked for benchmark-label top-k behavior, public/fallback dependence, and clean internal method support. This remains manual-audit research triage.</p>
        {table_html(high_impact_topk, ["slice", "n", "n_pos", "precision", "mean_score", "mean_confidence", "n_high_leakage", "mean_public_weight_fraction", "mean_fallback_weight_fraction", "claim_boundary"], 14)}
        {table_html(high_impact_leads, ["candidate_id", "manual_audit_priority", "stress_guarded_rank_global", "stress_guarded_final_review_score", "benchmark_label_status", "source_name", "hla_allele_4digit", "peptide", "public_qk_independence_note"], 18)}
        <h3>High-impact public/fallback ablation proxy</h3>
        {table_html(high_impact_ablation, ["candidate_id", "manual_audit_priority", "clean_only_proxy_score", "stress_guarded_clean_method_count", "public_weight_fraction", "fallback_weight_fraction", "top_internal_methods", "ablation_interpretation"], 12)}
        <h3>Reviewer kill audit</h3>
        <p class="muted">This table separates high-impact leads that survive automated leakage, near-neighbor, public/fallback dependence, and same-source negative-pressure checks from rows that require manual audit or are blocked from clean claims.</p>
        {table_html(reviewer_kill, ["candidate_id", "manual_audit_priority", "stress_guarded_rank_global", "stress_guarded_final_review_score", "source_name", "hla_allele_4digit", "peptide", "leakage_risk_level", "nearest_neighbor_similarity", "source_negatives_at_or_above", "public_weight_fraction", "fallback_weight_fraction", "reviewer_kill_disposition", "reviewer_kill_reason"], 20)}
        <h3>T1 translational readiness</h3>
        <p class="muted">T1 rows are pHLA assay-design/manual-review candidates. Patient/translational claims remain blocked until WT/gene/mutation/expression/clonality/patient/safety metadata are linked.</p>
        {table_html(t1_readiness, ["candidate_id", "stress_guarded_rank_global", "stress_guarded_final_review_score", "peptide", "hla_allele_4digit", "benchmark_evidence_score", "peptide_hla_assay_design_ready", "translational_metadata_score", "impact_readiness_tier", "recommended_next_action"], 12)}
        <h3>T1 impact upgrade plan</h3>
        <p class="muted">Execution layer: metadata intake, assay-design steps, and go/no-go criteria for moving from benchmark-clean pHLA candidates to stronger research evidence.</p>
        {table_html(t1_impact_assay, ["candidate_id", "assay_priority", "assay_id", "assay_or_review", "required_input", "current_status"], 24)}
        <h3>Nature-grade claim ladder</h3>
        <p class="muted">The upgrade path is explicit: supported framework claims now, blocked translational/clinical claims until metadata, public overlap audit, lockbox validation, and assay evidence are added.</p>
        {table_html(nature_claims, ["claim_level", "claim", "current_status", "evidence_now", "required_to_upgrade", "forbidden_overclaim"], 12)}
        <h3>Nature-grade validation gaps</h3>
        {table_html(nature_gaps, ["gap_id", "blocks_claim", "current_status", "minimum_fix", "impact_if_fixed"], 12)}
        <h3>Top BMA-ranked rows</h3>
        <p class="muted">Most current top rows still abstain because leakage/source-shift flags make clean benchmark claims invalid.</p>
        {table_html(top_bma, ["candidate_id", "barneo_bma_score", "patient_gated_bma_score", "bma_rank_global", "bma_confidence_score", "bma_confidence_bin", "bma_abstain", "bma_abstention_reason_primary", "selected_method_count", "selected_methods"], 25)}
        <h3>Non-abstain review prompts</h3>
        {table_html(clean_review, ["candidate_id", "barneo_bma_score", "patient_gated_bma_score", "bma_rank_global", "bma_confidence_score", "bma_confidence_bin", "selected_method_count", "selected_methods"], 25)}
        <h3>Abstention reason counts</h3>
        {table_html(abstention_summary, ["reason", "n"], 20)}
      </section>

      <section id="splits">
        <h2><span class="num">07</span>Split Robustness</h2>
        <p>All split sections stay visible even when metadata is weak. Source-heldout and HLA-heldout behavior are first-class outputs, not hidden failure modes.</p>
        {table_html(split_summary, ["split_contract", "n_metric_rows", "n_methods", "median_AUPRC", "median_top10_precision"], 20)}
        <h3>Source/HLA stress method positioning</h3>
        {table_html(stress_methods, ["method_name", "method_role", "mean_AUPRC", "source_heldout_median_delta_AUPRC", "hla_stress_median_delta_AUPRC", "korean_hla_median_delta_AUPRC", "stress_flags", "recommended_use"], 25)}
        <h3>Stress slice winners</h3>
        {table_html(stress_slices, ["split_contract", "split_group", "n_total_max", "n_pos_max", "overall_winner", "overall_winner_role", "overall_winner_AUPRC", "best_clean_method", "best_clean_AUPRC", "structure_lr_AUPRC", "clean_minus_anchor_AUPRC"], 45)}
        <h3>Korean-HLA top board</h3>
        {table_html(korean_hla_board[korean_hla_board["allele_rank"].le(3)] if not korean_hla_board.empty and "allele_rank" in korean_hla_board.columns else korean_hla_board, ["hla_allele", "method_name", "method_role", "best_contract", "n_total", "n_pos", "AUPRC", "top10_precision", "clean_claim_allowed"], 45)}
      </section>

      <section id="failure-aware">
        <h2><span class="num">08</span>Failure-Aware Reliability</h2>
        <p>The failure-aware layer applies observed source/HLA/leakage error modes to BAR-Neo-BMA. It separates clean non-abstain claims from practical manual review tiers.</p>
        {table_html(failure_aware.sort_values("failure_aware_score", ascending=False) if not failure_aware.empty and "failure_aware_score" in failure_aware.columns else failure_aware, ["candidate_id", "label", "source_name", "hla_allele_4digit", "base_patient_gated_bma_score", "failure_aware_score", "failure_aware_confidence", "failure_aware_review_tier", "failure_aware_abstain", "failure_aware_reason_primary"], 35)}
        {table_html(review_tier_summary, ["review_tier", "n"], 20)}
      </section>

      <section id="distribution-error">
        <h2><span class="num">09</span>Distribution Error Audit</h2>
        <p>This section links wrong predictions to training/external distribution, source prevalence, rare HLA support, low-prevalence stress groups, and leakage-risk regions.</p>
        <h3>Method vulnerability</h3>
        {table_html(distribution_vulnerability, ["method_name", "method_role", "n_scored", "source_score_prevalence_corr", "low_prevalence_high_ranked_negative_rate", "rare_hla_high_ranked_negative_rate", "rare_hla_missed_positive_rate", "external_missed_positive_rate", "distribution_vulnerability_score", "primary_distribution_issues"], 30)}
        <h3>Shift-linked error slices</h3>
        {table_html(shift_error_link, ["feature", "value", "train_n", "external_n", "train_positive_prevalence", "external_positive_prevalence", "distribution_shift_score", "max_error_pressure_score", "error_shift_link_score", "top_fp_methods", "top_fn_methods", "diagnosis"], 45)}
        <h3>Remediation plan</h3>
        {table_html(remediation_plan, ["priority", "failure_axis", "evidence", "recommended_action", "barneo_policy"], 12)}
      </section>

      <section id="contextual-bma">
        <h2><span class="num">10</span>Contextual BAR-Neo-BMA</h2>
        <p>Contextual BMA changes posterior expert weights by source, HLA support, low-prevalence stress, Korean-HLA focus, external/holdout status, and leakage risk. Public pretrained tools remain caveated support and are excluded from the clean contextual score view.</p>
        <h3>Top contextual candidates</h3>
        {table_html(contextual_scores, ["candidate_id", "source_name", "hla_allele_4digit", "contextual_bma_score", "clean_contextual_bma_score", "contextual_confidence_score", "contextual_abstain", "contextual_abstention_reason_primary", "selected_contextual_methods"], 35)}
        <h3>Context summary</h3>
        {table_html(contextual_summary, ["context_label", "n_candidates", "n_pos", "positive_prevalence", "mean_contextual_bma_score", "mean_clean_contextual_bma_score", "mean_confidence", "abstention_rate", "clean_claim_allowed", "caveated_public_used_rate"], 45)}
        <h3>Contextual method weights</h3>
        {table_html(contextual_weights, ["context_label", "method_name", "method_role", "method_family", "contextual_weight", "context_multiplier", "context_adjustment_reasons"], 55)}
      </section>

      <section id="challenge-pack">
        <h2><span class="num">11</span>Challenge Pack</h2>
        <p>The challenge pack converts current failures into next experiments: low-prevalence false-positive stress, rare-HLA rescue, high-score claim-blocked rows, external fragility, public/internal disagreement, Korean-HLA focus, and patient-gate metadata blockers.</p>
        <h3>Axis summary</h3>
        {table_html(challenge_summary, ["challenge_axis", "n_unique_candidates", "positive_prevalence", "mean_contextual_bma_score", "mean_clean_contextual_bma_score", "mean_confidence", "recommended_split_contract", "success_metric"], 30)}
        <h3>Top challenge rows</h3>
        {table_html(challenge_pack, ["global_challenge_rank", "challenge_axis", "candidate_id", "label", "source_name", "hla_allele_4digit", "contextual_bma_score", "clean_contextual_bma_score", "contextual_confidence_score", "challenge_reason", "recommended_next_action"], 55)}
        <h3>Next experiment plan</h3>
        {table_html(experiment_plan, ["priority", "experiment", "why", "success_metric", "expected_output"], 20)}
      </section>

      <section id="manual-queue">
        <h2><span class="num">12</span>Manual Review Queue</h2>
        <p>The manual queue is not a clean claim list. It separates high-score claim-blocked rows from rescue/watchlist rows and distribution challenge rows.</p>
        <h3>Priority Review Queue</h3>
        {table_html(manual_queue, ["candidate_id", "manual_review_priority_bin", "manual_review_priority_score", "failure_aware_review_tier", "label", "source_name", "hla_allele_4digit", "peptide", "base_patient_gated_bma_score", "failure_aware_score", "best_clean_internal_support", "best_caveated_public_support", "review_action_required"], 45)}
        <h3>Distribution Challenge Queue</h3>
        {table_html(challenge_queue, ["candidate_id", "challenge_axis", "manual_review_priority_bin", "manual_review_priority_score", "failure_aware_review_tier", "label", "source_name", "source_positive_prevalence", "hla_allele_4digit", "hla_allele_support_count", "peptide", "base_patient_gated_bma_score", "failure_aware_score", "review_action_required"], 60)}
      </section>

      <section id="public-audit">
        <h2><span class="num">13</span>Public Tool Overlap Audit</h2>
        <p><span class="warn">Public pretrained tools remain caveated.</span> Documentation-level provenance is not enough to call a public method a clean external baseline. Row-level candidate/peptide-HLA training-corpus overlap audit is required.</p>
        {table_html(public_audit[public_audit["uses_public_pretraining"].astype(bool)] if not public_audit.empty and "uses_public_pretraining" in public_audit.columns else public_audit, ["method_name", "method_role", "training_overlap_audit_status", "clean_comparator_allowed_after_audit", "reviewer_disposition", "caveat"], 30)}
        <h3>Overlap audit intake pack</h3>
        {table_html(public_intake, ["public_tool", "expected_drop_file", "template_file", "required_key", "clean_pass_rule", "current_clean_allowed", "current_disposition"], 30)}
        <h3>Row-Level Training Corpus Audit</h3>
        {table_html(public_training_row_summary, ["public_tool", "public_corpus_file", "n_candidate_overlaps", "candidate_overlap_fraction", "n_exact_peptide_hla", "n_exact_peptide", "n_near_peptide", "audit_status", "clean_comparator_allowed_after_row_audit"], 30)}
        <h3>Required Public Corpus Inputs</h3>
        {table_html(public_training_requirements, ["public_tool", "filename_tokens", "required_key", "minimum_columns", "clean_pass_rule", "drop_location"], 30)}
      </section>

      <section id="patient-gated">
        <h2><span class="num">14</span>PAAD/THCA Patient Gates</h2>
        <p><span class="warn">Current status: demo only.</span> Candidate rows lack enough disease timing, presentation, antigen expression, immune-context, and safety metadata for real patient-level confidence. All current patient-gated rows are research triage only and clinical_use=false.</p>
        <h3>Scenario gate matrix</h3>
        {table_html(patient_scenarios, ["scenario_id", "disease", "research_priority", "disease_context_gate", "presentation_gate", "antigen_gate", "immune_context_gate", "safety_gate", "scenario_gate_multiplier", "combination_strategy", "main_caveat"], 20)}
        <h3>Blocking metadata gaps</h3>
        {table_html(patient_blocked, ["field", "gate_group", "required_for", "metadata_status", "reviewer_safe_default_if_missing"], 60)}
        <h3>Demo gated queue</h3>
        {table_html(patient_queue, ["scenario_id", "candidate_id", "peptide", "hla_allele_4digit", "patient_gate_base_score", "model_confidence_gate", "scenario_gate_multiplier", "metadata_completion_cap", "demo_patient_gated_score", "patient_gated_confidence_bin", "abstention_reason_primary"], 45)}
      </section>

      <section id="caveats">
        <h2><span class="num">15</span>Caveats</h2>
        <div class="cards">
          <div class="card"><strong class="warn">Public tools:</strong><br>Public pretrained tools are caveated comparators until row-level training-corpus overlap audit is complete.</div>
          <div class="card"><strong class="warn">MHC class:</strong><br>Class I and Class II must not be pooled as a single predictor claim.</div>
          <div class="card"><strong class="warn">Clinical boundary:</strong><br>This is research triage, not clinical vaccine selection.</div>
          <div class="card"><strong class="warn">QK boundary:</strong><br>Allowed wording is bounded fallback/fusion component only.</div>
          <div class="card"><strong>Missing metadata:</strong><br><span class="muted">{html.escape(missing_text[:900])}</span></div>
          <div class="card"><strong>Manual review:</strong><br>High score and low confidence means investigate the row; it does not create a benchmark claim.</div>
        </div>
      </section>

      <section id="visuals">
        <h2><span class="num">16</span>Visual Dashboard</h2>
        <p>The separate visual dashboard contains 15 PNG panels for leaderboard ranking, source prevalence shift, abstention funnel, method vulnerability, challenge axes, contextual score-confidence, patient-gate metadata blockers, public-tool caveats, split heatmaps, anchor deltas, contextual weights, claim-safe top candidates, and stress-guarded ranking evidence.</p>
        <p><a href="clean_neobench_visual_dashboard_2026_05_10.html">Open CLEAN-NeoBench Visual Dashboard</a></p>
      </section>

      <section id="paths">
        <h2><span class="num">17</span>Sources + Paths</h2>
        <p class="path">Output root: {html.escape(str(output_root))}</p>
        <p class="path">Leaderboard: {html.escape(str(output_root / "clean_neobench_leaderboard.tsv"))}</p>
        <p class="path">Win/loss summary: {html.escape(str(output_root / "clean_neobench_winloss_method_summary.tsv"))}</p>
        <p class="path">Split winner board: {html.escape(str(output_root / "clean_neobench_split_winner_board.tsv"))}</p>
        <p class="path">Source/HLA stress audit: {html.escape(str(output_root / "clean_neobench_source_hla_stress_method_summary.tsv"))}</p>
        <p class="path">Korean-HLA method board: {html.escape(str(output_root / "clean_neobench_korean_hla_method_board.tsv"))}</p>
        <p class="path">Split metrics: {html.escape(str(output_root / "clean_neobench_split_metrics.tsv"))}</p>
        <p class="path">BMA candidate scores: {html.escape(str(output_root / "barneo_bma_candidate_scores.tsv"))}</p>
        <p class="path">Stress-guarded scores: {html.escape(str(output_root / "barneo_stress_guarded_candidate_scores.tsv"))}</p>
        <p class="path">High-impact leads: {html.escape(str(output_root / "barneo_high_impact_lead_candidates.tsv"))}</p>
        <p class="path">Reviewer kill audit: {html.escape(str(output_root / "barneo_high_impact_reviewer_kill_audit.tsv"))}</p>
        <p class="path">Top-pass reviewer evidence: {html.escape(str(output_root / "barneo_top_pass_reviewer_evidence.tsv"))}</p>
        <p class="path">High-impact lab handoff queue: {html.escape(str(output_root / "barneo_high_impact_lab_handoff_queue.tsv"))}</p>
        <p class="path">T1 translational readiness: {html.escape(str(output_root / "barneo_t1_translational_readiness.tsv"))}</p>
        <p class="path">T1 impact upgrade plan: {html.escape(str(output_root / "BAR_NEO_T1_IMPACT_UPGRADE_PLAN_KR.md"))}</p>
        <p class="path">T1 assay design matrix: {html.escape(str(output_root / "barneo_t1_assay_design_matrix.tsv"))}</p>
        <p class="path">BAR-Neo-NG scores: {html.escape(str(output_root / "barneo_ng_algorithm_scores.tsv"))}</p>
        <p class="path">BAR-Neo-NG T1 candidates: {html.escape(str(output_root / "barneo_ng_t1_candidates.tsv"))}</p>
        <p class="path">BAR-Neo-NG HTML: {html.escape(str(hub_root / "barneo_ng_algorithm_results_2026_05_10.html"))}</p>
        <p class="path">GA/RL bridge scores: {html.escape(str(output_root / "ga_rl_barneo_candidate_scores.tsv"))}</p>
        <p class="path">GA/RL bridge T1: {html.escape(str(output_root / "ga_rl_barneo_t1_confirmed.tsv"))}</p>
        <p class="path">GA/RL bridge unique T1 handoff: {html.escape(str(output_root / "ga_rl_barneo_t1_unique_candidates.tsv"))}</p>
        <p class="path">GA/RL high-GA blocked audit: {html.escape(str(output_root / "ga_rl_barneo_high_ga_blocked_audit.tsv"))}</p>
        <p class="path">GA/RL failure analysis: {html.escape(str(output_root / "ga_rl_barneo_failure_analysis.tsv"))}</p>
        <p class="path">GA/RL failure analysis KR: {html.escape(str(output_root / "GA_RL_BARNEO_FAILURE_ANALYSIS_KR.md"))}</p>
        <p class="path">GA/RL patient gate: {html.escape(str(output_root / "ga_rl_patient_gate_priority.tsv"))}</p>
        <p class="path">GA/RL patient gate summary: {html.escape(str(output_root / "ga_rl_patient_gate_priority_summary.tsv"))}</p>
        <p class="path">GA/RL patient gate KR: {html.escape(str(output_root / "GA_RL_BARNEO_PATIENT_GATE_KR.md"))}</p>
        <p class="path">GA/RL bridge HTML: {html.escape(str(hub_root / "ga_rl_barneo_bridge_results_2026_05_11.html"))}</p>
        <p class="path">Nature-grade evolution plan: {html.escape(str(output_root / "BAR_NEO_NATURE_GRADE_EVOLUTION_PLAN_KR.md"))}</p>
        <p class="path">Nature-grade claim ladder: {html.escape(str(output_root / "clean_neobench_nature_claim_ladder.tsv"))}</p>
        <p class="path">High-impact HTML: {html.escape(str(hub_root / "barneo_high_impact_candidate_pack_2026_05_10.html"))}</p>
        <p class="path">Top-pass reviewer evidence HTML: {html.escape(str(hub_root / "barneo_top_pass_reviewer_evidence_2026_05_10.html"))}</p>
        <p class="path">T1 translational readiness HTML: {html.escape(str(hub_root / "barneo_t1_translational_readiness_2026_05_10.html"))}</p>
        <p class="path">T1 impact upgrade HTML: {html.escape(str(hub_root / "barneo_t1_impact_upgrade_plan_2026_05_10.html"))}</p>
        <p class="path">Nature-grade evolution HTML: {html.escape(str(hub_root / "barneo_nature_grade_evolution_2026_05_10.html"))}</p>
        <p class="path">BMA method weights: {html.escape(str(output_root / "barneo_bma_method_weights.tsv"))}</p>
        <p class="path">BMA selector audit: {html.escape(str(output_root / "barneo_bma_selector_audit.tsv"))}</p>
        <p class="path">Distribution error audit: {html.escape(str(output_root / "clean_neobench_method_distribution_vulnerability.tsv"))}</p>
        <p class="path">Contextual BMA scores: {html.escape(str(output_root / "barneo_contextual_bma_candidate_scores.tsv"))}</p>
        <p class="path">Challenge pack: {html.escape(str(output_root / "clean_neobench_challenge_pack.tsv"))}</p>
        <p class="path">Public overlap intake: {html.escape(str(output_root / "clean_neobench_public_overlap_audit_intake.tsv"))}</p>
        <p class="path">Patient-gated demo: {html.escape(str(output_root / "patient_gated_clean_neo_candidate_queue.tsv"))}</p>
      </section>
    </main>
  </div>
</body>
</html>
"""
    page_path = hub_root / args.page_name
    page_path.write_text(html_text)
    update_manifest(
        output_root,
        "barneo_html_dossier",
        {
            "outputs": [str(page_path)],
            "warnings": ["HTML dossier is reviewer-safe research triage, not clinical vaccine selection or SOTA validation."],
        },
    )
    print(f"[barneo-dossier] wrote {page_path}")


if __name__ == "__main__":
    main()
