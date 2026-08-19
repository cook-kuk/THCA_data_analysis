#!/usr/bin/env python3
"""Build a physics-gate planning package for CROSS-Neo.

This layer sits after short MD and before expensive wetlab escalation.
It does not run new physics simulations; it ranks which candidates should
receive endpoint-energy, PMF, or FEP work and what each result would mean.
"""

from __future__ import annotations

import html
import json
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PKG = MD_OUT / "high_impact_decision_package"
LEAD = MD_OUT / "two_lead_impact_package"
PRECLIN = MD_OUT / "preclinical_validation_protocol"
FEEDBACK = MD_OUT / "assay_feedback_learner"
SCENARIO = MD_OUT / "assay_feedback_scenarios"
SIM = MD_OUT / "immunogenicity_simulation_escalation"
OUT = MD_OUT / "physics_gate_package"
FIG = MD_OUT / "figures"
HUB = REPO / "project/papers_hub_2026_05_04"
ASSET = HUB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_physics_gate.html"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
WEB_PAGE = WEB / "cross_neo_physics_gate.html"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def nfloat(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def esc(value: object) -> str:
    return html.escape("" if pd.isna(value) else str(value))


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return esc(value)


def table_html(df: pd.DataFrame, cols: list[str], n: int = 20) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    keep = [c for c in cols if c in df.columns]
    out = ["<table><thead><tr>"]
    out.extend(f"<th>{esc(c.replace('_', ' '))}</th>" for c in keep)
    out.append("</tr></thead><tbody>")
    for _, row in df.head(n).iterrows():
        out.append("<tr>")
        for c in keep:
            v = row.get(c, "")
            out.append(f"<td>{fmt(v) if isinstance(v, (int, float)) else esc(v)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


def md_table(df: pd.DataFrame, cols: list[str], n: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    keep = [c for c in cols if c in df.columns]
    return df[keep].head(n).to_markdown(index=False) if keep else "_No requested columns._"


def safe_copy(src: Path, dst: Path) -> None:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[physics-gate] skip permission-denied copy {dst}")


def load_inputs() -> dict[str, pd.DataFrame]:
    data = {
        "top13": read_tsv(PKG / "no_false_positive_top13_candidates.tsv"),
        "md_scores": read_tsv(MD_OUT / "md_evidence_scores.tsv"),
        "preclin": read_tsv(PRECLIN / "preclinical_readiness_scorecard.tsv"),
        "lead": read_tsv(LEAD / "two_lead_external_evidence.tsv"),
        "learner": read_tsv(FEEDBACK / "candidate_posterior_updates.tsv"),
        "scenarios": read_tsv(SCENARIO / "assay_feedback_scenario_outcomes.tsv"),
        "sim": read_tsv(SIM / "immunogenicity_simulation_matrix.tsv"),
        "batch": read_tsv(MD_OUT / "next_simulation_batch_recommendation.tsv"),
    }
    if data["top13"].empty or data["md_scores"].empty:
        raise FileNotFoundError("Required top13 or md_scores inputs are missing")
    return data


def build_physics_modules() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "module": "endpoint_mmgbsa_or_openmm_interaction_energy",
                "physics_question": "Does the mutant preserve a more favorable interface than WT/decoy over an ensemble?",
                "when_to_use": "After short MD exists; cheap ranking before PMF/FEP.",
                "input_requirement": "completed mutant/control trajectories or enough saved frames",
                "output": "trajectory endpoint interaction energies with bootstrap CI",
                "claim_boundary": "ranking diagnostic only; no activation or cytokine claim",
                "cost_class": "low",
            },
            {
                "module": "umbrella_or_steered_md_tcr_unbinding_pmf",
                "physics_question": "Does the mutant/TCR interface remain harder to unbind than WT/decoy?",
                "when_to_use": "Only if TCR-pMHC interface is already stable and multimer/TCR-binding decision depends on it.",
                "input_requirement": "stable TCR-pMHC model plus a defensible reaction coordinate",
                "output": "PMF / off-rate proxy / rupture path risk",
                "claim_boundary": "mechanistic proxy only; not cytokine or killing proof",
                "cost_class": "high",
            },
            {
                "module": "alchemical_fep_mutant_to_wt_delta_delta_g",
                "physics_question": "Is mutant more favorable than WT in a formally estimated free-energy sense?",
                "when_to_use": "Only for top flagship cases after WT mapping and stable 10ns controls.",
                "input_requirement": "clean WT mapping and stable starting structures",
                "output": "DeltaDeltaG with uncertainty",
                "claim_boundary": "does not prove immunogenicity; only supports specificity ranking",
                "cost_class": "very_high",
            },
            {
                "module": "structure_confidence_and_clash_crosscheck",
                "physics_question": "Is the starting pHLA/TCR-pMHC geometry physically plausible enough to justify expensive simulation?",
                "when_to_use": "Before any expensive physics module is launched.",
                "input_requirement": "modeled complex or template structure",
                "output": "clash / anchor / geometry sanity score",
                "claim_boundary": "presentation plausibility only",
                "cost_class": "very_low",
            },
        ]
    )


def build_candidate_physics_plan(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    top = data["top13"].copy()
    md = data["md_scores"].copy()
    pre = data["preclin"].copy()
    lead = data["lead"].copy()
    learner = data["learner"].copy()
    scenarios = data["scenarios"].copy()
    batch = data["batch"].copy()

    scenario_pressure = {}
    if not scenarios.empty:
        for peptide, sub in scenarios.groupby("peptide"):
            delta = pd.to_numeric(sub.get("posterior_delta", pd.Series(dtype=float)), errors="coerce").fillna(0)
            scenario_pressure[str(peptide)] = float(delta.abs().mean())

    rows = []
    for _, row in top.sort_values("priority_order").head(13).iterrows():
        peptide = str(row.get("peptide", ""))
        md_row = md[md["peptide"].astype(str).eq(peptide)] if not md.empty and "peptide" in md.columns else pd.DataFrame()
        pre_row = pre[pre["mutant_peptide"].astype(str).eq(peptide)] if not pre.empty and "mutant_peptide" in pre.columns else pd.DataFrame()
        lead_row = lead[lead["mutant_peptide"].astype(str).eq(peptide)] if not lead.empty and "mutant_peptide" in lead.columns else pd.DataFrame()
        learner_row = learner[learner["peptide"].astype(str).eq(peptide)] if not learner.empty and "peptide" in learner.columns else pd.DataFrame()
        batch_row = batch[batch["peptide"].astype(str).eq(peptide)] if not batch.empty and "peptide" in batch.columns else pd.DataFrame()
        paired = int(nfloat(row.get("paired_tcr_evidence_count"), 0))
        md_label = str(row.get("md_label", ""))
        md_score = nfloat(row.get("md_score"), nfloat(md_row.iloc[0].get("MD_evidence_score"), 0.0) if not md_row.empty else 0.0)
        priority = nfloat(row.get("priority_order"), 999)
        wetlab_priority = nfloat(row.get("wetlab_priority_score"), 0.0)
        preclin_ready = nfloat(pre_row.iloc[0].get("overall_readiness_score"), 0.0) if not pre_row.empty else 0.0
        learner_posterior = nfloat(learner_row.iloc[0].get("posterior_immunogenicity_readiness"), 0.0) if not learner_row.empty else preclin_ready
        scenario_delta = scenario_pressure.get(peptide, 0.0)
        wt_ready = bool(str(row.get("wt_or_decoy_ready", "")).strip()) or bool(not pre_row.empty)
        physics_priority_score = (
            0.28 * md_score
            + 0.18 * learner_posterior
            + 0.14 * wetlab_priority
            + 0.14 * min(paired / 20.0, 1.0)
            + 0.10 * min(scenario_delta / 0.25, 1.0)
            + 0.08 * (1.0 if wt_ready else 0.0)
            + 0.08 * (1.0 if "MD_VERY_STRONG" in md_label else 0.5 if "MD_STRONG" in md_label else 0.2)
        )
        if peptide in {"HMTEVVRHC", "GADGVGKSAL"}:
            physics_tier = "P0_FLAGSHIP_PHYSICS"
        elif physics_priority_score >= 0.60:
            physics_tier = "P1_PHYSICS_PRIORITY"
        elif physics_priority_score >= 0.40:
            physics_tier = "P2_PHYSICS_DIAGNOSTIC"
        else:
            physics_tier = "P3_NO_PHYSICS_NOW"

        if peptide == "HMTEVVRHC":
            module = "endpoint_mmgbsa_or_openmm_interaction_energy + umbrella_or_steered_md_tcr_unbinding_pmf + alchemical_fep_mutant_to_wt_delta_delta_g"
            next_step = "run endpoint energy now; hold PMF/FEP until WT/decoy control geometry is confirmed"
            cost_class = "high_to_very_high"
        elif peptide == "GADGVGKSAL":
            module = "endpoint_mmgbsa_or_openmm_interaction_energy + umbrella_or_steered_md_tcr_unbinding_pmf"
            next_step = "run endpoint energy now; PMF only if paired TCR model remains stable"
            cost_class = "high"
        elif physics_priority_score >= 0.60:
            module = "endpoint_mmgbsa_or_openmm_interaction_energy"
            next_step = "run cheap endpoint physics only; no expensive free-energy work yet"
            cost_class = "low"
        else:
            module = "structure_confidence_and_clash_crosscheck"
            next_step = "hold until assay or MD uncertainty justifies more physics budget"
            cost_class = "very_low"

        rows.append(
            {
                "row_id": row.get("row_id", ""),
                "peptide": peptide,
                "hla_4digit": row.get("hla_4digit", ""),
                "source_dataset": row.get("source_dataset", ""),
                "priority_order": priority,
                "paired_tcr_evidence_count": paired,
                "md_label": md_label,
                "md_score": md_score,
                "wt_or_decoy_ready": bool(wt_ready),
                "preclinical_readiness_score": preclin_ready,
                "learner_posterior_readiness": learner_posterior,
                "scenario_pressure": scenario_delta,
                "physics_priority_score": physics_priority_score,
                "physics_tier": physics_tier,
                "recommended_physics_module": module,
                "recommended_next_step": next_step,
                "cost_class": cost_class,
                "expected_runtime_class": "endpoint/low" if cost_class == "low" else "pmf/high" if cost_class == "high" else "fep/very_high" if cost_class == "very_high" else "diagnostic/very_low",
                "claim_boundary": "physics helps reject or rank; it does not prove immunogenicity",
            }
        )
    return pd.DataFrame(rows).sort_values(["physics_priority_score", "priority_order"], ascending=[False, True])


def build_budget(plan: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in plan.iterrows():
        if row["physics_tier"] == "P0_FLAGSHIP_PHYSICS":
            rows.append(
                {
                    "peptide": row["peptide"],
                    "hla_4digit": row["hla_4digit"],
                    "module": "endpoint_mmgbsa_or_openmm_interaction_energy",
                    "estimated_gpu_hours_low": 6.0,
                    "estimated_gpu_hours_high": 20.0,
                    "decision_use": "cheap ranking and consistency check",
                }
            )
            rows.append(
                {
                    "peptide": row["peptide"],
                    "hla_4digit": row["hla_4digit"],
                    "module": "umbrella_or_steered_md_tcr_unbinding_pmf",
                    "estimated_gpu_hours_low": 180.0,
                    "estimated_gpu_hours_high": 650.0,
                    "decision_use": "only if multimer/TCR-binding decision depends on it",
                }
            )
            rows.append(
                {
                    "peptide": row["peptide"],
                    "hla_4digit": row["hla_4digit"],
                    "module": "alchemical_fep_mutant_to_wt_delta_delta_g",
                    "estimated_gpu_hours_low": 250.0,
                    "estimated_gpu_hours_high": 900.0,
                    "decision_use": "WT specificity only, after the cheap gates pass",
                }
            )
        elif row["physics_tier"] == "P1_PHYSICS_PRIORITY":
            rows.append(
                {
                    "peptide": row["peptide"],
                    "hla_4digit": row["hla_4digit"],
                    "module": "endpoint_mmgbsa_or_openmm_interaction_energy",
                    "estimated_gpu_hours_low": 6.0,
                    "estimated_gpu_hours_high": 20.0,
                    "decision_use": "cheap ranking and consistency check",
                }
            )
    return pd.DataFrame(rows)


def build_rules() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("endpoint energy supports mutant > WT/decoy", "promote to PMF or FEP shortlist"),
            ("PMF supports stronger mutant unbinding barrier", "keep candidate alive for wetlab specificity testing"),
            ("FEP supports mutant > WT specificity", "use as highest-cost physics evidence before assay escalation"),
            ("WT or decoy looks as good as mutant", "hold mutant-specific claim and recheck specificity"),
            ("physics disagrees with MD and assay readiness", "treat as suspicion flag, not automatic rejection"),
            ("no paired TCR", "do not spend PMF/FEP on recognition claims"),
        ],
        columns=["observed_pattern", "physics_gate_response"],
    )


def build_transition_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ("endpoint_energy_pass", "Ranking support", "low", "pMHC or TCR-pMHC interface remains plausible"),
            ("PMF_pass", "Recognition plausibility", "high", "mutant interface survives unbinding pressure"),
            ("FEP_pass", "Specificity ranking", "very_high", "mutant wins over WT with uncertainty below margin"),
            ("WT_equally_good", "No specificity upgrade", "hold", "do not claim mutant-specific physics support"),
        ],
        columns=["physics_result", "claim_level", "action", "interpretation"],
    )


def make_figures(plan: pd.DataFrame, budget: pd.DataFrame, rules: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    if not plan.empty:
        top = plan.head(10)
        fig, ax = plt.subplots(figsize=(11.5, 5.0))
        ax.barh(np.arange(len(top)), top["physics_priority_score"], color="#4e79a7")
        ax.set_yticks(np.arange(len(top)), top["peptide"] + "\n" + top["physics_tier"])
        ax.invert_yaxis()
        ax.set_xlabel("physics priority score")
        ax.set_title("Physics gate priority")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIG / "fig_md64_physics_priority_rank.png", dpi=220)
        fig.savefig(FIG / "fig_md64_physics_priority_rank.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(11.5, 4.8))
        ax.axis("off")
        xs = np.linspace(0.08, 0.92, 4)
        labels = [
            ("endpoint energy", "cheap ranking"),
            ("PMF / umbrella", "off-rate proxy"),
            ("FEP", "specificity ranking"),
            ("wetlab", "activation/killing"),
        ]
        for i, (title, body) in enumerate(labels):
            ax.add_patch(plt.Rectangle((xs[i] - 0.085, 0.48), 0.17, 0.20, fc="#10243a", ec="#edf5ff", lw=1.2))
            ax.text(xs[i], 0.62, title, ha="center", color="#f2c46d", weight="bold", fontsize=10)
            ax.text(xs[i], 0.53, body, ha="center", color="white", fontsize=9)
            if i < len(labels) - 1:
                ax.annotate("", xy=(xs[i + 1] - 0.09, 0.58), xytext=(xs[i] + 0.09, 0.58), arrowprops=dict(arrowstyle="->", lw=2))
        ax.text(0.5, 0.84, "Physics gate ladder", ha="center", fontsize=15, weight="bold")
        ax.text(0.5, 0.23, "Expensive physics is reserved for candidates that already survived MD and claim-boundary checks.", ha="center", color="#9b1c31", fontsize=10)
        fig.savefig(FIG / "fig_md65_physics_gate_ladder.png", dpi=220, bbox_inches="tight")
        fig.savefig(FIG / "fig_md65_physics_gate_ladder.pdf", bbox_inches="tight")
        plt.close(fig)

    if not budget.empty:
        fig, ax = plt.subplots(figsize=(11.2, 4.9))
        budget = budget.sort_values("estimated_gpu_hours_high", ascending=False).head(10)
        ax.barh(np.arange(len(budget)), budget["estimated_gpu_hours_high"], color="#f28e2b")
        ax.set_yticks(np.arange(len(budget)), budget["peptide"] + "\n" + budget["module"])
        ax.invert_yaxis()
        ax.set_xlabel("estimated GPU hours high")
        ax.set_title("Physics budget by candidate/module")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIG / "fig_md66_physics_budget.png", dpi=220)
        fig.savefig(FIG / "fig_md66_physics_budget.pdf")
        plt.close(fig)

    counts = rules["physics_gate_response"].value_counts()
    fig, ax = plt.subplots(figsize=(10.4, 4.8))
    ax.barh(np.arange(len(counts)), counts.values, color="#59a14f")
    ax.set_yticks(np.arange(len(counts)), counts.index)
    ax.invert_yaxis()
    ax.set_xlabel("rule count")
    ax.set_title("Physics gate decision rules")
    fig.tight_layout()
    fig.savefig(FIG / "fig_md67_physics_gate_rules.png", dpi=220)
    fig.savefig(FIG / "fig_md67_physics_gate_rules.pdf")
    plt.close(fig)


def write_report(modules: pd.DataFrame, plan: pd.DataFrame, budget: pd.DataFrame, rules: pd.DataFrame) -> None:
    top2 = plan[plan["peptide"].isin(["HMTEVVRHC", "GADGVGKSAL"])]
    lines = [
        "# CROSS-Neo Physics Gate Package",
        "",
        "## Executive verdict",
        "",
        "This is the last expensive computational gate before wetlab specificity testing. Use cheap interface energies first, PMF/umbrella only when TCR-plausibility matters, and FEP only for the two flagship cases after WT mapping is clean.",
        "",
        "## Physics modules",
        "",
        md_table(modules, ["module", "physics_question", "when_to_use", "output", "claim_boundary", "cost_class"], 20),
        "",
        "## Candidate plan",
        "",
        md_table(plan, ["peptide", "hla_4digit", "physics_priority_score", "physics_tier", "recommended_physics_module", "recommended_next_step", "cost_class"], 20),
        "",
        "## Flagship physics candidates",
        "",
        md_table(top2, ["peptide", "hla_4digit", "physics_priority_score", "physics_tier", "recommended_physics_module", "recommended_next_step"], 10),
        "",
        "## Budget",
        "",
        md_table(budget, ["peptide", "hla_4digit", "module", "estimated_gpu_hours_low", "estimated_gpu_hours_high", "decision_use"], 20),
        "",
        "## Rules",
        "",
        md_table(rules, ["observed_pattern", "physics_gate_response"], 20),
        "",
        "## Claim boundary",
        "",
        "Physics supports rejection, ranking, and specificity sanity checks. It does not prove activation, killing, clinical utility, or universal immunogenicity.",
    ]
    (OUT / "PHYSICS_GATE_REPORT.md").write_text("\n".join(lines) + "\n")

    kr = [
        "# physics gate 요약",
        "",
        "이 패키지는 MD 다음 단계에서 추가로 돈이 많이 드는 physics simulation을 어디에 쓸지 정하는 계획입니다.",
        "",
        "- endpoint energy: 싸고 빠른 ranking",
        "- umbrella / PMF: TCR unbinding 또는 off-rate proxy",
        "- FEP: WT 대비 mutant specificity의 가장 비싼 검증",
        "",
        "두 flagship 후보만 FEP 급으로 올립니다: HMTEVVRHC, GADGVGKSAL.",
        "나머지는 cheap physics diagnostic에만 머무르게 합니다.",
        "",
        "중요: 이것도 immunogenicity proof가 아니라 physics gate입니다.",
    ]
    (OUT / "PHYSICS_GATE_ONE_PAGE_KR.md").write_text("\n".join(kr) + "\n")


def write_html(modules: pd.DataFrame, plan: pd.DataFrame, budget: pd.DataFrame, rules: pd.DataFrame, summary: dict) -> None:
    fig_cards = "\n".join(
        f"""
        <article>
          <img src="assets/cross_neo_md_audit/{name}.png" alt="{title}">
          <h3>{title}</h3>
          <p>{desc}</p>
          <a href="assets/cross_neo_md_audit/{name}.png">PNG</a>
          <a href="assets/cross_neo_md_audit/{name}.pdf">PDF</a>
        </article>
        """
        for name, title, desc in [
            ("fig_md64_physics_priority_rank", "Physics Priority", "Which candidates should receive expensive physics work first."),
            ("fig_md65_physics_gate_ladder", "Physics Ladder", "Endpoint energy, PMF, and FEP as escalating physics skepticism."),
            ("fig_md66_physics_budget", "Physics Budget", "Estimated GPU-hours for the expensive modules."),
            ("fig_md67_physics_gate_rules", "Physics Rules", "Rules that map physics outcomes to claim states."),
        ]
    )
    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo Physics Gate</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5;text-decoration:none}} header{{padding:40px 30px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1360px;margin:0 auto;padding:24px 30px 70px}} h1{{font-family:Georgia,serif;font-size:44px;margin:0;color:#fff8e8}}
h2{{font-family:Georgia,serif;color:#fff2d0}} .lead,.muted{{color:#aebdd1}} section{{border-bottom:1px solid #26364d;padding:22px 0}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin-top:18px}} .stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:12px}}
.stat b{{display:block;color:#f2c46d;font-size:26px}} table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{background:#13243a;color:#f2c46d}}
.grid2{{display:grid;grid-template-columns:repeat(2,1fr);gap:16px}} article{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:13px}} img{{width:100%;background:white;border-radius:6px}}
.warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:12px 14px}} .links{{display:flex;flex-wrap:wrap;gap:10px}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
@media(max-width:900px){{.grid2{{grid-template-columns:1fr}} h1{{font-size:32px}}}}
</style></head><body>
<header><div class="wrap">
<h1>CROSS-Neo Physics Gate</h1>
<p class="lead">Cheap interface-energy ranking first, PMF/umbrella only when TCR-plausibility matters, and FEP only for the two flagship cases after WT mapping is clean.</p>
<div class="stats">
<div class="stat"><b>{summary['n_candidates']}</b><span>candidates</span></div>
<div class="stat"><b>{summary['n_flagships']}</b><span>flagship physics cases</span></div>
<div class="stat"><b>{summary['n_physics_modules']}</b><span>modules</span></div>
<div class="stat"><b>{summary['n_budget_rows']}</b><span>budget rows</span></div>
</div></div></header>
<main class="wrap">
<section><h2>Boundary</h2>
<div class="warn">Physics is a skepticism filter, not immunogenicity proof. It ranks or rejects candidates before wetlab spending, but cannot prove activation or killing.</div>
<div class="links">
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_GATE_REPORT.md">Full report</a>
<a class="pill" href="assets/cross_neo_md_audit/PHYSICS_GATE_ONE_PAGE_KR.md">Korean brief</a>
<a class="pill" href="cross_neo_high_impact_decision.html">high-impact page</a>
<a class="pill" href="cross_neo_assay_feedback_loop.html">assay feedback loop</a>
<a class="pill" href="cross_neo_assay_feedback_scenarios.html">assay scenario simulator</a>
<a class="pill" href="cross_neo_physics_launch_sheet.html">physics launch sheet</a>
<a class="pill" href="cross_neo_endpoint_energy_gate.html">endpoint energy gate</a>
<a class="pill" href="cross_neo_physics_case_study_board.html">physics case study board</a>
<a class="pill" href="cross_neo_integrated_decision_matrix.html">integrated decision matrix</a>
<a class="pill" href="cross_neo_decision_atlas.html">decision atlas</a>
<a class="pill" href="cross_neo_executive_impact_console.html">executive impact console</a>
<a class="pill" href="cross_neo_flagship_execution_packet.html">flagship execution packet</a>
<a class="pill" href="cross_neo_flagship_command_center.html">flagship command center</a>
<a class="pill" href="cross_neo_flagship_decision_tower.html">flagship decision tower</a>
<a class="pill" href="assets/cross_neo_md_audit/physics_gate_candidate_plan.tsv">candidate plan TSV</a>
<a class="pill" href="assets/cross_neo_md_audit/physics_gate_budget.tsv">budget TSV</a>
<a class="pill" href="assets/cross_neo_md_audit/physics_gate_rules.tsv">rules TSV</a>
</div></section>
<section><h2>Figures</h2><div class="grid2">{fig_cards}</div></section>
<section><h2>Physics Modules</h2>{table_html(modules, ['module','physics_question','when_to_use','output','claim_boundary','cost_class'], 20)}</section>
<section><h2>Candidate Plan</h2>{table_html(plan, ['peptide','hla_4digit','physics_priority_score','physics_tier','recommended_physics_module','recommended_next_step','cost_class'], 20)}</section>
<section><h2>Budget</h2>{table_html(budget, ['peptide','hla_4digit','module','estimated_gpu_hours_low','estimated_gpu_hours_high','decision_use'], 20)}</section>
<section><h2>Rules</h2>{table_html(rules, ['observed_pattern','physics_gate_response'], 20)}</section>
</main></body></html>"""
    PAGE.write_text(html_text)


def deploy(paths: list[Path]) -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for path in paths:
        if path.exists():
            safe_copy(path, ASSET / path.name)
            safe_copy(path, WEB_ASSET / path.name)
    for fig in [
        "fig_md64_physics_priority_rank",
        "fig_md65_physics_gate_ladder",
        "fig_md66_physics_budget",
        "fig_md67_physics_gate_rules",
    ]:
        for ext in [".png", ".pdf"]:
            path = FIG / f"{fig}{ext}"
            if path.exists():
                safe_copy(path, ASSET / path.name)
                safe_copy(path, WEB_ASSET / path.name)
    WEB.mkdir(parents=True, exist_ok=True)
    safe_copy(PAGE, WEB_PAGE)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    modules = build_physics_modules()
    plan = build_candidate_physics_plan(data)
    budget = build_budget(plan)
    rules = build_rules()
    transition = build_transition_matrix()
    modules.to_csv(OUT / "physics_modules.tsv", sep="\t", index=False)
    plan.to_csv(OUT / "physics_gate_candidate_plan.tsv", sep="\t", index=False)
    budget.to_csv(OUT / "physics_gate_budget.tsv", sep="\t", index=False)
    rules.to_csv(OUT / "physics_gate_rules.tsv", sep="\t", index=False)
    transition.to_csv(OUT / "physics_gate_transition_matrix.tsv", sep="\t", index=False)
    make_figures(plan, budget, rules)
    write_report(modules, plan, budget, rules)
    summary = {
        "n_candidates": int(len(plan)),
        "n_flagships": int((plan["physics_tier"] == "P0_FLAGSHIP_PHYSICS").sum()),
        "n_physics_modules": int(len(modules)),
        "n_budget_rows": int(len(budget)),
        "figures": [
            "fig_md64_physics_priority_rank",
            "fig_md65_physics_gate_ladder",
            "fig_md66_physics_budget",
            "fig_md67_physics_gate_rules",
        ],
        "boundary": "physics planning and skepticism gating, not immunogenicity proof",
    }
    (OUT / "physics_gate_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_html(modules, plan, budget, rules, summary)
    deploy(
        [
            OUT / "physics_modules.tsv",
            OUT / "physics_gate_candidate_plan.tsv",
            OUT / "physics_gate_budget.tsv",
            OUT / "physics_gate_rules.tsv",
            OUT / "physics_gate_transition_matrix.tsv",
            OUT / "PHYSICS_GATE_REPORT.md",
            OUT / "PHYSICS_GATE_ONE_PAGE_KR.md",
            OUT / "physics_gate_summary.json",
        ]
    )
    print(json.dumps(summary, indent=2))
    print(f"[physics-gate] wrote {OUT}")
    print(f"[physics-gate] page {PAGE}")


if __name__ == "__main__":
    main()
