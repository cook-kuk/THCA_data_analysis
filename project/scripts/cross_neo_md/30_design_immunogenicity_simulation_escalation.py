#!/usr/bin/env python3
"""Design an immunogenicity-oriented simulation escalation layer.

The purpose is to add simulation work that is useful before immunogenicity
wetlab testing without overclaiming what simulation can prove.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PKG = MD_OUT / "high_impact_decision_package"
DIV = MD_OUT / "diverse_simulation_plan"
OUT = MD_OUT / "immunogenicity_simulation_escalation"
FIG = MD_OUT / "figures"


GPU_HOURS_PER_NS_LOW = 0.6
GPU_HOURS_PER_NS_HIGH = 1.2


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def clean_key(peptide: str, hla: str) -> str:
    return f"{peptide}|{hla}".replace("*", "").replace(":", "").replace("/", "_")


def has_text(value: object) -> bool:
    s = str(value).strip().lower()
    return bool(s) and s not in {"nan", "none", "null", "na"}


def build_existing_job_index() -> dict[str, list[dict]]:
    jobs = read_tsv(DIV / "diverse_simulation_job_plan.tsv")
    if jobs.empty:
        return {}
    idx: dict[str, list[dict]] = {}
    for _, row in jobs.iterrows():
        key = clean_key(str(row.get("peptide", "")), str(row.get("hla_4digit", "")))
        idx.setdefault(key, []).append(row.to_dict())
    return idx


def existing_status(candidate: pd.Series, scope: str, control: str = "") -> tuple[str, str]:
    key = clean_key(str(candidate.get("peptide", "")), str(candidate.get("hla_4digit", "")))
    matches = []
    for job in EXISTING_JOBS.get(key, []):
        if scope and str(job.get("complex_scope", "")) != scope:
            continue
        if control and control not in str(job.get("control_type", "")):
            continue
        matches.append(job)
    if not matches:
        return "not_manifested", "needs_structure_or_job_generation"
    if any(str(j.get("launch_group", "")) == "tier1_gad_full_replicates_controls" for j in matches):
        return "queued_or_ready_in_diverse_pack", "run_immediate_ready_batch.sh or queued watcher"
    if any("after_hmtevv_sync" in str(j.get("launch_group", "")) for j in matches):
        return "hold_after_hmtevv_sync", "finish/sync current HMTEVVRHC run before launch"
    if any(str(j.get("launch_group", "")) == "tier0_short_screens" for j in matches):
        return "ready_short_screen", "run_tier0_short_screens.sh"
    return "manifested", "see diverse_simulation_job_plan.tsv"


def add_row(rows: list[dict], candidate: pd.Series, **kwargs) -> None:
    total_ns = float(kwargs.get("replicates", 0)) * float(kwargs.get("ns_per_rep", 0))
    est_low = total_ns * GPU_HOURS_PER_NS_LOW
    est_high = total_ns * GPU_HOURS_PER_NS_HIGH
    if kwargs.get("simulation_module", "").startswith("alchemical_fep"):
        est_low, est_high = 250.0, 900.0
    if kwargs.get("simulation_module", "").startswith("umbrella"):
        est_low, est_high = 180.0, 650.0
    rows.append(
        {
            "priority_order": candidate.get("priority_order", ""),
            "row_id": candidate.get("row_id", ""),
            "peptide": candidate.get("peptide", ""),
            "hla_4digit": candidate.get("hla_4digit", ""),
            "source_dataset": candidate.get("source_dataset", ""),
            "label_binary": candidate.get("label_binary", candidate.get("actual_label", "")),
            "main_dl_score": candidate.get("main_dl_score", ""),
            "tcr_augmented_score_mean": candidate.get("tcr_augmented_score_mean", ""),
            "paired_tcr_evidence_count": candidate.get("paired_tcr_evidence_count", 0),
            "baker_structural_score": candidate.get("baker_structural_score", 0),
            "md_label": candidate.get("md_label", ""),
            "assay_readout_goal": kwargs["assay_readout_goal"],
            "simulation_module": kwargs["simulation_module"],
            "escalation_tier": kwargs["escalation_tier"],
            "complex_scope": kwargs["complex_scope"],
            "control_type": kwargs["control_type"],
            "replicates": kwargs.get("replicates", 0),
            "ns_per_rep": kwargs.get("ns_per_rep", 0),
            "total_requested_ns": total_ns,
            "estimated_gpu_hours_low": est_low,
            "estimated_gpu_hours_high": est_high,
            "launch_status": kwargs.get("launch_status", "not_manifested"),
            "blocking_dependency": kwargs.get("blocking_dependency", ""),
            "why_needed": kwargs["why_needed"],
            "go_no_go_signal": kwargs["go_no_go_signal"],
            "expected_output": kwargs["expected_output"],
            "claim_boundary": kwargs["claim_boundary"],
        }
    )


def build_matrix(top13: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for _, cand in top13.sort_values("priority_order").iterrows():
        rank = int(float(cand.get("priority_order", 999)))
        paired = float(cand.get("paired_tcr_evidence_count", 0) or 0) > 0
        structure_score = float(cand.get("baker_structural_score", 0) or 0)
        md_label = str(cand.get("md_label", "")).lower()
        wt_status = str(cand.get("wt_status", "")).lower()
        decoy_status = str(cand.get("anchor_preserved_decoy", "")).lower()
        wt_ready = has_text(wt_status) and "missing" not in wt_status and "blocked" not in wt_status
        decoy_ready = has_text(decoy_status) and "missing" not in decoy_status
        top4 = rank <= 4
        top2 = rank <= 2

        add_row(
            rows,
            cand,
            assay_readout_goal="HLA binding/stability prefilter",
            simulation_module="pmhc_structure_crosscheck_boltz_chai_af_or_template",
            escalation_tier="T0_cheap_structure",
            complex_scope="pMHC",
            control_type="mutant",
            replicates=0,
            ns_per_rep=0,
            launch_status="already_supported" if structure_score > 0 else "structure_model_needed",
            blocking_dependency="install_or_remote_tool_if_no_template" if structure_score <= 0 else "",
            why_needed="Before T-cell assay, confirm the peptide can be modeled in the HLA groove and is not an obvious structural failure.",
            go_no_go_signal="Reject or hold if peptide cannot sit in the groove or structure confidence/contact geometry is implausible.",
            expected_output="structure confidence, peptide groove pose, anchor contact sanity, static clash/contact score",
            claim_boundary="Presentation plausibility only; not recognition or immunogenicity.",
        )

        pmhc_status, pmhc_block = existing_status(cand, "pMHC", "mutant")
        add_row(
            rows,
            cand,
            assay_readout_goal="HLA binding/stability prefilter",
            simulation_module="openmm_pmhc_mutant_3x10ns",
            escalation_tier="T1_short_md",
            complex_scope="pMHC",
            control_type="mutant",
            replicates=3 if top4 else 1,
            ns_per_rep=10 if top4 else 5,
            launch_status=pmhc_status if top4 else "stage_after_top4_or_if_wetlab_slot_opens",
            blocking_dependency=pmhc_block if top4 else "run only after Tier-1 candidates are resolved",
            why_needed="Tests whether the candidate peptide remains anchored in MHC before spending PBMC/TCR assay resources.",
            go_no_go_signal="Promote only if peptide RMSD/drift plateau and pMHC native contacts remain stable across replicates.",
            expected_output="peptide RMSD, groove drift, anchor contacts, pMHC contact occupancy, native Q",
            claim_boundary="Supports peptide-HLA stability only; cannot prove immunogenic T-cell response.",
        )

        add_row(
            rows,
            cand,
            assay_readout_goal="Mutant specificity and WT cross-reactivity risk",
            simulation_module="openmm_wt_and_anchor_decoy_pmhc_3x10ns",
            escalation_tier="T1_specificity_controls",
            complex_scope="pMHC",
            control_type="wildtype_and_anchor_decoy",
            replicates=3 if top4 else 1,
            ns_per_rep=10 if top4 else 5,
            launch_status="ready_or_tentative" if wt_ready or decoy_ready else "blocked_wt_or_decoy_curation",
            blocking_dependency="manual WT confirmation and anchor-preserved decoy structure prep" if not (wt_ready and decoy_ready) else "structure prep then OpenMM",
            why_needed="Immunogenicity assays must show mutant response over WT/decoy controls, not just a stable peptide.",
            go_no_go_signal="Promote if mutant pMHC stability and TCR-facing exposure exceed WT/decoy; hold if WT is equally stable and TCR-facing.",
            expected_output="mutant-WT/decoy delta RMSD, delta contacts, anchor retention, exposed mutation proxy",
            claim_boundary="Specificity triage only; WT/decoy wetlab remains required.",
        )

        if paired:
            tcr_status, tcr_block = existing_status(cand, "TCR-pMHC", "mutant")
            add_row(
                rows,
                cand,
                assay_readout_goal="TCR binding or multimer follow-up",
                simulation_module="openmm_tcr_pmhc_mutant_3x10ns",
                escalation_tier="T1_tcr_recognition_md",
                complex_scope="TCR-pMHC",
                control_type="mutant",
                replicates=3 if top4 else 1,
                ns_per_rep=10 if top4 else 5,
                launch_status=tcr_status if top4 else "stage_after_top4_or_if_paired_tcr_clone_available",
                blocking_dependency=tcr_block if top4 else "requires paired TCR structure/template and assay relevance check",
                why_needed="Directly checks whether TCR-peptide and CDR3 contacts persist, which is closer to recognition than pMHC alone.",
                go_no_go_signal="Promote if CDR3-peptide contacts and mutation-site contacts persist; deprioritize if TCR mostly contacts MHC or interface falls apart.",
                expected_output="TCR-peptide contacts, CDR3 contacts, mutation-site occupancy, TCR-pMHC native Q, crossing-angle proxy",
                claim_boundary="Recognition plausibility only; cytokine/killing must be measured experimentally.",
            )
            add_row(
                rows,
                cand,
                assay_readout_goal="Mutant-specific TCR recognition and WT cross-reactivity",
                simulation_module="openmm_tcr_pmhc_wt_decoy_3x10ns",
                escalation_tier="T2_tcr_specificity_controls",
                complex_scope="TCR-pMHC",
                control_type="wildtype_and_anchor_decoy",
                replicates=3 if top2 else 1,
                ns_per_rep=10,
                launch_status="blocked_wt_or_decoy_tcr_structure" if not wt_ready else "structure_prep_needed",
                blocking_dependency="WT/decoy TCR-pMHC model and manual TCR relevance check",
                why_needed="A positive TCR interface is not enough; the mutant must beat WT/decoy controls to support wetlab specificity.",
                go_no_go_signal="Promote if mutant has higher persistent CDR3/mutation-site contacts than WT/decoy.",
                expected_output="TCR-peptide contact delta, CDR3 delta, interface stability delta, WT cross-reactivity risk",
                claim_boundary="Does not prove lack of self-reactivity; only flags relative structural risk.",
            )

        add_row(
            rows,
            cand,
            assay_readout_goal="Endpoint affinity ranking for HLA and TCR interface",
            simulation_module="endpoint_mmgbsa_or_openmm_interaction_energy",
            escalation_tier="T2_endpoint_energy",
            complex_scope="pMHC_and_TCR-pMHC_if_available",
            control_type="mutant_vs_wt_decoy",
            replicates=0,
            ns_per_rep=0,
            launch_status="analysis_after_10ns_trajectories",
            blocking_dependency="requires completed mutant/control trajectories",
            why_needed="Cheaply summarizes trajectory ensembles into relative interface-energy features before deciding on FEP/umbrella sampling.",
            go_no_go_signal="Promote only if energy/contact ranking is consistent with mutant > WT/decoy and replicate variance is acceptable.",
            expected_output="trajectory endpoint interaction energies, MM/GBSA if AmberTools exists, bootstrap CI",
            claim_boundary="Endpoint energies are noisy; use as ranking diagnostics only.",
        )

        if top2:
            add_row(
                rows,
                cand,
                assay_readout_goal="High-confidence mutant-WT specificity",
                simulation_module="alchemical_fep_mutant_to_wt_delta_delta_g",
                escalation_tier="T3_expensive_free_energy",
                complex_scope="pMHC_and_TCR-pMHC_if_available",
                control_type="mutant_to_wildtype",
                replicates=1,
                ns_per_rep=0,
                launch_status="hold_until_10ns_md_and_wt_confirmed",
                blocking_dependency="requires clean WT mapping, stable starting structures, and FEP toolchain",
                why_needed="Only worth running for the top cases if WT/decoy controls are central to the wetlab decision.",
                go_no_go_signal="Promote if mutant improves pMHC/TCR-pMHC free-energy proxies versus WT with uncertainty below decision margin.",
                expected_output="DeltaDeltaG with uncertainty for mutant vs WT in pMHC/TCR-pMHC context",
                claim_boundary="Free energy still does not simulate T-cell activation.",
            )
            if paired:
                add_row(
                    rows,
                    cand,
                    assay_readout_goal="TCR off-rate / multimer plausibility",
                    simulation_module="umbrella_or_steered_md_tcr_unbinding_pmf",
                    escalation_tier="T3_expensive_pmf",
                    complex_scope="TCR-pMHC",
                    control_type="mutant_vs_wt_decoy",
                    replicates=1,
                    ns_per_rep=0,
                    launch_status="hold_until_stable_tcr_pmhc_10ns",
                    blocking_dependency="requires stable TCR-pMHC interface and carefully chosen reaction coordinate",
                    why_needed="If wetlab multimer/TCR binding is decisive, a PMF/off-rate proxy can rank mutant vs WT/decoy.",
                    go_no_go_signal="Promote if mutant PMF/contact persistence exceeds WT/decoy and is consistent with multimer assay.",
                    expected_output="TCR-pMHC unbinding PMF proxy, interface rupture path, off-rate risk flag",
                    claim_boundary="PMF is a mechanistic proxy, not cytokine response.",
                )

        add_row(
            rows,
            cand,
            assay_readout_goal="Cytokine/ICS/killing assay interpretation",
            simulation_module="kinetic_proofreading_activation_proxy",
            escalation_tier="T4_diagnostic_model",
            complex_scope="TCR-pMHC",
            control_type="mutant_vs_wt_decoy",
            replicates=0,
            ns_per_rep=0,
            launch_status="analysis_only_after_tcr_binding_proxy",
            blocking_dependency="requires experimental or simulated TCR binding/off-rate proxies",
            why_needed="Maps structural/TCR binding proxies to expected assay risk classes without pretending to simulate cell biology.",
            go_no_go_signal="Use as a post-hoc interpretation layer only; never as the sole wetlab trigger.",
            expected_output="activation-risk class, sensitivity to koff/contact-persistence assumptions",
            claim_boundary="T-cell activation, IFN-gamma, and killing cannot be proven by simulation.",
        )

    return pd.DataFrame(rows)


def build_bridge() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "wetlab_assay": "HLA binding / stability assay",
                "simulation_support": "pMHC structure crosscheck; mutant pMHC 3x10ns; WT/decoy pMHC controls",
                "decision_use": "Cull peptides that do not remain in the groove before T-cell assays.",
                "not_supported_claim": "Does not prove TCR recognition.",
            },
            {
                "wetlab_assay": "pMHC multimer / tetramer with known TCR or T cells",
                "simulation_support": "TCR-pMHC 3x10ns; CDR3-peptide contacts; TCR unbinding PMF only for top cases",
                "decision_use": "Prioritize candidates whose TCR interface is persistent and peptide-facing.",
                "not_supported_claim": "Does not prove cytokine secretion or killing.",
            },
            {
                "wetlab_assay": "IFN-gamma ELISpot / ICS",
                "simulation_support": "TCR contact persistence, mutation-site contacts, WT/decoy specificity, kinetic proofreading proxy",
                "decision_use": "Interpret why a candidate might pass/fail activation assays.",
                "not_supported_claim": "Cannot simulate complete immune-cell activation.",
            },
            {
                "wetlab_assay": "Killing assay",
                "simulation_support": "Only upstream presentation/recognition plausibility plus WT/decoy risk",
                "decision_use": "Run only after binding and activation screen are positive.",
                "not_supported_claim": "No simulation here proves tumor killing.",
            },
        ]
    )


def build_go_no_go_rules() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "stage": "T0 structure",
                "go_rule": "Peptide modeled in groove; anchor contacts plausible; no severe clashes; confidence acceptable.",
                "no_go_or_hold": "Unsupported HLA/length, failed groove placement, heavy clash, or missing WT for specificity-critical case.",
            },
            {
                "stage": "T1 10ns pMHC MD",
                "go_rule": "Peptide RMSD/drift plateaus; native pMHC contacts and anchors persist in most replicates.",
                "no_go_or_hold": "Peptide exits groove, anchor contacts collapse, or replicate behavior is inconsistent.",
            },
            {
                "stage": "T1/T2 TCR-pMHC MD",
                "go_rule": "Persistent TCR-peptide/CDR3 contacts; mutation-site contact if biologically relevant; TCR not only MHC-focused.",
                "no_go_or_hold": "Interface falls apart, contacts are MHC-only, or mutation residue is buried away from TCR in all models.",
            },
            {
                "stage": "WT/decoy specificity",
                "go_rule": "Mutant shows stronger pMHC/TCR interface evidence than WT and anchor-preserved decoy.",
                "no_go_or_hold": "WT/decoy looks equally good or better; self-reactivity risk requires wetlab caution.",
            },
            {
                "stage": "T3 expensive FEP/PMF",
                "go_rule": "Only launch if 10ns mutant/control results are stable and wetlab decision depends on specificity ranking.",
                "no_go_or_hold": "Do not spend FEP/PMF budget on unstable, uncurated, or non-TCR-available cases.",
            },
        ]
    )


def summarize_budget(matrix: pd.DataFrame) -> pd.DataFrame:
    active = matrix[matrix["total_requested_ns"] > 0].copy()
    if active.empty:
        return pd.DataFrame()
    group_cols = ["peptide", "hla_4digit"]
    return (
        active.groupby(group_cols, dropna=False)
        .agg(
            planned_jobs=("simulation_module", "count"),
            requested_ns=("total_requested_ns", "sum"),
            gpu_hours_low=("estimated_gpu_hours_low", "sum"),
            gpu_hours_high=("estimated_gpu_hours_high", "sum"),
            blocked_jobs=("launch_status", lambda s: int(s.astype(str).str.contains("blocked|hold", case=False, regex=True).sum())),
        )
        .reset_index()
        .sort_values(["requested_ns", "planned_jobs"], ascending=False)
    )


def make_figures(matrix: pd.DataFrame, bridge: pd.DataFrame, budget: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    tier_order = ["T0_cheap_structure", "T1_short_md", "T1_specificity_controls", "T1_tcr_recognition_md", "T2_tcr_specificity_controls", "T2_endpoint_energy", "T3_expensive_free_energy", "T3_expensive_pmf", "T4_diagnostic_model"]
    tier_counts = matrix["escalation_tier"].value_counts().reindex(tier_order).dropna()
    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.bar(tier_counts.index, tier_counts.values, color=["#4777b4", "#59a14f", "#edc948", "#af7aa1", "#ff9da7", "#9c755f", "#e15759", "#f28e2b", "#76b7b2"][: len(tier_counts)])
    ax.set_ylabel("planned rows")
    ax.set_title("Immunogenicity simulation escalation ladder")
    ax.tick_params(axis="x", rotation=35, labelsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md31_immunogenicity_simulation_ladder.png", dpi=220)
    fig.savefig(FIG / "fig_md31_immunogenicity_simulation_ladder.pdf")
    plt.close(fig)

    if not budget.empty:
        top = budget.head(10).copy()
        labels = top["peptide"] + "\n" + top["hla_4digit"]
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.bar(labels, top["requested_ns"], color="#4e79a7")
        ax.set_ylabel("requested ns in planned MD modules")
        ax.set_title("Candidate simulation budget before immunogenicity testing")
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIG / "fig_md32_candidate_simulation_budget.png", dpi=220)
        fig.savefig(FIG / "fig_md32_candidate_simulation_budget.pdf")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 4.5))
    y = np.arange(len(bridge))
    ax.barh(y, [1] * len(bridge), color="#76b7b2")
    ax.set_yticks(y, bridge["wetlab_assay"])
    ax.set_xlim(0, 1.15)
    ax.set_xticks([])
    ax.set_title("Wetlab immunogenicity assay to simulation-support bridge")
    for i, row in bridge.iterrows():
        ax.text(0.03, i, row["simulation_support"], va="center", fontsize=8, color="black")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIG / "fig_md33_assay_to_simulation_bridge.png", dpi=220)
    fig.savefig(FIG / "fig_md33_assay_to_simulation_bridge.pdf")
    plt.close(fig)

    status = matrix["launch_status"].astype(str).value_counts().head(10)
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.barh(status.index[::-1], status.values[::-1], color="#f28e2b")
    ax.set_xlabel("planned modules")
    ax.set_title("Simulation launch/readiness status")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md34_immunogenicity_simulation_readiness.png", dpi=220)
    fig.savefig(FIG / "fig_md34_immunogenicity_simulation_readiness.pdf")
    plt.close(fig)


def write_commands(matrix: pd.DataFrame) -> None:
    lines = [
        "# Immunogenicity Simulation Escalation Commands",
        "",
        "These commands are launch scaffolds. They do not prove immunogenicity and should be run only when the current queue/resources are intentionally allocated.",
        "",
        "## Refresh status and reports",
        "",
        "```bash",
        "python project/scripts/cross_neo_md/13_watch_runpod_md_and_autoreport.py --once",
        "python project/scripts/cross_neo_md/30_design_immunogenicity_simulation_escalation.py",
        "python project/scripts/cross_neo_md/29_build_high_impact_decision_web.py",
        "```",
        "",
        "## Existing queued/ready OpenMM diverse simulation pack",
        "",
        "```bash",
        "cd project/results/cross_neo_md_audit_2026_05_10/diverse_simulation_plan",
        "bash run_immediate_ready_batch.sh",
        "bash run_after_hmtevv_sync.sh",
        "```",
        "",
        "## Structure-prep refresh before new WT/decoy/TCR jobs",
        "",
        "```bash",
        "python project/scripts/cross_neo_v2/13_prepare_tcr_structure_jobs.py",
        "python project/scripts/cross_neo_md/14_prepare_md_counterfactual_batch.py",
        "python project/scripts/cross_neo_md/17_prepare_counterfactual_structure_prep_jobs.py",
        "python project/scripts/cross_neo_md/25_prepare_diverse_md_launch_pack.py",
        "```",
        "",
        "## Endpoint energy / MMGBSA placeholder after trajectories finish",
        "",
        "```bash",
        "# If AmberTools/MMPBSA.py is installed, run per completed trajectory.",
        "# Otherwise use OpenMM interaction-energy decomposition from parsed frames.",
        "python project/scripts/cross_neo_md/03_analyze_pmhc_contacts.py",
        "python project/scripts/cross_neo_md/04_analyze_tcr_contacts.py",
        "python project/scripts/cross_neo_md/05_counterfactual_md_analysis.py",
        "```",
        "",
        "## Expensive hold-only modules",
        "",
        "- `alchemical_fep_mutant_to_wt_delta_delta_g`: run only for GADGVGKSAL/HMTEVVRHC after WT mapping and stable 10 ns controls.",
        "- `umbrella_or_steered_md_tcr_unbinding_pmf`: run only after stable TCR-pMHC 10 ns and if multimer/TCR-binding assay decision depends on it.",
        "",
    ]
    (OUT / "immunogenicity_simulation_run_commands.md").write_text("\n".join(lines))


def write_report(matrix: pd.DataFrame, bridge: pd.DataFrame, rules: pd.DataFrame, budget: pd.DataFrame) -> None:
    total_ns = matrix["total_requested_ns"].sum()
    active_jobs = int((matrix["total_requested_ns"] > 0).sum())
    blocked = int(matrix["launch_status"].astype(str).str.contains("blocked|hold", case=False, regex=True).sum())
    ready = int(matrix["launch_status"].astype(str).str.contains("ready|queued|manifested|supported", case=False, regex=True).sum())
    lines = [
        "# Immunogenicity Simulation Escalation Report",
        "",
        "## Executive verdict",
        "",
        "Add simulations as a staged pre-wetlab culling layer, not as immunogenicity proof. The best immediate value is WT/decoy-aware pMHC and TCR-pMHC stability for the strict Top-13, with expensive FEP/PMF held for the top two only.",
        "",
        "## Summary",
        "",
        f"- planned simulation/analysis modules: {len(matrix)}",
        f"- modules with explicit MD length: {active_jobs}",
        f"- requested short-MD trajectory length: {total_ns:.1f} ns",
        f"- ready/manifested/supported modules: {ready}",
        f"- blocked or hold modules: {blocked}",
        "",
        "## Highest-impact immediate actions",
        "",
        "1. Keep the completed/synced HMTEVVRHC 10 ns evidence locked in the reports and use it as the MD-strong case study.",
        "2. Run/monitor the already prepared diverse pack for GADGVGKSAL mutant TCR-pMHC replicates and same-HLA positive controls.",
        "3. Curate WT peptides and anchor-preserved decoys for Top-4 candidates before spending more GPU on specificity claims.",
        "4. Add endpoint interaction-energy/MMGBSA only after replicate trajectories exist.",
        "5. Reserve FEP or umbrella/steered MD for GADGVGKSAL and HMTEVVRHC only if WT/decoy 10 ns controls are clean.",
        "",
        "## Assay bridge",
        "",
        bridge.to_markdown(index=False),
        "",
        "## Go/no-go rules",
        "",
        rules.to_markdown(index=False),
        "",
        "## Candidate budget",
        "",
        budget.to_markdown(index=False) if not budget.empty else "No explicit-MD budget rows.",
        "",
        "## Claim boundary",
        "",
        "Simulation can support pMHC stability, TCR-interface plausibility, and mutant-over-WT/decoy specificity hypotheses. It cannot prove IFN-gamma release, T-cell activation, killing, clinical response, or universal immunogenicity.",
        "",
    ]
    (OUT / "IMMUNOGENICITY_SIMULATION_ESCALATION_REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    top13 = read_tsv(PKG / "no_false_positive_top13_candidates.tsv")
    if top13.empty:
        raise FileNotFoundError(PKG / "no_false_positive_top13_candidates.tsv")
    matrix = build_matrix(top13)
    bridge = build_bridge()
    rules = build_go_no_go_rules()
    budget = summarize_budget(matrix)
    matrix.to_csv(OUT / "immunogenicity_simulation_matrix.tsv", sep="\t", index=False)
    bridge.to_csv(OUT / "wetlab_assay_simulation_bridge.tsv", sep="\t", index=False)
    rules.to_csv(OUT / "simulation_go_no_go_rules.tsv", sep="\t", index=False)
    budget.to_csv(OUT / "candidate_simulation_budget.tsv", sep="\t", index=False)
    make_figures(matrix, bridge, budget)
    write_commands(matrix)
    write_report(matrix, bridge, rules, budget)
    summary = {
        "planned_modules": int(len(matrix)),
        "explicit_md_modules": int((matrix["total_requested_ns"] > 0).sum()),
        "total_requested_short_md_ns": float(matrix["total_requested_ns"].sum()),
        "blocked_or_hold_modules": int(matrix["launch_status"].astype(str).str.contains("blocked|hold", case=False, regex=True).sum()),
        "boundary": "simulation-support for immunogenicity testing, not immunogenicity proof",
    }
    (OUT / "immunogenicity_simulation_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


EXISTING_JOBS = build_existing_job_index()


if __name__ == "__main__":
    main()
