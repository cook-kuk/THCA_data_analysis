#!/usr/bin/env python3
"""Build a preclinical validation protocol package for CROSS-Neo leads."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
LEAD = MD_OUT / "two_lead_impact_package"
ASSAY = MD_OUT / "assay_ready_translation_package"
PKG = MD_OUT / "high_impact_decision_package"
OUT = MD_OUT / "preclinical_validation_protocol"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def nfloat(x: object, default: float = 0.0) -> float:
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    leads = read_tsv(LEAD / "two_lead_external_evidence.tsv")
    top13 = read_tsv(PKG / "no_false_positive_top13_candidates.tsv")
    assay = read_tsv(ASSAY / "assay_reagent_order_sheet.tsv")
    if leads.empty or top13.empty:
        raise FileNotFoundError("Required lead/top13 package is missing")
    return leads, top13, assay


def build_readiness_scorecard(leads: pd.DataFrame, top13: pd.DataFrame, assay: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, lead in leads.iterrows():
        top = top13[top13["peptide"].eq(lead["mutant_peptide"])]
        row = top.iloc[0] if not top.empty else pd.Series(dtype=object)
        reagent_rows = assay[assay["lead"].eq(lead["lead"])]
        wt_ready = bool(str(lead["wildtype_peptide"]).strip())
        p0_reagents = int(reagent_rows["priority"].astype(str).eq("P0").sum()) if not reagent_rows.empty else 0
        external = 1.0 if lead["external_support_level"] == "VERY_HIGH" else 0.85
        md = min(nfloat(lead["crossneo_md_score"]), 1.0)
        paired = min(nfloat(lead["paired_tcr_evidence_count"]) / 20.0, 1.0)
        controls = 1.0 if wt_ready and bool(str(lead["anchor_preserved_decoy"]).strip()) else 0.5
        reagent = min(p0_reagents / 4.0, 1.0)
        decision = 0.25 * external + 0.25 * md + 0.20 * paired + 0.20 * controls + 0.10 * reagent
        if decision >= 0.85:
            tier = "VALIDATION_P0_FLAGSHIP"
        elif decision >= 0.70:
            tier = "VALIDATION_P1_READY"
        else:
            tier = "VALIDATION_P2_CURATION"
        rows.append(
            {
                "lead": lead["lead"],
                "mutant_peptide": lead["mutant_peptide"],
                "hla": lead["hla"],
                "wildtype_peptide": lead["wildtype_peptide"],
                "decoy_peptide": lead["anchor_preserved_decoy"],
                "external_score": external,
                "md_score": md,
                "paired_tcr_score": paired,
                "control_readiness_score": controls,
                "reagent_readiness_score": reagent,
                "overall_readiness_score": decision,
                "validation_tier": tier,
                "recommended_next_action": "launch WT/decoy-controlled HLA stability + multimer/activation assays",
                "blocked_claim": "immunogenicity until activation exceeds WT/decoy controls",
                "source_row_id": row.get("row_id", ""),
            }
        )
    return pd.DataFrame(rows)


def build_validation_gates(leads: pd.DataFrame) -> pd.DataFrame:
    gate_rows = []
    gate_defs = [
        (
            1,
            "identity_and_reagent_qc",
            "Mutant, WT, and decoy peptide identity/purity confirmed; HLA context confirmed.",
            "Re-order or hold candidate until matched reagents are clean.",
            "No biological claim.",
        ),
        (
            2,
            "hla_presentation_gate",
            "Mutant pHLA stability/binding is measurable and assay background is controlled; WT/decoy included.",
            "Hold T-cell assays if mutant fails HLA presentation plausibility.",
            "Presentation plausibility.",
        ),
        (
            3,
            "tcr_or_multimer_gate",
            "Mutant pHLA binding/staining exceeds WT, decoy, and irrelevant pHLA at matched conditions.",
            "Downgrade to pMHC-only case if WT/decoy is equal or stronger.",
            "Recognition plausibility in tested TCR/cell context.",
        ),
        (
            4,
            "activation_gate",
            "IFN-gamma ELISpot/ICS response to mutant exceeds WT and decoy controls with replicate consistency.",
            "No immunogenicity claim if activation is absent or WT-like.",
            "Assay-specific immunogenicity.",
        ),
        (
            5,
            "functional_specificity_gate",
            "Killing assay is positive only in matched HLA/mutation target context and negative for WT/decoy controls.",
            "No cytotoxicity claim if target context is mismatched or controls are positive.",
            "Functional cytotoxicity in the tested system.",
        ),
    ]
    for _, lead in leads.iterrows():
        for order, gate, go, hold, claim in gate_defs:
            gate_rows.append(
                {
                    "lead": lead["lead"],
                    "mutant_peptide": lead["mutant_peptide"],
                    "hla": lead["hla"],
                    "gate_order": order,
                    "gate": gate,
                    "go_criterion": go,
                    "hold_or_no_go_action": hold,
                    "claim_allowed_if_pass": claim,
                    "claim_forbidden_if_fail": "Do not claim immunogenicity, clinical utility, or vaccine efficacy.",
                }
            )
    return pd.DataFrame(gate_rows)


def build_failure_actions(leads: pd.DataFrame) -> pd.DataFrame:
    failure_defs = [
        ("mutant_hla_binding_low", "Check peptide quality and HLA allele; run pMHC structure/MD only as diagnostic; deprioritize wetlab activation."),
        ("wt_equal_or_stronger", "Flag self-reactivity/cross-reactivity risk; do not claim mutant specificity; consider alternate TCR or lower priority."),
        ("decoy_positive", "Flag nonspecific recognition or motif artifact; require alternate decoy and orthogonal assay."),
        ("multimer_positive_activation_negative", "Treat as binding without productive activation; analyze TCR off-rate/kinetic proofreading proxy."),
        ("activation_positive_wt_positive", "Potential self-reactivity; do not advance without dose-response and safety context."),
        ("md_stable_wetlab_negative", "Use as false-positive explanation: structural plausibility is not sufficient for immunogenicity."),
        ("wetlab_positive_md_weak", "Use as model rescue case; re-check structure/template and run replicate/control MD if needed."),
    ]
    rows = []
    for _, lead in leads.iterrows():
        for failure, action in failure_defs:
            rows.append(
                {
                    "lead": lead["lead"],
                    "failure_mode": failure,
                    "interpretation": action,
                    "next_analysis": "update case audit, claim boundary, and threshold dashboard",
                }
            )
    return pd.DataFrame(rows)


def build_sample_material_request(leads: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, lead in leads.iterrows():
        rows += [
            {
                "lead": lead["lead"],
                "material": "HLA-matched donor PBMC or patient T cells",
                "requirement": f"Match {lead['hla']}; record donor HLA typing and cancer/context metadata.",
                "why_needed": "Activation assays require a cellular context; simulation cannot replace this.",
                "priority": "P0",
            },
            {
                "lead": lead["lead"],
                "material": "paired TCR clone or TCR-transduced reporter if available",
                "requirement": "Use mutant, WT, decoy, irrelevant peptide, and known positive pHLA controls.",
                "why_needed": "Separates cognate recognition from bulk PBMC noise.",
                "priority": "P0" if lead["mutant_peptide"] == "HMTEVVRHC" else "P1",
            },
            {
                "lead": lead["lead"],
                "material": "target cells for killing assay",
                "requirement": "Matched HLA, mutation expression, antigen processing/presentation evidence.",
                "why_needed": "Killing assay is only meaningful after binding and activation gates pass.",
                "priority": "P2_followup",
            },
        ]
    return pd.DataFrame(rows)


def make_figures(score: pd.DataFrame, gates: pd.DataFrame, failures: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    metrics = [
        "external_score",
        "md_score",
        "paired_tcr_score",
        "control_readiness_score",
        "reagent_readiness_score",
        "overall_readiness_score",
    ]
    heat = score.set_index("lead")[metrics]
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    im = ax.imshow(heat.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(heat.index)), heat.index)
    ax.set_xticks(np.arange(len(metrics)), [m.replace("_", "\n") for m in metrics], rotation=0, fontsize=8)
    ax.set_title("Preclinical validation readiness scorecard")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md42_preclinical_readiness_scorecard.png", dpi=220)
    fig.savefig(FIG / "fig_md42_preclinical_readiness_scorecard.pdf")
    plt.close(fig)

    gate_names = gates.drop_duplicates("gate_order").sort_values("gate_order")["gate"].tolist()
    fig, ax = plt.subplots(figsize=(12, 4.4))
    ax.axis("off")
    xs = np.linspace(0.08, 0.92, len(gate_names))
    for i, gate in enumerate(gate_names):
        ax.add_patch(plt.Circle((xs[i], 0.58), 0.065, color="#17324f", ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.58, str(i + 1), ha="center", va="center", color="#f2c46d", fontsize=16, weight="bold")
        ax.text(xs[i], 0.34, gate.replace("_", "\n"), ha="center", va="center", fontsize=8.5)
        if i < len(gate_names) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.08, 0.58), xytext=(xs[i] + 0.08, 0.58), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.84, "Preclinical validation gates: when claims can escalate", ha="center", fontsize=15, weight="bold")
    fig.savefig(FIG / "fig_md43_preclinical_validation_gates.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md43_preclinical_validation_gates.pdf", bbox_inches="tight")
    plt.close(fig)

    top_fail = failures["failure_mode"].value_counts()
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.barh(top_fail.index[::-1], top_fail.values[::-1], color="#e15759")
    ax.set_xlabel("lead-specific action rows")
    ax.set_title("Failure-mode action map")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md44_failure_mode_action_map.png", dpi=220)
    fig.savefig(FIG / "fig_md44_failure_mode_action_map.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.2, 4.4))
    ax.axis("off")
    states = [
        ("Computational\npriority", "Top-13 no-FP\nTCR + MD"),
        ("Assay-ready\ncontrols", "mutant + WT\n+ decoy"),
        ("Preclinical\ngates", "HLA -> TCR\n-> activation"),
        ("Claim\nboundary", "only escalate\nif controls pass"),
    ]
    xs = np.linspace(0.10, 0.90, len(states))
    for i, (title, body) in enumerate(states):
        ax.add_patch(plt.Rectangle((xs[i] - 0.09, 0.46), 0.18, 0.23, fc="#10243a", ec="#395170", lw=1.2))
        ax.text(xs[i], 0.62, title, ha="center", va="center", color="#f2c46d", fontsize=10, weight="bold")
        ax.text(xs[i], 0.51, body, ha="center", va="center", color="white", fontsize=9)
        if i < len(states) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.10, 0.57), xytext=(xs[i] + 0.10, 0.57), arrowprops=dict(arrowstyle="->", lw=2))
    ax.set_title("From computational candidate to preclinical validation protocol", fontsize=15, weight="bold")
    fig.savefig(FIG / "fig_md45_preclinical_translation_flow.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md45_preclinical_translation_flow.pdf", bbox_inches="tight")
    plt.close(fig)


def write_protocol(score: pd.DataFrame, gates: pd.DataFrame, failures: pd.DataFrame, materials: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo Preclinical Validation Protocol",
        "",
        "## Objective",
        "",
        "Convert the two flagship CROSS-Neo candidates into a claim-safe preclinical validation workflow. The protocol separates computational prioritization, structural plausibility, TCR recognition plausibility, and assay-specific immunogenicity.",
        "",
        "## Readiness Scorecard",
        "",
        score.to_markdown(index=False),
        "",
        "## Validation Gates",
        "",
        gates.to_markdown(index=False),
        "",
        "## Sample / Material Request",
        "",
        materials.to_markdown(index=False),
        "",
        "## Failure Mode Actions",
        "",
        failures.to_markdown(index=False),
        "",
        "## Claim Boundary",
        "",
        "A candidate becomes assay-specific immunogenic only after mutant activation exceeds WT and decoy controls. Structural MD, pHLA binding, and multimer binding are supportive layers, not substitutes for activation assays.",
    ]
    (OUT / "PRECLINICAL_VALIDATION_PROTOCOL.md").write_text("\n".join(lines) + "\n")

    collab = [
        "# Collaboration Brief: CROSS-Neo Two-Lead Validation",
        "",
        "We have two prioritized neoantigen candidates ready for controlled wetlab evaluation:",
        "",
        "- KRAS G12D `GADGVGKSAL / HLA-C*08:02`, with WT `GAGGVGKSAL` and decoy `GASGVKGADL`.",
        "- TP53 R175H `HMTEVVRHC / HLA-A*02:01`, with WT `HMTEVVRRC` and decoy `HMTHRVVEC`.",
        "",
        "The computational package includes strict current-label prioritization, paired TCR evidence, published external support, MD audit, reagent order sheets, plate maps, and go/no-go criteria. We are not claiming immunogenicity from computation. The requested experiment is to test whether mutant response exceeds WT and decoy controls in HLA stability, multimer/TCR binding, and IFN-gamma/ICS assays.",
        "",
        "Primary deliverable: determine whether either lead passes activation above WT/decoy controls under matched HLA context.",
    ]
    (OUT / "COLLABORATION_BRIEF.md").write_text("\n".join(collab) + "\n")

    kr = [
        "# Preclinical validation Korean one-page",
        "",
        "## 지금 상태",
        "",
        "두 후보는 이제 computational ranking을 넘어서 preclinical validation protocol까지 내려왔습니다.",
        "",
        "- KRAS G12D: `GADGVGKSAL / HLA-C*08:02`",
        "- TP53 R175H: `HMTEVVRHC / HLA-A*02:01`",
        "",
        "## GO 조건",
        "",
        "mutant가 WT와 decoy보다 HLA stability, TCR/multimer binding, IFN-gamma/ICS에서 강해야 합니다.",
        "",
        "## HOLD/NO-GO",
        "",
        "WT가 같거나 강하면 self-reactivity risk입니다. Decoy가 강하면 motif/nonspecific risk입니다. Activation이 없으면 immunogenicity claim은 금지입니다.",
    ]
    (OUT / "PRECLINICAL_VALIDATION_ONE_PAGE_KR.md").write_text("\n".join(kr) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    leads, top13, assay = load_inputs()
    score = build_readiness_scorecard(leads, top13, assay)
    gates = build_validation_gates(leads)
    failures = build_failure_actions(leads)
    materials = build_sample_material_request(leads)
    score.to_csv(OUT / "preclinical_readiness_scorecard.tsv", sep="\t", index=False)
    gates.to_csv(OUT / "preclinical_validation_gates.tsv", sep="\t", index=False)
    failures.to_csv(OUT / "preclinical_failure_mode_actions.tsv", sep="\t", index=False)
    materials.to_csv(OUT / "sample_material_request.tsv", sep="\t", index=False)
    make_figures(score, gates, failures)
    write_protocol(score, gates, failures, materials)
    summary = {
        "n_leads": int(len(score)),
        "n_gates": int(len(gates)),
        "n_failure_actions": int(len(failures)),
        "figures": [
            "fig_md42_preclinical_readiness_scorecard",
            "fig_md43_preclinical_validation_gates",
            "fig_md44_failure_mode_action_map",
            "fig_md45_preclinical_translation_flow",
        ],
        "boundary": "preclinical validation protocol, not immunogenicity proof",
    }
    (OUT / "preclinical_validation_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
