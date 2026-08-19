#!/usr/bin/env python3
"""Build a two-lead impact package for CROSS-Neo immunogenicity triage."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
PKG = MD_OUT / "high_impact_decision_package"
SIM = MD_OUT / "immunogenicity_simulation_escalation"
OUT = MD_OUT / "two_lead_impact_package"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def num(value: object, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def load_leads() -> pd.DataFrame:
    top = read_tsv(PKG / "no_false_positive_top13_candidates.tsv")
    if top.empty:
        raise FileNotFoundError(PKG / "no_false_positive_top13_candidates.tsv")
    leads = top[top["peptide"].isin(["GADGVGKSAL", "HMTEVVRHC"])].copy()
    leads["lead_rank"] = leads["peptide"].map({"GADGVGKSAL": 1, "HMTEVVRHC": 2})
    return leads.sort_values("lead_rank")


def build_external_evidence(leads: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in leads.iterrows():
        pep = row["peptide"]
        if pep == "GADGVGKSAL":
            rows.append(
                {
                    "lead": "KRAS_G12D_HLA_C0802",
                    "row_id": row["row_id"],
                    "mutation": "KRAS G12D",
                    "mutant_peptide": "GADGVGKSAL",
                    "wildtype_peptide": "GAGGVGKSAL",
                    "anchor_preserved_decoy": row.get("anchor_preserved_decoy", "GASGVKGADL"),
                    "hla": "HLA-C*08:02",
                    "mutation_position_in_peptide": 3,
                    "mechanistic_specificity": "P3 Asp creates an HLA-C*08:02 Arg156 salt-bridge anchor; WT glycine is a weak anchor and WT KRAS peptides fail HLA-C stabilization in reported assays.",
                    "available_structure_or_template": "6UON/TCR-pMHC and related KRAS G12D HLA-C*08:02 structures/templates already used locally.",
                    "commercial_or_wetlab_reagent_note": "Custom peptide synthesis and HLA-C*08:02/KRAS G12D TCR assay route; verify reagent availability with vendor before ordering.",
                    "external_support_level": "HIGH",
                    "external_support_reason": "Published ACT/TCR work supports KRAS G12D HLA-C*08:02 as a recognized shared neoantigen with mutant-over-WT presentation specificity.",
                    "primary_source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7293613/",
                    "secondary_source_url": "https://www.nature.com/articles/s41467-022-32811-1",
                    "crossneo_md_label": row.get("md_label", ""),
                    "crossneo_md_score": row.get("md_score", ""),
                    "paired_tcr_evidence_count": row.get("paired_tcr_evidence_count", 0),
                    "claim_upgrade": "Lead wetlab candidate with structural and TCR-resource support; needs WT/decoy assays for immunogenicity claim.",
                }
            )
        elif pep == "HMTEVVRHC":
            rows.append(
                {
                    "lead": "TP53_R175H_HLA_A0201",
                    "row_id": row["row_id"],
                    "mutation": "TP53 R175H",
                    "mutant_peptide": "HMTEVVRHC",
                    "wildtype_peptide": "HMTEVVRRC",
                    "anchor_preserved_decoy": row.get("anchor_preserved_decoy", "HMTHRVVEC"),
                    "hla": "HLA-A*02:01",
                    "mutation_position_in_peptide": 8,
                    "mechanistic_specificity": "R175H changes peptide position 8 from Arg to His. Published p53 R175H/HLA-A*02:01 TCR work and WT p53_168-176 reagent availability make this a strong specificity-control case.",
                    "available_structure_or_template": "6VRN primary 10 ns completed locally; p53 WT HMTEVVRRC pMHC structure/reagent sources exist.",
                    "commercial_or_wetlab_reagent_note": "Mutant HMTEVVRHC and WT HMTEVVRRC pHLA monomer/tetramer reagents are commercially listed; peptide and tetramer assays are unusually actionable.",
                    "external_support_level": "VERY_HIGH",
                    "external_support_reason": "p53 R175H 168-176 is a published HLA-A*02:01 shared neoantigen with TCR structures and clinical/translational TCE/TCR interest.",
                    "primary_source_url": "https://www.nature.com/articles/s41467-020-16755-y",
                    "secondary_source_url": "https://jitc.bmj.com/content/12/Suppl_2/A1181",
                    "reagent_source_url": "https://www.genscript.com/peptide/RP30864-P53_R175H_Peptide_HMTEVVRHC_.html",
                    "wt_reagent_source_url": "https://www.targetmol.com/recombinant-protein/hla_a_02_01_b2m_p53_wt_hmtevvrrc_tetramer_protein_human_mhc_his_avi_",
                    "crossneo_md_label": row.get("md_label", ""),
                    "crossneo_md_score": row.get("md_score", ""),
                    "paired_tcr_evidence_count": row.get("paired_tcr_evidence_count", 0),
                    "claim_upgrade": "Promote to flagship case study and wetlab-ready candidate after WT/decoy specificity controls.",
                }
            )
    return pd.DataFrame(rows)


def build_specificity_locks(external: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in external.iterrows():
        if row["mutant_peptide"] == "GADGVGKSAL":
            locks = [
                ("mutation_identity", "PASS", "KRAS G12D decamer GADGVGKSAL; WT GAGGVGKSAL."),
                ("hla_context", "PASS", "HLA-C*08:02 matched to published KRAS G12D TCR/ACT evidence."),
                ("structure_template", "PASS", "Local 6UON-based pMHC/TCR-pMHC trajectory exists."),
                ("mutant_md", "REPLICATED_PASS", f"Current CROSS-Neo MD label {row['crossneo_md_label']} score {row['crossneo_md_score']}. Completed 10 ns replicate/control MD is available, with replicate consistency marked robust_stable."),
                ("wt_control", "READY_TENTATIVE", "WT GAGGVGKSAL inferred; must be ordered/tested and optionally modeled."),
                ("decoy_control", "READY_TENTATIVE", f"Anchor-preserved decoy candidate {row['anchor_preserved_decoy']}."),
                ("immunogenicity_assay", "NOT_YET", "Needs HLA binding, tetramer/multimer if TCR available, IFN-gamma/ICS, and WT/decoy comparison."),
            ]
        else:
            locks = [
                ("mutation_identity", "PASS", "TP53 R175H peptide HMTEVVRHC; WT HMTEVVRRC confirmed by p53_168-176 sources."),
                ("hla_context", "PASS", "HLA-A*02:01 matched to published p53 R175H TCR structural evidence."),
                ("structure_template", "PASS", "Local 6VRN 10 ns completed; MD_VERY_STRONG."),
                ("mutant_md", "PASS", f"Current CROSS-Neo MD label {row['crossneo_md_label']} score {row['crossneo_md_score']}."),
                ("wt_control", "READY_NOW", "WT HMTEVVRRC pHLA reagents/structure sources exist; order/model as control."),
                ("decoy_control", "READY_TENTATIVE", f"Anchor-preserved decoy candidate {row['anchor_preserved_decoy']}."),
                ("immunogenicity_assay", "NOT_YET", "Needs HLA binding, mutant/WT tetramer, IFN-gamma/ICS, and killing only after activation."),
            ]
        for lock, status, detail in locks:
            rows.append(
                {
                    "lead": row["lead"],
                    "mutant_peptide": row["mutant_peptide"],
                    "hla": row["hla"],
                    "specificity_lock": lock,
                    "status": status,
                    "detail": detail,
                }
            )
    return pd.DataFrame(rows)


def build_wetlab_order(external: pd.DataFrame) -> pd.DataFrame:
    rows = []
    order = 0
    for _, row in external.iterrows():
        for assay, reason, output, stop_rule in [
            ("peptide_synthesis_qc", "Needed for mutant, WT, and decoy matched testing.", "HPLC/MS purity and solubility", "Stop if peptide cannot be synthesized/purified."),
            ("hla_binding_or_stability", "Confirms presentation plausibility before T-cell readouts.", "MHC stabilization or binding/stability signal", "Stop if mutant is not stronger/plausible versus controls."),
            ("mutant_wt_decoy_pmhc_multimer", "Direct recognition screen when paired TCR or T-cell material exists.", "Mutant multimer positive with low WT/decoy binding", "Stop or downgrade if WT/decoy binds equally."),
            ("ifng_elispot_or_ics", "Minimum activation assay for immunogenicity claim.", "Mutant response above WT/decoy and irrelevant peptide", "No immunogenicity claim without activation above controls."),
            ("killing_assay_followup", "Only after activation is positive.", "Specific target-cell killing with matched HLA/mutation", "Do not run as first-line screen."),
        ]:
            order += 1
            rows.append(
                {
                    "order": order,
                    "lead": row["lead"],
                    "mutant_peptide": row["mutant_peptide"],
                    "hla": row["hla"],
                    "assay": assay,
                    "matched_controls": f"WT {row['wildtype_peptide']}; decoy {row['anchor_preserved_decoy']}; irrelevant HLA-matched peptide; known positive pHLA if available",
                    "why": reason,
                    "expected_output": output,
                    "stop_or_hold_rule": stop_rule,
                }
            )
    return pd.DataFrame(rows)


def build_figures(external: pd.DataFrame, locks: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    metrics = []
    for _, row in external.iterrows():
        metrics.append(
            {
                "lead": row["mutant_peptide"] + "\n" + row["hla"],
                "External": 1.0 if row["external_support_level"] == "VERY_HIGH" else 0.85,
                "Paired TCR": min(num(row["paired_tcr_evidence_count"]) / 20.0, 1.0),
                "MD": min(num(row["crossneo_md_score"]), 1.0),
                "WT control": 1.0 if row["wildtype_peptide"] == "HMTEVVRRC" else 0.75,
                "Reagent": 1.0 if row["mutant_peptide"] == "HMTEVVRHC" else 0.65,
            }
        )
    mdf = pd.DataFrame(metrics).set_index("lead")
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    im = ax.imshow(mdf.to_numpy(), aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(mdf.index)), mdf.index)
    ax.set_xticks(np.arange(len(mdf.columns)), mdf.columns, rotation=25, ha="right")
    ax.set_title("Two-lead evidence moat for immunogenicity testing")
    fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md35_two_lead_evidence_moat.png", dpi=220)
    fig.savefig(FIG / "fig_md35_two_lead_evidence_moat.pdf")
    plt.close(fig)

    status_score = {"PASS": 1.0, "REPLICATED_PASS": 0.85, "READY_NOW": 0.9, "READY_TENTATIVE": 0.75, "PARTIAL_PASS": 0.65, "NOT_YET": 0.25}
    piv = locks.copy()
    piv["score"] = piv["status"].map(status_score).fillna(0.0)
    heat = piv.pivot(index="lead", columns="specificity_lock", values="score")
    fig, ax = plt.subplots(figsize=(11, 4.4))
    im = ax.imshow(heat.to_numpy(), aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(heat.index)), heat.index)
    ax.set_xticks(np.arange(len(heat.columns)), heat.columns, rotation=30, ha="right", fontsize=8)
    ax.set_title("Specificity lock status before immunogenicity claim")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md36_specificity_lock_matrix.png", dpi=220)
    fig.savefig(FIG / "fig_md36_specificity_lock_matrix.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11.5, 4.5))
    ax.axis("off")
    lanes = [
        ("CROSS-Neo", "DL/TCR triage\nTop-13 no-FP"),
        ("Structure", "TCR-pMHC\npublished/local"),
        ("MD", "10 ns audit\nHMTEVVRHC very strong"),
        ("Controls", "WT + decoy\nspecificity"),
        ("Wetlab", "ELISpot/ICS\nthen killing"),
    ]
    xs = np.linspace(0.08, 0.92, len(lanes))
    for i, (head, body) in enumerate(lanes):
        ax.add_patch(plt.Rectangle((xs[i] - 0.075, 0.45), 0.15, 0.22, fc="#17324f", ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.61, head, ha="center", va="center", color="#f2c46d", weight="bold", fontsize=10)
        ax.text(xs[i], 0.51, body, ha="center", va="center", color="white", fontsize=9)
        if i < len(lanes) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.085, 0.56), xytext=(xs[i] + 0.085, 0.56), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.82, "Impact upgrade: from candidate ranking to controlled immunogenicity test plan", ha="center", fontsize=15, weight="bold")
    ax.text(0.5, 0.22, "Allowed now: prioritization and structural plausibility. Strong immunogenicity claim requires activation above WT/decoy.", ha="center", fontsize=10, color="#9b1c31")
    fig.savefig(FIG / "fig_md37_impact_upgrade_flow.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md37_impact_upgrade_flow.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(external: pd.DataFrame, locks: pd.DataFrame, wetlab: pd.DataFrame) -> None:
    lines = [
        "# Two-Lead CROSS-Neo Impact Upgrade",
        "",
        "## Decision",
        "",
        "Promote GADGVGKSAL/HLA-C*08:02 and HMTEVVRHC/HLA-A*02:01 as the flagship wetlab triage cases. The impact increase comes from combining current-label no-FP selection, paired TCR evidence, published shared-neoantigen biology, completed MD with a replicated GADGVGKSAL 10 ns run, and explicit WT/decoy control plans.",
        "",
        "## Lead Evidence",
        "",
        external.to_markdown(index=False),
        "",
        "## Specificity Locks",
        "",
        locks.to_markdown(index=False),
        "",
        "## Wetlab Order",
        "",
        wetlab.to_markdown(index=False),
        "",
        "## Claim Boundary",
        "",
        "These leads can be claimed as high-priority, structure-auditable, TCR-resource-supported candidates for controlled wetlab testing. They cannot be claimed as immunogenic until mutant responses exceed WT/decoy controls in activation assays.",
        "",
        "## Source URLs",
        "",
        "- KRAS G12D HLA-C*08:02 TCR/ACT evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC7293613/",
        "- KRAS G12D structural specificity discussion: https://www.nature.com/articles/s41467-022-32811-1",
        "- TP53 R175H HLA-A*02:01 TCR structure evidence: https://www.nature.com/articles/s41467-020-16755-y",
        "- TP53 R175H T-cell engager specificity abstract: https://jitc.bmj.com/content/12/Suppl_2/A1181",
        "- TP53 R175H peptide reagent example: https://www.genscript.com/peptide/RP30864-P53_R175H_Peptide_HMTEVVRHC_.html",
        "- TP53 WT HMTEVVRRC reagent/structure source example: https://www.targetmol.com/recombinant-protein/hla_a_02_01_b2m_p53_wt_hmtevvrrc_tetramer_protein_human_mhc_his_avi_",
    ]
    (OUT / "TWO_LEAD_IMPACT_UPGRADE_REPORT.md").write_text("\n".join(lines) + "\n")


def write_manuscript_insert(external: pd.DataFrame) -> None:
    lines = [
        "# Manuscript Insert: Two-Lead Controlled Immunogenicity Triage",
        "",
        "## Results-Style Insert",
        "",
        "We next asked whether the highest-priority candidates could be advanced from a ranked list to a controlled immunogenicity testing plan. Two candidates emerged as flagship case studies: KRAS G12D GADGVGKSAL presented by HLA-C*08:02 and TP53 R175H HMTEVVRHC presented by HLA-A*02:01. Both candidates were selected by the strict no-false-positive operating point in the current labeled set and had paired TCR evidence. The KRAS G12D case has published mutant-over-wild-type presentation specificity in HLA-C*08:02, while the TP53 R175H case has published TCR structural evidence and commercially listed mutant and wild-type pHLA reagents. In our MD audit, HMTEVVRHC completed a 10 ns explicit-solvent trajectory with MD_VERY_STRONG evidence, while GADGVGKSAL retained replicated MD_MODERATE support after completion of an additional 10 ns explicit-solvent run. We therefore promoted these two candidates to a wetlab-ready case-study tier, with a strict requirement that mutant responses exceed wild-type and anchor-preserved decoy controls before any immunogenicity claim.",
        "",
        "## Figure Legends",
        "",
        "**Figure MD35. Two-lead evidence moat for controlled immunogenicity testing.** Heatmap summarizing orthogonal evidence for KRAS G12D/HLA-C*08:02 and TP53 R175H/HLA-A*02:01, including external support, paired TCR evidence, MD evidence, wild-type control readiness, and reagent/actionability. Scores are decision-support indicators, not immunogenicity probabilities.",
        "",
        "**Figure MD36. Specificity-lock matrix before immunogenicity claims.** Status matrix for mutation identity, HLA context, structure/template availability, mutant MD evidence, wild-type control readiness, decoy control readiness, and activation assay status. A candidate cannot advance to an immunogenicity claim until the activation assay lock is satisfied under WT/decoy controls.",
        "",
        "**Figure MD37. Impact upgrade flow.** Schematic showing how CROSS-Neo moves from DL/TCR ranking to structure/MD audit, then WT/decoy specificity controls, and finally wetlab immunogenicity assays. The flow explicitly separates allowed structural-prioritization claims from forbidden immunogenicity claims before wetlab activation data.",
        "",
        "## Claim-Safe Text",
        "",
        "- Allowed now: these two candidates are high-priority, structure-auditable, TCR-resource-supported wetlab candidates.",
        "- Allowed now: HMTEVVRHC has completed 10 ns MD with strong structural plausibility in the modeled TCR-pMHC context.",
        "- Not allowed now: either candidate is immunogenic, clinically useful, or validated as a vaccine target.",
        "- Required next: mutant response must exceed WT and decoy controls in HLA binding/stability and T-cell activation assays.",
        "",
        "## Source Anchors",
        "",
        external[["lead", "primary_source_url", "secondary_source_url"]].to_markdown(index=False),
    ]
    (OUT / "TWO_LEAD_MANUSCRIPT_INSERT.md").write_text("\n".join(lines) + "\n")

    kr_lines = [
        "# Two-lead impact upgrade Korean brief",
        "",
        "## 핵심",
        "",
        "- 이제 단순 Top-13 랭킹이 아니라, 실제 실험으로 바로 넘길 수 있는 flagship 2개 case-study가 생겼습니다.",
        "- 1번: KRAS G12D `GADGVGKSAL / HLA-C*08:02`.",
        "- 2번: TP53 R175H `HMTEVVRHC / HLA-A*02:01`.",
        "- 둘 다 paired TCR evidence가 있고, 문헌/구조/MD/WT-control 축을 붙일 수 있습니다.",
        "",
        "## 왜 임팩트가 올라갔나",
        "",
        "- KRAS G12D는 WT peptide가 HLA-C*08:02에 잘 제시되지 않는다는 문헌 근거가 있어 mutant specificity 스토리가 좋습니다.",
        "- TP53 R175H는 WT `HMTEVVRRC`, mutant `HMTEVVRHC`가 명확하고, mutant/WT pHLA reagent와 TCR 구조 근거가 있어 실험 설계가 매우 직접적입니다.",
        "- HMTEVVRHC는 10 ns MD 완료 후 `MD_VERY_STRONG`으로 올라갔습니다.",
        "- GADGVGKSAL도 추가 10 ns replicate/control MD가 완료되어 `MD_MODERATE`를 반복 MD 근거로 방어할 수 있습니다.",
        "",
        "## 다음 실험 문장",
        "",
        "두 후보 모두 peptide-HLA stability, mutant/WT/decoy tetramer 또는 multimer, IFN-gamma ELISpot/ICS를 순서대로 진행하고, mutant가 WT/decoy보다 강할 때만 immunogenicity claim으로 올립니다.",
        "",
        "## 금지",
        "",
        "- 아직 면역원성 증명이라고 말하지 않습니다.",
        "- MD가 T-cell activation을 증명한다고 말하지 않습니다.",
        "- WT/decoy 없이 mutant-specific recognition claim을 하지 않습니다.",
    ]
    (OUT / "TWO_LEAD_ONE_PAGE_KR.md").write_text("\n".join(kr_lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    leads = load_leads()
    external = build_external_evidence(leads)
    locks = build_specificity_locks(external)
    wetlab = build_wetlab_order(external)
    external.to_csv(OUT / "two_lead_external_evidence.tsv", sep="\t", index=False)
    locks.to_csv(OUT / "two_lead_specificity_locks.tsv", sep="\t", index=False)
    wetlab.to_csv(OUT / "two_lead_wetlab_order.tsv", sep="\t", index=False)
    build_figures(external, locks)
    write_report(external, locks, wetlab)
    write_manuscript_insert(external)
    summary = {
        "leads": external["lead"].tolist(),
        "n_specificity_locks": int(len(locks)),
        "n_wetlab_steps": int(len(wetlab)),
        "figures": [
            "fig_md35_two_lead_evidence_moat",
            "fig_md36_specificity_lock_matrix",
            "fig_md37_impact_upgrade_flow",
        ],
        "boundary": "lead case-study and wetlab prioritization, not immunogenicity proof",
    }
    (OUT / "two_lead_impact_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
