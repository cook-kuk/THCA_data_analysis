#!/usr/bin/env python3
"""Build an assay-ready translation package for the two CROSS-Neo lead cases."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
LEAD = MD_OUT / "two_lead_impact_package"
PKG = MD_OUT / "high_impact_decision_package"
OUT = MD_OUT / "assay_ready_translation_package"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def load_leads() -> pd.DataFrame:
    leads = read_tsv(LEAD / "two_lead_external_evidence.tsv")
    if leads.empty:
        raise FileNotFoundError(LEAD / "two_lead_external_evidence.tsv")
    return leads


def build_reagent_order_sheet(leads: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in leads.iterrows():
        lead = row["lead"]
        mutant = row["mutant_peptide"]
        wt = row["wildtype_peptide"]
        decoy = row["anchor_preserved_decoy"]
        hla = row["hla"]
        rows += [
            {
                "lead": lead,
                "item_type": "synthetic_peptide",
                "item_name": f"{mutant} mutant peptide",
                "sequence": mutant,
                "hla_context": hla,
                "minimum_spec": ">=90-95% purity by HPLC/MS; endotoxin-aware handling if cell assay",
                "purpose": "mutant test antigen",
                "priority": "P0",
                "ordering_note": "Use aliquoted lyophilized peptide; verify solubility and DMSO/water compatibility before assay.",
            },
            {
                "lead": lead,
                "item_type": "synthetic_peptide",
                "item_name": f"{wt} wild-type peptide",
                "sequence": wt,
                "hla_context": hla,
                "minimum_spec": "same purity/formulation as mutant peptide",
                "purpose": "self-reactivity and specificity control",
                "priority": "P0",
                "ordering_note": "Must be tested at matched molar concentrations against mutant peptide.",
            },
            {
                "lead": lead,
                "item_type": "synthetic_peptide",
                "item_name": f"{decoy} anchor-preserved decoy peptide",
                "sequence": decoy,
                "hla_context": hla,
                "minimum_spec": "same purity/formulation as mutant peptide",
                "purpose": "sequence-specificity control",
                "priority": "P0",
                "ordering_note": "Use as negative control; do not interpret if HLA binding is much weaker unless that is the intended stress test.",
            },
            {
                "lead": lead,
                "item_type": "pHLA_monomer_or_tetramer",
                "item_name": f"{mutant}/{hla} pHLA reagent",
                "sequence": mutant,
                "hla_context": hla,
                "minimum_spec": "biotinylated monomer and/or PE/APC tetramer where available",
                "purpose": "TCR/multimer binding screen",
                "priority": "P0" if mutant == "HMTEVVRHC" else "P1",
                "ordering_note": "For TP53 R175H, commercial mutant pHLA reagents are listed; for KRAS/HLA-C verify vendor availability or refold in-house.",
            },
            {
                "lead": lead,
                "item_type": "pHLA_monomer_or_tetramer",
                "item_name": f"{wt}/{hla} pHLA reagent",
                "sequence": wt,
                "hla_context": hla,
                "minimum_spec": "matched fluorophore/format to mutant pHLA",
                "purpose": "WT pHLA binding specificity",
                "priority": "P0",
                "ordering_note": "Required before claiming mutant-specific recognition.",
            },
        ]
    return pd.DataFrame(rows)


def build_plate_map(leads: pd.DataFrame) -> pd.DataFrame:
    rows = []
    wells = [f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)]
    concentrations = ["high", "mid", "low"]
    idx = 0
    for _, row in leads.iterrows():
        for assay in ["HLA_stability", "multimer_or_tetramer", "IFNg_ELISpot_or_ICS"]:
            for antigen_type, sequence in [
                ("mutant", row["mutant_peptide"]),
                ("wildtype", row["wildtype_peptide"]),
                ("decoy", row["anchor_preserved_decoy"]),
                ("irrelevant", "HLA_MATCHED_IRRELEVANT"),
            ]:
                for conc in concentrations:
                    for rep in [1, 2]:
                        if idx >= len(wells):
                            idx = 0
                        rows.append(
                            {
                                "plate": "plate_1" if assay != "IFNg_ELISpot_or_ICS" else "plate_2",
                                "well": wells[idx],
                                "lead": row["lead"],
                                "assay": assay,
                                "antigen_type": antigen_type,
                                "sequence": sequence,
                                "hla": row["hla"],
                                "concentration_level": conc,
                                "replicate": rep,
                                "readout": {
                                    "HLA_stability": "surface HLA or pHLA stability signal",
                                    "multimer_or_tetramer": "TCR+ multimer frequency or MFI",
                                    "IFNg_ELISpot_or_ICS": "spot count or IFN-gamma+ T-cell frequency",
                                }[assay],
                            }
                        )
                        idx += 1
    return pd.DataFrame(rows)


def build_go_no_go_thresholds(leads: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in leads.iterrows():
        rows += [
            {
                "lead": row["lead"],
                "gate": "HLA presentation",
                "go_threshold": "mutant HLA stability/binding is clearly above no-peptide and not worse than assay positive control; WT/decoy context recorded",
                "hold_or_no_go": "mutant fails HLA stabilization or assay is dominated by no-peptide/background signal",
                "claim_allowed_if_go": "presentation plausibility",
            },
            {
                "lead": row["lead"],
                "gate": "TCR/pHLA multimer",
                "go_threshold": "mutant pHLA binds TCR/T cells above WT, decoy, and irrelevant controls at matched concentration",
                "hold_or_no_go": "WT or decoy is equal/higher than mutant; nonspecific pHLA staining; no paired TCR/cell material",
                "claim_allowed_if_go": "recognition plausibility in tested TCR/cell context",
            },
            {
                "lead": row["lead"],
                "gate": "activation",
                "go_threshold": "mutant IFN-gamma/ICS response exceeds WT and decoy controls with replicate consistency",
                "hold_or_no_go": "activation is absent, WT-like, decoy-like, or only seen at nonphysiologic peptide dose",
                "claim_allowed_if_go": "assay-specific immunogenicity in the tested context",
            },
            {
                "lead": row["lead"],
                "gate": "killing",
                "go_threshold": "specific target-cell killing requires matched HLA, mutation expression, antigen processing/presentation, and WT/decoy negative controls",
                "hold_or_no_go": "do not run or claim killing before HLA and activation gates pass",
                "claim_allowed_if_go": "functional cytotoxicity in matched experimental system",
            },
        ]
    return pd.DataFrame(rows)


def build_reviewer_moat(leads: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "reviewer_attack": "The model was optimized on current labels.",
                "defense_asset": "Top-13 label-aware result is framed as a decision-support operating point; source-heldout/external validation is separated.",
                "where_addressed": "high-impact decision page; threshold optimizer; claim ladder",
            },
            {
                "reviewer_attack": "MD cannot prove immunogenicity.",
                "defense_asset": "MD is explicitly used only for structural plausibility and control selection; activation assays are required.",
                "where_addressed": "MD claim boundary; assay-ready go/no-go thresholds; fig_md37",
            },
            {
                "reviewer_attack": "TCR-aware prediction is not general because most datasets lack paired TCR.",
                "defense_asset": "TCR is a case-study/expert branch; main pMHC ranking is kept separate.",
                "where_addressed": "TCR extension decision report; two-lead package",
            },
            {
                "reviewer_attack": "WT cross-reactivity was not tested.",
                "defense_asset": "WT controls are now explicit: KRAS WT GAGGVGKSAL and TP53 WT HMTEVVRRC; no immunogenicity claim before WT/decoy gates.",
                "where_addressed": "specificity lock matrix; reagent order sheet; wetlab plate map",
            },
            {
                "reviewer_attack": "These are literature-known antigens, not novel discovery.",
                "defense_asset": "Use them as anchor case studies to validate the pipeline's triage logic, while keeping discovery claims for independent candidates after validation.",
                "where_addressed": "two-lead impact report; manuscript insert",
            },
        ]
    )


def make_figures(order: pd.DataFrame, thresholds: pd.DataFrame, moat: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    count = order.groupby(["lead", "item_type"]).size().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    count.plot(kind="bar", stacked=True, ax=ax, color=["#4e79a7", "#59a14f"])
    ax.set_ylabel("items")
    ax.set_title("Assay-ready reagent order sheet by lead")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md38_reagent_order_sheet.png", dpi=220)
    fig.savefig(FIG / "fig_md38_reagent_order_sheet.pdf")
    plt.close(fig)

    plate_counts = order.groupby(["lead", "purpose"]).size().reset_index(name="n")
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    y = np.arange(len(plate_counts))
    ax.barh(y, plate_counts["n"], color="#76b7b2")
    ax.set_yticks(y, plate_counts["lead"] + " | " + plate_counts["purpose"])
    ax.set_xlabel("reagent/control count")
    ax.set_title("Control coverage for two-lead wetlab package")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig_md39_control_coverage.png", dpi=220)
    fig.savefig(FIG / "fig_md39_control_coverage.pdf")
    plt.close(fig)

    gates = thresholds["gate"].drop_duplicates().tolist()
    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.axis("off")
    xs = np.linspace(0.10, 0.90, len(gates))
    for i, gate in enumerate(gates):
        ax.add_patch(plt.Circle((xs[i], 0.56), 0.07, color="#17324f", ec="#edf5ff", lw=1.2))
        ax.text(xs[i], 0.56, str(i + 1), ha="center", va="center", color="#f2c46d", fontsize=16, weight="bold")
        ax.text(xs[i], 0.33, gate, ha="center", va="center", fontsize=10)
        if i < len(gates) - 1:
            ax.annotate("", xy=(xs[i + 1] - 0.08, 0.56), xytext=(xs[i] + 0.08, 0.56), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.83, "Assay go/no-go ladder before immunogenicity claims", ha="center", fontsize=15, weight="bold")
    fig.savefig(FIG / "fig_md40_assay_go_no_go_ladder.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md40_assay_go_no_go_ladder.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.axis("off")
    for i, row in moat.iterrows():
        y = 0.88 - i * 0.17
        ax.add_patch(plt.Rectangle((0.04, y - 0.055), 0.92, 0.11, fc="#10243a", ec="#395170"))
        ax.text(0.06, y + 0.022, row["reviewer_attack"], color="#f2c46d", fontsize=9, weight="bold")
        ax.text(0.06, y - 0.025, row["defense_asset"], color="white", fontsize=8.5)
    ax.set_title("Reviewer-defense moat for assay-ready translation", fontsize=15, weight="bold")
    fig.savefig(FIG / "fig_md41_reviewer_defense_moat.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md41_reviewer_defense_moat.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(reagents: pd.DataFrame, plate: pd.DataFrame, thresholds: pd.DataFrame, moat: pd.DataFrame) -> None:
    lines = [
        "# Assay-Ready Translation Package",
        "",
        "## Decision",
        "",
        "The two lead candidates have been converted into an assay-ready package with mutant, WT, decoy, pHLA reagent, plate-map, and go/no-go criteria. This is the practical bridge from computational prioritization to controlled immunogenicity testing.",
        "",
        "## Reagent Order Sheet",
        "",
        reagents.to_markdown(index=False),
        "",
        "## Plate Map",
        "",
        plate.head(96).to_markdown(index=False),
        "",
        "## Go/No-Go Thresholds",
        "",
        thresholds.to_markdown(index=False),
        "",
        "## Reviewer Moat",
        "",
        moat.to_markdown(index=False),
        "",
        "## Claim Boundary",
        "",
        "This package makes the work assay-ready. It still does not prove immunogenicity until mutant responses beat WT and decoy controls in activation assays.",
    ]
    (OUT / "ASSAY_READY_TRANSLATION_REPORT.md").write_text("\n".join(lines) + "\n")

    kr = [
        "# Assay-ready translation Korean brief",
        "",
        "## 핵심",
        "",
        "이제 두 lead는 단순 후보가 아니라 바로 실험팀에 넘길 수 있는 형태입니다.",
        "",
        "- mutant peptide",
        "- WT peptide",
        "- anchor-preserved decoy",
        "- pHLA monomer/tetramer",
        "- plate map",
        "- go/no-go 기준",
        "- reviewer 방어 논리",
        "",
        "## 가장 중요한 claim boundary",
        "",
        "WT/decoy보다 mutant가 HLA binding, multimer/TCR binding, IFN-gamma/ICS에서 강해야만 immunogenicity claim으로 올립니다.",
    ]
    (OUT / "ASSAY_READY_ONE_PAGE_KR.md").write_text("\n".join(kr) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    leads = load_leads()
    reagents = build_reagent_order_sheet(leads)
    plate = build_plate_map(leads)
    thresholds = build_go_no_go_thresholds(leads)
    moat = build_reviewer_moat(leads)
    reagents.to_csv(OUT / "assay_reagent_order_sheet.tsv", sep="\t", index=False)
    plate.to_csv(OUT / "assay_plate_map.tsv", sep="\t", index=False)
    thresholds.to_csv(OUT / "assay_go_no_go_thresholds.tsv", sep="\t", index=False)
    moat.to_csv(OUT / "assay_reviewer_defense_moat.tsv", sep="\t", index=False)
    make_figures(reagents, thresholds, moat)
    write_report(reagents, plate, thresholds, moat)
    summary = {
        "n_reagent_rows": int(len(reagents)),
        "n_plate_rows": int(len(plate)),
        "n_go_no_go_rules": int(len(thresholds)),
        "figures": [
            "fig_md38_reagent_order_sheet",
            "fig_md39_control_coverage",
            "fig_md40_assay_go_no_go_ladder",
            "fig_md41_reviewer_defense_moat",
        ],
        "boundary": "assay readiness, not immunogenicity proof",
    }
    (OUT / "assay_ready_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
