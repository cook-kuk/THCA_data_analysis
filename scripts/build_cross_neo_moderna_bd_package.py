#!/usr/bin/env python3
"""Build an IP-safe Moderna/Merck BD package for CROSS-Neo KG-GA.

This creates non-confidential and NDA-stage materials for a prospective
validation/license discussion. It does not generate protected manuscript prose.
"""

from __future__ import annotations

import html
import json
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
GA = ROOT / "project/results/cross_neo_ga_rl_full_comparison_2026_05_10"
MASTER = GA / "ga_rl_master_algorithm_comparison.tsv"
PUBLIC = GA / "ga_rl_public_competitor_context_summary.tsv"
CLAIM = GA / "ga_rl_claim_boundary_matrix.tsv"
OUT = ROOT / "project/results/cross_neo_moderna_bd_package_2026_05_11"
FIG = OUT / "figures"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = HUB / "assets/cross_neo_moderna_bd_package"
LIVE_ASSET = LIVE / "assets/cross_neo_moderna_bd_package"
PAGE = HUB / "cross_neo_moderna_bd_package.html"
LIVE_PAGE = LIVE / "cross_neo_moderna_bd_package.html"

OFFICIAL_SOURCES = [
    {
        "source": "Merck/Moderna melanoma Phase 3 announcement",
        "url": "https://www.merck.com/news/merck-and-moderna-initiate-phase-3-study-evaluating-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-for-adjuvant-treatment-of-patients-with-resected-high-riskstage-iib-iv-melanom/",
        "relevance": "V940/mRNA-4157 Phase 3 adjuvant melanoma program; approximately 1,089 patients; RFS primary endpoint.",
    },
    {
        "source": "Merck/Moderna NSCLC Phase 3 announcement",
        "url": "https://www.merck.com/news/merck-and-moderna-initiate-phase-3-trial-evaluating-adjuvant-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-after-neoadjuvant-keytruda-and-chemotherapy-in-patients-with-certain-ty/",
        "relevance": "INTerpath-009 expansion; V940 described as mRNA coding for up to 34 algorithmically derived neoantigens.",
    },
]

CONTACT_ROUTES = [
    {
        "company": "Merck",
        "route": "Business development and licensing (BD&L)",
        "url": "https://www.merck.com/research/business-development-and-licensing/",
        "why": "Official Merck collaboration route; includes oncology, vaccines, and enabling/platform technologies.",
        "pre_nda_action": "Use non-confidential teaser only; request CDA/NDA and blinded validation discussion.",
        "confidentiality_note": "Do not send source code, feature list, candidate tables, or scoring weights before CDA/NDA.",
    },
    {
        "company": "Merck",
        "route": "Corporate contact fallback",
        "url": "https://www.merck.com/contact-us/",
        "why": "Fallback for routing to BD&L or oncology external innovation if web BD intake is not available.",
        "pre_nda_action": "Ask for the correct BD&L contact for individualized neoantigen vaccine selection technology.",
        "confidentiality_note": "General contact route; keep all technical details non-confidential.",
    },
    {
        "company": "Moderna",
        "route": "Official contact page",
        "url": "https://www.modernatx.com/en-US/contact-moderna",
        "why": "Official Moderna page lists corporate contact routes, HQ, media, investor, and general contact channels.",
        "pre_nda_action": "Use only to request routing or a warm BD contact; do not disclose internals through general inboxes.",
        "confidentiality_note": "No public BD-specific email was verified here; treat as routing only, not a confidential submission portal.",
    },
    {
        "company": "Moderna/Merck",
        "route": "Warm intro / conference BD route",
        "url": "N/A",
        "why": "Best practical route for a sensitive platform/selection-layer pitch when no public Moderna BD inbox is verified.",
        "pre_nda_action": "Send one-paragraph teaser and ask for CDA/NDA-enabled technical diligence.",
        "confidentiality_note": "Do not attach NDA-stage validation protocol unless CDA/NDA is in place.",
    },
]


AUTHOR_NAME = "Seungho Cook"
CONTACT_EMAIL = "kukshomr@snu.ac.kr"
CONTACT_HANDLE = "@Ho15421542"


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def esc(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return html.escape(str(value))


def fmt(value: object, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return esc(value)


def safe_copy(src: Path, dst: Path) -> None:
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    except PermissionError:
        print(f"[moderna-bd] skip permission-denied copy {dst}")


def table_html(df: pd.DataFrame, cols: list[str], max_rows: int = 50) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    keep = [c for c in cols if c in df.columns]
    rows = ["<table><thead><tr>"]
    rows += [f"<th>{esc(c.replace('_', ' '))}</th>" for c in keep]
    rows.append("</tr></thead><tbody>")
    for _, r in df.head(max_rows).iterrows():
        rows.append("<tr>")
        for c in keep:
            v = r[c]
            rows.append(f"<td>{fmt(v) if isinstance(v, (int, float, np.integer, np.floating)) else esc(v)}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def markdown_table(df: pd.DataFrame) -> str:
    return df.where(pd.notna(df), "NA").to_markdown(index=False)


def build_metrics_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    master = read_tsv(MASTER)
    public = read_tsv(PUBLIC)
    claim = read_tsv(CLAIM)

    same_board_models = [
        "KG_GA_evolved",
        "CROSS_integrated",
        "CROSS_claimsafe",
        "CROSS_stress",
        "CROSS_finetuned",
        "BigMHC_IM",
        "BigMHC_EL",
    ]
    same = master[master["algorithm"].isin(same_board_models)].copy()
    same = same[
        [
            "algorithm",
            "algorithm_family",
            "all_AUPRC",
            "all_AUROC",
            "frozen_validation_AUPRC",
            "low_medium_leakage_AUPRC",
            "v6_smoke_AUPRC",
            "primary_disposition",
        ]
    ].sort_values("all_AUPRC", ascending=False)
    for c in same.columns:
        if c.endswith("AUPRC") or c.endswith("AUROC"):
            same[c] = pd.to_numeric(same[c], errors="coerce")

    strict_models = [
        "BAR-Neo_confidence",
        "BAR-Neo",
        "RF_biophys",
        "MHCflurry_2.0_presentation",
        "KG_GA_evolved_controller",
    ]
    strict = master[master["algorithm"].isin(strict_models)].copy()
    strict = strict[
        [
            "algorithm",
            "algorithm_family",
            "claim_track",
            "neo_strict_AUPRC",
            "neo_strict_AUROC",
            "patients_evaluated",
            "primary_disposition",
        ]
    ].sort_values("neo_strict_AUPRC", ascending=False)
    for c in ["neo_strict_AUPRC", "neo_strict_AUROC", "patients_evaluated"]:
        strict[c] = pd.to_numeric(strict[c], errors="coerce")

    return same, strict, public, claim


def build_bd_tables(same: pd.DataFrame, strict: pd.DataFrame, public: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    readiness = pd.DataFrame(
        [
            {
                "stage": "0. Internal lock",
                "status": "done",
                "deliverable": "Frozen local KG-GA comparison package",
                "exit_criterion": "AUPRC 0.933 on 2,715-row board with explicit caveats",
                "share_level": "internal only",
            },
            {
                "stage": "1. IP hygiene",
                "status": "next",
                "deliverable": "Provisional patent / invention disclosure",
                "exit_criterion": "Claims cover KG-guided evolutionary neoantigen selection, audit gates, and assay-control generation",
                "share_level": "attorney / tech transfer",
            },
            {
                "stage": "2. Non-confidential teaser",
                "status": "ready",
                "deliverable": "One-page BD teaser with masked method details",
                "exit_criterion": "No code, no full feature list, no exact candidate tables",
                "share_level": "pre-NDA",
            },
            {
                "stage": "3. CDA/NDA",
                "status": "next",
                "deliverable": "Mutual NDA and blinded benchmark data-use terms",
                "exit_criterion": "Sponsor cannot use submitted scores to train without written license",
                "share_level": "legal",
            },
            {
                "stage": "4. Blinded validation",
                "status": "design-ready",
                "deliverable": "Frozen scoring container/API run on sponsor blinded candidates",
                "exit_criterion": "Pre-specified top-34 enrichment/AUPRC lift vs incumbent selection stack",
                "share_level": "NDA",
            },
            {
                "stage": "5. Option/license",
                "status": "conditional",
                "deliverable": "Field-limited oncology neoantigen selection license",
                "exit_criterion": "Payment tied to validation pass, clinical program use, and milestones",
                "share_level": "BD/legal",
            },
        ]
    )

    validation = pd.DataFrame(
        [
            {
                "module": "Blinded candidate input",
                "required_from_partner": "patient_id, mutant peptide, optional WT peptide, HLA, expression, variant context, current internal scores if available",
                "our_commitment": "Freeze parser and score generation before label unblinding",
                "pass_signal": ">=95% candidate parse success and reproducible score file hash",
            },
            {
                "module": "Primary endpoint",
                "required_from_partner": "binary or ordinal immunogenicity / selection / wetlab labels withheld until lock",
                "our_commitment": "Top-34 ranking and calibrated score per patient",
                "pass_signal": "AUPRC or top-34 hit-rate lift versus incumbent sponsor baseline",
            },
            {
                "module": "Slot-value endpoint",
                "required_from_partner": "per-patient vaccine slot budget, ideally 20-34",
                "our_commitment": "Rank candidates under a V940-like slot constraint",
                "pass_signal": "More validated candidates in fixed 34-slot payload without loss of HLA/source diversity",
            },
            {
                "module": "False-positive control",
                "required_from_partner": "decoys, WT peptides, benign/self-like candidates, failed assay records",
                "our_commitment": "Report claim-safe and production-priority scores separately",
                "pass_signal": "Lower false-positive burden at matched sensitivity",
            },
            {
                "module": "Explainability",
                "required_from_partner": "no extra confidential biology required",
                "our_commitment": "Per-candidate reason codes, uncertainty, and audit flags",
                "pass_signal": "Sponsor scientist can adjudicate top candidates without opening model internals",
            },
            {
                "module": "Wetlab translation",
                "required_from_partner": "assay format, HLA cell line/APC constraints, readout thresholds",
                "our_commitment": "Mutant/WT/decoy control map and top-N plate-ready export",
                "pass_signal": "Assay-ready list accepted by wetlab operations",
            },
        ]
    )

    claim_boundary = pd.DataFrame(
        [
            {
                "claim": "Allowed pre-NDA",
                "wording": "Retrospective CROSS-Neo prioritization engine with AUPRC 0.933 on a 2,715-candidate local board",
                "why": "Supported by same-board metrics",
                "do_not_say": "Clinically validated or Moderna-ready efficacy proof",
            },
            {
                "claim": "Allowed under NDA",
                "wording": "Frozen algorithm can be evaluated on blinded sponsor candidate pools for top-34 enrichment",
                "why": "Validation protocol is pre-specified",
                "do_not_say": "Will outperform sponsor stack before seeing blinded result",
            },
            {
                "claim": "Comparator boundary",
                "wording": "Same-board comparison beats BigMHC_IM/EL and fixed CROSS-Neo; public context includes RF_biophys/MHCflurry/PRIME/TransPHLA",
                "why": "Different benchmark tables have different leakage/overlap caveats",
                "do_not_say": "Universal SOTA over all public algorithms",
            },
            {
                "claim": "Strict reliability boundary",
                "wording": "BAR-Neo_confidence is still the strict no-overlap reliability winner in current local stack",
                "why": "Strict AUPRC 0.857 vs KG_GA_controller 0.624",
                "do_not_say": "KG-GA wins every split",
            },
            {
                "claim": "Wetlab boundary",
                "wording": "Assay-ready prioritization and control design",
                "why": "No prospective mutant-vs-WT/decoy readout yet",
                "do_not_say": "Immunogenicity proof",
            },
        ]
    )

    data_room = pd.DataFrame(
        [
            {
                "folder": "00_non_confidential",
                "item": "MODERNA_NONCONFIDENTIAL_TEASER.md",
                "share_timing": "pre-NDA",
                "contains": "Value proposition, masked metrics, validation ask, no algorithm internals",
            },
            {
                "folder": "01_validation_protocol",
                "item": "MODERNA_BLINDED_VALIDATION_PROTOCOL.md",
                "share_timing": "NDA",
                "contains": "Input schema, frozen scoring, endpoints, success criteria",
            },
            {
                "folder": "02_metrics",
                "item": "moderna_same_board_comparator_metrics.tsv",
                "share_timing": "NDA",
                "contains": "KG-GA vs CROSS/BigMHC same-board metrics",
            },
            {
                "folder": "02_metrics",
                "item": "moderna_strict_reliability_metrics.tsv",
                "share_timing": "NDA",
                "contains": "BAR-Neo, RF_biophys, MHCflurry, KG-GA-controller strict no-overlap context",
            },
            {
                "folder": "03_claim_boundary",
                "item": "moderna_claim_boundary.tsv",
                "share_timing": "pre-NDA summary / NDA detail",
                "contains": "Allowed claims and blocked claims",
            },
            {
                "folder": "04_deal",
                "item": "MODERNA_OPTION_LICENSE_TERM_SHEET_SCAFFOLD.md",
                "share_timing": "after technical interest",
                "contains": "Option/license terms scaffold; not legal advice",
            },
            {
                "folder": "05_ip",
                "item": "INVENTION_DISCLOSURE_SCAFFOLD.md",
                "share_timing": "attorney / tech transfer before outreach",
                "contains": "Provisional-patent disclosure skeleton and claim areas",
            },
            {
                "folder": "06_outreach",
                "item": "MODERNA_OUTREACH_EMAIL_AFTER_IP.md",
                "share_timing": "after IP hygiene and contact selection",
                "contains": "Short pre-NDA email draft that does not disclose internals",
            },
        ]
    )

    return readiness, validation, claim_boundary, data_room


def build_contact_routes() -> pd.DataFrame:
    return pd.DataFrame(CONTACT_ROUTES)


def write_markdown_docs(
    same: pd.DataFrame,
    strict: pd.DataFrame,
    public: pd.DataFrame,
    readiness: pd.DataFrame,
    validation: pd.DataFrame,
    claim_boundary: pd.DataFrame,
    data_room: pd.DataFrame,
    contact_routes: pd.DataFrame,
) -> None:
    kg = same[same["algorithm"].eq("KG_GA_evolved")].iloc[0]
    integrated = same[same["algorithm"].eq("CROSS_integrated")].iloc[0]
    big_im = same[same["algorithm"].eq("BigMHC_IM")].iloc[0]
    big_el = same[same["algorithm"].eq("BigMHC_EL")].iloc[0]
    bar = strict[strict["algorithm"].eq("BAR-Neo_confidence")].iloc[0]
    kg_controller = strict[strict["algorithm"].eq("KG_GA_evolved_controller")].iloc[0]

    teaser = [
        "# Non-Confidential Teaser: CROSS-Neo KG-GA Selection Layer",
        "",
        "## One-Line Positioning",
        "",
        "CROSS-Neo KG-GA is an algorithmic neoantigen-selection layer designed to improve which candidates enter a fixed-slot individualized mRNA vaccine payload, before manufacturing.",
        "",
        "## Why This Fits A V940-Like Workflow",
        "",
        "- Public Moderna/Merck materials describe V940/mRNA-4157 as an individualized neoantigen therapy encoding up to 34 algorithmically derived neoantigens.",
        "- The value gap is therefore not mRNA manufacturing; it is better candidate selection under a constrained payload budget.",
        "- CROSS-Neo outputs ranked candidates, uncertainty, reason codes, and wetlab-control-ready exports.",
        "",
        "## Masked Evidence Snapshot",
        "",
        f"- Same-board 2,715-candidate result: KG-GA AUPRC {kg['all_AUPRC']:.3f} / AUROC {kg['all_AUROC']:.3f}.",
        f"- Fixed integrated comparator: AUPRC {integrated['all_AUPRC']:.3f}.",
        f"- Public BigMHC context on the same board: IM AUPRC {big_im['all_AUPRC']:.3f}; EL AUPRC {big_el['all_AUPRC']:.3f}.",
        f"- Frozen validation-like subset: KG-GA AUPRC {kg['frozen_validation_AUPRC']:.3f}.",
        f"- Strict no-overlap reliability remains a separate lane: BAR-Neo_confidence AUPRC {bar['neo_strict_AUPRC']:.3f}; KG-GA-controller AUPRC {kg_controller['neo_strict_AUPRC']:.3f}.",
        "",
        "## Validation Ask",
        "",
        "Under CDA/NDA, provide a blinded candidate pool with sponsor labels withheld. We return frozen scores, top-34 ranking, uncertainty, and assay-control recommendations. Labels are unblinded only after score lock.",
        "",
        "## Pre-NDA Boundary",
        "",
        "No code, full feature list, candidate tables, training corpus details, or proprietary scoring weights should be shared before IP filing and CDA/NDA.",
        "",
        "## Official Context Sources",
        "",
    ]
    teaser += [f"- [{s['source']}]({s['url']}) — {s['relevance']}" for s in OFFICIAL_SOURCES]
    teaser += [
        "",
        "## Contact",
        "",
        f"- {AUTHOR_NAME}",
        f"- Email: {CONTACT_EMAIL}",
        f"- X/Twitter: {CONTACT_HANDLE}",
        "",
    ]
    (OUT / "MODERNA_NONCONFIDENTIAL_TEASER.md").write_text("\n".join(teaser) + "\n", encoding="utf-8")

    protocol = [
        "# Moderna/Merck Blinded Validation Protocol Scaffold",
        "",
        "## Objective",
        "",
        "Test whether CROSS-Neo KG-GA improves fixed-slot individualized neoantigen selection versus sponsor incumbent scoring on blinded candidate pools.",
        "",
        "## Data Flow",
        "",
        "1. Sponsor provides candidate table without outcome labels or incumbent rank if possible.",
        "2. CROSS-Neo team freezes parser, model version, config hash, and score columns before processing.",
        "3. CROSS-Neo returns per-candidate score, per-patient top-34 list, reason codes, uncertainty flags, and assay-control annotations.",
        "4. Sponsor unblinds labels and calculates pre-specified endpoints.",
        "5. Joint review decides whether to expand to wetlab or license option.",
        "",
        "## Minimum Input Schema",
        "",
        "| column | required | note |",
        "|---|---|---|",
        "| patient_id | yes | de-identified |",
        "| mutant_peptide | yes | 8-15 aa for MHC-I path; longer class-II can be separated |",
        "| hla_allele | yes | high-resolution preferred |",
        "| wt_peptide | preferred | enables mutant-vs-WT specificity filters |",
        "| expression / RNA support | preferred | optional but valuable |",
        "| variant context | preferred | SNV/indel/fusion/source metadata if available |",
        "| sponsor_score | optional | held out until after lock if possible |",
        "| wetlab_label | yes, withheld | binary/ordinal response, presentation, or selection outcome |",
        "",
        "## Success Criteria",
        "",
        markdown_table(validation),
        "",
        "## Statistical Readout",
        "",
        "- Primary: per-candidate AUPRC and per-patient top-34 hit rate versus incumbent baseline.",
        "- Secondary: AUROC, top-10/top-20/top-34 precision, calibration, false-positive burden at matched sensitivity.",
        "- Stratification: tumor type, HLA family, mutation class, expression availability, payload slot budget.",
        "- Hard rule: no model fitting, threshold tuning, or weight update on sponsor labels before primary analysis.",
        "",
        "## Decision Gates",
        "",
        markdown_table(readiness),
        "",
    ]
    (OUT / "MODERNA_BLINDED_VALIDATION_PROTOCOL.md").write_text("\n".join(protocol) + "\n", encoding="utf-8")

    terms = [
        "# Option/License Term Sheet Scaffold",
        "",
        "Not legal advice. This is a negotiation scaffold for counsel/BD review.",
        "",
        "## Proposed Structure",
        "",
        "- Stage 1: mutual CDA/NDA and blinded validation.",
        "- Stage 2: paid evaluation or sponsored validation project.",
        "- Stage 3: time-limited exclusive option in individualized oncology neoantigen selection.",
        "- Stage 4: field-limited license if validation passes.",
        "",
        "## Assets",
        "",
        "- KG-guided evolutionary neoantigen prioritization algorithm.",
        "- Claim-safe fallback scoring lane.",
        "- Audit flags, uncertainty, and reason-code layer.",
        "- Assay-control and top-N plate export layer.",
        "",
        "## Commercial Hooks",
        "",
        "- Upfront evaluation fee after NDA.",
        "- Option fee after blinded benchmark pass.",
        "- License fee for clinical-program use.",
        "- Milestones tied to clinical trial use, regulatory progress, and commercial launch.",
        "- Running royalty or per-patient scoring fee.",
        "",
        "## Guardrails",
        "",
        "- No sponsor training/use of submitted scores outside the agreed evaluation without license.",
        "- No disclosure of model internals before IP filing.",
        "- No exclusivity outside individualized oncology neoantigen selection unless separately priced.",
        "- Publications or public claims require mutual review.",
        "",
    ]
    (OUT / "MODERNA_OPTION_LICENSE_TERM_SHEET_SCAFFOLD.md").write_text("\n".join(terms) + "\n", encoding="utf-8")

    invention = [
        "# Invention Disclosure Scaffold",
        "",
        "Not legal advice. Use this as a structured intake for counsel or technology transfer before any external outreach.",
        "",
        "## Working Title",
        "",
        "Knowledge-graph-guided evolutionary selection of individualized cancer vaccine neoantigens with audit-controlled wetlab translation.",
        "",
        "## Problem",
        "",
        "Individualized mRNA vaccine payloads have limited antigen slots. Candidate pools can contain many plausible peptides, but existing presentation or immunogenicity scores do not fully solve false-positive burden, explainability, patient-level slot allocation, or wetlab-control design.",
        "",
        "## Technical Solution",
        "",
        "- Use a knowledge graph to encode biological, assay, HLA, presentation, immunogenicity, self-similarity, TCR, and translation constraints.",
        "- Use a genetic/evolutionary controller to search candidate scoring/selection policies.",
        "- Maintain separated production-priority and claim-safe scoring lanes.",
        "- Emit per-candidate reason codes, uncertainty/audit flags, and top-N wetlab-control maps.",
        "",
        "## Candidate Claim Areas",
        "",
        "- KG-guided evolutionary optimization for neoantigen ranking under a fixed payload slot budget.",
        "- Dual-lane prioritization: production/experiment-priority score plus claim-safe fallback score.",
        "- Blinded validation workflow with frozen score hashes before label unblinding.",
        "- Automated mutant/WT/decoy assay-control selection tied to candidate ranking.",
        "- Patient-level top-N payload selection with HLA/source diversity constraints and uncertainty flags.",
        "",
        "## Evidence Snapshot",
        "",
        f"- Same-board 2,715-candidate KG-GA: AUPRC {kg['all_AUPRC']:.3f}, AUROC {kg['all_AUROC']:.3f}.",
        f"- Frozen validation-like subset: AUPRC {kg['frozen_validation_AUPRC']:.3f}.",
        f"- Current strict no-overlap reliability caveat: BAR-Neo_confidence AUPRC {bar['neo_strict_AUPRC']:.3f} vs KG-GA-controller {kg_controller['neo_strict_AUPRC']:.3f}.",
        "",
        "## Do Not Disclose Before Filing",
        "",
        "- Source code, exact feature list, scoring weights, graph schema details, candidate-level score tables, training corpus joins, or endpoint labels.",
        "",
    ]
    (OUT / "INVENTION_DISCLOSURE_SCAFFOLD.md").write_text("\n".join(invention) + "\n", encoding="utf-8")

    outreach = [
        "# Moderna/Merck Outreach Email Draft After IP Hygiene",
        "",
        "Use only after invention disclosure/provisional filing and contact selection. Do not attach code, candidate tables, or full method details.",
        "",
        "## Subject Options",
        "",
        "- Blinded validation proposal: neoantigen selection layer for individualized mRNA vaccine payloads",
        "- NDA-stage benchmark proposal for top-34 neoantigen selection",
        "- IP-safe proposal: improving candidate selection before individualized mRNA vaccine manufacturing",
        "",
        "## Short Email",
        "",
        "Dear [Name],",
        "",
        "I am reaching out with an IP-protected neoantigen-selection technology that may be relevant to individualized mRNA vaccine workflows such as fixed-slot personalized neoantigen payload design.",
        "",
        "In a local 2,715-candidate benchmark, our knowledge-graph-guided evolutionary selector achieved AUPRC 0.933 / AUROC 0.943 and outperformed fixed internal and public presentation/immunogenicity comparators on that board. We are not claiming prospective clinical validation; the proposed next step is a blinded NDA-stage evaluation on sponsor-provided candidate pools, with labels withheld until score lock.",
        "",
        "The output is a per-patient top-N/top-34 ranking with uncertainty, reason codes, and assay-control recommendations. We can share a non-confidential one-page summary and, under CDA/NDA, a pre-specified validation protocol.",
        "",
        "Would your team be open to a brief discussion about a blinded benchmark collaboration or option-to-license evaluation?",
        "",
        "Best,",
        AUTHOR_NAME,
        CONTACT_EMAIL,
        f"X/Twitter: {CONTACT_HANDLE}",
        "",
        "## Attachment Rule",
        "",
        "Attach only `MODERNA_NONCONFIDENTIAL_TEASER.md` converted to PDF. Do not attach the full dossier, TSVs, figures with internals, source code, or candidate-level score exports pre-NDA.",
        "",
    ]
    (OUT / "MODERNA_OUTREACH_EMAIL_AFTER_IP.md").write_text("\n".join(outreach) + "\n", encoding="utf-8")

    overview = [
        "# CROSS-Neo Moderna-Ready BD Package",
        "",
        "## Executive Decision",
        "",
        "The package is now framed as a selection-layer collaboration, not a claim that CROSS-Neo replaces Moderna's mRNA platform.",
        "",
        "## Core Evidence",
        "",
        markdown_table(same),
        "",
        "## Strict Reliability Context",
        "",
        markdown_table(strict),
        "",
        "## Broader Public Competitor Context",
        "",
        markdown_table(public),
        "",
        "## Claim Boundary",
        "",
        markdown_table(claim_boundary),
        "",
        "## Data Room Index",
        "",
        markdown_table(data_room),
        "",
        "## Contact Routing",
        "",
        markdown_table(contact_routes),
        "",
    ]
    (OUT / "MODERNA_BD_PACKAGE_OVERVIEW.md").write_text("\n".join(overview) + "\n", encoding="utf-8")

    send_readme = [
        "# Pre-NDA Send Packet README",
        "",
        "This folder is the only packet intended for pre-NDA outreach.",
        "",
        "## Include",
        "",
        "- `MODERNA_NONCONFIDENTIAL_TEASER.md`",
        "- `MODERNA_OUTREACH_EMAIL_AFTER_IP.md`",
        "- `moderna_contact_routes.tsv`",
        "",
        "## Do Not Include Pre-NDA",
        "",
        "- Source code",
        "- Candidate-level score tables",
        "- Full feature lists or scoring weights",
        "- Training corpus joins",
        "- NDA validation protocol unless CDA/NDA is already signed",
        "- Any claims of prospective clinical efficacy or immunogenicity proof",
        "",
        "## Required Sequence",
        "",
        "1. File invention disclosure/provisional or receive counsel clearance.",
        "2. Identify the contact route or warm intro.",
        "3. Send the short outreach email with only the non-confidential teaser.",
        "4. Ask for CDA/NDA before sharing the blinded validation protocol.",
        "5. Run frozen blinded validation before discussing option/license economics.",
        "",
    ]
    (OUT / "PRE_NDA_SEND_PACKET_README.md").write_text("\n".join(send_readme) + "\n", encoding="utf-8")


def build_pre_nda_packet(contact_routes: pd.DataFrame) -> None:
    packet = OUT / "pre_nda_send_packet"
    if packet.exists():
        shutil.rmtree(packet)
    packet.mkdir(parents=True, exist_ok=True)
    for name in [
        "MODERNA_NONCONFIDENTIAL_TEASER.md",
        "MODERNA_OUTREACH_EMAIL_AFTER_IP.md",
        "PRE_NDA_SEND_PACKET_README.md",
    ]:
        shutil.copy2(OUT / name, packet / name)
    contact_routes.to_csv(packet / "moderna_contact_routes.tsv", sep="\t", index=False)
    archive_base = OUT / "pre_nda_send_packet"
    zip_path = OUT / "pre_nda_send_packet.zip"
    if zip_path.exists():
        zip_path.unlink()
    shutil.make_archive(str(archive_base), "zip", packet)


def build_figures(same: pd.DataFrame, strict: pd.DataFrame, readiness: pd.DataFrame, validation: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    plot = same.dropna(subset=["all_AUPRC"]).sort_values("all_AUPRC")
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    colors = ["#c7922b" if a == "KG_GA_evolved" else "#315f7d" for a in plot["algorithm"]]
    ax.barh(plot["algorithm"], plot["all_AUPRC"], color=colors)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Same-board AUPRC")
    ax.set_title("Selection-layer evidence: KG-GA vs current local comparators")
    for i, (_, r) in enumerate(plot.iterrows()):
        ax.text(float(r["all_AUPRC"]) + 0.01, i, f"{float(r['all_AUPRC']):.3f}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_same_board_auprc_positioning.png", dpi=220)
    fig.savefig(FIG / "fig1_same_board_auprc_positioning.pdf")
    plt.close(fig)

    strict_plot = strict.dropna(subset=["neo_strict_AUPRC"]).sort_values("neo_strict_AUPRC")
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    colors = ["#b24a62" if "KG_GA" in a else "#315f7d" for a in strict_plot["algorithm"]]
    ax.barh(strict_plot["algorithm"], strict_plot["neo_strict_AUPRC"], color=colors)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Strict no-overlap AUPRC")
    ax.set_title("Reliability context: keep BAR-Neo lane separate")
    for i, (_, r) in enumerate(strict_plot.iterrows()):
        ax.text(float(r["neo_strict_AUPRC"]) + 0.01, i, f"{float(r['neo_strict_AUPRC']):.3f}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_strict_reliability_boundary.png", dpi=220)
    fig.savefig(FIG / "fig2_strict_reliability_boundary.pdf")
    plt.close(fig)

    stages = readiness["stage"].tolist()
    status_score = readiness["status"].map({"done": 1.0, "ready": 0.8, "design-ready": 0.65, "next": 0.35, "conditional": 0.2}).fillna(0.2)
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(np.arange(len(stages)), status_score, marker="o", lw=2, color="#c7922b")
    ax.set_ylim(0, 1.05)
    ax.set_xticks(np.arange(len(stages)), stages, rotation=25, ha="right")
    ax.set_ylabel("Readiness")
    ax.set_title("BD readiness ladder: from local evidence to option/license")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_bd_readiness_ladder.png", dpi=220)
    fig.savefig(FIG / "fig3_bd_readiness_ladder.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.axis("off")
    boxes = [
        ("Sponsor blinded\ncandidate pool", 0.05, 0.55),
        ("Frozen CROSS-Neo\nKG-GA scoring", 0.29, 0.55),
        ("Top-34 payload\nranking + flags", 0.53, 0.55),
        ("Label unblind\nand endpoint test", 0.77, 0.55),
        ("Option/license\nif pass", 0.77, 0.18),
    ]
    for text, x, y in boxes:
        ax.add_patch(plt.Rectangle((x, y), 0.18, 0.22, facecolor="#102236", edgecolor="#c7922b", lw=1.5))
        ax.text(x + 0.09, y + 0.11, text, ha="center", va="center", color="white", fontsize=10, weight="bold")
    for x1, y1, x2, y2 in [(0.23, 0.66, 0.29, 0.66), (0.47, 0.66, 0.53, 0.66), (0.71, 0.66, 0.77, 0.66), (0.86, 0.55, 0.86, 0.40)]:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1.8, color="#c7922b"))
    ax.text(0.05, 0.08, "Hard rule: no model fitting on sponsor labels before primary analysis.", color="#cbd5e1", fontsize=11)
    ax.set_title("Blinded validation workflow", color="black", pad=12)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_blinded_validation_workflow.png", dpi=220)
    fig.savefig(FIG / "fig4_blinded_validation_workflow.pdf")
    plt.close(fig)


def build_html(
    same: pd.DataFrame,
    strict: pd.DataFrame,
    public: pd.DataFrame,
    readiness: pd.DataFrame,
    validation: pd.DataFrame,
    claim_boundary: pd.DataFrame,
    data_room: pd.DataFrame,
    contact_routes: pd.DataFrame,
) -> str:
    kg = same[same["algorithm"].eq("KG_GA_evolved")].iloc[0]
    cards = ""
    for fig, title, note in [
        ("fig1_same_board_auprc_positioning.png", "Same-Board Positioning", "KG-GA against CROSS and BigMHC comparators."),
        ("fig2_strict_reliability_boundary.png", "Strict Reliability Boundary", "BAR-Neo remains the strict no-overlap reliability lane."),
        ("fig3_bd_readiness_ladder.png", "BD Readiness Ladder", "What must happen before partner outreach and licensing."),
        ("fig4_blinded_validation_workflow.png", "Blinded Validation Workflow", "NDA-stage validation without leaking labels or internals."),
    ]:
        cards += f"""
        <article class="card">
          <img src="assets/cross_neo_moderna_bd_package/{fig}" alt="{esc(title)}">
          <h3>{esc(title)}</h3>
          <p>{esc(note)}</p>
        </article>
        """

    sources = "\n".join(
        [f"<li><a href='{esc(s['url'])}'>{esc(s['source'])}</a> — {esc(s['relevance'])}</li>" for s in OFFICIAL_SOURCES]
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo Moderna-Ready BD Package</title>
<style>
:root {{ --bg:#07111f; --panel:#101d2d; --ink:#edf5ff; --muted:#9fb0c5; --gold:#f2c46d; --line:#27384f; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter,Arial,sans-serif; line-height:1.5; }}
header {{ padding:48px 5vw 28px; border-bottom:1px solid var(--line); background:#0b1626; }}
.kicker {{ color:var(--gold); text-transform:uppercase; letter-spacing:.08em; font-size:12px; }}
h1 {{ font-size:clamp(30px,5vw,60px); margin:8px 0 10px; line-height:1.0; }}
.lead {{ color:var(--muted); max-width:1120px; font-size:18px; }}
.stats {{ display:grid; grid-template-columns:repeat(5,minmax(140px,1fr)); gap:12px; margin-top:24px; }}
.stat {{ background:var(--panel); border:1px solid var(--line); padding:14px; border-radius:8px; }}
.stat b {{ display:block; font-size:26px; color:var(--gold); }}
.stat span {{ color:var(--muted); font-size:12px; }}
main {{ padding:28px 5vw 60px; max-width:1500px; margin:0 auto; }}
section {{ margin:34px 0; }}
h2 {{ font-size:24px; margin-bottom:12px; }}
.note {{ border-left:4px solid var(--gold); background:#101d2d; padding:14px 16px; color:var(--muted); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:12px; }}
.card img {{ width:100%; background:#fff; border-radius:6px; }}
.card p {{ color:var(--muted); }}
table {{ width:100%; border-collapse:collapse; font-size:13px; background:var(--panel); }}
th, td {{ border:1px solid var(--line); padding:8px; vertical-align:top; }}
th {{ color:var(--gold); text-align:left; position:sticky; top:0; background:#101d2d; }}
.tablewrap {{ overflow:auto; border:1px solid var(--line); }}
a {{ color:var(--gold); }}
code {{ color:var(--gold); }}
@media (max-width:900px) {{ .stats {{ grid-template-columns:repeat(2,1fr); }} }}
</style></head><body>
<header>
<div class="kicker">CROSS-Neo · Moderna/Merck BD readiness · IP-safe package</div>
<h1>Position CROSS-Neo as the neoantigen selection layer.</h1>
<p class="lead">The pitch is not to replace an mRNA platform. The pitch is to improve candidate selection under a fixed individualized vaccine payload budget, then prove it on blinded sponsor data under NDA.</p>
<div class="stats">
<div class="stat"><b>{kg['all_AUPRC']:.3f}</b><span>same-board KG-GA AUPRC</span></div>
<div class="stat"><b>{kg['all_AUROC']:.3f}</b><span>same-board AUROC</span></div>
<div class="stat"><b>{kg['frozen_validation_AUPRC']:.3f}</b><span>frozen validation-like AUPRC</span></div>
<div class="stat"><b>34</b><span>V940-like payload slot frame</span></div>
<div class="stat"><b>NDA</b><span>required before internals</span></div>
</div>
</header><main>
<section><h2>Decision</h2>
<p class="note"><b>Ready for IP-safe BD packaging, not ready for unprotected raw outreach.</b> File invention disclosure/provisional first, send only the non-confidential teaser pre-NDA, then offer a blinded top-34 validation collaboration with option/license rights.</p></section>
<section><h2>Figure Board</h2><div class="grid">{cards}</div></section>
<section><h2>Same-Board Evidence</h2><div class="tablewrap">{table_html(same, ['algorithm','algorithm_family','all_AUPRC','all_AUROC','frozen_validation_AUPRC','low_medium_leakage_AUPRC','v6_smoke_AUPRC','primary_disposition'], 20)}</div></section>
<section><h2>Strict Reliability Context</h2><div class="tablewrap">{table_html(strict, ['algorithm','algorithm_family','claim_track','neo_strict_AUPRC','neo_strict_AUROC','patients_evaluated','primary_disposition'], 20)}</div></section>
<section><h2>Broader Public Competitor Context</h2><div class="tablewrap">{table_html(public, ['method','n_context_rows','total_scored','datasets_with_auroc','mean_context_AUROC','best_context_AUROC','min_context_AUROC','claim_use'], 20)}</div></section>
<section><h2>Validation Plan</h2><div class="tablewrap">{table_html(validation, ['module','required_from_partner','our_commitment','pass_signal'], 20)}</div></section>
<section><h2>BD Readiness</h2><div class="tablewrap">{table_html(readiness, ['stage','status','deliverable','exit_criterion','share_level'], 20)}</div></section>
<section><h2>Claim Boundary</h2><div class="tablewrap">{table_html(claim_boundary, ['claim','wording','why','do_not_say'], 20)}</div></section>
<section><h2>Data Room Index</h2><div class="tablewrap">{table_html(data_room, ['folder','item','share_timing','contains'], 20)}</div></section>
<section><h2>Contact Routing</h2><div class="tablewrap">{table_html(contact_routes, ['company','route','url','why','pre_nda_action','confidentiality_note'], 20)}</div></section>
<section><h2>Official Context Sources</h2><ul>{sources}</ul></section>
<section><h2>Contact</h2><p class="note"><b>{esc(AUTHOR_NAME)}</b><br>Email: <code>{esc(CONTACT_EMAIL)}</code><br>X/Twitter: <code>{esc(CONTACT_HANDLE)}</code></p></section>
<section><h2>Local Files</h2><ul>
<li><code>{esc(str(OUT))}</code></li>
<li><code>{esc(str(OUT / 'MODERNA_NONCONFIDENTIAL_TEASER.md'))}</code></li>
<li><code>{esc(str(OUT / 'MODERNA_BLINDED_VALIDATION_PROTOCOL.md'))}</code></li>
<li><code>{esc(str(OUT / 'MODERNA_OPTION_LICENSE_TERM_SHEET_SCAFFOLD.md'))}</code></li>
<li><code>{esc(str(OUT / 'pre_nda_send_packet.zip'))}</code></li>
</ul></section>
</main></body></html>"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ASSET.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET.mkdir(parents=True, exist_ok=True)

    same, strict, public, claim = build_metrics_tables()
    readiness, validation, claim_boundary, data_room = build_bd_tables(same, strict, public)
    contact_routes = build_contact_routes()

    same.to_csv(OUT / "moderna_same_board_comparator_metrics.tsv", sep="\t", index=False)
    strict.to_csv(OUT / "moderna_strict_reliability_metrics.tsv", sep="\t", index=False)
    public.to_csv(OUT / "moderna_public_competitor_context.tsv", sep="\t", index=False)
    readiness.to_csv(OUT / "moderna_bd_readiness_ladder.tsv", sep="\t", index=False)
    validation.to_csv(OUT / "moderna_blinded_validation_success_criteria.tsv", sep="\t", index=False)
    claim_boundary.to_csv(OUT / "moderna_claim_boundary.tsv", sep="\t", index=False)
    data_room.to_csv(OUT / "moderna_data_room_index.tsv", sep="\t", index=False)
    contact_routes.to_csv(OUT / "moderna_contact_routes.tsv", sep="\t", index=False)
    claim.to_csv(OUT / "source_ga_rl_claim_boundary_matrix.tsv", sep="\t", index=False)

    write_markdown_docs(same, strict, public, readiness, validation, claim_boundary, data_room, contact_routes)
    build_pre_nda_packet(contact_routes)
    build_figures(same, strict, readiness, validation)

    for src in FIG.glob("*"):
        if src.is_file():
            safe_copy(src, ASSET / src.name)
            safe_copy(src, LIVE_ASSET / src.name)

    for src in OUT.glob("*"):
        if src.is_file():
            safe_copy(src, ASSET / src.name)
            safe_copy(src, LIVE_ASSET / src.name)

    PAGE.write_text(build_html(same, strict, public, readiness, validation, claim_boundary, data_room, contact_routes), encoding="utf-8")
    safe_copy(PAGE, LIVE_PAGE)

    kg = same[same["algorithm"].eq("KG_GA_evolved")].iloc[0]
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": str(OUT),
        "html_path": str(PAGE),
        "live_html_path": str(LIVE_PAGE),
        "kg_ga_same_board_auprc": float(kg["all_AUPRC"]),
        "kg_ga_same_board_auroc": float(kg["all_AUROC"]),
        "kg_ga_frozen_validation_auprc": float(kg["frozen_validation_AUPRC"]),
        "package_status": "BD-ready scaffold; file IP/provisional before external outreach",
        "pre_nda_rule": "share teaser only; no code, full feature list, candidate tables, or weights",
        "pre_nda_send_packet_zip": str(OUT / "pre_nda_send_packet.zip"),
    }
    (OUT / "moderna_bd_package_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
