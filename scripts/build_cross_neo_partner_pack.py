#!/usr/bin/env python3
"""Build a partner-ready CROSS-Neo pack for external outreach.

This package is meant to sit above the executive brief: it combines the
selection-layer pitch, contact routing, latest SOTA watchlist, and a clean
send packet index in one shareable place.
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
MOD = ROOT / "project/results/cross_neo_moderna_bd_package_2026_05_11"
SOTA = ROOT / "project/results/cross_neo_latest_sota_watchlist_2026_05_11"
OUT = ROOT / "project/results/cross_neo_partner_pack_2026_05_11"
FIG = OUT / "figures"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = HUB / "assets/cross_neo_partner_pack"
LIVE_ASSET = LIVE / "assets/cross_neo_partner_pack"
PAGE = HUB / "cross_neo_partner_pack.html"
LIVE_PAGE = LIVE / "cross_neo_partner_pack.html"


def read_tsv(path: Path) -> pd.DataFrame:
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
        print(f"[partner-pack] skip permission-denied copy {dst}")


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


def build_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    master = read_tsv(GA / "ga_rl_master_algorithm_comparison.tsv")
    public = read_tsv(GA / "ga_rl_public_competitor_context_summary.tsv")
    strict = master[master["algorithm"].isin([
        "BAR-Neo_confidence", "BAR-Neo", "BAR-Neo_patient_gated", "BAR-Neo-X",
        "BAR-Neo-X_claim_safe", "RF_biophys", "MHCflurry_2.0_presentation", "KG_GA_evolved_controller",
    ])].copy()
    same = master[master["algorithm"].isin([
        "KG_GA_evolved", "CROSS_integrated", "CROSS_claimsafe", "CROSS_stress",
        "CROSS_BMA", "CROSS_finetuned", "BigMHC_IM", "BigMHC_EL",
    ])].copy()
    contact = read_tsv(MOD / "moderna_contact_routes.tsv")
    return same, strict, public, contact, master


def build_pack_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    exec_brief = pd.DataFrame(
        [
            {"section": "Same-board champion", "value": "KG_GA_evolved", "detail": "AUPRC 0.933 / AUROC 0.943 / top96 96/96"},
            {"section": "Frozen validation-like", "value": "KG_GA_evolved", "detail": "AUPRC 0.978 / AUROC 0.999"},
            {"section": "Low/medium leakage", "value": "KG_GA_evolved", "detail": "AUPRC 0.757 / top96 36/96"},
            {"section": "Strict no-overlap winner", "value": "BAR-Neo_confidence", "detail": "AUPRC 0.857 / AUROC 0.980"},
            {"section": "Latest literature watchlist", "value": "NetMHCpan-4.2", "detail": "plus ImmugenX, PepFore/NeoaPred, NeoGuider, TrambaHLApan, NeoTImmuML"},
        ]
    )
    send_packet = pd.DataFrame(
        [
            {"item": "Non-confidential teaser", "share_level": "pre-NDA", "path": str(MOD / "MODERNA_NONCONFIDENTIAL_TEASER.md")},
            {"item": "Outreach email draft", "share_level": "pre-NDA", "path": str(MOD / "MODERNA_OUTREACH_EMAIL_AFTER_IP.md")},
            {"item": "Impact memo", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_IMPACT_MEMO.md")},
            {"item": "Board memo", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_BOARD_MEMO.md")},
            {"item": "Competitive positioning", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_COMPETITIVE_POSITIONING.md")},
            {"item": "Next-step checklist", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_NEXT_STEP_CHECKLIST.md")},
            {"item": "Why now", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_WHY_NOW.md")},
            {"item": "Why sponsor", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_WHY_SPONSOR.md")},
            {"item": "Why team", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_WHY_TEAM.md")},
            {"item": "One-slide executive", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_ONE_SLIDE_EXECUTIVE.md")},
            {"item": "Killer summary", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_KILLER_SUMMARY.md")},
            {"item": "Decision snapshot", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_DECISION_SNAPSHOT.md")},
            {"item": "Cover letter", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_COVER_LETTER.md")},
            {"item": "Partner invitation", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_PARTNER_INVITATION.md")},
            {"item": "Board ask", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_BOARD_ASK.md")},
            {"item": "Pilot proposal", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_PILOT_PROPOSAL.md")},
            {"item": "Success metrics", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_SUCCESS_METRICS.md")},
            {"item": "Milestones", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_MILESTONES.md")},
            {"item": "Budget sketch", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_BUDGET_SKETCH.md")},
            {"item": "Mutual action plan", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_MUTUAL_ACTION_PLAN.md")},
            {"item": "Data room checklist", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_DATA_ROOM_CHECKLIST.md")},
            {"item": "Business case", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_BUSINESS_CASE.md")},
            {"item": "Approval path", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_APPROVAL_PATH.md")},
            {"item": "Integration checklist", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_INTEGRATION_CHECKLIST.md")},
            {"item": "Term sheet scaffold", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_TERM_SHEET_SCAFFOLD.md")},
            {"item": "ROI model", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_ROI_MODEL.md")},
            {"item": "Operating model", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_SPONSOR_OPERATING_MODEL.md")},
            {"item": "Investment memo", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_INTERNAL_INVESTMENT_MEMO.md")},
            {"item": "Pilot charter", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_PILOT_CHARTER.md")},
            {"item": "Readout template", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_EXEC_READOUT_TEMPLATE.md")},
            {"item": "Launch checklist", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_LAUNCH_CHECKLIST.md")},
            {"item": "RACI", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_RACI.md")},
            {"item": "Board packet", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_BOARD_PACKET.md")},
            {"item": "SOW scaffold", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_SOW_SCAFFOLD.md")},
            {"item": "Compliance checklist", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_COMPLIANCE_CHECKLIST.md")},
            {"item": "Exec decision memo", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_EXEC_DECISION_MEMO.md")},
            {"item": "Investment committee brief", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_INVESTMENT_COMMITTEE_BRIEF.md")},
            {"item": "Data room index", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_DATA_ROOM_INDEX.md")},
            {"item": "Platform guardrails", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_PLATFORM_GUARDRAILS.md")},
            {"item": "Executive Q&A", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_EXEC_QA.md")},
            {"item": "Stakeholder map", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_STAKEHOLDER_MAP.md")},
            {"item": "Sign-off memo", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_SIGNOFF_MEMO.md")},
            {"item": "Talking points", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_TALKING_POINTS.md")},
            {"item": "Launch deck outline", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_LAUNCH_DECK_OUTLINE.md")},
            {"item": "Board summary", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_BOARD_SUMMARY.md")},
            {"item": "Deal terms snapshot", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_DEAL_TERMS_SNAPSHOT.md")},
            {"item": "Risk mitigation", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_RISK_MITIGATION_ONE_PAGER.md")},
            {"item": "Appendix index", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_APPENDIX_INDEX.md")},
            {"item": "CEO brief", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_CEO_BRIEF.md")},
            {"item": "Deal snapshot", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_DEAL_SNAPSHOT.md")},
            {"item": "FAQ index", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_FAQ_INDEX.md")},
            {"item": "Sign-off tracker", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_SIGNOFF_TRACKER.md")},
            {"item": "Executive cover", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_EXECUTIVE_COVER.md")},
            {"item": "Value calculation", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_VALUE_CALCULATION.md")},
            {"item": "Fast-path plan", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_FAST_PATH_PLAN.md")},
            {"item": "Call script", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_CALL_SCRIPT.md")},
            {"item": "First meeting agenda", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_FIRST_MEETING_AGENDA.md")},
            {"item": "Due diligence questions", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_DD_QUESTIONS.md")},
            {"item": "Follow-up email", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_FOLLOWUP_EMAIL.md")},
            {"item": "FAQ / objection sheet", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_FAQ_OBJECTION_SHEET.md")},
            {"item": "One-page pitch", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_ONE_PAGE_PITCH.md")},
            {"item": "NDA checklist", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_NDA_CHECKLIST.md")},
            {"item": "Meeting notes template", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_MEETING_NOTES_TEMPLATE.md")},
            {"item": "Risk register", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_RISK_REGISTER.md")},
            {"item": "Contact tracker", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_CONTACT_TRACKER.md")},
            {"item": "Reply matrix", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_REPLY_HANDLING_MATRIX.md")},
            {"item": "Internal handoff", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_INTERNAL_HANDOFF_CHECKLIST.md")},
            {"item": "NDA ask template", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_NDA_ASK_TEMPLATE.md")},
            {"item": "Blinded data request", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_BLINDED_DATA_REQUEST.md")},
            {"item": "Score lock procedure", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_SCORE_LOCK_PROCEDURE.md")},
            {"item": "Post-NDA protocol", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_POST_NDA_PROTOCOL.md")},
            {"item": "Slide outline", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_SLIDE_OUTLINE.md")},
            {"item": "Internal review memo", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_INTERNAL_REVIEW_MEMO.md")},
            {"item": "Go/no-go criteria", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_GO_NO_GO_CRITERIA.md")},
            {"item": "Answer bank", "share_level": "NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_ANSWER_BANK.md")},
            {"item": "Deal memo", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_DEAL_MEMO.md")},
            {"item": "Value thesis", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_VALUE_THESIS.md")},
            {"item": "KPI scorecard", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_KPI_SCORECARD.md")},
            {"item": "Decision tree", "share_level": "pre-NDA", "path": str(ROOT / "project/results/cross_neo_partner_pack_2026_05_11/MODERNA_PARTNER_DECISION_TREE.md")},
            {"item": "Contact routes", "share_level": "routing only", "path": str(MOD / "moderna_contact_routes.tsv")},
            {"item": "Executive brief", "share_level": "internal/shareable summary", "path": str(ROOT / "project/results/cross_neo_exec_brief_2026_05_11/CROSS_NEO_EXEC_BRIEF_KR.md")},
            {"item": "Latest SOTA watchlist", "share_level": "internal/shareable summary", "path": str(SOTA / "latest_literature_sota_watchlist.tsv")},
        ]
    )
    claim = pd.DataFrame(
        [
            {"claim": "Allowed pre-NDA", "line": "Retrospective prioritization engine and validation ask", "blocked": "Efficacy or prospective proof"},
            {"claim": "Allowed under NDA", "line": "Blinded top-34 validation on sponsor data", "blocked": "Unblinded model internals before lock"},
            {"claim": "Allowed in BD", "line": "Selection layer for individualized mRNA vaccine payloads", "blocked": "Replacing Moderna's platform"},
            {"claim": "Not allowed", "line": "Universal SOTA claim over every public method", "blocked": "Different endpoints and datasets"},
        ]
    )
    return exec_brief, send_packet, claim


def build_figures(same: pd.DataFrame, strict: pd.DataFrame, public: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.style.use("default")

    same_plot = same.copy().sort_values("all_AUPRC")
    fig, ax = plt.subplots(figsize=(10.5, 5))
    ax.barh(same_plot["algorithm"], same_plot["all_AUPRC"], color=["#c7922b" if a == "KG_GA_evolved" else "#315f7d" for a in same_plot["algorithm"]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Same-board AUPRC")
    ax.set_title("Local selection-layer rank")
    fig.tight_layout()
    fig.savefig(FIG / "fig1_same_board_rank.png", dpi=220)
    fig.savefig(FIG / "fig1_same_board_rank.pdf")
    plt.close(fig)

    strict_plot = strict.copy().sort_values("neo_strict_AUPRC")
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.barh(strict_plot["algorithm"], strict_plot["neo_strict_AUPRC"], color=["#b24a62" if "KG_GA" in a else "#315f7d" for a in strict_plot["algorithm"]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Strict no-overlap AUPRC")
    ax.set_title("Reliability lane stays separate")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_strict_reliability.png", dpi=220)
    fig.savefig(FIG / "fig2_strict_reliability.pdf")
    plt.close(fig)

    pub = public.copy().sort_values("best_context_AUROC", ascending=False)
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.barh(pub["method"], pub["best_context_AUROC"], color=["#c7922b" if m == "RF_biophys" else "#315f7d" for m in pub["method"]])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Best context AUROC")
    ax.set_title("Broader public comparator context")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_public_context.png", dpi=220)
    fig.savefig(FIG / "fig3_public_context.pdf")
    plt.close(fig)


def build_html(same: pd.DataFrame, strict: pd.DataFrame, public: pd.DataFrame, contact: pd.DataFrame, exec_brief: pd.DataFrame, send_packet: pd.DataFrame, claim: pd.DataFrame) -> str:
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>CROSS-Neo Partner Pack</title>
<style>
:root {{ --bg:#07111f; --panel:#101d2d; --ink:#edf5ff; --muted:#9fb0c5; --gold:#f2c46d; --line:#27384f; --red:#ff7d8d; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:Inter,Arial,sans-serif; line-height:1.5; }}
header {{ padding:48px 5vw 28px; border-bottom:1px solid var(--line); background:#0b1626; }}
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
table {{ width:100%; border-collapse:collapse; font-size:13px; background:var(--panel); }}
th, td {{ border:1px solid var(--line); padding:8px; vertical-align:top; }}
th {{ color:var(--gold); text-align:left; position:sticky; top:0; background:#101d2d; }}
.tablewrap {{ overflow:auto; border:1px solid var(--line); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:16px; }}
.card {{ border:1px solid var(--line); background:var(--panel); border-radius:8px; padding:12px; }}
.card img {{ width:100%; background:#fff; border-radius:6px; }}
.card p {{ color:var(--muted); }}
a {{ color:var(--gold); }}
code {{ color:var(--gold); }}
.muted {{ color:var(--muted); }}
.bad {{ color:var(--red); }}
@media (max-width:900px) {{ .stats {{ grid-template-columns:repeat(2,1fr); }} }}
</style></head><body>
<header>
<h1>One partner page for selection, SOTA context, and outreach.</h1>
<p class="lead">This pack is the bridge between the technical comparison and the external conversation. It keeps same-board champion evidence, strict reliability boundaries, and latest literature watchlist separated so the pitch stays credible.</p>
<div class="stats">
<div class="stat"><b>0.933</b><span>KG-GA same-board AUPRC</span></div>
<div class="stat"><b>0.857</b><span>strict no-overlap winner</span></div>
<div class="stat"><b>7</b><span>latest SOTA literature groups</span></div>
<div class="stat"><b>pre-NDA</b><span>teaser only before NDA</span></div>
<div class="stat"><b>BD-ready</b><span>selection-layer framing</span></div>
</div>
</header><main>
<section><h2>Positioning</h2><p class="note"><b>We are not selling a replacement mRNA platform.</b> We are selling the selection layer that improves which candidates enter a fixed-slot individualized vaccine payload, with a frozen blinded validation path and clear claim boundary.</p></section>
<section><h2>Visual Board</h2><div class="grid">
<article class="card"><img src="assets/cross_neo_partner_pack/fig1_same_board_rank.png" alt="Same-board rank"><p>KG-GA vs local comparators.</p></article>
<article class="card"><img src="assets/cross_neo_partner_pack/fig2_strict_reliability.png" alt="Strict reliability"><p>BAR-Neo remains the strict lane winner.</p></article>
<article class="card"><img src="assets/cross_neo_partner_pack/fig3_public_context.png" alt="Public context"><p>Broader public comparator context.</p></article>
</div></section>
<section><h2>Executive Brief</h2><div class="tablewrap">{table_html(exec_brief, ['section','value','detail'], 20)}</div></section>
<section><h2>Same-Board Comparator Table</h2><div class="tablewrap">{table_html(same, ['algorithm','all_AUPRC','all_AUROC','frozen_validation_AUPRC','low_medium_leakage_AUPRC','v6_smoke_AUPRC','paperability_reward'], 20)}</div></section>
<section><h2>Strict Reliability Table</h2><div class="tablewrap">{table_html(strict, ['algorithm','neo_strict_AUPRC','neo_strict_AUROC'], 20)}</div></section>
<section><h2>Broader Public Context</h2><div class="tablewrap">{table_html(public, ['method','mean_context_AUROC','best_context_AUROC','claim_use'], 20)}</div></section>
<section><h2>Latest SOTA Watchlist</h2>
<p class="note">NetMHCpan-4.2, ImmugenX, PepFore/NeoaPred, NeoGuider, TrambaHLApan, and NeoTImmuML are the newer literature methods to watch. They are useful context, but their metrics are not merged into the local leaderboard because the datasets and endpoints differ.</p>
</section>
<section><h2>Send Packet</h2><div class="tablewrap">{table_html(send_packet, ['item','share_level','path'], 20)}</div></section>
<section><h2>Contact Routing</h2><div class="tablewrap">{table_html(contact, ['company','route','url','pre_nda_action','confidentiality_note'], 20)}</div></section>
<section><h2>Claim Boundary</h2><div class="tablewrap">{table_html(claim, ['claim','line','blocked'], 20)}</div></section>
<section><h2>What to Say</h2>
<div class="card">
<p><b>Say:</b> We have a frozen selection layer that improves candidate prioritization under a fixed-slot individualized mRNA vaccine payload.</p>
<p><b>Say:</b> We can run a blinded validation on sponsor candidates and return top-34 ranking plus uncertainty and assay-control annotations.</p>
<p><b>Do not say:</b> Prospective efficacy, universal SOTA, or replacement of the sponsor platform.</p>
</div></section>
<section><h2>Conversation Assets</h2>
<div class="tablewrap">{table_html(send_packet[send_packet['item'].isin(['Call script','First meeting agenda','Due diligence questions','Follow-up email'])], ['item','share_level','path'], 20)}</div></section>
</main></body></html>"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ASSET.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET.mkdir(parents=True, exist_ok=True)

    same, strict, public, contact, master = build_inputs()
    exec_brief, send_packet, claim = build_pack_tables()

    call_script = [
        "# Moderna/Merck Call Script",
        "",
        "Opening:",
        "Thanks for taking the time. I am reaching out with an IP-protected neoantigen selection layer for individualized mRNA vaccine workflows.",
        "",
        "One sentence value:",
        "In a 2,715-candidate local benchmark, the KG-guided GA selector improved same-board prioritization; the proposal is to test whether it improves blinded top-N candidate selection under a fixed payload budget.",
        "",
        "What to ask:",
        "- Who owns BD&L or external innovation for individualized cancer vaccine selection?",
        "- Would your team consider a CDA/NDA-enabled blinded validation on sponsor candidate pools?",
        "- What is the preferred format for a non-confidential teaser?",
        "- Is there a fixed-slot top-34 / top-30 / top-20 payload convention we should mirror?",
        "",
        "What not to disclose:",
        "No code, no candidate tables, no feature list, no scoring weights, no training joins, no prospective efficacy claims.",
        "",
        "Close:",
        "If this sounds relevant, I can send a one-page non-confidential teaser and then move to CDA/NDA for a blinded validation protocol.",
        "",
    ]
    (OUT / "MODERNA_CALL_SCRIPT.md").write_text("\n".join(call_script) + "\n", encoding="utf-8")

    agenda = [
        "# First Meeting Agenda",
        "",
        "1. Confirm the sponsor's individualized neoantigen workflow and payload slot budget.",
        "2. Explain the selection-layer value proposition in one slide or one paragraph.",
        "3. Review the pre-NDA teaser only.",
        "4. Align on blinded validation inputs, endpoints, and label-lock procedure.",
        "5. Confirm what can be sent under CDA/NDA and who signs off.",
        "6. Set next checkpoint: protocol exchange, frozen run, and decision timeline.",
        "",
    ]
    (OUT / "MODERNA_FIRST_MEETING_AGENDA.md").write_text("\n".join(agenda) + "\n", encoding="utf-8")

    ddq = [
        "# Due Diligence Questions",
        "",
        "- What is the exact payload slot budget per patient?",
        "- Do you want presentation-only, immunogenicity-only, or a dual-lane ranking?",
        "- Can we use blinded candidate pools with labels withheld until score lock?",
        "- Do you prefer top-20, top-30, top-34, or another fixed-slot analysis?",
        "- What false-positive controls and assay formats matter most?",
        "- Is the sponsor stack already using public baselines such as NetMHCpan, MHCflurry, BigMHC, PRIME, or TransPHLA?",
        "- What are the approval gates for moving from blinded validation to option/license discussion?",
        "",
    ]
    (OUT / "MODERNA_DD_QUESTIONS.md").write_text("\n".join(ddq) + "\n", encoding="utf-8")

    followup = [
        "# Follow-Up Email Draft",
        "",
        "Subject: Blinded validation proposal for individualized neoantigen selection",
        "",
        "Dear [Name],",
        "",
        "Thank you for the discussion. As a follow-up, I am attaching only the non-confidential teaser and routing materials.",
        "",
        "The key proposal remains the same: a frozen, blinded validation on sponsor candidate pools to test whether KG-guided neoantigen selection improves fixed-slot ranking under your current workflow.",
        "",
        "If useful, I can send a CDA/NDA-stage validation protocol and the full score lock procedure once your team confirms the appropriate contact and confidentiality path.",
        "",
        "Best,",
        "Seungho Cook",
        "kukshomr@snu.ac.kr",
        "",
    ]
    (OUT / "MODERNA_FOLLOWUP_EMAIL.md").write_text("\n".join(followup) + "\n", encoding="utf-8")

    faq = [
        "# FAQ / Objection Handling",
        "",
        "## Why not just use existing public predictors?",
        "",
        "Public predictors are useful baselines, but they are not optimized for our exact fixed-slot selection problem, blinded validation workflow, or dual-lane claim-safe gating.",
        "",
        "## Why should a partner care?",
        "",
        "Because the hard problem is not generating lots of candidate peptides; it is selecting the right small set under a payload budget with fewer false positives and a defensible audit trail.",
        "",
        "## Why NDA?",
        "",
        "Because the useful details are the frozen score logic, candidate-level rankings, and validation protocol. Those should not be disclosed before IP hygiene and confidentiality terms.",
        "",
        "## Why is this different from a standard immunogenicity model?",
        "",
        "It is a selection layer plus reliability lane, not a single monolithic predictor. The model is designed to route candidates, not just score them.",
        "",
        "## Why does the latest literature watchlist matter?",
        "",
        "It shows we are not ignoring newer methods such as NetMHCpan-4.2, ImmugenX, PepFore/NeoaPred, NeoGuider, TrambaHLApan, and NeoTImmuML. We are comparing them by task and keeping metrics separate when they are not directly comparable.",
        "",
        "## What do we ask for in the first partner call?",
        "",
        "A routing contact, a seat for the blinded validation discussion, and a clear NDA path.",
        "",
    ]
    (OUT / "MODERNA_FAQ_OBJECTION_SHEET.md").write_text("\n".join(faq) + "\n", encoding="utf-8")

    one_page = [
        "# One-Page Pitch",
        "",
        "## What We Do",
        "We improve which neoantigen candidates get into a fixed-slot individualized mRNA vaccine payload.",
        "",
        "## Why It Matters",
        "Payload space is limited; better ranking reduces false positives and makes wetlab validation more efficient.",
        "",
        "## Local Evidence",
        "- KG_GA_evolved same-board AUPRC 0.933 / AUROC 0.943.",
        "- Frozen validation-like AUPRC 0.978.",
        "- Public comparators on the same board include BigMHC_IM 0.672 and BigMHC_EL 0.648.",
        "",
        "## Reliability Boundary",
        "- BAR-Neo_confidence remains the strict no-overlap winner at AUPRC 0.857.",
        "",
        "## Why Partner",
        "- Frozen blinded validation on sponsor candidates.",
        "- Top-34 ranking under a V940-like slot budget.",
        "- Audit trail, uncertainty, and assay-control suggestions.",
        "",
        "## Ask",
        "A short intro, a CDA/NDA, and a blinded benchmark set.",
        "",
    ]
    (OUT / "MODERNA_ONE_PAGE_PITCH.md").write_text("\n".join(one_page) + "\n", encoding="utf-8")

    nda_checklist = [
        "# NDA Checklist",
        "",
        "- Confirm legal entity names for the CDA/NDA.",
        "- Confirm mutual confidentiality and permitted-use language.",
        "- Confirm whether score outputs are considered confidential data.",
        "- Confirm retention, destruction, and publication-review language.",
        "- Confirm whether sponsor can use outputs for internal benchmarking only or for model training.",
        "- Confirm dispute / audit rights for blinded validation runs.",
        "- Confirm whether the validation set can include decoys and WT peptides.",
        "",
    ]
    (OUT / "MODERNA_NDA_CHECKLIST.md").write_text("\n".join(nda_checklist) + "\n", encoding="utf-8")

    meeting_notes = [
        "# Meeting Notes Template",
        "",
        "## Attendees",
        "-",
        "",
        "## Sponsor Workflow",
        "-",
        "",
        "## Payload Slot Budget",
        "-",
        "",
        "## Current Baselines",
        "-",
        "",
        "## Validation Inputs",
        "-",
        "",
        "## NDA Path",
        "-",
        "",
        "## Action Items",
        "-",
        "",
    ]
    (OUT / "MODERNA_MEETING_NOTES_TEMPLATE.md").write_text("\n".join(meeting_notes) + "\n", encoding="utf-8")

    risk = [
        "# Risk Register",
        "",
        "| risk | mitigation |",
        "|---|---|",
        "| Overclaiming SOTA | Keep same-board, strict, and literature contexts separate |",
        "| Leaking internals pre-NDA | Use teaser-only send packet |",
        "| Sponsor wants model training use | Negotiate permitted-use and licensing upfront |",
        "| Different payload slot conventions | Ask in first call and mirror sponsor workflow |",
        "| Wetlab ambiguity | Require mutant/WT/decoy control mapping |",
        "| Output not adopted | Offer blinded validation and audit trail, not hype |",
        "",
    ]
    (OUT / "MODERNA_RISK_REGISTER.md").write_text("\n".join(risk) + "\n", encoding="utf-8")

    contact_tracker = [
        "# Contact Tracker",
        "",
        "| date | company | contact | route | status | next action |",
        "|---|---|---|---|---|---|",
        "| 2026-05-11 | Merck | TBD | BD&L | not sent | send teaser + ask for CDA/NDA contact |",
        "| 2026-05-11 | Moderna | TBD | routing page | not sent | request routing to BD or external innovation |",
        "",
    ]
    (OUT / "MODERNA_CONTACT_TRACKER.md").write_text("\n".join(contact_tracker) + "\n", encoding="utf-8")

    reply_matrix = [
        "# Reply Handling Matrix",
        "",
        "| reply type | meaning | response |",
        "|---|---|---|",
        "| interest | wants to learn more | send teaser, then NDA request |",
        "| intro request | wants to route internally | send one-page pitch and contact tracker |",
        "| NDA request | confidentiality path available | send NDA checklist + blinded validation ask |",
        "| technical questions | wants specifics | answer only at teaser level unless NDA signed |",
        "| no response | low priority | follow up once after 5-7 business days |",
        "| decline | not a fit now | thank them, preserve route, ask for future contact |",
        "",
    ]
    (OUT / "MODERNA_REPLY_HANDLING_MATRIX.md").write_text("\n".join(reply_matrix) + "\n", encoding="utf-8")

    handoff = [
        "# Internal Handoff Checklist",
        "",
        "- Confirm teaser version is current and does not include code or candidate tables.",
        "- Confirm NDA checklist is attached before any protocol details.",
        "- Confirm call notes were captured and entered into tracker.",
        "- Confirm follow-up email was sent or scheduled.",
        "- If NDA path exists, move to blinded validation protocol and legal review.",
        "- If partner requests deeper technical detail, stop and route through IP-approved channel.",
        "",
    ]
    (OUT / "MODERNA_INTERNAL_HANDOFF_CHECKLIST.md").write_text("\n".join(handoff) + "\n", encoding="utf-8")

    nda_ask = [
        "# NDA Ask Template",
        "",
        "We would like to continue the technical discussion under a CDA/NDA so we can share the blinded validation protocol and discuss sponsor-provided candidate pools.",
        "",
        "Requested terms:",
        "- Mutual confidentiality",
        "- Permitted use limited to evaluation of the selection-layer technology",
        "- No training or productization use of submitted scores without a separate license",
        "- Confidential treatment of score outputs and validation materials",
        "",
    ]
    (OUT / "MODERNA_NDA_ASK_TEMPLATE.md").write_text("\n".join(nda_ask) + "\n", encoding="utf-8")

    blinded_request = [
        "# Blinded Data Request",
        "",
        "Please provide a blinded candidate table with the following where available:",
        "",
        "- patient_id",
        "- mutant_peptide",
        "- optional WT peptide",
        "- HLA allele",
        "- expression / RNA support",
        "- variant context",
        "- payload slot budget / intended rank window",
        "- sponsor baseline score if it can remain blinded until lock",
        "- eventual labels withheld until score freeze",
        "",
        "Optional controls:",
        "- decoys",
        "- self-like peptides",
        "- known non-binders",
        "- wetlab-negative historical candidates",
        "",
        "Deliverable after freeze:",
        "- per-candidate score",
        "- per-patient top-N ranking",
        "- uncertainty flags",
        "- reason codes",
        "- assay-control annotations",
        "",
    ]
    (OUT / "MODERNA_BLINDED_DATA_REQUEST.md").write_text("\n".join(blinded_request) + "\n", encoding="utf-8")

    lock_proc = [
        "# Score Lock Procedure",
        "",
        "1. Freeze model version and config hash.",
        "2. Freeze input parser and feature schema.",
        "3. Generate scores on blinded sponsor data.",
        "4. Store outputs with timestamp and hash.",
        "5. Do not inspect labels until the lock is documented.",
        "6. After lock, compare against sponsor baseline and unblinded outcomes.",
        "",
        "Hard rule:",
        "No tuning, no calibration, no feature changes, no threshold updates after seeing labels.",
        "",
    ]
    (OUT / "MODERNA_SCORE_LOCK_PROCEDURE.md").write_text("\n".join(lock_proc) + "\n", encoding="utf-8")

    post_nda = [
        "# Post-NDA Validation Protocol",
        "",
        "## Phase 1: Freeze",
        "- receive blinded candidate file",
        "- confirm slot budget",
        "- confirm endpoint definitions",
        "",
        "## Phase 2: Run",
        "- execute frozen scorer",
        "- generate top-34 and optional top-20, top-30 views",
        "- emit uncertainty and reason codes",
        "",
        "## Phase 3: Unblind",
        "- compare to sponsor labels or wetlab outputs",
        "- calculate AUPRC, AUROC, precision@N, and false-positive burden",
        "",
        "## Phase 4: Decision",
        "- go to option/license if lift is meaningful",
        "- otherwise archive the run and keep route warm",
        "",
    ]
    (OUT / "MODERNA_POST_NDA_PROTOCOL.md").write_text("\n".join(post_nda) + "\n", encoding="utf-8")

    slide_outline = [
        "# Slide Outline",
        "",
        "1. Problem: fixed-slot individualized vaccine payloads are constrained.",
        "2. Solution: KG-guided selection layer with audit trail and claim-safe fallback.",
        "3. Local evidence: KG-GA same-board champion, frozen validation-like strength, strict reliability boundary.",
        "4. External context: newer literature methods are improving, but metrics are not directly comparable across datasets.",
        "5. Partner value: blinded top-34 validation and decision support before manufacturing.",
        "6. Ask: CDA/NDA and a blinded benchmark set.",
        "",
    ]
    (OUT / "MODERNA_SLIDE_OUTLINE.md").write_text("\n".join(slide_outline) + "\n", encoding="utf-8")

    internal_review = [
        "# Internal Review Memo",
        "",
        "## Summary",
        "This is a selection-layer technology with a frozen validation path, not a platform replacement.",
        "",
        "## Why It May Matter",
        "- Same-board local evidence is strong.",
        "- The external partner can test blinded top-N performance before committing to deeper diligence.",
        "- The output is operationally useful: rankings, uncertainty, and assay-control annotations.",
        "",
        "## What To Verify Internally",
        "- Is the slot budget aligned with our payload format?",
        "- Are the validation endpoints meaningful for the sponsor program?",
        "- Does the sponsor want a claim-safe lane and a production lane?",
        "- Can legal support CDA/NDA plus evaluation-only use?",
        "",
        "## Open Risks",
        "- Cross-dataset SOTA claims are not appropriate without reruns.",
        "- Prospective efficacy is not yet established.",
        "- Public predictors remain useful baselines and may be preferred in some workflows.",
        "",
    ]
    (OUT / "MODERNA_INTERNAL_REVIEW_MEMO.md").write_text("\n".join(internal_review) + "\n", encoding="utf-8")

    go_no_go = [
        "# Go / No-Go Criteria",
        "",
        "| criterion | go | no-go |",
        "|---|---|---|",
        "| NDA signed | yes | no |",
        "| Blinded data available | yes | no |",
        "| Score lock honored | yes | no |",
        "| Top-N lift over incumbent | meaningful and reproducible | flat or noisy |",
        "| Sponsor fit | individualized payload and selection pain exists | no relevant use case |",
        "| Legal path | evaluation and permitted use clear | training use / ownership unclear |",
        "| Wetlab path | controls and endpoint definitions exist | no actionable assay route |",
        "",
    ]
    (OUT / "MODERNA_GO_NO_GO_CRITERIA.md").write_text("\n".join(go_no_go) + "\n", encoding="utf-8")

    answer_bank = [
        "# Answer Bank",
        "",
        "Q: Are you claiming clinical superiority?",
        "A: No. The claim is bounded to retrospective prioritization and blinded validation.",
        "",
        "Q: Why not just use public baselines?",
        "A: We do compare to public baselines, but they are not optimized for the same fixed-slot selection problem or internal workflow.",
        "",
        "Q: What is your strongest result?",
        "A: KG_GA_evolved on the same 2,715-candidate board at AUPRC 0.933 / AUROC 0.943.",
        "",
        "Q: What is your strict reliability caveat?",
        "A: BAR-Neo_confidence remains the strict no-overlap winner.",
        "",
        "Q: What do you need next?",
        "A: A CDA/NDA, blinded data, and a sponsor contact who owns the evaluation path.",
        "",
    ]
    (OUT / "MODERNA_ANSWER_BANK.md").write_text("\n".join(answer_bank) + "\n", encoding="utf-8")

    executive_cover = [
        "# Executive Cover",
        "",
        "## Why this deserves attention now",
        "The sponsor already has an individualized mRNA vaccine engine. The leverage point is the selection layer: which candidates get into the limited payload window.",
        "",
        "## What we have",
        "- A frozen local champion on the 2,715-candidate board.",
        "- A claim-safe fallback lane.",
        "- A strict reliability lane kept separate.",
        "- A blinded validation protocol and a pre-NDA packet that can be sent immediately.",
        "",
        "## Why this is actionable",
        "- It can be validated on blinded sponsor candidates without changing manufacturing infrastructure.",
        "- It can be framed as evaluation first, license later.",
        "- It has a clean boundary between discovery, reliability, and literature watchlist context.",
        "",
    ]
    (OUT / "MODERNA_EXECUTIVE_COVER.md").write_text("\n".join(executive_cover) + "\n", encoding="utf-8")

    value_calc = [
        "# Value Calculation",
        "",
        "| driver | proxy | why it matters |",
        "|---|---|---|",
        "| fewer false positives entering synthesis | higher selection AUPRC | less wasted wetlab capacity |",
        "| faster shortlist generation | frozen top-N ranking | shorter review cycles |",
        "| better sponsor confidence | audit trail and uncertainty | easier internal signoff |",
        "| better assay design | control annotations | clearer wetlab interpretation |",
        "| easier adoption | selection layer, not platform replacement | lower integration friction |",
        "",
        "## Practical business claim",
        "The goal is not to replace existing manufacturing or public predictors. The goal is to improve the fraction of candidates that deserve expensive validation and payload inclusion.",
        "",
    ]
    (OUT / "MODERNA_VALUE_CALCULATION.md").write_text("\n".join(value_calc) + "\n", encoding="utf-8")

    fast_path = [
        "# Fast-Path Plan",
        "",
        "1. Send teaser + executive cover + contact route.",
        "2. Request one named owner for external innovation or BD&L.",
        "3. Ask for CDA/NDA and a blinded candidate pool.",
        "4. Run frozen score lock and top-N ranking.",
        "5. Compare against sponsor incumbent on pre-specified metrics.",
        "6. If signal is real, move to paid pilot or option discussion.",
        "",
    ]
    (OUT / "MODERNA_FAST_PATH_PLAN.md").write_text("\n".join(fast_path) + "\n", encoding="utf-8")

    deal_memo = [
        "# Deal Memo",
        "",
        "## Thesis",
        "A frozen neoantigen selection layer can improve individualized mRNA vaccine payload selection under a fixed slot budget and reduce the cost of false-positive candidate synthesis and testing.",
        "",
        "## Why It Is Commercially Interesting",
        "- The sponsor already has a manufacturing engine; the scarce resource is payload slot allocation and validation throughput.",
        "- Better selection reduces wasted synthesis, assay time, and downstream ambiguity.",
        "- The value proposition is easiest to validate in blinded top-N ranking, not in an open-ended model bake-off.",
        "",
        "## What Makes It Differentiated",
        "- Same-board local champion is KG_GA_evolved.",
        "- Claim-safe fallback is separated.",
        "- Strict reliability lane is separated.",
        "- Latest literature SOTA is explicitly tracked, not ignored.",
        "",
        "## Why This Could Be A Fast Win",
        "- Sponsor can evaluate it without changing manufacturing infrastructure.",
        "- A blinded validation can be staged before deeper integration.",
        "- One fixed-slot top-N decision can be enough to justify an option or pilot.",
        "",
        "## Ask",
        "CDA/NDA, blinded data, and one sponsor owner for evaluation.",
        "",
    ]
    (OUT / "MODERNA_DEAL_MEMO.md").write_text("\n".join(deal_memo) + "\n", encoding="utf-8")

    impact_memo = [
        "# Impact Memo",
        "",
        "## Bottom line",
        "This is a selection-layer technology for individualized neoantigen mRNA programs. The sponsor keeps its manufacturing engine; we improve which candidates deserve expensive validation and payload inclusion.",
        "",
        "## Why this is high-impact",
        "- Fixed-slot payloads are scarce.",
        "- False positives cost synthesis, assay time, and review time.",
        "- A blinded validation on sponsor data can settle whether the ranking gain is real without changing the sponsor's platform.",
        "",
        "## What makes it unusually strong",
        "- Same-board champion: KG_GA_evolved AUPRC 0.933 / AUROC 0.943.",
        "- Frozen validation-like AUPRC 0.978.",
        "- Strict reliability lane separated: BAR-Neo_confidence AUPRC 0.857.",
        "- Latest literature watchlist explicitly tracked, including NetMHCpan-4.2 and newer immunogenicity models.",
        "",
        "## Why a partner should care now",
        "- A short blinded pilot can produce an adoption signal fast.",
        "- If the sponsor already has a candidate pool, the marginal integration cost is low.",
        "- If the sponsor does not, the same framework still becomes a reusable scoring stack and decision layer.",
        "",
        "## Ask",
        "Send the teaser, route to the right owner, and ask for a CDA/NDA plus a blinded validation set.",
        "",
    ]
    (OUT / "MODERNA_IMPACT_MEMO.md").write_text("\n".join(impact_memo) + "\n", encoding="utf-8")

    board_memo = [
        "# Board Memo",
        "",
        "## Decision in one line",
        "This is a partner-ready selection layer for individualized neoantigen payloads, not a platform replacement and not a prospective efficacy claim.",
        "",
        "## What the board needs to know",
        "- Same-board local champion: KG_GA_evolved AUPRC 0.933 / AUROC 0.943.",
        "- Frozen validation-like strength: AUPRC 0.978.",
        "- Strict no-overlap reliability winner is separate: BAR-Neo_confidence AUPRC 0.857.",
        "- Latest literature SOTA is tracked separately, with NetMHCpan-4.2 and newer immunogenicity models in scope.",
        "",
        "## Why this matters commercially",
        "- Fixed payload slots are the bottleneck.",
        "- Better selection saves wetlab and synthesis spend.",
        "- The route to value is a blinded validation pilot, then a pilot/license decision.",
        "",
        "## Board ask",
        "Approve a pre-NDA outreach package, an NDA path, and a blinded validation discussion with the sponsor.",
        "",
    ]
    (OUT / "MODERNA_BOARD_MEMO.md").write_text("\n".join(board_memo) + "\n", encoding="utf-8")

    comp = [
        "# Competitive Positioning",
        "",
        "| dimension | our position | why it wins |",
        "|---|---|---|",
        "| discovery prioritization | champion | KG-GA leads same-board local ranking |",
        "| conservative fallback | claim-safe lane | keeps risky outputs separated |",
        "| reliability | separate strict lane | BAR-Neo remains benchmark winner there |",
        "| public baseline coverage | broad | BigMHC, RF_biophys, MHCflurry, PRIME, TransPHLA tracked |",
        "| latest SOTA awareness | current | NetMHCpan-4.2, ImmugenX, PepFore/NeoaPred, NeoGuider, TrambaHLApan, NeoTImmuML tracked |",
        "| commercial fit | selection layer | low-friction integration into sponsor workflow |",
        "",
    ]
    (OUT / "MODERNA_COMPETITIVE_POSITIONING.md").write_text("\n".join(comp) + "\n", encoding="utf-8")

    next_steps = [
        "# Next-Step Checklist",
        "",
        "- Send non-confidential teaser and executive cover.",
        "- Route to BD/technical owner.",
        "- Ask for CDA/NDA and blinded candidate set.",
        "- Confirm payload slot budget and endpoint definitions.",
        "- Freeze score lock procedure before seeing labels.",
        "- Run blinded validation and compare against incumbent.",
        "- If signal is meaningful, move to option/pilot term sheet.",
        "",
    ]
    (OUT / "MODERNA_NEXT_STEP_CHECKLIST.md").write_text("\n".join(next_steps) + "\n", encoding="utf-8")

    why_now = [
        "# Why Now",
        "",
        "- Individualized mRNA programs already exist and are maturing.",
        "- The limiting factor is not only manufacturing; it is candidate selection under a fixed payload budget.",
        "- That means a selection layer can create value without waiting for a platform rebuild.",
        "- The best moment to validate is before the sponsor locks in another generation of internal ranking logic.",
        "",
    ]
    (OUT / "MODERNA_WHY_NOW.md").write_text("\n".join(why_now) + "\n", encoding="utf-8")

    why_sponsor = [
        "# Why This Sponsor",
        "",
        "- They already operate in individualized oncology and can understand the value of better payload selection.",
        "- They have the manufacturing and clinical infrastructure to turn a selection win into an actionable program.",
        "- Their existing workflow makes a blinded validation on sponsor candidate pools practical.",
        "- A positive result can move from evaluation to option/license without needing a new platform introduction.",
        "",
    ]
    (OUT / "MODERNA_WHY_SPONSOR.md").write_text("\n".join(why_sponsor) + "\n", encoding="utf-8")

    why_team = [
        "# Why This Team",
        "",
        "- The local package already separates champion, claim-safe fallback, and strict reliability lanes.",
        "- The team has the discipline to keep claim boundaries explicit.",
        "- The team can hand over a blinded validation package instead of an overhyped pitch deck.",
        "- The team is already producing partner-facing artifacts instead of raw analysis only.",
        "",
    ]
    (OUT / "MODERNA_WHY_TEAM.md").write_text("\n".join(why_team) + "\n", encoding="utf-8")

    one_slide = [
        "# One-Slide Executive",
        "",
        "Title: Selection layer for individualized neoantigen payloads",
        "",
        "1. Problem: payload slots are scarce and candidate pools are noisy.",
        "2. Solution: KG-guided selection layer with claim-safe fallback and strict reliability separation.",
        "3. Evidence: same-board champion KG_GA_evolved AUPRC 0.933 / AUROC 0.943; frozen validation-like AUPRC 0.978.",
        "4. Boundary: BAR-Neo_confidence remains the strict no-overlap winner; this is a discovery/prioritization story, not a prospective efficacy claim.",
        "5. Ask: CDA/NDA plus a blinded validation set and a fixed-slot top-N evaluation.",
        "",
    ]
    (OUT / "MODERNA_ONE_SLIDE_EXECUTIVE.md").write_text("\n".join(one_slide) + "\n", encoding="utf-8")

    killer = [
        "# Killer Summary",
        "",
        "We help a sponsor choose the right few neoantigens from a noisy candidate pool for a fixed-slot individualized mRNA vaccine payload.",
        "The local board says the selection layer works: KG_GA_evolved leads same-board prioritization at AUPRC 0.933 / AUROC 0.943, while strict reliability remains separated in the BAR-Neo lane.",
        "The business move is simple: send a non-confidential teaser, ask for CDA/NDA, and run a blinded top-N validation before deeper integration.",
        "",
    ]
    (OUT / "MODERNA_KILLER_SUMMARY.md").write_text("\n".join(killer) + "\n", encoding="utf-8")

    decision_snapshot = [
        "# Decision Snapshot",
        "",
        "| question | answer |",
        "|---|---|",
        "| Is the local champion clear? | yes, KG_GA_evolved |",
        "| Is strict reliability lane separate? | yes, BAR-Neo_confidence |",
        "| Is latest literature tracked? | yes |",
        "| Is pre-NDA packet ready? | yes |",
        "| Is NDA path defined? | yes |",
        "| Is the pitch platform-replacement? | no, selection-layer only |",
        "| Is a blinded validation ready to run? | yes |",
        "",
    ]
    (OUT / "MODERNA_DECISION_SNAPSHOT.md").write_text("\n".join(decision_snapshot) + "\n", encoding="utf-8")

    cover_letter = [
        "# Cover Letter",
        "",
        "Dear [Name],",
        "",
        "I am sharing a non-confidential package for an IP-protected selection layer that improves fixed-slot individualized neoantigen payload selection.",
        "",
        "The local data are strong enough to justify a blinded evaluation discussion, but not enough to make prospective efficacy claims. That is why the package includes a teaser, a one-slide executive, a value thesis, and a pre-NDA send packet.",
        "",
        "If this is relevant to your team, I would like to move to CDA/NDA and a blinded validation set.",
        "",
        "Best,",
        "Seungho Cook",
        "kukshomr@snu.ac.kr",
        "",
    ]
    (OUT / "MODERNA_COVER_LETTER.md").write_text("\n".join(cover_letter) + "\n", encoding="utf-8")

    partner_inv = [
        "# Partner Invitation",
        "",
        "We invite your team to review a selection-layer technology for individualized neoantigen vaccine payloads.",
        "",
        "What we propose:",
        "- Send a non-confidential teaser first.",
        "- If relevant, move to CDA/NDA.",
        "- Run a blinded validation on sponsor candidate pools.",
        "- Compare against incumbent selection logic on pre-specified metrics.",
        "",
        "What you get:",
        "- A frozen, auditable scoring layer.",
        "- Claim-safe fallback and reliability separation.",
        "- A clear path to option/license only if the pilot shows signal.",
        "",
    ]
    (OUT / "MODERNA_PARTNER_INVITATION.md").write_text("\n".join(partner_inv) + "\n", encoding="utf-8")

    board_ask = [
        "# Board Ask",
        "",
        "Approve:",
        "1. Pre-NDA outreach to Merck/Moderna BD routing contacts.",
        "2. Immediate use of teaser-only materials outside the lab.",
        "3. A formal IP/provisional pathway before any disclosure of internals.",
        "4. A blinded validation plan as the default next technical step.",
        "5. A partner-facing storyline that frames the work as a selection layer, not a platform replacement.",
        "",
    ]
    (OUT / "MODERNA_BOARD_ASK.md").write_text("\n".join(board_ask) + "\n", encoding="utf-8")

    pilot_proposal = [
        "# Pilot Proposal",
        "",
        "## Objective",
        "Validate whether the selection layer improves fixed-slot individualized neoantigen payload ranking on blinded sponsor data.",
        "",
        "## Scope",
        "- Frozen scorer only.",
        "- Blinded candidate pool provided by sponsor.",
        "- Evaluate top-20, top-30, and top-34 cutoffs if sponsor uses those windows.",
        "- Compare against incumbent sponsor baseline and any public comparator already in use.",
        "",
        "## Deliverables",
        "- Score file with timestamp and hash.",
        "- Per-patient ranking table.",
        "- Uncertainty and reason code report.",
        "- Assay-control annotations.",
        "",
        "## Expected sponsor effort",
        "- One dataset handoff.",
        "- One technical owner.",
        "- One legal path for CDA/NDA.",
        "",
    ]
    (OUT / "MODERNA_PILOT_PROPOSAL.md").write_text("\n".join(pilot_proposal) + "\n", encoding="utf-8")

    success_metrics = [
        "# Success Metrics",
        "",
        "| metric | target | why it matters |",
        "|---|---|---|",
        "| top-N AUPRC lift | positive vs incumbent | better candidate ranking |",
        "| top-34 hit rate | higher than incumbent | better fixed-slot utility |",
        "| false-positive burden | lower at matched sensitivity | fewer wasted assays |",
        "| calibration | stable enough for review | easier internal adoption |",
        "| auditability | reason codes accepted | smoother governance |",
        "| turnaround | frozen run completed quickly | fast decision-making |",
        "",
    ]
    (OUT / "MODERNA_SUCCESS_METRICS.md").write_text("\n".join(success_metrics) + "\n", encoding="utf-8")

    milestones = [
        "# Milestones",
        "",
        "1. Teaser delivered.",
        "2. Contact routed to owner.",
        "3. NDA requested and signed.",
        "4. Blinded candidate data received.",
        "5. Score lock completed.",
        "6. Frozen validation run completed.",
        "7. Sponsor review of top-N lift and false-positive burden.",
        "8. Option/pilot discussion or archive.",
        "",
    ]
    (OUT / "MODERNA_MILESTONES.md").write_text("\n".join(milestones) + "\n", encoding="utf-8")

    budget_sketch = [
        "# Budget Sketch",
        "",
        "This is not a quote. It is a rough internal frame for a paid pilot or sponsored evaluation.",
        "",
        "| cost bucket | what it covers | comment |",
        "|---|---|---|",
        "| evaluation setup | data parsing, score lock, report generation | low integration cost |",
        "| technical review | meetings, parameter alignment, validation discussion | one owner at sponsor side helps |",
        "| pilot execution | blinded run on sponsor candidate pool | core value event |",
        "| legal / IP | CDA/NDA, permitted-use language, option terms | counsel-driven |",
        "| post-pilot | adoption / license / milestone work | only if validation passes |",
        "",
    ]
    (OUT / "MODERNA_BUDGET_SKETCH.md").write_text("\n".join(budget_sketch) + "\n", encoding="utf-8")

    mutual_action_plan = [
        "# Mutual Action Plan",
        "",
        "| step | owner | output | timing |",
        "|---|---|---|---|",
        "| 1. route | us | teaser + NDA ask sent | day 0 |",
        "| 2. contact confirmation | sponsor | named owner and legal path | day 1-5 |",
        "| 3. NDA execution | both | CDA/NDA in place | day 5-15 |",
        "| 4. data room setup | sponsor | blinded candidate pool + metadata | day 15-25 |",
        "| 5. score lock | us | frozen model hash + run log | day 25-30 |",
        "| 6. blinded validation | us | ranked output + audit trail | day 30-40 |",
        "| 7. readout | both | decision memo and next-step call | day 40-45 |",
        "",
    ]
    (OUT / "MODERNA_MUTUAL_ACTION_PLAN.md").write_text("\n".join(mutual_action_plan) + "\n", encoding="utf-8")

    data_room = [
        "# Data Room Checklist",
        "",
        "## Sponsor-provided items",
        "- Blinded candidate pool.",
        "- Optional label file held separately until score lock.",
        "- Slot budget policy per patient.",
        "- Any incumbent baseline output currently used.",
        "- Assay controls or wetlab annotations if available.",
        "",
        "## Our-provided items",
        "- Frozen score file with timestamp and hash.",
        "- Ranking table and top-N exports.",
        "- Uncertainty and reason-code report.",
        "- Audit trail summary.",
        "- Readout deck for internal review.",
        "",
        "## Access rule",
        "- No unblinded internals before NDA and score lock.",
        "- No feature list leakage outside the approved scope.",
        "",
    ]
    (OUT / "MODERNA_DATA_ROOM_CHECKLIST.md").write_text("\n".join(data_room) + "\n", encoding="utf-8")

    business_case = [
        "# Business Case",
        "",
        "This is the internal sponsor-facing rationale for why a partner should care now.",
        "",
        "| question | answer |",
        "|---|---|",
        "| What problem is solved? | Better selection under a fixed payload budget for individualized neoantigen workflows. |",
        "| What changes economically? | Fewer weak candidates enter the manufacturing and assay queue. |",
        "| Why now? | The workflow already exists; the selection layer can be evaluated without changing the platform. |",
        "| Why us? | Local same-board champion plus a strict reliability lane and a blinded validation path. |",
        "| What is the next decision? | NDA-backed pilot or archive. |",
        "",
    ]
    (OUT / "MODERNA_BUSINESS_CASE.md").write_text("\n".join(business_case) + "\n", encoding="utf-8")

    approval_path = [
        "# Approval Path",
        "",
        "1. Internal owner confirms relevance.",
        "2. Legal confirms NDA/CDA route.",
        "3. Technical sponsor accepts blinded validation scope.",
        "4. Partner allows a frozen score lock and top-N readout.",
        "5. Sponsor chooses archive, pilot, or term-sheet discussion.",
        "",
    ]
    (OUT / "MODERNA_APPROVAL_PATH.md").write_text("\n".join(approval_path) + "\n", encoding="utf-8")

    integration = [
        "# Integration Checklist",
        "",
        "## Sponsor side",
        "- Confirm candidate pool format.",
        "- Confirm fixed-slot payload convention.",
        "- Confirm label-hold / score-lock procedure.",
        "- Confirm who signs off on validation readout.",
        "",
        "## Our side",
        "- Freeze scorer hash.",
        "- Produce top-N output and uncertainty report.",
        "- Return audit trail and reason codes.",
        "- Keep internals out of the sponsor data room.",
        "",
    ]
    (OUT / "MODERNA_INTEGRATION_CHECKLIST.md").write_text("\n".join(integration) + "\n", encoding="utf-8")

    term_sheet = [
        "# Term Sheet Scaffold",
        "",
        "This is an NDA-stage skeleton only, not a final legal document.",
        "",
        "| section | placeholder |",
        "|---|---|",
        "| scope | blinded validation of a selection layer on sponsor candidate pools |",
        "| deliverable | frozen score file, ranked outputs, audit trail, readout memo |",
        "| decision gate | pilot continuation, option discussion, or archive |",
        "| data handling | sponsor data remains blinded until score lock |",
        "| compensation | pilot / evaluation fee or option-based structure, if applicable |",
        "| exclusions | no platform replacement claim, no unblinded internals pre-clearance |",
        "",
    ]
    (OUT / "MODERNA_TERM_SHEET_SCAFFOLD.md").write_text("\n".join(term_sheet) + "\n", encoding="utf-8")

    roi_model = [
        "# ROI Model",
        "",
        "This is a qualitative internal frame, not a commercial quote.",
        "",
        "| lever | direction | driver |",
        "|---|---|---|",
        "| assay burden | down | fewer weak candidates move forward |",
        "| manufacturing waste | down | better fixed-slot selection |",
        "| review time | down | frozen rankings and reason codes |",
        "| decision quality | up | clearer top-N lift under blinded validation |",
        "| integration cost | low | selection layer, not platform replacement |",
        "",
        "The value argument is strongest when the sponsor already has a candidate pool and wants to reduce false positives before expensive downstream work.",
        "",
    ]
    (OUT / "MODERNA_ROI_MODEL.md").write_text("\n".join(roi_model) + "\n", encoding="utf-8")

    operating_model = [
        "# Sponsor Operating Model",
        "",
        "## Roles",
        "- Business owner: approves confidentiality and pilot scope.",
        "- Technical owner: owns blinded validation inputs and readout.",
        "- Legal owner: routes CDA/NDA and any option language.",
        "- Data owner: handles candidate pool, labels, and lock procedure.",
        "",
        "## Cadence",
        "- One intake call.",
        "- One NDA check-in.",
        "- One data-room transfer.",
        "- One frozen run.",
        "- One decision readout.",
        "",
        "## Operating rule",
        "- Keep the model frozen during validation.",
        "- Keep sponsor labels blinded until lock.",
        "- Keep the next-step decision explicit: archive, pilot, or term-sheet discussion.",
        "",
    ]
    (OUT / "MODERNA_SPONSOR_OPERATING_MODEL.md").write_text("\n".join(operating_model) + "\n", encoding="utf-8")

    investment_memo = [
        "# Internal Investment Memo",
        "",
        "Why this is worth pursuing now:",
        "",
        "- The local evidence already shows a same-board champion and a strict reliability lane.",
        "- The external pitch is low-friction because it does not ask the sponsor to replace its platform.",
        "- The proposed next step is cheap relative to a full platform rewrite: a blinded validation and a score-lock readout.",
        "- The upside is leverage: better ranking, fewer false positives, and a cleaner internal decision path.",
        "",
        "Decision needed:",
        "- Approve pre-NDA outreach.",
        "- Approve NDA-stage validation framing.",
        "- Approve use of the partner pack as the external operating set.",
        "",
    ]
    (OUT / "MODERNA_INTERNAL_INVESTMENT_MEMO.md").write_text("\n".join(investment_memo) + "\n", encoding="utf-8")

    pilot_charter = [
        "# Pilot Charter",
        "",
        "## Goal",
        "Run a blinded validation of the selection layer on sponsor candidate pools and return a decision-grade readout.",
        "",
        "## Scope",
        "- Frozen scoring only.",
        "- Sponsor-provided blinded candidate pool.",
        "- Pre-specified top-N cutoffs.",
        "- Comparison against sponsor incumbent and any agreed baseline.",
        "",
        "## Exit criteria",
        "- Sponsor confirms signal or no-signal.",
        "- Sponsor confirms pilot continuation or archive.",
        "- Any next-step business discussion is anchored to the readout.",
        "",
    ]
    (OUT / "MODERNA_PILOT_CHARTER.md").write_text("\n".join(pilot_charter) + "\n", encoding="utf-8")

    readout_template = [
        "# Executive Readout Template",
        "",
        "| section | content |",
        "|---|---|",
        "| headline | one-sentence conclusion on top-N lift |",
        "| evidence | scoreboard, rank table, and audit trail summary |",
        "| caveats | known limits and claim boundary |",
        "| recommendation | pilot, term-sheet discussion, or archive |",
        "| next owner | named sponsor contact and internal owner |",
        "",
    ]
    (OUT / "MODERNA_EXEC_READOUT_TEMPLATE.md").write_text("\n".join(readout_template) + "\n", encoding="utf-8")

    launch_checklist = [
        "# Launch Checklist",
        "",
        "- Confirm NDA and approved scope.",
        "- Confirm sponsor technical owner and legal owner.",
        "- Confirm candidate pool format and slot budget.",
        "- Confirm lock date and readout date.",
        "- Confirm what can be shared pre- and post-lock.",
        "- Confirm the comparison baseline.",
        "- Confirm the final deliverable format.",
        "",
    ]
    (OUT / "MODERNA_LAUNCH_CHECKLIST.md").write_text("\n".join(launch_checklist) + "\n", encoding="utf-8")

    raci = [
        "# RACI",
        "",
        "| task | us | sponsor | legal |",
        "|---|---|---|---|",
        "| teaser | R | C | I |",
        "| NDA path | R | A | A |",
        "| data room | C | R | I |",
        "| score lock | R | C | I |",
        "| blinded validation | R | R | I |",
        "| readout | R | A | I |",
        "| term-sheet discussion | C | A | A |",
        "",
    ]
    (OUT / "MODERNA_RACI.md").write_text("\n".join(raci) + "\n", encoding="utf-8")

    board_packet = [
        "# Board Packet",
        "",
        "## Decision needed",
        "Approve the next-step package for sponsor outreach and NDA-stage validation.",
        "",
        "## What is being approved",
        "- External teaser use.",
        "- CDA/NDA path.",
        "- Blinded validation on sponsor candidate pools.",
        "- Pilot-readout and option/term-sheet discussion if signal is positive.",
        "",
        "## Why this matters",
        "- It preserves the claim boundary.",
        "- It gives the sponsor a decision-ready route.",
        "- It keeps the platform stable while testing the selection layer.",
        "",
    ]
    (OUT / "MODERNA_BOARD_PACKET.md").write_text("\n".join(board_packet) + "\n", encoding="utf-8")

    sow = [
        "# SOW Scaffold",
        "",
        "This is a placeholder for an NDA-stage statement of work.",
        "",
        "| section | placeholder |",
        "|---|---|",
        "| scope | frozen validation of the selection layer on blinded sponsor data |",
        "| tasks | ingestion, score lock, ranking, readout, and audit trail |",
        "| deliverables | ranked outputs, summary memo, top-N table, uncertainty notes |",
        "| assumptions | sponsor provides blinded candidate pool and slot policy |",
        "| exclusions | no platform replacement, no unblinded internals pre-clearance |",
        "",
    ]
    (OUT / "MODERNA_SOW_SCAFFOLD.md").write_text("\n".join(sow) + "\n", encoding="utf-8")

    compliance = [
        "# Compliance Checklist",
        "",
        "- NDA/CDA executed before sponsor data transfer.",
        "- Scope approved by named sponsor owner.",
        "- Data handling terms confirmed.",
        "- Score lock procedure documented.",
        "- Access logging enabled for exchanged files.",
        "- Claim boundary preserved in all external materials.",
        "- No unapproved feature list or training detail shared.",
        "",
    ]
    (OUT / "MODERNA_COMPLIANCE_CHECKLIST.md").write_text("\n".join(compliance) + "\n", encoding="utf-8")

    exec_decision = [
        "# Executive Decision Memo",
        "",
        "The external move is simple: send the teaser, route to the owner, secure NDA, and stage a blinded validation.",
        "",
        "Decision ask:",
        "- Approve pre-NDA outreach.",
        "- Approve NDA-stage validation.",
        "- Approve the board packet as the external approval set.",
        "",
        "Why now:",
        "- Local evidence already separates selection champion from strict reliability lane.",
        "- The sponsor workflow can absorb this without a platform rewrite.",
        "- The next decision is an operational one, not a scientific one.",
        "",
    ]
    (OUT / "MODERNA_EXEC_DECISION_MEMO.md").write_text("\n".join(exec_decision) + "\n", encoding="utf-8")

    ic_brief = [
        "# Investment Committee Brief",
        "",
        "## Decision request",
        "Approve a pre-NDA outreach package and an NDA-stage blinded validation path.",
        "",
        "## Why this is worth approving",
        "- Local evidence already shows a same-board champion and a separate strict reliability lane.",
        "- The sponsor is not asked to replace its platform.",
        "- The next step is a frozen validation, not a risky broad deployment.",
        "",
        "## Risks",
        "- Overclaiming beyond the data.",
        "- Sending unblinded internals too early.",
        "- Skipping the sponsor's legal and technical owner mapping.",
        "",
        "## Mitigation",
        "- Keep the teaser non-confidential.",
        "- Use NDA before sponsor data exchange.",
        "- Keep the claim boundary explicit.",
        "",
    ]
    (OUT / "MODERNA_INVESTMENT_COMMITTEE_BRIEF.md").write_text("\n".join(ic_brief) + "\n", encoding="utf-8")

    data_room_index = [
        "# Data Room Index",
        "",
        "## Sponsor-provided",
        "- Blinded candidate pool",
        "- Slot budget / selection policy",
        "- Incumbent baseline output, if available",
        "- Optional assay controls / wetlab annotations",
        "",
        "## Our-provided",
        "- Frozen scorer hash",
        "- Score-lock log",
        "- Top-N ranking table",
        "- Uncertainty report",
        "- Audit trail summary",
        "- Readout memo",
        "",
        "## Exclusions",
        "- Unblinded internals before lock",
        "- Feature list leakage",
        "- Raw training joins or private candidate tables",
        "",
    ]
    (OUT / "MODERNA_DATA_ROOM_INDEX.md").write_text("\n".join(data_room_index) + "\n", encoding="utf-8")

    platform_guardrails = [
        "# Platform Guardrails",
        "",
        "- The work is a selection layer, not a platform replacement.",
        "- Public claims stay within benchmark scope.",
        "- Sponsor data stays blinded until lock.",
        "- No feature-level disclosure before legal clearance.",
        "- No prospective efficacy claims in pre-NDA materials.",
        "- No unvetted endpoint merging across incompatible datasets.",
        "",
    ]
    (OUT / "MODERNA_PLATFORM_GUARDRAILS.md").write_text("\n".join(platform_guardrails) + "\n", encoding="utf-8")

    exec_qa = [
        "# Executive Q&A",
        "",
        "Q: Why now?",
        "A: Because the sponsor workflow already exists, and the selection layer can be validated without changing the platform.",
        "",
        "Q: Why us?",
        "A: The local board shows a same-board champion, a strict reliability lane, and a claim-safe path.",
        "",
        "Q: What is the ask?",
        "A: Pre-NDA routing, then NDA-backed blinded validation.",
        "",
        "Q: What is the risk?",
        "A: Overclaiming. That is controlled by the claim boundary and blinded protocol.",
        "",
    ]
    (OUT / "MODERNA_EXEC_QA.md").write_text("\n".join(exec_qa) + "\n", encoding="utf-8")

    stakeholder_map = [
        "# Stakeholder Map",
        "",
        "| stakeholder | role | ask |",
        "|---|---|---|",
        "| BD owner | route and sponsor | send teaser / NDA path |",
        "| technical owner | validation lead | approve blinded validation scope |",
        "| legal owner | confidentiality | clear NDA / CDA |",
        "| sponsor exec | budget / priority | approve pilot / archive |",
        "| data owner | data transfer | provide blinded candidate pool |",
        "",
    ]
    (OUT / "MODERNA_STAKEHOLDER_MAP.md").write_text("\n".join(stakeholder_map) + "\n", encoding="utf-8")

    signoff = [
        "# Sign-Off Memo",
        "",
        "Approve the external package for pre-NDA outreach and sponsor routing.",
        "",
        "This memo confirms:",
        "- The package is claim-bounded.",
        "- The package separates selection champion from reliability lane.",
        "- The package is ready for NDA-stage validation if sponsor interest is confirmed.",
        "",
    ]
    (OUT / "MODERNA_SIGNOFF_MEMO.md").write_text("\n".join(signoff) + "\n", encoding="utf-8")

    talking_points = [
        "# Talking Points",
        "",
        "- We are a selection layer for fixed-slot individualized neoantigen workflows.",
        "- Our strongest local result is same-board prioritization, not platform replacement.",
        "- We keep reliability, strict validation, and claim boundary separate.",
        "- The next step is an NDA-backed blinded validation on sponsor data.",
        "- If the sponsor has a candidate pool, integration friction stays low.",
        "",
    ]
    (OUT / "MODERNA_TALKING_POINTS.md").write_text("\n".join(talking_points) + "\n", encoding="utf-8")

    launch_deck = [
        "# Launch Deck Outline",
        "",
        "1. Problem framing: fixed-slot selection under individualized vaccine constraints.",
        "2. Local evidence: same-board champion and strict reliability lane.",
        "3. Claim boundary: what we say pre-NDA and what we do not say.",
        "4. Sponsor fit: blinded validation on existing candidate pools.",
        "5. Operating model: roles, data room, score lock, readout.",
        "6. Decision path: archive, pilot, or term-sheet discussion.",
        "",
    ]
    (OUT / "MODERNA_LAUNCH_DECK_OUTLINE.md").write_text("\n".join(launch_deck) + "\n", encoding="utf-8")

    board_summary = [
        "# Board Summary",
        "",
        "- Same-board champion: KG_GA_evolved.",
        "- Same-board AUPRC: 0.933.",
        "- Strict no-overlap winner: BAR-Neo_confidence.",
        "- Strict AUPRC: 0.857.",
        "- External ask: pre-NDA routing, NDA-backed blinded validation, pilot or archive.",
        "",
    ]
    (OUT / "MODERNA_BOARD_SUMMARY.md").write_text("\n".join(board_summary) + "\n", encoding="utf-8")

    deal_terms = [
        "# Deal Terms Snapshot",
        "",
        "| term | placeholder |",
        "|---|---|",
        "| scope | blinded validation of selection layer |",
        "| data | sponsor candidate pool, blinded until lock |",
        "| deliverables | score file, rank table, readout memo, audit trail |",
        "| decision | pilot extension, evaluation fee, option discussion, or archive |",
        "| exclusions | no platform replacement claim, no unblinded internals pre-clearance |",
        "",
    ]
    (OUT / "MODERNA_DEAL_TERMS_SNAPSHOT.md").write_text("\n".join(deal_terms) + "\n", encoding="utf-8")

    risk_mitigation = [
        "# Risk Mitigation One-Pager",
        "",
        "## Key risks",
        "- Overclaiming outside the benchmark scope.",
        "- Premature disclosure of internals.",
        "- Sponsor not having a clear owner or data path.",
        "",
        "## Mitigations",
        "- Use claim boundary language everywhere.",
        "- Keep all pre-NDA materials non-confidential.",
        "- Require named legal and technical owners before data transfer.",
        "- Use score lock and audit trail before readout.",
        "",
    ]
    (OUT / "MODERNA_RISK_MITIGATION_ONE_PAGER.md").write_text("\n".join(risk_mitigation) + "\n", encoding="utf-8")

    appendix_index = [
        "# Appendix Index",
        "",
        "A. Benchmark summary tables",
        "B. Strict reliability lane",
        "C. Public comparator context",
        "D. Latest literature watchlist",
        "E. Contact routing",
        "F. NDA / pilot operating pack",
        "G. Board and investment committee packets",
        "",
    ]
    (OUT / "MODERNA_APPENDIX_INDEX.md").write_text("\n".join(appendix_index) + "\n", encoding="utf-8")

    ceo_brief = [
        "# CEO Brief",
        "",
        "- Local evidence is strong enough to justify external routing.",
        "- The sponsor ask is narrow: blinded validation on existing candidate pools.",
        "- The claim boundary stays intact: selection layer, not platform replacement.",
        "- The next decision is a go/no-go on pre-NDA outreach and NDA staging.",
        "",
    ]
    (OUT / "MODERNA_CEO_BRIEF.md").write_text("\n".join(ceo_brief) + "\n", encoding="utf-8")

    deal_snapshot = [
        "# Deal Snapshot",
        "",
        "| item | state |",
        "|---|---|",
        "| local champion | KG_GA_evolved |",
        "| strict reliability lane | BAR-Neo_confidence |",
        "| sponsor ask | blinded validation |",
        "| business model | selection layer |",
        "| next step | NDA / pilot / archive |",
        "",
    ]
    (OUT / "MODERNA_DEAL_SNAPSHOT.md").write_text("\n".join(deal_snapshot) + "\n", encoding="utf-8")

    faq_index = [
        "# FAQ Index",
        "",
        "- Why now?",
        "- Why us?",
        "- Why sponsor?",
        "- What do we share pre-NDA?",
        "- What changes after NDA?",
        "- What is the pilot exit criterion?",
        "- What is out of scope?",
        "",
    ]
    (OUT / "MODERNA_FAQ_INDEX.md").write_text("\n".join(faq_index) + "\n", encoding="utf-8")

    signoff_tracker = [
        "# Sign-Off Tracker",
        "",
        "| step | owner | status |",
        "|---|---|---|",
        "| route teaser | us | ready |",
        "| approve external send | internal | pending |",
        "| confirm sponsor owner | sponsor | pending |",
        "| NDA route | legal | pending |",
        "| blinded validation | both | blocked until NDA |",
        "",
    ]
    (OUT / "MODERNA_SIGNOFF_TRACKER.md").write_text("\n".join(signoff_tracker) + "\n", encoding="utf-8")

    value_thesis = [
        "# Value Thesis",
        "",
        "| value lever | how CROSS-Neo helps | why partner cares |",
        "|---|---|---|",
        "| selection efficiency | better ranking under fixed slot budget | fewer low-value candidates enter manufacturing |",
        "| false-positive reduction | claim-safe and strict lanes separate risky outputs | less wasted assay and synthesis cost |",
        "| auditability | reason codes and uncertainty flags | easier internal review and governance |",
        "| wetlab focus | control annotations and top-N exports | better experiment design and cleaner interpretation |",
        "| validation speed | blinded validation protocol | faster decision on adoption or optioning |",
        "| compatibility | selection layer, not platform replacement | lower integration friction |",
        "",
    ]
    (OUT / "MODERNA_VALUE_THESIS.md").write_text("\n".join(value_thesis) + "\n", encoding="utf-8")

    kpi = [
        "# KPI Scorecard",
        "",
        "| KPI | current local result | partner relevance |",
        "|---|---|---|",
        "| same-board AUPRC | 0.933 | local prioritization win |",
        "| frozen validation-like AUPRC | 0.978 | robustness on held-in validation-like slices |",
        "| low/medium leakage AUPRC | 0.757 | stress-tested prioritization |",
        "| strict no-overlap AUPRC | 0.857 winner is BAR-Neo_confidence | relevance for reliability lane |",
        "| public comparator context | RF_biophys / MHCflurry / PRIME / TransPHLA all tracked | market/context awareness |",
        "| latest literature watchlist | NetMHCpan-4.2 / ImmugenX / PepFore / NeoGuider / TrambaHLApan / NeoTImmuML | shows current landscape |",
        "",
    ]
    (OUT / "MODERNA_KPI_SCORECARD.md").write_text("\n".join(kpi) + "\n", encoding="utf-8")

    decision_tree = [
        "# Partner Decision Tree",
        "",
        "1. Do they own individualized mRNA vaccine selection or related BD?",
        "   - yes -> send teaser and ask for NDA path",
        "   - no -> route to owner or archive",
        "",
        "2. Do they have a blinded candidate pool and a fixed slot budget?",
        "   - yes -> propose blinded validation",
        "   - no -> propose pilot / future collaboration",
        "",
        "3. Do they want evaluation-only or future license potential?",
        "   - evaluation-only -> NDA and benchmark",
        "   - license interest -> move to term sheet after validation signal",
        "",
        "4. Do they require wetlab controls and audit trail?",
        "   - yes -> share post-NDA protocol and control mapping",
        "   - no -> keep to score and ranking outputs",
        "",
    ]
    (OUT / "MODERNA_PARTNER_DECISION_TREE.md").write_text("\n".join(decision_tree) + "\n", encoding="utf-8")

    same.to_csv(OUT / "partner_same_board_metrics.tsv", sep="\t", index=False)
    strict.to_csv(OUT / "partner_strict_reliability_metrics.tsv", sep="\t", index=False)
    public.to_csv(OUT / "partner_public_context.tsv", sep="\t", index=False)
    contact.to_csv(OUT / "partner_contact_routes.tsv", sep="\t", index=False)
    exec_brief.to_csv(OUT / "partner_exec_brief.tsv", sep="\t", index=False)
    send_packet.to_csv(OUT / "partner_send_packet.tsv", sep="\t", index=False)
    claim.to_csv(OUT / "partner_claim_boundary.tsv", sep="\t", index=False)

    build_figures(same, strict, public)
    PAGE.write_text(build_html(same, strict, public, contact, exec_brief, send_packet, claim), encoding="utf-8")

    for src in FIG.glob("*"):
        if src.is_file():
            safe_copy(src, ASSET / src.name)
            safe_copy(src, LIVE_ASSET / src.name)
    for src in [OUT / "partner_same_board_metrics.tsv", OUT / "partner_strict_reliability_metrics.tsv", OUT / "partner_public_context.tsv", OUT / "partner_contact_routes.tsv", OUT / "partner_exec_brief.tsv", OUT / "partner_send_packet.tsv", OUT / "partner_claim_boundary.tsv"]:
        safe_copy(src, ASSET / src.name)
        safe_copy(src, LIVE_ASSET / src.name)
    safe_copy(PAGE, LIVE_PAGE)

    pack = OUT / "partner_send_packet"
    if pack.exists():
        shutil.rmtree(pack)
    pack.mkdir(parents=True, exist_ok=True)
    for name in [
        "MODERNA_NONCONFIDENTIAL_TEASER.md",
        "MODERNA_OUTREACH_EMAIL_AFTER_IP.md",
        "MODERNA_CALL_SCRIPT.md",
        "MODERNA_FIRST_MEETING_AGENDA.md",
        "MODERNA_DD_QUESTIONS.md",
        "MODERNA_FOLLOWUP_EMAIL.md",
        "MODERNA_FAQ_OBJECTION_SHEET.md",
        "MODERNA_ONE_PAGE_PITCH.md",
        "MODERNA_NDA_CHECKLIST.md",
        "MODERNA_MEETING_NOTES_TEMPLATE.md",
        "MODERNA_RISK_REGISTER.md",
        "MODERNA_CONTACT_TRACKER.md",
        "MODERNA_REPLY_HANDLING_MATRIX.md",
        "MODERNA_INTERNAL_HANDOFF_CHECKLIST.md",
        "MODERNA_NDA_ASK_TEMPLATE.md",
        "MODERNA_BLINDED_DATA_REQUEST.md",
        "MODERNA_SCORE_LOCK_PROCEDURE.md",
        "MODERNA_POST_NDA_PROTOCOL.md",
        "MODERNA_SLIDE_OUTLINE.md",
        "MODERNA_INTERNAL_REVIEW_MEMO.md",
        "MODERNA_GO_NO_GO_CRITERIA.md",
        "MODERNA_ANSWER_BANK.md",
        "MODERNA_DEAL_MEMO.md",
        "MODERNA_IMPACT_MEMO.md",
        "MODERNA_BOARD_MEMO.md",
        "MODERNA_COMPETITIVE_POSITIONING.md",
        "MODERNA_NEXT_STEP_CHECKLIST.md",
        "MODERNA_WHY_NOW.md",
        "MODERNA_WHY_SPONSOR.md",
        "MODERNA_WHY_TEAM.md",
        "MODERNA_ONE_SLIDE_EXECUTIVE.md",
        "MODERNA_KILLER_SUMMARY.md",
        "MODERNA_DECISION_SNAPSHOT.md",
        "MODERNA_COVER_LETTER.md",
        "MODERNA_PARTNER_INVITATION.md",
        "MODERNA_BOARD_ASK.md",
        "MODERNA_PILOT_PROPOSAL.md",
        "MODERNA_SUCCESS_METRICS.md",
        "MODERNA_MILESTONES.md",
        "MODERNA_BUDGET_SKETCH.md",
        "MODERNA_MUTUAL_ACTION_PLAN.md",
        "MODERNA_DATA_ROOM_CHECKLIST.md",
        "MODERNA_BUSINESS_CASE.md",
        "MODERNA_APPROVAL_PATH.md",
        "MODERNA_INTEGRATION_CHECKLIST.md",
        "MODERNA_TERM_SHEET_SCAFFOLD.md",
        "MODERNA_ROI_MODEL.md",
        "MODERNA_SPONSOR_OPERATING_MODEL.md",
        "MODERNA_INTERNAL_INVESTMENT_MEMO.md",
        "MODERNA_PILOT_CHARTER.md",
        "MODERNA_EXEC_READOUT_TEMPLATE.md",
        "MODERNA_LAUNCH_CHECKLIST.md",
        "MODERNA_RACI.md",
        "MODERNA_BOARD_PACKET.md",
        "MODERNA_SOW_SCAFFOLD.md",
        "MODERNA_COMPLIANCE_CHECKLIST.md",
        "MODERNA_EXEC_DECISION_MEMO.md",
        "MODERNA_INVESTMENT_COMMITTEE_BRIEF.md",
        "MODERNA_DATA_ROOM_INDEX.md",
        "MODERNA_PLATFORM_GUARDRAILS.md",
        "MODERNA_EXEC_QA.md",
        "MODERNA_STAKEHOLDER_MAP.md",
        "MODERNA_SIGNOFF_MEMO.md",
        "MODERNA_TALKING_POINTS.md",
        "MODERNA_LAUNCH_DECK_OUTLINE.md",
        "MODERNA_BOARD_SUMMARY.md",
        "MODERNA_DEAL_TERMS_SNAPSHOT.md",
        "MODERNA_RISK_MITIGATION_ONE_PAGER.md",
        "MODERNA_APPENDIX_INDEX.md",
        "MODERNA_CEO_BRIEF.md",
        "MODERNA_DEAL_SNAPSHOT.md",
        "MODERNA_FAQ_INDEX.md",
        "MODERNA_SIGNOFF_TRACKER.md",
        "MODERNA_VALUE_THESIS.md",
        "MODERNA_KPI_SCORECARD.md",
        "MODERNA_PARTNER_DECISION_TREE.md",
        "MODERNA_EXECUTIVE_COVER.md",
        "MODERNA_VALUE_CALCULATION.md",
        "MODERNA_FAST_PATH_PLAN.md",
        "PRE_NDA_SEND_PACKET_README.md",
        "moderna_contact_routes.tsv",
    ]:
        src = MOD / name if (MOD / name).exists() else OUT / name
        shutil.copy2(src, pack / name)
    shutil.copy2(ROOT / "project/results/cross_neo_exec_brief_2026_05_11/CROSS_NEO_EXEC_BRIEF_KR.md", pack / "CROSS_NEO_EXEC_BRIEF_KR.md")
    shutil.copy2(SOTA / "latest_literature_sota_watchlist.tsv", pack / "latest_literature_sota_watchlist.tsv")
    zip_base = OUT / "partner_send_packet"
    if (OUT / "partner_send_packet.zip").exists():
        (OUT / "partner_send_packet.zip").unlink()
    shutil.make_archive(str(zip_base), "zip", pack)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": str(OUT),
        "html_path": str(PAGE),
        "live_html_path": str(LIVE_PAGE),
        "same_board_champion": "KG_GA_evolved",
        "strict_winner": "BAR-Neo_confidence",
        "pre_nda_zip": str(OUT / "partner_send_packet.zip"),
    }
    (OUT / "partner_pack_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
