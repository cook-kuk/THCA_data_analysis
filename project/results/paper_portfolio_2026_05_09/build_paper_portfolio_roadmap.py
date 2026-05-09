#!/usr/bin/env python3
"""Build a portfolio-level paper decision roadmap from current local state.

This is deliberately a decision scaffold, not manuscript prose. It avoids
voice-protected sections and fixes the claim boundary for each candidate.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "project/results/paper_portfolio_2026_05_09"
HUB = ROOT / "project/papers_hub_2026_05_04"


ROWS = [
    {
        "rank": 1,
        "paper_id": "Paper 1",
        "track": "DM1 molecular dark matter",
        "decision_band": "SUBMIT_NOW",
        "current_grade": "A",
        "venue_read": "Realistic primary: Cell Reports Medicine / JCI Insight. Stretch: Nature Communications only with stronger extension or prospective cohort. Fallback: npj Precision Oncology.",
        "one_line": "8-gene DM1/DM2 molecular dark-matter subtype with fusion enrichment, epigenetic thyroid-lineage silencing, clinical/meta support, and reviewer defense through v18.",
        "why_it_can_publish": "Most complete package: manuscript v8 scaffold, Fig 1-8 map, Q1-Q14 reviewer surface, v13/v14 mechanism locks, deconvolution rollup, methods dossier.",
        "main_blocker": "Author-only voice-protected manuscript sections and no wet-lab causality. Bundang/prospective cohort would raise ceiling but is not required for fallback.",
        "next_action": "Finish author-voice sections, freeze main claims, submit preprint/manuscript before opening new analyses.",
        "do_not_claim": "Do not claim causal methylation mechanism or treatment response restoration.",
        "source_paths": "project/manuscript_v8/02_outline.md; project/HANDOFF_2026_05_09.md; project/papers_hub_2026_05_04/paper1_reviewer_defense_dashboard.html",
    },
    {
        "rank": 2,
        "paper_id": "Paper 1B",
        "track": "BRAF-cPTC DM2 / DM1xTERT escape",
        "decision_band": "MERGE_OR_STANDALONE_DECISION",
        "current_grade": "A-",
        "venue_read": "Nature Communications stretch if separated and framed cleanly; otherwise strengthens Paper 1. Fallback: CRM / JCI Insight.",
        "one_line": "In BRAF-driven cPTC, DM2 appears to be the adverse RAI/progression arm, while DM1 is usually immune-rich/protective except for a DM1 x TERT aggressive escape state.",
        "why_it_can_publish": "Strong survival and external Landa replication: DM1 x TERT pooled HR about 19; BRAF-like DM2 PFI HR about 4.7; RAI-refractory signal directionally coherent.",
        "main_blocker": "Potential cannibalization/conflict with Paper 1; Lee 2024 lacks survival/driver/TERT fields; image rescue is not load-bearing.",
        "next_action": "Make a binary decision: fold into Paper 1 as high-impact extension, or hold for standalone after Paper 1 locks.",
        "do_not_claim": "Do not say generic DM1 is aggressive in BRAF-cPTC; the escape state is DM1 x TERT.",
        "source_paths": "project/results/p2_braf_nature_sprint_2026_05_09/_synthesis/NATURE_TIER_SYNTHESIS_V2_CODEX_TAKEOVER.md; project/papers_hub_2026_05_04/paper1_paper2_braf_axis_takeover_dossier.html",
    },
    {
        "rank": 3,
        "paper_id": "Paper 11",
        "track": "Pan-cancer portable DM1 axis",
        "decision_band": "NEXT_MAJOR_AFTER_P1",
        "current_grade": "A-",
        "venue_read": "Nature Communications primary candidate. Fallback: Genome Medicine / npj Precision Oncology.",
        "one_line": "Thyroid-derived DM1/dedifferentiation architecture ports across 32 TCGA lineages, is prognostic in multiple cancers, and links to MYC/NAMPT dependency plus MAPK-inhibitor PRISM sensitivity.",
        "why_it_can_publish": "Large scale: 10,978-11,069 TCGA samples, 32-33 lineages, 10 prognostic lineages, DepMap/PRISM actionability, ICI reserve, epigenetic proxy.",
        "main_blocker": "Manuscript not started; functional validation absent; full pan-cancer methylation fetch blocked and replaced by RNA epigenetic-machinery proxy.",
        "next_action": "Do not start during Paper 1 voice marathon unless making only a handoff/dossier. After Paper 1 submission, build Nat Commun 6-figure story.",
        "do_not_claim": "Do not overstate full causal methylation without HM450/wet-lab layer.",
        "source_paths": "project/results/paper11_pancancer/PAPER11_OUTLINE.md; project/papers_hub_2026_05_04/paper11_pancancer.html",
    },
    {
        "rank": 4,
        "paper_id": "Paper 4",
        "track": "Korean GD / Pan-Asian HLA",
        "decision_band": "UPGRADE_IF_DATA_ARRIVES",
        "current_grade": "B+ now; A if de novo cohort arrives",
        "venue_read": "Nature Communications only with Korean adult GD germline NGS case-control cohort. Fallback now: Genes and Immunity / HLA / Scientific Reports.",
        "one_line": "Pan-Asian AITD/GD HLA architecture centered on DPB1*05:01, B*46:01, A*02:07, and C*01:02, with specificity and ancestry/resource layers.",
        "why_it_can_publish": "Meta/resource story is coherent; DPB1*05:01 pooled OR about 2.1; B*46:01 and A*02:07 signals; C*01:02 GD-specific focus.",
        "main_blocker": "No de novo Korean adult GD NGS cases + matched controls; without that, this is synthesis/resource rather than original Korean association.",
        "next_action": "Keep as high-upside backlog. Prioritize cohort acquisition over more computation.",
        "do_not_claim": "Do not call current data a definitive Korean GD association discovery.",
        "source_paths": "project/results/hla_two_paper_synthesis_2026_05_09/HLA_TWO_PAPER_NCOMM_SYNTHESIS.md; project/papers_hub_2026_05_04/hla_ncomm_upgrade_2026_05_09.html",
    },
    {
        "rank": 5,
        "paper_id": "Paper 2",
        "track": "H&E to DM1 pathology projection",
        "decision_band": "VIABLE_PILOT",
        "current_grade": "B+",
        "venue_read": "npj Digital Medicine / JCI Insight. CRM possible only if framed modestly and external H&E improves. Fallback: Paper 1 supplement.",
        "one_line": "Foundation-model pathology classifier recovers part of the DM1 axis from H&E, with N=59 final AUC around 0.75 and UNI LOTO about 0.85, but BRAF/cPTC image rescue fails.",
        "why_it_can_publish": "Architecture-bound closure was overturned; CLAM/ViT-L/UNI signal exists; audit dossier and confound case study make the limitations usable.",
        "main_blocker": "Small N, calibration issues, TSS/center confounding, no external Korean H&E validation, BRAF/cPTC generic image rescue negative.",
        "next_action": "Keep as pilot/methods paper; do not let it block Paper 1. External H&E or K2 slides are the only real upgrade.",
        "do_not_claim": "Do not claim clinical deployment or robust BRAF/cPTC image biology.",
        "source_paths": "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/PAPER2_OUTLINE_2026_05_08.md; project/papers_hub_2026_05_04/paper2_image_dm1_audit_dossier.html",
    },
    {
        "rank": 6,
        "paper_id": "Paper 2B",
        "track": "HT-overlap PTC immune/HLA context",
        "decision_band": "BIOLOGY_STRONG_ALLELE_WEAK",
        "current_grade": "B+ biology; C allele",
        "venue_read": "JCI Insight / Frontiers in Immunology / Cancer Immunology Research style if kept expression/spatial. Allele association not ready.",
        "one_line": "HT-overlap/immune-active PTC has HLA-II, IFN-gamma, TLS, APC, and spatial niche biology; allele-level DRB1*04:05 remains exploratory.",
        "why_it_can_publish": "Expression-module and spatial biology are strong: HLA-II/IFNG/TLS/APC layers replicate; negative controls against HLA-G and somatic HLA loss are useful.",
        "main_blocker": "Explicit HT-PTC allele test is n=9 vs n=9; FDR-negative. Needs germline NGS HLA cohort with HT labels.",
        "next_action": "Use as Paper 2 biology reserve or separate immune-context paper only after Paper 1/Paper 2 scope is locked.",
        "do_not_claim": "Do not claim HLA allele risk/outcome prediction from cancer cohorts.",
        "source_paths": "project/results/hla_two_paper_synthesis_2026_05_09/HLA_TWO_PAPER_NCOMM_SYNTHESIS.md; project/papers_hub_2026_05_04/paper2_hla.html",
    },
    {
        "rank": 7,
        "paper_id": "Paper 3",
        "track": "ICI vulnerability in molecularly dark thyroid cancer",
        "decision_band": "FROZEN_DESIGN",
        "current_grade": "B design; not executable now",
        "venue_read": "JCI Insight / CRM possible after Track B; not currently submission-ready.",
        "one_line": "HLA loss, neoantigen architecture, and immune ecotypes define ICI vulnerability/readiness, not thyroid ICI response prediction.",
        "why_it_can_publish": "Track A registry is detailed; Track B-lite and Paper 11 ICI reserve show directionally useful immune-response transfer signals.",
        "main_blocker": "Track B explicitly blocked until Paper 1 bioRxiv + Paper 2 closures + explicit author command; no robust thyroid ICI-treated cohort.",
        "next_action": "Keep frozen. Use only as future plan or Paper 11 reserve.",
        "do_not_claim": "Do not claim thyroid ICI response prediction.",
        "source_paths": "project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md; project/papers_hub_2026_05_04/paper3_ici_status_dossier_2026_05_09.html",
    },
    {
        "rank": 8,
        "paper_id": "Paper 9",
        "track": "Synthetic lethality / perturbation vulnerability",
        "decision_band": "WETLAB_REQUIRED",
        "current_grade": "C+ now; A only with experiments",
        "venue_read": "After wet-lab: JCI Insight / Cancer Research / CRM. Now: internal dossier only.",
        "one_line": "Public DepMap/PRISM and Paper9 metabolic screens prioritize vulnerability candidates: MAPK-inhibitor sensitivity, MYC/NAMPT dependency, and SLC1A5 > GLUD1 > GLS wet-lab queue.",
        "why_it_can_publish": "Strong hypothesis generator: PRISM FDR hits 11 with 7 canonical MAPK; AZD-0364 FDR 5.88e-7; MYC/NAMPT dependency; glutamine-axis prioritization.",
        "main_blocker": "Not validated synthetic lethality. Strict thyroid model coverage underpowered; drug-response claims need matched models and rescue.",
        "next_action": "Keep as actionability box/reviewer reserve until wet-lab perturbation exists.",
        "do_not_claim": "Do not say validated synthetic lethality, clinical recommendation, or patient-selection biomarker.",
        "source_paths": "project/results/p_synthetic_lethality_2026_05_09/SYNTHETIC_LETHALITY_BRIEF_KR.md; project/papers_hub_2026_05_04/paper1_synthetic_lethality_dossier.html",
    },
    {
        "rank": 9,
        "paper_id": "Paper 10/12",
        "track": "Atlas / network architecture companion",
        "decision_band": "ABSORB_INTO_P11",
        "current_grade": "B- standalone; useful support",
        "venue_read": "Best as Paper 11 figures/supplement. Standalone only as resource/methods note if expanded.",
        "one_line": "Cell-line atlas and 21-gene network/module architecture explain the pan-cancer DM1 axis and candidate modules.",
        "why_it_can_publish": "Useful figure inventory and module network exist, but the standalone question is weaker than Paper 11.",
        "main_blocker": "Not enough independent biological claim beyond Paper 11; would fragment the portfolio.",
        "next_action": "Use as Paper 11 support, not separate manuscript during marathon.",
        "do_not_claim": "Do not split it into a separate paper unless Paper 11 scope demands a companion resource.",
        "source_paths": "project/results/paper10_atlas/; project/results/paper12_network/summary.json",
    },
    {
        "rank": 10,
        "paper_id": "Neo/CROSS-Neo",
        "track": "Neoantigen prediction / vaccine platform",
        "decision_band": "SEPARATE_UNIVERSE",
        "current_grade": "B; promising but not thyroid paper",
        "venue_read": "Bioinformatics / Genome Medicine / methods-resource venue. Nat Biomed Eng only after much stronger external validation.",
        "one_line": "Post-TESLA neoantigen ranking/abstention framework with source-aware audits and vaccine triage dossiers.",
        "why_it_can_publish": "Large curation/platform surface and useful decision dossiers exist; source leakage controls are being built.",
        "main_blocker": "Source-heldout performance remains weak/heterogeneous; thyroid vaccine translation is speculative; could distract from Paper 1.",
        "next_action": "Keep completely separate from THCA Paper 1-4 unless only citing as future platform.",
        "do_not_claim": "Do not imply thyroid vaccine readiness from PAAD/PDAC vaccine evidence.",
        "source_paths": "project/results/cross_neo_v1_lockdown/locked_anchor_summary.md; project/results/p_cancer_vaccine_strategy_2026_05_09/SUMMARY.md",
    },
    {
        "rank": 11,
        "paper_id": "Spatial CNV / PantheonOS side tracks",
        "track": "Spatial ecosystem / CNV / TLS tooling",
        "decision_band": "PARK",
        "current_grade": "C+ as paper; useful methods/demo",
        "venue_read": "Methods/demo or supplementary support only unless a direct biological claim is hardened.",
        "one_line": "Spatial and agentic tooling generated useful caveats and niche/TLS observations, but most are not Paper 1 blocking.",
        "why_it_can_publish": "Some pages and figures are polished; v16-v18 spatial autopsies are strong reviewer caveats.",
        "main_blocker": "Core positive biological claim is not strong enough for a standalone paper; spatial MAPK-panel anti-correlation failed.",
        "next_action": "Use as caveat/support pages; do not spend manuscript time now.",
        "do_not_claim": "Do not convert failed spatial MAPK-panel result into positive mechanism evidence.",
        "source_paths": "project/papers_hub_2026_05_04/paper1_spatial_signal_decomposition_v17.html; project/papers_hub_2026_05_04/paper1_spatial_lag_pockets_v18.html",
    },
]


def write_tsv(path: Path) -> None:
    fields = list(ROWS[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(ROWS)


def md_table(rows: list[dict[str, object]]) -> str:
    fields = ["rank", "paper_id", "track", "decision_band", "current_grade", "venue_read", "next_action"]
    out = ["| " + " | ".join(fields) + " |", "| " + " | ".join(["---"] * len(fields)) + " |"]
    for row in rows:
        vals = [str(row[f]).replace("|", "/") for f in fields]
        out.append("| " + " | ".join(vals) + " |")
    return "\n".join(out)


def write_md(path: Path) -> None:
    submit_now = [r for r in ROWS if r["decision_band"] == "SUBMIT_NOW"]
    next_major = [r for r in ROWS if r["decision_band"] in {"MERGE_OR_STANDALONE_DECISION", "NEXT_MAJOR_AFTER_P1"}]
    gated = [r for r in ROWS if r["decision_band"] in {"UPGRADE_IF_DATA_ARRIVES", "WETLAB_REQUIRED", "FROZEN_DESIGN"}]
    park = [r for r in ROWS if r["decision_band"] in {"ABSORB_INTO_P11", "SEPARATE_UNIVERSE", "PARK"}]

    text = f"""# Paper Portfolio Roadmap — 2026-05-09

## Bottom line

지금 전체 포트폴리오의 정답은 **Paper 1을 먼저 ship**하는 것입니다. 그 다음 높은 upside는 **BRAF-cPTC DM2 / DM1xTERT axis**와 **Paper 11 pan-cancer portable DM1**입니다. Paper 4는 Korean adult GD NGS가 오면 바로 급상승하고, Paper 9는 wet-lab 전까지 synthetic lethality가 아니라 perturbation-priority dossier입니다.

## Immediate order

1. **Paper 1:** submit-now line. 새 분석보다 author-voice manuscript completion이 병목입니다.
2. **Paper 1B / BRAF axis:** Paper 1에 흡수할지 standalone으로 잠글지 결정해야 합니다.
3. **Paper 11:** Paper 1 이후 Nature Communications 후보로 크게 갑니다.
4. **Paper 2 image:** pilot paper로 살리되 Paper 1을 막지 않습니다.
5. **Paper 4 / Paper 9:** 각각 cohort / wet-lab data가 오기 전까지는 대기입니다.

## Submit / Upgrade / Park

### Submit now

{md_table(submit_now)}

### High-upside next

{md_table(next_major)}

### Gated but valuable

{md_table(gated)}

### Park / absorb

{md_table(park)}

## Full matrix

{md_table(ROWS)}

## Hard boundaries

- Voice-protected manuscript sections remain author-keyboard only.
- Paper 9 must not be called validated synthetic lethality without wet-lab perturbation and rescue.
- Paper 2 HT/HLA must not become an allele-association paper from n=9 vs n=9.
- Paper 3 must not claim thyroid ICI response prediction without thyroid ICI-treated data.
- Paper 10/12 should support Paper 11 unless a stronger standalone resource claim is created.
"""
    path.write_text(text, encoding="utf-8")


def esc(s: object) -> str:
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def row_html(row: dict[str, object]) -> str:
    cols = ["rank", "paper_id", "track", "decision_band", "current_grade", "venue_read", "one_line", "main_blocker", "next_action", "do_not_claim"]
    return "<tr>" + "".join(f"<td>{esc(row[c])}</td>" for c in cols) + "</tr>"


def write_html(path: Path) -> None:
    cards = [
        ("Submit now", "Paper 1", "Finish author-voice manuscript; stop opening new analyses."),
        ("Highest upside", "Paper 11", "Pan-cancer DM1 has Nat Commun shape after Paper 1."),
        ("Decision fork", "BRAF axis", "Merge into Paper 1 or hold as standalone."),
        ("Conditional leap", "Paper 4", "Korean GD NGS turns synthesis into original discovery."),
        ("Do not overclaim", "Paper 9", "Perturbation priority, not validated synthetic lethality."),
    ]
    card_html = "\n".join(
        f"<div class=\"card\"><span>{esc(k)}</span><b>{esc(v)}</b><p>{esc(t)}</p></div>" for k, v, t in cards
    )
    table = "\n".join(row_html(r) for r in ROWS)
    source_cards = "\n".join(
        f"<li><strong>{esc(r['paper_id'])}</strong>: <code>{esc(r['source_paths'])}</code></li>" for r in ROWS
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Paper Portfolio Roadmap — 2026-05-09</title>
  <style>
    :root {{
      --bg:#101217; --panel:#181d25; --panel2:#202735; --ink:#eee8dc; --muted:#aeb7c5;
      --line:#343d4c; --gold:#d9ad5f; --red:#f07167; --green:#8fd694; --blue:#8bbcff;
      --mono:"JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      --serif:"Cormorant Garamond", Georgia, serif; --sans:Inter, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    * {{ box-sizing:border-box; }}
    body {{ margin:0; background:var(--bg); color:var(--ink); font-family:var(--sans); line-height:1.5; }}
    .hero {{ padding:56px max(26px,6vw) 30px; border-bottom:1px solid var(--line); background:#121720; }}
    .kicker {{ color:var(--gold); font:700 12px var(--mono); text-transform:uppercase; }}
    h1 {{ font-family:var(--serif); font-size:clamp(44px,7vw,88px); line-height:.95; margin:12px 0 12px; letter-spacing:0; }}
    .lead {{ max-width:1050px; color:#d8dde8; font-size:19px; }}
    .cards {{ display:grid; grid-template-columns:repeat(5,minmax(0,1fr)); gap:12px; margin-top:24px; }}
    .card {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:15px; min-height:135px; }}
    .card span {{ display:block; color:var(--gold); font:700 11px var(--mono); text-transform:uppercase; }}
    .card b {{ display:block; font-size:22px; margin:8px 0; }}
    .card p {{ margin:0; color:var(--muted); font-size:13px; }}
    main {{ max-width:1560px; margin:0 auto; padding:28px 24px 70px; }}
    h2 {{ font-family:var(--serif); font-size:36px; margin:34px 0 12px; letter-spacing:0; }}
    .tl {{ background:var(--panel); border:1px solid var(--line); border-left:5px solid var(--gold); border-radius:8px; padding:16px 18px; font-size:18px; }}
    .table-wrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:8px; }}
    table {{ width:100%; border-collapse:collapse; min-width:1500px; background:var(--panel); font-size:13px; }}
    th, td {{ padding:10px 11px; border-bottom:1px solid var(--line); vertical-align:top; }}
    th {{ color:var(--gold); font:700 11px var(--mono); text-transform:uppercase; background:#151922; text-align:left; }}
    td:nth-child(4) {{ color:#f1d392; font-weight:700; }}
    td:nth-child(5) {{ color:#b8e6be; font-weight:700; }}
    code {{ font-family:var(--mono); color:#f1d392; font-size:.9em; }}
    .cols {{ display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:14px; }}
    .box {{ background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:16px; }}
    .box b {{ color:var(--gold); }}
    li {{ margin:8px 0; color:var(--muted); }}
    @media(max-width:1000px) {{ .cards,.cols {{ grid-template-columns:1fr; }} table {{ min-width:1200px; }} }}
  </style>
</head>
<body>
  <header class="hero">
    <div class="kicker">Portfolio decision map · 2026-05-09 · factual scaffold only</div>
    <h1>Where Each Paper Goes</h1>
    <p class="lead">A compact ranking of every active manuscript possibility: what can ship now, what needs one decisive dataset, what should be absorbed, and what must not be overclaimed.</p>
    <div class="cards">{card_html}</div>
  </header>
  <main>
    <section>
      <h2>Bottom Line</h2>
      <div class="tl"><strong>Decision:</strong> Paper 1 ships first. The next real upside is BRAF-axis triage and Paper 11 pan-cancer DM1. Paper 4 waits for Korean GD NGS. Paper 9 waits for wet-lab. Paper 2 image survives as a pilot, not the load-bearing main story.</div>
    </section>
    <section>
      <h2>Decision Matrix</h2>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Rank</th><th>Paper</th><th>Track</th><th>Band</th><th>Grade</th><th>Venue read</th><th>One-line</th><th>Main blocker</th><th>Next action</th><th>Do not claim</th></tr></thead>
          <tbody>{table}</tbody>
        </table>
      </div>
    </section>
    <section>
      <h2>Hard Boundaries</h2>
      <div class="cols">
        <div class="box"><b>Voice</b><p>Hook, Aim, Discussion mechanism, Limitations, Cover para 1, and Q9 remain author-keyboard only.</p></div>
        <div class="box"><b>Overclaim</b><p>Synthetic lethality, HLA allele association, and thyroid ICI response prediction are not validated by current data.</p></div>
        <div class="box"><b>Marathon</b><p>New work should be Paper-1-blocking only unless explicitly opened by the author.</p></div>
      </div>
    </section>
    <section>
      <h2>Source Paths</h2>
      <ul>{source_cards}</ul>
    </section>
  </main>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "paper_portfolio_decision_matrix.tsv")
    (OUT / "paper_portfolio_decision_matrix.json").write_text(json.dumps(ROWS, indent=2, ensure_ascii=False), encoding="utf-8")
    write_md(OUT / "PAPER_PORTFOLIO_ROADMAP_KR.md")
    write_html(HUB / "paper_portfolio_roadmap_2026_05_09.html")
    print("wrote", OUT / "paper_portfolio_decision_matrix.tsv")
    print("wrote", OUT / "PAPER_PORTFOLIO_ROADMAP_KR.md")
    print("wrote", HUB / "paper_portfolio_roadmap_2026_05_09.html")


if __name__ == "__main__":
    main()
