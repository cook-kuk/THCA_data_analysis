#!/usr/bin/env python3
"""THYRAI rebrand orchestrator.

Produces the THYRAI landing page, pipeline page, 4 platform pages,
publications / press / team pages, and accompanying assets. Idempotent.

Visual system inspired by modern biotech platform conventions.
Not affiliated with any drug-discovery AI company.
"""
from __future__ import annotations

import argparse
import html
import json
import logging
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable

# ------------- constants / paths -------------
PROJECT_ROOT = Path("/opt/thyroid-dash/project")
HTML_ROOT = PROJECT_ROOT / "reports/html"
PAGES = HTML_ROOT / "pages"
ASSETS = HTML_ROOT / "assets"
LOG_DIR = PROJECT_ROOT / "logs"
LOG_FILE = LOG_DIR / "rebrand_insilico.log"

YEAR = datetime.now().year

LEGAL_FOOTER = (
    "Visual system inspired by modern biotech platform conventions. "
    "Not affiliated with any drug-discovery AI company."
)

PLATFORM_TILES = [
    {
        "code": "OMEGA-THCA",
        "title": "Multi-omics biomarker discovery",
        "desc": "Transcriptomic, methylation, and immune signals from 1,509 TCGA and external GEO samples, reduced to 2,773 replicated biomarkers and eight druggable targets.",
        "href": "/reports/html/pages/platform_omega.html",
    },
    {
        "code": "CHEM-THCA",
        "title": "Chemistry triage and compound ranking",
        "desc": "ChEMBL-anchored compound sets for each target, scored by pChEMBL, mechanism-of-action clarity, and literature anchoring. Outputs are hit lists, never approvals.",
        "href": "/reports/html/pages/platform_chem.html",
    },
    {
        "code": "PANEL-THCA",
        "title": "Clinical-grade 54-gene panel benchmark",
        "desc": "A compact 54-gene LogReg panel with AUC 1.00 TCGA internal and 0.97–0.89 on four external GEO cohorts, benchmarked against the 112-gene reference.",
        "href": "/reports/html/pages/platform_panel.html",
    },
    {
        "code": "CLINIC-THCA",
        "title": "Bethesda-indeterminate triage simulator",
        "desc": "A simulation layer that projects panel performance into Bethesda III/IV cytology workflows — strictly for computational triage and research review.",
        "href": "/reports/html/pages/platform_clinic.html",
    },
]

PLATFORM_DETAIL = {
    "omega": {
        "index": 1,
        "code": "OMEGA-THCA",
        "title": "Discovery layer across multi-omics evidence",
        "summary": [
            "OMEGA-THCA ingests TCGA and four external GEO cohorts and unifies them into a cross-cohort biomarker panel.",
            "Replication is enforced: every one of the 2,773 retained biomarkers must clear fold-change, FDR, and external-cohort direction checks.",
            "Output feeds directly into the target shortlist consumed by CHEM-THCA and PANEL-THCA.",
        ],
        "capabilities": [
            ("COHORTS", "1,509 samples across TCGA-THCA and four external GEO cohorts."),
            ("REPLICATION", "Effect direction and magnitude checked on held-out cohorts."),
            ("OMICS", "Expression, methylation, and immune deconvolution in one plane."),
            ("OUTPUT", "Ranked biomarker table feeding panel and target-selection stages."),
        ],
        "primary_cta": ("EXPLORE DATASETS", "/reports/html/pages/02_datasets.html"),
        "figures": [
            {"title": "PCA embedding — cross-cohort batch structure", "src": "../figs_interactive/embeddings_pca.html"},
            {"title": "Biomarker replication scatter — TCGA vs external", "src": "../figs_interactive/biomarker_replication_scatter.html"},
        ],
        "reports": [
            ("Analysis summary v3", "../../analysis_summary_v3.md"),
            ("Biomarker analysis", "../../biomarker_analysis.md"),
            ("External validation v3", "../../external_validation_v3.md"),
        ],
    },
    "chem": {
        "index": 2,
        "code": "CHEM-THCA",
        "title": "Chemistry triage for shortlisted targets",
        "summary": [
            "CHEM-THCA pulls ChEMBL bioactivities for each shortlisted target and ranks compounds by pChEMBL, MoA clarity, and literature anchoring.",
            "Cross-references Connectivity-Map reversal signatures to surface compounds that oppose the tumor expression state.",
            "All outputs are computational hit lists; no preclinical or clinical claims are made.",
        ],
        "capabilities": [
            ("COMPOUND SETS", "ChEMBL-anchored with pChEMBL and MoA filters."),
            ("REVERSAL", "CMap score overlay of candidate perturbagens."),
            ("NETWORK", "Protein-protein interaction context around each target."),
            ("OUTPUT", "Ranked compound shortlist per target, published in full."),
        ],
        "primary_cta": ("OPEN DRUG DISCOVERY", "/reports/html/pages/15_drug_discovery.html"),
        "figures": [
            {"title": "CMap top reversers — perturbagen score overlay", "src": "../figs_interactive/cmap_top_reversers.html"},
            {"title": "PPI network — novel targets (QUBO view)", "src": "../figs_interactive/ppi_network_novel_qubo.html"},
        ],
        "reports": [
            ("Biomarker-to-drug report", "../../biomarker_to_drug_report.md"),
            ("Analysis summary v3", "../../analysis_summary_v3.md"),
            ("Next steps v3", "../../next_steps_v3.md"),
        ],
    },
    "panel": {
        "index": 3,
        "code": "PANEL-THCA",
        "title": "Clinical-grade 54-gene panel benchmark",
        "summary": [
            "PANEL-THCA benchmarks a compact 54-gene LogReg panel against the 112-gene reference on TCGA and external GEO cohorts.",
            "The compact panel reaches AUC 1.00 on TCGA internal cross-validation and 0.97–0.89 on four external cohorts.",
            "Released as the reference triage panel for indeterminate-cytology research simulations.",
        ],
        "capabilities": [
            ("SIZE", "54 genes — lean footprint for targeted panels."),
            ("BENCHMARK", "Compared against 112-gene reference on five cohorts."),
            ("PERFORMANCE", "AUC 1.00 internal, 0.97–0.89 external."),
            ("REPRODUCIBILITY", "Full coefficients and risk score exported to TSV."),
        ],
        "primary_cta": ("OPEN PANEL COMPARISON", "/reports/html/pages/08_panel_comparison.html"),
        "figures": [
            {"title": "Panel grouped metrics — compact vs reference", "src": "../figs_interactive/panel_grouped_metrics.html"},
            {"title": "Panel marginal gain — gene-count curve", "src": "../figs_interactive/panel_marginal_gain.html"},
        ],
        "reports": [
            ("ML report v2", "../../ml_report_v2.md"),
            ("Analysis summary v3", "../../analysis_summary_v3.md"),
            ("Paper outline v3", "../../paper_outline_v3.md"),
        ],
    },
    "clinic": {
        "index": 4,
        "code": "CLINIC-THCA",
        "title": "Bethesda-indeterminate triage simulation",
        "summary": [
            "CLINIC-THCA projects panel behavior into Bethesda III and IV cytology workflows as a research-only simulation.",
            "All simulation outputs are explicitly computational and do not constitute clinical guidance.",
            "Built to help clinical collaborators reason about how the compact panel would behave in a triage slot.",
        ],
        "capabilities": [
            ("COHORT MODEL", "Simulated Bethesda III/IV mix anchored to public priors."),
            ("PANEL INPUT", "Uses PANEL-THCA coefficients directly."),
            ("DECISION VIEW", "Operating-point exploration for triage thresholds."),
            ("SCOPE", "Research simulation only — not medical guidance."),
        ],
        "primary_cta": ("OPEN BETHESDA SIMULATOR", "/reports/html/pages/23_bethesda_sim.html"),
        "figures": [
            {"title": "Bethesda simulation overview", "src": "../figs_interactive/panel_cost_utility.html"},
            {"title": "Panel grouped metrics", "src": "../figs_interactive/panel_grouped_metrics.html"},
        ],
        "reports": [
            ("Analysis summary v3", "../../analysis_summary_v3.md"),
            ("Next steps v3", "../../next_steps_v3.md"),
            ("MOBILE UX checklist", "../../MOBILE_UX_CHECKLIST.md"),
        ],
    },
}


# ------------- helpers -------------

def setup_logging() -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    lg = logging.getLogger("rebrand")
    lg.setLevel(logging.INFO)
    for h in list(lg.handlers):
        lg.removeHandler(h)
    fh = logging.FileHandler(LOG_FILE, mode="a")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    lg.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(logging.Formatter("%(message)s"))
    lg.addHandler(sh)
    return lg


def esc(s: str) -> str:
    return html.escape(s or "", quote=True)


# ------------- shared HTML fragments -------------

def head_block(title: str, extra_css: str = "", rel_base: str = "") -> str:
    # rel_base is the relative path prefix to reach HTML_ROOT from the page.
    # index.html uses "", pages/* uses "../".
    return f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">
  <meta name=\"color-scheme\" content=\"dark\">
  <meta name=\"theme-color\" content=\"#000000\">
  <title>{esc(title)}</title>
  <meta name=\"description\" content=\"THYRAI — computational biomarker-to-compound triage for thyroid cancer. Research dashboard only.\">
  <link rel=\"icon\" href=\"{rel_base}assets/img/favicon.svg\" type=\"image/svg+xml\">
  <link rel=\"preload\" href=\"{rel_base}assets/fonts/inter-400.woff2\" as=\"font\" type=\"font/woff2\" crossorigin>
  <link rel=\"preload\" href=\"{rel_base}assets/fonts/inter-600.woff2\" as=\"font\" type=\"font/woff2\" crossorigin>
  <link rel=\"preload\" href=\"{rel_base}assets/fonts/jetbrains-mono-400.woff2\" as=\"font\" type=\"font/woff2\" crossorigin>
  <link rel=\"stylesheet\" href=\"{rel_base}assets/css/brand.css\">
  {extra_css}
</head>
<body class=\"thyrai\">
"""


def topnav(rel_base: str, active: str = "") -> str:
    def acls(name: str) -> str:
        return ' aria-current="page"' if active == name else ""
    return f"""<nav class=\"topnav\" aria-label=\"Primary\">
  <div class=\"topnav__inner\">
    <a class=\"topnav__brand\" href=\"{rel_base}index.html\" aria-label=\"THYRAI home\">
      <img src=\"{rel_base}assets/img/brand/thyrai_wordmark.svg\" alt=\"THYRAI\">
    </a>
    <ul class=\"topnav__links\">
      <li class=\"topnav__item\"><a class=\"topnav__link\"{acls('platform')} href=\"{rel_base}pages/platform_omega.html\">Platform</a>
        <ul class=\"topnav__menu\" role=\"menu\">
          <li><span class=\"topnav__menu-label\">Platforms</span></li>
          <li><a href=\"{rel_base}pages/platform_omega.html\">OMEGA-THCA</a></li>
          <li><a href=\"{rel_base}pages/platform_chem.html\">CHEM-THCA</a></li>
          <li><a href=\"{rel_base}pages/platform_panel.html\">PANEL-THCA</a></li>
          <li><a href=\"{rel_base}pages/platform_clinic.html\">CLINIC-THCA</a></li>
          <li><span class=\"topnav__menu-label\" style=\"margin-top:8px\">Research detail</span></li>
          <li><a href=\"{rel_base}pages/01_overview.html\">01 · Overview</a></li>
          <li><a href=\"{rel_base}pages/15_drug_discovery.html\">15 · Drug discovery</a></li>
          <li><a href=\"{rel_base}pages/17_biomarker_insights.html\">17 · Biomarkers</a></li>
          <li><a href=\"{rel_base}pages/23_bethesda_sim.html\">23 · Bethesda sim</a></li>
          <li><a href=\"{rel_base}pages/99_glossary.html\">99 · Glossary</a></li>
        </ul>
      </li>
      <li><a class=\"topnav__link\"{acls('pipeline')} href=\"{rel_base}pages/pipeline.html\">Pipeline</a></li>
      <li><a class=\"topnav__link\"{acls('publications')} href=\"{rel_base}pages/publications.html\">Publications</a></li>
      <li><a class=\"topnav__link\"{acls('team')} href=\"{rel_base}pages/team.html\">Team</a></li>
      <li><a class=\"topnav__link\"{acls('contact')} href=\"mailto:kukshomr@gmail.com?subject=THYRAI%20contact\">Contact</a></li>
    </ul>
    <a class=\"topnav__cta\" href=\"mailto:kukshomr@gmail.com?subject=THYRAI%20access%20request\">Request Access</a>
    <button class=\"topnav__hamburger\" type=\"button\" aria-label=\"Toggle navigation\" aria-expanded=\"false\"><span></span></button>
  </div>
</nav>
"""


def footer(rel_base: str, include_legacy_link: bool = True) -> str:
    legacy = (
        f'<li><a href="{rel_base}index_v2_legacy.html">v2 Research mode</a></li>'
        if include_legacy_link else ""
    )
    return f"""<footer class=\"site-footer\">
  <div class=\"container\">
    <div class=\"site-footer__grid\">
      <div class=\"site-footer__col site-footer__brand\">
        <img src=\"{rel_base}assets/img/brand/thyrai_wordmark.svg\" alt=\"THYRAI\">
        <p>Research dashboard for computational thyroid-cancer triage. Computational triage only. Not medical advice.</p>
      </div>
      <div class=\"site-footer__col\">
        <h4>Platform</h4>
        <ul>
          <li><a href=\"{rel_base}pages/platform_omega.html\">OMEGA-THCA</a></li>
          <li><a href=\"{rel_base}pages/platform_chem.html\">CHEM-THCA</a></li>
          <li><a href=\"{rel_base}pages/platform_panel.html\">PANEL-THCA</a></li>
          <li><a href=\"{rel_base}pages/platform_clinic.html\">CLINIC-THCA</a></li>
        </ul>
      </div>
      <div class=\"site-footer__col\">
        <h4>Research</h4>
        <ul>
          <li><a href=\"{rel_base}pages/pipeline.html\">Pipeline</a></li>
          <li><a href=\"{rel_base}pages/publications.html\">Publications</a></li>
          <li><a href=\"{rel_base}pages/press.html\">Press</a></li>
          <li><a href=\"{rel_base}pages/99_glossary.html\">Glossary</a></li>
          {legacy}
        </ul>
      </div>
      <div class=\"site-footer__col\">
        <h4>Contact</h4>
        <ul>
          <li><a href=\"mailto:kukshomr@gmail.com\">kukshomr@gmail.com</a></li>
          <li><a href=\"{rel_base}pages/team.html\">Team</a></li>
          <li><a href=\"mailto:kukshomr@gmail.com?subject=THYRAI%20licensing\">Licensing</a></li>
        </ul>
      </div>
    </div>
    <div class=\"site-footer__meta\">
      <span>&copy; {YEAR} THYRAI Research</span>
      <span class=\"site-footer__disclaimer\">{LEGAL_FOOTER}</span>
    </div>
  </div>
</footer>
<script src=\"{rel_base}assets/js/brand.js\" defer></script>
</body>
</html>
"""


# ------------- latest updates from markdown -------------

def latest_markdown_updates(limit: int = 3) -> list[dict]:
    reports_dir = PROJECT_ROOT / "reports"
    md_files = [p for p in reports_dir.glob("*.md") if p.is_file()]
    md_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    picked = md_files[:limit]
    out = []
    for p in picked:
        title = p.stem.replace("_", " ").title()
        excerpt = ""
        try:
            txt = p.read_text(errors="ignore")
            # strip headings/blank/frontmatter
            cleaned_lines = []
            for line in txt.splitlines():
                s = line.strip()
                if not s:
                    continue
                if s.startswith(("#", "-", "|", ">", "*", "`")):
                    continue
                cleaned_lines.append(s)
                if len(" ".join(cleaned_lines)) > 200:
                    break
            excerpt = " ".join(cleaned_lines)[:140]
            if len(excerpt) == 140:
                excerpt += "…"
        except Exception:
            excerpt = ""
        mtime = datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d")
        out.append({
            "title": title,
            "excerpt": excerpt,
            "date": mtime,
            "href": f"../{p.name}",
            "abs_href": f"/reports/{p.name}",
        })
    return out


# ------------- page builders -------------

def build_index(log: logging.Logger) -> None:
    updates = latest_markdown_updates(3)
    update_cards = "\n".join(
        f"""<a class=\"update-card\" href=\"/reports/{esc(u['href'].split('/')[-1])}\">
  <div class=\"update-card__date\">{esc(u['date'])}</div>
  <div class=\"update-card__title\">{esc(u['title'])}</div>
  <p class=\"update-card__excerpt\">{esc(u['excerpt'])}</p>
</a>"""
        for u in updates
    )

    platform_tiles_html = "\n".join(
        f"""<a class=\"platform-tile\" href=\"{esc(t['href'])}\">
  <span>
    <span class=\"section-label\">Platform · {i+1}/4</span>
    <span class=\"platform-tile__codename\">{esc(t['code'])}</span>
  </span>
  <span>
    <p class=\"platform-tile__desc\">{esc(t['desc'])}</p>
    <span class=\"platform-tile__arrow\">Explore &rarr;</span>
  </span>
</a>"""
        for i, t in enumerate(PLATFORM_TILES)
    )

    press_logos = "\n".join(
        f'<img src="assets/img/brand/press_{n}.svg" alt="{label}" loading="lazy">'
        for n, label in [
            ("tcga", "TCGA"), ("geo", "GEO"), ("chembl", "ChEMBL"),
            ("clinvar", "ClinVar"), ("pubmed", "PubMed"),
        ]
    )

    # Mini pipeline 8x6 preview — first-column filled for each target; CYP1B1 also stage 2
    with open(ASSETS / "data/pipeline.json") as f:
        pipeline = json.load(f)
    mini_cells = []
    for _row in pipeline[:8]:
        # Only show 5 cells of the first row as a representative strip
        pass
    # Build a 6-column strip summarizing how many targets sit in each stage
    stage_counts = {"target_id": 0, "hit_to_lead": 0}
    for p in pipeline:
        stage_counts[p["stage"]] = stage_counts.get(p["stage"], 0) + 1
    mini_strip_cells = "\n".join(
        f'<div class="mini-pipeline__cell {"mini-pipeline__cell--filled-1" if i == 0 else ("mini-pipeline__cell--filled-2" if (i == 1 and stage_counts.get("hit_to_lead", 0) > 0) else "")}"></div>'
        for i in range(6)
    )

    head = head_block("THYRAI — Generative AI for thyroid cancer", rel_base="")

    body = f"""{topnav(rel_base='', active='home')}

<section class=\"hero\">
  <div class=\"container hero__inner\">
    <span class=\"eyebrow\">Generative AI for thyroid cancer</span>
    <h1>From 2,773 public multi-omics biomarkers to drug-discoverable targets.</h1>
    <p class=\"hero__lede\">THYRAI is an open research dashboard that turns replicated thyroid-cancer biomarker evidence into an eight-target, fifteen-compound computational pipeline — released with every table, figure, and report behind it.</p>
    <div class=\"hero__ctas\">
      <a class=\"btn-primary\" href=\"pages/pipeline.html\">Explore pipeline</a>
      <a class=\"btn-ghost\" href=\"pages/platform_omega.html\">Platform overview</a>
    </div>
  </div>
</section>

<section aria-label=\"Featured in\">
  <div class=\"container\">
    <div class=\"press-bar\">
      <span class=\"press-bar__label\">Featured in</span>
      <div class=\"press-bar__logos\">
        {press_logos}
      </div>
    </div>
  </div>
</section>

<section class=\"section\" aria-labelledby=\"platform-h\">
  <div class=\"container\">
    <span class=\"section-label\">Platform</span>
    <h2 id=\"platform-h\">Four computational layers, one thyroid-cancer program.</h2>
    <div class=\"platform-grid\" style=\"margin-top:32px\">
      {platform_tiles_html}
    </div>
  </div>
</section>

<section class=\"section section--tight\" aria-labelledby=\"pipeline-preview-h\">
  <div class=\"container\">
    <span class=\"eyebrow\">Pipeline preview</span>
    <h2 id=\"pipeline-preview-h\">15 compounds across 8 targets.</h2>
    <div class=\"mini-pipeline\" aria-hidden=\"true\">{mini_strip_cells}</div>
    <a class=\"btn-ghost\" href=\"pages/pipeline.html\">See full pipeline &rarr;</a>
  </div>
</section>

<section aria-labelledby=\"stats-h\">
  <div class=\"container\">
    <div class=\"key-stats\">
      <div class=\"key-stats__item\"><div class=\"key-stats__value\">1,509</div><div class=\"key-stats__label\">Samples</div></div>
      <div class=\"key-stats__item\"><div class=\"key-stats__value\">2,773</div><div class=\"key-stats__label\">Biomarkers</div></div>
      <div class=\"key-stats__item\"><div class=\"key-stats__value\">8</div><div class=\"key-stats__label\">Targets</div></div>
      <div class=\"key-stats__item\"><div class=\"key-stats__value\">15</div><div class=\"key-stats__label\">Compounds</div></div>
    </div>
  </div>
</section>

<section class=\"pull-quote\" aria-label=\"Authority quote\">
  <div class=\"container\">
    <p class=\"pull-quote__text\">Classical LogReg on a 54-gene panel reaches AUC = 1.00 on TCGA internal CV and 0.97–0.89 on external GEO cohorts.</p>
    <div class=\"pull-quote__attr\">— THCA v2 Technical Report · 2026</div>
  </div>
</section>

<section class=\"section\" aria-labelledby=\"updates-h\">
  <div class=\"container\">
    <span class=\"eyebrow\">Latest updates</span>
    <h2 id=\"updates-h\">Recent reports from the research thread.</h2>
    <div class=\"updates-grid\" style=\"margin-top:32px\">
      {update_cards}
    </div>
  </div>
</section>

{footer(rel_base='')}"""

    (HTML_ROOT / "index.html").write_text(head + body)
    log.info("wrote index.html")


# ------------- pipeline page -------------

STAGE_ORDER = ["target_id", "hit_to_lead", "lead_opt", "ind", "phase1", "phase2"]
STAGE_CLASS = {
    "target_id": "pipeline-bar--target",
    "hit_to_lead": "pipeline-bar--hit",
    "lead_opt": "pipeline-bar--lead",
    "ind": "pipeline-bar--ind",
    "phase1": "pipeline-bar--phase1",
    "phase2": "pipeline-bar--phase2",
}


def build_pipeline_page(log: logging.Logger) -> None:
    with open(ASSETS / "data/pipeline.json") as f:
        pipeline = json.load(f)

    legend_cells = [
        ("Target ID", ""),
        ("Hit-to-lead", ""),
        ("Lead opt", ""),
        ("IND-enabling", ""),
        ("Phase 1", "n/a — computational triage"),
        ("Phase 2", "n/a — computational triage"),
    ]
    legend_html = "".join(
        f"<div class=\"pipeline-legend__cell\">{esc(l)}{('<small>' + esc(s) + '</small>') if s else ''}</div>"
        for l, s in legend_cells
    )

    rows_html = []
    for p in pipeline:
        target_stage_idx = STAGE_ORDER.index(p["stage"]) if p["stage"] in STAGE_ORDER else 0
        bars = []
        for i, s in enumerate(STAGE_ORDER):
            if i <= target_stage_idx:
                cls = STAGE_CLASS[s]
            else:
                cls = "pipeline-bar--empty"
            bars.append(f'<div class="pipeline-bar {cls}" aria-label="{esc(s.replace("_"," "))}"></div>')
        bars_html = "".join(bars)
        pchembl_txt = f"{p['best_pchembl']:.2f}" if p.get('best_pchembl') else "—"
        tooltip = (
            f"Compounds {p['compounds']} · best pChEMBL {pchembl_txt} · literature {p['literature_count']}"
        )
        mailto = (
            f"mailto:kukshomr@gmail.com?subject=THYRAI%20licensing%20{esc(p['target'])}"
        )
        rows_html.append(f"""<div class=\"pipeline-row\" data-target=\"{esc(p['target'])}\" tabindex=\"0\">
  <div class=\"pipeline-row__target\">
    <span class=\"pipeline-row__target-symbol\">{esc(p['target'])}</span>
    <span class=\"pipeline-row__target-indication\">{esc(p['indication'])}</span>
  </div>
  <div class=\"pipeline-row__bars\">{bars_html}</div>
  <a class=\"pipeline-row__license\" href=\"{mailto}\">
    <svg viewBox=\"0 0 20 20\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"1.5\" aria-hidden=\"true\"><path d=\"M2 3h2l1.5 10h10L18 5H5\"/><circle cx=\"8\" cy=\"16\" r=\"1.2\"/><circle cx=\"15\" cy=\"16\" r=\"1.2\"/></svg>
    Available for licensing
  </a>
  <div class=\"pipeline-tooltip\" role=\"tooltip\">{esc(tooltip)}</div>
</div>""")

    # Target detail cards
    detail_cards = []
    for p in pipeline:
        classification = p['classification'].replace('_', ' ').title()
        best_pchembl = f"{p['best_pchembl']:.2f}" if p.get('best_pchembl') else "—"
        top_compound = p.get('top_compound') or "—"
        raw_rationale = p.get('rationale') or ""
        if raw_rationale:
            cleaned = re.sub(r"\*\*(.+?)\*\*", r"\1", raw_rationale)
            cleaned = cleaned.strip()
        else:
            cleaned = (
                f"{p['target']} sits in the {classification.lower()} tier of the THYRAI shortlist "
                f"for {p['indication'].lower()}. Evidence summary: {p['compounds']} ChEMBL compounds, "
                f"{p['literature_count']} anchoring literature records, best pChEMBL {best_pchembl}. "
                "Released as a computational hypothesis, not a preclinical claim."
            )
        # Trim to 2-3 sentences
        sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+", cleaned) if s.strip()]
        trimmed = " ".join(sentences[:3])
        detail_cards.append(f"""<div class=\"target-card\">
  <span class=\"section-label\">{esc(classification)}</span>
  <h3>{esc(p['target'])}</h3>
  <p class=\"target-card__rationale\">{esc(trimmed)}</p>
  <div class=\"target-card__stats\">
    <div><span class=\"target-card__stat-label\">Compounds</span><span class=\"target-card__stat-value\">{esc(str(p['compounds']))}</span></div>
    <div><span class=\"target-card__stat-label\">Best pChEMBL</span><span class=\"target-card__stat-value\">{esc(best_pchembl)}</span></div>
    <div><span class=\"target-card__stat-label\">Literature</span><span class=\"target-card__stat-value\">{esc(str(p['literature_count']))}</span></div>
  </div>
  <a class=\"target-card__link\" href=\"15_drug_discovery.html\">Explore compound space &rarr;</a>
</div>""")

    head = head_block("THYRAI — Therapeutic pipeline", rel_base="../")
    body = f"""{topnav(rel_base='../', active='pipeline')}

<header class=\"page-header\">
  <div class=\"container\">
    <span class=\"eyebrow\">Therapeutic pipeline</span>
    <h1>Computational discovery pipeline</h1>
    <p class=\"page-header__subtitle\">Eight targets surfaced from 2,773 replicated biomarkers. Every bar below reflects computational evidence only. No compound has entered preclinical validation, IND-enabling studies, or clinical trials.</p>
  </div>
</header>

<section class=\"section\">
  <div class=\"container\">
    <div class=\"pipeline-legend\" role=\"list\">{legend_html}</div>
    <div class=\"pipeline-grid\">
      {''.join(rows_html)}
    </div>
  </div>
</section>

<section class=\"section\" aria-labelledby=\"target-detail-h\">
  <div class=\"container\">
    <span class=\"eyebrow\">Target rationale</span>
    <h2 id=\"target-detail-h\">Eight thyroid-cancer targets · classification, evidence, next step.</h2>
    <div class=\"target-cards\" style=\"margin-top:32px\">
      {''.join(detail_cards)}
    </div>
  </div>
</section>

<section>
  <div class=\"container\">
    <div class=\"disclaimer-band\" role=\"note\">
      <strong>Research scope</strong>
      All stages reflect computational evidence only. No compound has entered preclinical validation, IND-enabling studies, or clinical trials. This pipeline view is a research dashboard, not a drug-development claim.
    </div>
  </div>
</section>

{footer(rel_base='../')}"""

    (PAGES / "pipeline.html").write_text(head + body)
    log.info("wrote pages/pipeline.html")


# ------------- platform pages -------------

def build_platform_page(key: str, log: logging.Logger) -> None:
    cfg = PLATFORM_DETAIL[key]
    title = f"THYRAI — {cfg['code']}"

    cap_items = "\n".join(
        f"<li><span class=\"capability-list__label\">{esc(label)}</span><p>{esc(text)}</p></li>"
        for label, text in cfg["capabilities"]
    )
    fig_items = "\n".join(
        f"""<figure>
  <figcaption class=\"section-label\">{esc(f['title'])}</figcaption>
  <iframe class=\"fig-embed\" data-src=\"{esc(f['src'])}\" title=\"{esc(f['title'])}\" loading=\"lazy\"></iframe>
</figure>"""
        for f in cfg["figures"]
    )
    report_items = "\n".join(
        f"<li><a href=\"{esc(href)}\"><span>{esc(label)}</span><span class=\"report-list__date\">MD &rarr;</span></a></li>"
        for label, href in cfg["reports"]
    )
    summary_html = "".join(f"<p class=\"hero__lede\">{esc(s)}</p>" for s in cfg["summary"])

    head = head_block(title, rel_base="../")
    body = f"""{topnav(rel_base='../', active='platform')}

<section class=\"hero hero--md\">
  <div class=\"container hero__inner\">
    <span class=\"eyebrow\">Platform · {cfg['index']}/4</span>
    <div class=\"platform-tile__codename\" style=\"margin-bottom:24px\">{esc(cfg['code'])}</div>
    <h1>{esc(cfg['title'])}</h1>
    {summary_html}
    <div class=\"hero__ctas\"><a class=\"btn-primary\" href=\"..{esc(cfg['primary_cta'][1])[len('/reports/html'):] if cfg['primary_cta'][1].startswith('/reports/html') else cfg['primary_cta'][1]}\">{esc(cfg['primary_cta'][0])}</a><a class=\"btn-ghost\" href=\"pipeline.html\">See pipeline</a></div>
  </div>
</section>

<section class=\"section section--tight\" aria-labelledby=\"cap-h\">
  <div class=\"container\">
    <span class=\"section-label\">Capabilities</span>
    <h2 id=\"cap-h\">What this layer does.</h2>
    <ul class=\"capability-list\">{cap_items}</ul>
  </div>
</section>

<section class=\"section\" aria-labelledby=\"fig-h\">
  <div class=\"container\">
    <span class=\"section-label\">Embedded figures</span>
    <h2 id=\"fig-h\">Interactive outputs from this layer.</h2>
    {fig_items}
  </div>
</section>

<section class=\"section section--tight\" aria-labelledby=\"rep-h\">
  <div class=\"container\">
    <span class=\"section-label\">Related reports</span>
    <h2 id=\"rep-h\">Read the full research thread.</h2>
    <ul class=\"report-list\">{report_items}</ul>
  </div>
</section>

{footer(rel_base='../')}"""

    # Fix the primary CTA path if it pointed at /reports/html/...
    cta_path = cfg["primary_cta"][1]
    if cta_path.startswith("/reports/html/"):
        cta_fixed = cta_path[len("/reports/html/"):]  # e.g. pages/02_datasets.html
        # Relative from pages/platform_x.html to pages/02_datasets.html is same directory
        if cta_fixed.startswith("pages/"):
            cta_fixed = cta_fixed[len("pages/"):]
        body = body.replace(
            f"..{cta_path[len('/reports/html'):] if cta_path.startswith('/reports/html') else cta_path}",
            cta_fixed,
        )

    (PAGES / f"platform_{key}.html").write_text(head + body)
    log.info("wrote pages/platform_%s.html", key)


# ------------- publications / press / team -------------

def build_publications_page(log: logging.Logger) -> None:
    reports_dir = PROJECT_ROOT / "reports"
    md_files = sorted(reports_dir.glob("*.md"), key=lambda p: p.stat().st_mtime, reverse=True)

    tech = []
    analysis = []
    next_ = []
    for p in md_files:
        stem = p.stem.lower()
        if "next_steps" in stem or "paper_outline" in stem:
            next_.append(p)
        elif "report" in stem or "audit" in stem or "ml_report" in stem or "debug" in stem or "validation" in stem or "qa" in stem:
            tech.append(p)
        else:
            analysis.append(p)

    def render_group(title: str, items: list) -> str:
        if not items:
            return ""
        li = []
        for p in items:
            date = datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d")
            nice_title = p.stem.replace("_", " ").title()
            excerpt_line = ""
            try:
                for line in p.read_text(errors="ignore").splitlines():
                    s = line.strip()
                    if not s:
                        continue
                    if s.startswith(("#", "-", "|", ">", "*", "`")):
                        continue
                    excerpt_line = s[:160]
                    break
            except Exception:
                pass
            li.append(
                f"<li><a href=\"../../{esc(p.name)}\"><span>{esc(date)} &middot; {esc(nice_title)}<br><span style=\"color:var(--text-muted);font-size:13px\">{esc(excerpt_line)}</span></span><span class=\"report-list__date\">MD &rarr;</span></a></li>"
            )
        return f"""<section class=\"section section--tight\">
  <div class=\"container\">
    <span class=\"section-label\">{esc(title)}</span>
    <h2>{esc(title)}</h2>
    <ul class=\"report-list\">{''.join(li)}</ul>
  </div>
</section>"""

    head = head_block("THYRAI — Publications", rel_base="../")
    body = f"""{topnav(rel_base='../', active='publications')}

<header class=\"page-header\">
  <div class=\"container\">
    <span class=\"eyebrow\">Publications</span>
    <h1>Research reports, summaries, and next-step memos.</h1>
    <p class=\"page-header__subtitle\">Auto-generated index of the markdown research thread under <code>reports/</code>. Sorted by most recent first.</p>
  </div>
</header>

{render_group('Technical Reports', tech)}
{render_group('Analysis Summaries', analysis)}
{render_group('Next Steps', next_)}

{footer(rel_base='../')}"""
    (PAGES / "publications.html").write_text(head + body)
    log.info("wrote pages/publications.html")


def build_press_page(log: logging.Logger) -> None:
    slots = "\n".join(
        f"""<div class=\"slot-card\">
  <div class=\"section-label\">Slot {i+1}</div>
  <p style=\"color:var(--text-muted);font-family:var(--font-sans);text-transform:none;letter-spacing:0;font-size:14px;margin-top:8px\">Add citation — featured media, invited talk, or partner mention.</p>
</div>"""
        for i in range(5)
    )
    head = head_block("THYRAI — Press", rel_base="../")
    body = f"""{topnav(rel_base='../', active='press')}

<header class=\"page-header\">
  <div class=\"container\">
    <span class=\"eyebrow\">Press</span>
    <h1>Featured citations, media coverage, and invited talks.</h1>
    <p class=\"page-header__subtitle\">Press slots will populate as coverage lands. Every entry links back to the underlying report.</p>
  </div>
</header>

<section class=\"section\">
  <div class=\"container\">
    <div class=\"platform-grid\">{slots}</div>
  </div>
</section>

{footer(rel_base='../')}"""
    (PAGES / "press.html").write_text(head + body)
    log.info("wrote pages/press.html")


def build_team_page(log: logging.Logger) -> None:
    advisor_slots = "\n".join(
        f"""<div class=\"slot-card\">
  <div class=\"section-label\">Advisor · {i+1}</div>
  <p style=\"color:var(--text-muted);font-family:var(--font-sans);text-transform:none;letter-spacing:0;font-size:14px;margin-top:8px\">Name &middot; affiliation &middot; role. Placeholder slot.</p>
</div>"""
        for i in range(3)
    )
    head = head_block("THYRAI — Team", rel_base="../")
    body = f"""{topnav(rel_base='../', active='team')}

<header class=\"page-header\">
  <div class=\"container\">
    <span class=\"eyebrow\">Team</span>
    <h1>Lead researcher and advisors.</h1>
    <p class=\"page-header__subtitle\">Small independent research group. Contact the lead to discuss data access, collaboration, or licensing.</p>
  </div>
</header>

<section class=\"section\">
  <div class=\"container\">
    <div class=\"platform-grid\">
      <div class=\"card-flat\">
        <span class=\"section-label\">Lead</span>
        <h3 style=\"font-family:var(--font-mono);font-weight:500;font-size:26px;margin-top:8px\">THYRAI Principal</h3>
        <p>Lead on biomarker discovery, panel benchmarking, and pipeline triage. Correspondence address below.</p>
        <a class=\"btn-ghost\" href=\"mailto:kukshomr@gmail.com\">kukshomr@gmail.com</a>
      </div>
      {advisor_slots}
    </div>
  </div>
</section>

{footer(rel_base='../')}"""
    (PAGES / "team.html").write_text(head + body)
    log.info("wrote pages/team.html")


# ------------- main -------------

def run_pipeline_generator(log: logging.Logger) -> None:
    script = PROJECT_ROOT / "scripts/generate_pipeline_data.py"
    python = PROJECT_ROOT / ".venv/bin/python"
    r = subprocess.run([str(python), str(script)], capture_output=True, text=True)
    log.info("generate_pipeline_data.py rc=%d", r.returncode)
    log.info(r.stdout.strip())
    if r.stderr.strip():
        log.info("stderr: %s", r.stderr.strip())
    if r.returncode != 0:
        raise SystemExit(f"pipeline data generator failed: {r.stderr}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=str(PROJECT_ROOT))
    ap.add_argument("--brand-name", default="THYRAI")
    ap.add_argument("--accent", default="#F5A623")
    args = ap.parse_args()

    log = setup_logging()
    log.info("=== THYRAI rebrand start brand=%s accent=%s ===", args.brand_name, args.accent)
    log.info("project_root=%s", args.project_root)

    # Ensure dirs
    for d in [
        HTML_ROOT / "assets/img/brand",
        HTML_ROOT / "assets/data",
        HTML_ROOT / "assets/css",
        HTML_ROOT / "assets/js",
        HTML_ROOT / "assets/fonts",
        PAGES,
    ]:
        d.mkdir(parents=True, exist_ok=True)

    # Step 1: pipeline data
    run_pipeline_generator(log)

    # Step 2: build all pages
    build_index(log)
    build_pipeline_page(log)
    for k in ["omega", "chem", "panel", "clinic"]:
        build_platform_page(k, log)
    build_publications_page(log)
    build_press_page(log)
    build_team_page(log)

    # Step 3: overwrite main.css to defer to brand.css + legacy. To stay CSP-safe and preserve
    # existing 01-28 pages that reference main.css, main.css now just imports its legacy content.
    main_css = HTML_ROOT / "assets/css/main.css"
    main_css_legacy = HTML_ROOT / "assets/css/main_v2_legacy.css"
    if main_css_legacy.exists():
        main_css.write_text(
            "/* THYRAI main.css — legacy styles preserved for pages 01-28; brand styles live in brand.css */\n"
            "@import url('./main_v2_legacy.css');\n"
        )
        log.info("rewired main.css -> main_v2_legacy.css import")

    log.info("=== THYRAI rebrand complete ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
