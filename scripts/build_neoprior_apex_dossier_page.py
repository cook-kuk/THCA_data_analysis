#!/usr/bin/env python3
"""Build an editor/advisor-facing NEO-PRIOR apex dossier page.

The page is a static HTML hub for the high-impact strategy:
NEO-PRIOR standard + NeoBench-Vax benchmark + CROSS-Neo implementation +
product-assisted triage. It is designed as a showable link, not a new
scientific claim. It keeps the conservative boundary explicit throughout.
"""

from __future__ import annotations

import html
import shutil
from pathlib import Path
from textwrap import dedent


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
RESULT = ROOT / "project/results/cross_neo_two_paper_strategy_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
ASSET_DIR = HUB / "assets/cross_neo_apex"
OUT = HUB / "neoprior_neobench_crossneo_apex_dossier.html"


def read(path: Path, fallback: str = "") -> str:
    if path.exists():
        return path.read_text(encoding="utf-8", errors="replace")
    return fallback


def copy_assets() -> list[tuple[str, str]]:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    candidates = [
        RESULT / "impact_apex/figures/fig_apex_one_slide_launch_story.png",
        RESULT / "impact_apex/figures/fig_apex_neobench_flywheel.png",
        RESULT / "impact_apex/figures/fig_apex_90_day_roadmap.png",
        RESULT / "impact_apex/figures/fig_apex_claim_ladder.png",
        RESULT / "impact_ultra/figures/fig_ultra_submission_routes.png",
        RESULT / "impact_ultra/figures/fig_ultra_neobench_public_infrastructure.png",
        RESULT / "impact_max/figures/fig_editor_one_screen_story.png",
        RESULT / "impact_max/figures/fig_neoprior_neobench_crossneo_ecosystem.png",
        RESULT / "figures/two_paper_impact_map.png",
    ]
    copied: list[tuple[str, str]] = []
    for src in candidates:
        if src.exists():
            dst = ASSET_DIR / src.name
            shutil.copy2(src, dst)
            copied.append((src.name, f"assets/cross_neo_apex/{src.name}"))
    return copied


def md_excerpt(path: Path, max_chars: int = 1400) -> str:
    text = read(path, "Missing local file.")
    text = text.replace("# ", "").replace("## ", "")
    if len(text) > max_chars:
        text = text[:max_chars].rstrip() + "\n..."
    return html.escape(text)


def make_html(assets: list[tuple[str, str]]) -> str:
    asset_cards = "\n".join(
        f'<figure><img src="{html.escape(rel)}" alt="{html.escape(name)}"><figcaption>{html.escape(name)}</figcaption></figure>'
        for name, rel in assets
    )
    apex_memo = md_excerpt(RESULT / "impact_apex/APEX_IMPACT_LAUNCH_MEMO.md")
    abstracts = md_excerpt(RESULT / "impact_apex/JOURNAL_SPECIFIC_ABSTRACTS.md", 1900)
    ultra = md_excerpt(RESULT / "impact_ultra/ULTRA_IMPACT_EDITORIAL_DOSSIER.md", 1800)
    checklist = md_excerpt(RESULT / "impact_ultra/NEO_PRIOR_TWO_PAGE_CHECKLIST.md", 1700)
    demo = md_excerpt(RESULT / "impact_apex/WEDNESDAY_DEMO_EXECUTIVE_SCRIPT.md", 1200)

    return dedent(
        f"""
        <!doctype html>
        <html lang="en">
        <head>
          <meta charset="utf-8" />
          <meta name="viewport" content="width=device-width, initial-scale=1" />
          <title>NEO-PRIOR / NeoBench-Vax / CROSS-Neo Apex Dossier</title>
          <style>
            :root {{
              --bg: #070b13;
              --panel: #101827;
              --panel2: #0f1f33;
              --ink: #f8fafc;
              --muted: #a7b2c7;
              --line: rgba(148, 163, 184, .25);
              --cyan: #67e8f9;
              --green: #86efac;
              --amber: #fde68a;
              --rose: #fda4af;
              --blue: #93c5fd;
            }}
            * {{ box-sizing: border-box; }}
            html {{ scroll-behavior: smooth; }}
            body {{
              margin: 0;
              background:
                radial-gradient(circle at 20% 0%, rgba(14, 116, 144, .34), transparent 36%),
                linear-gradient(180deg, #07101e 0%, var(--bg) 48%, #05070d 100%);
              color: var(--ink);
              font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
              line-height: 1.55;
            }}
            a {{ color: var(--cyan); text-decoration: none; }}
            .hero {{
              min-height: 84vh;
              display: grid;
              align-items: end;
              padding: 72px min(7vw, 96px) 56px;
              border-bottom: 1px solid var(--line);
            }}
            .kicker {{
              color: var(--cyan);
              text-transform: uppercase;
              letter-spacing: .13em;
              font-weight: 800;
              font-size: 12px;
            }}
            h1 {{
              max-width: 1080px;
              font-size: clamp(42px, 6vw, 88px);
              line-height: .95;
              letter-spacing: 0;
              margin: 16px 0 18px;
            }}
            .lead {{
              max-width: 920px;
              font-size: clamp(18px, 2vw, 25px);
              color: var(--muted);
              margin: 0;
            }}
            .strip {{
              display: grid;
              grid-template-columns: repeat(4, minmax(0, 1fr));
              gap: 12px;
              max-width: 1140px;
              margin-top: 34px;
            }}
            .stat {{
              border: 1px solid var(--line);
              background: rgba(16, 24, 39, .74);
              padding: 16px;
              border-radius: 8px;
            }}
            .stat b {{
              display: block;
              font-size: 22px;
              color: #fff;
              margin-bottom: 4px;
            }}
            .layout {{
              display: grid;
              grid-template-columns: 280px minmax(0, 1fr);
              gap: 32px;
              padding: 34px min(7vw, 96px) 80px;
            }}
            nav {{
              position: sticky;
              top: 18px;
              align-self: start;
              border: 1px solid var(--line);
              border-radius: 8px;
              background: rgba(16, 24, 39, .78);
              padding: 18px;
              max-height: calc(100vh - 36px);
              overflow: auto;
            }}
            nav b {{ display: block; margin-bottom: 10px; }}
            nav a {{ display: block; padding: 7px 0; color: var(--muted); font-size: 14px; }}
            section {{
              border-top: 1px solid var(--line);
              padding: 42px 0;
            }}
            section:first-child {{ border-top: 0; padding-top: 0; }}
            h2 {{
              font-size: clamp(26px, 3vw, 38px);
              margin: 0 0 18px;
              letter-spacing: 0;
            }}
            h3 {{ margin: 0 0 10px; font-size: 19px; }}
            p {{ color: var(--muted); }}
            .grid3 {{
              display: grid;
              grid-template-columns: repeat(3, minmax(0, 1fr));
              gap: 14px;
            }}
            .grid2 {{
              display: grid;
              grid-template-columns: repeat(2, minmax(0, 1fr));
              gap: 14px;
            }}
            .card {{
              border: 1px solid var(--line);
              background: rgba(16, 24, 39, .78);
              border-radius: 8px;
              padding: 20px;
            }}
            .card strong {{ color: #fff; }}
            .good {{ color: var(--green); }}
            .warn {{ color: var(--amber); }}
            .danger {{ color: var(--rose); }}
            table {{
              width: 100%;
              border-collapse: collapse;
              font-size: 14px;
              overflow: hidden;
              border: 1px solid var(--line);
              border-radius: 8px;
            }}
            th, td {{
              border-bottom: 1px solid var(--line);
              padding: 12px;
              text-align: left;
              vertical-align: top;
            }}
            th {{ color: var(--cyan); background: rgba(15, 31, 51, .95); }}
            tr:last-child td {{ border-bottom: 0; }}
            .figgrid {{
              display: grid;
              grid-template-columns: repeat(2, minmax(0, 1fr));
              gap: 16px;
            }}
            figure {{
              margin: 0;
              border: 1px solid var(--line);
              background: #fff;
              border-radius: 8px;
              overflow: hidden;
            }}
            figure img {{
              display: block;
              width: 100%;
              height: auto;
              background: #fff;
            }}
            figcaption {{
              color: #334155;
              font-size: 12px;
              padding: 8px 10px;
              border-top: 1px solid #e2e8f0;
              background: #f8fafc;
            }}
            pre {{
              white-space: pre-wrap;
              border: 1px solid var(--line);
              background: rgba(15, 31, 51, .70);
              border-radius: 8px;
              padding: 18px;
              color: #dbeafe;
              overflow: auto;
              font-size: 13px;
            }}
            .decision {{
              border-left: 5px solid var(--amber);
              background: rgba(253, 230, 138, .08);
            }}
            .footer {{
              color: #64748b;
              font-size: 13px;
              padding-top: 32px;
            }}
            @media (max-width: 980px) {{
              .layout {{ grid-template-columns: 1fr; }}
              nav {{ position: static; }}
              .strip, .grid3, .grid2, .figgrid {{ grid-template-columns: 1fr; }}
            }}
          </style>
        </head>
        <body>
          <header class="hero">
            <div>
              <div class="kicker">NEO-PRIOR · NeoBench-Vax · CROSS-Neo</div>
              <h1>Auditable AI for the cancer vaccine shortlist.</h1>
              <p class="lead">A field-standard launch package for deciding whether neoantigen AI rankings are trustworthy enough to guide assay, manufacturing and clinical prioritization.</p>
              <div class="strip">
                <div class="stat"><b>Standard</b>NEO-PRIOR reporting contract</div>
                <div class="stat"><b>Benchmark</b>NeoBench-Vax public infrastructure</div>
                <div class="stat"><b>Method</b>CROSS-Neo first implementation</div>
                <div class="stat"><b>Boundary</b>Internal/locked retrospective only</div>
              </div>
            </div>
          </header>
          <div class="layout">
            <nav>
              <b>Dossier</b>
              <a href="#tl-dr">TL;DR</a>
              <a href="#why-now">Why Now</a>
              <a href="#architecture">Three-Layer Architecture</a>
              <a href="#evidence">Evidence Snapshot</a>
              <a href="#figures">Figures</a>
              <a href="#papers">Two Papers</a>
              <a href="#product">Business Demo</a>
              <a href="#claim">Claim Boundary</a>
              <a href="#send">What To Send</a>
              <a href="#sources">Sources</a>
            </nav>
            <main>
              <section id="tl-dr">
                <h2>TL;DR</h2>
                <div class="card decision">
                  <h3>Editor-facing thesis</h3>
                  <p><strong>Personalized cancer vaccines need auditable AI for the shortlist that decides what gets manufactured.</strong></p>
                  <p>This is stronger than a single-model claim because it creates a reporting standard, a benchmark infrastructure, and a first implementation under strict claim boundaries.</p>
                </div>
              </section>

              <section id="why-now">
                <h2>Why Now</h2>
                <div class="grid3">
                  <div class="card"><h3>Clinical pressure</h3><p>Personalized vaccine teams must choose a small number of candidates from a larger mutation-derived pool.</p></div>
                  <div class="card"><h3>Benchmark gap</h3><p>Published predictors mix tasks, public overlaps, HLA shortcuts, source bias, negative-label ambiguity and AUROC-heavy reporting.</p></div>
                  <div class="card"><h3>Actionable standard</h3><p>NEO-PRIOR specifies what evidence any model should report before its shortlist is treated as vaccine-ready.</p></div>
                </div>
              </section>

              <section id="architecture">
                <h2>Three-Layer Architecture</h2>
                <table>
                  <tr><th>Layer</th><th>Role</th><th>High-impact value</th></tr>
                  <tr><td><strong>NEO-PRIOR</strong></td><td>Reporting and evaluation standard</td><td>Defines the field rulebook rather than competing only as another predictor.</td></tr>
                  <tr><td><strong>NeoBench-Vax</strong></td><td>Public benchmark infrastructure</td><td>Creates reusable split manifests, audits, model cards and clean/product-assisted leaderboard tracks.</td></tr>
                  <tr><td><strong>CROSS-Neo</strong></td><td>First implementation</td><td>Demonstrates the standard with conservative internal/locked retrospective evaluation.</td></tr>
                  <tr><td><strong>Product triage</strong></td><td>Business-facing audit queue</td><td>Supports Wednesday demo while separating product-assisted utility from clean scientific claims.</td></tr>
                </table>
              </section>

              <section id="evidence">
                <h2>Evidence Snapshot</h2>
                <table>
                  <tr><th>Result family</th><th>Internal result</th><th>Interpretation</th></tr>
                  <tr><td>v0 clean C_counterfactual</td><td>HLA-stratified AUPRC 0.388; top10 precision 0.400</td><td>Best clean locked v0 baseline, but not enough for broad external claims.</td></tr>
                  <tr><td>v0 fixed late fusion</td><td>HLA-stratified AUPRC 0.473; top10 precision 0.500</td><td>Positive complementarity signal, but weights must remain fold-safe and not test-tuned.</td></tr>
                  <tr><td>Product-assisted internal demo</td><td>Strict set n=89; positives=21; product score AUPRC 0.582; top10 precision 0.700</td><td>Useful business triage signal, not clean scientific superiority.</td></tr>
                  <tr><td>Source-heldout stress</td><td>NEPdb/TESLA top-k often collapses</td><td>Central reason for OOD, abstention and source-stress reporting.</td></tr>
                </table>
              </section>

              <section id="figures">
                <h2>Launch Figures</h2>
                <div class="figgrid">
                  {asset_cards}
                </div>
              </section>

              <section id="papers">
                <h2>Two-Paper Strategy</h2>
                <div class="grid2">
                  <div class="card"><h3>Paper 1</h3><p><strong>Vaccine-ready AI for personalized cancer immunotherapy.</strong><br>Perspective/Comment route: Nature Medicine, Nature Cancer, Nature Reviews Clinical Oncology, Nature Reviews Immunology.</p></div>
                  <div class="card"><h3>Paper 2</h3><p><strong>Benchmarking vaccine-ready AI systems for neoantigen prioritization.</strong><br>Benchmark/system route: Nature Biomedical Engineering, Nature Machine Intelligence, Cell Reports Medicine, Patterns.</p></div>
                </div>
                <h3>Abstract Pack</h3>
                <pre>{abstracts}</pre>
              </section>

              <section id="product">
                <h2>Business Demo Script</h2>
                <pre>{demo}</pre>
              </section>

              <section id="claim">
                <h2>Claim Boundary</h2>
                <div class="grid3">
                  <div class="card"><h3 class="good">Can Say</h3><p>Internal locked retrospective evaluation; auditable triage; reporting standard; benchmark infrastructure; first implementation.</p></div>
                  <div class="card"><h3 class="warn">Separate</h3><p>Clean comparator track versus product-assisted track. Public predictor scores are disclosed and not used for clean superiority claims.</p></div>
                  <div class="card"><h3 class="danger">Do Not Say</h3><p>External validation, quantum advantage, universal superiority, clinical decision support approval, or definitive biological negatives.</p></div>
                </div>
              </section>

              <section id="send">
                <h2>What To Send</h2>
                <table>
                  <tr><th>Audience</th><th>Send</th><th>Purpose</th></tr>
                  <tr><td>Nature Medicine / Nature Cancer editor</td><td>Presubmission inquiry + one-slide story + NEO-PRIOR checklist</td><td>Frame as clinical AI governance for personalized vaccines.</td></tr>
                  <tr><td>Nature Biomedical Engineering / NMI editor</td><td>Benchmark README + NeoBench-Vax infrastructure figure + claim ladder</td><td>Frame as computational health system and evaluation contract.</td></tr>
                  <tr><td>Senior vaccine advisor</td><td>Consensus-style statement + advisor ask + clinical value chain</td><td>Strengthen clinical legitimacy.</td></tr>
                  <tr><td>Business partner</td><td>Wednesday demo script + product-assisted queue + landing page</td><td>Show practical triage without overclaiming.</td></tr>
                </table>
              </section>

              <section id="sources">
                <h2>Local Source Excerpts</h2>
                <h3>Apex Memo</h3>
                <pre>{apex_memo}</pre>
                <h3>Ultra Editorial Dossier</h3>
                <pre>{ultra}</pre>
                <h3>NEO-PRIOR Checklist</h3>
                <pre>{checklist}</pre>
                <p class="footer">Generated from local strategy outputs under {html.escape(str(RESULT))}. This page is a strategic dossier, not a clinical validation claim.</p>
              </section>
            </main>
          </div>
        </body>
        </html>
        """
    ).strip() + "\n"


def main() -> None:
    HUB.mkdir(parents=True, exist_ok=True)
    assets = copy_assets()
    OUT.write_text(make_html(assets), encoding="utf-8")
    print(str(OUT))
    print(f"assets={len(assets)}")


if __name__ == "__main__":
    main()
