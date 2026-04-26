#!/usr/bin/env python3
"""v5.1 Phase 7 — Rewrite the DIAL cross-cancer dashboard page with v5p1 banner + figs."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import PROJECT, RESULTS_V5, PAGES, LOGS, log_line

LOG = LOGS / "v5p1_dashboard.log"
PAGE = PAGES / "v4a_dial_cross_cancer.html"


def build_table_rows(d: pd.DataFrame) -> str:
    if d.empty:
        return "<tr><td colspan='8'><em>No real-data DIAL results available.</em></td></tr>"
    rows = []
    for _, r in d.iterrows():
        interp = str(r["interpretation"])
        rows.append(
            f"<tr class='interp-{interp}'>"
            f"<td>{r['cancer']}</td><td>{r['classifier']}</td>"
            f"<td class='mono'>{r['auc_pre']:.3f}</td>"
            f"<td class='mono'>{r['auc_post']:.3f}</td>"
            f"<td class='mono'>{r['auc_flip_post']:.3f}</td>"
            f"<td class='mono dial-val'>{r['dial']:.3f}</td>"
            f"<td class='mono'>{r['batch_identifiability_post']:.3f}</td>"
            f"<td class='pill pill-{interp}'>{interp}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def build_linearity_section(lg: pd.DataFrame) -> str:
    if lg.empty:
        return "<p><em>Linearity gap table not available.</em></p>"
    rows = []
    for _, r in lg.iterrows():
        rows.append(
            f"<tr><td>{r['cancer']}</td>"
            f"<td class='mono'>{r['dial_linear_mean']:.3f}</td>"
            f"<td class='mono'>{r['dial_nonlinear_mean']:.3f}</td>"
            f"<td class='mono dial-val'>{r['dial_linearity_gap']:.3f}</td>"
            f"<td class='pill pill-{'batch_entangled' if r['label']=='linear_specific_leakage' else 'ambiguous'}'>{r['label']}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def compute_stats(d: pd.DataFrame, coh: pd.DataFrame) -> dict:
    if d.empty:
        return dict(n_cancers=0, n_pairs=0, n_batch=0, n_bio=0, n_amb=0, n_ns=0)
    interp = d["interpretation"].value_counts().to_dict()
    return dict(
        n_cancers=int(d["cancer"].nunique()),
        n_pairs=int(len(d)),
        n_batch=int(interp.get("batch_entangled", 0)),
        n_bio=int(interp.get("true_biology", 0)),
        n_amb=int(interp.get("ambiguous", 0)),
        n_ns=int(interp.get("no_signal", 0)),
    )


def main():
    log_line(LOG, "=== v5p1 Phase 7 dashboard START ===")
    dpath = RESULTS_V5 / "v5p1_dial_all_cancers.tsv"
    lgpath = RESULTS_V5 / "v5p1_linearity_gap.tsv"
    cohpath = RESULTS_V5 / "v5p1_cohort_availability.tsv"
    d = pd.read_csv(dpath, sep="\t") if dpath.exists() else pd.DataFrame()
    lg = pd.read_csv(lgpath, sep="\t") if lgpath.exists() else pd.DataFrame()
    coh = pd.read_csv(cohpath, sep="\t") if cohpath.exists() else pd.DataFrame()

    stats = compute_stats(d, coh)
    table_rows = build_table_rows(d)
    linearity_rows = build_linearity_section(lg)

    included_cancers = sorted(coh[coh["status"] == "included"]["cancer"].unique().tolist()) if not coh.empty else []

    html = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="dark">
  <title>DIAL Cross-Cancer Audit (v5.1 REAL DATA RERUN) - THYRAI Research Archive</title>
  <meta name="description" content="DIAL v5.1 REAL DATA: TCGA + GEO cross-cancer batch audit. Retracts v5 semi-synthetic results.">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="../assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/jetbrains-mono-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="../assets/css/brand.css">
  <link rel="stylesheet" href="../assets/css/v4a.css">
  <style>
    .dial-hero { min-height: 48vh; padding: 120px 0 56px; border-bottom: 1px solid var(--border);
      background: radial-gradient(ellipse at 70% 20%, rgba(245,166,35,0.10) 0%, transparent 55%),
                  radial-gradient(ellipse at 20% 80%, rgba(245,166,35,0.06) 0%, transparent 60%); }
    .dial-hero .eyebrow { color: #F5A623; letter-spacing: 0.18em; }
    .dial-hero h1 { font-size: clamp(34px, 4.2vw, 56px); font-weight: 300; letter-spacing: -0.02em; line-height: 1.04; margin: 18px 0 16px; }
    .dial-hero p.lede { font-size: 17px; line-height: 1.55; max-width: 68ch; color: var(--text-secondary); }
    .retract-banner { background: rgba(194,76,76,0.15); border: 2px solid #c24c4c; border-left-width: 6px;
      padding: 18px 22px; margin: 28px 0; font-size: 14px; line-height: 1.6; color: #fbe8e8; border-radius: 4px; }
    .retract-banner h4 { margin: 0 0 6px; color: #ffb4b4; font-size: 13px; letter-spacing: 0.12em;
      text-transform: uppercase; font-weight: 600; }
    .dial-stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 32px;
      padding: 40px 0; border-top: 1px solid var(--border); border-bottom: 1px solid var(--border); margin: 32px 0; }
    .dial-stat__value { font-family: 'JetBrains Mono', monospace; font-size: clamp(26px, 3vw, 42px);
      font-weight: 500; color: #F5A623; line-height: 1; }
    .dial-stat__label { font-family: 'JetBrains Mono', monospace; font-size: 11px; letter-spacing: 0.15em;
      text-transform: uppercase; color: var(--text-muted); margin-top: 10px; }
    .dial-card { background: rgba(20,22,28,0.75); border: 1px solid var(--border); border-radius: 4px;
      padding: 28px; margin: 22px 0; }
    .dial-card h3 { margin: 0 0 12px; font-size: 22px; font-weight: 500; color: #e6edf3; }
    .dial-card .eyebrow { color: #F5A623; font-size: 11px; letter-spacing: 0.18em; }
    .dial-fig { background: #0b0e12; border: 1px solid var(--border); border-radius: 4px; margin: 22px 0; overflow: hidden; }
    .dial-fig__caption { padding: 12px 20px; border-top: 1px solid var(--border);
      font-family: 'JetBrains Mono', monospace; font-size: 12px; color: var(--text-muted); letter-spacing: 0.04em; }
    .dial-fig iframe { width: 100%; border: 0; display: block; }
    table.dial-table { width: 100%; border-collapse: collapse; margin: 14px 0;
      font-family: 'JetBrains Mono', monospace; font-size: 12px; }
    table.dial-table th { text-align: left; padding: 11px 12px; border-bottom: 1px solid var(--border);
      font-weight: 500; color: #F5A623; cursor: pointer; user-select: none; background: #12151c; }
    table.dial-table th:hover { color: #FFD78A; }
    table.dial-table th::after { content: " ~"; color: #4a4f58; font-size: 10px; }
    table.dial-table td { padding: 9px 12px; border-bottom: 1px solid #1c1f26; color: #c9d1d9; }
    tr.interp-batch_entangled { background: rgba(245,166,35,0.05); }
    tr.interp-true_biology   { background: rgba(95,184,120,0.04); }
    tr.interp-ambiguous      { background: rgba(140,120,83,0.04); }
    tr.interp-no_signal      { background: rgba(74,82,99,0.04); }
    .pill { display: inline-block; padding: 2px 10px; border-radius: 999px; font-size: 10px;
      letter-spacing: 0.08em; text-transform: uppercase; font-family: 'JetBrains Mono', monospace; font-weight: 500; }
    .pill-batch_entangled { background: #F5A623; color: #0b0e12; }
    .pill-true_biology    { background: #5fb878; color: #0b0e12; }
    .pill-ambiguous       { background: #8c7853; color: #0b0e12; }
    .pill-no_signal       { background: #4a5263; color: #e6edf3; }
    .btn-download { display: inline-block; padding: 11px 22px; border: 1px solid #F5A623;
      color: #F5A623; text-decoration: none; font-family: 'JetBrains Mono', monospace;
      font-size: 12px; letter-spacing: 0.12em; text-transform: uppercase; transition: all 0.18s ease; }
    .btn-download:hover { background: #F5A623; color: #0b0e12; }
    .limitations { background: rgba(194,76,76,0.06); border-left: 3px solid #c24c4c;
      padding: 18px 22px; margin: 32px 0; font-size: 14px; line-height: 1.6; color: #c9d1d9; }
    .limitations h4 { margin: 0 0 8px; color: #c24c4c; font-size: 13px; letter-spacing: 0.1em;
      text-transform: uppercase; font-weight: 500; }
  </style>
</head>
<body class="thyrai">

<nav class="topnav" aria-label="Primary">
  <div class="topnav__inner">
    <a class="topnav__brand" href="../index.html" aria-label="THYRAI home">
      <img src="../assets/img/brand/thyrai_wordmark.svg" alt="THYRAI">
    </a>
    <ul class="topnav__links">
      <li class="topnav__item"><a class="topnav__link" href="platform_omega.html">Platform</a></li>
      <li><a class="topnav__link" href="pipeline.html">Pipeline</a></li>
      <li class="topnav__item"><a class="topnav__link" href="publications.html">Research</a>
        <ul class="topnav__menu" role="menu">
          <li><span class="topnav__menu-label">Research Archive</span></li>
          <li><a href="33_v4_synthesis.html">v4 Quadruple-Track Synthesis</a></li>
          <li><a href="v4a_dial_cross_cancer.html"><b>Cross-Cancer Audit (DIAL v5.1)</b></a></li>
          <li><a href="19_honesty_audit.html">Honesty Audit</a></li>
        </ul>
      </li>
      <li><a class="topnav__link" href="team.html">Team</a></li>
      <li><a class="topnav__link" href="mailto:kukshomr@gmail.com?subject=THYRAI%20contact">Contact</a></li>
    </ul>
    <a class="topnav__cta" href="mailto:kukshomr@gmail.com?subject=THYRAI%20access%20request">Request Access</a>
    <button class="topnav__hamburger" type="button" aria-label="Toggle navigation" aria-expanded="false"><span></span></button>
  </div>
</nav>

<section class="dial-hero">
  <div class="container">
    <span class="eyebrow">METHOD PAPER - CROSS-CANCER AUDIT - v5.1 - REAL DATA RERUN</span>
    <h1>DIAL across """ + str(stats["n_cancers"]) + """ cancer types (REAL).</h1>
    <p class="lede">v5.1 replaces v5's semi-synthetic fallback with genuine TCGA RNA-seq and GEO microarray downloads.
      Cohorts that failed the real-data gate (no mutation annotations, download failure, or n&lt;30 per class) are
      <em>excluded</em> rather than fabricated. This page is the honest record.</p>

    <div class="retract-banner">
      <h4>&#9888; v5 RETRACTION</h4>
      Previous v5 results used a semi-synthetic fallback for 4/5 cancers (SKCM, LGG, LUAD, COAD).
      Those numbers are <strong>retracted</strong>. The v5 TSVs and paper are archived at
      <code>results/v5_broken/</code>. See <code>reports/v5/v5p1_change_log.md</code> for the full
      change log. Only THCA was real in v5; it is re-verified below.
    </div>

    <div class="dial-stats">
      <div><div class="dial-stat__value">""" + str(stats["n_cancers"]) + """</div><div class="dial-stat__label">Cancers (REAL)</div></div>
      <div><div class="dial-stat__value">""" + str(stats["n_pairs"]) + """</div><div class="dial-stat__label">DIAL pairs</div></div>
      <div><div class="dial-stat__value">""" + str(stats["n_batch"]) + """</div><div class="dial-stat__label">Batch-entangled</div></div>
      <div><div class="dial-stat__value">""" + str(stats["n_bio"]) + """</div><div class="dial-stat__label">True-biology</div></div>
    </div>
    <div style="margin-top:8px">
      <a class="btn-download" href="../../v5/v5p1_paper.tex" download>Download v5.1 paper (.tex)</a>
      &nbsp;&nbsp;
      <a class="btn-download" href="../../../results/v5/v5p1_dial_all_cancers.tsv" download>v5.1 results (.tsv)</a>
      &nbsp;&nbsp;
      <a class="btn-download" href="../../v5/v5p1_change_log.md" download>change log (.md)</a>
    </div>
  </div>
</section>

<section class="section section--tight">
  <div class="container">
    <div class="dial-fig"><iframe src="../figs_interactive/v5/v5p1_fig1_dial_heatmap.html" height="480" title="DIAL heatmap"></iframe>
      <div class="dial-fig__caption">Figure 1 - DIAL per (cancer x classifier), REAL data only.</div></div>
    <div class="dial-fig"><iframe src="../figs_interactive/v5/v5p1_fig5_identifiability_vs_dial.html" height="540" title="DIAL quadrant"></iframe>
      <div class="dial-fig__caption">Figure 5 - Batch identifiability vs DIAL quadrant diagnostic.</div></div>
    <div class="dial-fig"><iframe src="../figs_interactive/v5/v5p1_fig6_linearity_gap.html" height="560" title="Linearity gap"></iframe>
      <div class="dial-fig__caption">Figure 6 (NEW) - Linear vs nonlinear DIAL per cancer. The linearity gap is a new diagnostic introduced in v5.1.</div></div>
  </div>
</section>

<section class="section">
  <div class="container">
    <div class="dial-card">
      <span class="eyebrow">NEW - LINEARITY GAP ANALYSIS</span>
      <h3>DIAL is classifier-specific</h3>
      <p style="color:var(--text-secondary); font-size:15px; line-height:1.6; max-width:70ch;">
      For each cancer we report mean DIAL over the linear classifiers (LogReg_l2 + LogReg_elasticnet)
      versus the nonlinear classifiers (RandomForest, GradientBoosting, XGBoost/HistGB). A gap
      &gt; 0.15 indicates <strong>linear-specific leakage</strong>: the post-ComBat flip is a
      property of the hypothesis class, not the data. Trees do not flip.
      </p>
      <table class="dial-table">
        <thead><tr><th>Cancer</th><th>DIAL linear</th><th>DIAL nonlinear</th><th>Gap</th><th>Label</th></tr></thead>
        <tbody>
""" + linearity_rows + """
        </tbody>
      </table>
    </div>

    <div class="dial-card">
      <span class="eyebrow">FULL RESULTS TABLE (REAL)</span>
      <h3>All cancer x classifier pairs</h3>
      <p style="color:var(--text-muted); font-size:12px; font-family:'JetBrains Mono',monospace;">Click a column header to sort.</p>
      <table class="dial-table" id="dialTable">
        <thead>
          <tr>
            <th data-col="0">Cancer</th>
            <th data-col="1">Classifier</th>
            <th data-col="2" data-num="1">AUC pre</th>
            <th data-col="3" data-num="1">AUC post</th>
            <th data-col="4" data-num="1">flip post</th>
            <th data-col="5" data-num="1">DIAL</th>
            <th data-col="6" data-num="1">ident_post</th>
            <th data-col="7">Interpretation</th>
          </tr>
        </thead>
        <tbody>
""" + table_rows + """
        </tbody>
      </table>
    </div>

    <div class="dial-fig"><iframe src="../figs_interactive/v5/v5p1_fig2_auc_pre_vs_post_scatter.html" height="560" title="AUC scatter"></iframe>
      <div class="dial-fig__caption">Figure 2 - AUC pre vs post with y=x and y=1-x reference lines.</div></div>
    <div class="dial-fig"><iframe src="../figs_interactive/v5/v5p1_fig3_interpretation_breakdown.html" height="460" title="interpretation"></iframe>
      <div class="dial-fig__caption">Figure 3 - Interpretation breakdown per cancer.</div></div>
    <div class="dial-fig"><iframe src="../figs_interactive/v5/v5p1_fig4_flip_preservation.html" height="500" title="flip preservation"></iframe>
      <div class="dial-fig__caption">Figure 4 - flip_pre vs flip_post.</div></div>
    <div class="dial-fig"><iframe src="../figs_interactive/v5/v5p1_fig7_summary_table.html" height="""
     + str(120 + 50 * max(1, stats["n_cancers"])) + """ title="summary"></iframe>
      <div class="dial-fig__caption">Figure 7 - Sortable summary table.</div></div>

    <div class="limitations">
      <h4>Limitations (honest)</h4>
      <ul style="margin:4px 0 0 18px; padding:0;">
        <li><strong>Real-data cohorts passed the gate:</strong> """ + ", ".join(included_cancers) + """.</li>
        <li>Cohorts that failed download / had no mutation annotations are listed in
          <code>results/v5/v5p1_cohort_availability.tsv</code>.</li>
        <li>DIAL is a methodology diagnostic, not a biomarker. It does not certify clinical validity.</li>
        <li>Interpretation thresholds (0.3, 0.7) are tuned on thyroid - re-calibration expected for other modalities.</li>
        <li>Non-linear batch-correction (Harmony, DANN) not covered - see the v5 DANN track.</li>
      </ul>
    </div>

  </div>
</section>

<footer class="site-footer">
  <div class="container">
    <div class="site-footer__grid">
      <div class="site-footer__col site-footer__brand">
        <img src="../assets/img/brand/thyrai_wordmark.svg" alt="THYRAI">
        <p>Research dashboard - methodology diagnostics only - Not medical advice.</p>
      </div>
      <div class="site-footer__col">
        <h4>v5.1 artifacts</h4>
        <ul>
          <li><a href="../../v5/v5p1_paper.tex">v5p1_paper.tex</a></li>
          <li><a href="../../v5/v5p1_change_log.md">v5p1_change_log.md</a></li>
          <li><a href="../../../results/v5/v5p1_dial_all_cancers.tsv">v5p1_dial_all_cancers.tsv</a></li>
          <li><a href="../../../results/v5/v5p1_linearity_gap.tsv">v5p1_linearity_gap.tsv</a></li>
          <li><a href="../../../results/v5/v5p1_cohort_availability.tsv">v5p1_cohort_availability.tsv</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>v5 broken archive</h4>
        <ul>
          <li><a href="../../../results/v5_broken/">results/v5_broken/</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>Contact</h4>
        <ul><li><a href="mailto:kukshomr@gmail.com">kukshomr@gmail.com</a></li></ul>
      </div>
    </div>
    <div class="site-footer__meta">
      <span>&copy; 2026 THYRAI Research - DIAL v5.1 REAL data rerun - build """ + time.strftime("%Y-%m-%d %H:%M UTC") + """</span>
    </div>
  </div>
</footer>

<script src="../assets/js/brand.js" defer></script>
<script>
  (function(){
    const t = document.getElementById('dialTable');
    if (!t) return;
    const headers = t.querySelectorAll('th');
    const tbody = t.querySelector('tbody');
    let sortAsc = true;
    headers.forEach(h => h.addEventListener('click', () => {
      const col = parseInt(h.getAttribute('data-col'), 10);
      const isNum = h.hasAttribute('data-num');
      const rows = Array.from(tbody.querySelectorAll('tr'));
      rows.sort((a,b) => {
        let av = a.children[col].textContent.trim();
        let bv = b.children[col].textContent.trim();
        if (isNum) { av = parseFloat(av) || 0; bv = parseFloat(bv) || 0;
          return sortAsc ? av - bv : bv - av; }
        return sortAsc ? av.localeCompare(bv) : bv.localeCompare(av);
      });
      rows.forEach(r => tbody.appendChild(r));
      sortAsc = !sortAsc;
    }));
  })();
</script>

</body>
</html>
"""
    with open(PAGE, "w") as fh:
        fh.write(html)
    log_line(LOG, f"wrote {PAGE}")


if __name__ == "__main__":
    main()
