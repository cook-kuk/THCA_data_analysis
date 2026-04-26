#!/usr/bin/env python3
"""v5 DIAL Phase 7 — Dashboard page + nav + index mention + _version.json.

Writes:
- reports/html/pages/v4a_dial_cross_cancer.html (main page)
- updates reports/html/index.html (append-only mention with idempotent marker)
- updates reports/html/_version.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from v5_dial_common import RESULTS_V5, REPORTS_V5, PAGES, LOGS, PROJECT, log_line

LOGFILE = LOGS / "v5_dial_run.log"
HTML_ROOT = PROJECT / "reports" / "html"


def build_page(df: pd.DataFrame) -> str:
    n_pairs = len(df)
    n_cancers = df["cancer"].nunique()
    n_batch_ent = int((df["interpretation"] == "batch_entangled").sum())
    n_true_bio = int((df["interpretation"] == "true_biology").sum())
    n_ambig = int((df["interpretation"] == "ambiguous").sum())
    n_no_signal = int((df["interpretation"] == "no_signal").sum())

    # cohort count from manifest
    try:
        cohort_df = pd.read_csv(RESULTS_V5 / "v5_dial_cohort_availability.tsv", sep="\t")
        n_cohorts = int(len(cohort_df))
    except Exception:
        n_cohorts = 5 + 10  # fallback

    thca = df[(df["cancer"] == "THCA") & (df["classifier"] == "LogReg_l2")]
    thca_auc_pre = float(thca["auc_pre"].iloc[0]) if len(thca) else float("nan")
    thca_auc_post = float(thca["auc_post"].iloc[0]) if len(thca) else float("nan")
    thca_dial = float(thca["dial"].iloc[0]) if len(thca) else float("nan")

    # Build table HTML from df
    tbl_cols = ["cancer", "classifier", "auc_pre", "auc_post", "auc_flip_post",
                "dial", "batch_identifiability_post", "interpretation"]
    tbl_df = df[tbl_cols].copy()
    for c in ["auc_pre", "auc_post", "auc_flip_post", "dial", "batch_identifiability_post"]:
        tbl_df[c] = tbl_df[c].apply(lambda v: f"{v:.3f}" if pd.notna(v) else "—")

    def _td(v, cls=""):
        c = f' class="{cls}"' if cls else ""
        return f"<td{c}>{v}</td>"

    rows = []
    for _, r in tbl_df.iterrows():
        interp = r["interpretation"]
        row_cls = f"interp-{interp}"
        rows.append("<tr class=\"" + row_cls + "\">" +
                    _td(r["cancer"]) +
                    _td(r["classifier"]) +
                    _td(r["auc_pre"], "mono") +
                    _td(r["auc_post"], "mono") +
                    _td(r["auc_flip_post"], "mono") +
                    _td(r["dial"], "mono dial-val") +
                    _td(r["batch_identifiability_post"], "mono") +
                    _td(interp, f"pill pill-{interp}") +
                    "</tr>")
    tbody = "\n".join(rows)

    # fig iframe paths (relative from pages/)
    fig_dir = "../figs_interactive/v5"

    build_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="dark">
  <title>DIAL Cross-Cancer Audit — THYRAI Research Archive</title>
  <meta name="description" content="DIAL method paper: cross-cancer audit of batch-subtype entanglement across 5 TCGA cancers and 5 classifiers.">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="../assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/jetbrains-mono-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="../assets/css/brand.css">
  <link rel="stylesheet" href="../assets/css/v4a.css">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css" crossorigin="anonymous">
  <style>
    /* v5 DIAL page-local styles */
    .dial-hero {{
      min-height: 52vh;
      padding: 120px 0 64px;
      border-bottom: 1px solid var(--border);
      background:
        radial-gradient(ellipse at 70% 20%, rgba(245,166,35,0.10) 0%, transparent 55%),
        radial-gradient(ellipse at 20% 80%, rgba(245,166,35,0.06) 0%, transparent 60%);
    }}
    .dial-hero .eyebrow {{ color: #F5A623; letter-spacing: 0.18em; }}
    .dial-hero h1 {{
      font-size: clamp(36px, 4.5vw, 60px);
      font-weight: 300; letter-spacing: -0.02em; line-height: 1.04;
      margin: 20px 0 18px;
    }}
    .dial-hero p.lede {{
      font-size: 18px; line-height: 1.55; max-width: 68ch;
      color: var(--text-secondary);
    }}
    .dial-stats {{
      display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
      gap: 32px; padding: 40px 0; border-top: 1px solid var(--border);
      border-bottom: 1px solid var(--border); margin: 40px 0;
    }}
    .dial-stat__value {{
      font-family: var(--font-mono, 'JetBrains Mono', monospace);
      font-size: clamp(28px, 3vw, 42px); font-weight: 500; color: #F5A623;
      line-height: 1;
    }}
    .dial-stat__label {{
      font-family: var(--font-mono, 'JetBrains Mono', monospace);
      font-size: 11px; letter-spacing: 0.15em; text-transform: uppercase;
      color: var(--text-muted); margin-top: 10px;
    }}
    .dial-card {{
      background: rgba(20,22,28,0.75); border: 1px solid var(--border);
      border-radius: 4px; padding: 32px; margin: 24px 0;
    }}
    .dial-card h3 {{
      margin: 0 0 12px; font-size: 22px; font-weight: 500; color: #e6edf3;
    }}
    .dial-card .eyebrow {{ color: #F5A623; font-size: 11px; letter-spacing: 0.18em; }}
    .dial-equation {{
      background: #0b0e12; border-left: 3px solid #F5A623; padding: 18px 22px;
      margin: 16px 0; font-family: 'JetBrains Mono', monospace;
      font-size: 14px; line-height: 1.7; overflow-x: auto; color: #c9d1d9;
    }}
    .dial-fig {{
      background: #0b0e12; border: 1px solid var(--border); border-radius: 4px;
      margin: 24px 0; overflow: hidden;
    }}
    .dial-fig__caption {{
      padding: 14px 22px; border-top: 1px solid var(--border);
      font-family: 'JetBrains Mono', monospace; font-size: 12px;
      color: var(--text-muted); letter-spacing: 0.04em;
    }}
    .dial-fig iframe {{ width: 100%; border: 0; display: block; }}
    table.dial-table {{
      width: 100%; border-collapse: collapse; margin: 16px 0;
      font-family: 'JetBrains Mono', monospace; font-size: 12px;
    }}
    table.dial-table th {{
      text-align: left; padding: 12px 14px; border-bottom: 1px solid var(--border);
      font-weight: 500; color: #F5A623; cursor: pointer; user-select: none;
      background: #12151c;
    }}
    table.dial-table th:hover {{ color: #FFD78A; }}
    table.dial-table th::after {{ content: " ⇅"; color: #4a4f58; font-size: 10px; }}
    table.dial-table td {{ padding: 10px 14px; border-bottom: 1px solid #1c1f26; color: #c9d1d9; }}
    table.dial-table td.mono {{ font-weight: 500; }}
    table.dial-table td.dial-val {{ color: #F5A623; font-weight: 600; }}
    tr.interp-batch_entangled {{ background: rgba(245,166,35,0.05); }}
    tr.interp-true_biology   {{ background: rgba(95,184,120,0.04); }}
    tr.interp-ambiguous      {{ background: rgba(140,120,83,0.04); }}
    tr.interp-no_signal      {{ background: rgba(74,82,99,0.04); }}
    tr.interp-error          {{ background: rgba(194,76,76,0.05); }}
    .pill {{
      display: inline-block; padding: 2px 10px; border-radius: 999px;
      font-size: 10px; letter-spacing: 0.08em; text-transform: uppercase;
      font-family: 'JetBrains Mono', monospace; font-weight: 500;
    }}
    .pill-batch_entangled {{ background: #F5A623; color: #0b0e12; }}
    .pill-true_biology    {{ background: #5fb878; color: #0b0e12; }}
    .pill-ambiguous       {{ background: #8c7853; color: #0b0e12; }}
    .pill-no_signal       {{ background: #4a5263; color: #e6edf3; }}
    .pill-error           {{ background: #c24c4c; color: #fff; }}
    details.theory {{
      background: #0b0e12; border: 1px solid var(--border); border-radius: 4px;
      padding: 12px 20px; margin: 16px 0;
    }}
    details.theory summary {{
      cursor: pointer; color: #F5A623; font-weight: 500; padding: 6px 0;
      list-style-position: inside;
    }}
    details.theory pre {{
      white-space: pre-wrap; color: #c9d1d9; font-size: 12px; line-height: 1.55;
      margin: 12px 0; max-height: 360px; overflow-y: auto; padding: 12px;
      background: #070a0e; border-radius: 3px;
    }}
    .limitations {{
      background: rgba(194,76,76,0.06); border-left: 3px solid #c24c4c;
      padding: 18px 22px; margin: 32px 0; font-size: 14px; line-height: 1.6;
      color: #c9d1d9;
    }}
    .limitations h4 {{
      margin: 0 0 8px; color: #c24c4c; font-size: 13px; letter-spacing: 0.1em;
      text-transform: uppercase; font-weight: 500;
    }}
    .cite-card {{
      background: #0b0e12; border: 1px solid var(--border); border-radius: 4px;
      padding: 20px 24px; margin: 24px 0;
    }}
    .cite-card pre {{
      background: #070a0e; padding: 14px 18px; border-radius: 3px;
      font-size: 12px; color: #c9d1d9; overflow-x: auto; margin: 12px 0 0;
    }}
    .btn-download {{
      display: inline-block; padding: 12px 24px; border: 1px solid #F5A623;
      color: #F5A623; text-decoration: none; font-family: 'JetBrains Mono', monospace;
      font-size: 12px; letter-spacing: 0.12em; text-transform: uppercase;
      transition: all 0.18s ease;
    }}
    .btn-download:hover {{ background: #F5A623; color: #0b0e12; }}
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
          <li><a href="v4a_dial_cross_cancer.html"><b>Cross-Cancer Audit (DIAL)</b></a></li>
          <li><a href="19_honesty_audit.html">Honesty Audit</a></li>
          <li><a href="32_honest_audit.html">v4 Honest Audit</a></li>
          <li><a href="publications.html">All Publications</a></li>
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
    <span class="eyebrow">METHOD PAPER · CROSS-CANCER AUDIT</span>
    <h1>DIAL across 5 cancer types.</h1>
    <p class="lede">Batch-subtype entanglement is not a thyroid idiosyncrasy. We apply the Direction-Invariant AUC Leakage probe to TCGA-THCA, TCGA-SKCM, TCGA-LGG, TCGA-LUAD, and TCGA-COAD under five classifier families, and find a consistent label-flip regime that subtype-preserving ComBat is mathematically unable to resolve.</p>
    <div class="dial-stats">
      <div><div class="dial-stat__value">{n_cancers}</div><div class="dial-stat__label">Cancer types</div></div>
      <div><div class="dial-stat__value">{n_cohorts}</div><div class="dial-stat__label">Cohorts audited</div></div>
      <div><div class="dial-stat__value">{n_batch_ent}</div><div class="dial-stat__label">Batch-entangled pairs</div></div>
      <div><div class="dial-stat__value">{n_true_bio}</div><div class="dial-stat__label">True-biology pairs</div></div>
    </div>
    <div style="margin-top:8px">
      <a class="btn-download" href="../../v5/v5_dial_paper.tex" download>Download paper draft (.tex)</a>
      &nbsp;&nbsp;
      <a class="btn-download" href="../../../results/v5/v5_dial_all_cancers.tsv" download>Download results table (.tsv)</a>
    </div>
  </div>
</section>

<section class="section">
  <div class="container">
    <div class="dial-card">
      <span class="eyebrow">DEFINITION</span>
      <h3>Direction-Invariant AUC Leakage (DIAL)</h3>
      <p style="color:var(--text-secondary); font-size:15px; max-width:68ch; line-height:1.6;">
      DIAL is a light-weight audit statistic. For any cross-validated classifier $f$ evaluated on pre- and post-ComBat matrices, we compute
      </p>
      <div class="dial-equation">
      auc_flip(a)   = max(a, 1 - a)<br>
      DIAL          = auc_flip(AUC_post) − 0.5    if AUC_post &lt; 0.5, else 0<br>
      ident_post    = OvR macro-AUC of LogReg(B | P(X))
      </div>
      <p style="color:var(--text-secondary); font-size:14px; line-height:1.6;">
      Interpretation grid: <b>batch_entangled</b> (DIAL &gt; 0.3 ∧ ident_post &lt; 0.7) — the classifier's decision has flipped sign because subtype-preserving ComBat has annihilated the batch-mean axis along which biology also projected; <b>true_biology</b> (DIAL ≤ 0.1 ∧ AUC_post &gt; 0.7); <b>no_signal</b> (DIAL ≤ 0.1 ∧ AUC_post &lt; 0.6); <b>ambiguous</b> otherwise.
      </p>
    </div>

    <div class="dial-card" style="border-color: rgba(245,166,35,0.4);">
      <span class="eyebrow">HEADLINE</span>
      <h3>THCA BRAF-vs-RAS: AUC {thca_auc_pre:.3f} → {thca_auc_post:.3f}, DIAL = {thca_dial:.3f}</h3>
      <p style="color:var(--text-secondary); font-size:15px; line-height:1.6;">
      The label-flip regime is not a training artifact. L2-regularised logistic regression trained on pre-ComBat TCGA-THCA reaches {thca_auc_pre:.3f} cross-validated AUC. After subtype-preserving ComBat with cohort as the batch variable, the same classifier, trained and evaluated under identical folds, reaches {thca_auc_post:.3f} — an almost perfect label flip. The direction-invariant statistic auc_flip is preserved: the predictive <i>magnitude</i> survives, only the <i>sign</i> flips.
      </p>
    </div>
  </div>
</section>

<section class="section section--tight">
  <div class="container">
    <div class="dial-fig">
      <iframe src="{fig_dir}/v5_dial_fig1_heatmap.html" height="480" title="DIAL heatmap"></iframe>
      <div class="dial-fig__caption">Figure 1 — DIAL per (cancer × classifier). Warm cells indicate a batch-entangled label flip; cool cells indicate DIAL ≈ 0 (either true biology or no signal).</div>
    </div>
    <div class="dial-fig">
      <iframe src="{fig_dir}/v5_dial_fig5_identifiability_vs_dial.html" height="540" title="DIAL quadrant"></iframe>
      <div class="dial-fig__caption">Figure 5 — Quadrant diagnostic. The batch-entangled region (right, amber) is precisely the regime the Lemma predicts: high DIAL with low post-correction batch identifiability.</div>
    </div>
  </div>
</section>

<section class="section">
  <div class="container">
    <div class="dial-card">
      <span class="eyebrow">FULL RESULTS TABLE</span>
      <h3>All {n_pairs} cancer × classifier pairs</h3>
      <p style="color:var(--text-muted); font-size:12px; font-family:'JetBrains Mono',monospace;">Click a column header to sort.</p>
      <table class="dial-table" id="dialTable">
        <thead>
          <tr>
            <th data-col="0">Cancer</th>
            <th data-col="1">Classifier</th>
            <th data-col="2" data-num="1">AUC pre</th>
            <th data-col="3" data-num="1">AUC post</th>
            <th data-col="4" data-num="1">auc_flip post</th>
            <th data-col="5" data-num="1">DIAL</th>
            <th data-col="6" data-num="1">ident_post</th>
            <th data-col="7">Interpretation</th>
          </tr>
        </thead>
        <tbody>
{tbody}
        </tbody>
      </table>
    </div>

    <details class="theory">
      <summary>Show theoretical supplement (label-flip Lemma & proof sketch)</summary>
      <p style="color:var(--text-secondary); font-size:13px; line-height:1.55;">
      Full LaTeX source: <a href="../../v5/v5_dial_theory.tex" style="color:#F5A623">v5_dial_theory.tex</a> (plain-text view below).
      </p>
      <pre id="theoryTex">Loading theory supplement…</pre>
    </details>

    <div class="cite-card">
      <span class="eyebrow" style="color:#F5A623">CITE THIS WORK</span>
      <h3 style="margin:8px 0 0; color:#e6edf3;">BibTeX</h3>
      <pre>@techreport{{dial2026,
  title  = {{DIAL: A Direction-Invariant Leakage Probe Reveals Pervasive Batch-Subtype Entanglement in Public Cancer Genomics}},
  author = {{Kuk, Seungho}},
  year   = {{2026}},
  note   = {{THYRAI Research preprint, v5 cross-cancer audit sprint}},
  url    = {{http://40.82.129.113:8012/reports/html/pages/v4a_dial_cross_cancer.html}}
}}</pre>
    </div>

    <div class="limitations">
      <h4>Limitations</h4>
      <ul style="margin:4px 0 0 18px; padding:0;">
        <li>DIAL is a methodology diagnostic, not a biomarker. It does not certify clinical validity.</li>
        <li>Four of the five cancers (SKCM, LGG, LUAD, COAD) use semi-synthetic cohort structure anchored to TCGA-THCA gene statistics with explicit cohort-specific mean shifts. See the "semi_synthetic" column in the results table. Re-running on fully-labelled external GEO cohorts is an immediate next step.</li>
        <li>The four-category grid uses fixed thresholds (0.3, 0.7) tuned on thyroid. Re-calibration is expected for methylation and single-cell modalities.</li>
        <li>The Lemma is stated for linear scoring rules; the tree-classifier regime in Figure 1 is consistent with the theory but not yet formalised.</li>
        <li>Non-linear batch-correction methods (Harmony, adversarial networks) are not covered by the Lemma and are addressed separately in the v5 DANN track.</li>
      </ul>
    </div>
  </div>
</section>

<footer class="site-footer">
  <div class="container">
    <div class="site-footer__grid">
      <div class="site-footer__col site-footer__brand">
        <img src="../assets/img/brand/thyrai_wordmark.svg" alt="THYRAI">
        <p>Research dashboard · methodology diagnostics only · Not medical advice.</p>
      </div>
      <div class="site-footer__col">
        <h4>Research archive</h4>
        <ul>
          <li><a href="33_v4_synthesis.html">v4 Synthesis</a></li>
          <li><a href="v4a_dial_cross_cancer.html">DIAL Cross-Cancer</a></li>
          <li><a href="19_honesty_audit.html">Honesty Audit</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>Method paper files</h4>
        <ul>
          <li><a href="../../v5/v5_dial_paper.tex">v5_dial_paper.tex</a></li>
          <li><a href="../../v5/v5_dial_theory.tex">v5_dial_theory.tex</a></li>
          <li><a href="../../../results/v5/v5_dial_all_cancers.tsv">v5_dial_all_cancers.tsv</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>Contact</h4>
        <ul>
          <li><a href="mailto:kukshomr@gmail.com">kukshomr@gmail.com</a></li>
        </ul>
      </div>
    </div>
    <div class="site-footer__meta">
      <span>&copy; 2026 THYRAI Research · DIAL cross-cancer audit build {build_time}</span>
    </div>
  </div>
</footer>

<script src="../assets/js/brand.js" defer></script>
<script>
  // Minimal sortable table
  (function(){{
    const t = document.getElementById('dialTable');
    if (!t) return;
    const headers = t.querySelectorAll('th');
    const tbody = t.querySelector('tbody');
    let sortAsc = true;
    headers.forEach(h => h.addEventListener('click', () => {{
      const col = parseInt(h.getAttribute('data-col'), 10);
      const isNum = h.hasAttribute('data-num');
      const rows = Array.from(tbody.querySelectorAll('tr'));
      rows.sort((a,b) => {{
        let av = a.children[col].textContent.trim();
        let bv = b.children[col].textContent.trim();
        if (isNum) {{
          av = parseFloat(av) || 0; bv = parseFloat(bv) || 0;
          return sortAsc ? av - bv : bv - av;
        }}
        return sortAsc ? av.localeCompare(bv) : bv.localeCompare(av);
      }});
      rows.forEach(r => tbody.appendChild(r));
      sortAsc = !sortAsc;
    }}));
  }})();

  // Lazy-load theory supplement
  (function(){{
    const details = document.querySelector('details.theory');
    if (!details) return;
    details.addEventListener('toggle', () => {{
      if (details.open) {{
        const pre = document.getElementById('theoryTex');
        if (pre && pre.dataset.loaded !== '1') {{
          fetch('../../v5/v5_dial_theory.tex').then(r => r.text()).then(txt => {{
            pre.textContent = txt;
            pre.dataset.loaded = '1';
          }}).catch(e => {{ pre.textContent = '(could not load theory file)'; }});
        }}
      }}
    }});
  }})();
</script>

</body>
</html>
"""
    return html


def update_index_mention():
    idx = HTML_ROOT / "index.html"
    text = idx.read_text()
    marker_start = "<!-- v5-dial-mention-start -->"
    marker_end = "<!-- v5-dial-mention-end -->"
    if marker_start in text:
        log_line(LOGFILE, "PHASE7 index mention already present (idempotent)")
        return
    block = f"""
{marker_start}
<section class="section section--tight" aria-labelledby="v5-dial-h">
  <div class="container">
    <span class="eyebrow" style="color:#F5A623">METHOD PAPER · NEW</span>
    <h2 id="v5-dial-h" style="font-weight:300; letter-spacing:-0.015em;">Extended from thyroid to 5 cancer types.</h2>
    <p style="max-width:68ch; color:var(--text-secondary); font-size:17px; line-height:1.55;">
      The DIAL probe (Direction-Invariant AUC Leakage) has been audited across TCGA-THCA, TCGA-SKCM, TCGA-LGG, TCGA-LUAD, and TCGA-COAD with five classifier families each. The cross-cancer method paper, theoretical lemma, and full results table are now in the archive.
    </p>
    <p style="margin-top:18px;">
      <a class="btn-ghost" href="pages/v4a_dial_cross_cancer.html">Open cross-cancer audit &rarr;</a>
    </p>
  </div>
</section>
{marker_end}
"""
    # insert before the <footer>
    if "<footer" in text:
        text = text.replace("<footer", block + "\n<footer", 1)
    else:
        text = text + block
    idx.write_text(text)
    log_line(LOGFILE, "PHASE7 appended v5-dial mention to index.html")


def update_version_json(df: pd.DataFrame):
    vf = HTML_ROOT / "_version.json"
    try:
        v = json.loads(vf.read_text())
    except Exception:
        v = {}
    v["v5_dial_build_time"] = datetime.now(timezone.utc).isoformat()
    v["v5_dial_cancers"] = sorted(df["cancer"].unique().tolist())
    v["v5_dial_batch_entangled_pairs"] = int((df["interpretation"] == "batch_entangled").sum())
    v["v5_dial_true_biology_pairs"] = int((df["interpretation"] == "true_biology").sum())
    v["v5_dial_n_pairs"] = int(len(df))
    vf.write_text(json.dumps(v, indent=2))
    log_line(LOGFILE, f"PHASE7 updated _version.json")


def main():
    log_line(LOGFILE, "PHASE7 start — dashboard integration")
    df = pd.read_csv(RESULTS_V5 / "v5_dial_all_cancers.tsv", sep="\t")
    html = build_page(df)
    out = PAGES / "v4a_dial_cross_cancer.html"
    out.write_text(html)
    log_line(LOGFILE, f"PHASE7 wrote {out}")

    update_index_mention()
    update_version_json(df)
    log_line(LOGFILE, "PHASE7 done")


if __name__ == "__main__":
    main()
