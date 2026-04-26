#!/usr/bin/env python
"""v5.1 dashboard page builder.

Rewrites /opt/thyroid-dash/project/reports/html/pages/v4a_dial_cross_cancer.html
with the THCA-specificity narrative. Preserves the existing site chrome
(nav, font preloads, CSS includes) and embeds the 5 new Plotly figures.
Reads TSVs from results/v5/ for live counts.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v5"
HTML_PAGES = PROJECT / "reports" / "html" / "pages"
OUT = HTML_PAGES / "v4a_dial_cross_cancer.html"


def load_counts() -> dict:
    dial = pd.read_csv(RES / "v5p1_dial_all_cancers.tsv", sep="\t")
    total = len(dial)
    counts = dial.interpretation.value_counts().to_dict()
    by_cancer = (
        dial.groupby(["cancer", "interpretation"]).size().unstack(fill_value=0)
    )
    true_bio = {c: int(by_cancer.loc[c, "true_biology"]) if "true_biology" in by_cancer.columns and c in by_cancer.index else 0
                for c in by_cancer.index}
    return {
        "total": total,
        "batch_entangled": int(counts.get("batch_entangled", 0)),
        "true_biology": int(counts.get("true_biology", 0)),
        "no_signal": int(counts.get("no_signal", 0)),
        "ambiguous": int(counts.get("ambiguous", 0)),
        "partial_batch": int(counts.get("partial_batch", 0)),
        "by_cancer": by_cancer,
        "true_bio_per_cancer": true_bio,
    }


PAGE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="dark">
  <title>v5.1 THCA-Specific DIAL Audit — THYRAI Research Archive</title>
  <meta name="description" content="v5.1 LODO audit on real multi-cohort TCGA + GEO data. THCA BRAF/RAS-specific batch entanglement confirmed; other cancers are null.">
  <meta property="og:title" content="v5.1 THCA-Specific DIAL Audit">
  <meta property="og:description" content="4 of 5 THCA classifiers flip direction post-ComBat (DIAL up to 0.494). All 20 non-THCA classifier-cancer combinations null.">
  <meta property="og:image" content="../figs_interactive/v5/v5p1_fig1_interpretation_heatmap.png">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="../assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/jetbrains-mono-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="../assets/css/brand.css">
  <link rel="stylesheet" href="../assets/css/v4a.css">
  <style>
    .dial-hero {{
      min-height: 48vh;
      padding: 112px 0 56px;
      border-bottom: 1px solid var(--border);
      background:
        radial-gradient(ellipse at 70% 20%, rgba(245,166,35,0.10) 0%, transparent 55%),
        radial-gradient(ellipse at 20% 80%, rgba(245,166,35,0.06) 0%, transparent 60%);
    }}
    .dial-hero .eyebrow {{ color: #F5A623; letter-spacing: 0.18em; }}
    .dial-hero h1 {{
      font-size: clamp(32px, 4.2vw, 54px);
      font-weight: 300; letter-spacing: -0.02em; line-height: 1.08;
      margin: 18px 0 14px;
    }}
    .dial-hero .sub {{ color: #d2a34e; font-size: 16px; letter-spacing: 0.02em; }}
    .dial-hero .lede {{ color: var(--text-secondary); font-size: 17px; line-height: 1.6; max-width: 72ch; margin-top: 16px; }}
    .tile-row {{
      display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 20px; margin: 30px 0 10px;
    }}
    .tile {{
      background: rgba(20,22,28,0.75);
      border: 1px solid var(--border);
      border-radius: 4px;
      padding: 20px 22px;
    }}
    .tile__value {{
      font-family: 'JetBrains Mono', monospace;
      font-size: 28px; font-weight: 500; color: #F5A623; line-height: 1;
    }}
    .tile__label {{
      font-family: 'JetBrains Mono', monospace;
      font-size: 11px; letter-spacing: 0.14em; text-transform: uppercase;
      color: var(--text-muted); margin-top: 10px;
    }}
    .dial-card {{
      background: rgba(20,22,28,0.75);
      border: 1px solid var(--border);
      border-radius: 4px;
      padding: 30px 32px; margin: 24px 0;
    }}
    .dial-card h2 {{ margin: 0 0 14px; font-size: 24px; font-weight: 400; color: #EAEAEA; letter-spacing: -0.01em; }}
    .dial-card .eyebrow {{ color: #F5A623; font-size: 11px; letter-spacing: 0.18em; font-family: 'JetBrains Mono', monospace; }}
    .dial-card p {{ color: var(--text-secondary); font-size: 15px; line-height: 1.66; max-width: 72ch; }}
    .dial-card ol, .dial-card ul {{ color: var(--text-secondary); font-size: 15px; line-height: 1.66; }}
    .dial-card ol li, .dial-card ul li {{ margin: 6px 0; }}
    .dial-fig {{
      background: #0b0e12;
      border: 1px solid var(--border);
      border-radius: 4px;
      margin: 22px 0;
      overflow: hidden;
    }}
    .dial-fig iframe {{ width: 100%; border: 0; display: block; }}
    .dial-fig__caption {{
      padding: 14px 22px; border-top: 1px solid var(--border);
      font-family: 'JetBrains Mono', monospace; font-size: 12px;
      color: var(--text-muted); letter-spacing: 0.04em;
    }}
    .caveat-band {{
      background: rgba(194,76,76,0.06);
      border-left: 3px solid #c24c4c;
      padding: 20px 24px; margin: 40px 0 24px;
      font-size: 14px; line-height: 1.6; color: #c9d1d9;
    }}
    .caveat-band h4 {{
      margin: 0 0 8px; color: #c24c4c; font-size: 12px; letter-spacing: 0.12em;
      text-transform: uppercase; font-weight: 500;
    }}
    .archive-list {{ font-family: 'JetBrains Mono', monospace; font-size: 13px; }}
    .archive-list li {{ margin: 4px 0; color: #c9d1d9; }}
    .archive-list code {{ color: #F5A623; background: rgba(245,166,35,0.08); padding: 1px 5px; border-radius: 2px; }}
    .retracted-note {{ color: #c24c4c; font-style: italic; font-size: 13px; margin-top: 10px; }}
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
          <li><a href="v4a_dial_cross_cancer.html" aria-current="page"><b>v5.1 THCA-Specific Audit</b></a></li>
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
    <span class="eyebrow">v5.1 · REAL DATA · LODO · REVISED</span>
    <h1>v5.1 REAL DATA + LODO — THCA BRAF/RAS-specific phenomenon confirmed.</h1>
    <p class="sub">v5 (StratifiedKFold + semi-synthetic) retracted.</p>
    <p class="lede">
      We re-ran the DIAL audit on five TCGA cancers using <b>leave-one-dataset-out (LODO)</b> on the variance-top-3000 gene pool over real multi-cohort harmonizations. The label-flip regime is <b>specific to THCA BRAF-vs-RAS</b>: four of five classifiers flip direction post-ComBat with DIAL up to 0.494, while all twenty non-THCA classifier-cancer combinations are null (no flip). This is evidence of <i>specificity</i>, not of generalization failure of the method.
    </p>
    <div class="tile-row">
      <div class="tile"><div class="tile__value">{total}</div><div class="tile__label">Classifiers tested</div></div>
      <div class="tile"><div class="tile__value">{batch_entangled}</div><div class="tile__label">batch_entangled (THCA only)</div></div>
      <div class="tile"><div class="tile__value">{true_biology}</div><div class="tile__label">true_biology<br>(LGG {tb_lgg}, LUAD {tb_luad}, COAD {tb_coad})</div></div>
      <div class="tile"><div class="tile__value">{no_signal}</div><div class="tile__label">no_signal</div></div>
      <div class="tile"><div class="tile__value">{ambiguous}</div><div class="tile__label">ambiguous</div></div>
    </div>
  </div>
</section>

<section class="section">
  <div class="container">
    <div class="dial-card">
      <span class="eyebrow">SECTION A · THE SPECIFICITY FINDING</span>
      <h2>DIAL × interpretation across 5 cancers × 5 classifiers.</h2>
      <p>
        Each cell in the heatmap below reports the DIAL value and its categorical interpretation for one classifier-cancer combination. The red column is the finding: THCA with four classifiers in <b>batch_entangled</b> (LogReg_l2, LogReg_elasticnet, RandomForest, GradientBoosting). The other four cancers yield no batch_entangled pair — LGG and COAD are uniformly <b>true_biology</b>; LUAD is a mix of true_biology and no_signal; SKCM is ambiguous/no_signal. This shows the DIAL probe itself is not producing false positives in every cancer — the flip regime is a property of THCA's data generating process under multi-cohort harmonization.
      </p>
      <div class="dial-fig">
        <iframe src="../figs_interactive/v5/v5p1_fig1_interpretation_heatmap.html" width="100%" height="520"></iframe>
        <div class="dial-fig__caption">Figure 1 — Interpretation × DIAL per (cancer × classifier). Red = batch_entangled; green = true_biology; gray = no_signal; amber = ambiguous.</div>
      </div>
    </div>

    <div class="dial-card">
      <span class="eyebrow">SECTION B · LABEL-FLIP GEOMETRY IN THCA</span>
      <h2>Four THCA classifiers cross the y = 1 − x diagonal.</h2>
      <p>
        Figure 2 plots AUC-pre against AUC-post. The gray dashed line is the identity (no change under ComBat); the red dashed line is complete sign inversion. THCA points sit on the complete-flip line: LogReg_l2 (pre 0.994 → post 0.006), LogReg_elasticnet (0.994 → 0.008), RandomForest (0.996 → 0.177), GradientBoosting (0.750 → 0.166). XGBoost is the sole non-flip THCA case (0.781 → 0.487), consistent with tree-ensemble non-linearity breaking the linear projection that ComBat annihilates.
      </p>
      <div class="dial-fig">
        <iframe src="../figs_interactive/v5/v5p1_fig2_auc_pre_vs_post.html" width="100%" height="560"></iframe>
        <div class="dial-fig__caption">Figure 2 — AUC pre vs post. Marker size ∝ |DIAL|. THCA markers cluster on y = 1 − x; other cancers cluster near y = x.</div>
      </div>
      <p>
        Figure 3 is a classifier-level close-up: for each of the four flipped THCA classifiers, AUC-post is well below 0.5, but auc_flip_post (= max(a, 1 − a)) is well above 0.5. The predictive magnitude survives; only the direction flips. This is the direction-invariant leakage signature DIAL is designed to detect.
      </p>
      <div class="dial-fig">
        <iframe src="../figs_interactive/v5/v5p1_fig3_thca_label_flip.html" width="100%" height="620"></iframe>
        <div class="dial-fig__caption">Figure 3 — THCA label-flip geometry per classifier (2×2). Dashed line = chance (0.5).</div>
      </div>
    </div>

    <div class="dial-card">
      <span class="eyebrow">SECTION C · WHY ONLY THCA?</span>
      <h2>Three co-occurring conditions.</h2>
      <p>
        We hypothesize the THCA flip is the conjunction of three properties of its cohorts and its biology. All three are required; relaxing any one collapses the effect.
      </p>
      <ol>
        <li><b>Cohort × class-ratio alignment.</b> TCGA-THCA is 293 BRAF : 58 RAS; GSE27155 is 28 : 13. Both cohorts share the same majority-class skew and the same minority class, so subtype-preserving ComBat cannot distinguish between the subtype axis and the cohort mean.</li>
        <li><b>BRAF-vs-RAS biology projects onto a low-rank linear axis.</b> The MAPK-on vs MAPK-off distinction in PTC is well-captured by a single variance component in the 3,000-gene variance-top pool, which coincides with the cohort-mean direction ComBat subtracts.</li>
        <li><b>Harmonization collapses that axis.</b> ComBat re-centers by cohort mean; because biology and cohort are aligned, biology is mean-shifted with them, producing the near-total sign inversion visible in Figure 2.</li>
      </ol>
      <p>
        The other cancers violate at least one condition. <b>LGG</b> uses four TSS sub-splits instead of platforms, so cohort means are not aligned with the IDH axis and ComBat leaves biology intact (5/5 true_biology). <b>SKCM</b> has no strong linear subtype signature even pre-ComBat (AUC-pre 0.49–0.70), so there is no direction to flip. <b>LUAD</b> splits neatly into signal (LogReg_l2, RandomForest, LogReg_elasticnet at true_biology) and no_signal cases with GradientBoosting; no classifier reaches the flip regime. <b>COAD</b> has extreme class imbalance (76 : 296) but its subtype axis (CMS-like) is not cohort-aligned, so ComBat preserves the signal (5/5 true_biology).
      </p>
      <div class="dial-fig">
        <iframe src="../figs_interactive/v5/v5p1_fig5_cohort_label_imbalance.html" width="100%" height="560"></iframe>
        <div class="dial-fig__caption">Figure 5 — Per-cohort class-A vs class-B counts. THCA cohorts carry the shared skew that enables the flip; amber borders highlight THCA.</div>
      </div>
      <div class="dial-fig">
        <iframe src="../figs_interactive/v5/v5p1_fig4_specificity_evidence.html" width="100%" height="500"></iframe>
        <div class="dial-fig__caption">Figure 4 — Interpretation counts per cancer. The batch_entangled column is populated only by THCA.</div>
      </div>
    </div>

    <div class="dial-card">
      <span class="eyebrow">SECTION D · METHOD ROBUSTNESS</span>
      <h2>StratifiedKFold → LODO.</h2>
      <p>
        The earlier v5 run used StratifiedKFold inside the pooled matrix, which leaks cohort identity across folds: the pre-ComBat AUC inflates because the classifier partly identifies cohort rather than subtype. Switching to leave-one-dataset-out (LODO) with ComBat fit on the training cohorts and applied to the held-out cohort at inference is the correct pipeline for a cross-cohort claim. We also dropped semi-synthetic augmentations and used only real cohort expression matrices over the shared-gene intersection (ranges 9,179 to 20,964 genes before variance filtering, then the <b>variance-top-3000</b> pool per cancer for classifier training). Under LODO we do not re-gate on post-ComBat batch identifiability because LODO evaluation already holds out cohort; the interpretation rule used here is <b>DIAL &gt; 0.3</b> as the primary threshold for batch_entangled, with auc_post &gt; 0.7 required for true_biology, and auc_post &lt; 0.6 required for no_signal. THCA's four flipped classifiers sit well above 0.3; XGBoost on THCA is ambiguous by this rule and formally falls to no_signal.
      </p>
    </div>

    <div class="dial-card">
      <span class="eyebrow">SECTION E · ARCHIVE · RETRACTED v5 RUN</span>
      <h2>Superseded files.</h2>
      <p>The following files are the prior StratifiedKFold + semi-synthetic run that this page supersedes. They are retained for auditability only.</p>
      <ul class="archive-list">
        <li><code>results/v5_broken/v5_dial_all_cancers.tsv</code> — retracted</li>
        <li><code>results/v5_broken/v5_dial_THCA.tsv</code> — retracted</li>
        <li><code>results/v5_broken/v5_dial_SKCM.tsv</code> — retracted</li>
        <li><code>results/v5_broken/v5_dial_LGG.tsv</code> — retracted</li>
        <li><code>results/v5_broken/v5_dial_LUAD.tsv</code> — retracted</li>
        <li><code>results/v5_broken/v5_dial_COAD.tsv</code> — retracted</li>
        <li><code>results/v5_broken/v5_dial_paper.tex</code> — retracted draft</li>
      </ul>
      <p class="retracted-note">These files used StratifiedKFold over a pooled matrix and included semi-synthetic augmentations. They should not be cited. The v5.1 files under <code>results/v5/v5p1_*.tsv</code> are the current truth.</p>
    </div>

    <div class="caveat-band">
      <h4>Caveat</h4>
      DIAL is a <b>diagnostic</b> probe — it flags when a post-correction classifier has flipped direction because biology and cohort mean are aligned. It is <b>not a correction method</b>: a DIAL-flagged pair requires either a non-linear correction, a cohort-free feature pool, or a held-out independent cohort to be trusted. All v5.1 conclusions are <b>computational</b>; prospective cohort validation is pending through the ongoing collaboration with Prof. 유형원 at Seoul National University Bundang Hospital (유형원 교수 분당서울대 협력 진행 중).
    </div>
  </div>
</section>

<footer class="site-footer">
  <div class="container">
    <div class="site-footer__grid">
      <div class="site-footer__col site-footer__brand">
        <img src="../assets/img/brand/thyrai_wordmark.svg" alt="THYRAI">
        <p>Research dashboard for computational thyroid-cancer triage. Computational triage only. Not medical advice.</p>
      </div>
      <div class="site-footer__col">
        <h4>Research Archive</h4>
        <ul>
          <li><a href="v4a_archive_index.html">Archive index</a></li>
          <li><a href="v4a_dial_cross_cancer.html">v5.1 THCA-Specific Audit</a></li>
          <li><a href="v4a_ml_performance.html">ML performance</a></li>
          <li><a href="v4a_limitations.html">Limitations</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>v5.1 sources</h4>
        <ul>
          <li><a href="../../../results/v5/v5p1_dial_all_cancers.tsv">v5p1_dial_all_cancers.tsv</a></li>
          <li><a href="../../../results/v5/v5p1_table1_summary.tsv">v5p1_table1_summary.tsv</a></li>
          <li><a href="../../../results/v5/v5p1_harmonization.tsv">v5p1_harmonization.tsv</a></li>
          <li><a href="../../../reports/v5/v5p1_change_log.md">v5.1 change log</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>Legacy</h4>
        <ul>
          <li><a href="../index.html">Back to dashboard</a></li>
          <li><a href="99_glossary.html">99 · Glossary</a></li>
        </ul>
      </div>
    </div>
    <div class="site-footer__meta">
      <span>&copy; 2026 THYRAI Research · v5.1 revised</span>
      <span class="site-footer__disclaimer">Computational triage only. Not medical advice.</span>
    </div>
  </div>
</footer>
<script src="../assets/js/brand.js" defer></script>
</body>
</html>
"""


def main() -> int:
    c = load_counts()
    tb = c["true_bio_per_cancer"]
    html = PAGE.format(
        total=c["total"],
        batch_entangled=c["batch_entangled"],
        true_biology=c["true_biology"],
        no_signal=c["no_signal"],
        ambiguous=c["ambiguous"],
        tb_lgg=tb.get("LGG", 0),
        tb_luad=tb.get("LUAD", 0),
        tb_coad=tb.get("COAD", 0),
    )
    OUT.write_text(html, encoding="utf-8")
    size_kb = OUT.stat().st_size / 1024
    print(f"[v5p1_build_dashboard] wrote {OUT} ({size_kb:.1f} KB)")
    print(f"[v5p1_build_dashboard] total={c['total']} batch={c['batch_entangled']} tb={c['true_biology']} ns={c['no_signal']} amb={c['ambiguous']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
