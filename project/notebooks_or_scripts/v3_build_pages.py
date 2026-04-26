"""v3 STEP 8 — Build HTML pages 19-24 + index.html banner + markdown summaries.

Pages 19 (honesty audit), 20 (panel size curve), 21 (three class),
      22 (survival), 23 (bethesda sim), 24 (multimodal).
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v3_common import FIGS, HTML, PAGES, RESULTS_TABLES

INDEX = HTML / "index.html"
VERSION_JSON = HTML / "_version.json"

# ---------- page template ----------
PAGE_TMPL = """<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="dark light">
  <meta name="theme-color" content="#07132b">
  <title>{page_title} | THCA Multi-Omics Dashboard (v3)</title>
  <meta name="description" content="{page_subtitle}">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="stylesheet" href="../assets/css/main.css">
  <style>
    body{{background:linear-gradient(180deg,#07132b 0%,#0a1a3a 100%);color:#e2e8f0;
         font-family:Inter,system-ui,sans-serif;margin:0;padding:0}}
    .v3-badge{{display:inline-block;padding:4px 10px;border-radius:999px;
               background:linear-gradient(135deg,#14b8a6,#06b6d4);color:#041b1f;
               font-weight:700;font-size:.7rem;letter-spacing:.14em;text-transform:uppercase}}
    .v3-hero{{padding:36px 32px;border-radius:20px;margin:22px;
             background:linear-gradient(135deg,rgba(6,17,38,.9),rgba(11,26,58,.7));
             border:1px solid rgba(94,234,212,.25)}}
    .v3-hero h1{{font-size:clamp(1.6rem,1.1rem+1.5vw,2.4rem);margin:8px 0 6px;color:#fff}}
    .v3-hero p{{color:#94a3b8;max-width:780px;margin:0}}
    .v3-kpi{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:18px 22px}}
    .v3-kpi div{{padding:14px 16px;border-radius:12px;border:1px solid rgba(94,234,212,.25);
                background:linear-gradient(135deg,rgba(20,184,166,.12),rgba(6,182,212,.05))}}
    .v3-kpi strong{{display:block;font-size:1.4rem;color:#e2fff9;line-height:1}}
    .v3-kpi span{{display:block;font-size:.72rem;color:#94a3b8;margin-top:4px;
                 text-transform:uppercase;letter-spacing:.12em}}
    .v3-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(440px,1fr));
              gap:16px;margin:22px}}
    .v3-card{{padding:16px;border-radius:16px;border:1px solid rgba(148,163,184,.15);
             background:rgba(11,26,58,.45)}}
    .v3-card h3{{margin:0 0 6px;font-size:1.05rem;color:#e2e8f0}}
    .v3-card p.note{{margin:0 0 10px;font-size:.82rem;color:#94a3b8}}
    .v3-card iframe{{width:100%;height:520px;border:0;border-radius:10px;
                     background:rgba(7,19,43,.55)}}
    .v3-footer{{margin:24px 22px;padding:16px;border-radius:12px;
                background:rgba(11,26,58,.45);border:1px solid rgba(148,163,184,.12);
                color:#94a3b8;font-size:.88rem;display:flex;justify-content:space-between;
                gap:12px;flex-wrap:wrap;align-items:center}}
    .v3-footer a.btn{{padding:6px 12px;border-radius:8px;background:rgba(94,234,212,.14);
                      color:#5eead4;text-decoration:none;border:1px solid rgba(94,234,212,.3)}}
    .topnav{{position:sticky;top:0;z-index:5;background:rgba(7,19,43,.85);backdrop-filter:blur(10px);
            border-bottom:1px solid rgba(148,163,184,.12)}}
    .topnav-inner{{display:flex;align-items:center;justify-content:space-between;
                   padding:10px 22px;gap:12px}}
    .brand{{color:#5eead4;font-weight:700;letter-spacing:.12em;text-transform:uppercase;
           text-decoration:none}}
    .nav-links{{display:flex;gap:8px;flex-wrap:wrap;font-size:.72rem}}
    .nav-links a{{color:#94a3b8;text-decoration:none;padding:4px 8px;border-radius:6px}}
    .nav-links a:hover{{color:#e2e8f0;background:rgba(148,163,184,.08)}}
    .nav-links a.active{{color:#5eead4;background:rgba(94,234,212,.1)}}
  </style>
</head>
<body class="dark">
  <nav class="topnav" aria-label="Main">
    <div class="topnav-inner">
      <a class="brand" href="../index.html">THYROID DASH / v3</a>
      <div class="nav-links">
        <a href="../index.html">HOME</a>
        <a href="14_caveats.html">CAVEATS</a>
        <a href="19_honesty_audit.html"{active_19}>19 HONESTY</a>
        <a href="20_panel_size.html"{active_20}>20 PANEL-k</a>
        <a href="21_three_class.html"{active_21}>21 3-CLASS</a>
        <a href="22_survival.html"{active_22}>22 SURVIVAL</a>
        <a href="23_bethesda_sim.html"{active_23}>23 BETHESDA</a>
        <a href="24_multimodal.html"{active_24}>24 MULTIMODAL</a>
      </div>
    </div>
  </nav>
  <main>
    <section class="v3-hero">
      <span class="v3-badge">v3 · {build_time}</span>
      <h1>{page_title}</h1>
      <p>{page_subtitle}</p>
    </section>
    {kpi_section}
    <section class="v3-grid">
      {cards_html}
    </section>
    <div class="v3-footer">
      <div>
        Decision-support prototype · exploratory · retrospective computational triage.
        See <a href="14_caveats.html" style="color:#5eead4">page 14 (caveats)</a>
        and <a href="19_honesty_audit.html" style="color:#5eead4">page 19 (honesty audit)</a>.
      </div>
      <a class="btn" href="../index.html">Return home</a>
    </div>
  </main>
</body>
</html>
"""


def kpi_html(kpis):
    if not kpis:
        return ""
    items = "".join(f"<div><strong>{v}</strong><span>{k}</span></div>"
                    for k, v in kpis)
    return f'<section class="v3-kpi">{items}</section>'


def card_html(title, note, fig_src):
    return (f'<div class="v3-card"><h3>{title}</h3>'
            f'<p class="note">{note}</p>'
            f'<iframe loading="lazy" src="{fig_src}" '
            f'title="{title}"></iframe></div>')


def render_page(page_filename, page_title, page_subtitle, cards, kpis=None,
                 active_key: str = "19"):
    active_map = {k: "" for k in ["19", "20", "21", "22", "23", "24"]}
    active_map[active_key] = ' class="active"'
    html = PAGE_TMPL.format(
        page_title=page_title,
        page_subtitle=page_subtitle,
        build_time=datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC"),
        kpi_section=kpi_html(kpis or []),
        cards_html="\n      ".join(card_html(*c) for c in cards),
        active_19=active_map["19"], active_20=active_map["20"],
        active_21=active_map["21"], active_22=active_map["22"],
        active_23=active_map["23"], active_24=active_map["24"],
    )
    out = PAGES / page_filename
    out.write_text(html, encoding="utf-8")
    return out


def safe_read_tsv(path):
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def load_leakage_k10():
    df = safe_read_tsv(RESULTS_TABLES / "v3_leakage_curve.tsv")
    if df.empty:
        return None
    r = df[df["k_removed"] == 10]
    if r.empty:
        return None
    return float(r.iloc[0]["mean_auc"])


def build_page_19():
    # KPIs from honesty audit
    kpis = []
    lev = safe_read_tsv(RESULTS_TABLES / "v3_leakage_curve.tsv")
    if not lev.empty:
        r0 = lev[lev.k_removed == 0]
        r10 = lev[lev.k_removed == 10]
        if not r0.empty:
            kpis.append(("TCGA AUC (k=0 removed)", f"{float(r0.iloc[0]['mean_auc']):.3f}"))
        if not r10.empty:
            kpis.append(("TCGA AUC (k=10 removed)", f"{float(r10.iloc[0]['mean_auc']):.3f}"))
    mapk = safe_read_tsv(RESULTS_TABLES / "v3_mapk_ablation.tsv")
    if not mapk.empty:
        for _, row in mapk.iterrows():
            kpis.append((row["variant"], f"{float(row['mean_cv_auc']):.3f}"))
    perm = safe_read_tsv(RESULTS_TABLES / "v3_permutation_null.tsv")
    if not perm.empty:
        kpis.append(("Permutation p", f"{float(perm.iloc[0]['empirical_p']):.4f}"))

    cards = [
        ("Leakage curve",
         "Bootstrap TCGA AUC as top-k |Cohen's d| genes are removed from TierA67_clean.",
         "../figs_interactive/v3_leakage_curve.html"),
        ("MAPK-output ablation",
         "Full TierA67_clean vs TierA67_clean minus MAPK_OUTPUT (DUSP4/5/6, SPRY1/2/4, ETV4/5, FOSL1, PHLDA1).",
         "../figs_interactive/v3_mapk_ablation.html"),
        ("Permutation null",
         "1000 label shuffles; empirical p for observed AUC.",
         "../figs_interactive/v3_permutation_null.html"),
        ("Dataset identifiability",
         "OvR LogReg on TierA67_clean features predicting which cohort a tumor came from (batch leakage).",
         "../figs_interactive/v3_dataset_identifiability.html"),
        ("LODO AUC",
         "Leave-one-dataset-out: train TCGA-only vs TCGA+others, test external.",
         "../figs_interactive/v3_lodo.html"),
        ("Calibration (TCGA → GSE27155)",
         "Reliability diagram with isotonic and Platt recalibration on GSE27155 internal 5-fold.",
         "../figs_interactive/v3_calibration_gse27155.html"),
        ("Decision curve (GSE27155)",
         "Net benefit vs threshold for isotonic-recalibrated model vs treat-all/treat-none.",
         "../figs_interactive/v3_decision_curve_gse27155.html"),
    ]
    return render_page("19_honesty_audit.html",
                        "Honesty Audit",
                        "TCGA AUC drops once genes that double as BRAF-vs-RAS label-definers are removed. Six orthogonal leakage tests.",
                        cards, kpis=kpis, active_key="19")


def build_page_20():
    best = {}
    try:
        best = json.loads((RESULTS_TABLES / "v3_panel_best.json").read_text())
    except Exception:
        pass

    df = safe_read_tsv(RESULTS_TABLES / "v3_panel_size_curve.tsv")
    kpis = []
    if best:
        kpis.append(("Best strategy", best.get("strategy", "n/a")))
        kpis.append(("Best k", str(best.get("k", "n/a"))))
        kpis.append(("GSE27155 AUC", f"{best.get('gse27155_auc', float('nan')):.3f}"))
        kpis.append(("TCGA AUC", f"{best.get('tcga_auc', float('nan')):.3f}"))

    cards = [
        ("Panel-k curve (3 strategies, 3 cohorts)",
         "In-fold feature selection: univariate |d|, TDS16 union, and QUBO-neal. Dashed = external AUCs (isotonic recalibrated).",
         "../figs_interactive/v3_panel_size_curve.html"),
        ("NPV at Se≥0.95 heatmap (GSE27155)",
         "Rule-out performance: how pure the negatives are when we catch ≥95% of true positives.",
         "../figs_interactive/v3_panel_size_npv_heatmap.html"),
    ]
    return render_page("20_panel_size.html",
                        "Panel Size Curve",
                        "Exploratory sweep of panel sizes k ∈ {4…60} and three feature-selection strategies, evaluated on TCGA CV + 2 external cohorts.",
                        cards, kpis=kpis, active_key="20")


def build_page_21():
    kpis = []
    df = safe_read_tsv(RESULTS_TABLES / "v3_three_class_metrics.tsv")
    if not df.empty:
        macro = df[(df.model == "logreg") & (df.class_name == "MACRO")]
        if not macro.empty:
            kpis.append(("LogReg macro AUC", f"{float(macro.iloc[0]['ovr_auc']):.3f}"))
        gbm = df[(df.model == "gb") & (df.class_name == "MACRO")]
        if not gbm.empty:
            kpis.append(("GradBoost macro AUC", f"{float(gbm.iloc[0]['ovr_auc']):.3f}"))
        for cls in ["BRAF_like", "RAS_like", "dedifferentiated"]:
            r = df[(df.model == "logreg") & (df.class_name == cls)]
            if not r.empty:
                kpis.append((f"{cls} n", str(int(r.iloc[0]['n']))))

    cards = [
        ("3-class OvR ROC (LogReg)",
         "BRAF_like vs RAS_like vs dedifferentiated. Small-n subclasses flagged in legend.",
         "../figs_interactive/v3_three_class_roc.html"),
        ("Confusion matrix (LogReg, 5-fold CV)",
         "Where the model confuses RAS_like with dedifferentiated.",
         "../figs_interactive/v3_three_class_confusion.html"),
        ("SHAP top-15 per class (GradBoost)",
         "Mean|SHAP| on 200-sample stratified holdout. Features = TierA67_clean minus MAPK_OUTPUT.",
         "../figs_interactive/v3_three_class_shap.html"),
    ]
    return render_page("21_three_class.html",
                        "Three-Class Classifier",
                        "BRAF_like / RAS_like / dedifferentiated — multinomial LogReg + GradientBoosting, with SHAP interpretability.",
                        cards, kpis=kpis, active_key="21")


def build_page_22():
    kpis = []
    df = safe_read_tsv(RESULTS_TABLES / "v3_survival_cox.tsv")
    if not df.empty:
        mv = df[df["mode"] == "multivariate"]
        kpis.append(("Multivariate rows", str(len(mv))))
        try:
            n_total = int(df["n"].max())
            kpis.append(("Max n (Cox)", str(n_total)))
        except Exception:
            pass

    cards = [
        ("OS by molecular subtype",
         "BRAF_like vs RAS_like vs dedifferentiated.",
         "../figs_interactive/v3_km_subtype.html"),
        ("OS by TDS16 tertile",
         "Thyroid differentiation score tertile (low/mid/high).",
         "../figs_interactive/v3_km_tds_tertile.html"),
        ("OS by dediff-proxy tertile",
         "Dedifferentiation proxy score tertile.",
         "../figs_interactive/v3_km_dediff_tertile.html"),
        ("Multivariate Cox forest",
         "Hazard ratios adjusting for age + stage + sex.",
         "../figs_interactive/v3_cox_forest.html"),
    ]
    return render_page("22_survival.html",
                        "Survival (TCGA-THCA)",
                        "KM + log-rank pairwise + univariate/multivariate Cox. TCGA-THCA has few OS events — interpret with caution.",
                        cards, kpis=kpis, active_key="22")


def build_page_23():
    kpis = []
    surg = safe_read_tsv(RESULTS_TABLES / "v3_bethesda_surgery_reduction.tsv")
    if not surg.empty:
        for _, r in surg.iterrows():
            kpis.append((f"prev={r['prev']:.0%} reduction",
                         f"{r['reduction_vs_treat_all'] * 100:.0f}%"))

    cards = [
        ("Operating curves across prevalence",
         "Se / Sp / NPV / unnecessary-surgery rate vs decision threshold, one trace per prevalence.",
         "../figs_interactive/v3_bethesda_operating.html"),
        ("Decision curve (DCA)",
         "Net benefit vs threshold across prevalences, with treat-all/treat-none reference.",
         "../figs_interactive/v3_bethesda_decision_curve.html"),
        ("Surgery reduction at Se≥0.95",
         "Reduction in surgery rate vs treat-all at sensitivity ≥0.95, by prevalence.",
         "../figs_interactive/v3_bethesda_surgery_reduction.html"),
    ]
    return render_page("23_bethesda_sim.html",
                        "Bethesda Prevalence Simulation",
                        "Synthetic cohort (N=10,000 per prevalence, σ=0.3 Gaussian noise). Decision-support prototype only — not clinically validated.",
                        cards, kpis=kpis, active_key="23")


def build_page_24():
    kpis = []
    try:
        fusion = json.loads((RESULTS_TABLES / "v3_fusion_results.json").read_text())
        kpis.append(("Fusion status", fusion.get("status", "n/a")))
    except Exception:
        pass
    try:
        scna = json.loads((RESULTS_TABLES / "v3_scna_results.json").read_text())
        kpis.append(("SCNA status", scna.get("status", "n/a")))
    except Exception:
        pass
    meth = safe_read_tsv(RESULTS_TABLES / "v3_methylation_results.tsv")
    if not meth.empty:
        r = meth.iloc[0]
        kpis.append((f"GSE97466 n={int(r['n_samples'])}",
                     f"AUC={r['mean_cv_auc']:.3f}"
                     if not pd.isna(r['mean_cv_auc']) else "AUC=NA"))

    cards = [
        ("Methylation (GSE97466) — within-cohort",
         "Tumor-vs-normal CV AUC on beta top-5000. No TCGA methylation available "
         "locally, so inter-cohort LODO is infeasible.",
         "../figs_interactive/v3_multimodal_methylation.html"),
    ]
    return render_page("24_multimodal.html",
                        "Multimodal (methylation only)",
                        "Fusion callset and SCNA not locally available — skipped stubs. Methylation analysis limited to within-cohort (GSE97466).",
                        cards, kpis=kpis, active_key="24")


# ---------- index.html banner ----------
BANNER_START = '<!-- v3-banner:start -->'
BANNER_END = '<!-- v3-banner:end -->'


def append_banner_to_index():
    auc_k10 = load_leakage_k10()
    if auc_k10 is None:
        auc_text = "(AUC pending)"
    else:
        auc_text = f"{auc_k10:.3f}"

    banner = (f'{BANNER_START}\n'
              f'<div class="v3-banner" role="status" '
              f'style="padding:16px 28px;margin:0;background:linear-gradient(90deg,#0f766e,#14b8a6);'
              f'color:#03241f;font-weight:600;font-size:.95rem;border-bottom:1px solid rgba(3,36,31,.3)">'
              f'🔬 <strong>v3 update posted.</strong> Honesty audit shows TierA67_clean TCGA AUC drops '
              f'from 1.00 to {auc_text} after removing MAPK-output circular features. '
              f'<a href="pages/19_honesty_audit.html" '
              f'style="color:#03241f;text-decoration:underline">See page 19</a>.'
              f'</div>\n'
              f'{BANNER_END}\n')

    html = INDEX.read_text(encoding="utf-8")
    # Remove previous banner if present
    html = re.sub(
        re.escape(BANNER_START) + r".*?" + re.escape(BANNER_END) + r"\n?",
        "",
        html, flags=re.DOTALL,
    )
    # Insert just after <body...>
    m = re.search(r'<body[^>]*>', html)
    if m is None:
        # fallback: prepend
        html = banner + html
    else:
        idx = m.end()
        html = html[:idx] + "\n" + banner + html[idx:]
    INDEX.write_text(html, encoding="utf-8")


# ---------- version bump ----------
def update_version():
    try:
        v = json.loads(VERSION_JSON.read_text())
    except Exception:
        v = {}
    v["v3_build_time"] = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    v["v3_pages"] = ["19_honesty_audit.html", "20_panel_size.html",
                      "21_three_class.html", "22_survival.html",
                      "23_bethesda_sim.html", "24_multimodal.html"]
    VERSION_JSON.write_text(json.dumps(v, indent=2))


def run():
    outs = [
        build_page_19(),
        build_page_20(),
        build_page_21(),
        build_page_22(),
        build_page_23(),
        build_page_24(),
    ]
    append_banner_to_index()
    update_version()
    return {"pages": [str(o) for o in outs]}


def main():
    print(json.dumps(run(), indent=2, default=str))


if __name__ == "__main__":
    main()
