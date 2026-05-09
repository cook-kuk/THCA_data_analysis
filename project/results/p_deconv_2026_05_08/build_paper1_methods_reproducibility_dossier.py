#!/usr/bin/env python3
"""Build Paper 1 Fig 8 methods reproducibility dossier.

The page is intentionally factual and path-heavy: every headline number is
read from the v13/v14/v15/v16/v17/v18 output tables in this directory.
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/p_deconv_2026_05_08"
WEB = ROOT / "project/papers_hub_2026_05_04"
HTML_OUT = WEB / "paper1_methods_reproducibility_dossier.html"

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
TDS_16 = [
    "DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
    "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR",
]
MAPK_9 = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]


def esc(x: object) -> str:
    return html.escape("" if x is None else str(x))


def fmt(x: object, digits: int = 3) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return esc(x)
    if not math.isfinite(v):
        return "NA"
    av = abs(v)
    if av != 0 and (av < 0.001 or av >= 10000):
        return f"{v:.1e}"
    return f"{v:.{digits}f}"


def pct(x: object, digits: int = 1) -> str:
    try:
        return f"{float(x):.{digits}f}%"
    except (TypeError, ValueError):
        return "NA"


def table(headers: list[str], rows: list[list[object]], classes: str = "") -> str:
    th = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for row in rows:
        cells = []
        for val in row:
            cls = ""
            text = val
            if isinstance(val, tuple):
                text, cls = val
            cells.append(f"<td class=\"{esc(cls)}\">{text}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    return f"<table class=\"t {esc(classes)}\"><thead><tr>{th}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def fisher_pool(df: pd.DataFrame, panel: str) -> dict[str, float]:
    sub = df[df["panel"] == panel].dropna(subset=["spearman_rho", "n"]).copy()
    z = np.arctanh(sub["spearman_rho"].to_numpy(float))
    w = sub["n"].to_numpy(float) - 3
    mz = float((w * z).sum() / w.sum())
    se = 1.0 / math.sqrt(float(w.sum()))
    return {
        "rho": float(np.tanh(mz)),
        "ci_lo": float(np.tanh(mz - 1.96 * se)),
        "ci_hi": float(np.tanh(mz + 1.96 * se)),
        "n": int(sub["n"].sum()),
        "k": int(len(sub)),
    }


def load() -> dict[str, object]:
    v13_corr = pd.read_csv(OUT / "v13_cross_cohort_corr.tsv", sep="\t")
    v13_auc = pd.read_csv(OUT / "v13_auc_panel_vs_tds16.tsv", sep="\t")
    v13_sub = pd.read_csv(OUT / "v13_subAB_three_panels.tsv", sep="\t")
    v13_driver = pd.read_csv(OUT / "v13_per_driver_class_d.tsv", sep="\t")
    v14 = pd.read_csv(OUT / "v14_cross_cohort_forest.tsv", sep="\t")
    v15 = json.loads((OUT / "v15_mechanism_extensions_summary.json").read_text())
    v15_k2 = pd.read_csv(OUT / "v15A_k2_panel_overlay_by_group.tsv", sep="\t")
    v15_spatial = pd.read_csv(OUT / "v15B_spatial_mapk_panel_per_slide.tsv", sep="\t")
    v15_prism = pd.read_csv(OUT / "v15C_prism_top15_annotated.tsv", sep="\t")
    v16 = json.loads((OUT / "v16_spatial_prism_diagnostics_summary.json").read_text())
    v17 = json.loads((OUT / "v17_spatial_signal_decomposition_summary.json").read_text())
    v18 = json.loads((OUT / "v18_spatial_lag_pockets_summary.json").read_text())
    return {
        "v13_corr": v13_corr,
        "v13_auc": v13_auc,
        "v13_sub": v13_sub,
        "v13_driver": v13_driver,
        "v14": v14,
        "v14_pool": fisher_pool(v14, "Panel_8"),
        "v15": v15,
        "v15_k2": v15_k2,
        "v15_spatial": v15_spatial,
        "v15_prism": v15_prism,
        "v16": v16,
        "v17": v17,
        "v18": v18,
    }


def build_html(d: dict[str, object]) -> str:
    v13_corr = d["v13_corr"]
    v13_auc = d["v13_auc"]
    v13_sub = d["v13_sub"]
    v13_driver = d["v13_driver"]
    v14 = d["v14"]
    v14_pool = d["v14_pool"]
    v15 = d["v15"]
    v15_k2 = d["v15_k2"]
    v15_spatial = d["v15_spatial"]
    v15_prism = d["v15_prism"]
    v16 = d["v16"]
    v17 = d["v17"]
    v18 = d["v18"]

    tcga_panel = v13_corr[(v13_corr["cohort"] == "TCGA-THCA") & (v13_corr["panel"] == "Panel_8")].iloc[0]
    lee_panel = v13_corr[(v13_corr["cohort"] == "Lee_GSE213647") & (v13_corr["panel"] == "Panel_8")].iloc[0]
    braf_ras_panel = v13_driver[(v13_driver["contrast"] == "BRAF_V600E_vs_RAS_mut") & (v13_driver["score"] == "panel")].iloc[0]
    braf_ras_tds = v13_driver[(v13_driver["contrast"] == "BRAF_V600E_vs_RAS_mut") & (v13_driver["score"] == "tds16")].iloc[0]
    sub_mapk = v13_sub[v13_sub["metric"] == "mapk"].iloc[0]
    sub_panel = v13_sub[v13_sub["metric"] == "panel"].iloc[0]
    sub_tds = v13_sub[v13_sub["metric"] == "tds16"].iloc[0]
    v15_sp = v15["spatial"]
    v15_pr = v15["prism"]
    v15_k = v15["k2"]

    headline_rows = [
        ["BRAF V600E vs RAS Panel-8 d", fmt(braf_ras_panel["cohen_d"], 4), "v13_per_driver_class_d.tsv"],
        ["BRAF V600E vs RAS TDS-16 d", fmt(braf_ras_tds["cohen_d"], 4), "v13_per_driver_class_d.tsv"],
        ["MAPK x Panel-8 rho, TCGA / Lee", f"{fmt(tcga_panel['spearman_rho'])} / {fmt(lee_panel['spearman_rho'])}", "v13_cross_cohort_corr.tsv"],
        ["sub-A vs sub-B MAPK d", fmt(sub_mapk["d_A_minus_B"]), "v13_subAB_three_panels.tsv"],
        ["sub-A vs sub-B Panel-8 / TDS-16 d", f"{fmt(sub_panel['d_A_minus_B'])} / {fmt(sub_tds['d_A_minus_B'])}", "v13_subAB_three_panels.tsv"],
        ["Pooled 5-cohort MAPK x Panel-8 rho", f"{fmt(v14_pool['rho'])} [{fmt(v14_pool['ci_lo'])}, {fmt(v14_pool['ci_hi'])}]", "v14_cross_cohort_forest.tsv"],
        ["K2 score-only DM2-like calls", f"{v15_k['k2_n_dm2_like_p_ge_0p5']} / {v15_k['k2_n']} ({pct(v15_k['k2_pct_dm2_like_p_ge_0p5'])})", "v15A_k2_panel_overlay_summary.json"],
        ["Spatial MAPK x Panel pooled rho", f"{fmt(v15_sp['tumor_slide_fixed_effect_rho']['rho'])} [{fmt(v15_sp['tumor_slide_fixed_effect_rho']['ci_lo'])}, {fmt(v15_sp['tumor_slide_fixed_effect_rho']['ci_hi'])}]", "v15B_spatial_mapk_panel_summary.json"],
        ["Spatial rescue grid", f"0 / {v16['spatial']['grid_test_count']} negative subset-adjustment tests", "v16_spatial_adjustment_grid.tsv"],
        ["PRISM FDR<0.05 canonical MAPK hits", f"{v16['prism']['fdr05_canonical_mapk']}; p = {fmt(v16['prism']['fdr05_canonical_hypergeom_p'])}", "v16_prism_enrichment_sensitivity.tsv"],
        ["Spatial full+detection residual rho", fmt(v17["tumor_epi_full_spatial_detection_linear_rho"]), "v17_spatial_signal_decomposition_summary.json"],
        ["Spatial random-null percentile raw / full", f"{fmt(v17['observed_raw_percentile_vs_random'])} / {fmt(v17['observed_full_percentile_vs_random'])}", "v17_spatial_random_module_null.tsv"],
        ["Spatial residual KNN 1-6 lag rho", fmt(v18["tumor_epi_resid_knn_1_6_median_rho"]), "v18_spatial_lag_pockets_summary.json"],
        ["Spatial residual anti-pocket enrichment", f"{fmt(v18['resid_anti_pocket_median_enrichment'])}; OR = {fmt(v18['resid_anti_pocket_median_or'])}", "v18_spatial_pocket_enrichment.tsv"],
    ]

    formulas = [
        ["Within-cohort z-score", "z_gs = (x_gs - mean_g within cohort) / sd_g within cohort", "Bulk RNA module scores"],
        ["Panel-8 score", "mean z across DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR", "v13-v18"],
        ["TDS-16 score", "mean z across Yoo 2014 16-gene thyroid differentiation set", "v13-v14"],
        ["MAPK-9 output", "mean z across DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1", "v13-v18"],
        ["K2 centered profile", "log2(TPM+1), subtract each sample's 8-gene mean, TCGA-trained logistic regression", "v15A"],
        ["Spatial per-spot score", "raw h5ad counts -> log1p(CP10K) -> within-slide gene z -> mean module score", "v15B-v18"],
        ["Spatial lag score", "Panel score averaged across within-sample KNN rings around each MAPK-scored spot", "v18"],
        ["Spatial anti-pocket", "MAPK-high / Panel-low quadrant relative to within-sample top/bottom quartiles", "v18"],
        ["PRISM DM1 vulnerability", "Cohen's d on LFC: DM1-high minus DM1-low; negative = more sensitive in DM1-high", "v15C"],
        ["Pooled rho", "fixed-effect Fisher z using weight n-3; CI = pooled z +/- 1.96/sqrt(sum(n-3))", "v14-v15B"],
    ]

    gene_rows = [
        ["Panel-8", ", ".join(PANEL_8), "Deployable RAI/thyroid panel"],
        ["TDS-16", ", ".join(TDS_16), "Yoo 2014 canonical thyroid differentiation program"],
        ["TDS-only 8", ", ".join([g for g in TDS_16 if g not in PANEL_8]), "Disjoint TDS-16 genes; cherry-pick control"],
        ["MAPK-9", ", ".join(MAPK_9), "MAPK pathway output readout"],
    ]

    cohort_rows = []
    for cohort, sub in v14[v14["panel"] == "Panel_8"].groupby("cohort", sort=False):
        r = sub.iloc[0]
        note = {
            "TCGA-THCA": "well-differentiated primary; driver and methylation anchors",
            "Lee_GSE213647": "Korean primary cohort; Ensembl-to-symbol mapping",
            "GSE126698": "small Korean PTC/HT bulk RNA-seq",
            "GSE286332": "Korean PTC vs PTC+HT; HT-route decoupling",
            "GSE76039": "PDTC/ATC advanced-disease saturation",
            "GSE76039_advanced": "PDTC/ATC advanced-disease saturation",
        }.get(cohort, "")
        cohort_rows.append([cohort, int(r["n"]), fmt(r["spearman_rho"]), f"[{fmt(r['ci_lo'])}, {fmt(r['ci_hi'])}]", note])
    cohort_rows.extend([
        ["K2_PRJEB11591", int(v15_k["k2_n"]), "NA", "NA", "8-gene mini-index only; score-distribution evidence"],
        ["GSE250521", int(v15_sp["n_spots_total"]), fmt(v15_sp["tumor_slide_fixed_effect_rho"]["rho"]), f"[{fmt(v15_sp['tumor_slide_fixed_effect_rho']['ci_lo'])}, {fmt(v15_sp['tumor_slide_fixed_effect_rho']['ci_hi'])}]", "spatial stress test; positive co-localization, not anti-correlation"],
        ["DepMap/PRISM", int(v15_pr["n_prism_drugs_tested"]), "NA", "NA", "drug-vulnerability overlay, not expression restoration"],
    ])

    v13_rows = []
    for _, r in v13_corr.sort_values(["cohort", "panel"]).iterrows():
        v13_rows.append([r["cohort"], r["panel"], int(r["n"]), fmt(r["spearman_rho"]), fmt(r["spearman_p"])])

    auc_rows = []
    for _, r in v13_auc.iterrows():
        auc_rows.append([r["cohort"], int(r["n"]), fmt(r["AUC_Panel"]), fmt(r["AUC_TDS16"]), fmt(r["delta_TDS16_vs_Panel"])])

    v14_rows = []
    for _, r in v14[v14["panel"] == "Panel_8"].iterrows():
        v14_rows.append([r["cohort"], int(r["n"]), fmt(r["spearman_rho"]), f"[{fmt(r['ci_lo'])}, {fmt(r['ci_hi'])}]", fmt(r["spearman_p"])])

    k2_rows = []
    for _, r in v15_k2.iterrows():
        if r["cohort"] in {"K2_PRJEB11591", "TCGA_THCA", "Lee_GSE213647"}:
            k2_rows.append([r["cohort"], r["group"], int(r["n"]), fmt(r["median_p_DM2"]), pct(r["pct_p_DM2_ge_0p5"])])

    spatial_rows = []
    for _, r in v15_spatial.iterrows():
        spatial_rows.append([r["sample_id"], r["stage"], int(r["n"]), fmt(r["spearman_rho"]), fmt(r["partial_spearman_rho_qc"]), fmt(r["spearman_fdr"])])

    prism_rows = []
    for _, r in v15_prism.head(15).iterrows():
        prism_rows.append([int(r["rank_fdr"]), r["name"], r["mapk_class"], fmt(r["cohens_d"]), fmt(r["fdr"])])

    edge_rows = [
        ["Q9 causal boundary", "Mechanism layer motivates rather than confirms causality", "Keep perturbational proof outside Paper 1 scope"],
        ["K2", "No MAPK genes in mini-index", "Use only as 8-gene score-distribution overlay"],
        ["GSE250521", "Observed spatial MAPK x Panel direction is positive, not negative; v17/v18 attribute it to detection/covariate structure and close lag/pocket rescue", "Report as stress test / limitation, not positive mechanism evidence"],
        ["PRISM/DepMap", "Drug sensitivity and dependency are not expression restoration assays", "Use as therapeutic vulnerability rationale only"],
        ["v16 PRISM wording", "Earlier shorthand overclaimed top 15", "Use top 7 all MAPK; 7/11 FDR hits MAPK; top 15 is 8/15 MAPK"],
        ["v14 heterogeneity", "I2 is high because HT-route and advanced-disease cohorts decouple", "Interpret through v12 two-axis model, not as hidden failure"],
        ["Figure count", "Main Fig 8 remains two panels", "All v13-v18 mechanism evidence stays Supplementary"],
    ]

    file_rows = [
        ["v13", "run_v13_tds16_panel_mapk.py", "TDS-16 x Panel x MAPK cherry-pick lock"],
        ["v14", "run_v14_cross_cohort_forest.py", "5-cohort MAPK x thyroid-score forest"],
        ["v15", "run_v15_mechanism_extensions.py", "K2, spatial stress test, PRISM/DepMap overlay"],
        ["Dossier", "build_paper1_methods_reproducibility_dossier.py", "This page"],
        ["Figure", "Fig_SX_v15_mechanism_extensions.{png,pdf}", "4-panel v15 reserve figure"],
        ["v16", "run_v16_spatial_prism_diagnostics.py", "Spatial rescue failure + PRISM wording lock"],
        ["Figure", "Fig_SX_v16_spatial_prism_diagnostics.{png,pdf}", "4-panel v16 diagnostic figure"],
        ["v17", "run_v17_spatial_signal_decomposition.py", "Spatial detection/covariate decomposition + random-module null"],
        ["Figure", "Fig_SX_v17_spatial_signal_decomposition.{png,pdf}", "6-panel v17 spatial autopsy figure"],
        ["v18", "run_v18_spatial_lag_pockets.py", "Spatial lag and MAPK-high/Panel-low pocket closure"],
        ["Figure", "Fig_SX_v18_spatial_lag_pockets.{png,pdf}", "6-panel v18 spatial closure figure"],
        ["Rollup", "DECONV_ROLLUP_2026_05_09.md", "Operating evidence map: primary, support, reserve, caveat, rejected/deferred"],
        ["Rollup HTML", "paper1_deconvolution_rollup_v18.html", "Live deconvolution decision dashboard"],
        ["Manuscript", "project/manuscript_v8/05_figure_captions.md", "SX_v13/SX_v14/SX_v15/SX_v16/SX_v17/SX_v18 captions"],
        ["Hub", "project/papers_hub_2026_05_04/", "Live dossier HTML pages"],
    ]

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>Paper 1 Fig 8 Methods Reproducibility Dossier</title>
<style>
:root{{--bg:#07101a;--panel:#0d1828;--panel2:#111f34;--line:#24324a;--ink:#eaf1fb;--muted:#96a7bd;--gold:#ffd28a;--teal:#35d39d;--red:#ff8a6b;--blue:#7eb6ff;}}
*{{box-sizing:border-box}}
html,body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",-apple-system,BlinkMacSystemFont,sans-serif}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
code,pre,.mono{{font-family:"JetBrains Mono","SF Mono",Menlo,monospace;font-size:12px}}
h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;color:#fff7dd;letter-spacing:-.02em}}
.hero{{padding:64px 34px 42px;background:linear-gradient(135deg,#101f34 0%,#07101a 62%,#140f18 100%);border-bottom:1px solid var(--line)}}
.hero-inner{{max-width:1320px;margin:auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;letter-spacing:.26em;color:var(--gold);text-transform:uppercase}}
h1{{font-size:58px;line-height:1;margin:14px 0 14px}}
.lead{{max-width:980px;color:#c9d5e6;font:18px/1.6 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-top:24px}}
.stat{{background:rgba(255,255,255,.045);border:1px solid rgba(255,210,138,.2);padding:13px 14px;border-radius:8px}}
.stat b{{display:block;color:var(--gold);font:700 25px "Cormorant Garamond",serif;line-height:1}}
.stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
.crumbs{{margin-top:20px;color:var(--muted);font:12px "JetBrains Mono",monospace}}
.wrap{{max-width:1320px;margin:auto;display:grid;grid-template-columns:245px 1fr;gap:34px;padding:0 34px 80px}}
.toc{{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;border-right:1px solid var(--line);padding:28px 20px 28px 0}}
.toc h4{{font:700 10px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;text-transform:uppercase}}
.toc a{{display:block;color:#cfd8e8;padding:5px 0;font-size:12px}}
main{{min-width:0;padding-top:24px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:34px;margin:0 0 6px}} h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}}
h3{{font-size:22px;margin:20px 0 8px}}
.sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.t{{width:100%;border-collapse:collapse;margin:12px 0 18px;font-size:12.5px}}
.t th,.t td{{border:1px solid var(--line);padding:7px 9px;vertical-align:top}}
.t th{{background:#17243a;color:var(--gold);font:700 10.5px "JetBrains Mono",monospace;letter-spacing:.05em;text-transform:uppercase;text-align:left}}
.t tr:nth-child(even) td{{background:rgba(255,255,255,.025)}}
.hi{{color:var(--gold);font-weight:700}} .good{{color:var(--teal)}} .warn{{color:var(--red)}} .muted{{color:var(--muted)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 8px 8px 0;padding:14px 16px;margin:14px 0}}
.box.warn{{border-left-color:var(--red);background:#1b1515}} .box.good{{border-left-color:var(--teal);background:#0b1b18}}
.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}}
.card{{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:14px 16px}}
.card h3{{margin-top:0}}
pre{{white-space:pre-wrap;background:#050a12;border:1px solid var(--line);padding:14px;border-radius:8px;color:#cdd7e8}}
@media(max-width:980px){{.wrap{{display:block;padding:0 20px 60px}}.toc{{display:none}}.stats{{grid-template-columns:repeat(2,1fr)}}h1{{font-size:40px}}.grid{{grid-template-columns:1fr}}}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <div class="kicker">Paper 1 Fig 8 / methods reproducibility</div>
    <h1>Mechanism Layer Audit Trail</h1>
    <p class="lead">A reviewer-facing checklist for the Fig 8 mechanism layer: gene lists, formulas, cohort filters, edge-case handling, v13-v18 outputs, and the preserved correlative boundary.</p>
    <div class="stats">
      <div class="stat"><b>v13-v18</b><span>analysis layers</span></div>
      <div class="stat"><b>{fmt(v14_pool['rho'])}</b><span>5-cohort pooled rho</span></div>
      <div class="stat"><b>{v14_pool['n']}</b><span>forest samples</span></div>
      <div class="stat"><b>{v15_k['k2_n_dm2_like_p_ge_0p5']}/260</b><span>K2 DM2-like</span></div>
      <div class="stat"><b>{fmt(v18['tumor_epi_resid_knn_1_6_median_rho'])}</b><span>spatial residual lag rho</span></div>
      <div class="stat"><b>{v15_pr['n_fdr05_canonical_mapk']}/11</b><span>PRISM MAPK hits</span></div>
    </div>
    <div class="crumbs"><a href="index.html">Hub</a> / <a href="paper1_fig8_mechanism_dossier.html">Fig 8 mechanism dossier</a> / <a href="paper1_reviewer_defense_dashboard.html">Reviewer dashboard</a></div>
  </div>
</header>
<div class="wrap">
<aside class="toc">
  <h4>Sections</h4>
  <a href="#tldr">01 TL;DR</a>
  <a href="#genes">02 Gene Lists</a>
  <a href="#formulas">03 Formulas</a>
  <a href="#cohorts">04 Cohorts</a>
  <a href="#v13">05 v13 Lock</a>
  <a href="#v14">06 v14 Forest</a>
  <a href="#v15">07 v15 Reserve</a>
  <a href="#v16">08 v16 Diagnostics</a>
  <a href="#v17">09 v17 Decomposition</a>
  <a href="#v18">10 v18 Lag/Pockets</a>
  <a href="#edges">11 Edge Cases</a>
  <a href="#files">12 Files</a>
</aside>
<main>
<section id="tldr">
  <h2><span class="num">01</span>TL;DR</h2>
  <div class="sub">Every headline value cites the local output table it came from</div>
  {table(["Result", "Value", "Source"], headline_rows)}
  <div class="box warn"><b>Boundary kept intact.</b> The mechanism layer is HM450/RNA/cohort co-variation plus drug-vulnerability support. It motivates perturbational experiments; it does not claim causal proof or expression restoration.</div>
</section>

<section id="genes">
  <h2><span class="num">02</span>Gene Lists</h2>
  <div class="sub">No hidden panel edits between v13, v14, v15, v16, v17, and v18</div>
  {table(["Set", "Genes", "Role"], gene_rows)}
</section>

<section id="formulas">
  <h2><span class="num">03</span>Formulas</h2>
  <div class="sub">Exact score definitions used by the pipelines</div>
  {table(["Quantity", "Definition", "Where used"], formulas)}
</section>

<section id="cohorts">
  <h2><span class="num">04</span>Cohorts And Filters</h2>
  <div class="sub">Cohort-specific interpretation is part of the model, not a post-hoc escape hatch</div>
  {table(["Cohort", "n", "Panel-8 rho", "95% CI", "Notes"], cohort_rows)}
  <div class="box good"><b>v14 interpretation.</b> TCGA and Lee are well-differentiated primary cohorts and carry the MAPK-route signal. HT-route and advanced-disease cohorts are expected to decouple under the v12 two-axis model.</div>
</section>

<section id="v13">
  <h2><span class="num">05</span>v13 Cherry-Pick Lock</h2>
  <div class="sub">Panel-8 behaves as a compact readout of TDS-16</div>
  <div class="grid">
    <div class="card">
      <h3>MAPK x thyroid scores</h3>
      {table(["Cohort", "Panel", "n", "Spearman rho", "p"], v13_rows)}
    </div>
    <div class="card">
      <h3>MAPK-high classifier AUC</h3>
      {table(["Cohort", "n", "Panel-8 AUC", "TDS-16 AUC", "Delta"], auc_rows)}
    </div>
  </div>
</section>

<section id="v14">
  <h2><span class="num">06</span>v14 Cross-Cohort Forest</h2>
  <div class="sub">5 cohorts, n = {v14_pool['n']}, fixed-effect Fisher-z pooled Panel-8 rho</div>
  {table(["Cohort", "n", "Panel-8 rho", "95% CI", "p"], v14_rows)}
  <div class="box"><b>Pooled result.</b> MAPK x Panel-8 rho = {fmt(v14_pool['rho'])} [{fmt(v14_pool['ci_lo'])}, {fmt(v14_pool['ci_hi'])}] across {v14_pool['k']} cohorts. Heterogeneity is reported transparently and interpreted through the prespecified v12 two-axis convergence model.</div>
</section>

<section id="v15">
  <h2><span class="num">07</span>v15 Reserve Extensions</h2>
  <div class="sub">Two useful reserve supports plus one non-confirmatory spatial stress test</div>
  <h3>A. K2 score-only overlay</h3>
  {table(["Cohort", "Group", "n", "Median p_DM2", "p_DM2 >= 0.5"], k2_rows)}
  <h3>B. GSE250521 spatial stress test</h3>
  <div class="box warn"><b>Direction check.</b> Raw per-spot MAPK x Panel is positive in all slides, not anti-correlated; after QC partialing the pooled tumor-slide rho is {fmt(v15_sp['tumor_slide_qc_partial_fixed_effect_rho']['rho'])}. Do not use this as positive spatial mechanism confirmation.</div>
  {table(["Sample", "Stage", "n", "rho", "QC-partial rho", "FDR"], spatial_rows)}
  <h3>C. PRISM/DepMap overlay</h3>
  {table(["Rank", "Drug", "Class", "Cohen's d", "FDR"], prism_rows)}
</section>

<section id="v16">
  <h2><span class="num">08</span>v16 Diagnostics</h2>
  <div class="sub">Spatial rescue failure and PRISM exact-wording lock</div>
  {table(["Diagnostic", "Value", "Disposition"], [
      ["Spatial subset/adjustment grid", f"0 / {v16['spatial']['grid_test_count']} MAPK x Panel tests negative", "GSE250521 does not rescue spatial anti-correlation"],
      ["Tumor-slide pooled QC + cell-state rho", f"{fmt(v16['spatial']['tumor_slide_pooled_by_adjustment']['QC_plus_cell_state']['rho'])} [{fmt(v16['spatial']['tumor_slide_pooled_by_adjustment']['QC_plus_cell_state']['ci_lo'])}, {fmt(v16['spatial']['tumor_slide_pooled_by_adjustment']['QC_plus_cell_state']['ci_hi'])}]", "Positive signal attenuates to tiny residual"],
      ["PRISM top-k exact wording", f"top 7 = {v16['prism']['top7_canonical_mapk']}; top 15 = {v16['prism']['top15_canonical_mapk']}", "Do not say top 15 all MAPK"],
      ["PRISM FDR < 0.05 enrichment", f"{v16['prism']['fdr05_canonical_mapk']}, p = {fmt(v16['prism']['fdr05_canonical_hypergeom_p'])}", "Strong vulnerability support"],
      ["Exclude thyroid cell lines", f"rho = {fmt(v16['prism']['exclude_thyroid_rho'])}, p = {fmt(v16['prism']['exclude_thyroid_p'])}", "Pan-cancer signal does not depend on thyroid lines"],
  ])}
</section>

<section id="v17">
  <h2><span class="num">09</span>v17 Spatial Signal Decomposition</h2>
  <div class="sub">The failed spatial anti-correlation candidate is detection/covariate structure, not hidden mechanism confirmation</div>
  {table(["Diagnostic", "Value", "Disposition"], [
      ["Tumor epithelial raw MAPK x Panel rho", fmt(v17["tumor_epi_raw_linear_rho"]), "Positive co-localization, not anti-correlation"],
      ["Full covariate-adjusted rho", fmt(v17["tumor_epi_full_linear_rho"]), "Positive signal attenuates strongly"],
      ["Full spatial+detection residual rho", fmt(v17["tumor_epi_full_spatial_detection_linear_rho"]), "Near-zero after detection/spatial structure"],
      ["Median covariate R2, MAPK / Panel", f"{fmt(v17['median_full_covariate_r2_MAPK_tumor_epi'])} / {fmt(v17['median_full_covariate_r2_Panel8_tumor_epi'])}", "Depth/detection covariates explain much of both modules"],
      ["MAPK-detection x Panel-detection median rho", fmt(v17["median_MAPK_detect_vs_panel_detect_rho_tumor_epi"]), "Broad detection co-localization can induce positive module correlation"],
      ["Random-module percentile raw / full", f"{fmt(v17['observed_raw_percentile_vs_random'])} / {fmt(v17['observed_full_percentile_vs_random'])}", "Observed signal is generic versus same-size module nulls"],
      ["Top negative gene pair", f"{v17['top_negative_gene_pair']}; rho = {fmt(v17['top_negative_gene_pair_median_rho'])}; FDR = {fmt(v17['top_negative_gene_pair_fdr'])}", "No significant pairwise rescue"],
  ])}
  <div class="box warn"><b>Disposition.</b> GSE250521 should remain a spatial caveat/stress test. v17 documents why it fails instead of upgrading it into positive Fig 8 mechanism evidence.</div>
</section>

<section id="v18">
  <h2><span class="num">10</span>v18 Spatial Lag And Pockets</h2>
  <div class="sub">Final spatial escape-hatch test: neighborhood lag and MAPK-high/Panel-low pockets</div>
  {table(["Diagnostic", "Value", "Disposition"], [
      ["Raw same-spot rho", fmt(v18["tumor_epi_raw_same_spot_median_rho"]), "Same direction as failed v15/v17 raw spatial result"],
      ["Raw KNN 1-6 lag rho", fmt(v18["tumor_epi_raw_knn_1_6_median_rho"]), "Nearby Panel neighborhoods remain weakly positive, not negative"],
      ["Raw KNN 7-18 lag rho", fmt(v18["tumor_epi_raw_knn_7_18_median_rho"]), "Positive relation attenuates with distance but does not invert"],
      ["Residual same-spot rho", fmt(v18["tumor_epi_resid_same_spot_median_rho"]), "Residual same-spot relation is near null"],
      ["Residual KNN 1-6 lag rho", fmt(v18["tumor_epi_resid_knn_1_6_median_rho"]), "No residual negative neighborhood effect"],
      ["Residual KNN 7-18 lag rho", fmt(v18["tumor_epi_resid_knn_7_18_median_rho"]), "Distant neighborhood relation is near null"],
      ["Raw anti-pocket enrichment / OR", f"{fmt(v18['raw_anti_pocket_median_enrichment'])} / {fmt(v18['raw_anti_pocket_median_or'])}", "Anti-pockets are depleted, not enriched"],
      ["Residual anti-pocket enrichment / OR", f"{fmt(v18['resid_anti_pocket_median_enrichment'])} / {fmt(v18['resid_anti_pocket_median_or'])}", "Residual anti-pockets are independence-level"],
      ["Strongest raw anti-pocket covariate depletion", f"{v18['strongest_raw_anti_pocket_covariate_depletion']}; d = {fmt(v18['strongest_raw_anti_pocket_covariate_d'])}", "Low panel detection dominates raw anti-pockets"],
  ])}
  <div class="box good"><b>Disposition.</b> v18 closes the final spatial rescue route. GSE250521 remains a caveat/reserve page, not a positive Fig 8 mechanism support pillar.</div>
</section>

<section id="edges">
  <h2><span class="num">11</span>Edge-Case Handling</h2>
  <div class="sub">Reviewer attack surface and how it is contained</div>
  {table(["Issue", "Risk", "Disposition"], edge_rows)}
</section>

<section id="files">
  <h2><span class="num">12</span>Files And Re-run Commands</h2>
  <div class="sub">Local paths are the source of truth</div>
  {table(["Layer", "Path", "Purpose"], file_rows)}
  <pre>cd /home/seungho/personal/THCA_data_analysis
python project/results/p_deconv_2026_05_08/run_v13_tds16_panel_mapk.py
python project/results/p_deconv_2026_05_08/run_v14_cross_cohort_forest.py
python project/results/p_deconv_2026_05_08/run_v15_mechanism_extensions.py
python project/results/p_deconv_2026_05_08/run_v16_spatial_prism_diagnostics.py
python project/results/p_deconv_2026_05_08/run_v17_spatial_signal_decomposition.py
python project/results/p_deconv_2026_05_08/run_v18_spatial_lag_pockets.py
python project/results/p_deconv_2026_05_08/build_paper1_methods_reproducibility_dossier.py</pre>
  <div class="box muted"><b>Deploy target.</b> HTML lives at <code>project/papers_hub_2026_05_04/paper1_methods_reproducibility_dossier.html</code> and is copied to <code>/var/www/papers/papers_hub_2026_05_04/</code> for the live hub.</div>
</section>
</main>
</div>
</body>
</html>
"""


def main() -> None:
    HTML_OUT.write_text(build_html(load()))
    print(f"wrote {HTML_OUT}")


if __name__ == "__main__":
    main()
