#!/usr/bin/env python3
from __future__ import annotations

import html
from pathlib import Path

import pandas as pd

from v17p2_common import FIG, RPT, TAB


def read_json(path: Path) -> dict:
    import json
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def fig_card(name: str, title: str) -> str:
    return f"<div class='card'><h3>{html.escape(title)}</h3><iframe src='../../results/v17p2/figs/{name}' loading='lazy'></iframe><a href='../../results/v17p2/figs/{name}' target='_blank' rel='noopener'>새 창</a></div>"


def html_table(df: pd.DataFrame, tid: str) -> str:
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    body = []
    for _, row in df.iterrows():
        body.append("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in row.tolist()) + "</tr>")
    return f"<table id='{tid}' class='display compact stripe'><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def wrap(title: str, body: str) -> str:
    return f"""<!doctype html>
<html lang='ko'>
<head>
  <meta charset='UTF-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>{html.escape(title)}</title>
  <link rel='preconnect' href='https://fonts.googleapis.com'>
  <link rel='preconnect' href='https://fonts.gstatic.com' crossorigin>
  <link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700&display=swap' rel='stylesheet'>
  <link rel='stylesheet' href='https://cdn.jsdelivr.net/npm/datatables.net-dt@1.13.8/css/jquery.dataTables.min.css'>
  <style>
    :root{{--bg:#07132b;--panel:#0f172a;--fg:#e2e8f0;--accent:#f59e0b;}}
    body{{margin:0;background:linear-gradient(180deg,#020617,#07132b);color:var(--fg);font-family:Inter,'Noto Sans KR',sans-serif}}
    .wrap{{max-width:1380px;margin:0 auto;padding:28px}}
    .hero,.panel{{background:rgba(15,23,42,.85);border:1px solid rgba(255,255,255,.08);border-radius:24px;padding:24px}}
    .hero h1{{font-size:48px;line-height:1.05;margin:8px 0 12px}}
    .muted{{color:#94a3b8}}
    .grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px}}
    .metric{{background:rgba(255,255,255,.04);border-radius:18px;padding:16px;border:1px solid rgba(255,255,255,.08)}}
    .metric b{{display:block;font-size:30px}}
    .section{{margin-top:18px}}
    .cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:14px}}
    .card{{background:#081120;border:1px solid rgba(255,255,255,.08);border-radius:18px;padding:14px}}
    .card iframe{{width:100%;height:340px;border:0;background:#000;border-radius:12px}}
    a{{color:#7dd3fc}}
    table.dataTable{{width:100%!important}}
  </style>
</head>
<body>
<div class='wrap'>{body}</div>
<script src='https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js'></script>
<script src='https://cdn.jsdelivr.net/npm/datatables.net@1.13.8/js/jquery.dataTables.min.js'></script>
<script>document.querySelectorAll('table.display').forEach(el=>new DataTable(el,{{pageLength:10,order:[]}}));</script>
</body></html>"""


def main() -> None:
    robust = read_json(TAB / "robustness_summary.json")
    clinical = read_json(TAB / "clinical_summary.json")
    biology = read_json(TAB / "biology_summary.json")
    external = read_json(TAB / "external_summary.json")
    dial = read_json(TAB / "dial_method_summary.json")
    hero = f"""
    <section class='hero'>
      <div class='muted'>v17 Phase 2 / 5-Layer Validation</div>
      <h1>DM1/DM2 5-Layer Validation</h1>
      <p>Dark Matter subtype를 단순 clustering 결과가 아니라, 통계적 강건성·임상 분화·생물학적 기전·외부 검증·방법론 audit 다섯 층으로 다시 검증한 한국어 리포트다.</p>
      <div class='grid'>
        <div class='metric'><b>{robust.get('median_persistence','NA')}</b><span>Bootstrap persistence</span></div>
        <div class='metric'><b>{clinical.get('age_p_welch','NA')}</b><span>Age p-value</span></div>
        <div class='metric'><b>{external.get('cohorts_dia_auc_gt_0_85','NA')}</b><span>Cohorts DIA-AUC > 0.85</span></div>
        <div class='metric'><b>{dial.get('combat_threshold_lambda_ident_lt_0_5','NA')}</b><span>ComBat threshold</span></div>
      </div>
    </section>
    """
    summary_df = pd.DataFrame([
        {"Layer": "1 Robustness", "Key metric": "mean off-diagonal ARI", "Value": robust.get("mean_ari_offdiag")},
        {"Layer": "2 Clinical", "Key metric": "OS events", "Value": clinical.get("os_events")},
        {"Layer": "3 Biology", "Key metric": "Hallmark FDR<0.05", "Value": biology.get("hallmark_fdr_lt_0_05")},
        {"Layer": "4 External", "Key metric": "robust cohorts", "Value": external.get("cohorts_dia_auc_gt_0_85")},
        {"Layer": "5 DIAL", "Key metric": "pancancer tested", "Value": dial.get("pancancer_tested")},
    ])
    body = hero + f"""
    <section class='panel section'>
      <h2>핵심 메트릭</h2>
      {html_table(summary_df, 'summary')}
    </section>
    <section class='panel section'><h2>Layer 1 — Robustness</h2><div class='cards'>
      {fig_card('robustness_method_consistency.html','5-method consistency')}
      {fig_card('bootstrap_consensus.html','Bootstrap consensus')}
      {fig_card('permutation_null_distribution.html','Permutation null')}
      {fig_card('feature_ablation_curve.html','Feature ablation')}
    </div></section>
    <section class='panel section'><h2>Layer 2 — Clinical</h2><div class='cards'>
      {fig_card('km_curves_panel.html','Kaplan-Meier')}
      {fig_card('cox_forest_plot.html','Cox forest')}
      {fig_card('age_violin_with_pvalues.html','Age separation')}
      {fig_card('stage_mosaic.html','Stage distribution')}
    </div></section>
    <section class='panel section'><h2>Layer 3 — Biology</h2><div class='cards'>
      {fig_card('gsea_top_pathways_bar.html','GSEA top pathways')}
      {fig_card('hallmark_heatmap_per_sample.html','Hallmark heatmap')}
      {fig_card('immune_composition_stacked.html','Immune proxy')}
      {fig_card('mapk_activity_violin.html','MAPK activity')}
      {fig_card('thyroid_diff_score_violin.html','Thyroid differentiation')}
    </div></section>
    <section class='panel section'><h2>Layer 4 — External Validation</h2><div class='cards'>
      {fig_card('external_4cohort_forest.html','External forest')}
      {fig_card('external_confusion_matrices.html','Classifier comparison')}
      {fig_card('scrna_dm_signature_umap.html','scRNA projection')}
      {fig_card('scrna_per_patient_heterogeneity.html','Per-patient heterogeneity')}
    </div></section>
    <section class='panel section'><h2>Layer 5 — DIAL Method Audit</h2><div class='cards'>
      {fig_card('dial_combat_threshold.html','ComBat threshold')}
      {fig_card('dial_pancancer_forest.html','Pan-cancer feasibility')}
      {fig_card('dial_vs_alternatives_radar.html','Alternatives comparison')}
      {fig_card('dial_failure_mode_simulation.html','Failure-mode simulation')}
    </div></section>
    """
    (RPT / "index.html").write_text(wrap("v17 Phase 2", body), encoding="utf-8")


if __name__ == "__main__":
    main()

