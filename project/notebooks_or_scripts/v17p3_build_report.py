#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from v17p3_common import FIG, RPT, TAB


def j(name: str) -> dict:
    p = TAB / name
    return json.loads(p.read_text()) if p.exists() else {}


def card(fig: str, title: str) -> str:
    return f"<div class='card'><h3>{title}</h3><iframe src='../../results/v17p3/figs/{fig}' loading='lazy'></iframe></div>"


def main() -> None:
    f1, f2, f3 = j("F1_summary.json"), j("F2_summary.json"), j("F3_summary.json")
    a1, a2, a3 = j("A1_summary.json"), j("A2_summary.json"), j("A3_summary.json")
    a4, a5, a6 = j("A4_summary.json"), j("A5_summary.json"), j("A6_summary.json")
    body = f"""<!doctype html><html lang='ko'><head><meta charset='UTF-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>v17 Phase 3</title>
<link rel='preconnect' href='https://fonts.googleapis.com'><link rel='preconnect' href='https://fonts.gstatic.com' crossorigin>
<link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&family=Noto+Sans+KR:wght@400;500;700&display=swap' rel='stylesheet'>
<style>body{{font-family:Inter,'Noto Sans KR',sans-serif;background:#04101f;color:#e2e8f0;margin:0}} .wrap{{max-width:1400px;margin:0 auto;padding:24px}} .hero,.panel{{background:#0b1629;border:1px solid rgba(255,255,255,.08);border-radius:24px;padding:24px;margin-bottom:18px}} .grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px}} .metric{{background:rgba(255,255,255,.04);padding:14px;border-radius:16px}} .metric b{{display:block;font-size:28px}} .cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:14px}} .card{{background:#071120;border-radius:18px;padding:12px;border:1px solid rgba(255,255,255,.08)}} iframe{{width:100%;height:340px;border:0;border-radius:12px;background:#000}}</style></head><body><div class='wrap'>
<section class='hero'><div>v17 Phase 3</div><h1>DM1/DM2 — From Discovery to Clinical Actionability</h1>
<div class='grid'>
<div class='metric'><b>{f1.get('robust_cohort_count_dia_auc_gt_0_85','NA')}</b><span>robust cohorts</span></div>
<div class='metric'><b>{f2.get('hallmark_fdr_lt_0_01','NA')}</b><span>Hallmark FDR&lt;0.01</span></div>
<div class='metric'><b>{f3.get('n_endpoints_p_lt_0_05','NA')}</b><span>clinical endpoints</span></div>
<div class='metric'><b>{a3.get('n_drugs_fdr_lt_0_1','NA')}</b><span>selective drugs</span></div>
<div class='metric'><b>{a6.get('pdtc_auc_using_low_rai','NA')}</b><span>RAI validation AUC</span></div>
</div></section>
<section class='panel'><h2>Fix</h2><div class='cards'>
{card('F1_5cohort_forest.html','F1 External recovery')}
{card('F2_gsea_top_pathways_dual_dir.html','F2 Proper GSEA')}
{card('F3_endpoint_forest.html','F3 Clinical reframe')}
</div></section>
<section class='panel'><h2>Amplify</h2><div class='cards'>
{card('A1_dominance_per_patient_bar.html','A1 scRNA heterogeneity')}
{card('A2_dm_score_by_driver.html','A2 DM in BRAF/RAS-positive')}
{card('A3_drug_volcano_dm1_vs_dm2.html','A3 Drug response')}
{card('A4_pancancer_heatmap.html','A4 Pan-cancer transfer')}
{card('A5_immune_landscape.html','A5 Genomic/immune')}
{card('A6_rai_score_distribution.html','A6 TF/RAI')}
</div></section></div></body></html>"""
    (RPT / "index.html").write_text(body, encoding="utf-8")


if __name__ == "__main__":
    main()

