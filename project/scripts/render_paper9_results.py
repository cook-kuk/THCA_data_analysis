#!/usr/bin/env python3
"""Render Paper 9 first-pass results md → HTML + Paper 9 plan & results → PDF.

PDF build via WeasyPrint. NO Unifont in font-family (per memory weasyprint_unifont_gotcha).
"""
from pathlib import Path

import markdown
from weasyprint import HTML

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
HUB = ROOT / "project/papers_hub_2026_05_04"

# 1) results md → HTML
RES_MD = ROOT / "project/reports/2026_05_04_paper9_sl_first_pass_results.md"
RES_HTML = HUB / "paper9_first_pass.html"
res_body = markdown.markdown(
    RES_MD.read_text(),
    extensions=["extra", "tables", "toc", "sane_lists"],
)

shell_results = """<!DOCTYPE html>
<html lang="ko"><head>
<meta charset="UTF-8" /><meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Paper 9 — First-Pass Results | Synthetic-lethal vulnerability map</title>
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;600&family=Noto+Sans+KR:wght@400;500;700&display=swap" rel="stylesheet" />
<link rel="stylesheet" href="style.css" />
<style>
.plan-doc h1{font-family:'Cormorant Garamond',serif;color:var(--accent);font-size:36px;line-height:1.15;margin:18px 0 10px}
.plan-doc h2{font-family:'Cormorant Garamond',serif;color:var(--ink);font-size:26px;margin:32px 0 10px;border-bottom:2px solid var(--rule);padding-bottom:6px}
.plan-doc h3{font-family:'Cormorant Garamond',serif;color:var(--ink-soft);font-size:19px;margin:22px 0 8px}
.plan-doc h4{font-family:'JetBrains Mono',monospace;font-size:13px;color:var(--accent);text-transform:uppercase;margin:18px 0 6px}
.plan-doc p{margin:0 0 12px;line-height:1.7;font-size:15.5px}
.plan-doc strong{color:var(--accent)}
.plan-doc em{color:var(--ink-soft);font-style:italic}
.plan-doc ul,.plan-doc ol{margin:0 0 14px 24px}
.plan-doc li{margin-bottom:5px;line-height:1.65}
.plan-doc table{border-collapse:collapse;margin:14px 0;width:100%;font-size:13.5px}
.plan-doc th,.plan-doc td{border:1px solid var(--rule);padding:7px 11px;text-align:left;vertical-align:top}
.plan-doc th{background:#F4EFE6;font-family:'JetBrains Mono',monospace;font-size:11px;text-transform:uppercase}
.plan-doc code{font-family:'JetBrains Mono',monospace;font-size:13px;background:#F4EFE6;padding:1px 6px;border-radius:3px}
.plan-doc pre{background:#F4EFE6;padding:10px 14px;border-radius:4px;overflow-x:auto;font-size:13px;font-family:'JetBrains Mono',monospace;margin:10px 0}
.plan-doc hr{border:0;border-top:1px solid var(--rule);margin:30px 0}
.figure-strip{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:24px 0}
.figure-strip figure{margin:0}
.figure-strip img{width:100%;border:1px solid var(--rule);border-radius:4px;background:#fff;padding:6px}
.figure-strip figcaption{font-size:12px;color:var(--ink-soft);font-style:italic;margin-top:5px}
.print-btn{position:fixed;top:16px;right:16px;z-index:9999;padding:10px 16px;background:var(--accent);color:#fff;border:none;border-radius:6px;font-family:'JetBrains Mono',monospace;font-size:13px;font-weight:600;cursor:pointer;box-shadow:0 2px 10px rgba(0,0,0,.18)}
.print-btn:hover{background:#000}
@media print{
  .print-btn,.topnav{display:none}
  body{background:#fff;color:#000}
  .container{max-width:100%;padding:0;margin:0}
  table,figure{page-break-inside:avoid}
  h1,h2,h3{page-break-after:avoid}
  a{color:#000;text-decoration:none}
}
</style>
</head><body>
<button class="print-btn" onclick="window.print()" title="이 페이지를 PDF로 인쇄/저장">📄 PDF로 저장 (Ctrl+P)</button>
<div class="container">
<nav class="topnav">
  <a href="index.html">← Portfolio hub</a>
  <span class="center">Paper 9 — First-Pass Results</span>
  <a href="portfolio_paper9.html">Paper 9 detail →</a>
</nav>
<div class="figure-strip">
  <figure><img src="../results/paper9_sl_first_pass/figures/F2_ccle_dm1_target_corr.png" alt="F2 CCLE DM1 corr"/><figcaption>F2 — CCLE DM1 × candidate target Pearson r (n=13)</figcaption></figure>
  <figure><img src="../results/paper9_sl_first_pass/figures/F3_dge_candidate_ranking.png" alt="F3 DGE ranking"/><figcaption>F3 — DM1-vs-DM2 DGE candidate ranking (TCGA-THCA bulk, q-corrected)</figcaption></figure>
  <figure><img src="../results/paper9_sl_first_pass/figures/F6_tcga_target_survival.png" alt="F6 Cox HR forest"/><figcaption>F6 — TCGA-THCA univariate Cox HR per candidate (whole + DM1-high)</figcaption></figure>
</div>
<div class="plan-doc">
__BODY__
</div>
</div>
</body></html>"""

RES_HTML.write_text(shell_results.replace("__BODY__", res_body))
print(f"Wrote {RES_HTML} ({RES_HTML.stat().st_size:,} bytes)")

# 2) PDF build for plan + results + portfolio detail
PLAN_HTML = HUB / "paper9_plan.html"
DETAIL_HTML = HUB / "portfolio_paper9.html"

for html_path in [PLAN_HTML, RES_HTML, DETAIL_HTML]:
    pdf_path = html_path.with_suffix(".pdf")
    print(f"Building PDF: {html_path.name} -> {pdf_path.name}")
    HTML(filename=str(html_path), base_url=str(HUB)).write_pdf(str(pdf_path))
    print(f"  {pdf_path.stat().st_size/1024:.1f} KB")

print("Done.")
