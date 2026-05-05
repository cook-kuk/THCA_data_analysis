#!/usr/bin/env python3
"""Render Paper 9 strategic plan md → styled HTML inside papers_hub_2026_05_04/."""
from pathlib import Path

import markdown

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SRC_MD = ROOT / "project/reports/2026_05_04_paper9_synthetic_lethality_ruppin_style_plan.md"
OUT_HTML = ROOT / "project/papers_hub_2026_05_04/paper9_plan.html"

md_text = SRC_MD.read_text()
body_html = markdown.markdown(
    md_text,
    extensions=["extra", "tables", "toc", "sane_lists"],
)

shell = """<!DOCTYPE html>
<html lang="ko"><head>
<meta charset="UTF-8" /><meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Paper 9 Strategic Plan | Ruppin-style synthetic-lethal vulnerability map</title>
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;600&family=Noto+Sans+KR:wght@400;500;700&display=swap" rel="stylesheet" />
<link rel="stylesheet" href="style.css" />
<style>
.plan-doc h1{font-family:'Cormorant Garamond',serif;color:var(--accent);font-size:36px;line-height:1.15;margin:18px 0 10px;letter-spacing:-0.5px}
.plan-doc h2{font-family:'Cormorant Garamond',serif;color:var(--ink);font-size:26px;margin:32px 0 10px;border-bottom:2px solid var(--rule);padding-bottom:6px}
.plan-doc h3{font-family:'Cormorant Garamond',serif;color:var(--ink-soft);font-size:19px;margin:22px 0 8px}
.plan-doc h4{font-family:'JetBrains Mono',monospace;font-size:13px;color:var(--accent);text-transform:uppercase;letter-spacing:.05em;margin:18px 0 6px}
.plan-doc p{margin:0 0 12px 0;line-height:1.7;font-size:15.5px}
.plan-doc strong{color:var(--accent)}
.plan-doc em{color:var(--ink-soft);font-style:italic}
.plan-doc ul,.plan-doc ol{margin:0 0 14px 24px}
.plan-doc li{margin-bottom:5px;line-height:1.65}
.plan-doc table{border-collapse:collapse;margin:14px 0;width:100%;font-size:14px}
.plan-doc th,.plan-doc td{border:1px solid var(--rule);padding:8px 12px;text-align:left;vertical-align:top}
.plan-doc th{background:#F4EFE6;font-family:'JetBrains Mono',monospace;font-size:11px;text-transform:uppercase;letter-spacing:.04em}
.plan-doc code{font-family:'JetBrains Mono',monospace;font-size:13px;background:#F4EFE6;padding:1px 6px;border-radius:3px;color:var(--ink)}
.plan-doc hr{border:0;border-top:1px solid var(--rule);margin:30px 0}
.plan-doc blockquote{border-left:4px solid var(--accent);padding:8px 14px;background:#FCF8EE;margin:12px 0;font-style:italic;color:var(--ink-soft)}
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
  <span class="center">Paper 9 — Strategic Plan</span>
  <a href="portfolio_paper9.html">Paper 9 detail →</a>
</nav>
<div class="plan-doc">
__BODY__
</div>
</div>
</body></html>"""

OUT_HTML.write_text(shell.replace("__BODY__", body_html))
print(f"Wrote {OUT_HTML}  ({len(shell):,} -> {OUT_HTML.stat().st_size:,} bytes)")
