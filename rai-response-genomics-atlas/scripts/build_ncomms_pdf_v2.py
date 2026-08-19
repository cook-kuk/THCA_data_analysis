#!/usr/bin/env python3
"""Render NComms-style manuscript markdown to journal-formatted PDF with figures embedded.

Builds both English and Korean v2 + cover letter.

Korean uses Noto Sans/Serif CJK KR fonts (per weasyprint Unifont gotcha — DO NOT
include 'Unifont' in font-family).
"""
from __future__ import annotations
from pathlib import Path
import sys
import re

import markdown as md_lib
from weasyprint import HTML

ROOT = Path(__file__).resolve().parent.parent
MANU = ROOT / "manuscript"


COMMON_BASE_CSS = """
@page {
  size: A4;
  margin: 18mm 16mm 22mm 16mm;
  @bottom-center { content: counter(page) " / " counter(pages); font-size: 8.5pt; color: #555; }
}
body { hyphens: auto; text-align: justify; }

h1, h2, h3, h4 { color: #0e2a4a; line-height: 1.25; }
h1 { font-size: 17pt; margin: 0 0 6pt 0; font-weight: 700; }
h2 { font-size: 12.5pt; margin: 14pt 0 4pt 0; font-weight: 700; border-bottom: 0.6pt solid #b0bcd0; padding-bottom: 2pt; }
h3 { font-size: 11pt; margin: 10pt 0 3pt 0; font-weight: 600; color: #1f4f88; }
h4 { font-size: 10.5pt; margin: 6pt 0 2pt 0; font-weight: 600; }

p { margin: 4pt 0; }
em { font-style: italic; }
strong { font-weight: 700; color: #0e2a4a; }

/* author block */
h1 + p { color: #444; font-size: 10pt; margin-top: 2pt; }
h1 + p + p { font-size: 9.2pt; color: #555; font-style: italic; }

/* abstract emphasis */
h2 + p:first-of-type { font-size: 10pt; padding: 6pt 8pt; background: #f4f6fa; border-left: 2pt solid #1f4f88; }

/* references hanging indent */
ol { padding-left: 0; counter-reset: ref; }
ol li { list-style: none; counter-increment: ref; padding-left: 2.6em; text-indent: -2.6em; margin: 1.5pt 0; font-size: 9pt; line-height: 1.35; }
ol li::before { content: counter(ref) ". "; font-weight: 700; color: #0e2a4a; }

table { width: 100%; border-collapse: collapse; margin: 6pt 0 10pt 0; font-size: 9pt; }
th, td { border: 0.4pt solid #aab; padding: 3pt 5pt; vertical-align: top; }
th { background: #e9eff8; color: #0e2a4a; }

code { background: #eef; padding: 0 3pt; border-radius: 2pt; font-size: 8.8pt; font-family: "Liberation Mono","Courier New",monospace; }
hr { border: 0; border-top: 0.6pt solid #b0bcd0; margin: 12pt 0; }
ul { padding-left: 1.4em; }
ul li { margin: 1.5pt 0; }

img { max-width: 100%; height: auto; display: block; margin: 6pt auto 4pt auto; border: 0.3pt solid #d4dae8; }
img + p, img + em { font-size: 9pt; color: #555; font-style: italic; text-align: center; margin-top: 1pt; }

p, li, h3, h4 { orphans: 3; widows: 3; }
h2, h3 { break-after: avoid; page-break-after: avoid; }

/* keep figures with their captions */
img { page-break-inside: avoid; }
"""


EN_CSS = COMMON_BASE_CSS + """
html { font-family: "Liberation Serif","Times New Roman","Nimbus Roman",serif; font-size: 10pt; line-height: 1.5; color: #1c1c1c; }
h1, h2, h3, h4 { font-family: "Liberation Sans","Helvetica","Arial",sans-serif; }
@page {
  @top-left { content: "Cook & Yu · RAI differentiation panel"; font-family: "Liberation Sans","Helvetica","Arial",sans-serif; font-size: 8pt; color: #555; }
  @top-right { content: "Nature Communications · v2 draft · 2026-05-21"; font-family: "Liberation Sans","Helvetica","Arial",sans-serif; font-size: 8pt; color: #555; }
}
"""

KR_CSS = COMMON_BASE_CSS + """
html { font-family: "Noto Serif CJK KR","Liberation Serif",serif; font-size: 10pt; line-height: 1.55; color: #1c1c1c; word-break: keep-all; }
h1, h2, h3, h4 { font-family: "Noto Sans CJK KR","Liberation Sans",sans-serif; }
body { text-align: left; }   /* CJK looks better left-aligned than full justify */
@page {
  @top-left { content: "Cook & Yu · 갑상선암 RAI 분화 패널"; font-family: "Noto Sans CJK KR"; font-size: 8pt; color: #555; }
  @top-right { content: "Nature Communications · 한글본 v2 · 2026-05-21"; font-family: "Noto Sans CJK KR"; font-size: 8pt; color: #555; }
}
"""

COVER_CSS = EN_CSS  # cover letter uses English style


def build(md_path: Path, css: str, base_url: Path, out_html: Path, out_pdf: Path):
    md_text = md_path.read_text()
    html_body = md_lib.markdown(md_text, extensions=["tables", "fenced_code", "sane_lists"])
    html = f"""<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8"/>
  <title>{md_path.stem}</title>
  <style>{css}</style>
</head>
<body>
{html_body}
</body>
</html>
"""
    out_html.write_text(html)
    HTML(string=html, base_url=str(base_url)).write_pdf(out_pdf)
    print(f"wrote {out_html.name}  ({out_html.stat().st_size:,} B)")
    print(f"wrote {out_pdf.name}  ({out_pdf.stat().st_size:,} B)")


def main():
    # English v2
    build(MANU / "ncomms_rai_atlas_v2.md", EN_CSS, MANU,
          MANU / "ncomms_rai_atlas_v2.html", MANU / "ncomms_rai_atlas_v2.pdf")
    # Korean v2
    build(MANU / "ncomms_rai_atlas_v2_kr.md", KR_CSS, MANU,
          MANU / "ncomms_rai_atlas_v2_kr.html", MANU / "ncomms_rai_atlas_v2_kr.pdf")
    # Korean FINAL (v3)
    final_kr = MANU / "ncomms_rai_atlas_FINAL_kr.md"
    if final_kr.exists():
        build(final_kr, KR_CSS, MANU,
              MANU / "ncomms_rai_atlas_FINAL_kr.html", MANU / "ncomms_rai_atlas_FINAL_kr.pdf")
    # Cover letter
    build(MANU / "cover_letter_v1.md", COVER_CSS, MANU,
          MANU / "cover_letter_v1.html", MANU / "cover_letter_v1.pdf")


if __name__ == "__main__":
    main()
