#!/usr/bin/env python3
"""Render the NComms-style manuscript markdown to a journal-formatted PDF.

Uses weasyprint with a Nature-style CSS (sans-serif headings, serif body,
numbered references with superscript callouts, figure legends).
"""
from __future__ import annotations
from pathlib import Path
import re
import sys

try:
    import markdown as md_lib
except ImportError:
    print("installing markdown..."); import subprocess; subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet", "markdown"])
    import markdown as md_lib

from weasyprint import HTML, CSS

ROOT = Path(__file__).resolve().parent.parent
MD_SRC = ROOT / "manuscript" / "ncomms_rai_atlas_v1.md"
OUT_HTML = ROOT / "manuscript" / "ncomms_rai_atlas_v1.html"
OUT_PDF  = ROOT / "manuscript" / "ncomms_rai_atlas_v1.pdf"

NATURE_CSS = """
@page {
  size: A4;
  margin: 18mm 16mm 22mm 16mm;
  @top-left { content: "Cook & Yu  ·  RAI differentiation panel"; font-family: "Liberation Sans","Helvetica","Arial",sans-serif; font-size: 8pt; color: #555; }
  @top-right { content: "Nature Communications  ·  v1 draft  ·  2026-05-21"; font-family: "Liberation Sans","Helvetica","Arial",sans-serif; font-size: 8pt; color: #555; }
  @bottom-center { content: counter(page) " / " counter(pages); font-family: "Liberation Sans","Helvetica","Arial",sans-serif; font-size: 8.5pt; color: #555; }
}
html { font-family: "Liberation Serif","Times New Roman","Nimbus Roman",serif; font-size: 10pt; line-height: 1.5; color: #1c1c1c; }
body { hyphens: auto; text-align: justify; }

h1, h2, h3, h4 { font-family: "Liberation Sans","Helvetica","Arial",sans-serif; color: #0e2a4a; line-height: 1.25; }
h1 { font-size: 17pt; margin: 0 0 6pt 0; font-weight: 700; }
h2 { font-size: 12.5pt; margin: 14pt 0 4pt 0; font-weight: 700; border-bottom: 0.6pt solid #b0bcd0; padding-bottom: 2pt; }
h3 { font-size: 11pt; margin: 10pt 0 3pt 0; font-weight: 600; color: #1f4f88; }
h4 { font-size: 10.5pt; margin: 6pt 0 2pt 0; font-weight: 600; }

p { margin: 4pt 0; }
em { font-style: italic; }
strong { font-weight: 700; color: #0e2a4a; }

/* author block */
h1 + p { font-family: "Liberation Sans","Helvetica","Arial",sans-serif; color: #444; font-size: 10pt; margin-top: 2pt; }
h1 + p + p { font-size: 9.2pt; color: #555; font-style: italic; }

/* abstract emphasis */
h2 + p:first-of-type { font-size: 10pt; padding: 6pt 8pt; background: #f4f6fa; border-left: 2pt solid #1f4f88; }

/* references hanging indent */
ol { padding-left: 0; counter-reset: ref; }
ol li { list-style: none; counter-increment: ref; padding-left: 2.6em; text-indent: -2.6em; margin: 1.5pt 0; font-size: 9pt; line-height: 1.35; }
ol li::before { content: counter(ref) ". "; font-weight: 700; color: #0e2a4a; }

/* figure legends */
strong:first-child { font-family: "Liberation Sans","Helvetica","Arial",sans-serif; color: #0e2a4a; }

/* tables */
table { width: 100%; border-collapse: collapse; margin: 6pt 0 10pt 0; font-size: 9pt; }
th, td { border: 0.4pt solid #aab; padding: 3pt 5pt; vertical-align: top; }
th { background: #e9eff8; color: #0e2a4a; }

code { background: #eef; padding: 0 3pt; border-radius: 2pt; font-size: 8.8pt; font-family: "Liberation Mono","Courier New",monospace; }

hr { border: 0; border-top: 0.6pt solid #b0bcd0; margin: 12pt 0; }

/* List items in main text */
ul { padding-left: 1.4em; }
ul li { margin: 1.5pt 0; }

/* Subtle drop-cap-free Introduction first paragraph */
h2 + p { margin-top: 4pt; }

/* prevent orphans */
p, li, h3, h4 { orphans: 3; widows: 3; }
h2, h3 { break-after: avoid; page-break-after: avoid; }
"""


def main():
    md_text = MD_SRC.read_text()

    # Auto-convert markdown superscript citations like "outcome¹⁻³" — actually those are already unicode
    # Build raw HTML
    html_body = md_lib.markdown(md_text, extensions=["tables", "fenced_code", "sane_lists"])

    # Add minimal wrapper
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <title>RAI atlas — NComms draft</title>
  <style>{NATURE_CSS}</style>
</head>
<body>
{html_body}
</body>
</html>
"""
    OUT_HTML.write_text(html)
    print(f"wrote {OUT_HTML}  ({len(html):,} bytes)")

    # PDF
    HTML(string=html, base_url=str(ROOT)).write_pdf(OUT_PDF)
    print(f"wrote {OUT_PDF}  ({OUT_PDF.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
