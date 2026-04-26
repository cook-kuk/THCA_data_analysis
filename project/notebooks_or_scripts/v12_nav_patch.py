#!/usr/bin/env python3
"""Inject Timeline / Reasoning / v12 Biology links into the top-nav of pages
that are missing them. Idempotent: skips pages that already have the v12 link."""
import re
from pathlib import Path

PAGES = [
    "platform_omega.html",
    "platform_chem.html",
    "platform_chem_showcase.html",
    "platform_panel.html",
    "platform_clinic.html",
    "v9_interactive_dial.html",
    "v11_novel_pathways.html",
]
ROOT = Path("/opt/thyroid-dash/project/reports/html/pages")

# Pipeline anchor matches both "../pages/pipeline.html" (nested) and "pipeline.html" (flat).
PIPELINE_RX = re.compile(
    r'(<li><a class="topnav__link" href="(?:\.\./pages/)?pipeline\.html">Pipeline</a></li>)',
    re.MULTILINE,
)

def patch(path: Path):
    txt = path.read_text(encoding="utf-8")
    if "v12_biological_evidence" in txt:
        return "skip-already-present"
    m = PIPELINE_RX.search(txt)
    if not m:
        return "skip-no-pipeline-anchor"
    # Detect prefix used: nested pages use "../pages/", flat use ""
    prefix = "../pages/" if "../pages/pipeline.html" in m.group(0) else ""
    inject = (
        f'      <li><a class="topnav__link" href="{prefix}v4a_dial_cross_cancer.html" style="color:#F5A623">v5.1 Audit</a></li>\n'
        f'      <li><a class="topnav__link" href="{prefix}research_timeline.html">Timeline</a></li>\n'
        f'      <li><a class="topnav__link" href="{prefix}research_reasoning.html" style="color:#6AB04C">Reasoning</a></li>\n'
        f'      <li><a class="topnav__link" href="{prefix}v12_biological_evidence.html" style="color:#6AB04C">v12 Biology</a></li>'
    )
    new = txt[:m.end()] + "\n" + inject + txt[m.end():]
    path.write_text(new, encoding="utf-8")
    return "patched"

for fname in PAGES:
    p = ROOT / fname
    if not p.exists():
        print(f"  miss   {fname}")
        continue
    res = patch(p)
    print(f"  {res:8} {fname}")
