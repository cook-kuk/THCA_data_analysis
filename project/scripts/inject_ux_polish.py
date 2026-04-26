#!/usr/bin/env python3
"""Inject power-user UX features (table export, print-to-PDF, guided tour,
threshold compare mode) into every applicable THCA dashboard HTML page.

Idempotent: re-running makes no changes when markers already exist.

Targets (relative to ``reports/html``):

  - table-export.js   : pages containing ``<table>`` or ``dx-table`` (static
    HTML OR async-rendered via renderGridTable / gene-explorer / view_*).
  - print.css + print-button.js : only for the report-style pages listed in
    ``PRINT_ALLOW``. ``print.css`` is linked with ``media="print"``.
  - tour.css + tour.js: every page. The tour adds a "? 투어" button to the
    existing ``.nav-actions`` cluster and auto-runs once per browser.
  - threshold-compare.js : only on ``07_ml_baseline.html``.

Does NOT modify Python pipeline, templates, or the HTML schema — inserts
are localized ``<link>`` / ``<script>`` tags with stable marker comments.

Usage::

    python scripts/inject_ux_polish.py
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
HTML_ROOT = ROOT / "reports" / "html"
PAGES_DIR = HTML_ROOT / "pages"
INDEX = HTML_ROOT / "index.html"

PRINT_ALLOW = {
    "01_overview.html",
    "07_ml_baseline.html",
    "14_caveats.html",
    "15_drug_discovery.html",
    "16_quantum.html",
    "17_biomarker_insights.html",
    "view_investor.html",
    "view_researcher.html",
}

TABLE_MARKERS = ("<table", "dx-table", "renderGridTable", "dx-table-wrap")


def asset_prefix(page_path: Path) -> str:
    """Return the correct relative prefix ('../' for /pages, '' for root)."""
    return "../" if page_path.parent.name == "pages" else ""


def has_table(text: str) -> bool:
    return any(m in text for m in TABLE_MARKERS)


def inject_once(html: str, marker: str, payload: str, *, before: str) -> tuple[str, bool]:
    """Insert ``payload`` once immediately before ``before`` if ``marker`` is
    absent. Returns (new_html, changed)."""
    if marker in html:
        return html, False
    idx = html.rfind(before)
    if idx < 0:
        return html, False
    return html[:idx] + payload + html[idx:], True


def process_page(path: Path) -> list[str]:
    text = path.read_text()
    original = text
    changes: list[str] = []
    page_name = path.name
    prefix = asset_prefix(path)

    # ---- 1. CSV table export (conditional on table presence) ----
    if has_table(text):
        marker = "<!-- ux:table-export -->"
        payload = (
            f'\n  {marker}<script src="{prefix}assets/js/table-export.js" defer></script>'
        )
        text, changed = inject_once(text, marker, payload, before="</body>")
        if changed:
            changes.append("table-export")

    # ---- 2. Tour (every page) ----
    tour_css_marker = "<!-- ux:tour-css -->"
    tour_css = (
        f'\n  {tour_css_marker}<link rel="stylesheet" href="{prefix}assets/css/tour.css">'
    )
    text, changed = inject_once(text, tour_css_marker, tour_css, before="</head>")
    if changed:
        changes.append("tour.css")

    tour_js_marker = "<!-- ux:tour-js -->"
    tour_js = (
        f'\n  {tour_js_marker}<script src="{prefix}assets/js/tour.js" defer></script>'
    )
    text, changed = inject_once(text, tour_js_marker, tour_js, before="</body>")
    if changed:
        changes.append("tour.js")

    # ---- 3. Print stylesheet + button (allowlist pages) ----
    if page_name in PRINT_ALLOW:
        print_css_marker = "<!-- ux:print-css -->"
        print_css = (
            f'\n  {print_css_marker}<link rel="stylesheet" media="print" '
            f'href="{prefix}assets/css/print.css">'
        )
        text, changed = inject_once(text, print_css_marker, print_css, before="</head>")
        if changed:
            changes.append("print.css")

        print_js_marker = "<!-- ux:print-js -->"
        print_js = (
            f'\n  {print_js_marker}<script src="{prefix}assets/js/print-button.js" defer></script>'
        )
        text, changed = inject_once(text, print_js_marker, print_js, before="</body>")
        if changed:
            changes.append("print-button.js")

    # ---- 4. Threshold compare mode (07 only) ----
    if page_name == "07_ml_baseline.html":
        tc_marker = "<!-- ux:threshold-compare -->"
        tc_payload = (
            f'\n  {tc_marker}<script src="{prefix}assets/js/threshold-compare.js" defer></script>'
        )
        text, changed = inject_once(text, tc_marker, tc_payload, before="</body>")
        if changed:
            changes.append("threshold-compare.js")

    if text != original:
        path.write_text(text)
    return changes


def iter_pages() -> Iterable[Path]:
    yield INDEX
    for p in sorted(PAGES_DIR.glob("*.html")):
        yield p


def main() -> None:
    summary = {
        "table-export": 0,
        "tour.css": 0,
        "tour.js": 0,
        "print.css": 0,
        "print-button.js": 0,
        "threshold-compare.js": 0,
    }
    total = 0
    unchanged = 0
    for page in iter_pages():
        if not page.exists():
            continue
        total += 1
        changes = process_page(page)
        if changes:
            print(f"patched {page.relative_to(HTML_ROOT)}: {', '.join(changes)}")
            for c in changes:
                summary[c] = summary.get(c, 0) + 1
        else:
            unchanged += 1

    print("")
    print("summary:")
    for k, v in summary.items():
        print(f"  {k:<22} {v:>3} pages")
    print(f"  pages visited          {total:>3}")
    print(f"  pages already current  {unchanged:>3}")


if __name__ == "__main__":
    main()
