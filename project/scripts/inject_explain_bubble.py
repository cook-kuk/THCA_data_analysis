#!/usr/bin/env python3
"""Inject explain-bubble CSS/JS link + data-explain-id attributes across all pages.
Idempotent via marker comments.

Targets:
  1. <link> to assets/css/explain-bubble.css in <head>
  2. <script> to assets/js/explain-bubble.js before </body>
  3. data-explain-id on h2.section-title, iframe, .figure-wrap, .card, etc.
  4. Ensure matching entry exists in explain_content.json (placeholder if missing)
"""
from __future__ import annotations
import re
import json
from pathlib import Path

ROOT = Path("/opt/thyroid-dash/project/reports/html")
PAGES_DIR = ROOT / "pages"
CONTENT_JSON = ROOT / "assets/data/explain_content.json"

MARKER_ASSETS = "<!-- explain-bubble-assets -->"
MARKER_IDS = "<!-- explain-ids-injected -->"


def load_content() -> dict:
    if CONTENT_JSON.exists():
        return json.loads(CONTENT_JSON.read_text(encoding="utf-8"))
    return {}


def save_content(d: dict) -> None:
    CONTENT_JSON.write_text(json.dumps(d, ensure_ascii=False, indent=2), encoding="utf-8")


def relative_assets_path(page_path: Path) -> str:
    """Return relative path prefix from a page to assets/. Index is 'assets/', pages/* is '../assets/'."""
    if page_path.name == "index.html" and page_path.parent == ROOT:
        return "assets/"
    return "../assets/"


def ensure_assets(html: str, page_path: Path) -> tuple[str, bool]:
    """Add explain-bubble CSS link + JS script to <head> / before </body>. Idempotent."""
    if MARKER_ASSETS in html:
        return html, False
    prefix = relative_assets_path(page_path)
    css_link = f'<link rel="stylesheet" href="{prefix}css/explain-bubble.css">'
    js_script = f'<script src="{prefix}js/explain-bubble.js" defer></script>'
    marker_block_head = f"{MARKER_ASSETS}\n{css_link}\n"
    marker_block_body = f"{MARKER_ASSETS}\n{js_script}\n"

    changed = False
    # Add CSS before </head>
    if "</head>" in html and css_link not in html:
        html = html.replace("</head>", f"{marker_block_head}</head>", 1)
        changed = True
    # Add JS before </body>
    if "</body>" in html and js_script not in html:
        html = html.replace("</body>", f"{marker_block_body}</body>", 1)
        changed = True
    return html, changed


# Regex patterns for target elements
# Match <h2 ... id="xxx" ...> — ANY h2 with an id
PAT_H2_SECTION = re.compile(
    r'(<h2(?![^>]*data-explain-id)[^>]*\bid="([^"]+)"[^>]*)(>)',
    re.IGNORECASE,
)
# Match <section ... id="xxx" ...> — ANY section with id
PAT_SECTION_ID = re.compile(
    r'(<section(?![^>]*data-explain-id)[^>]*\bid="([^"]+)"[^>]*)(>)',
    re.IGNORECASE,
)
# Skip id values that look like chrome (not content sections)
SKIP_IDS = {
    "main-content", "page-hero-title", "command-palette", "figure-modal",
    "figure-modal-title", "td-drawer", "td-drawer-backdrop", "td-menu-toggle",
    "palette-input", "palette-results", "figure-modal-frame", "figure-modal-open",
    "figure-modal-png", "figure-modal-tsv", "figure-modal-copy", "figure-modal-fullscreen",
    "figure-modal-title", "landing-title", "q-hero-title", "bm-hero-title",
    "biz-hero-title", "sec-ml-nav", "sec-ml-task", "sec-ml-curves", "sec-ml-summary",
    "sec-threshold-table", "sec-threshold-explorer",
}


def has_data_explain(tag_open: str) -> bool:
    return "data-explain-id" in tag_open


def inject_explain_ids(html: str, page_slug: str, content: dict) -> tuple[str, int, list[str]]:
    """Inject data-explain-id on section headings. Return new html, n_injected, list of ids."""
    injected_ids: list[str] = []

    def h2_repl(m: re.Match) -> str:
        tag_open, sec_id, close = m.group(1), m.group(2), m.group(3)
        if has_data_explain(tag_open) or sec_id in SKIP_IDS:
            return m.group(0)
        full_id = f"{page_slug}-{sec_id}"
        injected_ids.append(full_id)
        return f'{tag_open} data-explain-id="{full_id}"{close}'

    def section_repl(m: re.Match) -> str:
        tag_open, sec_id, close = m.group(1), m.group(2), m.group(3)
        if has_data_explain(tag_open) or sec_id in SKIP_IDS:
            return m.group(0)
        full_id = f"{page_slug}-{sec_id}"
        # Check if matching h2 id already injected (avoid dup)
        if any(full_id == x for x in injected_ids):
            return m.group(0)
        injected_ids.append(full_id)
        return f'{tag_open} data-explain-id="{full_id}"{close}'

    # Pass 1: h2 with id
    html2 = PAT_H2_SECTION.sub(h2_repl, html)
    # Pass 2: section with id (for pages where h2 has no id but section does)
    html2 = PAT_SECTION_ID.sub(section_repl, html2)
    # Pass 3: h2 WITHOUT id → auto-generate id from text
    PAT_H2_NOID = re.compile(
        r'(<h2(?![^>]*\bid=)(?![^>]*data-explain-id)[^>]*)(>)(.*?)(</h2>)',
        re.IGNORECASE | re.DOTALL,
    )
    counter = [0]
    def h2_noid_repl(m: re.Match) -> str:
        tag_open, close, inner, close_tag = m.group(1), m.group(2), m.group(3), m.group(4)
        # derive slug from visible text
        text = re.sub(r"<[^>]+>", "", inner).strip()
        slug = re.sub(r"[^\w가-힣\-]+", "-", text)[:40].strip("-").lower()
        if not slug:
            counter[0] += 1
            slug = f"h2-{counter[0]}"
        full_id = f"{page_slug}-auto-{slug}"
        if full_id in SKIP_IDS or any(full_id == x for x in injected_ids):
            return m.group(0)
        injected_ids.append(full_id)
        return f'{tag_open} data-explain-id="{full_id}"{close}{inner}{close_tag}'
    html2 = PAT_H2_NOID.sub(h2_noid_repl, html2)
    return html2, len(injected_ids), injected_ids


def ensure_content_entry(content: dict, idx: str, title: str = "") -> bool:
    """Add placeholder if missing. Return True if added."""
    if idx in content:
        return False
    content[idx] = {
        "title": title or idx,
        "tagline": "설명이 곧 추가됩니다",
        "body": f"자동 생성된 플레이스홀더 entry입니다. `explain_content.json`에서 `{idx}` 항목을 업데이트하세요.",
        "tags": ["placeholder"],
    }
    return True


def extract_section_title(html: str, sec_id: str) -> str:
    """Grab the text content of the h2 inside this section."""
    pat = re.compile(
        rf'<h2[^>]*id="{re.escape(sec_id)}"[^>]*>(.*?)</h2>',
        re.IGNORECASE | re.DOTALL,
    )
    m = pat.search(html)
    if not m:
        return ""
    text = re.sub(r"<[^>]+>", "", m.group(1)).strip()
    return text


def process_page(page_path: Path, content: dict, stats: dict) -> None:
    try:
        html = page_path.read_text(encoding="utf-8")
    except Exception as e:
        stats["errors"].append(f"{page_path.name}: read error {e}")
        return
    if "<html" not in html.lower():
        return

    original = html
    page_slug = page_path.stem  # e.g. "07_ml_baseline"
    changed_any = False

    # 1) Assets
    html, changed_assets = ensure_assets(html, page_path)
    if changed_assets:
        changed_any = True
        stats["assets_added"] += 1

    # 2) data-explain-id injection
    html, n_inj, ids = inject_explain_ids(html, page_slug, content)
    if n_inj > 0:
        changed_any = True
        stats["ids_added"] += n_inj
        # Add content entries
        for idx in ids:
            sec_id = idx.split("-", 1)[1]
            title_hint = extract_section_title(html, sec_id) or idx
            if ensure_content_entry(content, idx, title_hint):
                stats["content_entries_added"] += 1
        # Drop idempotency marker at top of <body>
        if MARKER_IDS not in html:
            html = html.replace("<body", f"{MARKER_IDS}\n<body", 1)

    if changed_any:
        page_path.write_text(html, encoding="utf-8")
        stats["pages_touched"] += 1


def main() -> None:
    pages = [ROOT / "index.html"] + sorted(PAGES_DIR.glob("*.html"))
    content = load_content()
    stats = {
        "pages_touched": 0,
        "assets_added": 0,
        "ids_added": 0,
        "content_entries_added": 0,
        "errors": [],
    }
    for p in pages:
        process_page(p, content, stats)
    save_content(content)

    print(f"Pages scanned           : {len(pages)}")
    print(f"Pages touched           : {stats['pages_touched']}")
    print(f"Assets links added      : {stats['assets_added']}")
    print(f"data-explain-id added   : {stats['ids_added']}")
    print(f"Content entries added   : {stats['content_entries_added']}")
    print(f"Total content entries   : {len(content)}")
    if stats["errors"]:
        print("Errors:")
        for e in stats["errors"]:
            print(f"  - {e}")


if __name__ == "__main__":
    main()
