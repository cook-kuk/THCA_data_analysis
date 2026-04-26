#!/usr/bin/env python3
"""In-place patch: add 투자자 + 연구자 persona links to topnav, drawer, and commandIndex on every page.

Idempotent: re-running does not duplicate links. Links are inserted BEFORE the 용어 (glossary)
nav link, matching the pattern previously used by inject_drug_navlink.py.
"""
import json
import re
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/reports/html")
pages = list((ROOT / "pages").glob("*.html")) + [ROOT / "index.html"]

PERSONA_LINKS = [
    ("투자자", "view_investor.html"),
    ("연구자", "view_researcher.html"),
]


def process_page(p: Path):
    if not p.exists():
        return
    text = p.read_text()
    in_pages = p.parent.name == "pages"
    page_name = p.name

    # Compute hrefs relative to this page location.
    def href_for(target):
        if in_pages:
            return target  # same directory
        else:
            return f"pages/{target}"

    active_map = {t: (page_name == t) for _, t in PERSONA_LINKS}

    # ---- 1. TOPNAV injection ----
    # Match the 용어 topnav link anchor and insert persona links BEFORE it, if not already present.
    def build_topnav_links():
        parts = []
        for label, target in PERSONA_LINKS:
            href = href_for(target)
            cls = "nav-link active" if active_map[target] else "nav-link"
            parts.append(f'<a class="{cls}" href="{href}" role="menuitem">{label}</a>')
        return "".join(parts)

    topnav_links_html = build_topnav_links()

    # Only inject if not already there (idempotent check for BOTH persona targets)
    has_investor_nav = ('href="view_investor.html"' in text) or ('href="pages/view_investor.html"' in text)
    has_researcher_nav = ('href="view_researcher.html"' in text) or ('href="pages/view_researcher.html"' in text)

    new_text = text
    if not (has_investor_nav and has_researcher_nav):
        # Pattern candidates for 용어 topnav link (role=menuitem variant)
        # We match the full anchor so we can place our links BEFORE it.
        glossary_topnav_patterns = [
            # With role=menuitem
            re.compile(r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?99_glossary\.html" role="menuitem">용어</a>)', re.UNICODE),
            # Without role=menuitem
            re.compile(r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?99_glossary\.html">용어</a>)', re.UNICODE),
        ]

        def topnav_sub(m):
            return topnav_links_html + m.group(1)

        for pat in glossary_topnav_patterns:
            replaced, n = pat.subn(topnav_sub, new_text, count=1)
            if n:
                new_text = replaced
                break

    # ---- 2. DRAWER injection ----
    # Drawer anchors use class "" or "active" (no role) and sit inside #td-drawer.
    has_investor_drawer = re.search(r'<a class="[^"]*" href="(?:pages/)?view_investor\.html">투자자</a>', new_text) is not None
    has_researcher_drawer = re.search(r'<a class="[^"]*" href="(?:pages/)?view_researcher\.html">연구자</a>', new_text) is not None

    if not (has_investor_drawer and has_researcher_drawer):
        drawer_parts = []
        for label, target in PERSONA_LINKS:
            href = href_for(target)
            cls = "active" if active_map[target] else ""
            drawer_parts.append(f'<a class="{cls}" href="{href}">{label}</a>')
        drawer_links_html = "".join(drawer_parts)

        glossary_drawer_patterns = [
            re.compile(r'(<a class="[^"]*" href="(?:\.\./)?(?:pages/)?99_glossary\.html">용어</a>)', re.UNICODE),
        ]

        def drawer_sub(m):
            return drawer_links_html + m.group(1)

        for pat in glossary_drawer_patterns:
            replaced, n = pat.subn(drawer_sub, new_text, count=1)
            if n:
                new_text = replaced
                break

    # ---- 3. commandIndex injection ----
    # Find ThyroidDash.commandIndex = [...] and prepend persona entries if absent.
    # We match the array as JSON, parse, modify, write back.
    cmd_idx_pat = re.compile(
        r'(window\.ThyroidDash\.commandIndex\s*=\s*)(\[[\s\S]*?\])(\s*;)',
        re.UNICODE,
    )
    m_cmd = cmd_idx_pat.search(new_text)
    if m_cmd:
        arr_text = m_cmd.group(2)
        try:
            arr = json.loads(arr_text)
        except Exception:
            arr = None
        if isinstance(arr, list):
            existing_hrefs = {e.get("href", "") for e in arr if isinstance(e, dict)}
            persona_entries = []
            for label, target in PERSONA_LINKS:
                href = href_for(target)
                if href not in existing_hrefs and target not in existing_hrefs:
                    persona_entries.append({"label": label, "href": href, "tags": "page persona"})
            if persona_entries:
                # insert persona entries right after the first entry (index 0, usually 홈)
                new_arr = [arr[0]] + persona_entries + arr[1:] if arr else persona_entries
                replacement = m_cmd.group(1) + json.dumps(new_arr, ensure_ascii=False) + m_cmd.group(3)
                new_text = new_text[:m_cmd.start()] + replacement + new_text[m_cmd.end():]

    if new_text != text:
        p.write_text(new_text)
        print(f"patched: {p.name}")
    else:
        print(f"no change: {p.name}")


if __name__ == "__main__":
    for pp in pages:
        process_page(pp)
