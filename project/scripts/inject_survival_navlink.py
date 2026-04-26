#!/usr/bin/env python3
"""In-place idempotent patch: add a '생존' (survival) nav link pointing to
21_survival.html on every dashboard HTML page.

Insertion rule: AFTER the '해석' (18_pathway_immune_meth.html) link. If
that link is absent, the page is skipped. Idempotent — won't double-patch.
"""
import json
import re
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/reports/html")
TARGET = "21_survival.html"
LABEL = "생존"

pages = list((ROOT / "pages").glob("*.html")) + [ROOT / "index.html"]


def process(p: Path):
    if not p.exists():
        return
    text = p.read_text()
    in_pages = p.parent.name == "pages"
    is_self = p.name == TARGET
    href = TARGET if in_pages else f"pages/{TARGET}"

    if TARGET in text and f">{LABEL}</a>" in text:
        print(f"skip (already has link): {p.name}")
        return

    new_text = text

    # --- 1. TOPNAV: after 해석 link
    pat_nav_role = re.compile(
        r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?18_pathway_immune_meth\.html" role="menuitem">해석</a>)',
        re.UNICODE,
    )
    pat_nav_plain = re.compile(
        r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?18_pathway_immune_meth\.html">해석</a>)',
        re.UNICODE,
    )

    def sub_nav(m):
        cls = "nav-link active" if is_self else "nav-link"
        source = m.group(1)
        if 'href="pages/18_pathway_immune_meth.html"' in source:
            hh = f"pages/{TARGET}"
        elif 'href="../pages/18_pathway_immune_meth.html"' in source:
            hh = f"../pages/{TARGET}"
        elif 'href="../18_pathway_immune_meth.html"' in source:
            hh = f"../{TARGET}"
        else:
            hh = href
        if ' role="menuitem"' in source:
            return source + f'<a class="{cls}" href="{hh}" role="menuitem">{LABEL}</a>'
        return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'

    replaced, n = pat_nav_role.subn(sub_nav, new_text, count=1)
    if n:
        new_text = replaced
    else:
        replaced, n = pat_nav_plain.subn(sub_nav, new_text, count=1)
        if n:
            new_text = replaced

    # --- 2. DRAWER
    pat_drawer = re.compile(
        r'(<a class="[^"]*" href="(?:\.\./)?(?:pages/)?18_pathway_immune_meth\.html">해석</a>)',
        re.UNICODE,
    )

    def sub_drawer(m):
        cls = "active" if is_self else ""
        source = m.group(1)
        if 'href="pages/18_pathway_immune_meth.html"' in source:
            hh = f"pages/{TARGET}"
        elif 'href="../18_pathway_immune_meth.html"' in source:
            hh = f"../{TARGET}"
        else:
            hh = href
        return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'

    replaced, n = pat_drawer.subn(sub_drawer, new_text, count=1)
    if n:
        new_text = replaced

    # --- 3. commandIndex
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
            candidate_hrefs = {TARGET, f"pages/{TARGET}", f"../{TARGET}"}
            if not (existing_hrefs & candidate_hrefs):
                new_entry = {
                    "label": "생존",
                    "href": href,
                    "tags": "page survival kaplan-meier cox prognosis",
                }
                idx_after = None
                for i, e in enumerate(arr):
                    if not isinstance(e, dict):
                        continue
                    if "18_pathway_immune_meth.html" in e.get("href", "") or e.get("label") == "해석":
                        idx_after = i
                        break
                new_arr = (arr + [new_entry]) if idx_after is None else arr[: idx_after + 1] + [new_entry] + arr[idx_after + 1:]
                replacement = m_cmd.group(1) + json.dumps(new_arr, ensure_ascii=False) + m_cmd.group(3)
                new_text = new_text[:m_cmd.start()] + replacement + new_text[m_cmd.end():]

    if new_text != text:
        p.write_text(new_text)
        print(f"patched: {p.name}")
    else:
        print(f"no change: {p.name}")


if __name__ == "__main__":
    for pp in pages:
        process(pp)
