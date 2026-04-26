#!/usr/bin/env python3
"""Idempotent: inject a '리뷰어' nav link pointing to v3_reviewer_guide.html.

Inserts immediately after the '통계' (20_meta_power.html) nav link in the
topnav / drawer / commandIndex on every dashboard HTML page. Safe to re-run:
if the target link already exists on a page, that page is left untouched.

Usage:
    python3 /opt/thyroid-dash/project/scripts/inject_reviewer_navlink.py
"""
import json
import re
from pathlib import Path

ROOT = Path("/opt/thyroid-dash/project/reports/html")
TARGET = "v3_reviewer_guide.html"
LABEL = "리뷰어"
ANCHOR_FILE = "20_meta_power.html"
ANCHOR_LABEL = "통계"

pages = list((ROOT / "pages").glob("*.html")) + [ROOT / "index.html"]


def process(p: Path) -> None:
    if not p.exists():
        return
    text = p.read_text()
    in_pages = p.parent.name == "pages"
    is_self = p.name == TARGET

    # Idempotency: if a link with our target+label already exists, skip.
    if TARGET in text and f">{LABEL}</a>" in text:
        print(f"skip (already has link): {p.name}")
        return

    # Build candidate hrefs for this page context.
    if in_pages:
        self_href = TARGET
    else:
        self_href = f"pages/{TARGET}"

    new_text = text
    patched = False

    # --- TOPNAV with role=menuitem anchored on '통계'
    pat_nav_role = re.compile(
        r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?' + re.escape(ANCHOR_FILE) +
        r'" role="menuitem">' + ANCHOR_LABEL + r'</a>)',
        re.UNICODE,
    )
    pat_nav_plain = re.compile(
        r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?' + re.escape(ANCHOR_FILE) +
        r'">' + ANCHOR_LABEL + r'</a>)',
        re.UNICODE,
    )

    def sub_nav(m):
        cls = "nav-link active" if is_self else "nav-link"
        source = m.group(1)
        if 'href="pages/' in source:
            hh = f"pages/{TARGET}"
        elif 'href="../pages/' in source:
            hh = f"../pages/{TARGET}"
        elif 'href="../' in source:
            hh = f"../{TARGET}"
        else:
            hh = self_href
        if ' role="menuitem"' in source:
            return source + f'<a class="{cls}" href="{hh}" role="menuitem">{LABEL}</a>'
        return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'

    replaced, n = pat_nav_role.subn(sub_nav, new_text, count=1)
    if n:
        new_text = replaced
        patched = True
    else:
        replaced, n = pat_nav_plain.subn(sub_nav, new_text, count=1)
        if n:
            new_text = replaced
            patched = True

    # --- DRAWER (class without nav-link prefix, e.g. class="" or class="active")
    pat_drawer = re.compile(
        r'(<a class="[^"]*" href="(?:\.\./)?(?:pages/)?' + re.escape(ANCHOR_FILE) +
        r'">' + ANCHOR_LABEL + r'</a>)',
        re.UNICODE,
    )

    def sub_drawer(m):
        cls = "active" if is_self else ""
        source = m.group(1)
        # Skip if this was already handled by the topnav regex (class starts with nav-link)
        if 'class="nav-link' in source:
            return source
        if 'href="pages/' in source:
            hh = f"pages/{TARGET}"
        elif 'href="../' in source:
            hh = f"../{TARGET}"
        else:
            hh = self_href
        return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'

    replaced, n = pat_drawer.subn(sub_drawer, new_text, count=1)
    if n and replaced != new_text:
        new_text = replaced
        patched = True

    # --- commandIndex array
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
                    "label": LABEL,
                    "href": self_href,
                    "tags": "page reviewer guide v3 honesty audit",
                }
                idx_after = None
                for i, e in enumerate(arr):
                    if not isinstance(e, dict):
                        continue
                    if ANCHOR_FILE in e.get("href", "") or e.get("label") == ANCHOR_LABEL:
                        idx_after = i
                        break
                if idx_after is None:
                    new_arr = arr + [new_entry]
                else:
                    new_arr = arr[: idx_after + 1] + [new_entry] + arr[idx_after + 1:]
                replacement = (
                    m_cmd.group(1)
                    + json.dumps(new_arr, ensure_ascii=False)
                    + m_cmd.group(3)
                )
                new_text = new_text[: m_cmd.start()] + replacement + new_text[m_cmd.end():]
                patched = True

    if patched and new_text != text:
        p.write_text(new_text)
        print(f"patched: {p.name}")
    else:
        print(f"no change: {p.name}")


if __name__ == "__main__":
    count = 0
    for pp in pages:
        before = pp.read_text() if pp.exists() else ""
        process(pp)
        after = pp.read_text() if pp.exists() else ""
        if before != after:
            count += 1
    print(f"---\nreviewer nav link injected into {count} page(s)")
