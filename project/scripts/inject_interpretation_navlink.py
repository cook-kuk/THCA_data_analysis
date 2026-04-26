#!/usr/bin/env python3
"""In-place idempotent patch: add a '해석' (interpretation) nav link pointing to
18_pathway_immune_meth.html on every dashboard HTML page.

Insertion rule: AFTER the 마커 (17_biomarker_insights.html) link in topnav,
drawer, and commandIndex. If the 마커 link is absent, the page is skipped
(the page is expected to carry it; run inject_biomarker_navlink.py first if
ever needed).

Safe to re-run (idempotent): checks for existing '해석' reference before
inserting.
"""
import json
import re
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/reports/html")
TARGET = "18_pathway_immune_meth.html"
LABEL = "해석"

pages = list((ROOT / "pages").glob("*.html")) + [ROOT / "index.html"]


def process(p: Path):
    if not p.exists():
        return
    text = p.read_text()
    in_pages = p.parent.name == "pages"
    is_self = p.name == TARGET
    href = TARGET if in_pages else f"pages/{TARGET}"

    # Idempotency: if any reference to our target page already exists AND
    # label '해석' is present, do nothing.
    has_target = TARGET in text
    has_label_nav = f'>{LABEL}</a>' in text and TARGET in text

    new_text = text

    # --- 1. TOPNAV: insert after 마커 link
    if not has_label_nav:
        # role=menuitem variant
        pat_nav_role = re.compile(
            r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?17_biomarker_insights\.html" role="menuitem">마커</a>)',
            re.UNICODE,
        )
        pat_nav_plain = re.compile(
            r'(<a class="nav-link[^"]*" href="(?:\.\./)?(?:pages/)?17_biomarker_insights\.html">마커</a>)',
            re.UNICODE,
        )

        def sub_nav(m):
            cls = "nav-link active" if is_self else "nav-link"
            # Determine href flavour from the matched link (keep consistent ./ vs pages/)
            source = m.group(1)
            if 'href="pages/17_biomarker_insights.html"' in source:
                hh = f"pages/{TARGET}"
            elif 'href="../pages/17_biomarker_insights.html"' in source:
                hh = f"../pages/{TARGET}"
            elif 'href="../17_biomarker_insights.html"' in source:
                hh = f"../{TARGET}"
            else:
                hh = href
            if ' role="menuitem"' in source:
                return source + f'<a class="{cls}" href="{hh}" role="menuitem">{LABEL}</a>'
            else:
                return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'

        replaced, n = pat_nav_role.subn(sub_nav, new_text, count=1)
        if n:
            new_text = replaced
        else:
            replaced, n = pat_nav_plain.subn(sub_nav, new_text, count=1)
            if n:
                new_text = replaced

    # --- 2. DRAWER: pattern is <a class="" href="...">마커</a>
    has_drawer_label = re.search(r'<a class="[^"]*" href="[^"]*18_pathway_immune_meth\.html">해석</a>', new_text) is not None
    if not has_drawer_label:
        pat_drawer = re.compile(
            r'(<a class="[^"]*" href="(?:\.\./)?(?:pages/)?17_biomarker_insights\.html">마커</a>)',
            re.UNICODE,
        )

        def sub_drawer(m):
            cls = "active" if is_self else ""
            source = m.group(1)
            if 'href="pages/17_biomarker_insights.html"' in source:
                hh = f"pages/{TARGET}"
            elif 'href="../17_biomarker_insights.html"' in source:
                hh = f"../{TARGET}"
            else:
                hh = href
            return source + f'<a class="{cls}" href="{hh}">{LABEL}</a>'

        replaced, n = pat_drawer.subn(sub_drawer, new_text, count=1)
        if n:
            new_text = replaced

    # --- 3. commandIndex — add entry if not present
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
                cmd_href = href
                new_entry = {
                    "label": "해석",
                    "href": cmd_href,
                    "tags": "page interpretation pathway immune methylation",
                }
                # Place it immediately after the 바이오마커/마커 entry if we can find it;
                # else append to the end.
                idx_after = None
                for i, e in enumerate(arr):
                    if not isinstance(e, dict):
                        continue
                    if "17_biomarker_insights.html" in e.get("href", "") or e.get("label") in ("바이오마커 인사이트", "마커"):
                        idx_after = i
                        break
                if idx_after is None:
                    new_arr = arr + [new_entry]
                else:
                    new_arr = arr[: idx_after + 1] + [new_entry] + arr[idx_after + 1:]
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
