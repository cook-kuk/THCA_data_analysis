#!/usr/bin/env python3
"""In-place patch: add a nav link to the drug discovery page on every HTML."""
import re
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/reports/html")
pages = list((ROOT / "pages").glob("*.html")) + [ROOT / "index.html"]

def process_page(p):
    if not p.exists():
        return
    text = p.read_text()
    if "15_drug_discovery.html" in text:
        print(f"skip (already has link): {p.name}")
        return
    in_pages = p.parent.name == "pages"
    self_href = "15_drug_discovery.html" if in_pages else "pages/15_drug_discovery.html"
    is_self = p.name == "15_drug_discovery.html"

    state = {"patched": False}
    # Pattern A: topnav with role=menuitem ending on 14_caveats
    patA = re.compile(r'(<a class="nav-link[^"]*" href="(?:\.\./)?pages?/?14_caveats\.html" role="menuitem">주의점</a>)', re.UNICODE)
    patA2 = re.compile(r'(<a class="nav-link[^"]*" href="14_caveats\.html" role="menuitem">주의점</a>)', re.UNICODE)
    # Pattern B: drawer menu, short anchor, ending on 14_caveats
    patB = re.compile(r'(<a class="" href="14_caveats\.html">주의점</a>)', re.UNICODE)
    patB2 = re.compile(r'(<a class="" href="pages/14_caveats\.html">주의점</a>)', re.UNICODE)
    # Pattern C: hand-written topnav in new page itself (won't match since we skip)
    patC = re.compile(r'(<a class="nav-link" href="14_caveats\.html">주의점</a>)', re.UNICODE)
    patC2 = re.compile(r'(<a class="nav-link" href="pages/14_caveats\.html">주의점</a>)', re.UNICODE)
    # Home page index.html: nav-link with pages/
    patD = re.compile(r'(<a class="nav-link[^"]*" href="pages/14_caveats\.html"[^>]*>주의점</a>)', re.UNICODE)

    def insert(match):
        state["patched"] = True
        cls = "nav-link active" if is_self else "nav-link"
        link = f'<a class="{cls}" href="{self_href}" role="menuitem">드럭</a>'
        return match.group(1) + link

    def insert_simple(match):
        state["patched"] = True
        cls = "active" if is_self else ""
        link = f'<a class="{cls}" href="{self_href}">드럭</a>'
        return match.group(1) + link

    def insert_navlink_nomenuitem(match):
        state["patched"] = True
        cls = "nav-link active" if is_self else "nav-link"
        link = f'<a class="{cls}" href="{self_href}">드럭</a>'
        return match.group(1) + link

    new_text = text
    for pat in (patA, patA2, patD):
        new_text = pat.sub(insert, new_text)
    for pat in (patB, patB2):
        new_text = pat.sub(insert_simple, new_text)
    for pat in (patC, patC2):
        new_text = pat.sub(insert_navlink_nomenuitem, new_text)

    if state["patched"] and new_text != text:
        p.write_text(new_text)
        print(f"patched: {p.name}")
    else:
        print(f"not patched (no marker): {p.name}")


for pp in pages:
    process_page(pp)
