#!/usr/bin/env python3
"""In-place patch: add a '통계' nav link to the meta/power/reproducibility
page on every dashboard HTML page. Idempotent — skips pages that already
contain the 20_meta_power link.

Inserts immediately after the '해석' (18_pathway_immune_meth) nav link so
the ordering in the topnav stays monotonic: 마커 → 해석 → 통계 → 투자자.
"""
import re
from pathlib import Path

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/reports/html")
pages = list((ROOT / "pages").glob("*.html")) + [ROOT / "index.html"]


def process_page(p: Path) -> None:
    if not p.exists():
        return
    text = p.read_text()
    if "20_meta_power.html" in text:
        print(f"skip (already has link): {p.name}")
        return

    in_pages = p.parent.name == "pages"
    self_href = "20_meta_power.html" if in_pages else "pages/20_meta_power.html"
    is_self = p.name == "20_meta_power.html"

    state = {"patched": False}

    # Pattern A: topnav inside /pages/ with role=menuitem ending on 18_pathway_immune_meth
    patA = re.compile(
        r'(<a class="nav-link[^"]*" href="18_pathway_immune_meth\.html" role="menuitem">해석</a>)',
        re.UNICODE,
    )
    # Pattern B: topnav outside /pages/ (home page) using pages/ prefix
    patB = re.compile(
        r'(<a class="nav-link[^"]*" href="pages/18_pathway_immune_meth\.html"[^>]*>해석</a>)',
        re.UNICODE,
    )
    # Pattern C: drawer menu short anchors inside /pages/
    patC = re.compile(
        r'(<a class="" href="18_pathway_immune_meth\.html">해석</a>)',
        re.UNICODE,
    )
    # Pattern D: drawer menu with pages/ prefix
    patD = re.compile(
        r'(<a class="" href="pages/18_pathway_immune_meth\.html">해석</a>)',
        re.UNICODE,
    )
    # Pattern E: hand-written nav-link without role=menuitem
    patE = re.compile(
        r'(<a class="nav-link" href="18_pathway_immune_meth\.html">해석</a>)',
        re.UNICODE,
    )
    # Pattern F (fallback): pages with a short nav that ends at 16_quantum.html
    # (the quantum page itself) — append after the 16_quantum nav entry so the
    # topnav still surfaces the link even without a 해석 anchor.
    patF = re.compile(
        r'(<a class="nav-link(?:\s+active)?" href="(?:\.\./)?(?:pages/)?16_quantum\.html" role="menuitem">양자</a>)',
        re.UNICODE,
    )

    def insert_menuitem(match):
        state["patched"] = True
        cls = "nav-link active" if is_self else "nav-link"
        link = f'<a class="{cls}" href="{self_href}" role="menuitem">통계</a>'
        return match.group(1) + link

    def insert_simple(match):
        state["patched"] = True
        cls = "active" if is_self else ""
        link = f'<a class="{cls}" href="{self_href}">통계</a>'
        return match.group(1) + link

    def insert_navlink(match):
        state["patched"] = True
        cls = "nav-link active" if is_self else "nav-link"
        link = f'<a class="{cls}" href="{self_href}">통계</a>'
        return match.group(1) + link

    new_text = text
    for pat in (patA, patB):
        new_text = pat.sub(insert_menuitem, new_text)
    for pat in (patC, patD):
        new_text = pat.sub(insert_simple, new_text)
    new_text = patE.sub(insert_navlink, new_text)
    if not state["patched"]:
        new_text = patF.sub(insert_menuitem, new_text)

    if state["patched"] and new_text != text:
        p.write_text(new_text)
        print(f"patched: {p.name}")
    else:
        print(f"not patched (no anchor found): {p.name}")


if __name__ == "__main__":
    for pp in pages:
        process_page(pp)
