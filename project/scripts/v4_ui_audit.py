#!/usr/bin/env python3
"""
v4 UI consistency audit for the THYRAI / THCA dashboard.

Scans every page under reports/html/pages/ plus the two top-level entries, but
NEVER touches index.html, pages/29_*, pages/30_*, pages/31_*, pages/32_*,
pages/33_*, or any v4_* file (owned by the parallel v4 sprint agent).

Writes:
  logs/v4_ui_audit.md
  logs/v4_ui_audit_stats.json
  logs/v4_ui_audit_screenshots/<page>.png  (6 mobile shots)

Run after scripts/v4_ui_consistency_fix.py so the audit reflects post-fix state.
"""
from __future__ import annotations

import glob
import html as html_mod
import json
import os
import re
import sys
from pathlib import Path
from urllib.parse import urljoin, urlparse

import urllib.request
from playwright.sync_api import sync_playwright

PROJECT = Path("/opt/thyroid-dash/project")
ROOT = PROJECT / "reports/html"
PAGES_DIR = ROOT / "pages"
FIGS_INT = ROOT / "figs_interactive"
FIGS_PRE = ROOT / "figs_preview"
LOGS = PROJECT / "logs"
SHOTS = LOGS / "v4_ui_audit_screenshots"
BASE = "http://127.0.0.1:8012"

V4_PREFIXES = ("29_", "30_", "31_", "32_", "33_", "v4_")


def is_v4_zone(path: Path) -> bool:
    n = path.name
    if n == "index.html" and path.parent == ROOT:
        return True
    return any(n.startswith(p) for p in V4_PREFIXES)


# ---------------------------------------------------------------------------
# STATIC SCAN
# ---------------------------------------------------------------------------
def scan_page(path: Path) -> dict:
    html = path.read_text(encoding="utf-8", errors="replace")
    body_class_m = re.search(r"<body\s[^>]*class=\"([^\"]*)\"", html)
    body_class = body_class_m.group(1) if body_class_m else ""
    has_main = 'href="../assets/css/main.css"' in html or 'href="assets/css/main.css"' in html
    has_brand = 'href="../assets/css/brand.css"' in html or 'href="assets/css/brand.css"' in html
    has_hamburger = ("td-menu-toggle" in html) or ("topnav__hamburger" in html)
    has_triage = "Computational triage only" in html
    has_medical = "Not medical advice" in html
    has_notaffil = "Not affiliated with any drug-discovery" in html

    # nav links
    nav_links: list[tuple[str, str]] = []
    for m in re.finditer(r"<a[^>]+href=\"([^\"]+)\"[^>]*>([^<]*)</a>", html):
        href = m.group(1)
        text = m.group(2).strip()
        # restrict to nav region heuristics: near 'topnav' or 'nav-links' within 2000 chars
        pre = html[max(0, m.start() - 3000): m.start()].lower()
        if "<nav" in pre and "</nav>" not in pre.rsplit("<nav", 1)[-1]:
            nav_links.append((href, text))

    # duplicate nav link detection: same href appears >=2x in nav region
    dup_nav: dict[str, int] = {}
    for href, _ in nav_links:
        dup_nav[href] = dup_nav.get(href, 0) + 1
    dup_nav = {k: v for k, v in dup_nav.items() if v > 1 and not k.startswith(("mailto:", "#", "http"))}

    # figures
    iframes = []
    for m in re.finditer(r"<iframe\b([^>]*)>", html, re.IGNORECASE):
        attrs = m.group(1)
        src_m = re.search(r'(?<![\w-])src\s*=\s*["\']([^"\']+)["\']', attrs)
        data_src_m = re.search(r'data-src\s*=\s*["\']([^"\']+)["\']', attrs)
        src = src_m.group(1) if src_m else (data_src_m.group(1) if data_src_m else None)
        has_title = re.search(r"\btitle\s*=", attrs, re.IGNORECASE) is not None
        has_lazy = re.search(r"\bloading\s*=", attrs, re.IGNORECASE) is not None
        iframes.append({"src": src, "has_title": has_title, "has_lazy": has_lazy})

    imgs = []
    for m in re.finditer(r"<img\b([^>]*)>", html, re.IGNORECASE):
        attrs = m.group(1)
        src_m = re.search(r'(?<![\w-])src\s*=\s*["\']([^"\']+)["\']', attrs)
        src = src_m.group(1) if src_m else None
        has_lazy = re.search(r"\bloading\s*=", attrs, re.IGNORECASE) is not None
        imgs.append({"src": src, "has_lazy": has_lazy})

    # Also record lazy-loaded iframe placeholders (e.g. <div data-src="..."> that
    # the page's JS upgrades to an iframe at scroll-into-view time).
    extra_data_src: list[str] = []
    for m in re.finditer(r'data-src\s*=\s*["\']([^"\']+\.html)["\']', html):
        extra_data_src.append(m.group(1))

    return {
        "file": str(path),
        "name": path.name,
        "body_class": body_class,
        "has_main_css": has_main,
        "has_brand_css": has_brand,
        "has_hamburger": has_hamburger,
        "has_triage": has_triage,
        "has_medical": has_medical,
        "has_notaffil": has_notaffil,
        "nav_links_count": len(nav_links),
        "dup_nav": dup_nav,
        "iframes": iframes,
        "imgs": imgs,
        "extra_data_src": extra_data_src,
    }


def resolve_local(page_file: Path, src: str) -> Path | None:
    if not src or src.startswith(("http://", "https://", "//", "mailto:", "data:", "#", "javascript:")):
        return None
    if "+" in src or "${" in src or "{{" in src:
        return None
    src = src.split("?", 1)[0].split("#", 1)[0]
    if src.startswith("/"):
        cand = PROJECT / src.lstrip("/")
        if cand.exists():
            return cand
        stripped = src.replace("/reports/html/", "", 1).lstrip("/")
        return (ROOT / stripped)
    return (page_file.parent / src).resolve()


def check_figures(scan: list[dict]) -> dict:
    missing = []
    oversized = []
    untitled = []
    referenced: set[Path] = set()
    for pg in scan:
        p = Path(pg["file"])
        for frame in pg["iframes"]:
            src = frame["src"]  # may come from src= OR data-src=
            if not frame["has_title"]:
                untitled.append({"page": pg["name"], "src": src})
            resolved = resolve_local(p, src) if src else None
            if resolved is None:
                continue
            if not resolved.exists():
                missing.append({"page": pg["name"], "src": src})
                continue
            referenced.add(resolved.resolve())
            try:
                sz = resolved.stat().st_size
                if sz > 5 * 1024 * 1024:
                    oversized.append({"file": str(resolved), "bytes": sz, "page": pg["name"]})
                if sz == 0:
                    missing.append({"page": pg["name"], "src": src, "reason": "0 bytes"})
            except OSError:
                pass
        for img in pg["imgs"]:
            src = img["src"]
            resolved = resolve_local(p, src) if src else None
            if resolved is None:
                continue
            if not resolved.exists():
                missing.append({"page": pg["name"], "src": src})
                continue
            referenced.add(resolved.resolve())
            try:
                sz = resolved.stat().st_size
                if sz > 5 * 1024 * 1024:
                    oversized.append({"file": str(resolved), "bytes": sz, "page": pg["name"]})
            except OSError:
                pass
        # lazy-loaded iframes via data-src on other tags
        for src in pg.get("extra_data_src", []):
            resolved = resolve_local(p, src)
            if resolved is None:
                continue
            if not resolved.exists():
                missing.append({"page": pg["name"], "src": src, "kind": "data-src"})
                continue
            referenced.add(resolved.resolve())
            try:
                sz = resolved.stat().st_size
                if sz > 5 * 1024 * 1024:
                    oversized.append({"file": str(resolved), "bytes": sz, "page": pg["name"]})
            except OSError:
                pass

    # orphans: present in figs_interactive/figs_preview but never referenced.
    # SVG/PNG/TSV/CSV files are typically download companions for an HTML
    # figure with the same stem; treat those as referenced when the sibling
    # .html is referenced (so they don't flood the orphan list).
    referenced_stems: set[str] = {p.stem for p in referenced}
    orphans = []
    for d in (FIGS_INT, FIGS_PRE):
        for f in d.glob("*"):
            if not f.is_file():
                continue
            if f.suffix in (".tsv", ".csv"):
                continue
            if f.suffix in (".svg", ".png") and f.stem in referenced_stems:
                continue
            if f.resolve() not in referenced:
                orphans.append(str(f))

    return {"missing": missing, "oversized": oversized, "untitled_iframes": untitled, "orphans": orphans}


# ---------------------------------------------------------------------------
# LINK CRAWL (Playwright)
# ---------------------------------------------------------------------------
def link_crawl(scan: list[dict]) -> list[dict]:
    """HEAD-check every href/src on every page under localhost."""
    to_check: set[tuple[str, str]] = set()  # (page_name, url)
    urls_by_page: dict[str, set[str]] = {}

    def norm(base_url: str, href: str) -> str | None:
        if not href or href.startswith(("mailto:", "#", "javascript:", "data:", "tel:")):
            return None
        u = urljoin(base_url, href)
        if not u.startswith("http://127.0.0.1:8012"):
            # skip external
            if u.startswith(("http://", "https://")):
                return None
            return u
        return u

    for pg in scan:
        name = pg["name"]
        if name == "index.html":
            page_url = f"{BASE}/reports/html/index.html"
        else:
            page_url = f"{BASE}/reports/html/pages/{name}"
        html = Path(pg["file"]).read_text(encoding="utf-8", errors="replace")
        urls: set[str] = set()
        for m in re.finditer(r'(?:href|src)\s*=\s*["\']([^"\']+)["\']', html):
            href = m.group(1)
            # skip template strings
            if "+" in href or "${" in href or "{{" in href:
                continue
            u = norm(page_url, href)
            if u:
                urls.add(u)
        urls_by_page[name] = urls

    # HEAD check
    results: list[dict] = []
    cache: dict[str, int] = {}
    for name, urls in urls_by_page.items():
        for u in sorted(urls):
            code = cache.get(u)
            if code is None:
                try:
                    req = urllib.request.Request(u, method="HEAD")
                    with urllib.request.urlopen(req, timeout=6) as r:
                        code = r.getcode()
                except urllib.error.HTTPError as e:
                    code = e.code
                except Exception:
                    # fall back to GET since some static servers reject HEAD
                    try:
                        req = urllib.request.Request(u, method="GET")
                        with urllib.request.urlopen(req, timeout=6) as r:
                            code = r.getcode()
                    except urllib.error.HTTPError as e:
                        code = e.code
                    except Exception:
                        code = -1
                cache[u] = code
            if code and code >= 400:
                results.append({"page": name, "url": u, "status": code})
    return results


# ---------------------------------------------------------------------------
# MOBILE VIEWPORT (Playwright)
# ---------------------------------------------------------------------------
MOBILE_PAGES = [
    ("index.html", "reports/html/index.html"),
    ("pipeline.html", "reports/html/pages/pipeline.html"),
    ("17_biomarker_insights.html", "reports/html/pages/17_biomarker_insights.html"),
    ("10_gene_explorer.html", "reports/html/pages/10_gene_explorer.html"),
    ("07_ml_baseline.html", "reports/html/pages/07_ml_baseline.html"),
    ("99_glossary.html", "reports/html/pages/99_glossary.html"),
]


def mobile_audit() -> list[dict]:
    SHOTS.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    with sync_playwright() as p:
        browser = p.chromium.launch()
        context = browser.new_context(viewport={"width": 390, "height": 844})
        page = context.new_page()
        console_errors: list[str] = []
        page.on("pageerror", lambda e: console_errors.append(str(e)))
        page.on("console", lambda msg: console_errors.append(f"{msg.type}:{msg.text}") if msg.type == "error" else None)
        for label, path in MOBILE_PAGES:
            errs_before = len(console_errors)
            try:
                page.goto(f"{BASE}/{path}", wait_until="domcontentloaded", timeout=15000)
                page.wait_for_load_state("networkidle", timeout=8000)
            except Exception as e:
                results.append({"page": label, "error": f"navigation: {e}"})
                continue
            overflow = page.evaluate(
                "() => ({scrollW: document.documentElement.scrollWidth, clientW: document.documentElement.clientWidth})"
            )
            has_overflow = overflow["scrollW"] > overflow["clientW"] + 1
            hamburger_visible = False
            try:
                btn = page.query_selector("#td-menu-toggle, .topnav__hamburger, .menu-toggle")
                if btn:
                    hamburger_visible = btn.is_visible()
            except Exception:
                pass
            shot_path = SHOTS / f"{label.replace('/', '_')}.png"
            try:
                page.screenshot(path=str(shot_path), full_page=False)
            except Exception:
                pass
            results.append({
                "page": label,
                "scrollWidth": overflow["scrollW"],
                "clientWidth": overflow["clientW"],
                "horizontal_overflow": has_overflow,
                "hamburger_visible": hamburger_visible,
                "console_errors": console_errors[errs_before:][:10],
                "screenshot": str(shot_path) if shot_path.exists() else None,
            })
        browser.close()
    return results


# ---------------------------------------------------------------------------
# REPORT BUILDER
# ---------------------------------------------------------------------------
def build_report(scan, figure_findings, link_404s, mobile, total_fix_stats) -> str:
    lines: list[str] = []
    lines.append("# v4 UI Consistency Audit")
    lines.append("")
    lines.append(f"- Pages scanned: {len(scan)}")
    lines.append(f"- Total broken links (local 404): {len(link_404s)}")
    lines.append(f"- Missing referenced figures: {len(figure_findings['missing'])}")
    lines.append(f"- Oversized (>5MB) figures: {len(figure_findings['oversized'])}")
    lines.append(f"- Orphan figs (on disk, never referenced): {len(figure_findings['orphans'])}")
    lines.append(f"- Mobile viewports checked: {len(mobile)}")
    lines.append("")
    lines.append("## 1. Navigation coherence")
    thyrai = [s for s in scan if s["body_class"] == "thyrai"]
    v2 = [s for s in scan if s["body_class"] == "dark"]
    light = [s for s in scan if s["body_class"] == "light"]
    empty = [s for s in scan if not s["body_class"]]
    lines.append("")
    lines.append(f"- THYRAI chrome (body.thyrai): {len(thyrai)} pages")
    lines.append(f"- v2 dark chrome (body.dark): {len(v2)} pages")
    lines.append(f"- v2 light chrome (body.light): {len(light)} pages")
    lines.append(f"- Pages missing body class entirely: {len(empty)}  -> {[s['name'] for s in empty]}")
    lines.append("")
    no_hamb = [s["name"] for s in scan if not s["has_hamburger"] and s["body_class"] in ("dark", "light")]
    lines.append(f"- v2 pages without `#td-menu-toggle` or `.menu-toggle` (hamburger): {len(no_hamb)}")
    if no_hamb:
        lines.append("  - " + ", ".join(no_hamb))
    dup_examples = [(s["name"], s["dup_nav"]) for s in scan if s["dup_nav"]]
    lines.append(f"- Pages with duplicate nav `href` (same URL linked ≥2×): {len(dup_examples)}")
    for name, dups in dup_examples:
        lines.append(f"  - `{name}`: {dups}")
    lines.append("")
    lines.append("## 2. Theme inconsistency")
    lines.append("")
    lines.append(f"- Old pages loading brand.css but NOT main.css: see auto-fix table")
    bad = [s["name"] for s in scan if s["has_brand_css"] and not s["has_main_css"] and not s["name"].startswith(("pipeline", "platform_", "press.", "publications", "team.", "v3_reviewer", "view_"))]
    lines.append(f"- v2/v3 pages with brand.css but no main.css: {bad}")
    orphan_theme = [s["name"] for s in scan if not s["has_main_css"] and not s["has_brand_css"]]
    lines.append(f"- Pages with NEITHER main.css nor brand.css (fully standalone): {orphan_theme}")
    lines.append("")
    lines.append("## 3. Figure placement")
    lines.append("")
    lines.append(f"- Missing file references ({len(figure_findings['missing'])}):")
    for m in figure_findings["missing"]:
        lines.append(f"  - `{m['page']}` -> `{m.get('src')}` {m.get('reason','')}")
    lines.append(f"- Oversized figures (>5MB):")
    for m in figure_findings["oversized"]:
        lines.append(f"  - `{m['file']}` ({m['bytes']/1024/1024:.1f} MB)")
    lines.append(f"- Orphans (disk but not referenced, first 30):")
    for o in figure_findings["orphans"][:30]:
        lines.append(f"  - `{o}`")
    if len(figure_findings["orphans"]) > 30:
        lines.append(f"  - ... and {len(figure_findings['orphans']) - 30} more")
    lines.append("")
    lines.append("## 4. Link integrity (local 404s)")
    lines.append("")
    if not link_404s:
        lines.append("- No local 404s found.")
    else:
        for l in link_404s:
            lines.append(f"- `{l['page']}` -> {l['url']} (status {l['status']})")
    lines.append("")
    lines.append("## 5. Mobile viewport (390×844)")
    lines.append("")
    for r in mobile:
        lines.append(f"### {r['page']}")
        if "error" in r:
            lines.append(f"- ERROR: {r['error']}")
            continue
        lines.append(f"- scrollWidth: {r['scrollWidth']}, clientWidth: {r['clientWidth']}")
        lines.append(f"- horizontal_overflow: **{r['horizontal_overflow']}**")
        lines.append(f"- hamburger visible: {r['hamburger_visible']}")
        if r["console_errors"]:
            lines.append(f"- console errors ({len(r['console_errors'])}):")
            for err in r["console_errors"][:5]:
                lines.append(f"  - `{err}`")
        if r.get("screenshot"):
            lines.append(f"- screenshot: `{r['screenshot']}`")
    lines.append("")
    lines.append("## 6. Disclaimer audit")
    lines.append("")
    missing_triage = [s["name"] for s in scan if not s["has_triage"]]
    missing_medical = [s["name"] for s in scan if not s["has_medical"]]
    missing_notaffil = [s["name"] for s in scan if not s["has_notaffil"] and s["body_class"] == "thyrai"]
    lines.append(f"- THYRAI pages missing \"Not affiliated with any drug-discovery\" footer: {missing_notaffil}")
    lines.append(f"- Pages missing \"Computational triage only\": {len(missing_triage)} of {len(scan)}")
    lines.append(f"- Pages missing \"Not medical advice\": {len(missing_medical)} of {len(scan)}")
    lines.append("")
    lines.append("## 7. Version / build badge")
    version_path = ROOT / "_version.json"
    version_json = json.loads(version_path.read_text()) if version_path.exists() else {}
    actual = len(list(PAGES_DIR.glob("*.html")))
    lines.append(f"- `_version.json.total_pages`: {version_json.get('total_pages')}")
    lines.append(f"- actual pages/*.html count: {actual}")
    if version_json.get("total_pages") != actual:
        lines.append(f"- **MISMATCH**: _version.json is stale.")
    lines.append("")
    lines.append("## Auto-fixes applied")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps(total_fix_stats, indent=2))
    lines.append("```")
    lines.append("")
    lines.append("## Issues requiring user decision")
    lines.append("")
    lines.append("- Pages 25-28 are fully standalone (no shared chrome, no main.css, no topnav) — recompose into v2 chrome or accept as-is?")
    lines.append("- v2 dark pages still use Korean nav labels (개요, 데이터셋, ...); THYRAI pages use English — unify?")
    lines.append("- v3 pages 19-24 have their own mini-nav that cross-links only to other v3 pages — merge into main nav?")
    lines.append("- 99_glossary / view_investor / view_researcher have interleaved duplicate nav runs (`투자자 연구자 투자자 연구자`) — remove the duplicate pair?")
    lines.append("- OLD pages (01-28) do not carry the THYRAI \"Not affiliated\" disclaimer — should they?")
    return "\n".join(lines) + "\n"


def main() -> int:
    LOGS.mkdir(exist_ok=True)
    SHOTS.mkdir(exist_ok=True)

    # collect pages to scan
    page_paths = []
    for p in sorted(glob.glob(str(PAGES_DIR / "*.html"))):
        P = Path(p)
        if is_v4_zone(P):
            continue
        page_paths.append(P)
    page_paths.append(ROOT / "index.html")  # for scanning only, not editing

    scan = [scan_page(p) for p in page_paths]
    figure_findings = check_figures(scan)
    link_404s = link_crawl(scan)

    # load auto-fix stats (pre-computed by v4_ui_consistency_fix.py earlier)
    fix_log = LOGS / "v4_ui_fix_stats.json"
    total_fix_stats = json.loads(fix_log.read_text()) if fix_log.exists() else {"note": "not run yet"}

    mobile = mobile_audit()

    report = build_report(scan, figure_findings, link_404s, mobile, total_fix_stats)
    (LOGS / "v4_ui_audit.md").write_text(report, encoding="utf-8")

    stats = {
        "pages_scanned": len(scan),
        "broken_links": len(link_404s),
        "missing_figures": len(figure_findings["missing"]),
        "oversized_figures": len(figure_findings["oversized"]),
        "orphan_figures": len(figure_findings["orphans"]),
        "iframes_without_title": sum(1 for s in scan for f in s["iframes"] if not f["has_title"]),
        "mobile_overflow_pages": [r["page"] for r in mobile if r.get("horizontal_overflow")],
        "mobile_no_hamburger": [r["page"] for r in mobile if r.get("hamburger_visible") is False],
        "theme_counts": {
            "thyrai": sum(1 for s in scan if s["body_class"] == "thyrai"),
            "dark": sum(1 for s in scan if s["body_class"] == "dark"),
            "light": sum(1 for s in scan if s["body_class"] == "light"),
            "empty": sum(1 for s in scan if not s["body_class"]),
        },
        "fix_stats": total_fix_stats,
        "link_404s": link_404s,
        "version_total_pages": json.loads((ROOT / "_version.json").read_text()).get("total_pages") if (ROOT / "_version.json").exists() else None,
        "actual_total_pages": len(list(PAGES_DIR.glob("*.html"))),
    }
    (LOGS / "v4_ui_audit_stats.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps({"ok": True, **{k: v for k, v in stats.items() if k != "link_404s"}}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
