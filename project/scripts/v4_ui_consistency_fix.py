#!/usr/bin/env python3
"""
v4 UI consistency auto-fixer for the THYRAI / THCA dashboard.

Idempotent and safe to re-run. Only performs the "obvious, non-controversial" fixes
identified in the v4 UI audit (see logs/v4_ui_audit.md for full findings):

  1. Drop exact duplicate consecutive <a class="nav-link"> tags that point to the
     same href (e.g. the duplicated 투자자 / 연구자 links introduced by an older
     injector script).
  2. For any old-style dashboard page (01-28, 99, v3_reviewer_guide, view_*) that
     loads brand.css but NOT main.css, inject main.css (without removing brand).
     These pages were built for the v2 dark palette and break without main.css.
  3. Add missing title="..." to <iframe> tags. Value is derived from the
     surrounding <h3>/<h4>/card context; falls back to "Interactive figure".
  4. Add loading="lazy" to <img> and <iframe> below the first 40 lines
     ("below-fold"), skipping any tag that already has loading= set or lives in
     the hero / <header>.
  5. Comment out <img> / <iframe> tags whose local src file is missing, wrapping
     them in an HTML comment noting the missing path.

STRICT PARALLEL-SAFETY: this script skips any file whose basename starts with
``29_`` / ``30_`` / ``31_`` / ``32_`` / ``33_`` or ``v4_`` and never touches
index.html. The v4 sprint agent owns those files.

Usage:
    python3 scripts/v4_ui_consistency_fix.py [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path("/opt/thyroid-dash/project/reports/html")
PAGES_DIR = ROOT / "pages"
FIGS_INT = ROOT / "figs_interactive"
FIGS_PRE = ROOT / "figs_preview"

V4_PREFIXES = ("29_", "30_", "31_", "32_", "33_", "v4_")


def is_v4_zone(path: Path) -> bool:
    name = path.name
    if name == "index.html" and path.parent == ROOT:
        return True
    return any(name.startswith(p) for p in V4_PREFIXES)


# ---------------------------------------------------------------------------
# Fix 1: drop exact duplicate consecutive nav anchors
# ---------------------------------------------------------------------------
# We look at ANY pair of <a ...> tags that sit inside the top-nav region and
# point to the same href with the same visible text, only whitespace between.
# The v2/v3 nav-injector scripts sometimes added both a full `class="nav-link"`
# version AND a bare `class=""` version right after it (see 99_glossary,
# view_investor, view_researcher where "투자자"/"연구자" appear twice).
A_ANY_RE = re.compile(
    r'<a\s+class="[^"]*"\s+href="([^"]+)"[^>]*>([^<]*)</a>'
)


def dedupe_consecutive_nav_links(html: str) -> tuple[str, int]:
    """Remove the FIRST of any adjacent <a> pair that has the same href AND the
    same visible text with only whitespace between them (strict "exact
    duplicate consecutive" definition from the audit spec).

    Note: interleaved duplicate nav runs (e.g. inv/res/inv/res) are intentionally
    NOT auto-removed here — those are reported in the findings doc for user
    review because their removal boundaries are ambiguous.
    """
    matches = list(A_ANY_RE.finditer(html))
    if len(matches) < 2:
        return html, 0
    to_remove: list[tuple[int, int]] = []
    i = 0
    while i < len(matches) - 1:
        a, b = matches[i], matches[i + 1]
        between = html[a.end(): b.start()]
        if between.strip() == "" and a.group(1) == b.group(1) and a.group(2).strip() == b.group(2).strip():
            to_remove.append((a.start(), a.end()))
            i += 2  # don't re-use b
        else:
            i += 1
    if not to_remove:
        return html, 0
    out = html
    for s, e in reversed(to_remove):
        m = re.match(r"\s*", out[e:])
        trail = m.end() if m else 0
        out = out[:s] + out[e + trail:]
    return out, len(to_remove)


# ---------------------------------------------------------------------------
# Fix 2: ensure old pages that somehow only load brand.css also load main.css
# ---------------------------------------------------------------------------
OLD_PAGE_RE = re.compile(r"^(0[1-9]|1[0-9]|2[0-8]|99_|v3_reviewer|view_)")


def is_old_page(name: str) -> bool:
    return bool(OLD_PAGE_RE.match(name))


def ensure_main_css_on_old_pages(html: str, page_name: str) -> tuple[str, int]:
    if not is_old_page(page_name):
        return html, 0
    has_brand = 'href="../assets/css/brand.css"' in html
    has_main = 'href="../assets/css/main.css"' in html
    if has_brand and not has_main:
        # insert main.css just before brand.css link
        new_link = '<link rel="stylesheet" href="../assets/css/main.css">\n  '
        html = html.replace(
            '<link rel="stylesheet" href="../assets/css/brand.css">',
            new_link + '<link rel="stylesheet" href="../assets/css/brand.css">',
            1,
        )
        return html, 1
    return html, 0


# ---------------------------------------------------------------------------
# Fix 3: add missing title="" to iframes
# ---------------------------------------------------------------------------
IFRAME_OPEN_RE = re.compile(r"<iframe\b([^>]*)>", re.IGNORECASE)


def _derive_iframe_title(html: str, start_idx: int, src: str | None) -> str:
    """Look back ~800 chars for an <h3>/<h4>/<strong> or card heading."""
    window = html[max(0, start_idx - 800):start_idx]
    for tag in ("h3", "h4", "h2", "strong", "figcaption"):
        m = None
        for m in re.finditer(fr"<{tag}[^>]*>(.*?)</{tag}>", window, re.IGNORECASE | re.DOTALL):
            pass
        if m:
            txt = re.sub(r"<[^>]+>", "", m.group(1)).strip()
            if txt:
                return txt[:120]
    if src:
        base = src.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        base = base.replace("_", " ").replace("-", " ").strip()
        if base:
            return f"Interactive figure: {base}"
    return "Interactive figure"


def add_iframe_titles(html: str) -> tuple[str, int]:
    added = 0
    out_parts: list[str] = []
    last = 0
    for m in IFRAME_OPEN_RE.finditer(html):
        attrs = m.group(1)
        # already has title?
        if re.search(r"\btitle\s*=", attrs, re.IGNORECASE):
            continue
        src_m = re.search(r'(?<![\w-])src\s*=\s*["\']([^"\']*)["\']', attrs, re.IGNORECASE)
        src = src_m.group(1) if src_m else None
        title = _derive_iframe_title(html, m.start(), src)
        title_esc = (
            title.replace("&", "&amp;").replace('"', "&quot;")
        )
        new_attrs = attrs.rstrip() + f' title="{title_esc}"'
        out_parts.append(html[last: m.start()])
        out_parts.append(f"<iframe{new_attrs}>")
        last = m.end()
        added += 1
    if added == 0:
        return html, 0
    out_parts.append(html[last:])
    return "".join(out_parts), added


# ---------------------------------------------------------------------------
# Fix 4: add loading="lazy" to below-fold img/iframe
# ---------------------------------------------------------------------------
IMG_RE = re.compile(r"<img\b([^>]*)>", re.IGNORECASE)


def add_lazy_loading(html: str) -> tuple[str, int]:
    """Add loading="lazy" to <img> and <iframe> tags that do NOT have a loading
    attribute yet AND are not inside the first 40 lines (hero band).
    """
    lines = html.splitlines(keepends=True)
    if len(lines) <= 40:
        return html, 0
    head = "".join(lines[:40])
    tail = "".join(lines[40:])
    added = 0

    def _inject(tag_re: re.Pattern, text: str) -> tuple[str, int]:
        out: list[str] = []
        last = 0
        n_added = 0
        for m in tag_re.finditer(text):
            attrs = m.group(1)
            if re.search(r"\bloading\s*=", attrs, re.IGNORECASE):
                continue
            # skip if inside hero/header (best effort: look at last 200 chars before tag)
            ctx = text[max(0, m.start() - 200): m.start()].lower()
            if 'class="hero' in ctx or "<header" in ctx:
                continue
            tagname = "img" if tag_re is IMG_RE else "iframe"
            new_attrs = attrs.rstrip() + ' loading="lazy"'
            out.append(text[last: m.start()])
            # self-closing for img, not for iframe
            if tagname == "img":
                out.append(f"<img{new_attrs}>")
            else:
                out.append(f"<iframe{new_attrs}>")
            last = m.end()
            n_added += 1
        if n_added == 0:
            return text, 0
        out.append(text[last:])
        return "".join(out), n_added

    tail, a1 = _inject(IMG_RE, tail)
    tail, a2 = _inject(IFRAME_OPEN_RE, tail)
    added = a1 + a2
    if added == 0:
        return html, 0
    return head + tail, added


# ---------------------------------------------------------------------------
# Fix 4b: repair the specific wrong-filename links in v3_reviewer_guide.
# The reviewer-guide page was generated with draft filenames that never
# landed; the actual pages shipped as 20_panel_size.html and 23_bethesda_sim.html.
# This is a targeted, safe rename — both target pages exist on disk.
# ---------------------------------------------------------------------------
V3_REVIEWER_GUIDE_FIXES = {
    "20_panel_size_main.html": "20_panel_size.html",
    "23_bethesda_triage.html": "23_bethesda_sim.html",
}


def fix_v3_reviewer_guide_links(html: str, page_name: str) -> tuple[str, int]:
    if page_name != "v3_reviewer_guide.html":
        return html, 0
    n = 0
    for bad, good in V3_REVIEWER_GUIDE_FIXES.items():
        if bad in html:
            new = html.replace(f'href="{bad}"', f'href="{good}"')
            if new != html:
                n += html.count(f'href="{bad}"') - new.count(f'href="{bad}"')
                html = new
    return html, n


# ---------------------------------------------------------------------------
# Fix 5: comment out img/iframe referencing missing local files
# ---------------------------------------------------------------------------
LOCAL_SRC_RE = re.compile(
    r"<(?P<tag>img|iframe)\b(?P<attrs>[^>]*?)>", re.IGNORECASE
)


def _resolve_local(page_file: Path, src: str) -> Path | None:
    """Return absolute Path for a local src, or None if src is remote / anchor
    or is clearly a JS template placeholder (contains `+`, `${`, or `"+`).
    """
    if not src or src.startswith(("http://", "https://", "//", "mailto:", "data:", "#", "javascript:")):
        return None
    # Skip JS-template srcs like `'+src+'` or `${path}`
    if "+" in src or "${" in src or "{{" in src:
        return None
    src = src.split("?", 1)[0].split("#", 1)[0]
    if src.startswith("/"):
        # served from site root; map "/reports/html/..." -> ROOT/..
        # site root is /, but local file root is ROOT's parent's parent
        # the server rewrites / -> /reports/html/, so absolute paths are relative to site root
        candidate = Path("/opt/thyroid-dash/project") / src.lstrip("/")
        if "reports/html" in str(candidate) and candidate.exists():
            return candidate
        # try stripping /reports/html/ prefix and rooting at ROOT
        stripped = src.replace("/reports/html/", "", 1).lstrip("/")
        return (ROOT / stripped)
    # relative
    return (page_file.parent / src).resolve()


def comment_missing_media(html: str, page_file: Path) -> tuple[str, int, list[str]]:
    missing = 0
    missing_list: list[str] = []
    out: list[str] = []
    last = 0
    for m in LOCAL_SRC_RE.finditer(html):
        attrs = m.group("attrs")
        src_m = re.search(r'(?<![\w-])src\s*=\s*["\']([^"\']+)["\']', attrs)
        if not src_m:
            continue
        src = src_m.group(1)
        resolved = _resolve_local(page_file, src)
        if resolved is None:
            continue
        if resolved.exists():
            continue
        # missing — comment out
        out.append(html[last: m.start()])
        tag_text = m.group(0)
        out.append(f"<!-- v4-audit: missing file {src} -->\n<!-- {tag_text} -->")
        last = m.end()
        missing += 1
        missing_list.append(src)
    if missing == 0:
        return html, 0, []
    out.append(html[last:])
    return "".join(out), missing, missing_list


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def process_file(path: Path, dry_run: bool) -> dict:
    original = path.read_text(encoding="utf-8")
    html = original
    stats = {"file": str(path), "fixes": {}}

    html, n = dedupe_consecutive_nav_links(html)
    if n:
        stats["fixes"]["dedupe_nav"] = n
    html, n = ensure_main_css_on_old_pages(html, path.name)
    if n:
        stats["fixes"]["add_main_css"] = n
    html, n = add_iframe_titles(html)
    if n:
        stats["fixes"]["iframe_title"] = n
    html, n = add_lazy_loading(html)
    if n:
        stats["fixes"]["lazy_loading"] = n
    html, n = fix_v3_reviewer_guide_links(html, path.name)
    if n:
        stats["fixes"]["v3_reviewer_relink"] = n
    html, n, missing = comment_missing_media(html, path)
    if n:
        stats["fixes"]["commented_missing"] = n
        stats["missing_files"] = missing

    if html != original and not dry_run:
        path.write_text(html, encoding="utf-8")
        stats["written"] = True
    else:
        stats["written"] = False
    return stats


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not PAGES_DIR.is_dir():
        print(f"pages dir not found: {PAGES_DIR}", file=sys.stderr)
        return 2

    results: list[dict] = []
    for path in sorted(PAGES_DIR.glob("*.html")):
        if is_v4_zone(path):
            continue
        results.append(process_file(path, args.dry_run))

    total = {k: 0 for k in (
        "dedupe_nav", "add_main_css", "iframe_title", "lazy_loading",
        "v3_reviewer_relink", "commented_missing",
    )}
    files_written = 0
    for r in results:
        if r.get("written"):
            files_written += 1
        for k, v in r.get("fixes", {}).items():
            total[k] = total.get(k, 0) + v

    summary = {
        "dry_run": args.dry_run,
        "files_considered": len(results),
        "files_written": files_written,
        "total_fixes": total,
    }
    # persist for the audit script
    if not args.dry_run:
        log_dir = Path("/opt/thyroid-dash/project/logs")
        log_dir.mkdir(exist_ok=True)
        (log_dir / "v4_ui_fix_stats.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
