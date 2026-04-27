#!/usr/bin/env python3
"""Convert all relevant .md, .json, and .tsv files to nicely-styled PDFs AND HTML viewers.

For each input file we now generate two artifacts in pdf_exports:
  - <name>.pdf   — printable PDF (weasyprint)
  - <name>.html  — standalone HTML viewer with toolbar (download original, download PDF, scroll search)

When dashboards link to the .html viewer, clicking shows a pretty table instead of triggering
a raw file download. The toolbar still exposes both download buttons.
"""
from __future__ import annotations
import html as html_lib
import json
import subprocess
import sys
from pathlib import Path

import weasyprint

PROJECT = Path("/opt/thyroid-dash/project")
OUT = PROJECT / "submission" / "npj" / "pdf_exports"
OUT.mkdir(parents=True, exist_ok=True)

SOURCE_DIRS = [
    PROJECT / "submission" / "npj",
    PROJECT / "outreach",
    PROJECT / "results" / "v17_korean",
    PROJECT / "results" / "v17_realfix",
    PROJECT / "results" / "v17_ultimate",
]

CSS = """
@page { size: A4; margin: 1.8cm 1.8cm 2.2cm 1.8cm; @bottom-center { content: counter(page) " / " counter(pages); font-size: 9pt; color: #6b7280; } }
* { box-sizing: border-box; }
html, body { font-family: "Noto Sans CJK KR", "NanumGothic", "Nanum Gothic", "Malgun Gothic", "Apple SD Gothic Neo", "WenQuanYi Zen Hei", "DejaVu Sans", sans-serif; line-height: 1.55; font-size: 10.5pt; color: #1f2937; }
body { max-width: none; padding: 0; margin: 0; }
h1 { font-size: 20pt; color: #0f172a; border-bottom: 3px solid #2563eb; padding-bottom: 6pt; margin: 0 0 14pt 0; line-height: 1.25; }
h2 { font-size: 14.5pt; color: #0f172a; border-bottom: 1px solid #e5e7eb; padding-bottom: 3pt; margin-top: 16pt; line-height: 1.3; page-break-after: avoid; }
h3 { font-size: 12pt; color: #1e293b; margin-top: 11pt; line-height: 1.3; page-break-after: avoid; }
h4 { font-size: 11pt; color: #334155; margin-top: 9pt; }
p { margin: 6pt 0; }
code { background: #f1f5f9; padding: 1pt 4pt; border-radius: 3pt; font-size: 9.2pt; font-family: "DejaVu Sans Mono", "Courier New", monospace; word-break: break-all; }
pre { background: #f8fafc; padding: 9pt 11pt; border-radius: 5pt; border-left: 3px solid #94a3b8; overflow-x: auto; font-size: 9pt; line-height: 1.45; white-space: pre-wrap; word-wrap: break-word; }
pre code { background: transparent; padding: 0; }
table { border-collapse: collapse; width: 100%; margin: 8pt 0; font-size: 9.7pt; page-break-inside: avoid; }
th, td { border: 1px solid #d1d5db; padding: 4pt 8pt; text-align: left; vertical-align: top; }
th { background: #f3f4f6; font-weight: 600; color: #0f172a; }
blockquote { border-left: 3px solid #2563eb; margin: 8pt 0; padding: 4pt 12pt; background: #eff6ff; color: #1e3a8a; font-style: italic; }
a { color: #2563eb; text-decoration: none; word-break: break-all; }
ul, ol { padding-left: 22pt; margin: 6pt 0; }
li { margin: 2pt 0; }
hr { border: none; border-top: 1px solid #e5e7eb; margin: 12pt 0; }
strong, b { color: #0f172a; font-weight: 700; }
em, i { color: #334155; }
img { max-width: 100%; height: auto; }
.title-banner { background: linear-gradient(135deg, #2563eb 0%, #4f46e5 100%); color: white; padding: 14pt 18pt; border-radius: 8pt; margin-bottom: 14pt; }
.title-banner .doc-name { font-size: 16pt; font-weight: 700; margin: 0; }
.title-banner .doc-meta { font-size: 9.5pt; opacity: 0.85; margin-top: 4pt; }
.json-table th { background: #e0e7ff; color: #1e3a8a; }
.json-table.nested { margin: 0; }
.json-key { font-family: "DejaVu Sans Mono", monospace; color: #6d28d9; font-size: 9.2pt; word-break: break-all; }
.json-num { color: #b91c1c; font-variant-numeric: tabular-nums; text-align: right; }
.json-str { color: #166534; }
.json-bool { color: #c2410c; font-weight: 600; }
.json-null { color: #94a3b8; font-style: italic; }
.json-list-scalar { font-family: "DejaVu Sans Mono", monospace; font-size: 9.2pt; color: #1e293b; }
.tsv-table { font-size: 9.5pt; }
.tsv-table th { background: #dbeafe; color: #1e3a8a; font-size: 9.4pt; }
.tsv-table td.num { font-variant-numeric: tabular-nums; text-align: right; font-family: "DejaVu Sans Mono", monospace; font-size: 9.2pt; }
.tsv-table td.metric-cell { background: #fef3c7; color: #92400e; font-weight: 700; }
.tsv-table th.metric-col { background: #fde68a; color: #78350f; }
.tsv-summary { background: #f8fafc; border: 1px solid #cbd5e1; padding: 6pt 10pt; border-radius: 5pt; margin: 8pt 0; font-size: 9.5pt; color: #334155; }
.tsv-summary b { color: #0f172a; }
"""


def html_escape(s) -> str:
    return html_lib.escape(str(s))


VIEWER_CSS = CSS + """
.toolbar { position: sticky; top: 0; background: white; border-bottom: 2px solid #e5e7eb; padding: 10px 14px; margin-bottom: 12px; display: flex; align-items: center; gap: 10px; flex-wrap: wrap; z-index: 10; box-shadow: 0 2px 6px rgba(0,0,0,0.04); }
.toolbar a, .toolbar button { background: #2563eb; color: white; padding: 7px 14px; border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 13px; border: none; cursor: pointer; font-family: inherit; }
.toolbar a:hover, .toolbar button:hover { background: #1d4ed8; }
.toolbar a.secondary { background: white; color: #1e3a8a; border: 1.5px solid #1e3a8a; }
.toolbar a.secondary:hover { background: #eff6ff; }
.toolbar .breadcrumb { color: #475569; font-size: 13px; flex-grow: 1; padding-left: 6px; word-break: break-all; }
.toolbar input[type="search"] { padding: 6px 10px; border: 1px solid #cbd5e1; border-radius: 6px; font-size: 13px; min-width: 180px; }
.viewer-content { padding: 0 24px 40px; max-width: 1280px; margin: 0 auto; }
@media print { .toolbar { display: none; } .viewer-content { padding: 0; } }
"""

VIEWER_SEARCH_JS = """<script>
(function() {
  const input = document.getElementById('row-search');
  if (!input) return;
  input.addEventListener('input', function() {
    const q = this.value.toLowerCase();
    document.querySelectorAll('table tbody tr').forEach(tr => {
      tr.style.display = (!q || tr.textContent.toLowerCase().includes(q)) ? '' : 'none';
    });
  });
})();
</script>"""


def viewer_html(title: str, body: str, src_url: str, pdf_url: str, kind: str, rel: str,
                show_search: bool = False) -> str:
    search = (
        '<input id="row-search" type="search" placeholder="filter rows…">'
        if show_search else ""
    )
    toolbar = (
        f'<div class="toolbar">'
        f'<span class="breadcrumb"><b>{html_escape(kind)}</b> · {html_escape(rel)}</span>'
        f'{search}'
        f'<a href="{html_escape(src_url)}" download>📥 Download original</a>'
        f'<a class="secondary" href="{html_escape(pdf_url)}" target="_blank">📕 PDF</a>'
        f'</div>'
    )
    return (
        "<!DOCTYPE html><html lang='ko'><head><meta charset='utf-8'>"
        f"<title>{html_escape(title)}</title>"
        f"<style>{VIEWER_CSS}</style></head><body>"
        f"{toolbar}"
        f'<div class="viewer-content">{body}</div>'
        f"{VIEWER_SEARCH_JS if show_search else ''}"
        "</body></html>"
    )


def write_viewer(out_html: Path, title: str, body: str, src_path: Path, out_pdf: Path, kind: str,
                 show_search: bool = False):
    rel = src_path.relative_to(PROJECT)
    # Relative URLs from out_html dir → src_path and out_pdf
    src_url = relpath_url(out_html, src_path)
    pdf_url = relpath_url(out_html, out_pdf)
    html = viewer_html(title, body, src_url, pdf_url, kind, str(rel), show_search=show_search)
    out_html.write_text(html, encoding="utf-8")


def relpath_url(from_file: Path, to_file: Path) -> str:
    import os
    return os.path.relpath(to_file, from_file.parent)


def render_value(v, depth=0):
    if v is None:
        return '<span class="json-null">null</span>'
    if isinstance(v, bool):
        return f'<span class="json-bool">{str(v).lower()}</span>'
    if isinstance(v, (int, float)):
        if isinstance(v, float):
            return f'<span class="json-num">{v:.6g}</span>'
        return f'<span class="json-num">{v}</span>'
    if isinstance(v, str):
        return f'<span class="json-str">{html_escape(v)}</span>'
    if isinstance(v, list):
        return render_list(v, depth + 1)
    if isinstance(v, dict):
        return render_dict(v, depth + 1)
    return html_escape(repr(v))


def render_list(lst, depth=0):
    if not lst:
        return '<span class="json-null">[ ]</span>'
    if all(not isinstance(x, (dict, list)) for x in lst):
        items = ", ".join(render_value(x, depth) for x in lst)
        return f'<span class="json-list-scalar">[ {items} ]</span>'
    if all(isinstance(x, dict) for x in lst):
        all_keys = []
        seen = set()
        for d in lst:
            for k in d.keys():
                if k not in seen:
                    seen.add(k)
                    all_keys.append(k)
        rows = []
        rows.append("<tr>" + "".join(f'<th>{html_escape(k)}</th>' for k in all_keys) + "</tr>")
        for d in lst:
            rows.append(
                "<tr>"
                + "".join(f"<td>{render_value(d.get(k), depth)}</td>" for k in all_keys)
                + "</tr>"
            )
        return f'<table class="json-table nested">{"".join(rows)}</table>'
    items = "".join(f"<li>{render_value(x, depth)}</li>" for x in lst)
    return f"<ol>{items}</ol>"


def render_dict(d, depth=0):
    if not d:
        return '<span class="json-null">{ }</span>'
    rows = []
    for k, v in d.items():
        rows.append(
            f'<tr><td class="json-key">{html_escape(k)}</td>'
            f"<td>{render_value(v, depth)}</td></tr>"
        )
    return (
        '<table class="json-table nested">'
        '<tr><th style="width:30%">key</th><th>value</th></tr>'
        f'{"".join(rows)}</table>'
    )


def json_to_pdf(json_path: Path, out_pdf: Path):
    rel = json_path.relative_to(PROJECT)
    title = json_path.stem.replace("_", " ")
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception as e:
        print(f"  ! json parse failed on {json_path.name}: {e}")
        return False

    if isinstance(data, list):
        body = render_list(data)
    elif isinstance(data, dict):
        body = render_dict(data)
    else:
        body = f"<pre>{html_escape(json.dumps(data, indent=2, ensure_ascii=False))}</pre>"

    banner = (
        f'<div class="title-banner">'
        f'<div class="doc-name">📊 {title}</div>'
        f'<div class="doc-meta">{rel} · JSON 데이터 · 2026-04-27</div>'
        f'</div>'
    )
    full_html = (
        f"<!DOCTYPE html><html lang='ko'><head><meta charset='utf-8'>"
        f"<title>{title}</title><style>{CSS}</style></head>"
        f"<body>{banner}{body}</body></html>"
    )
    try:
        weasyprint.HTML(string=full_html, base_url=str(json_path.parent)).write_pdf(str(out_pdf))
        write_viewer(out_pdf.with_suffix(".html"), title, banner + body,
                     json_path, out_pdf, "JSON")
        return True
    except Exception as e:
        print(f"  ! weasyprint failed on {json_path.name}: {e}")
        return False


METRIC_KEYS = {"auc", "p_value", "p", "pvalue", "p-value", "hr", "ci_lower", "ci_upper",
               "auc_ci_lower", "auc_ci_upper", "n", "n_dm1", "n_dm2", "p_dm2",
               "concordance", "stat", "logrank_p", "cox_hr"}


def is_metric_col(name: str) -> bool:
    return name.strip().lower() in METRIC_KEYS or any(
        k in name.strip().lower() for k in ("auc", "p_value", "p_val", "_pval", "_hr", "logrank", "cohen")
    )


def looks_numeric(v: str) -> bool:
    if v is None or v == "":
        return False
    try:
        float(v)
        return True
    except ValueError:
        return False


def fmt_num(v: str) -> str:
    try:
        f = float(v)
        if abs(f) >= 10000 or (0 < abs(f) < 0.0001):
            return f"{f:.3e}"
        if f == int(f):
            return str(int(f))
        return f"{f:.4g}"
    except ValueError:
        return v


def tsv_to_pdf(tsv_path: Path, out_pdf: Path):
    rel = tsv_path.relative_to(PROJECT)
    title = tsv_path.stem.replace("_", " ")
    try:
        text = tsv_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        print(f"  ! tsv read failed on {tsv_path.name}: {e}")
        return False
    lines = text.rstrip("\n").split("\n")
    if not lines:
        return False
    sep = "\t" if "\t" in lines[0] else ","
    headers = lines[0].split(sep)
    data = [ln.split(sep) for ln in lines[1:] if ln.strip()]
    n_rows = len(data)

    metric_idx = {i for i, h in enumerate(headers) if is_metric_col(h)}

    # render table
    head_cells = "".join(
        f'<th class="metric-col">{html_escape(h)}</th>' if i in metric_idx else f'<th>{html_escape(h)}</th>'
        for i, h in enumerate(headers)
    )
    body_rows = []
    for row in data[:300]:  # cap at 300 rows for PDF readability
        cells = []
        for i, v in enumerate(row):
            css = []
            if looks_numeric(v):
                css.append("num")
                v_disp = fmt_num(v)
            else:
                v_disp = v
            if i in metric_idx and looks_numeric(v):
                css.append("metric-cell")
            cls = f' class="{" ".join(css)}"' if css else ""
            cells.append(f"<td{cls}>{html_escape(v_disp)}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")

    truncated_note = ""
    if n_rows > 300:
        truncated_note = (
            f'<div class="tsv-summary">⚠ Showing first 300 of <b>{n_rows}</b> rows '
            f'(full TSV: <code>{rel}</code>)</div>'
        )

    summary = (
        f'<div class="tsv-summary"><b>📊 {n_rows}</b> rows × <b>{len(headers)}</b> cols. '
        f'highlight columns: {", ".join(html_escape(headers[i]) for i in sorted(metric_idx)) or "(none auto-detected)"}</div>'
    )

    table = (
        f'<table class="json-table tsv-table">'
        f'<thead><tr>{head_cells}</tr></thead>'
        f'<tbody>{"".join(body_rows)}</tbody></table>'
    )

    banner = (
        f'<div class="title-banner">'
        f'<div class="doc-name">📑 {title}</div>'
        f'<div class="doc-meta">{rel} · TSV 데이터 · 2026-04-27</div>'
        f'</div>'
    )
    full_html = (
        f"<!DOCTYPE html><html lang='ko'><head><meta charset='utf-8'>"
        f"<title>{title}</title><style>{CSS}</style></head>"
        f"<body>{banner}{summary}{truncated_note}{table}</body></html>"
    )
    try:
        weasyprint.HTML(string=full_html, base_url=str(tsv_path.parent)).write_pdf(str(out_pdf))
        write_viewer(out_pdf.with_suffix(".html"), title,
                     banner + summary + truncated_note + table,
                     tsv_path, out_pdf, "TSV", show_search=(n_rows > 10))
        return True
    except Exception as e:
        print(f"  ! weasyprint failed on {tsv_path.name}: {e}")
        return False


def md_to_pdf(md_path: Path, out_pdf: Path):
    rel = md_path.relative_to(PROJECT)
    title = md_path.stem.replace("_", " ")
    try:
        body = subprocess.check_output(
            ["pandoc", str(md_path),
             "-f", "markdown+yaml_metadata_block+pipe_tables+raw_html",
             "-t", "html5", "--no-highlight"],
            stderr=subprocess.PIPE,
        ).decode("utf-8")
    except subprocess.CalledProcessError as e:
        print(f"  ! pandoc failed on {md_path.name}: {e.stderr.decode('utf-8')[:200]}")
        return False

    banner = (
        f'<div class="title-banner">'
        f'<div class="doc-name">{title}</div>'
        f'<div class="doc-meta">{rel} · 2026-04-27</div>'
        f'</div>'
    )
    full_html = (
        f"<!DOCTYPE html><html lang='ko'><head><meta charset='utf-8'>"
        f"<title>{title}</title><style>{CSS}</style></head>"
        f"<body>{banner}{body}</body></html>"
    )
    try:
        weasyprint.HTML(string=full_html, base_url=str(md_path.parent)).write_pdf(str(out_pdf))
        write_viewer(out_pdf.with_suffix(".html"), title, banner + body,
                     md_path, out_pdf, "Markdown")
        return True
    except Exception as e:
        print(f"  ! weasyprint failed on {md_path.name}: {e}")
        return False


def main():
    md_files = []
    json_files = []
    tsv_files = []
    for d in SOURCE_DIRS:
        if not d.exists():
            continue
        md_files.extend(d.rglob("*.md"))
        json_files.extend(d.rglob("*.json"))
        tsv_files.extend(d.rglob("*.tsv"))
    md_files = sorted(set(md_files))
    json_files = sorted(set(json_files))
    tsv_files = sorted(set(tsv_files))
    # Skip massive auto-generated files (e.g. >2 MB)
    json_files = [j for j in json_files if j.stat().st_size < 2_000_000]
    tsv_files = [t for t in tsv_files if t.stat().st_size < 1_500_000]
    print(f"  {len(md_files)} .md, {len(json_files)} .json, {len(tsv_files)} .tsv found")

    ok = 0
    fail = 0
    for md in md_files:
        rel = md.relative_to(PROJECT)
        out_pdf = OUT / rel.with_suffix(".pdf")
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        if md_to_pdf(md, out_pdf):
            ok += 1
            print(f"  ✓ md   {rel}")
        else:
            fail += 1

    for jp in json_files:
        rel = jp.relative_to(PROJECT)
        out_pdf = OUT / rel.with_suffix(".pdf")
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        if json_to_pdf(jp, out_pdf):
            ok += 1
            print(f"  ✓ json {rel}")
        else:
            fail += 1

    for tp in tsv_files:
        rel = tp.relative_to(PROJECT)
        out_pdf = OUT / rel.with_suffix(".pdf")
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        if tsv_to_pdf(tp, out_pdf):
            ok += 1
            print(f"  ✓ tsv  {rel}")
        else:
            fail += 1

    print(f"\n  ✓ {ok} converted, {fail} failed")
    print(f"  output dir: {OUT}")


if __name__ == "__main__":
    sys.exit(main())
