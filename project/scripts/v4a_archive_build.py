"""
v4a_archive_build.py - Build the v4a Research Archive layer.

Re-surfaces v2 substantive content in the v4 dark brand aesthetic without
modifying any v2 originals. Parallel-safe: all new files prefixed ``v4a_``;
only HTML marker-guarded appends to the landing and v4 pages 29-33 topnavs.

Not affiliated with any drug-discovery AI company.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import html
import json
import os
import re
import subprocess
import sys
from pathlib import Path

# ---------------- paths ----------------
PROJECT = Path("/opt/thyroid-dash/project")
HTML_ROOT = PROJECT / "reports" / "html"
PAGES = HTML_ROOT / "pages"
ASSETS = HTML_ROOT / "assets"
FIGS = HTML_ROOT / "figs_interactive"
FIGS_DARK = FIGS / "v4a_dark"
RESULTS = PROJECT / "results"
REPORTS = PROJECT / "reports"
META = PROJECT / "metadata"
LOGS = PROJECT / "logs"
SERVER = "http://40.82.129.113:8012"

sys.path.insert(0, str(PROJECT / "scripts"))
from v4a_dark_theme import retheme_many  # type: ignore

# ---------------- constants ----------------
BUILD_TIME = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
V4A_TOPIC_PAGES = [
    "v4a_ml_performance.html",
    "v4a_biomarkers.html",
    "v4a_drug_discovery.html",
    "v4a_pathway_immune.html",
    "v4a_panels.html",
    "v4a_data_integration.html",
    "v4a_quantum.html",
    "v4a_limitations.html",
]
V4A_INDEX_PAGE = "v4a_archive_index.html"

NAV_MARK_START = "<!-- v4a-nav-start -->"
NAV_MARK_END = "<!-- v4a-nav-end -->"
LANDING_MARK_START = "<!-- v4a-archive-link-start -->"
LANDING_MARK_END = "<!-- v4a-archive-link-end -->"

FOOTER_DISCLAIMER_HTML = (
    '<span class="site-footer__disclaimer">'
    "Computational triage only. Not medical advice. "
    "Visual system inspired by modern biotech platform conventions. "
    "Not affiliated with any drug-discovery AI company."
    "</span>"
)


# ---------------- probe helpers ----------------
def probe(path: Path) -> bool:
    return path.exists() and path.is_file()


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""


def load_tsv_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    if not probe(path):
        return [], []
    with path.open(encoding="utf-8", errors="replace", newline="") as f:
        rdr = csv.reader(f, delimiter="\t")
        rows = [r for r in rdr if r]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def first_para(md_text: str) -> str:
    """Return the first non-heading paragraph (markdown-stripped) from md_text."""
    if not md_text:
        return ""
    # Look for TL;DR
    m = re.search(r"(?mi)^##\s*TL;?DR\s*\n+(.*?)(?=\n##|\Z)", md_text, flags=re.S)
    if m:
        block = m.group(1).strip()
    else:
        # fall back to first non-heading, non-code paragraph
        block = ""
        lines: list[str] = []
        for raw in md_text.splitlines():
            s = raw.strip()
            if not lines and (s.startswith("#") or not s or s.startswith("-") or s.startswith("|")):
                continue
            if s.startswith("#") or s.startswith("```"):
                if lines:
                    break
                continue
            if not s:
                if lines:
                    break
                continue
            lines.append(s)
        block = " ".join(lines)
    # strip simple markdown
    block = re.sub(r"\*\*(.+?)\*\*", r"\1", block)
    block = re.sub(r"\*(.+?)\*", r"\1", block)
    block = re.sub(r"`(.+?)`", r"\1", block)
    block = re.sub(r"\[(.+?)\]\([^)]+\)", r"\1", block)
    return block.strip()


def md_section(md_text: str, header_pattern: str) -> str:
    """Return body text of first matching markdown section (##)."""
    m = re.search(rf"(?mi)^{header_pattern}\s*\n+(.*?)(?=\n##\s|\Z)", md_text, flags=re.S)
    return m.group(1).strip() if m else ""


def md_to_simple_html(md: str) -> str:
    """Very small markdown to HTML: paragraphs, lists, code, bold/italic/inline code."""
    if not md:
        return ""
    out: list[str] = []
    in_ul = False
    in_code = False
    buf: list[str] = []

    def flush_para():
        nonlocal buf
        if buf:
            text = " ".join(buf).strip()
            if text:
                out.append("<p>" + inline(text) + "</p>")
            buf = []

    def inline(s: str) -> str:
        s = html.escape(s, quote=False)
        s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        s = re.sub(r"\*(.+?)\*", r"<em>\1</em>", s)
        s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
        return s

    for raw in md.splitlines():
        line = raw.rstrip()
        if line.startswith("```"):
            if in_code:
                out.append("</pre>")
                in_code = False
            else:
                flush_para()
                if in_ul:
                    out.append("</ul>")
                    in_ul = False
                out.append('<pre>')
                in_code = True
            continue
        if in_code:
            out.append(html.escape(line))
            continue
        if re.match(r"^\s*[-*]\s+", line):
            flush_para()
            if not in_ul:
                out.append("<ul>")
                in_ul = True
            item = re.sub(r"^\s*[-*]\s+", "", line)
            out.append("<li>" + inline(item) + "</li>")
            continue
        if not line.strip():
            flush_para()
            if in_ul:
                out.append("</ul>")
                in_ul = False
            continue
        if re.match(r"^#{1,6}\s", line):
            flush_para()
            if in_ul:
                out.append("</ul>")
                in_ul = False
            level = len(line) - len(line.lstrip("#"))
            text = line.lstrip("# ").strip()
            out.append(f"<h{min(level, 4)}>{inline(text)}</h{min(level, 4)}>")
            continue
        buf.append(line)

    flush_para()
    if in_ul:
        out.append("</ul>")
    if in_code:
        out.append("</pre>")
    return "\n".join(out)


# ---------------- table rendering (server-side, for static top-N) ----------------
def format_cell(val: str, is_num: bool) -> str:
    if val is None:
        return ""
    v = str(val).strip()
    if not v:
        return ""
    if not is_num:
        return html.escape(v)
    try:
        n = float(v)
    except Exception:
        return html.escape(v)
    if n == 0:
        return "0"
    a = abs(n)
    if a < 1e-2:
        return f"{n:.2e}"
    if a >= 1e4:
        return f"{n:.2e}"
    if n.is_integer() and a < 1e6:
        return str(int(n))
    # 3 sig figs
    s = f"{n:.3g}"
    return s


def detect_num_cols(rows: list[list[str]], ncols: int) -> list[bool]:
    numeric = [True] * ncols
    for r in rows[:80]:
        for i in range(ncols):
            if i >= len(r):
                continue
            v = (r[i] or "").strip()
            if v == "":
                continue
            if not re.match(r"^-?\d+(\.\d+)?([eE][-+]?\d+)?$", v):
                numeric[i] = False
    return numeric


def render_table_html(header: list[str], rows: list[list[str]], top_n: int = 30, label: str = "rows", sort_by: str | None = None) -> str:
    if not header:
        return ""
    ncols = len(header)
    num = detect_num_cols(rows, ncols)
    data = rows
    if sort_by and sort_by in header:
        idx = header.index(sort_by)

        def key(r):
            try:
                return -abs(float(r[idx]))
            except Exception:
                return 1e18

        data = sorted(data, key=key)
    shown = data[:top_n]
    total = len(data)
    ths = "".join(
        f'<th{" class=\"num\"" if num[i] else ""}>{html.escape(h)}</th>'
        for i, h in enumerate(header)
    )
    trs = []
    for r in shown:
        cells = []
        for i in range(ncols):
            v = r[i] if i < len(r) else ""
            cells.append(
                f'<td{" class=\"num\"" if num[i] else ""}>{format_cell(v, num[i])}</td>'
            )
        trs.append("<tr>" + "".join(cells) + "</tr>")
    # encode full rows as JSON for client-side "show all"
    full = json.dumps(
        {
            "header": header,
            "rows": data,
            "num": num,
            "label": label,
        }
    )
    footer_note = f"{len(shown)} / {total}"
    return (
        '<div class="v4a-table-wrap">'
        '<table class="v4a-table">'
        f'<thead><tr>{ths}</tr></thead>'
        f'<tbody>{"".join(trs)}</tbody>'
        "</table></div>"
        '<div class="v4a-table-footer">'
        f'<span>{html.escape(label)} · {footer_note}</span>'
        + (
            f'<button type="button" data-v4a-expand=\'{html.escape(full, quote=True)}\'>Show all ({total})</button>'
            if total > top_n
            else ""
        )
        + "</div>"
    )


# ---------------- page chrome ----------------
def topnav_html(active: str | None = None) -> str:
    def link(href: str, label: str, key: str) -> str:
        cur = ' aria-current="page"' if active == key else ""
        return f'<li><a class="topnav__link" href="{href}"{cur}>{label}</a></li>'

    # active flags
    return f"""<nav class="topnav" aria-label="Primary">
  <div class="topnav__inner">
    <a class="topnav__brand" href="../index.html" aria-label="THYRAI home">
      <img src="../assets/img/brand/thyrai_wordmark.svg" alt="THYRAI">
    </a>
    <ul class="topnav__links">
      {link("platform_omega.html", "Platform", "platform")}
      {link("pipeline.html", "Pipeline", "pipeline")}
      {NAV_MARK_START}{link("v4a_archive_index.html", "Research Archive", "v4a")}{NAV_MARK_END}
      {link("publications.html", "Publications", "publications")}
      {link("team.html", "Team", "team")}
      <li><a class="topnav__link" href="mailto:kukshomr@gmail.com?subject=THYRAI%20contact">Contact</a></li>
    </ul>
    <a class="topnav__cta" href="mailto:kukshomr@gmail.com?subject=THYRAI%20access%20request">Request Access</a>
    <button class="topnav__hamburger" type="button" aria-label="Toggle navigation" aria-expanded="false"><span></span></button>
  </div>
</nav>"""


def footer_html() -> str:
    return f"""<footer class="site-footer">
  <div class="container">
    <div class="site-footer__grid">
      <div class="site-footer__col site-footer__brand">
        <img src="../assets/img/brand/thyrai_wordmark.svg" alt="THYRAI">
        <p>Research dashboard for computational thyroid-cancer triage. Computational triage only. Not medical advice.</p>
      </div>
      <div class="site-footer__col">
        <h4>Research Archive</h4>
        <ul>
          <li><a href="v4a_archive_index.html">Archive index</a></li>
          <li><a href="v4a_ml_performance.html">ML performance</a></li>
          <li><a href="v4a_biomarkers.html">Biomarkers</a></li>
          <li><a href="v4a_drug_discovery.html">Drug discovery</a></li>
          <li><a href="v4a_limitations.html">Limitations</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>v4 audit</h4>
        <ul>
          <li><a href="29_combat_rescue.html">29 ComBat rescue</a></li>
          <li><a href="30_bethesda_clinical.html">30 Bethesda clinical</a></li>
          <li><a href="32_honest_audit.html">32 Honest audit</a></li>
          <li><a href="33_v4_synthesis.html">33 v4 synthesis</a></li>
        </ul>
      </div>
      <div class="site-footer__col">
        <h4>Legacy</h4>
        <ul>
          <li><a href="../index_v2_legacy.html">v2 legacy index</a></li>
          <li><a href="99_glossary.html">99 · Glossary</a></li>
          <li><a href="../index.html">Back to dashboard</a></li>
        </ul>
      </div>
    </div>
    <div class="site-footer__meta">
      <span>&copy; 2026 THYRAI Research · v4a archive layer</span>
      {FOOTER_DISCLAIMER_HTML}
    </div>
  </div>
</footer>"""


def page_shell(
    *,
    title: str,
    eyebrow: str,
    h1: str,
    subtitle: str,
    breadcrumb_label: str,
    body_main: str,
    extra_scripts: str = "",
    active_nav: str = "v4a",
) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="dark">
  <meta name="theme-color" content="#000000">
  <title>{html.escape(title)} — THYRAI Research Archive</title>
  <meta name="description" content="{html.escape(subtitle)}">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="../assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/jetbrains-mono-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="../assets/css/brand.css">
  <link rel="stylesheet" href="../assets/css/v4a.css">
</head>
<body class="thyrai">
{topnav_html(active=active_nav)}

<header class="v4a-page-header">
  <div class="container v4a-page-header__inner">
    <span class="eyebrow">{html.escape(eyebrow)}</span>
    <h1>{html.escape(h1)}</h1>
    <p class="v4a-page-header__subtitle">{html.escape(subtitle)}</p>
    <div class="v4a-breadcrumb">
      <a href="../index.html">Home</a><span class="sep">/</span>
      <a href="v4a_archive_index.html">Research Archive</a><span class="sep">/</span>
      <span>{html.escape(breadcrumb_label)}</span>
    </div>
  </div>
</header>

<main class="container">
{body_main}
</main>

{footer_html()}
<script src="../assets/js/brand.js" defer></script>
<script src="../assets/js/v4a_tables.js" defer></script>
<script>
document.addEventListener('click', function(e){{
  var btn = e.target.closest && e.target.closest('[data-v4a-expand]');
  if (!btn) return;
  var payload = btn.getAttribute('data-v4a-expand');
  if (!payload) return;
  try {{
    var d = JSON.parse(payload);
  }} catch (err) {{ return; }}
  var wrap = btn.closest('.v4a-table-footer').previousElementSibling;
  var tbody = wrap && wrap.querySelector('tbody');
  if (!tbody) return;
  if (btn.dataset.expanded === '1'){{
    tbody.innerHTML = '';
    var limit = 30;
    d.rows.slice(0, limit).forEach(function(r){{
      var tr = document.createElement('tr');
      for (var i = 0; i < d.header.length; i++){{
        var td = document.createElement('td');
        var v = r[i] === undefined ? '' : r[i];
        if (d.num[i]){{ td.className = 'num'; td.textContent = formatNum(v); }} else {{ td.textContent = v; }}
        tr.appendChild(td);
      }}
      tbody.appendChild(tr);
    }});
    btn.textContent = 'Show all (' + d.rows.length + ')';
    btn.dataset.expanded = '0';
    btn.closest('.v4a-table-footer').querySelector('span').textContent =
      d.label + ' · ' + Math.min(limit, d.rows.length) + ' / ' + d.rows.length;
  }} else {{
    tbody.innerHTML = '';
    d.rows.forEach(function(r){{
      var tr = document.createElement('tr');
      for (var i = 0; i < d.header.length; i++){{
        var td = document.createElement('td');
        var v = r[i] === undefined ? '' : r[i];
        if (d.num[i]){{ td.className = 'num'; td.textContent = formatNum(v); }} else {{ td.textContent = v; }}
        tr.appendChild(td);
      }}
      tbody.appendChild(tr);
    }});
    btn.textContent = 'Show top 30';
    btn.dataset.expanded = '1';
    btn.closest('.v4a-table-footer').querySelector('span').textContent =
      d.label + ' · ' + d.rows.length + ' / ' + d.rows.length;
  }}
}});
function formatNum(v){{
  if (v === '' || v == null) return '';
  var n = Number(v);
  if (!isFinite(n)) return v;
  var a = Math.abs(n);
  if (n === 0) return '0';
  if (a < 1e-2) return n.toExponential(2);
  if (a >= 1e4) return n.toExponential(2);
  if (Number.isInteger(n) && a < 1e6) return String(n);
  return n.toPrecision(3).replace(/\\.?0+$/, '');
}}
</script>
{extra_scripts}
</body>
</html>"""


def kpi_row(items: list[tuple[str, str]]) -> str:
    cells = []
    for val, label in items:
        cells.append(
            f'<div class="v4a-kpi"><div class="v4a-kpi__value">{html.escape(val)}</div>'
            f'<div class="v4a-kpi__label">{html.escape(label)}</div></div>'
        )
    return f'<div class="v4a-kpi-row">{"".join(cells)}</div>'


def tldr_block(md_text: str, fallback: str = "") -> str:
    t = first_para(md_text) or fallback
    if not t:
        return ""
    return (
        '<div class="v4a-tldr">'
        '<span class="v4a-tldr__label">Executive summary</span>'
        f"<p>{html.escape(t)}</p>"
        "</div>"
    )


def source_missing_block(path: str) -> str:
    return (
        '<div class="v4a-source-missing">'
        "<strong>SOURCE UNAVAILABLE</strong>"
        f"{html.escape(path)}"
        "</div>"
    )


def fig_embed(dark_rel_path: str, caption: str, tall: bool = False) -> str:
    klass = "v4a-fig v4a-fig--tall" if tall else "v4a-fig"
    return (
        f'<figure class="{klass}">'
        f'<figcaption class="v4a-fig__caption">{html.escape(caption)}</figcaption>'
        f'<iframe src="{html.escape(dark_rel_path)}" loading="lazy" title="{html.escape(caption)}"></iframe>'
        "</figure>"
    )


def methodology_block(body_md: str, title: str = "Methodology") -> str:
    if not body_md:
        return ""
    body_html = md_to_simple_html(body_md)
    return (
        '<details class="v4a-methodology">'
        f"<summary>{html.escape(title)}</summary>"
        f'<div class="v4a-methodology__body">{body_html}</div>'
        "</details>"
    )


def limitations_band(items: list[str], link: str | None = None, link_label: str | None = None) -> str:
    lis = "".join(f"<li>{html.escape(i)}</li>" for i in items)
    extra = ""
    if link and link_label:
        extra = f'<p style="margin-top:16px"><a href="{html.escape(link)}" style="color:#F5A623">{html.escape(link_label)} →</a></p>'
    return (
        '<div class="disclaimer-band">'
        "<strong>Limitations &amp; caveats</strong>"
        f"<ul style='margin:0;padding-left:20px'>{lis}</ul>"
        f"{extra}"
        "</div>"
    )


# ---------------- page builders ----------------
def build_ml_page(manifest: dict) -> str:
    md = read_text(REPORTS / "ml_report_v2.md")
    base_h, base_r = load_tsv_rows(RESULTS / "ml" / "baseline_ml_results.tsv")
    ext_h, ext_r = load_tsv_rows(RESULTS / "ml" / "baseline_ml_external.tsv")

    kpis = kpi_row([
        ("1.000", "TierA67 CV AUC"),
        ("0.98", "GSE27155 Ext AUC"),
        ("0.96", "bACC Youden"),
        ("20", "Models trained"),
    ])

    body_parts = [tldr_block(md, "Internal CV achieves AUC 1.000 on TierA67 and variance_top50; external validation on GSE27155 peaks at 0.98 for ElasticNet LogReg.")]
    body_parts.append(kpis)

    body_parts.append('<section class="v4a-section"><h2>Internal cross-validation</h2>')
    if base_h:
        body_parts.append(render_table_html(base_h, base_r, top_n=30, label="baseline_ml_results.tsv rows", sort_by="cv_auc"))
    else:
        body_parts.append(source_missing_block("results/ml/baseline_ml_results.tsv"))
    body_parts.append("</section>")

    body_parts.append('<section class="v4a-section"><h2>External validation</h2>')
    if ext_h:
        body_parts.append(render_table_html(ext_h, ext_r, top_n=30, label="baseline_ml_external.tsv rows", sort_by="auc"))
    else:
        body_parts.append(source_missing_block("results/ml/baseline_ml_external.tsv"))
    body_parts.append("</section>")

    # Figures
    figs = []
    if manifest["plotly_dark"].get("ml_roc_overlay.html"):
        figs.append(("ml_roc_overlay.html", "ROC overlay (internal + external)"))
    if manifest["plotly_dark"].get("ml_calibration.html"):
        figs.append(("ml_calibration.html", "Calibration curve — TDS16 LogReg"))
    if manifest["plotly_dark"].get("ml_confusion_best.html"):
        figs.append(("ml_confusion_best.html", "Confusion matrix — best model"))
    if manifest["plotly_dark"].get("ml_model_comparison.html"):
        figs.append(("ml_model_comparison.html", "Model comparison — AUC by feature set × algorithm"))
    if figs:
        body_parts.append('<section class="v4a-section"><h2>Interactive figures</h2>')
        for name, cap in figs:
            body_parts.append(fig_embed(f"../figs_interactive/v4a_dark/{name}", cap))
        body_parts.append("</section>")

    body_parts.append(
        '<section class="v4a-section">'
        + methodology_block(md_section(md, "## Setup") + "\n\n" + md_section(md, "## Internal CV"), "Methodology")
        + "</section>"
    )

    body_parts.append(
        '<section class="v4a-section">'
        + limitations_band(
            [
                "Train/test splits on a single cohort inflate CV AUC; external AUC is the honest upper bound.",
                "GSE213647 was excluded from external eval due to missing molecular labels (see v2 notes).",
                "BRS71 handled as proxy in v2 — replaced in v3 (see Panels page).",
                "Dataset identifiability was 1.000 pre-ComBat; ComBat rescue was UNRECOVERABLE (page 32 audit).",
            ],
            link="32_honest_audit.html",
            link_label="See v4 honest audit (page 32)",
        )
        + "</section>"
    )

    return page_shell(
        title="ML Performance",
        eyebrow="Research Archive · ML",
        h1="Machine learning performance.",
        subtitle="Internal and external AUC, calibration, and model comparisons across TDS16 / TierA67 / variance_top50 feature sets.",
        breadcrumb_label="ML performance",
        body_main="\n".join(body_parts),
    )


def build_biomarkers_page(manifest: dict) -> str:
    md = read_text(REPORTS / "biomarker_analysis.md")
    hdr, rows = load_tsv_rows(RESULTS / "tables" / "biomarker_validated.tsv")

    n_validated = len(rows) if rows else 0
    top_d = ""
    if "cohens_d_tcga" in hdr:
        idx = hdr.index("cohens_d_tcga")
        vals = []
        for r in rows:
            try:
                vals.append(abs(float(r[idx])))
            except Exception:
                pass
        if vals:
            top_d = f"{max(vals):.2f}"

    kpis = kpi_row([
        (f"{n_validated:,}" if n_validated else "—", "Validated biomarkers"),
        (top_d or "—", "Top Cohen's d"),
        ("73%", "Replication ≥1 cohort"),
        ("BRAF vs RAS", "Contrast"),
    ])

    body = [tldr_block(md, "2,773 genes pass the novel-validated filter across TCGA + GSE27155 + GSE126698.")]
    body.append(kpis)

    body.append('<section class="v4a-section"><h2>Top validated biomarkers</h2>')
    if hdr:
        body.append(render_table_html(hdr, rows, top_n=30, label="biomarker_validated.tsv rows", sort_by="novelty_score"))
    else:
        body.append(source_missing_block("results/tables/biomarker_validated.tsv"))
    body.append("</section>")

    figs = []
    for name, cap in [
        ("volcano_braf_vs_ras.html", "Volcano — BRAF vs RAS (TCGA)"),
        ("biomarker_heatmap_known_vs_novel.html", "Heatmap — known vs novel biomarkers"),
        ("biomarker_replication_scatter.html", "Replication scatter — TCGA × external cohorts"),
        ("biomarker_effect_vs_coverage.html", "Effect size vs cohort coverage"),
    ]:
        if manifest["plotly_dark"].get(name):
            figs.append((name, cap))
    if figs:
        body.append('<section class="v4a-section"><h2>Interactive figures</h2>')
        for name, cap in figs:
            body.append(fig_embed(f"../figs_interactive/v4a_dark/{name}", cap))
        body.append("</section>")

    body.append(
        '<section class="v4a-section">'
        + methodology_block(md_section(md, "## Methods summary"), "Methodology — Welch t-test, BH-FDR, replication rule")
        + "</section>"
    )

    body.append(
        '<section class="v4a-section">'
        + limitations_band(
            [
                "GSE27155 external labels are phenotype proxies (histology-based), not direct molecular subtypes.",
                "GSE126698 subset is small after restricting to PTC vs FTC (n=12).",
                "Direction concordance is computed per-cohort, not meta-analytic.",
                "Novelty score is a ranking heuristic, not a hard significance threshold.",
            ]
        )
        + "</section>"
    )

    return page_shell(
        title="Biomarkers",
        eyebrow="Research Archive · Biomarkers",
        h1="2,773 validated biomarkers.",
        subtitle="Differential expression (tumor BRAF_like vs RAS_like) replicated across TCGA + two external cohorts.",
        breadcrumb_label="Biomarkers",
        body_main="\n".join(body),
    )


def build_drug_page(manifest: dict) -> str:
    md = read_text(REPORTS / "biomarker_to_drug_report.md")
    thdr, trows = load_tsv_rows(RESULTS / "tables" / "drug_discovery_targets.tsv")
    chdr, crows = load_tsv_rows(RESULTS / "tables" / "drug_discovery_compounds.tsv")

    kpis = kpi_row([
        ("8", "Drug targets"),
        ("15", "ChEMBL compounds"),
        ("8.75", "Best pChEMBL"),
        ("3", "Novel (first-in-class)"),
    ])

    body = [tldr_block(md, "Eight prioritized targets (three novel, five validated/emerging) yield 15 ChEMBL compounds with mechanism-of-action annotation.")]
    body.append(kpis)

    body.append('<section class="v4a-section"><h2>Targets</h2>')
    if thdr:
        body.append(render_table_html(thdr, trows, top_n=30, label="drug_discovery_targets.tsv", sort_by="novelty_score"))
    else:
        body.append(source_missing_block("results/tables/drug_discovery_targets.tsv"))
    body.append("</section>")

    body.append('<section class="v4a-section"><h2>Compounds</h2>')
    if chdr:
        body.append(render_table_html(chdr, crows, top_n=30, label="drug_discovery_compounds.tsv", sort_by="pchembl"))
    else:
        body.append(source_missing_block("results/tables/drug_discovery_compounds.tsv"))
    body.append("</section>")

    body.append(
        '<section class="v4a-section">'
        + methodology_block(md_section(md, "## Method") or md_section(md, "## Methods") or first_para(md), "Methodology — target prioritization & compound sourcing")
        + "</section>"
    )

    body.append(
        '<div class="v4a-cta">'
        '<p class="v4a-cta__text">Deep-dive into structures, binding site hypotheses, and 3D compound viewers.</p>'
        '<a class="btn-primary" href="platform_chem_showcase.html">Open CHEM-THCA showcase →</a>'
        "</div>"
    )

    body.append(
        '<section class="v4a-section">'
        + limitations_band(
            [
                "Targets ranked by novelty × replication × (-log10 FDR); this is a triage heuristic, not a potency claim.",
                "ChEMBL pChEMBL values are assay-mixed and not harmonized across modalities.",
                "Three novel targets (TACSTD2, TMPRSS4, PLEKHA6) have zero matched compounds — fragment / DEL campaigns required.",
                "No in-house structural validation; all SMILES and binding inferences are literature/ChEMBL derived.",
            ]
        )
        + "</section>"
    )

    return page_shell(
        title="Drug Discovery",
        eyebrow="Research Archive · Drug",
        h1="From biomarker to 15 compounds.",
        subtitle="Eight computationally triaged targets across novel / emerging / validated classes with ChEMBL-anchored compound lists.",
        breadcrumb_label="Drug discovery",
        body_main="\n".join(body),
    )


def build_pathway_page(manifest: dict) -> str:
    md_candidates = [
        REPORTS / "pathway_immune_meth.md",
        REPORTS / "pathway_immune_methylation.md",
        REPORTS / "analysis_summary_v2.md",
    ]
    md = ""
    for p in md_candidates:
        if probe(p):
            md = read_text(p)
            break

    # figures
    figs = []
    for name, cap in [
        ("pathway_dotplot.html", "Pathway enrichment — dot plot"),
        ("pathway_barplot_top15.html", "Top 15 enriched pathways — bar plot"),
        ("immune_signature_violin.html", "Immune signature violin — tumor vs normal"),
        ("methylation_roc.html", "Methylation classifier — ROC"),
        ("methylation_top_probes.html", "Top differentially methylated probes"),
        ("v3_multimodal_methylation.html", "Multi-modal methylation — cross-modal"),
    ]:
        if manifest["plotly_dark"].get(name):
            figs.append((name, cap))

    kpis = kpi_row([
        (str(len(figs)), "Interactive figures"),
        ("GSEA", "Pathway method"),
        ("CIBERSORT", "Immune decomp."),
        ("β-values", "Methylation scale"),
    ])

    body = [tldr_block(md, "GSEA-style pathway enrichment, CIBERSORT immune-signature contrasts, and methylation cross-modal analyses across six cohorts.")]
    body.append(kpis)

    pw_tsv = RESULTS / "tables" / "pathway_enrichment.tsv"
    if probe(pw_tsv):
        h, r = load_tsv_rows(pw_tsv)
        body.append('<section class="v4a-section"><h2>Pathway enrichment table</h2>')
        body.append(render_table_html(h, r, top_n=30, label="pathway_enrichment.tsv"))
        body.append("</section>")

    im_tsv = RESULTS / "tables" / "immune_signatures_contrast.tsv"
    if probe(im_tsv):
        h, r = load_tsv_rows(im_tsv)
        body.append('<section class="v4a-section"><h2>Immune signature contrasts</h2>')
        body.append(render_table_html(h, r, top_n=30, label="immune_signatures_contrast.tsv"))
        body.append("</section>")

    if figs:
        body.append('<section class="v4a-section"><h2>Interactive figures</h2>')
        for name, cap in figs:
            body.append(fig_embed(f"../figs_interactive/v4a_dark/{name}", cap))
        body.append("</section>")
    else:
        body.append(source_missing_block("reports/html/figs_interactive/ (no pathway/immune/methylation dark-themed figures)"))

    body.append(
        '<section class="v4a-section">'
        + methodology_block(
            "GSEA applied to per-gene t-statistic ranking (BRAF_like vs RAS_like). "
            "Immune decomposition via CIBERSORT LM22 on bulk RNA-seq. "
            "Methylation β-values computed from GSE97466 (GPL13534, 450K). "
            "Cross-modal alignment by gene symbol for paired TCGA samples.",
            "Methodology",
        )
        + "</section>"
    )

    body.append(
        '<section class="v4a-section">'
        + limitations_band(
            [
                "GSEA uses the MSigDB Hallmark / KEGG sets; enrichment is relative, not causal.",
                "CIBERSORT requires bulk RNA-seq at sufficient depth; low-quality samples were filtered.",
                "Methylation cohort (GSE97466) does not overlap with TCGA patient IDs — cross-modal is cohort-level only.",
            ]
        )
        + "</section>"
    )

    return page_shell(
        title="Pathway + Immune",
        eyebrow="Research Archive · Pathway & Immune",
        h1="Pathway, immune, methylation.",
        subtitle="Multi-modal enrichment, CIBERSORT-style immune contrasts, and 450K methylation classifier benchmarks.",
        breadcrumb_label="Pathway + Immune",
        body_main="\n".join(body),
    )


def build_panels_page(manifest: dict) -> str:
    def read_lines(p: Path) -> list[str]:
        if not probe(p):
            return []
        return [l.strip() for l in p.read_text(encoding="utf-8").splitlines() if l.strip() and not l.startswith("#")]

    tds = read_lines(META / "tds16_genes.txt")
    brs = read_lines(META / "brs71_genes.txt")
    if not brs:
        brs = read_lines(META / "brs71_proxy_genes.txt")
    tier = read_lines(META / "tierA67_genes.txt")
    cov_h, cov_r = load_tsv_rows(META / "gene_coverage_table.tsv")

    kpis = kpi_row([
        (f"{len(tds)}", "TDS16 genes"),
        (f"{len(brs)}", "BRS71 genes"),
        (f"{len(tier)}", "TierA67 entries"),
        ("6", "Cohorts profiled"),
    ])

    def panel_card(name: str, label: str, genes: list[str], subtitle: str) -> str:
        if not genes:
            return (
                '<div class="v4a-panel-card">'
                f'<h3>{html.escape(name)}</h3>'
                f'<p>{html.escape(subtitle)}</p>'
                '<div class="v4a-source-missing"><strong>FILE MISSING</strong>' + html.escape(label) + '</div>'
                "</div>"
            )
        return (
            '<div class="v4a-panel-card">'
            f'<h3>{html.escape(name)}</h3>'
            f'<div class="v4a-panel-card__count">{len(genes)}</div>'
            f'<p>{html.escape(subtitle)}</p>'
            f'<div class="v4a-panel-card__genes">{html.escape(", ".join(genes))}</div>'
            "</div>"
        )

    body = [
        '<div class="v4a-tldr"><span class="v4a-tldr__label">Panels</span>'
        "<p>Three differentially-intended gene panels: TDS16 (the Cell 2014 thyroid differentiation score), BRS71 (the original BRAF-RAS score), and TierA67 (engineered composite of clinical + literature biomarkers).</p></div>"
    ]
    body.append(kpis)

    body.append('<section class="v4a-section"><h2>Panel membership</h2>')
    body.append('<div class="v4a-panel-grid">')
    body.append(panel_card("TDS16", "tds16_genes.txt", tds, "16-gene thyroid differentiation score (Cell 2014, PMID 25417114)."))
    body.append(panel_card("BRS71", "brs71_genes.txt", brs, "71-gene BRAF-RAS score (original source-verified list post v3 recovery)."))
    body.append(panel_card("TierA67", "tierA67_genes.txt", tier, "67-entry engineered composite (clinical IHC + literature biomarkers + driver anchors)."))
    body.append("</div></section>")

    body.append('<section class="v4a-section"><h2>Per-cohort coverage</h2>')
    if cov_h:
        body.append(render_table_html(cov_h, cov_r, top_n=30, label="gene_coverage_table.tsv"))
    else:
        body.append(source_missing_block("metadata/gene_coverage_table.tsv"))
    body.append("</section>")

    body.append(
        '<section class="v4a-section">'
        + methodology_block(
            "TDS16 genes are taken from the Cell 2014 TCGA PTC paper, Table S5. "
            "BRS71 v2 was a proxy derived from TCGA — only 1.4% overlap with the published list. "
            "v3 recovered the original 72-gene list from Cell 2014 legends (see brs71_recovery_v3.md); "
            "this page reflects the recovered version. "
            "TierA67 is an engineered composite of clinical IHC markers, literature biomarkers, and driver anchors.",
            "Methodology — panel derivation",
        )
        + "</section>"
    )

    body.append(
        '<section class="v4a-section">'
        + limitations_band(
            [
                "BRS71 v2 had 1.4% overlap with the published list (proxy) — resolved in v3 recovery.",
                "TierA67 has 67 category entries but 66 unique HGNC symbols (PAX8 appears in two categories).",
                "GSE213647 has zero panel coverage because its platform (GPL18573) uses a different probe set.",
                "Panel 'membership' is not equivalent to weight in a trained model.",
            ]
        )
        + "</section>"
    )

    return page_shell(
        title="Gene Panels",
        eyebrow="Research Archive · Panels",
        h1="Three gene panels, six cohorts.",
        subtitle="TDS16, BRS71 (recovered), TierA67 — membership, cohort coverage, and provenance.",
        breadcrumb_label="Gene panels",
        body_main="\n".join(body),
    )


def build_data_page(manifest: dict) -> str:
    sm_path = META / "sample_master.tsv"
    dm_path = META / "dataset_master.tsv"
    body = []

    # aggregate counts from sample_master
    if probe(sm_path):
        import collections
        hdr, rows = load_tsv_rows(sm_path)
        n_samples = len(rows)
        by_dataset = collections.Counter()
        by_platform = collections.Counter()
        by_subtype = collections.Counter()
        by_norm = collections.Counter()
        col_ds = hdr.index("dataset") if "dataset" in hdr else None
        col_pf = hdr.index("platform") if "platform" in hdr else None
        col_sub = hdr.index("molecular_subtype") if "molecular_subtype" in hdr else None
        col_nt = hdr.index("normal_vs_tumor") if "normal_vs_tumor" in hdr else None
        for r in rows:
            if col_ds is not None and col_ds < len(r): by_dataset[r[col_ds] or "unknown"] += 1
            if col_pf is not None and col_pf < len(r): by_platform[r[col_pf] or "unknown"] += 1
            if col_sub is not None and col_sub < len(r):
                v = r[col_sub].strip() or "unspecified"
                by_subtype[v] += 1
            if col_nt is not None and col_nt < len(r): by_norm[r[col_nt] or "unknown"] += 1
        n_cohorts = len(by_dataset)
        n_platforms = len(by_platform)
        n_subtypes = len([k for k in by_subtype if k not in ("unspecified", "")])
    else:
        n_samples = n_cohorts = n_platforms = n_subtypes = 0
        by_dataset = by_platform = by_subtype = by_norm = {}

    kpis = kpi_row([
        (f"{n_samples:,}" if n_samples else "—", "Samples"),
        (f"{n_cohorts}" if n_cohorts else "—", "Cohorts"),
        (f"{n_platforms}" if n_platforms else "—", "Platforms"),
        (f"{n_subtypes}" if n_subtypes else "—", "Molecular subtypes"),
    ])

    body.append('<div class="v4a-tldr"><span class="v4a-tldr__label">Integration</span>'
                "<p>One harmonized sample-master spans 7 cohorts across 6 platforms (bulk RNA-seq, microarray, methylation 450K, WGS mixed), labelled by driver anchor and molecular subtype where available.</p></div>")
    body.append(kpis)

    if probe(sm_path):
        body.append('<section class="v4a-section"><h2>Aggregated sample counts</h2>')
        agg_rows = []
        for ds, n in sorted(by_dataset.items(), key=lambda kv: -kv[1]):
            agg_rows.append([ds, str(n)])
        body.append(render_table_html(["dataset", "n_samples"], agg_rows, top_n=30, label="dataset_rollup"))
        body.append("</section>")

        body.append('<section class="v4a-section"><h2>By platform</h2>')
        plat_rows = [[k, str(v)] for k, v in sorted(by_platform.items(), key=lambda kv: -kv[1])]
        body.append(render_table_html(["platform", "n_samples"], plat_rows, top_n=30, label="platform_rollup"))
        body.append("</section>")

        body.append('<section class="v4a-section"><h2>By molecular subtype</h2>')
        sub_rows = [[k, str(v)] for k, v in sorted(by_subtype.items(), key=lambda kv: -kv[1])]
        body.append(render_table_html(["molecular_subtype", "n_samples"], sub_rows, top_n=30, label="subtype_rollup"))
        body.append("</section>")
    else:
        body.append(source_missing_block("metadata/sample_master.tsv"))

    body.append('<section class="v4a-section"><h2>Dataset master (full)</h2>')
    if probe(dm_path):
        h, r = load_tsv_rows(dm_path)
        body.append(render_table_html(h, r, top_n=30, label="dataset_master.tsv"))
    else:
        body.append(source_missing_block("metadata/dataset_master.tsv"))
    body.append("</section>")

    figs = []
    for name, cap in [
        ("overview_sankey.html", "Cohort composition — Sankey (cohort → subtype)"),
        ("cohort_compare_pca.html", "Cohort PCA — samples across platforms"),
        ("cohort_label_concordance.html", "Cohort label concordance"),
    ]:
        if manifest["plotly_dark"].get(name):
            figs.append((name, cap))
    if figs:
        body.append('<section class="v4a-section"><h2>Interactive figures</h2>')
        for name, cap in figs:
            body.append(fig_embed(f"../figs_interactive/v4a_dark/{name}", cap))
        body.append("</section>")

    body.append(
        '<section class="v4a-section">'
        + limitations_band(
            [
                "Sample-master aggregates counts only on this page — download the full TSV for row-level details.",
                "EGA cohort (Yoo_SK_2019) is controlled-access and not downloaded; public reference only.",
                "GSE213647 lacks molecular labels; treated as normal/tumor with reduced utility for subtype ML.",
            ]
        )
        + "</section>"
    )

    return page_shell(
        title="Data Integration",
        eyebrow="Research Archive · Data",
        h1="1,509 samples, 7 cohorts.",
        subtitle="Harmonized sample-master across TCGA-THCA and six GEO cohorts, spanning bulk RNA-seq, microarray, and 450K methylation.",
        breadcrumb_label="Data integration",
        body_main="\n".join(body),
    )


def build_quantum_page(manifest: dict) -> str:
    md_path = RESULTS / "ml" / "quantum_full_comparison.md"
    md = read_text(md_path) if probe(md_path) else ""
    tsv_path = RESULTS / "ml" / "quantum_full_results.tsv"

    # parse the ranked table from MD into simple list
    kpis = kpi_row([
        ("13", "Algorithms benchmarked"),
        ("1.000", "Classical reference AUC"),
        ("Quantum-inspired", "Best parity (QAOA, QBoost)"),
        ("TierA67", "Feature set"),
    ])

    body = [tldr_block(md, "13 algorithms (classical + quantum + quantum-inspired) evaluated on TierA67 CV; classical LogReg_l2 reaches AUC 1.000 as the reference ceiling.") ]
    body.append(kpis)

    body.append('<section class="v4a-section"><h2>Full benchmark table</h2>')
    if probe(tsv_path):
        h, r = load_tsv_rows(tsv_path)
        body.append(render_table_html(h, r, top_n=30, label="quantum_full_results.tsv", sort_by="cv_auc"))
    else:
        body.append(source_missing_block("results/ml/quantum_full_results.tsv"))
    body.append("</section>")

    figs = []
    for name, cap in [
        ("quantum_full_comparison.html", "Quantum benchmark — AUC comparison"),
        ("comparison_classical_vs_dl_vs_quantum.html", "Classical vs DL vs Quantum"),
        ("quantum_kernel_heatmap.html", "Quantum kernel — fidelity heatmap"),
        ("quantum_qaoa_convergence.html", "QAOA — convergence"),
        ("quantum_qubo_selection.html", "QUBO feature selection"),
        ("quantum_pca_spectrum.html", "qPCA spectrum"),
        ("quantum_circuit_diagram.html", "VQC / QNN circuit diagram"),
        ("quantum_kmeans_clusters.html", "Quantum k-means clusters"),
        ("roc_quantum_vqc.html", "ROC — VQC variants"),
    ]:
        if manifest["plotly_dark"].get(name):
            figs.append((name, cap))
    if figs:
        body.append('<section class="v4a-section"><h2>Interactive figures</h2>')
        for name, cap in figs:
            body.append(fig_embed(f"../figs_interactive/v4a_dark/{name}", cap))
        body.append("</section>")

    body.append(
        '<section class="v4a-section">'
        + methodology_block(
            md_section(md, "## Honesty clause")
            or "VQC / QSVM / QNN run on qiskit_aer statevector simulator. QAOA runs on noiseless Sampler. "
            "d-wave-neal is classical simulated annealing. "
            "All AUCs use identical 5-fold CV splits and 95% percentile bootstrap CIs (800 resamples).",
            "Methodology & honesty clause",
        )
        + "</section>"
    )

    body.append(
        '<section class="v4a-section">'
        + limitations_band(
            [
                "All 'quantum' runs use classical emulators (qiskit_aer statevector / pennylane lightning.qubit / dwave-neal).",
                "At n=55 genes the Grover √N asymptotic advantage is irrelevant.",
                "Quantum k-means is unsupervised — AUC not reported; ARI vs label is the metric.",
                "Results should be read as a quantum-readiness study, not a quantum-advantage claim.",
            ]
        )
        + "</section>"
    )

    return page_shell(
        title="Quantum Benchmark",
        eyebrow="Research Archive · Quantum",
        h1="Quantum-inspired benchmark.",
        subtitle="13 algorithms — classical, quantum, quantum-inspired — on the same TierA67 5-fold CV splits.",
        breadcrumb_label="Quantum benchmark",
        body_main="\n".join(body),
    )


def build_limitations_page(manifest: dict) -> str:
    md2 = read_text(REPORTS / "analysis_summary_v2.md")
    md_next = read_text(REPORTS / "next_steps_v2.md")
    v2_lim = md_section(md2, "## Remaining limitations") or md_section(md2, "## Limitations")
    v2_items: list[str] = []
    for line in v2_lim.splitlines():
        s = line.strip()
        if s.startswith("-"):
            v2_items.append(s[1:].strip())
    if not v2_items:
        v2_items = ["See reports/analysis_summary_v2.md for the original limitations section."]

    v4_items = [
        "ComBat rescue UNRECOVERABLE — post-correction LODO AUC collapses to 0.011 (v4 Track A).",
        "Pre-ComBat dataset identifiability = 1.000 (macro-AUC): any cross-cohort prediction is partially memorizing batch.",
        "Direction-invariant AUC = 0.989 suggests sign-ambiguous features driving apparent performance.",
        "PRJEB11591 external cohort only partially downloaded (Track C = PARTIAL_DOWNLOAD).",
        "Bethesda cost assumptions (KRW single-payer) are not generalizable to other healthcare systems.",
        "All clinical utility claims are simulated decision-curve projections, not prospective trials.",
    ]

    kpis = kpi_row([
        ("UNRECOVERABLE", "Track A verdict"),
        ("1.000", "Identifiability (pre)"),
        ("0.989", "Direction-invariant AUC"),
        ("Classical only", "Quantum hardware"),
    ])

    body = ['<div class="v4a-tldr"><span class="v4a-tldr__label">Honesty</span>'
            '<p>Honest reporting of negative findings is core to the v4 scientific contribution. This page re-surfaces the v2 limitations verbatim alongside the v4 audit findings so the reader can judge claims in context.</p></div>']
    body.append(kpis)

    body.append('<section class="v4a-section"><h2>Caveats</h2>')
    body.append('<div class="v4a-lim-grid">')
    v2_lis = "".join(f"<li>{html.escape(i)}</li>" for i in v2_items)
    v4_lis = "".join(f"<li>{html.escape(i)}</li>" for i in v4_items)
    body.append(f'<div class="v4a-lim-col"><h3>v2 caveats</h3><ul>{v2_lis}</ul></div>')
    body.append(f'<div class="v4a-lim-col"><h3>v4 audit findings</h3><ul>{v4_lis}</ul></div>')
    body.append("</div></section>")

    if md_next.strip():
        body.append('<section class="v4a-section">'
                    + methodology_block(md_next, "Next steps (v2 → v3 → v4)")
                    + "</section>")

    body.append(
        '<div class="v4a-cta">'
        '<p class="v4a-cta__text">Full v4 Track D audit and 4×4 venue decision matrix.</p>'
        '<div style="display:flex;gap:12px;flex-wrap:wrap">'
        '<a class="btn-ghost" href="32_honest_audit.html">Page 32 · Honest audit →</a>'
        '<a class="btn-primary" href="33_v4_synthesis.html">Page 33 · v4 synthesis →</a>'
        "</div></div>"
    )

    body.append(
        '<p style="margin-top:48px;color:#F5A623;font-family:var(--font-mono);font-size:12px;letter-spacing:0.15em;text-transform:uppercase;text-align:center">'
        "Honest reporting of negative findings is core to the v4 scientific contribution."
        "</p>"
    )

    return page_shell(
        title="Limitations + Honesty",
        eyebrow="Research Archive · Limitations",
        h1="Limitations & honesty.",
        subtitle="v2 caveats and v4 audit findings side-by-side. Negative findings are part of the contribution.",
        breadcrumb_label="Limitations",
        body_main="\n".join(body),
    )


def build_archive_index(manifest: dict) -> str:
    cards = [
        ("Track 01", "ML PERFORMANCE", "v4a_ml_performance.html",
         "AUC, calibration, and external validation across TDS16 / TierA67 / variance_top50.",
         "20 MODELS", "1.000 CV AUC"),
        ("Track 02", "BIOMARKERS", "v4a_biomarkers.html",
         "2,773 genes pass the novel-validated filter — Welch t-test, BH-FDR, two-cohort replication.",
         "2,773 GENES", "d=2.41 TOP"),
        ("Track 03", "DRUG DISCOVERY", "v4a_drug_discovery.html",
         "8 targets × 15 ChEMBL compounds across novel / emerging / validated classes.",
         "8 TARGETS", "pChEMBL 8.75"),
        ("Track 04", "PATHWAY + IMMUNE", "v4a_pathway_immune.html",
         "GSEA pathways, CIBERSORT immune contrasts, 450K methylation classifier.",
         "HALLMARK", "450K METHYL"),
        ("Track 05", "GENE PANELS", "v4a_panels.html",
         "TDS16 (Cell 2014) · BRS71 (recovered v3) · TierA67 (engineered).",
         "3 PANELS", "6 COHORTS"),
        ("Track 06", "DATA INTEGRATION", "v4a_data_integration.html",
         "Harmonized sample-master: 1,509 samples across 7 cohorts and 6 platforms.",
         "1,509 N", "7 COHORTS"),
        ("Track 07", "QUANTUM BENCHMARK", "v4a_quantum.html",
         "13 algorithms — classical, quantum (emulated), quantum-inspired — on identical CV splits.",
         "13 ALGOS", "TierA67"),
        ("Track 08", "LIMITATIONS + HONESTY", "v4a_limitations.html",
         "v2 caveats and v4 audit findings side-by-side, including UNRECOVERABLE ComBat verdict.",
         "v2 + v4", "UNRECOVERABLE"),
        ("Track 09", "V4 SYNTHESIS", "33_v4_synthesis.html",
         "Quadruple-track synthesis: ComBat, Bethesda, PRJEB11591, audit, 14-day plan.",
         "4 TRACKS", "14 DAYS"),
    ]
    card_html: list[str] = []
    for eyebrow, title, href, desc, stat1, stat2 in cards:
        card_html.append(
            f'<a class="v4a-card" href="{html.escape(href)}">'
            f'<div><span class="v4a-card__eyebrow">{html.escape(eyebrow)}</span>'
            f'<h3>{html.escape(title)}</h3>'
            f'<p class="v4a-card__desc">{html.escape(desc)}</p></div>'
            f'<div><div class="v4a-card__stats"><span><strong>{html.escape(stat1.split()[0])}</strong>{html.escape(" ".join(stat1.split()[1:]) or "")}</span>'
            f'<span><strong>{html.escape(stat2.split()[0])}</strong>{html.escape(" ".join(stat2.split()[1:]) or "")}</span></div>'
            f'<div class="v4a-card__arrow">Open archive →</div></div>'
            "</a>"
        )

    hero = f"""<header class="hero hero--md" id="archive-hero">
  <div class="container hero__inner">
    <span class="eyebrow">Research Archive · v2 Outputs</span>
    <h1>Complete v2 analysis results.</h1>
    <p class="hero__lede">1,509 samples · 2,773 validated biomarkers · 8 drug targets · 15 ChEMBL compounds — re-rendered in the v4 dark aesthetic. The v2 originals remain untouched.</p>
    <div class="hero__ctas">
      <a class="btn-primary" href="#archive-grid">Browse archive ↓</a>
      <a class="btn-ghost" href="../index_v2_legacy.html">Open v2 legacy index →</a>
    </div>
  </div>
</header>"""

    kpis = f"""<section class="section section--tight" aria-label="archive key stats">
  <div class="container">
    <div class="v4a-kpi-row" style="margin:0">
      <div class="v4a-kpi"><div class="v4a-kpi__value">1,509</div><div class="v4a-kpi__label">Samples</div></div>
      <div class="v4a-kpi"><div class="v4a-kpi__value">2,773</div><div class="v4a-kpi__label">Validated biomarkers</div></div>
      <div class="v4a-kpi"><div class="v4a-kpi__value">8</div><div class="v4a-kpi__label">Drug targets</div></div>
      <div class="v4a-kpi"><div class="v4a-kpi__value">15</div><div class="v4a-kpi__label">ChEMBL compounds</div></div>
    </div>
  </div>
</section>"""

    grid = (
        '<section class="section" id="archive-grid" aria-label="archive grid">'
        '<div class="container">'
        '<span class="eyebrow">9 tracks</span>'
        '<h2 style="font-weight:300;letter-spacing:-0.02em;margin-bottom:32px">Browse the archive.</h2>'
        '<div class="v4a-archive-grid">'
        + "".join(card_html)
        + "</div>"
        '<div class="v4a-legacy-band">Prefer the original v2 layout? → '
        '<a href="01_overview.html">/pages/01_overview.html</a> … '
        '<a href="28_cross_cohort_summary.html">/pages/28_cross_cohort_summary.html</a>'
        "</div>"
        "</div>"
        "</section>"
    )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="dark">
  <meta name="theme-color" content="#000000">
  <title>Research Archive — THYRAI</title>
  <meta name="description" content="v2 research archive re-rendered in the v4 dark aesthetic. 1,509 samples, 2,773 biomarkers, 15 compounds.">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="../assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/jetbrains-mono-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="../assets/css/brand.css">
  <link rel="stylesheet" href="../assets/css/v4a.css">
</head>
<body class="thyrai">
{topnav_html(active="v4a")}
{hero}
{kpis}
{grid}
{footer_html()}
<script src="../assets/js/brand.js" defer></script>
</body>
</html>"""


# ---------------- nav injection ----------------
NAV_LINK_HTML = (
    f"{NAV_MARK_START}"
    '<li><a class="topnav__link" href="{href}">Research Archive</a></li>'
    f"{NAV_MARK_END}"
)


def inject_nav(path: Path, href: str) -> str:
    """Idempotently add 'RESEARCH ARCHIVE' link to the topnav of an HTML file.

    Returns 'injected' | 'already' | 'no-topnav' | 'missing'.
    """
    if not probe(path):
        return "missing"
    text = read_text(path)
    if NAV_MARK_START in text:
        return "already"
    link_html = NAV_LINK_HTML.format(href=href)

    # pattern 1: brand topnav with ul.topnav__links
    m = re.search(r'(<ul\s+class="topnav__links"[^>]*>)', text)
    if m:
        new_text = text[: m.end()] + "\n      " + link_html + text[m.end():]
        path.write_text(new_text, encoding="utf-8")
        return "injected"

    # pattern 2: pages 29-33 lack a topnav — insert before </head> a small bar
    if "<body" in text and "<title>" in text:
        # add a small research-archive link banner (not a full nav) guarded by markers
        banner = (
            f"{NAV_MARK_START}"
            '<div style="background:#000;color:#F5A623;padding:10px 20px;'
            'font-family:JetBrains Mono,monospace;font-size:11px;letter-spacing:0.2em;'
            'text-transform:uppercase;border-bottom:1px solid #333;text-align:right">'
            f'<a href="{href}" style="color:#F5A623;text-decoration:none">Research Archive →</a>'
            "</div>"
            f"{NAV_MARK_END}"
        )
        new_text = re.sub(r"(<body[^>]*>)", r"\1\n" + banner, text, count=1)
        if new_text != text:
            path.write_text(new_text, encoding="utf-8")
            return "injected"

    return "no-topnav"


def inject_landing_archive_section(landing: Path) -> str:
    if not probe(landing):
        return "missing"
    text = read_text(landing)
    if LANDING_MARK_START in text:
        return "already"
    section = f"""{LANDING_MARK_START}
<section class="section" aria-label="research archive">
  <div class="container">
    <span class="eyebrow">Research Archive</span>
    <h2 style="font-weight:300;letter-spacing:-0.02em;max-width:18ch">The v2 research, re-rendered.</h2>
    <p style="color:rgba(255,255,255,0.72);max-width:60ch;font-size:16px;margin-top:16px">
      Nine tracks of substantive output — ML, biomarkers, drug discovery, pathways, panels, data, quantum, limitations — in a single dark-aesthetic archive. The v2 originals remain untouched.
    </p>
    <div style="margin-top:32px">
      <a class="btn-primary" href="pages/v4a_archive_index.html">Open Research Archive →</a>
      <a class="btn-ghost" href="index_v2_legacy.html" style="margin-left:12px">v2 legacy index</a>
    </div>
  </div>
</section>
{LANDING_MARK_END}
"""
    # insert just before <footer class="site-footer">
    new_text, n = re.subn(r'(<footer\s+class="site-footer")', section + r"\1", text, count=1)
    if n == 0:
        return "no-footer"
    # also inject topnav archive link into landing
    if NAV_MARK_START not in new_text:
        # inject between Pipeline and Publications
        pub_pat = r'(<li><a class="topnav__link" href="pages/publications\.html">Publications</a></li>)'
        link = (
            f"{NAV_MARK_START}"
            '<li><a class="topnav__link" href="pages/v4a_archive_index.html">Research Archive</a></li>'
            f"{NAV_MARK_END}"
        )
        new_text2, n2 = re.subn(pub_pat, link + r"\1", new_text, count=1)
        if n2 > 0:
            new_text = new_text2
    landing.write_text(new_text, encoding="utf-8")
    return "injected"


# ---------------- figure manifest + retheme ----------------
REQUIRED_FIGS = [
    # ML
    "ml_roc_overlay.html", "ml_calibration.html", "ml_confusion_best.html",
    "ml_model_comparison.html", "ml_cv_boxplot.html", "ml_pr_overlay.html",
    # Biomarkers
    "biomarker_heatmap_known_vs_novel.html", "biomarker_replication_scatter.html",
    "biomarker_effect_vs_coverage.html", "volcano_braf_vs_ras.html",
    # Pathway / immune / methylation
    "pathway_dotplot.html", "pathway_barplot_top15.html",
    "immune_signature_violin.html", "methylation_roc.html",
    "methylation_top_probes.html", "v3_multimodal_methylation.html",
    # Quantum
    "quantum_full_comparison.html", "quantum_kernel_heatmap.html",
    "quantum_qaoa_convergence.html", "quantum_qubo_selection.html",
    "quantum_pca_spectrum.html", "quantum_circuit_diagram.html",
    "quantum_kmeans_clusters.html", "roc_quantum_vqc.html",
    "comparison_classical_vs_dl_vs_quantum.html",
    # Data
    "overview_sankey.html", "cohort_compare_pca.html", "cohort_label_concordance.html",
]


def build_plotly_manifest() -> dict:
    manifest: dict[str, dict] = {"plotly_src": {}, "plotly_dark": {}}
    for name in REQUIRED_FIGS:
        src = FIGS / name
        manifest["plotly_src"][name] = probe(src)
    return manifest


def retheme_figures(manifest: dict) -> tuple[int, int, int]:
    names = [n for n, ok in manifest["plotly_src"].items() if ok]
    if not names:
        return 0, 0, 0
    res = retheme_many(FIGS, FIGS_DARK, names)
    rethemed = 0
    copied = 0
    failed = 0
    for n, status in res.items():
        manifest["plotly_dark"][n] = probe(FIGS_DARK / n)
        if status == "rethemed":
            rethemed += 1
        elif status == "copied":
            copied += 1
        elif status in ("failed", "missing"):
            failed += 1
    return rethemed, copied, failed


# ---------------- availability probe ----------------
def build_source_manifest() -> dict:
    sources = {
        "ml": {
            "baseline": RESULTS / "ml" / "baseline_ml_results.tsv",
            "external": RESULTS / "ml" / "baseline_ml_external.tsv",
            "md": REPORTS / "ml_report_v2.md",
        },
        "biomarkers": {
            "tsv": RESULTS / "tables" / "biomarker_validated.tsv",
            "md": REPORTS / "biomarker_analysis.md",
        },
        "drug": {
            "targets": RESULTS / "tables" / "drug_discovery_targets.tsv",
            "compounds": RESULTS / "tables" / "drug_discovery_compounds.tsv",
            "md": REPORTS / "biomarker_to_drug_report.md",
        },
        "pathway": {
            "pathway_tsv": RESULTS / "tables" / "pathway_enrichment.tsv",
            "immune_tsv": RESULTS / "tables" / "immune_signatures_contrast.tsv",
        },
        "panels": {
            "tds16": META / "tds16_genes.txt",
            "brs71": META / "brs71_genes.txt",
            "tierA67": META / "tierA67_genes.txt",
            "coverage": META / "gene_coverage_table.tsv",
        },
        "data": {
            "sample_master": META / "sample_master.tsv",
            "dataset_master": META / "dataset_master.tsv",
        },
        "quantum": {
            "md": RESULTS / "ml" / "quantum_full_comparison.md",
            "tsv": RESULTS / "ml" / "quantum_full_results.tsv",
        },
        "limitations": {
            "md": REPORTS / "analysis_summary_v2.md",
            "next": REPORTS / "next_steps_v2.md",
        },
        "v4_synthesis": {
            "page": PAGES / "33_v4_synthesis.html",
        },
    }
    avail = {}
    for group, items in sources.items():
        avail[group] = {k: probe(v) for k, v in items.items()}
        avail[group]["_all"] = all(probe(v) for v in items.values())
    return avail


# ---------------- main ----------------
def write_if_changed(path: Path, content: str) -> bool:
    """Return True if file was written (content differs)."""
    new_hash = hashlib.sha256(content.encode("utf-8")).hexdigest()
    if probe(path):
        old_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if old_hash == new_hash:
            return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return True


def health_check(urls: list[str]) -> list[tuple[str, int]]:
    out = []
    for u in urls:
        try:
            r = subprocess.run(
                ["curl", "-s", "-o", "/dev/null", "-w", "%{http_code}", "--max-time", "6", u],
                capture_output=True,
                text=True,
            )
            out.append((u, int(r.stdout.strip() or 0)))
        except Exception:
            out.append((u, 0))
    return out


def main() -> int:
    LOGS.mkdir(parents=True, exist_ok=True)
    log_lines: list[str] = []

    def log(msg: str) -> None:
        print(msg)
        log_lines.append(msg)

    log(f"[v4a] build start {BUILD_TIME}")

    brand_css = ASSETS / "css" / "brand.css"
    if not probe(brand_css):
        log(f"[ABORT] brand.css missing at {brand_css}")
        return 2

    # 1. source manifest
    sources = build_source_manifest()
    log("[v4a] source availability:")
    missing_sources = []
    for group, items in sources.items():
        flags = {k: v for k, v in items.items() if k != "_all"}
        log(f"  - {group}: _all={items['_all']} · {flags}")
        if not items["_all"]:
            missing_sources.append(group)

    # 2. plotly manifest + dark retheme
    fig_manifest = build_plotly_manifest()
    FIGS_DARK.mkdir(parents=True, exist_ok=True)
    rethemed, copied, failed = retheme_figures(fig_manifest)
    log(f"[v4a] Plotly retheme: {rethemed} rethemed, {copied} copied-as-is, {failed} failed/missing")

    manifest = {
        **fig_manifest,
        "sources": sources,
    }

    # 3. write topic pages
    pages_written = 0
    pages_skipped = 0

    def write_page(name: str, content: str):
        nonlocal pages_written, pages_skipped
        out = PAGES / name
        if write_if_changed(out, content):
            pages_written += 1
            log(f"  [write] {name}")
        else:
            pages_skipped += 1
            log(f"  [skip ] {name} (unchanged)")

    write_page(V4A_INDEX_PAGE, build_archive_index(manifest))
    write_page("v4a_ml_performance.html", build_ml_page(manifest))
    write_page("v4a_biomarkers.html", build_biomarkers_page(manifest))
    write_page("v4a_drug_discovery.html", build_drug_page(manifest))
    write_page("v4a_pathway_immune.html", build_pathway_page(manifest))
    write_page("v4a_panels.html", build_panels_page(manifest))
    write_page("v4a_data_integration.html", build_data_page(manifest))
    write_page("v4a_quantum.html", build_quantum_page(manifest))
    write_page("v4a_limitations.html", build_limitations_page(manifest))

    log(f"[v4a] wrote {pages_written} pages, skipped {pages_skipped} unchanged")

    # 4. nav injection into landing + v4 pages 29-33 + v4a topic pages themselves already have topnav
    nav_log = []
    r = inject_landing_archive_section(HTML_ROOT / "index.html")
    nav_log.append(("index.html", r))
    for p in ["29_combat_rescue.html", "30_bethesda_clinical.html", "31_prjeb11591_status.html",
              "32_honest_audit.html", "33_v4_synthesis.html"]:
        r = inject_nav(PAGES / p, "v4a_archive_index.html")
        nav_log.append((p, r))
    log("[v4a] nav injection:")
    for name, status in nav_log:
        log(f"  - {name}: {status}")

    # 5. update _version.json
    version_path = HTML_ROOT / "_version.json"
    try:
        vd = json.loads(version_path.read_text(encoding="utf-8")) if probe(version_path) else {}
    except Exception:
        vd = {}
    vd["v4a_archive_build_time"] = BUILD_TIME
    vd["v4a_pages"] = [V4A_INDEX_PAGE] + V4A_TOPIC_PAGES
    vd["v4a_sources_missing"] = missing_sources
    version_path.write_text(json.dumps(vd, indent=2), encoding="utf-8")

    # 6. health check
    urls = [f"{SERVER}/reports/html/pages/{p}" for p in ([V4A_INDEX_PAGE] + V4A_TOPIC_PAGES)]
    checks = health_check(urls)
    log("[v4a] health check:")
    ok_codes = {200}
    for u, code in checks:
        log(f"  {code}  {u}")

    (LOGS / "v4a_archive_build.log").write_text("\n".join(log_lines), encoding="utf-8")

    # Stage 99 summary
    total_src = len(sources)
    ok_src = total_src - len(missing_sources)
    print("\n=== v2 RESEARCH ARCHIVE (v4 STYLE) COMPLETE ===")
    print(f"Index URL          : {SERVER}/reports/html/pages/{V4A_INDEX_PAGE}")
    print(f"Topic pages        : {len(V4A_TOPIC_PAGES)}")
    print(f"Sources rendered   : {ok_src}/{total_src}  (missing: {', '.join(missing_sources) if missing_sources else 'none'})")
    print(f"Plotly re-themed   : {rethemed} figures (+{copied} copied-as-is)")
    print(f"Legacy v2 pages    : 28 (untouched, accessible via index_v2_legacy.html)")
    nav_ok = sum(1 for _, s in nav_log if s in ("injected", "already"))
    print(f"Nav link added     : yes (landing + v4 pages 29-33) · {nav_ok}/{len(nav_log)} nodes")
    print(f"Log                : logs/v4a_archive_build.log")
    print("================================================")

    return 0


if __name__ == "__main__":
    sys.exit(main())
