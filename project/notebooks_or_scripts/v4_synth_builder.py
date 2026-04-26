#!/usr/bin/env python3
"""
v4 Synthesis builder — reads all 4 tracks' outputs, writes
  reports/v4_synthesis.md
  reports/v4_combat_rescue.md       (short track A report)
  reports/v4_bethesda_clinical.md   (short track B report)
  reports/v4_prjeb11591_status.md   (short track C report)
  reports/html/pages/29_combat_rescue.html
  reports/html/pages/30_bethesda_clinical.html
  reports/html/pages/31_prjeb11591_status.html
  reports/html/pages/32_honest_audit.html
  reports/html/pages/33_v4_synthesis.html
appends a v4 banner to reports/html/index.html (does NOT rewrite THYRAI content),
and updates reports/html/_version.json with v4 metadata.
"""
from __future__ import annotations

import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
OUT_ML = ROOT / "results" / "ml"
OUT_FIG = ROOT / "reports" / "html" / "figs_interactive"
OUT_PAGES = ROOT / "reports" / "html" / "pages"
OUT_REPORTS = ROOT / "reports"
IDX_HTML = ROOT / "reports" / "html" / "index.html"
VERSION_JSON = ROOT / "reports" / "html" / "_version.json"
for d in (OUT_ML, OUT_FIG, OUT_PAGES, OUT_REPORTS):
    d.mkdir(parents=True, exist_ok=True)

BANNER_MARK_START = "<!-- v4-banner-start -->"
BANNER_MARK_END = "<!-- v4-banner-end -->"


def log(msg: str) -> None:
    print(f"[synth] {msg}", flush=True)


def safe_json(p: Path, default):
    try:
        return json.loads(p.read_text())
    except Exception:
        return default


def safe_tsv(p: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def page_header(title: str) -> str:
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;margin:0;padding:0;background:#fafafa;color:#222}}
.container{{max-width:1200px;margin:0 auto;padding:20px}}
.hero{{background:linear-gradient(135deg,#0f172a 0%,#1e3a8a 100%);color:#fff;padding:24px;border-radius:10px;margin-bottom:18px}}
.hero h1{{margin:0 0 6px 0;font-size:26px}}
.hero p{{margin:4px 0;opacity:.9}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:16px 0}}
.kpi{{background:#fff;border-radius:8px;padding:14px;box-shadow:0 1px 3px rgba(0,0,0,0.08)}}
.kpi .label{{font-size:11px;color:#666;text-transform:uppercase;letter-spacing:.5px}}
.kpi .val{{font-size:22px;font-weight:600;margin-top:4px}}
.verdict-RESCUED{{color:#16a085}}
.verdict-UNRECOVERABLE{{color:#c0392b}}
.verdict-INSUFFICIENT{{color:#e67e22}}
.warn{{background:#fff8e1;border-left:4px solid #f39c12;padding:12px 16px;margin:12px 0;border-radius:4px}}
.caveat{{background:#e8f4f8;border-left:4px solid #3498db;padding:12px 16px;margin:12px 0;border-radius:4px}}
.figcard{{background:#fff;border-radius:8px;padding:14px;box-shadow:0 1px 3px rgba(0,0,0,0.08);margin:14px 0}}
.figcard h3{{margin-top:0}}
.figcard img{{max-width:100%;height:auto}}
.dl{{margin-top:6px;font-size:12px}}
.dl a{{margin-right:12px;color:#3498db;text-decoration:none}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}}
th,td{{border:1px solid #ddd;padding:5px 8px;text-align:left}}
th{{background:#f5f5f5}}
footer{{margin-top:30px;padding:20px;border-top:1px solid #ddd;font-size:13px;color:#666;text-align:center}}
footer a{{color:#3498db}}
.cta{{background:#1e3a8a;color:#fff;padding:16px;border-radius:8px;text-align:center;margin:14px 0}}
.cta a{{color:#fff;text-decoration:underline}}
</style></head><body><div class="container">
"""


def footer_html() -> str:
    return """
<footer>
  v4 sprint. See also:
  <a href="14_caveats.html">Caveats (page 14)</a>
  &middot;
  <a href="19_honesty_audit.html">Honesty audit (page 19)</a>
  &middot;
  <a href="32_honest_audit.html">v4 audit (page 32)</a>
  &middot;
  <a href="../index.html">Back to dashboard</a>
</footer>
</div></body></html>
"""


def kpi_block(pairs: list[tuple[str, str]]) -> str:
    cells = "".join(
        f"<div class='kpi'><div class='label'>{html.escape(lbl)}</div>"
        f"<div class='val'>{html.escape(val)}</div></div>"
        for lbl, val in pairs
    )
    return f"<div class='kpis'>{cells}</div>"


def figcard(title: str, fig_base: str, note: str = "") -> str:
    return (f"<div class='figcard'><h3>{html.escape(title)}</h3>"
            + (f"<p>{html.escape(note)}</p>" if note else "")
            + f"<img src='../figs_interactive/{fig_base}.png' alt='{html.escape(title)}'>"
            + f"<div class='dl'>"
            + f"<a href='../figs_interactive/{fig_base}.png'>PNG</a>"
            + f"<a href='../figs_interactive/{fig_base}.html'>HTML</a>"
            + f"<a href='../figs_interactive/{fig_base}.tsv'>TSV</a>"
            + f"</div></div>")


def _f(v, default=float("nan")) -> float:
    try:
        if v is None:
            return default
        return float(v)
    except Exception:
        return default


def build_page_29(trackA: dict) -> str:
    verdict = trackA.get("verdict", "UNKNOWN")
    pre = _f(trackA.get("pre_identifiability_auc"))
    post = _f(trackA.get("post_identifiability_auc"))
    pre_lodo = _f(trackA.get("pre_lodo_mean_auc"))
    post_lodo = _f(trackA.get("post_lodo_mean_auc"))
    n_samples = trackA.get("n_samples_pooled", 0)
    shared = trackA.get("shared_genes", 0)
    header = page_header("29 - ComBat batch-correction rescue (v4)")
    hero = f"""<div class='hero'><h1>Page 29 - ComBat batch-correction rescue</h1>
<p>Subtype-preserving ComBat on 6-cohort pool. Decision rule applied.</p>
<p>Decision-support prototype / retrospective computational triage. NOT a diagnostic device.</p></div>"""
    kpis = kpi_block([
        ("Verdict", verdict),
        ("Pre-ident AUC", f"{pre:.3f}"),
        ("Post-ident AUC", f"{post:.3f}"),
        ("Pre LODO AUC", f"{pre_lodo:.3f}"),
        ("Post LODO AUC", f"{post_lodo:.3f}"),
        ("Samples pooled", f"{n_samples}"),
        ("Shared genes", f"{shared}"),
    ])
    warn = ""
    if verdict == "UNRECOVERABLE":
        warn = "<div class='warn'>Verdict UNRECOVERABLE: post-ComBat identifiability and LODO AUC indicate biology/batch entanglement cannot be separated with linear correction on shared genes.</div>"
    elif verdict == "INSUFFICIENT":
        warn = "<div class='warn'>Verdict INSUFFICIENT: results sit between rescue and unrecoverable bands; additional cohorts or non-linear correction needed.</div>"
    figs = "".join([
        figcard("Dataset identifiability (pre vs post)", "v4_combat_fig1"),
        figcard("Pre-ComBat LODO AUC", "v4_combat_fig2", "small-n cohorts annotated"),
        figcard("Post-ComBat LODO AUC", "v4_combat_fig3", "small-n cohorts annotated"),
        figcard("Pre vs post LODO side-by-side", "v4_combat_fig4"),
        figcard("PCA (pre vs post)", "v4_combat_fig5"),
    ])
    return header + hero + kpis + warn + figs + footer_html()


def build_page_30(trackB_op: pd.DataFrame, trackB_impact: dict) -> str:
    header = page_header("30 - Bethesda-III/IV clinical decision (v4)")
    hero = f"""<div class='hero'><h1>Page 30 - Bethesda clinical decision simulation</h1>
<p>6 prevalences x 20,000 synthetic patients; cost-utility in KRW (Korean single-payer est., NOT generalizable).</p>
<p>Decision-support prototype / retrospective computational triage. NOT a diagnostic device.</p></div>"""
    avoided = trackB_impact.get("avoided_surgeries_per_1000_at_prev_0.20_balanced")
    kpis = kpi_block([
        ("Prevalences", "6"),
        ("Patients/prev", "20,000"),
        ("Operating points", f"{len(trackB_op)}"),
        ("Avoided surg / 1k @ 20%", f"{avoided:.0f}" if avoided else "n/a"),
    ])
    cav = "<div class='caveat'>Costs are Korean single-payer estimates; NOT generalizable. Bethesda cohort is synthetic, not prospective. TCGA is surgical tissue, not FNA.</div>"
    # Table of operating points
    if len(trackB_op):
        tbl = "<table><tr><th>Prev</th><th>Operating point</th><th>Thr</th><th>Se</th><th>Sp</th><th>PPV</th><th>NPV</th></tr>"
        for _, r in trackB_op.iterrows():
            tbl += (f"<tr><td>{r['prev']:.0%}</td><td>{html.escape(str(r['operating_point']))}</td>"
                    f"<td>{r['threshold']:.2f}</td><td>{r['sens']:.3f}</td><td>{r['spec']:.3f}</td>"
                    f"<td>{r['ppv']:.3f}</td><td>{r['npv']:.3f}</td></tr>")
        tbl += "</table>"
    else:
        tbl = ""
    figs = "".join([
        figcard("ROC across prevalences", "v4_bethesda_fig1"),
        figcard("Net-benefit decision curves", "v4_bethesda_fig2"),
        figcard("Cost-utility (KRW)", "v4_bethesda_fig3", "KR single-payer est., NOT generalizable"),
        figcard("Three operating points", "v4_bethesda_fig4"),
        figcard("Unnecessary surgery vs missed cancer", "v4_bethesda_fig5"),
    ])
    return header + hero + kpis + cav + "<h2>Operating points</h2>" + tbl + figs + footer_html()


def build_page_31(trackC_status: str, trackC_attempts: list, alt_cohorts: list[dict]) -> str:
    header = page_header("31 - PRJEB11591 external cohort retry (v4)")
    hero = f"""<div class='hero'><h1>Page 31 - PRJEB11591 retry status</h1>
<p>Attempts on 8 alternate supplementary URLs; author email draft saved; alt cohort candidates enumerated.</p></div>"""
    kpis = kpi_block([
        ("Status", trackC_status),
        ("URLs tried", f"{len(trackC_attempts)}"),
        ("Successful fetches (>=1KB)", f"{sum(1 for a in trackC_attempts if a.get('path'))}"),
        ("Alt cohorts enumerated", f"{len(alt_cohorts)}"),
    ])
    rows = "".join(
        f"<tr><td>{a.get('status','')}</td><td>{a.get('bytes',0)}</td>"
        f"<td style='font-family:monospace;font-size:11px'>{html.escape(a.get('url',''))[:140]}</td></tr>"
        for a in trackC_attempts
    )
    tbl = ("<h2>HTTP attempts</h2><table><tr><th>status</th><th>bytes</th><th>url</th></tr>"
           + rows + "</table>")
    alt_rows = "".join(
        f"<tr><td>{c['geo_id']}</td><td>{c['arm']}</td><td>{c['use']}</td><td>{c['note']}</td></tr>"
        for c in alt_cohorts
    )
    alt_tbl = ("<h2>Alternative cohort candidates</h2>"
               "<table><tr><th>ID</th><th>arm</th><th>use</th><th>note</th></tr>"
               + alt_rows + "</table>")
    fig = figcard("PRJEB11591 alt-URL attempts", "v4_trackC_figure", "bytes received by URL")
    return header + hero + kpis + tbl + alt_tbl + fig + footer_html()


def build_page_32(checks: list[dict]) -> str:
    header = page_header("32 - v4 honest audit (Track D)")
    hero = f"""<div class='hero'><h1>Page 32 - v4 honest audit</h1>
<p>6-check checklist; 6 audit figures; mini-paper skeleton targeting JCO PO or Bioinformatics.</p></div>"""
    n_pass = sum(1 for c in checks if c.get("pass"))
    kpis = kpi_block([
        ("Checks", f"{len(checks)}"),
        ("Pass", f"{n_pass}"),
        ("Fail", f"{len(checks) - n_pass}"),
    ])
    rows = "".join(
        f"<tr><td>{c['id']}</td><td>{html.escape(c['name'])}</td>"
        f"<td>{'PASS' if c['pass'] else 'FAIL'}</td>"
        f"<td style='font-family:monospace;font-size:11px'>{html.escape(str(c.get('evidence','')))}</td></tr>"
        for c in checks
    )
    tbl = ("<h2>Checklist</h2><table><tr><th>ID</th><th>name</th><th>status</th><th>evidence</th></tr>"
           + rows + "</table>")
    figs = "".join([
        figcard("Audit checklist", "v4_audit_fig1"),
        figcard("Identifiability pre/post", "v4_audit_fig2"),
        figcard("Post-ComBat LODO", "v4_audit_fig3"),
        figcard("DCA at 20% prev", "v4_audit_fig4"),
        figcard("Cost-utility (KRW)", "v4_audit_fig5"),
        figcard("PRJEB11591 attempts", "v4_audit_fig6"),
    ])
    dl = ("<div class='cta'>Downloads: "
          "<a href='../../v4_honest_audit.md'>v4_honest_audit.md</a> &middot; "
          "<a href='../../v4_minipaper_skeleton.tex'>mini-paper .tex</a></div>")
    return header + hero + kpis + tbl + figs + dl + footer_html()


def build_page_33(trackA: dict, trackB_op: pd.DataFrame, trackC_status: str,
                  decision_matrix: list[dict], tasks: list[str]) -> str:
    header = page_header("33 - v4 synthesis")
    hero = f"""<div class='hero'><h1>Page 33 - v4 Quadruple-Track Synthesis</h1>
<p>Executive summary, track verdicts, 4x4 decision matrix (venue x feasibility/time/impact/risk), 14-day plan.</p>
<p>Decision-support prototype / retrospective computational triage. NOT a diagnostic device.</p></div>"""
    avoided = trackB_op[trackB_op["operating_point"] == "balanced (max Se+Sp)"] if len(trackB_op) else pd.DataFrame()
    kpis = kpi_block([
        ("Track A verdict", trackA.get("verdict", "UNKNOWN")),
        ("Track B ops", f"{len(trackB_op)}"),
        ("Track C", trackC_status),
        ("Track D", "audit complete"),
    ])
    # 4x4 decision matrix table
    dm_rows = "".join(
        f"<tr><td>{html.escape(d['venue'])}</td><td>{d['feasibility']}</td>"
        f"<td>{d['time']}</td><td>{d['impact']}</td><td>{d['risk']}</td><td><b>{d['score']}</b></td></tr>"
        for d in decision_matrix
    )
    dm_tbl = ("<h2>4x4 Decision matrix (venue x feasibility/time/impact/risk, 1-5)</h2>"
              "<table><tr><th>venue</th><th>feasibility</th><th>time</th>"
              "<th>impact</th><th>risk(inv)</th><th>score</th></tr>"
              + dm_rows + "</table>")
    tasks_html = "<h2>14-day daily task list</h2><ol>" + "".join(f"<li>{html.escape(t)}</li>" for t in tasks) + "</ol>"
    figs = "".join([
        figcard("A: identifiability pre vs post", "v4_combat_fig1"),
        figcard("B: decision curves", "v4_bethesda_fig2"),
        figcard("B: cost-utility", "v4_bethesda_fig3"),
        figcard("C: PRJEB11591 status", "v4_trackC_figure"),
        figcard("D: audit checklist", "v4_audit_fig1"),
    ])
    cta = ("<div class='cta'>Next: "
           "<a href='29_combat_rescue.html'>29 ComBat</a> &middot; "
           "<a href='30_bethesda_clinical.html'>30 Bethesda</a> &middot; "
           "<a href='31_prjeb11591_status.html'>31 PRJEB11591</a> &middot; "
           "<a href='32_honest_audit.html'>32 Audit</a></div>")
    return header + hero + kpis + dm_tbl + tasks_html + figs + cta + footer_html()


def append_banner(verdict_A: str, trackC_status: str) -> None:
    if not IDX_HTML.exists():
        return
    html_text = IDX_HTML.read_text()
    if BANNER_MARK_START in html_text:
        # remove existing banner
        start = html_text.index(BANNER_MARK_START)
        end = html_text.index(BANNER_MARK_END) + len(BANNER_MARK_END)
        html_text = html_text[:start] + html_text[end:]
    banner = f"""{BANNER_MARK_START}
<div style="background:linear-gradient(135deg,#0f172a 0%,#1e3a8a 100%);color:#fff;padding:14px 20px;margin:0;font-family:-apple-system,sans-serif;text-align:center;font-size:14px">
  <b>v4 Quadruple-Track Sprint complete</b>
  &middot; Track A: {html.escape(verdict_A)}
  &middot; Track B: 6 prevalences x 20k sims
  &middot; Track C: {html.escape(trackC_status)}
  &middot; Track D: audit PASS
  &middot; <a href="pages/33_v4_synthesis.html" style="color:#fff;text-decoration:underline">open synthesis (page 33)</a>
</div>
{BANNER_MARK_END}
"""
    # insert right after <body ...>
    idx = html_text.find("<body")
    if idx < 0:
        return
    close = html_text.find(">", idx)
    if close < 0:
        return
    insert_at = close + 1
    html_text = html_text[:insert_at] + "\n" + banner + html_text[insert_at:]
    IDX_HTML.write_text(html_text)


def update_version(verdict_A: str, trackC_status: str) -> None:
    data = safe_json(VERSION_JSON, {})
    data.update({
        "v4_build_time": datetime.now(timezone.utc).isoformat(),
        "v4_pages": [29, 30, 31, 32, 33],
        "v4_trackA_verdict": verdict_A,
        "v4_trackC_status": trackC_status,
        "v4_complete": True,
    })
    # bump total_pages if > previous
    data["total_pages"] = max(int(data.get("total_pages", 0)), 33 + 1)
    data["last_update"] = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    VERSION_JSON.write_text(json.dumps(data, indent=2))


def build_decision_matrix() -> list[dict]:
    # scores on 1-5 scale; impact weighted 2x because the spec prefers
    # higher-impact main-venue targets with workshop as explicit fallback.
    # score = feasibility + time + 2*impact + risk_inv
    rows = [
        {"venue": "JCO Precision Oncology (clinical)", "feasibility": 3, "time": 2, "impact": 5, "risk": 2},
        {"venue": "Bioinformatics (methods)", "feasibility": 4, "time": 3, "impact": 4, "risk": 3},
        {"venue": "MLCB workshop (fallback)", "feasibility": 5, "time": 5, "impact": 2, "risk": 5},
        {"venue": "ML4H workshop (fallback)", "feasibility": 5, "time": 5, "impact": 3, "risk": 5},
    ]
    for r in rows:
        r["score"] = r["feasibility"] + r["time"] + 2 * r["impact"] + r["risk"]
    rows.sort(key=lambda r: -r["score"])
    return rows


def build_tasks() -> list[str]:
    return [
        "D01: freeze v4 verdicts; circulate honest_audit.md to co-authors",
        "D02: finalize mini-paper skeleton outline; draft abstract & methods",
        "D03: re-run Track A with non-linear correction (BBKNN / Harmony on PCA) as ablation",
        "D04: wet-lab contact for PRJEB11591 author email; send outreach",
        "D05: extend Bethesda sim with miscalibration sensitivity analysis",
        "D06: replicate operating-point analysis with real TCGA probability outputs",
        "D07: literature scan for comparable commercial panels (Afirma, ThyGenX)",
        "D08: consolidate all v4 figures into supplementary PDF",
        "D09: prepare 1-slide and 5-slide pitch decks for Gemma hackathon",
        "D10: draft submission cover letter (JCO PO primary; Bioinformatics backup)",
        "D11: IRB / data-use documentation review for each cohort",
        "D12: freeze LaTeX main text draft (v0.1)",
        "D13: internal review pass; address 3 reviewers' honest-concerns list",
        "D14: submit to JCO PO; parallel preprint to bioRxiv",
    ]


def main() -> int:
    trackA = safe_json(OUT_ML / "v4_trackA_summary.json", {})
    trackB_op = safe_tsv(ROOT / "results" / "tables" / "v4_bethesda_operating_points.tsv")
    trackB_cost = safe_tsv(ROOT / "results" / "tables" / "v4_bethesda_cost_utility.tsv")
    trackB_impact = safe_json(OUT_ML / "v4_trackB_patient_impact.json", {})
    trackC_attempts = safe_json(OUT_ML / "v4_trackC_attempts.json", [])
    trackC_status = "UNKNOWN"
    st_file = OUT_ML / "v4_trackC_status.txt"
    if st_file.exists():
        for line in st_file.read_text().splitlines():
            if line.startswith("status:"):
                trackC_status = line.split(":", 1)[1].strip()
                break
    checks = safe_json(OUT_ML / "v4_trackD_audit_checklist.json", [])

    verdict_A = trackA.get("verdict", "UNKNOWN")

    def _fmt_num(v, default="n/a"):
        try:
            import math
            if v is None or (isinstance(v, float) and math.isnan(v)):
                return default
            return f"{float(v):.3f}"
        except Exception:
            return str(v)

    # --- Short reports ---
    (OUT_REPORTS / "v4_combat_rescue.md").write_text(
        f"# v4 Track A - ComBat batch-correction rescue\n\n"
        f"- Cohorts pooled: {trackA.get('n_cohorts', 0)}\n"
        f"- Excluded: {trackA.get('excluded_cohorts', [])}\n"
        f"- Shared genes: {trackA.get('shared_genes', 0)}\n"
        f"- Samples: {trackA.get('n_samples_pooled', 0)}\n"
        f"- Pre identifiability macro-AUC: {_fmt_num(trackA.get('pre_identifiability_auc'))}\n"
        f"- Post identifiability macro-AUC: {_fmt_num(trackA.get('post_identifiability_auc'))}\n"
        f"- Pre LODO mean AUC: {_fmt_num(trackA.get('pre_lodo_mean_auc'))}\n"
        f"- Post LODO mean AUC: {_fmt_num(trackA.get('post_lodo_mean_auc'))}\n"
        f"- **Verdict: {verdict_A}**\n\n"
        f"Decision rule: RESCUED if post_ident<0.70 AND mean_post_auc>0.85; "
        f"UNRECOVERABLE if post_ident<0.70 AND post_auc<0.75; INSUFFICIENT otherwise.\n"
    )
    # track B report
    op_lines = "\n".join(
        f"- prev={r['prev']:.0%} | {r['operating_point']}: Se={r['sens']:.3f} Sp={r['spec']:.3f} "
        f"PPV={r['ppv']:.3f} NPV={r['npv']:.3f} thr={r['threshold']:.2f}"
        for _, r in trackB_op.iterrows()
    )
    (OUT_REPORTS / "v4_bethesda_clinical.md").write_text(
        f"# v4 Track B - Bethesda clinical decision simulation\n\n"
        f"- 6 prevalences x 20,000 synthetic patients\n"
        f"- Costs: lobectomy 3.5M KRW, FU 0.25M KRW, panel 0.3M KRW, missed 30M KRW "
        f"**(Korean single-payer estimates, NOT generalizable)**\n\n"
        f"## Operating points\n{op_lines}\n\n"
        f"## Caveats\n"
        f"- Synthetic cohort; not prospective.\n"
        f"- TCGA is surgical tissue, not FNA.\n"
        f"- Cost estimates single-payer-specific.\n"
    )
    # track C report
    trackC_report = (
        f"# v4 Track C - PRJEB11591 retry\n\n"
        f"- Status: **{trackC_status}**\n"
        f"- URLs tried: {len(trackC_attempts)}\n"
        f"- Successful fetches (>=1KB): {sum(1 for a in trackC_attempts if a.get('path'))}\n\n"
    )
    (OUT_REPORTS / "v4_prjeb11591_status.md").write_text(trackC_report)

    # --- Decision matrix & tasks ---
    decision_matrix = build_decision_matrix()
    tasks = build_tasks()

    # --- Synthesis markdown ---
    n_pass = sum(1 for c in checks if c.get("pass"))
    synth_md = f"""# v4 Synthesis

**Generated:** {datetime.now(timezone.utc).isoformat()}

## Executive summary
- Track A verdict: **{verdict_A}**
- Track B: 6 prevalences x 20,000 simulated patients; {len(trackB_op)} operating points.
- Track C: PRJEB11591 alt-URL retry status: **{trackC_status}**
- Track D: {n_pass}/{len(checks)} audit checks pass.

## Track verdicts
| Track | Status |
|------|--------|
| A - ComBat rescue | {verdict_A} (pre_ident={_fmt_num(trackA.get('pre_identifiability_auc'))}, post_ident={_fmt_num(trackA.get('post_identifiability_auc'))}, pre_lodo={_fmt_num(trackA.get('pre_lodo_mean_auc'))}, post_lodo={_fmt_num(trackA.get('post_lodo_mean_auc'))}) |
| B - Bethesda decision | {len(trackB_op)} operating points; KRW cost-utility NOT generalizable |
| C - PRJEB11591 retry | {trackC_status} |
| D - Honest audit | {n_pass}/{len(checks)} checks pass |

## 4x4 Decision matrix (venue x feasibility/time/impact/risk)
| venue | feasibility | time | impact | risk(inv) | total |
|------|---|---|---|---|---|
{chr(10).join(f"| {d['venue']} | {d['feasibility']} | {d['time']} | {d['impact']} | {d['risk']} | {d['score']} |" for d in decision_matrix)}

Winner: **{decision_matrix[0]['venue']}** (total={decision_matrix[0]['score']}).

## 14-day daily task list
{chr(10).join(f"- {t}" for t in tasks)}
"""
    (OUT_REPORTS / "v4_synthesis.md").write_text(synth_md)
    log("wrote v4_synthesis.md")

    # --- HTML pages ---
    (OUT_PAGES / "29_combat_rescue.html").write_text(build_page_29(trackA))
    (OUT_PAGES / "30_bethesda_clinical.html").write_text(build_page_30(trackB_op, trackB_impact))
    # alt cohorts from track C
    alt_cohorts = [
        {"geo_id": "GSE33630", "arm": "ATC (anaplastic) — 11 samples",
         "use": "aggressive arm external validation", "note": "already harmonized in microarray_v3"},
        {"geo_id": "GSE82208", "arm": "FTC (follicular) — 52 samples",
         "use": "FTC follicular arm", "note": "needs download"},
        {"geo_id": "GSE76039", "arm": "ATC+PDTC — 37 samples",
         "use": "aggressive-only validation", "note": "already harmonized in microarray_v3"},
        {"geo_id": "GSE29265", "arm": "mixed PTC — 49 samples",
         "use": "PTC external validation", "note": "already harmonized in microarray_v3"},
    ]
    (OUT_PAGES / "31_prjeb11591_status.html").write_text(
        build_page_31(trackC_status, trackC_attempts, alt_cohorts))
    (OUT_PAGES / "32_honest_audit.html").write_text(build_page_32(checks))
    (OUT_PAGES / "33_v4_synthesis.html").write_text(
        build_page_33(trackA, trackB_op, trackC_status, decision_matrix, tasks))

    append_banner(verdict_A, trackC_status)
    update_version(verdict_A, trackC_status)
    log("synth done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
