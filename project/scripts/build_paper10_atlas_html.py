#!/usr/bin/env python3
"""Build paper10_atlas.html — rich visualization-first companion page.

Loads all P9 + P10 TSVs and embeds:
  - 17 figures
  - 9 sortable HTML tables
  - thick prose descriptions per section
"""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P9 = ROOT / "project/results/paper9_sl_first_pass"
P10 = ROOT / "project/results/paper10_atlas"
OUT_HTML = ROOT / "project/papers_hub_2026_05_04/paper10_atlas.html"

# ----- load everything -----
candidates = pd.read_csv(P9 / "candidate_targets.tsv", sep="\t")
ccle_state = pd.read_csv(P9 / "ccle_dm1_state.tsv", sep="\t")
ccle_corr = pd.read_csv(P9 / "ccle_dm1_target_corr.tsv", sep="\t")
dge = pd.read_csv(P9 / "dge_candidate_ranking.tsv", sep="\t")
cox_full = pd.read_csv(P9 / "tcga_thca_target_cox.tsv", sep="\t")
cox_dm1 = pd.read_csv(P9 / "tcga_thca_dm1_target_cox.tsv", sep="\t")
class_summary = pd.read_csv(P10 / "class_summary.tsv", sep="\t")
summary_p9 = json.loads((P9 / "summary.json").read_text())

# Class color tokens
CLASS_COLORS = {
    "NAD_salvage": "#7B1F2A", "NAD_de_novo_SR_partner": "#A04451",
    "JAK_STAT": "#34547A", "JAK_STAT_upstream": "#5E80B0",
    "Epigenetic": "#B8893C", "SFK": "#3C6B4F",
    "Ion_channel": "#3F7A8A", "DDR": "#962E2E",
    "Metabolism": "#8B5E00", "Metabolism_glycolysis": "#A5722B",
    "Metabolism_master": "#5C3D00", "TROP2_ADC_non_SL": "#5E5E5E",
    "Lineage_TF_anchor": "#0F1A2E",
}

def class_chip(c: str) -> str:
    color = CLASS_COLORS.get(c, "#888")
    return (f'<span class="cls-chip" style="background:{color}1a;'
            f'color:{color};border:1px solid {color}55">{c}</span>')


def fmt_p(p):
    if pd.isna(p):
        return ""
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_num(x, n=3):
    if pd.isna(x):
        return ""
    return f"{x:.{n}f}"


def render_table(df: pd.DataFrame, table_id: str, columns: list[tuple]) -> str:
    """columns: list of (key, label, fmt). fmt = 'cls' | 'num:N' | 'pval' | 'str'"""
    head = "<thead><tr>" + "".join(
        f'<th onclick="sortTable(\'{table_id}\',{i})">{lbl} ⇅</th>'
        for i, (_, lbl, _) in enumerate(columns)
    ) + "</tr></thead>"
    rows_html = []
    for _, r in df.iterrows():
        cells = []
        for key, _, fmt in columns:
            v = r.get(key, "")
            if pd.isna(v):
                cells.append("<td></td>")
                continue
            if fmt == "cls":
                cells.append(f"<td>{class_chip(str(v))}</td>")
            elif fmt == "pval":
                cells.append(f'<td data-sort="{v}">{fmt_p(v)}</td>')
            elif fmt.startswith("num:"):
                n = int(fmt.split(":")[1])
                cells.append(f'<td data-sort="{v}">{fmt_num(v, n)}</td>')
            elif fmt == "int":
                cells.append(f'<td data-sort="{v}">{int(v)}</td>')
            elif fmt == "bold":
                cells.append(f"<td><strong>{v}</strong></td>")
            else:
                cells.append(f"<td>{v}</td>")
        rows_html.append("<tr>" + "".join(cells) + "</tr>")
    body = "<tbody>" + "".join(rows_html) + "</tbody>"
    return (f'<table class="atlas-table" id="{table_id}">{head}{body}</table>')


# ====== Build tables ======
T1 = render_table(
    candidates,
    "T1_panel",
    [("gene", "Gene", "bold"), ("class", "Class", "cls")],
)

T2 = render_table(
    ccle_state[["line", "name", "oncotree", "DM1_like_score",
                 "TF_collapse_score_z", "DM1_high"]],
    "T2_ccle_state",
    [
        ("line", "CCLE line", "str"),
        ("name", "Short name", "str"),
        ("oncotree", "Oncotree", "str"),
        ("DM1_like_score", "DM1_like score", "num:3"),
        ("TF_collapse_score_z", "TF collapse z", "num:3"),
        ("DM1_high", "DM1 high?", "int"),
    ],
)

T3 = render_table(
    ccle_corr.assign(present=ccle_corr["pearson_r"].notna()).query("present"),
    "T3_ccle_corr",
    [
        ("gene", "Gene", "bold"),
        ("class", "Class", "cls"),
        ("pearson_r", "Pearson r", "num:3"),
        ("p", "p", "pval"),
        ("mean_dm1_high", "mean (DM1-high)", "num:2"),
        ("mean_dm1_low", "mean (DM1-low)", "num:2"),
        ("delta", "Δ (high − low)", "num:2"),
    ],
)

dge_view = dge.copy()
dge_view["dm1_minus_dm2"] = -dge_view["log2FC_DM2vDM1"]
T4 = render_table(
    dge_view.sort_values("dm1_minus_dm2", ascending=False)[
        ["gene", "class", "dm1_minus_dm2", "log2FC_DM2vDM1", "t", "q_BH",
         "mean_DM1", "mean_DM2"]
    ],
    "T4_dge",
    [
        ("gene", "Gene", "bold"),
        ("class", "Class", "cls"),
        ("dm1_minus_dm2", "log2(DM1/DM2)", "num:2"),
        ("log2FC_DM2vDM1", "log2(DM2/DM1)", "num:2"),
        ("t", "t-stat", "num:2"),
        ("q_BH", "q (BH)", "pval"),
        ("mean_DM1", "mean DM1 expr", "num:2"),
        ("mean_DM2", "mean DM2 expr", "num:2"),
    ],
)

T5 = render_table(
    cox_full.sort_values("p"),
    "T5_cox_full",
    [
        ("gene", "Gene", "bold"),
        ("class", "Class", "cls"),
        ("n", "n", "int"),
        ("hr", "HR", "num:3"),
        ("hr_lo", "HR low 95%", "num:3"),
        ("hr_hi", "HR high 95%", "num:3"),
        ("p", "p", "pval"),
    ],
)

T6 = render_table(
    cox_dm1.sort_values("p"),
    "T6_cox_dm1",
    [
        ("gene", "Gene", "bold"),
        ("class", "Class", "cls"),
        ("n", "n", "int"),
        ("hr", "HR (DM1-high)", "num:3"),
        ("hr_lo", "HR low 95%", "num:3"),
        ("hr_hi", "HR high 95%", "num:3"),
        ("p", "p", "pval"),
    ],
)

T7 = render_table(
    class_summary.sort_values("dge_mean_log2_dm1_up", ascending=False),
    "T7_class_summary",
    [
        ("class", "Class", "cls"),
        ("n_genes", "n genes", "int"),
        ("dge_mean_log2_dm1_up", "DGE mean log2(DM1/DM2)", "num:2"),
        ("dge_mean_t", "DGE mean t", "num:2"),
        ("ccle_n_measured", "CCLE measured", "int"),
        ("ccle_mean_r", "CCLE mean r", "num:3"),
        ("cox_n", "Cox n", "int"),
        ("cox_mean_log_hr", "Cox mean log HR", "num:3"),
        ("cox_min_p", "Cox min p", "pval"),
    ],
)

# T8 — drug landscape (manual structured)
drug_rows = [
    ("NAD salvage", "NAMPT (target) ↔ NAPRT (SR partner)",
     "NAMPT inhibitor",
     "FK866/APO866 · KPT-9274 · OT-82",
     "Phase I-II"),
    ("JAK/STAT upstream receptors", "OSMR · IL6R",
     "Anti-receptor + JAK/TYK2 inhibitor",
     "tocilizumab · ruxolitinib · deucravacitinib",
     "Approved (other indications)"),
    ("JAK/STAT", "JAK1 · JAK2 · TYK2 · STAT3",
     "Selective JAK / STAT3 inhibitor",
     "ruxolitinib · TTI-101 · napabucasin (caveat)",
     "Approved + clinical"),
    ("Epigenetic", "DNMT1/3A/3B · HDAC1/2/6 · KDM1A · EZH2",
     "DNMT inhibitor + HDACi + LSD1i + EZH2i",
     "decitabine · guadecitabine · vorinostat · GSK-2879552 · tazemetostat",
     "Approved + clinical"),
    ("SFK (SRC-family)", "LYN · SRC · FYN · YES1",
     "SRC-family inhibitor",
     "dasatinib · bosutinib · saracatinib",
     "Approved (CML/Ph+ALL)"),
    ("Ion channel", "KCNN4 (KCa3.1)",
     "Calcium-activated K+ channel blocker",
     "senicapoc (clinical) · TRAM-34 (tool)",
     "Phase II (other indications)"),
    ("DDR / replication stress", "ATR · CHEK1/2 · PARP1/2",
     "ATRi + PARPi (combination)",
     "ceralasertib · olaparib · talazoparib · niraparib",
     "Approved (BRCA) + clinical"),
    ("Metabolism", "GLS · LDHA · HK2 · IDH2 · MYC",
     "Glutaminase + glycolysis + OXPHOS inhibitor",
     "CB-839 · IACS-010759 · 2-DG · enasidenib (IDH2)",
     "Phase I-II + approved (IDH2)"),
    ("TROP2 (carved-out non-SL)", "TACSTD2",
     "Antibody-drug conjugate (ADC)",
     "sacituzumab govitecan (Trodelvy) · datopotamab deruxtecan",
     "Approved (TNBC, urothelial)"),
]
T8_rows = "".join(
    f'<tr><td><strong>{a}</strong></td><td>{b}</td><td>{c}</td>'
    f'<td>{d}</td><td>{e}</td></tr>'
    for a, b, c, d, e in drug_rows
)
T8 = (f'<table class="atlas-table" id="T8_drug">'
       f'<thead><tr><th>Class</th><th>Top targets</th><th>Drug class</th>'
       f'<th>Lead compounds</th><th>Stage</th></tr></thead>'
       f'<tbody>{T8_rows}</tbody></table>')

# T9 — deferred registry (manual structured)
deferred_rows = [
    ("D1", "DepMap CRISPR genome-wide × DM1 state",
     "Disk budget + matrix size", "Post-marathon"),
    ("D2", "SL / SR pair inference", "Requires D1 first", "After D1"),
    ("D3", "PRISM 19Q4 / 24Q2 IC50 stratification",
     "Raw cache on /opt/thyroid-dash, not on this dev box",
     "After Pod re-mount or re-pull"),
    ("D4", "GSE76039 per-target gene replication",
     "DM1 score on disk; per-gene matrix not pre-extracted",
     "1-day work post-marathon"),
    ("D5", "GDSC / CTRP cross-validation",
     "Out of first-pass scope", "Post D1"),
    ("D6", "Wet-lab Tier 1 (siRNA / drug viability)",
     "Requires collaborator wet-lab access",
     "Conditional on Paper 1 in print"),
    ("D7", "PDX / organoid validation",
     "Tier 3 — beyond first paper", "Follow-up"),
    ("D8", "Korean K2 (PRJEB11591) replication",
     "Korean pivot is Paper 4 territory",
     "Post Paper 1 + 2 + 4"),
]
T9_rows = "".join(
    f'<tr><td><strong>{a}</strong></td><td>{b}</td><td>{c}</td><td>{d}</td></tr>'
    for a, b, c, d in deferred_rows
)
T9 = (f'<table class="atlas-table" id="T9_deferred">'
       f'<thead><tr><th>#</th><th>Step</th><th>Why deferred</th>'
       f'<th>Earliest</th></tr></thead><tbody>{T9_rows}</tbody></table>')

# ====== Page shell ======
HTML = f"""<!DOCTYPE html>
<html lang="ko"><head>
<meta charset="UTF-8" /><meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Paper 10 — DM1 Synthetic-Lethal Vulnerability Atlas</title>
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;600&family=Noto+Sans+KR:wght@400;500;700&display=swap" rel="stylesheet" />
<link rel="stylesheet" href="style.css" />
<style>
.atlas-hero{{
  background:linear-gradient(135deg,#FAEEED 0%,#FCF8EE 50%,#E5EEF0 100%);
  border:1px solid var(--rule);border-left:8px solid #7B1F2A;
  padding:36px 44px;margin-bottom:30px;border-radius:6px
}}
.atlas-hero h1{{font-family:'Cormorant Garamond',serif;font-size:48px;line-height:1.05;color:#7B1F2A;letter-spacing:-1px}}
.atlas-hero h1 .num{{color:#0F1A2E;margin-right:14px;font-family:'JetBrains Mono',monospace;font-size:36px}}
.atlas-hero .sub{{color:#3A4658;font-size:17px;margin-top:14px;font-style:italic;line-height:1.55;max-width:920px}}
.atlas-hero .meta{{display:flex;gap:14px;margin-top:18px;font-family:'JetBrains Mono',monospace;font-size:12px;color:#3A4658;flex-wrap:wrap}}
.atlas-hero .meta span{{padding:4px 12px;background:rgba(0,0,0,0.04);border-radius:3px}}
.atlas-hero .meta strong{{color:#7B1F2A}}
.atlas-toc{{display:grid;grid-template-columns:repeat(3,1fr);gap:8px 16px;margin:18px 0 24px;background:#fff;border:1px solid var(--rule);padding:18px 22px;border-radius:5px}}
.atlas-toc a{{font-family:'JetBrains Mono',monospace;font-size:11.5px;color:var(--accent);text-decoration:none;padding:3px 0}}
.atlas-toc a:hover{{text-decoration:underline}}
section.atlas-sec{{margin:36px 0;padding-top:14px;border-top:1px dotted var(--rule)}}
section.atlas-sec h2{{font-family:'Cormorant Garamond',serif;font-size:30px;color:#0F1A2E;margin-bottom:6px;letter-spacing:-0.3px}}
section.atlas-sec .lead{{color:#3A4658;font-size:14.5px;font-style:italic;margin-bottom:14px;max-width:920px;line-height:1.6}}
section.atlas-sec .body{{font-size:15.5px;line-height:1.7;margin-bottom:14px}}
.atlas-figure{{margin:18px 0;padding:14px;background:#fff;border:1px solid var(--rule);border-radius:5px}}
.atlas-figure img{{width:100%;display:block;border:1px solid var(--rule-soft);border-radius:3px}}
.atlas-figure .cap{{margin-top:10px;font-size:13.5px;color:#3A4658;line-height:1.55}}
.atlas-figure .cap b{{color:#7B1F2A;font-family:'JetBrains Mono',monospace;font-size:11px;text-transform:uppercase;letter-spacing:.05em}}
.atlas-figure-pair{{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin:18px 0}}
.atlas-figure-trio{{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:18px 0}}
.atlas-table{{width:100%;border-collapse:collapse;font-size:12.5px;margin:12px 0;background:#fff}}
.atlas-table th,.atlas-table td{{border:1px solid var(--rule);padding:6px 10px;text-align:left;vertical-align:top}}
.atlas-table th{{background:#F4EFE6;font-family:'JetBrains Mono',monospace;font-size:10.5px;text-transform:uppercase;letter-spacing:.04em;cursor:pointer;user-select:none;white-space:nowrap}}
.atlas-table th:hover{{background:#EBE1C9}}
.atlas-table tbody tr:nth-child(odd){{background:#FCFAF5}}
.atlas-table tbody tr:hover{{background:#F4EFE6}}
.cls-chip{{display:inline-block;padding:2px 8px;border-radius:999px;font-family:'JetBrains Mono',monospace;font-size:10px;font-weight:700;letter-spacing:.04em}}
.kpi-strip{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:14px 0 22px}}
.kpi-card{{background:#fff;border:1px solid var(--rule);border-top:4px solid #3F7A8A;padding:14px 16px;border-radius:4px}}
.kpi-card.red{{border-top-color:#7B1F2A}}
.kpi-card.green{{border-top-color:#3C6B4F}}
.kpi-card.gold{{border-top-color:#B8893C}}
.kpi-card .k{{font-family:'JetBrains Mono',monospace;font-size:10px;color:#3A4658;text-transform:uppercase;letter-spacing:.06em}}
.kpi-card .v{{font-family:'Cormorant Garamond',serif;font-size:30px;color:#0F1A2E;font-weight:600;margin-top:4px}}
.kpi-card .s{{font-size:12.5px;color:#3A4658;margin-top:4px;line-height:1.45}}
.callout{{padding:14px 18px;margin:18px 0;border-left:4px solid #3F7A8A;background:#F4FAFB;font-size:14px;line-height:1.65}}
.callout.warn{{border-left-color:#B8893C;background:#FCF8EE}}
.callout.danger{{border-left-color:#7B1F2A;background:#FAEEED}}
.callout strong{{color:#7B1F2A}}
.print-btn{{position:fixed;top:16px;right:16px;z-index:9999;padding:10px 16px;background:#7B1F2A;color:#fff;border:none;border-radius:6px;font-family:'JetBrains Mono',monospace;font-size:13px;font-weight:600;cursor:pointer;box-shadow:0 2px 10px rgba(0,0,0,.18)}}
.print-btn:hover{{background:#000}}
@media (max-width:900px){{
  .atlas-toc{{grid-template-columns:1fr}}
  .kpi-strip{{grid-template-columns:1fr 1fr}}
  .atlas-figure-pair,.atlas-figure-trio{{grid-template-columns:1fr}}
}}
@media print{{
  .print-btn,.topnav,.atlas-toc{{display:none}}
  body{{background:#fff;color:#000}}
  section.atlas-sec{{page-break-inside:avoid}}
  .atlas-figure,.atlas-table{{page-break-inside:avoid}}
  .atlas-hero{{page-break-after:avoid}}
}}
</style>
</head><body>
<button class="print-btn" onclick="window.print()" title="이 페이지를 PDF로 인쇄/저장">📄 PDF로 저장 (Ctrl+P)</button>
<div class="container">
<nav class="topnav">
  <a href="index.html">← Portfolio hub</a>
  <span class="center">Paper 10 — DM1 SL Vulnerability Atlas</span>
  <a href="portfolio_paper9.html">Paper 9 detail →</a>
</nav>

<header class="atlas-hero">
  <h1><span class="num">10</span>DM1 Synthetic-Lethal Vulnerability Atlas</h1>
  <div class="sub">
    Visualization-first companion to Paper 9. Same DM1 anchor, same triangulation —
    expanded into <strong>17 atlas-grade figures</strong>, <strong>9 sortable tables</strong>,
    and per-class drug-landscape annotations.
    Designed for live inspection at <code style="background:rgba(0,0,0,0.06);padding:2px 8px;border-radius:3px">http://40.82.129.113/papers_hub_2026_05_04/paper10_atlas.html</code>.
  </div>
  <div class="meta">
    <span>Author: <strong>Seungho Cook</strong></span>
    <span>Anchor: <strong>Paper 1 DM1 / RAI-lineage axis</strong></span>
    <span>Date: 2026-05-04</span>
    <span>Status: <strong>strategic + first-pass</strong></span>
    <span>Earliest full-execution: 2026-06-16</span>
  </div>
</header>

<div class="kpi-strip">
  <div class="kpi-card red">
    <div class="k">Candidate panel</div>
    <div class="v">36 <span style="font-size:14px;color:#3A4658">genes</span></div>
    <div class="s">across 8 + 1 carved-out target classes</div>
  </div>
  <div class="kpi-card">
    <div class="k">DGE coverage</div>
    <div class="v">36/36</div>
    <div class="s">all candidates measured (TCGA-THCA bulk, 51,711 genes)</div>
  </div>
  <div class="kpi-card green">
    <div class="k">CCLE coverage</div>
    <div class="v">8/36</div>
    <div class="s">measurable in local 4007-gene thyroid panel (n=13)</div>
  </div>
  <div class="kpi-card gold">
    <div class="k">TCGA Cox</div>
    <div class="v">35/36</div>
    <div class="s">univariate Cox fit; n=560, 18 OS events (3.2%)</div>
  </div>
</div>

<div class="atlas-toc">
  <a href="#sec01">§ 01 Concept (F01)</a>
  <a href="#sec02">§ 02 CCLE state (F02 · T2)</a>
  <a href="#sec03">§ 03 CCLE corr per class (F03 · T3)</a>
  <a href="#sec04">§ 04 DGE volcano (F04–F06 · T4)</a>
  <a href="#sec05">§ 05 TCGA Cox forest (F07–F08 · T5–T6)</a>
  <a href="#sec06">§ 06 Triangulation (F09 · F11 · F12 · T7)</a>
  <a href="#sec07">§ 07 Top-5 KM curves (F10)</a>
  <a href="#sec08">§ 08 Lineage TF anchor (F13)</a>
  <a href="#sec09">§ 09 Drug landscape (F16 · T8)</a>
  <a href="#sec10">§ 10 Deferred registry (F17 · T9)</a>
  <a href="#sec11">§ 11 Provenance (F15)</a>
  <a href="#sec12">§ 12 Paper 9 ↔ 10 (F14)</a>
</div>

<section class="atlas-sec" id="sec01">
  <h2>§ 01 — Concept</h2>
  <div class="lead">DM1 = lineage-collapsed thyroid state (FOXE1 / NKX2-1 / PAX8 / HHEX ↓) with compensatory rewiring (STAT3 / AP-1 / DNMT / SFK / cytokine receptor / replication stress ↑). The Atlas asks: which compensatory nodes constitute druggable synthetic-lethal vulnerabilities?</div>
  <div class="body">
    <p>Paper 10 sits downstream of Paper 1 (which defines the DM1 state) and parallels Paper 9 (which lays out the strategic plan). The Atlas adds <strong>visualization-first triangulation</strong> across three local data layers — CCLE n=13, DM1-vs-DM2 DGE on TCGA-THCA bulk, and TCGA-THCA univariate Cox — so any reviewer can read effect size, direction, and provenance at a glance.</p>
    <p>The conceptual scaffolding is Ruppin-style (ISLE / SELECT) but does <em>not</em> lift their text, signatures, or pair tables. We cite their methods as prior art and run our own state-anchored analysis. The state anchor is <strong>transcriptional</strong>, not mutation-based — this is the core thyroid-specific adaptation.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F01_concept.png" alt="F01 concept"/>
    <div class="cap"><b>F01 — Concept.</b> DM1 lineage collapse → narrowed survival corridor → druggable SL/SR nodes. Triangulation across CCLE n=13, DGE 51,711 genes, and TCGA-THCA Cox n=560.</div>
  </div>
</section>

<section class="atlas-sec" id="sec02">
  <h2>§ 02 — CCLE thyroid: DM1_like state distribution (n=13)</h2>
  <div class="lead">Within-sample z-score over the 8 RAI genes (TPO/DIO1/TSHR/PAX8/TG/FOXE1/NKX2-1/SLC5A5), sign-flipped → DM1_like. Median split partitions 6 DM1-high vs 7 DM1-low lines.</div>
  <div class="body">
    <p>The CCLE n=13 panel is too small for any statistical claim about a single SL pair, but it provides a reproducible <strong>directional sanity check</strong> for the DM1 axis. Lines that are anaplastic / poorly differentiated (e.g., ML-1, MB1, TT, S117) cluster on the DM1-high side; the well-differentiated follicular lines (FTC133) and the BRAF V600E papillary line (BCPAP) cluster on the DM1-low side.</p>
    <p>The TF-collapse score (FOXE1/NKX2-1/PAX8/HHEX z-mean) tracks closely with DM1_like, as expected — these are not orthogonal axes but two views of the same lineage-collapse phenomenon.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F02_ccle_dm1_distribution.png" alt="F02 CCLE distribution"/>
    <div class="cap"><b>F02 — CCLE DM1 state.</b> 13 thyroid lines ranked by DM1_like score. Red = DM1-high (above median), blue = DM1-low. Gold dashed line = median split.</div>
  </div>
  <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:14px 0 6px">Table T2 — CCLE line state (sortable)</h3>
  {T2}
</section>

<section class="atlas-sec" id="sec03">
  <h2>§ 03 — CCLE candidate-target correlation, per class</h2>
  <div class="lead">Pearson r between each candidate's CCLE expression and the DM1_like score. Positive r = candidate up in DM1-like lines.</div>
  <div class="body">
    <p>Only 8 of 36 candidates are present in the local CCLE thyroid expression file (a thyroid-only 4007-gene matrix). For these 8, the directional signal points exactly where Paper 9's hypothesis predicts: <strong>NAMPT, TACSTD2, KCNN4, LYN positive</strong> (DM1-up); <strong>OSMR, PAX8, NKX2-1, FOXE1 negative</strong> (DM1-down — three of those four are lineage TF anchors and serve as the floor sanity check).</p>
    <p>Because n=13, all candidate p-values are nominal (only the lineage anchor FOXE1 reaches p&lt;1e-4). The full DepMap CRISPR contrast is the gating analysis for any SL claim — it is deferred (D1 in the registry).</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F03_ccle_corr_class_facet.png" alt="F03 CCLE corr facet"/>
    <div class="cap"><b>F03 — CCLE Pearson r per target class.</b> Facets by class. Limited to genes measurable in the local 4007-gene thyroid panel.</div>
  </div>
  <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:14px 0 6px">Table T3 — CCLE Pearson r (sortable, only present-in-CCLE genes)</h3>
  {T3}
</section>

<section class="atlas-sec" id="sec04">
  <h2>§ 04 — DM1 vs DM2 DGE: full candidate volcano + class summary</h2>
  <div class="lead">All 36 candidates measurable in the existing TCGA-THCA bulk DGE (51,711 genes). Volcano colored by class; per-class boxplot summarizes coherence.</div>
  <div class="body">
    <p>This is the <strong>strongest layer</strong> of the first-pass triangulation. Every candidate clears measurement; many clear q&lt;0.05 with substantial effect sizes. The picture that emerges is <strong>not</strong> "JAK/STAT core (STAT3/JAK1/JAK2) is the SL window" — STAT3 itself is essentially flat (q NS) at bulk. Instead, the JAK/STAT signal lives at <strong>upstream cytokine receptors</strong> (OSMR log2FC=−1.44 q=2e-15; IL6R log2FC=−1.34 q=5e-14). This is a meaningful re-framing for any drug-class strategy: anti-IL6R (tocilizumab-class) or anti-OSMR / dual-JAK becomes the candidate, not napabucasin.</p>
    <p>SFK is a <strong>coherent class signal</strong> — LYN, SRC, FYN all DM1-up at q&lt;<<0.05; YES1 alone goes the other way (DM2-up, q=0.011) which is consistent with thyroid-specific YES1 biology and supports treating YES1 as a class-internal control.</p>
    <p>KCNN4 (q=4e-32, log2FC=−4.00) and TACSTD2 (q=8e-24, log2FC=−4.14) are the two single-gene strongest DM1 markers — KCNN4 is in the SL track, TACSTD2 is the carved-out non-SL ADC track.</p>
    <p>Lineage TFs (FOXE1, HHEX, PAX8) are all DM2-enriched as expected — sanity floor confirmed.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F04_dge_volcano.png" alt="F04 volcano"/>
    <div class="cap"><b>F04 — DGE volcano.</b> All 36 candidates labeled; classes colored. q=0.05 cutoff dashed. Right side = DM1-enriched.</div>
  </div>
  <div class="atlas-figure-pair">
    <div class="atlas-figure">
      <img src="../results/paper10_atlas/figures/F05_dge_class_box.png" alt="F05 class box"/>
      <div class="cap"><b>F05 — Per-class DGE distribution.</b> Box per class (sorted by mean log2(DM1/DM2)). Coherence within class = candidates likely co-regulated.</div>
    </div>
    <div class="atlas-figure">
      <img src="../results/paper10_atlas/figures/F06_top50_dge.png" alt="F06 top 50"/>
      <div class="cap"><b>F06 — Whole-DGE top-50 DM1-up vs DM2-up.</b> Context for the candidate panel; shows the candidates are not an arbitrary cherry-pick.</div>
    </div>
  </div>
  <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:14px 0 6px">Table T4 — Full DGE candidate ranking (sortable)</h3>
  {T4}
</section>

<section class="atlas-sec" id="sec05">
  <h2>§ 05 — TCGA-THCA univariate Cox HR per candidate (whole + DM1-high)</h2>
  <div class="lead">Lifelines CoxPHFitter on UCSC pancan log2(norm+1) expression × OS. Whole cohort n=560, ~18 events; DM1-high subset n=280.</div>
  <div class="body">
    <p>TCGA-THCA's well-known limitation: very low OS event rate (~3.2%). Cox is therefore <strong>underpowered</strong> for most genes — the absence of significance is mostly a power statement, not evidence of irrelevance. Effect sizes from this layer are best treated as <em>directional priors</em> for downstream replication (GSE76039, Korean K2, Bundang).</p>
    <p>That said, three nominally significant hits emerge in the whole cohort: <strong>TYK2 HR=0.15 p=0.006</strong> (low TYK2 → worse OS — re-frame as protective modifier, not a target); <strong>ATR HR=6.88 p=0.024</strong> (high ATR → worse OS — fits replication-stress addiction hypothesis); <strong>CHEK2 HR=0.48 p=0.040</strong> (DDR loss → worse OS, direction-consistent with ATR finding). Together ATR + CHEK2 motivate the ATRi+PARPi combination as a state-conditional vulnerability even though the bulk DGE direction is flat.</p>
    <p>The DM1-high subset (n=280) is even more underpowered; no gene reaches p&lt;0.05. Top trends are MYC (p=0.05), HHEX (p=0.10), LDHA (p=0.10), GLS (p=0.15) — all consistent with a metabolic-axis story but each individual signal is weak.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F07_cox_forest_full.png" alt="F07 Cox forest"/>
    <div class="cap"><b>F07 — Whole-cohort Cox HR forest.</b> All 35 fittable candidates, sorted by HR. Class-colored. ✱ marks p&lt;0.05.</div>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F08_cox_whole_vs_dm1high.png" alt="F08 cox compare"/>
    <div class="cap"><b>F08 — Cox HR consistency.</b> log HR in whole TCGA-THCA vs DM1-high subset. Diagonal = invariant; deviation = state-conditional behavior.</div>
  </div>
  <div class="atlas-figure-pair">
    <div class="atlas-figure">
      <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:0 0 8px">Table T5 — whole TCGA-THCA Cox</h3>
      {T5}
    </div>
    <div class="atlas-figure">
      <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:0 0 8px">Table T6 — DM1-high subset Cox</h3>
      {T6}
    </div>
  </div>
</section>

<section class="atlas-sec" id="sec06">
  <h2>§ 06 — Triangulation: target class × evidence layer</h2>
  <div class="lead">Per-class summary across DGE / CCLE / Cox. Heatmap (F09) standardizes within-layer; scoreboard (F11) lays out verdicts; radar (F12) shows class-level shape.</div>
  <div class="body">
    <p>The triangulation matrix (F09) collapses 36 genes × 3 evidence layers into 12 classes × 3 z-scores. <strong>Ion_channel</strong>, <strong>JAK/STAT_upstream</strong>, <strong>SFK</strong>, <strong>TROP2</strong>, and <strong>NAD_salvage</strong> have positive DGE × positive CCLE — co-directional across the two transcriptomic layers. <strong>DDR</strong>, by contrast, is flat at bulk DGE but lights up in Cox — that is the state-conditional fingerprint we are looking for.</p>
    <p>The scoreboard (F11) names verdicts: KCNN4 (lead candidate), OSMR/IL6R (receptor-blockade SL), SFK (consistent class signal), DDR (state-conditional), MYC + metabolic (SL pair), NAMPT (needs DepMap), TROP2 (carved-out ADC).</p>
    <p>The radar (F12) shows the per-class shape: ion-channel is sharp DGE+ with weak Cox; DDR is sharp Cox with weak DGE; SFK is balanced.</p>
  </div>
  <div class="atlas-figure-trio">
    <div class="atlas-figure">
      <img src="../results/paper10_atlas/figures/F09_triangulation_matrix.png" alt="F09 matrix"/>
      <div class="cap"><b>F09 — Triangulation matrix.</b> Class × layer; cell text = raw value, color = within-layer z.</div>
    </div>
    <div class="atlas-figure">
      <img src="../results/paper10_atlas/figures/F11_triangulation_scoreboard.png" alt="F11 scoreboard"/>
      <div class="cap"><b>F11 — Scoreboard.</b> Provisional Paper 10 priority list with class verdicts and drug pointers.</div>
    </div>
    <div class="atlas-figure">
      <img src="../results/paper10_atlas/figures/F12_class_radar.png" alt="F12 radar"/>
      <div class="cap"><b>F12 — Class triangulation radar.</b> Per-axis min-max normalized; shape captures class strength profile across layers.</div>
    </div>
  </div>
  <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:14px 0 6px">Table T7 — Class summary (sortable)</h3>
  {T7}
</section>

<section class="atlas-sec" id="sec07">
  <h2>§ 07 — Top-5 Cox candidates: TCGA-THCA Kaplan-Meier</h2>
  <div class="lead">Top 5 candidates by whole-cohort Cox p-value, stratified by median expression. Logrank p reported.</div>
  <div class="body">
    <p>Median-split KM is more interpretable than continuous Cox for clinical readers. The top-5 panel (TYK2, ATR, CHEK2, YES1, HDAC1) recapitulates the Cox direction: high-ATR cohort separates worse than low-ATR cohort, etc. Logrank p values largely track Cox p but are typically slightly less significant due to the dichotomization.</p>
    <p>Note the y-axis range (0.6–1.0) — very few events occur, so visual separation between curves is small in absolute terms even when statistically suggestive.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F10_top5_KM.png" alt="F10 KM curves"/>
    <div class="cap"><b>F10 — Top-5 Kaplan-Meier panel.</b> TCGA-THCA OS, median-split per gene. Logrank p annotated on each panel.</div>
  </div>
</section>

<section class="atlas-sec" id="sec08">
  <h2>§ 08 — Lineage TF anchor sanity floor</h2>
  <div class="lead">FOXE1 / NKX2-1 / PAX8 / HHEX 4-panel: must be DM2-enriched and CCLE-negative-r if our DM1 axis is real.</div>
  <div class="body">
    <p>This is the primary check that the entire Atlas is built on a meaningful state. All four lineage TFs are DM2-enriched in DGE (FOXE1 log2FC=+1.04 q=8e-8; HHEX +0.73 q=1e-4; PAX8 +0.64 q=2e-4; NKX2-1 not reaching nominal q) and three of four (FOXE1, NKX2-1, PAX8) show negative Pearson r vs DM1_like in CCLE. <strong>FOXE1 r=−0.91 p=1.5e-5 in n=13</strong> is the cleanest single-gene calibration we have.</p>
    <p>If a future iteration ever shows lineage TFs going DM1-up, the entire DM1 axis is suspect — keep this sanity floor as a regression test.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F13_lineage_tf_anchor.png" alt="F13 lineage TF"/>
    <div class="cap"><b>F13 — Lineage TF anchor 4-panel.</b> Each TF: blue bar = CCLE r, red bar = DGE log2(DM1/DM2). All four TFs lineage-collapsed in DM1.</div>
  </div>
</section>

<section class="atlas-sec" id="sec09">
  <h2>§ 09 — Drug-class landscape</h2>
  <div class="lead">Per-class lead compounds and clinical stage. The intent is translational: Paper 10 picks classes where existing approved or clinical-stage agents already exist, lowering downstream cost.</div>
  <div class="body">
    <p>Approved or late-clinical-stage agents already exist for every major target class on the priority list — this is by design, not by chance. Composition-of-matter IP is mostly third-party; the Paper 10 IP angle is <strong>method-of-treatment</strong>: "method of selecting thyroid cancer patients for [target-class] therapy by measuring the 8-gene DM1 readout + TF-collapse score."</p>
    <p>The TROP2 row is preserved here as a reminder that <strong>ADC vulnerability is not synthetic lethality</strong>; the entry is included because TROP2 is the strongest DM1-surface marker after KCNN4 and is already being prosecuted in the Phase B drug platform.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F16_drug_landscape.png" alt="F16 drug landscape"/>
    <div class="cap"><b>F16 — Drug landscape figure.</b> Same content as Table T8 below in figure form for slide reuse.</div>
  </div>
  <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:14px 0 6px">Table T8 — Per-class drug landscape</h3>
  {T8}
</section>

<section class="atlas-sec" id="sec10">
  <h2>§ 10 — Deferred work registry (honest scope)</h2>
  <div class="lead">Eight deferred steps. None of these run before Paper 1 in print + marathon close (≥ 2026-06-13) + Yu sign-off.</div>
  <div class="callout warn">
    <strong>Scope rule.</strong> The Atlas is the maximum honest first-pass possible from local data without breaking the marathon-mode rule of "no new analysis unless paper-blocking for Paper 1". Every claim above is replicable with the script <code>p9_sl_first_pass.py</code> + <code>p10_atlas_figures.py</code>. Nothing in the deferred registry has been fudged or back-filled.
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F17_deferred_registry.png" alt="F17 deferred"/>
    <div class="cap"><b>F17 — Deferred registry figure.</b> 8 deferred steps with reason + earliest-start gate.</div>
  </div>
  <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;letter-spacing:.05em;margin:14px 0 6px">Table T9 — Deferred work registry</h3>
  {T9}
</section>

<section class="atlas-sec" id="sec11">
  <h2>§ 11 — Data provenance</h2>
  <div class="lead">Every cell in this Atlas traces back to a local file. No new download was performed for Paper 10.</div>
  <div class="body">
    <p>Three data layers feed all 17 figures and 9 tables:</p>
    <ul>
      <li><strong>CCLE thyroid (n=13)</strong> — local 4007-gene panel under <code>project/results/v8p1_rigor/f_ccle_validation/</code>, originally generated for the v8p1 BRS validation track.</li>
      <li><strong>DM1-vs-DM2 DGE (51,711 genes)</strong> — pre-existing TCGA-THCA bulk DGE under <code>project/results/dark_matter_phase2/web/data/dge_dm1_vs_dm2.tsv</code>; computed during the dark-matter phase 2 work.</li>
      <li><strong>TCGA pancan expression (11,069 samples)</strong> — UCSC Xena log2(norm+1) under <code>project/data/raw/TCGA_pancan/pancan_geneExp.gz</code>; merged with <code>project_external_st/results/extra/s_tcga_thca_scored.tsv</code> for DM1-state and OS metadata.</li>
    </ul>
    <p>Disk impact this session: 0 new downloads. All Atlas output written under <code>project/results/paper10_atlas/figures/</code> (~2.0 MB total).</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F15_data_provenance.png" alt="F15 provenance"/>
    <div class="cap"><b>F15 — Provenance diagram.</b> Three sources upstream; eight deferred steps downstream; Atlas in the middle.</div>
  </div>
</section>

<section class="atlas-sec" id="sec12">
  <h2>§ 12 — Paper 9 ↔ Paper 10 positioning</h2>
  <div class="lead">Same DM1 anchor, same triangulation. Different communication mode. Both gated on Paper 1 in print + marathon close.</div>
  <div class="body">
    <p>Paper 9 = strategic plan + first-pass results, optimized for written argument. 14 sections of text, deferred-work logic, IP angle, 12-week execution plan. Best for reviewer / advisor consumption when sequence-of-reasoning matters.</p>
    <p>Paper 10 = visualization-first ATLAS, optimized for inspection. 17 figures + 9 sortable tables + thick captions. Best for in-meeting reference, slide reuse, and quick lookup.</p>
    <p>Both papers point at the same execution plan: full DepMap CRISPR + PRISM + GSE76039 replication start ≥ 2026-06-16.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper10_atlas/figures/F14_paper9_vs_paper10.png" alt="F14 positioning"/>
    <div class="cap"><b>F14 — Paper 9 ↔ Paper 10.</b> Both anchored on Paper 1 DM1/RAI-lineage axis. Paper 9 feeds Paper 10's figure inventory.</div>
  </div>
</section>

<footer style="margin-top:48px;padding-top:24px;border-top:2px solid #7B1F2A;font-size:13px;color:#3A4658">
  <p><strong>Reproducibility.</strong> All figures regenerable via <code>source project/.venv/bin/activate &amp;&amp; python project/scripts/p10_atlas_figures.py</code>; tables regenerable via <code>python project/scripts/build_paper10_atlas_html.py</code>.</p>
  <p><strong>What is NOT claimed.</strong> No clinical efficacy. No SL-pair claim (deferred). No ISLE/SELECT equivalence (concept-only borrowing). No replacement of driver-mutation profiling. TROP2 explicitly NOT framed as SL.</p>
  <p><strong>Scope guard.</strong> Voice-protected prose not authored. Paper 1 / 2 / 3 / 4 manuscripts not touched. No new data downloaded. No GPU. No wet-lab.</p>
  <p style="margin-top:14px;font-style:italic;color:#7B1F2A">Hub URL: <code style="background:rgba(0,0,0,0.06);padding:2px 8px;border-radius:3px">http://40.82.129.113/papers_hub_2026_05_04/paper10_atlas.html</code></p>
</footer>

</div>

<script>
function sortTable(tableId, col) {{
  const tbl = document.getElementById(tableId);
  const tbody = tbl.querySelector('tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const dir = tbl.dataset.sortDir === 'asc' && tbl.dataset.sortCol == col ? 'desc' : 'asc';
  rows.sort((a, b) => {{
    const ac = a.children[col];
    const bc = b.children[col];
    const av = ac.dataset.sort !== undefined ? parseFloat(ac.dataset.sort) : ac.textContent.trim();
    const bv = bc.dataset.sort !== undefined ? parseFloat(bc.dataset.sort) : bc.textContent.trim();
    if (typeof av === 'number' && typeof bv === 'number' && !isNaN(av) && !isNaN(bv)) {{
      return dir === 'asc' ? av - bv : bv - av;
    }}
    return dir === 'asc'
      ? String(av).localeCompare(String(bv))
      : String(bv).localeCompare(String(av));
  }});
  rows.forEach(r => tbody.appendChild(r));
  tbl.dataset.sortDir = dir;
  tbl.dataset.sortCol = col;
}}
</script>
</body></html>"""

OUT_HTML.write_text(HTML)
print(f"Wrote {OUT_HTML} ({OUT_HTML.stat().st_size:,} bytes)")
