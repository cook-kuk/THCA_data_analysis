#!/usr/bin/env python3
"""Build paper11_pancancer.html + paper12_network.html — rich atlas pages."""
from pathlib import Path
import json

import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
HUB = ROOT / "project/papers_hub_2026_05_04"
P11 = ROOT / "project/results/paper11_pancancer"
P12 = ROOT / "project/results/paper12_network"

CLASS_COLORS = {
    "NAD_salvage": "#7B1F2A", "NAD_de_novo_SR_partner": "#A04451",
    "JAK_STAT": "#34547A", "JAK_STAT_upstream": "#5E80B0",
    "Epigenetic": "#B8893C", "SFK": "#3C6B4F",
    "Ion_channel": "#3F7A8A", "DDR": "#962E2E",
    "Metabolism": "#8B5E00", "Metabolism_glycolysis": "#A5722B",
    "Metabolism_master": "#5C3D00", "TROP2_ADC_non_SL": "#5E5E5E",
    "Lineage_TF_anchor": "#0F1A2E",
}


def cls_chip(c):
    color = CLASS_COLORS.get(c, "#888")
    return (f'<span style="display:inline-block;padding:2px 8px;'
            f'border-radius:999px;background:{color}1a;color:{color};'
            f'border:1px solid {color}55;font-family:\'JetBrains Mono\',monospace;'
            f'font-size:10px;font-weight:700">{c}</span>')


# ===== Atlas page CSS shell =====
COMMON_CSS = """<style>
.atlas-hero{background:linear-gradient(135deg,#E5EEF0 0%,#FCF8EE 50%,#FAEEED 100%);border:1px solid var(--rule);border-left:8px solid #3F7A8A;padding:36px 44px;margin-bottom:30px;border-radius:6px}
.atlas-hero.acc-red{border-left-color:#7B1F2A;background:linear-gradient(135deg,#FAEEED 0%,#FCF8EE 50%,#E5EEF0 100%)}
.atlas-hero h1{font-family:'Cormorant Garamond',serif;font-size:46px;line-height:1.05;color:#3F7A8A;letter-spacing:-1px}
.atlas-hero.acc-red h1{color:#7B1F2A}
.atlas-hero h1 .num{color:#0F1A2E;margin-right:14px;font-family:'JetBrains Mono',monospace;font-size:34px}
.atlas-hero .sub{color:#3A4658;font-size:16px;margin-top:14px;font-style:italic;line-height:1.55;max-width:920px}
.atlas-hero .meta{display:flex;gap:14px;margin-top:18px;font-family:'JetBrains Mono',monospace;font-size:12px;color:#3A4658;flex-wrap:wrap}
.atlas-hero .meta span{padding:4px 12px;background:rgba(0,0,0,0.04);border-radius:3px}
.atlas-toc{display:grid;grid-template-columns:repeat(3,1fr);gap:8px 16px;margin:18px 0 24px;background:#fff;border:1px solid var(--rule);padding:18px 22px;border-radius:5px}
.atlas-toc a{font-family:'JetBrains Mono',monospace;font-size:11.5px;color:var(--accent);text-decoration:none;padding:3px 0}
.atlas-toc a:hover{text-decoration:underline}
section.atlas-sec{margin:36px 0;padding-top:14px;border-top:1px dotted var(--rule)}
section.atlas-sec h2{font-family:'Cormorant Garamond',serif;font-size:30px;color:#0F1A2E;margin-bottom:6px;letter-spacing:-0.3px}
section.atlas-sec .lead{color:#3A4658;font-size:14.5px;font-style:italic;margin-bottom:14px;max-width:920px;line-height:1.6}
section.atlas-sec .body{font-size:15.5px;line-height:1.7;margin-bottom:14px}
.atlas-figure{margin:18px 0;padding:14px;background:#fff;border:1px solid var(--rule);border-radius:5px}
.atlas-figure img{width:100%;display:block;border:1px solid var(--rule-soft);border-radius:3px}
.atlas-figure .cap{margin-top:10px;font-size:13.5px;color:#3A4658;line-height:1.55}
.atlas-figure .cap b{color:#7B1F2A;font-family:'JetBrains Mono',monospace;font-size:11px;text-transform:uppercase;letter-spacing:.05em}
.atlas-figure-pair{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin:18px 0}
.atlas-table{width:100%;border-collapse:collapse;font-size:12.5px;margin:12px 0;background:#fff}
.atlas-table th,.atlas-table td{border:1px solid var(--rule);padding:6px 10px;text-align:left;vertical-align:top}
.atlas-table th{background:#F4EFE6;font-family:'JetBrains Mono',monospace;font-size:10.5px;text-transform:uppercase;letter-spacing:.04em;cursor:pointer;user-select:none}
.atlas-table th:hover{background:#EBE1C9}
.atlas-table tbody tr:nth-child(odd){background:#FCFAF5}
.atlas-table tbody tr:hover{background:#F4EFE6}
.kpi-strip{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:14px 0 22px}
.kpi-card{background:#fff;border:1px solid var(--rule);border-top:4px solid #3F7A8A;padding:14px 16px;border-radius:4px}
.kpi-card.red{border-top-color:#7B1F2A}
.kpi-card.green{border-top-color:#3C6B4F}
.kpi-card.gold{border-top-color:#B8893C}
.kpi-card .k{font-family:'JetBrains Mono',monospace;font-size:10px;color:#3A4658;text-transform:uppercase}
.kpi-card .v{font-family:'Cormorant Garamond',serif;font-size:30px;color:#0F1A2E;font-weight:600;margin-top:4px}
.kpi-card .s{font-size:12.5px;color:#3A4658;margin-top:4px;line-height:1.45}
.callout{padding:14px 18px;margin:18px 0;border-left:4px solid #3F7A8A;background:#F4FAFB;font-size:14px;line-height:1.65}
.print-btn{position:fixed;top:16px;right:16px;z-index:9999;padding:10px 16px;background:#3F7A8A;color:#fff;border:none;border-radius:6px;font-family:'JetBrains Mono',monospace;font-size:13px;font-weight:600;cursor:pointer;box-shadow:0 2px 10px rgba(0,0,0,.18)}
.print-btn:hover{background:#000}
@media (max-width:900px){.atlas-toc{grid-template-columns:1fr}.kpi-strip{grid-template-columns:1fr 1fr}.atlas-figure-pair{grid-template-columns:1fr}}
@media print{.print-btn,.topnav,.atlas-toc{display:none}body{background:#fff;color:#000}.atlas-figure,.atlas-table{page-break-inside:avoid}}
</style>"""

SORT_JS = """<script>
function sortTable(tableId, col){
  const tbl=document.getElementById(tableId);
  const tbody=tbl.querySelector('tbody');
  const rows=Array.from(tbody.querySelectorAll('tr'));
  const dir=tbl.dataset.sortDir==='asc'&&tbl.dataset.sortCol==col?'desc':'asc';
  rows.sort((a,b)=>{
    const ac=a.children[col],bc=b.children[col];
    const av=ac.dataset.sort!==undefined?parseFloat(ac.dataset.sort):ac.textContent.trim();
    const bv=bc.dataset.sort!==undefined?parseFloat(bc.dataset.sort):bc.textContent.trim();
    if(typeof av==='number'&&typeof bv==='number'&&!isNaN(av)&&!isNaN(bv))return dir==='asc'?av-bv:bv-av;
    return dir==='asc'?String(av).localeCompare(String(bv)):String(bv).localeCompare(String(av));
  });
  rows.forEach(r=>tbody.appendChild(r));
  tbl.dataset.sortDir=dir;tbl.dataset.sortCol=col;
}
</script>"""


# ============================================================
#  PAPER 11
# ============================================================
print("Building paper11_pancancer.html …")
ls = pd.read_csv(P11 / "pancan_lineage_stats.tsv", sep="\t")
corr = pd.read_csv(P11 / "lineage_candidate_corr.tsv", sep="\t").rename(
    columns={"Unnamed: 0": "lineage"}
) if (P11 / "lineage_candidate_corr.tsv").exists() else pd.DataFrame()
sm11 = json.loads((P11 / "summary.json").read_text())

# Build T11_lineage table
t11_rows = []
for _, r in ls.iterrows():
    is_thca = r["lineage"] == "thyroid carcinoma"
    bg = ' style="background:#FAEEED;font-weight:bold"' if is_thca else ''
    t11_rows.append(
        f'<tr{bg}><td>{r["lineage"]}</td>'
        f'<td data-sort="{r["n"]}">{int(r["n"])}</td>'
        f'<td data-sort="{r["median_DM1"]}">{r["median_DM1"]:+.3f}</td>'
        f'<td data-sort="{r["mean_DM1"]}">{r["mean_DM1"]:+.3f}</td>'
        f'<td data-sort="{r["pct_above_thca_med"]}">'
        f'{r["pct_above_thca_med"]*100:.1f}%</td>'
        f'<td data-sort="{r["pct_above_pan_med"]}">'
        f'{r["pct_above_pan_med"]*100:.1f}%</td></tr>'
    )
T11 = (
    '<table class="atlas-table" id="T11_lineage">'
    '<thead><tr>'
    + "".join(f'<th onclick="sortTable(\'T11_lineage\',{i})">{l} ⇅</th>'
              for i, l in enumerate([
                  "Lineage", "n", "Median DM1", "Mean DM1",
                  "% > THCA median", "% > pan-cancer median"]))
    + '</tr></thead><tbody>'
    + "".join(t11_rows)
    + "</tbody></table>"
)


HTML11 = f"""<!DOCTYPE html>
<html lang="ko"><head>
<meta charset="UTF-8" /><meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Paper 11 — Pan-Cancer DM1 Transfer Atlas</title>
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;600&family=Noto+Sans+KR:wght@400;500;700&display=swap" rel="stylesheet" />
<link rel="stylesheet" href="style.css" />
{COMMON_CSS}
</head><body>
<button class="print-btn" onclick="window.print()" title="이 페이지를 PDF로 인쇄/저장">📄 PDF로 저장 (Ctrl+P)</button>
<div class="container">
<nav class="topnav">
  <a href="index.html">← Portfolio hub</a>
  <span class="center">Paper 11 — Pan-Cancer DM1 Transfer</span>
  <a href="paper10_atlas.html">Paper 10 ATLAS →</a>
</nav>

<header class="atlas-hero">
  <h1><span class="num">11</span>Pan-Cancer DM1 Transfer Atlas</h1>
  <div class="sub">
    Score the DM1_like axis (8-gene RAI within-sample-z, sign-flipped) across all <strong>{sm11['pancan_n_samples']:,}</strong> TCGA pancan samples · <strong>{sm11['n_lineages']}</strong> lineages.
    Where else outside thyroid does this DM1 phenotype appear, and is the candidate-target panel portable to those non-thyroid contexts?
  </div>
  <div class="meta">
    <span>Author: <strong>Seungho Cook</strong></span>
    <span>Anchor: <strong>Paper 1 DM1 axis · 8-gene RAI signature</strong></span>
    <span>Substrate: <strong>TCGA pancan log2(norm+1)</strong></span>
    <span>Date: 2026-05-04</span>
  </div>
</header>

<div class="kpi-strip">
  <div class="kpi-card red"><div class="k">Pancan samples scored</div><div class="v">{sm11['pancan_n_samples']:,}</div><div class="s">across {sm11['n_lineages']} TCGA primary diseases</div></div>
  <div class="kpi-card"><div class="k">THCA reference median</div><div class="v">{sm11['thca_median_dm1']:+.2f}</div><div class="s">DM1_like z-score (sign-flipped RAI_8 mean z)</div></div>
  <div class="kpi-card green"><div class="k">Top 5 non-THCA DM1+</div><div class="v">5</div><div class="s">{', '.join(sm11['top5_non_thca_dm1'][:3])}…</div></div>
  <div class="kpi-card gold"><div class="k">Top portability</div><div class="v">{len(sm11['top_portability'])}</div><div class="s">non-thyroid lineages with DM1+ AND TF-collapse</div></div>
</div>

<div class="atlas-toc">
  <a href="#sec01">§ 01 Lineage DM1 distribution (F11.01)</a>
  <a href="#sec02">§ 02 % above THCA median (F11.02)</a>
  <a href="#sec03">§ 03 Sample size × DM1 (F11.03)</a>
  <a href="#sec04">§ 04 Cross-lineage candidate heatmap (F11.04)</a>
  <a href="#sec05">§ 05 DM1 vs TF-collapse (F11.05)</a>
  <a href="#sec06">§ 06 Portability ranking (F11.06)</a>
  <a href="#sec07">§ 07 Lineage table (T11)</a>
</div>

<section class="atlas-sec" id="sec01">
  <h2>§ 01 — Per-lineage DM1_like distribution</h2>
  <div class="lead">Boxplot per TCGA primary disease, sorted by median DM1_like. THCA highlighted in dark red.</div>
  <div class="body">
    <p>The 8-gene RAI signature is a thyroid-developmental signature, so by construction THCA samples score high on DM1_like. The interesting question is which non-thyroid lineages have a tail or shoulder of high DM1_like samples, since those would be candidates for SL strategies developed against the DM1 axis.</p>
    <p>Reading the figure: sarcomas, gliomas, thymic and adrenal lineages can show elevated DM1_like in subsets — these are partial-lineage-loss cancers, plausible portability targets. Differentiated tissue cancers (e.g., prostate, well-differentiated PTC reference) sit at the other end.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper11_pancancer/figures/F11_01_lineage_dm1_box.png" alt="F11.01"/>
    <div class="cap"><b>F11.01 — Pan-cancer DM1_like distribution per lineage.</b> Sorted by lineage median; THCA in dark red. Dashed line = THCA median.</div>
  </div>
</section>

<section class="atlas-sec" id="sec02">
  <h2>§ 02 — % of lineage above THCA median</h2>
  <div class="lead">Threshold = THCA median DM1_like. Lineages with high prevalence may share the DM1 vulnerability state.</div>
  <div class="atlas-figure">
    <img src="../results/paper11_pancancer/figures/F11_02_lineage_pct_above_thca.png" alt="F11.02"/>
    <div class="cap"><b>F11.02 — Cross-cancer DM1-like prevalence.</b> Bar color: red ≥ 40%, blue otherwise. THCA flagged.</div>
  </div>
</section>

<section class="atlas-sec" id="sec03">
  <h2>§ 03 — Lineage sample size vs median DM1_like</h2>
  <div class="lead">Big lineages (n>500) with median DM1_like outside [-0.1, +0.1] = high statistical power for portability test.</div>
  <div class="atlas-figure">
    <img src="../results/paper11_pancancer/figures/F11_03_lineage_n_vs_dm1.png" alt="F11.03"/>
    <div class="cap"><b>F11.03 — Lineage n × median DM1.</b> Labels for outliers and THCA reference.</div>
  </div>
</section>

<section class="atlas-sec" id="sec04">
  <h2>§ 04 — Cross-lineage candidate × DM1 correlation heatmap</h2>
  <div class="lead">Top 5 non-THCA DM1+ lineages + THCA reference. Per-lineage Pearson r between each candidate and DM1_like.</div>
  <div class="body">
    <p>If the SL strategy is portable, candidates should show <strong>concordant Pearson direction across multiple lineages</strong>. This is a stricter test than just "DM1+ in the lineage" — it asks whether the same compensatory rewiring (KCNN4 ↑, OSMR ↑, SFK ↑, etc.) is recovered.</p>
    <p>Class-color tick labels link back to Paper 9/10 candidate panel.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper11_pancancer/figures/F11_04_lineage_candidate_heatmap.png" alt="F11.04"/>
    <div class="cap"><b>F11.04 — Lineage × candidate correlation heatmap.</b> Red = candidate up in DM1, blue = down. Tick label colors = candidate class.</div>
  </div>
</section>

<section class="atlas-sec" id="sec05">
  <h2>§ 05 — DM1_like vs TF-collapse, per lineage</h2>
  <div class="lead">Two views of the same lineage-collapse phenomenon. A lineage with high DM1_like AND low TF-collapse is the strongest portability candidate.</div>
  <div class="atlas-figure">
    <img src="../results/paper11_pancancer/figures/F11_05_dm1_vs_tfcollapse_lineage.png" alt="F11.05"/>
    <div class="cap"><b>F11.05 — Lineage median: DM1_like and TF-collapse, side by side.</b></div>
  </div>
</section>

<section class="atlas-sec" id="sec06">
  <h2>§ 06 — Portability ranking</h2>
  <div class="lead">portability_score = median DM1_like − median TF-collapse. High = both axes signal lineage collapse → strongest candidate for SL transfer.</div>
  <div class="atlas-figure">
    <img src="../results/paper11_pancancer/figures/F11_06_portability.png" alt="F11.06"/>
    <div class="cap"><b>F11.06 — Top non-thyroid lineages by SL portability score.</b></div>
  </div>
</section>

<section class="atlas-sec" id="sec07">
  <h2>§ 07 — Lineage table T11 (sortable)</h2>
  <div class="lead">All {sm11['n_lineages']} TCGA primary diseases with DM1_like + TF-collapse stats. THCA row highlighted.</div>
  {T11}
</section>

<div class="callout">
<strong>What this Atlas does NOT claim.</strong> Pan-cancer DM1+ subsets do not automatically inherit thyroid SL biology. Portability requires (a) shared compensatory rewiring direction (F11.04) and (b) state-conditional dependency confirmation in DepMap — the latter is deferred. Treat this as a <em>hypothesis-generating ranking</em>, not as a clinical recommendation.
</div>

<footer style="margin-top:48px;padding-top:24px;border-top:2px solid #3F7A8A;font-size:13px;color:#3A4658">
  <p><strong>Reproducibility.</strong> Figures regenerable via <code>python project/scripts/p11_p12_extended_atlas.py</code>; HTML via <code>python project/scripts/build_paper11_12_html.py</code>.</p>
  <p><strong>Hub URL.</strong> <code style="background:rgba(0,0,0,0.06);padding:2px 8px;border-radius:3px">http://40.82.129.113/papers_hub_2026_05_04/paper11_pancancer.html</code></p>
</footer>
</div>
{SORT_JS}
</body></html>"""

(HUB / "paper11_pancancer.html").write_text(HTML11)
print(f"  Wrote paper11_pancancer.html "
      f"({(HUB / 'paper11_pancancer.html').stat().st_size:,} bytes)")


# ============================================================
#  PAPER 12
# ============================================================
print("Building paper12_network.html …")
sm12 = json.loads((P12 / "summary.json").read_text())
mod_summary = pd.read_csv(P12 / "module_summary.tsv", sep="\t")
modules = pd.read_csv(P12 / "candidate_modules.tsv", sep="\t")
edges = pd.read_csv(P12 / "edges_top.tsv", sep="\t")

# T12_modules
t12_mod_rows = []
for _, r in mod_summary.iterrows():
    members = modules[modules["module"] == r["module"]]
    chips = " ".join(cls_chip(c) for c in members["class"].unique())
    t12_mod_rows.append(
        f'<tr><td><strong>M{int(r["module"])}</strong></td>'
        f'<td>{int(r["n_genes"])}</td>'
        f'<td data-sort="{r["mean_intra_corr"]}">'
        f'{r["mean_intra_corr"]:+.3f}</td>'
        f'<td>{", ".join(members["gene"])}</td>'
        f'<td>{chips}</td></tr>'
    )
T12_mod = (
    '<table class="atlas-table" id="T12_mod">'
    '<thead><tr>'
    + "".join(f'<th onclick="sortTable(\'T12_mod\',{i})">{l} ⇅</th>'
              for i, l in enumerate(
                  ["Module", "n genes", "Mean intra r", "Members", "Classes"]))
    + '</tr></thead><tbody>' + "".join(t12_mod_rows)
    + "</tbody></table>"
)

# T12_edges (top 30)
t12_edge_rows = []
top_edges = edges.reindex(
    edges["r"].abs().sort_values(ascending=False).index
).head(30)

# Need candidate class lookup
CANDIDATES = dict(zip(modules["gene"], modules["class"]))

for i, (_, e) in enumerate(top_edges.iterrows(), 1):
    c1 = CANDIDATES.get(e["g1"], "")
    c2 = CANDIDATES.get(e["g2"], "")
    t12_edge_rows.append(
        f'<tr><td>{i}</td><td><strong>{e["g1"]}</strong> {cls_chip(c1)}</td>'
        f'<td><strong>{e["g2"]}</strong> {cls_chip(c2)}</td>'
        f'<td data-sort="{e["r"]}">{e["r"]:+.3f}</td></tr>'
    )
T12_edge = (
    '<table class="atlas-table" id="T12_edge">'
    '<thead><tr>'
    + "".join(f'<th onclick="sortTable(\'T12_edge\',{i})">{l} ⇅</th>'
              for i, l in enumerate(
                  ["#", "Gene 1", "Gene 2", "Pearson r"]))
    + '</tr></thead><tbody>' + "".join(t12_edge_rows)
    + "</tbody></table>"
)


HTML12 = f"""<!DOCTYPE html>
<html lang="ko"><head>
<meta charset="UTF-8" /><meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Paper 12 — Candidate Co-expression Network</title>
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;600&family=Noto+Sans+KR:wght@400;500;700&display=swap" rel="stylesheet" />
<link rel="stylesheet" href="style.css" />
{COMMON_CSS}
</head><body>
<button class="print-btn" onclick="window.print()" title="이 페이지를 PDF로 인쇄/저장">📄 PDF로 저장 (Ctrl+P)</button>
<div class="container">
<nav class="topnav">
  <a href="index.html">← Portfolio hub</a>
  <span class="center">Paper 12 — Candidate Co-expression Network</span>
  <a href="paper11_pancancer.html">Paper 11 →</a>
</nav>

<header class="atlas-hero acc-red">
  <h1><span class="num">12</span>Candidate Co-expression Network</h1>
  <div class="sub">
    Gene × gene Pearson correlation among the <strong>{sm12['n_candidates_in_network']}</strong> Paper 9/10 candidate genes in TCGA-THCA (n={sm12['n_thca_samples']}).
    Identifies <strong>{sm12['n_modules']} co-expression modules</strong> and the network's <strong>hub genes</strong>.
    Bridges the candidate panel to a regulatory-module reading.
  </div>
  <div class="meta">
    <span>Author: <strong>Seungho Cook</strong></span>
    <span>Substrate: <strong>TCGA-THCA n={sm12['n_thca_samples']}</strong></span>
    <span>Network edges (|r|≥0.4): <strong>{sm12['n_top_edges']}</strong></span>
    <span>Top hubs: <strong>{', '.join(sm12['top_hub_genes'])}</strong></span>
  </div>
</header>

<div class="kpi-strip">
  <div class="kpi-card red"><div class="k">Candidates in network</div><div class="v">{sm12['n_candidates_in_network']}</div><div class="s">of 36 panel genes (those measurable in pancan)</div></div>
  <div class="kpi-card"><div class="k">Co-expression modules</div><div class="v">{sm12['n_modules']}</div><div class="s">via average-linkage clustering on |1−r|</div></div>
  <div class="kpi-card green"><div class="k">Strong edges (|r|≥0.4)</div><div class="v">{sm12['n_top_edges']}</div><div class="s">edges in the network graph</div></div>
  <div class="kpi-card gold"><div class="k">Top 5 hub genes</div><div class="v" style="font-size:18px;line-height:1.3">{'<br>'.join(sm12['top_hub_genes'])}</div><div class="s">highest degree centrality</div></div>
</div>

<div class="atlas-toc">
  <a href="#sec01">§ 01 Co-expression matrix (F12.01)</a>
  <a href="#sec02">§ 02 Module summary (F12.02 · T12.M)</a>
  <a href="#sec03">§ 03 Network graph (F12.03)</a>
  <a href="#sec04">§ 04 Top edges (F12.04 · T12.E)</a>
  <a href="#sec05">§ 05 Module × Module (F12.05)</a>
  <a href="#sec06">§ 06 Class-mean vs DM1 (F12.06)</a>
  <a href="#sec07">§ 07 DM1-high vs DM1-low Δr (F12.07)</a>
  <a href="#sec08">§ 08 Hub genes (F12.08)</a>
</div>

<section class="atlas-sec" id="sec01">
  <h2>§ 01 — Candidate × candidate co-expression matrix</h2>
  <div class="lead">Pearson r in TCGA-THCA, hierarchically clustered, class-color-coded ticks. Black lines mark module boundaries.</div>
  <div class="body">
    <p>This is the workhorse figure. Reading the matrix:</p>
    <ul>
      <li>Lineage TFs (FOXE1, NKX2-1, PAX8, HHEX) cluster together with strongly negative correlation against most other classes — the lineage-collapse signature is the dominant axis of variation in the candidate panel.</li>
      <li>SFK members (LYN, SRC, FYN) cluster together; YES1 sits separately, consistent with its DGE-DM2-up direction.</li>
      <li>JAK/STAT-upstream cytokine receptors (OSMR, IL6R) cluster with SFK and inflammatory companions, consistent with cytokine-driven SFK activation.</li>
      <li>DDR (PARP1/2, ATR, CHEK1/2) sits in a separate, more loosely coordinated module — consistent with replication-stress addiction being state-conditional rather than co-induced with the inflammatory program.</li>
    </ul>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper12_network/figures/F12_01_corr_heatmap.png" alt="F12.01"/>
    <div class="cap"><b>F12.01 — Co-expression matrix.</b> Pearson r in TCGA-THCA n={sm12['n_thca_samples']}. Class-color ticks. Black lines = module boundaries.</div>
  </div>
</section>

<section class="atlas-sec" id="sec02">
  <h2>§ 02 — Module summary</h2>
  <div class="lead">{sm12['n_modules']} modules from average-linkage clustering on |1 − r| distance.</div>
  <div class="atlas-figure">
    <img src="../results/paper12_network/figures/F12_02_module_summary.png" alt="F12.02"/>
    <div class="cap"><b>F12.02 — Module table figure.</b> Same content as Table T12.M below in figure form.</div>
  </div>
  <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;margin:14px 0 6px">Table T12.M — Module composition (sortable)</h3>
  {T12_mod}
</section>

<section class="atlas-sec" id="sec03">
  <h2>§ 03 — Network graph (MDS layout)</h2>
  <div class="lead">Nodes = candidates; edges = strong correlations (|r|≥0.4); positions from MDS on |1−r| distance.</div>
  <div class="atlas-figure">
    <img src="../results/paper12_network/figures/F12_03_network.png" alt="F12.03"/>
    <div class="cap"><b>F12.03 — Network graph.</b> Red edges = positive correlation; blue = anti. Nodes colored by candidate class.</div>
  </div>
</section>

<section class="atlas-sec" id="sec04">
  <h2>§ 04 — Top-30 edges</h2>
  <div class="lead">|r| ranked edges between candidate genes in TCGA-THCA.</div>
  <div class="atlas-figure-pair">
    <div class="atlas-figure">
      <img src="../results/paper12_network/figures/F12_04_top_edges.png" alt="F12.04"/>
      <div class="cap"><b>F12.04 — Top-30 edges (figure).</b></div>
    </div>
    <div>
      <h3 style="font-family:'JetBrains Mono',monospace;font-size:13px;color:#7B1F2A;text-transform:uppercase;margin:0 0 8px">Table T12.E — Top edges (sortable)</h3>
      {T12_edge}
    </div>
  </div>
</section>

<section class="atlas-sec" id="sec05">
  <h2>§ 05 — Module × Module mean correlation</h2>
  <div class="lead">Diagonal = intra-module mean r; off-diagonal = inter-module. Color: red = positive, blue = negative.</div>
  <div class="atlas-figure">
    <img src="../results/paper12_network/figures/F12_05_module_corr.png" alt="F12.05"/>
    <div class="cap"><b>F12.05 — Module pair correlation.</b> Lineage-TF module typically anti-correlated to inflammatory modules.</div>
  </div>
</section>

<section class="atlas-sec" id="sec06">
  <h2>§ 06 — Class-mean expression vs DM1_like</h2>
  <div class="lead">Per-class mean expression × DM1_like in TCGA-THCA. Pearson r reported per panel.</div>
  <div class="atlas-figure">
    <img src="../results/paper12_network/figures/F12_06_class_vs_dm1.png" alt="F12.06"/>
    <div class="cap"><b>F12.06 — DDR / SFK / JAK-upstream class-mean × DM1_like.</b></div>
  </div>
</section>

<section class="atlas-sec" id="sec07">
  <h2>§ 07 — Δ co-expression: DM1-high vs DM1-low subsets</h2>
  <div class="lead">Difference in correlation matrix between DM1-high and DM1-low halves. Identifies <em>state-conditional</em> co-regulation.</div>
  <div class="body">
    <p>Edges that are strong in DM1-high but absent in DM1-low (green cells in the heatmap) are candidate <strong>state-conditional regulatory wiring</strong> — exactly the kind of co-dependency the synthetic-lethality framework targets. Edges that are strong in both subsets (~0 Δ) reflect housekeeping-level co-regulation and are less interesting for SL.</p>
  </div>
  <div class="atlas-figure">
    <img src="../results/paper12_network/figures/F12_07_corr_diff.png" alt="F12.07"/>
    <div class="cap"><b>F12.07 — Δ Pearson r (DM1-high − DM1-low).</b> Green = stronger in DM1-high; pink = stronger in DM1-low.</div>
  </div>
</section>

<section class="atlas-sec" id="sec08">
  <h2>§ 08 — Hub genes (degree centrality)</h2>
  <div class="lead">Number of |r|≥0.4 edges per node. High-degree candidates are network hubs and prime targets for combinatorial knockdown.</div>
  <div class="atlas-figure">
    <img src="../results/paper12_network/figures/F12_08_hub_genes.png" alt="F12.08"/>
    <div class="cap"><b>F12.08 — Degree centrality.</b> Bar color = candidate class. Top hub: {sm12['top_hub_genes'][-1]}.</div>
  </div>
</section>

<div class="callout">
<strong>Bridge to Paper 9 / 10 priorities.</strong> The hub genes from this network analysis (top 5: {", ".join(sm12['top_hub_genes'])}) should be cross-referenced with the triangulation scoreboard (Paper 10 F11). Hub genes that are also DGE-up + Cox-significant + clinically tractable are the strongest combinatorial-target candidates.
</div>

<footer style="margin-top:48px;padding-top:24px;border-top:2px solid #7B1F2A;font-size:13px;color:#3A4658">
  <p><strong>Reproducibility.</strong> Figures regenerable via <code>python project/scripts/p11_p12_extended_atlas.py</code>; HTML via <code>python project/scripts/build_paper11_12_html.py</code>.</p>
  <p><strong>Hub URL.</strong> <code style="background:rgba(0,0,0,0.06);padding:2px 8px;border-radius:3px">http://40.82.129.113/papers_hub_2026_05_04/paper12_network.html</code></p>
</footer>
</div>
{SORT_JS}
</body></html>"""

(HUB / "paper12_network.html").write_text(HTML12)
print(f"  Wrote paper12_network.html "
      f"({(HUB / 'paper12_network.html').stat().st_size:,} bytes)")
print("Done.")
