#!/usr/bin/env python3
from __future__ import annotations

import html
import shutil
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROP = ROOT / "kthyro_ctc_emt_samsung_proposal"
NOYU = PROP / "outputs/in_silico_paper_no_yu"
PILOT = ROOT / "kthyro_public_pilot/results"
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = "kthyro_no_yu_in_silico"
HUB_ASSET = HUB / "assets" / ASSET
LIVE_ASSET = LIVE / "assets" / ASSET
HTML = HUB / "kthyro_no_yu_in_silico_ctc_emt_dossier.html"
LIVE_HTML = LIVE / "kthyro_no_yu_in_silico_ctc_emt_dossier.html"
LOCAL_HTML = NOYU / "reports/kthyro_no_yu_in_silico_ctc_emt_dossier.html"
ZIP_PATH = NOYU / "kthyro_no_yu_in_silico_paper_packet_2026_05_09.zip"


def esc(x) -> str:
    return html.escape(str(x), quote=True)


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def fmt(x, nd=3):
    try:
        if pd.isna(x):
            return ""
        if abs(float(x)) >= 100:
            return f"{float(x):,.0f}"
        return f"{float(x):.{nd}g}"
    except Exception:
        return str(x)


def df_table(df: pd.DataFrame, cols: list[str] | None = None, max_rows: int = 12) -> str:
    if df.empty:
        return "<p class='muted'>No table available.</p>"
    if cols:
        cols = [c for c in cols if c in df.columns]
        df = df[cols]
    df = df.head(max_rows)
    out = ["<table class='t'><thead><tr>"]
    for c in df.columns:
        out.append(f"<th>{esc(c)}</th>")
    out.append("</tr></thead><tbody>")
    for _, row in df.iterrows():
        out.append("<tr>")
        for c in df.columns:
            v = row[c]
            cls = " class='num'" if isinstance(v, (int, float)) or str(v).replace(".", "", 1).replace("-", "", 1).isdigit() else ""
            out.append(f"<td{cls}>{esc(fmt(v) if cls else v)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def markdown_excerpt(path: Path, max_chars=2400) -> str:
    if not path.exists():
        return ""
    txt = path.read_text(encoding="utf-8")
    if len(txt) > max_chars:
        txt = txt[:max_chars].rstrip() + "\n\n..."
    return f"<pre>{esc(txt)}</pre>"


def copy_assets() -> None:
    HUB_ASSET.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET.mkdir(parents=True, exist_ok=True)
    files: list[Path] = []
    files += sorted((NOYU / "figures").glob("*.png"))
    files += sorted((NOYU / "visual_qc").glob("*.png"))
    files += sorted((NOYU / "source_figures").glob("*.png"))
    files += sorted((NOYU / "tables").glob("*.tsv"))
    files += sorted((NOYU / "reports").glob("*.md"))
    extra = [
        PILOT / "extra_analyses/figures/spatial_niche_adjacency_enrichment_heatmap.png",
        PILOT / "extra_analyses/figures/spatial_condition_and_adjacency_trends.png",
        PILOT / "extra_analyses/figures/spatial_hotspot_colocalization_heatmap.png",
        PILOT / "extra_analyses/figures/spatial_interface_distance_z_by_condition.png",
        PILOT / "extra_analyses/figures/tcga_clinical_stage_axis_associations.png",
        PILOT / "extra_analyses/figures/tcga_axis_covariate_correlation_heatmap.png",
        PILOT / "extra_analyses/figures/external_exact_kthyro_meta_direction_consistency.png",
        PILOT / "extra_analyses/figures/gse151179_rai_direct_label_reanalysis.png",
        PILOT / "figures/proposal/spatial_representative_maps_dark.png",
        PILOT / "figures/proposal/tcga_vulnerability_axes_heatmap_dark.png",
        PILOT / "figures/proposal/tcga_therapeutic_quadrant_dark.png",
    ]
    files += [p for p in extra if p.exists()]
    with zipfile.ZipFile(ZIP_PATH, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in files:
            if p.exists():
                zf.write(p, arcname=p.name)
    files.append(ZIP_PATH)
    seen = set()
    for p in files:
        if not p.exists() or p.name in seen:
            continue
        seen.add(p.name)
        shutil.copy2(p, HUB_ASSET / p.name)
        shutil.copy2(p, LIVE_ASSET / p.name)


def fig_card(img: str, title: str, caption: str = "", cls: str = "") -> str:
    return f"""
    <article class="fig-card {cls}">
      <div class="fig-head"><b>{esc(title)}</b></div>
      <a href="assets/{ASSET}/{esc(img)}"><img src="assets/{ASSET}/{esc(img)}" alt="{esc(title)}" loading="lazy" /></a>
      <p>{esc(caption)}</p>
    </article>
    """


def link_card(label: str, fn: str, note: str = "") -> str:
    return f"<a class='link-card' href='assets/{ASSET}/{esc(fn)}'><b>{esc(label)}</b><span>{esc(note or fn)}</span></a>"


def metric_cards(metrics: list[tuple[str, str, str]]) -> str:
    return "<div class='metric-grid'>" + "".join(
        f"<div class='metric'><b>{esc(v)}</b><span>{esc(k)}</span><p>{esc(note)}</p></div>"
        for k, v, note in metrics
    ) + "</div>"


def build_html() -> str:
    ev = read_tsv(NOYU / "tables/in_silico_evidence_matrix.tsv")
    pub = read_tsv(NOYU / "tables/publishability_decision_table.tsv")
    gap = read_tsv(NOYU / "tables/serial_ctc_ngs_dataset_gap_map.tsv")
    markers = read_tsv(NOYU / "tables/ctc_emt_marker_prior_table.tsv")
    braf = read_tsv(PILOT / "extra_analyses/tables/tcga_braf_only_label_distribution.tsv")
    variance = read_tsv(PILOT / "extra_analyses/tables/tcga_axis_variance_models.tsv")
    clinical = read_tsv(PILOT / "extra_analyses/tables/tcga_clinical_axis_associations_expanded.tsv")
    survival = read_tsv(PILOT / "extra_analyses/tables/tcga_survival_exploratory_axis_associations.tsv")
    spatial = read_tsv(PILOT / "tables/spatial_coherence_statistics.tsv")
    adjacency = read_tsv(PILOT / "extra_analyses/tables/spatial_niche_adjacency_enrichment.tsv")
    hotspot = read_tsv(PILOT / "extra_analyses/tables/spatial_hotspot_colocalization_tests.tsv")
    interface = read_tsv(PILOT / "extra_analyses/tables/spatial_interface_distance_tests.tsv")
    scrna = read_tsv(PILOT / "extra_analyses/tables/scrna_module_celltype_attribution.tsv")
    prot = read_tsv(PILOT / "extra_analyses/tables/bulk_proteomics_kthyro_module_contrasts.tsv")
    prot_trend = read_tsv(PILOT / "extra_analyses/tables/bulk_proteomics_kthyro_dediff_trends.tsv")

    braf_n = int(braf["n"].sum()) if not braf.empty else "NA"
    braf_largest = float(braf["fraction_of_braf"].max() * 100) if not braf.empty else 0
    niche = spatial[spatial["label_set"].eq("niche_label")].dropna(subset=["z"]) if not spatial.empty else pd.DataFrame()
    slides = niche["sample_id"].nunique() if not niche.empty else "NA"
    z_gt2 = int((niche["z"] > 2).sum()) if not niche.empty else "NA"
    z_med = float(niche["z"].median()) if not niche.empty else 0

    driver_only = variance[variance["model"].eq("driver_only")].copy() if not variance.empty else pd.DataFrame()
    if not driver_only.empty:
        driver_only["r_squared_pct"] = driver_only["r_squared"] * 100
        driver_only["unexplained_pct"] = driver_only["unexplained_fraction"] * 100
    prot_atc = prot[prot["contrast"].eq("ATC_vs_PTC")].copy() if not prot.empty else pd.DataFrame()
    if not prot_atc.empty:
        prot_atc = prot_atc.sort_values("cohens_d_group_a_minus_b", key=lambda s: s.abs(), ascending=False)
    if not adjacency.empty and "median_log2_enrichment_proxy" in adjacency.columns:
        adjacency = adjacency.sort_values("median_log2_enrichment_proxy", ascending=False)
    elif not adjacency.empty and "log2_enrichment_proxy" in adjacency.columns:
        adjacency = adjacency.sort_values("log2_enrichment_proxy", ascending=False)
    if not clinical.empty and "fdr_q" in clinical.columns:
        clinical = clinical.sort_values("fdr_q")

    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width,initial-scale=1.0" />
<title>No-Yu-data In Silico CTC-EMT Dossier</title>
<style>
:root{{
  --bg:#080d18; --panel:#101a2c; --panel2:#16213a; --line:#263954;
  --text:#e9f0fb; --muted:#9aa9bd; --gold:#ffd28a; --cyan:#5bdcff;
  --teal:#35d39d; --purple:#b38cff; --orange:#ff9b66; --warn:#ff8066;
}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--text);font-family:Inter,"Noto Sans KR",-apple-system,sans-serif;font-size:14px;line-height:1.55}}
a{{color:var(--teal);text-decoration:none}} a:hover{{color:var(--gold);text-decoration:underline}}
h1,h2,h3{{font-family:"Cormorant Garamond","Noto Serif KR",serif;color:#fff8e7;letter-spacing:-.4px}}
code,pre,.mono{{font-family:"JetBrains Mono","SF Mono",Menlo,monospace}}
.hero{{padding:68px 34px 48px;border-bottom:1px solid var(--line);background:radial-gradient(ellipse at 18% 10%,rgba(91,220,255,.18),transparent 34%),radial-gradient(ellipse at 82% 8%,rgba(179,140,255,.14),transparent 30%),linear-gradient(135deg,#0b1020,#10182b 58%,#050811)}}
.hero-inner{{max-width:1380px;margin:0 auto}}
.kicker{{font:800 10px/1.4 "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.22em;text-transform:uppercase;margin-bottom:16px}}
h1{{font-size:64px;line-height:1;margin:0 0 16px}} h1 em{{color:var(--gold);font-style:italic}}
.lead{{font-size:18px;max-width:1050px;color:#cfdae8;font-family:"Noto Serif KR",serif;line-height:1.72}}
.actions{{display:grid;grid-template-columns:1.2fr 1fr 1fr;gap:14px;max-width:1060px;margin:24px 0}}
.actions a{{display:block;border:1px solid rgba(255,210,138,.38);border-radius:14px;padding:18px 20px;background:linear-gradient(135deg,rgba(91,220,255,.16),rgba(255,210,138,.08));text-decoration:none}}
.actions a b{{display:block;color:#fff8e7;font-size:23px;font-family:"Cormorant Garamond","Noto Serif KR",serif;margin-bottom:4px}}
.actions a span{{font:10px/1.45 "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.06em;text-transform:uppercase}}
.urlbox{{max-width:1060px;border:1px solid rgba(91,220,255,.3);border-left:4px solid var(--cyan);background:#091321;border-radius:0 10px 10px 0;padding:10px 13px;margin:10px 0;font:12px/1.45 "JetBrains Mono",monospace;color:#cfe9ff;overflow-wrap:anywhere}}
.wrap{{max-width:1380px;margin:0 auto;display:grid;grid-template-columns:250px 1fr;gap:38px;padding:0 32px 78px}}
.toc{{position:sticky;top:0;align-self:start;max-height:100vh;overflow:auto;padding:28px 0;border-right:1px solid var(--line)}}
.toc h4{{font:800 10px/1.4 "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;text-transform:uppercase;margin:0 0 12px}}
.toc ol{{list-style:none;counter-reset:t;margin:0;padding:0}} .toc li{{counter-increment:t;margin:5px 0;font-size:12px}}
.toc li:before{{content:counter(t,decimal-leading-zero) "  ";font-family:"JetBrains Mono",monospace;color:var(--gold);font-size:10px}}
.toc a{{color:#d5dfec}}
main{{padding-top:24px;min-width:0}} section{{padding:30px 0;border-bottom:1px solid var(--line)}}
h2{{font-size:36px;margin:0 0 6px}} h2 .num{{font:400 14px/1 "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:12px}}
.sub{{font:11px/1.45 "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.09em;text-transform:uppercase;margin-bottom:18px}}
.metric-grid{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin:20px 0}}
.metric{{background:rgba(255,255,255,.045);border:1px solid rgba(255,210,138,.18);border-radius:12px;padding:13px}}
.metric b{{display:block;color:var(--gold);font:800 28px/1 "Cormorant Garamond",serif;margin-bottom:7px}}
.metric span{{display:block;color:var(--muted);font:800 9px/1.35 "JetBrains Mono",monospace;letter-spacing:.08em;text-transform:uppercase}}
.metric p{{margin:6px 0 0;color:#d5dfec;font-size:12px}}
.fig-grid{{display:grid;grid-template-columns:repeat(2,1fr);gap:18px;margin-top:16px}}
.fig-grid.wide{{grid-template-columns:1fr}}
.fig-card{{background:#0e1728;border:1px solid var(--line);border-radius:13px;overflow:hidden;box-shadow:0 14px 38px rgba(0,0,0,.24)}}
.fig-card.feature{{border-color:rgba(91,220,255,.45)}}
.fig-head{{padding:12px 14px;background:#121f35;border-bottom:1px solid var(--line)}} .fig-head b{{color:#fff8e7}}
.fig-card img{{display:block;width:100%;height:auto;background:#050811;border-bottom:1px solid var(--line)}}
.fig-card p{{margin:0;padding:12px 14px 14px;color:#d5dfec;font-size:13px}}
.grid{{display:grid;gap:14px}} .grid-2{{grid-template-columns:1fr 1fr}} .grid-3{{grid-template-columns:repeat(3,1fr)}}
.box{{background:#0e1728;border:1px solid var(--line);border-left:3px solid var(--gold);border-radius:0 9px 9px 0;padding:14px 17px;margin:14px 0}}
.box.warn{{border-left-color:var(--warn);background:#1d1410}} .box.good{{border-left-color:var(--teal);background:#0e1d18}}
.box h3{{font:800 12px/1.4 "JetBrains Mono",monospace;letter-spacing:.1em;text-transform:uppercase;color:var(--gold);margin:0 0 8px}}
table.t{{border-collapse:collapse;width:100%;font-size:12.5px;margin:10px 0 16px}} .t th,.t td{{border:1px solid var(--line);padding:8px 9px;text-align:left;vertical-align:top}}
.t th{{background:#16213a;color:var(--gold);font:800 10.5px/1.3 "JetBrains Mono",monospace;letter-spacing:.04em;text-transform:uppercase}}
.t tr:nth-child(even) td{{background:rgba(255,255,255,.025)}} .t td.num{{font-family:"JetBrains Mono",monospace;text-align:right;color:var(--gold)}}
.links{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}} .link-card{{display:block;background:#101a2c;border:1px solid var(--line);border-radius:8px;padding:12px;min-height:78px;text-decoration:none}}
.link-card b{{display:block;color:var(--gold);font:800 10px/1.4 "JetBrains Mono",monospace;letter-spacing:.08em;text-transform:uppercase;margin-bottom:4px}}
.link-card span{{color:#d5dfec;font-size:12px;overflow-wrap:anywhere}}
pre{{white-space:pre-wrap;background:#0a111d;border:1px solid var(--line);border-radius:8px;padding:14px;color:#dce7f4;font-size:12px;overflow:auto;max-height:520px}}
.muted{{color:var(--muted)}} footer{{padding:26px 32px;border-top:1px solid var(--line);text-align:center;color:var(--muted);font:11px/1.4 "JetBrains Mono",monospace}}
@media(max-width:1100px){{.wrap{{grid-template-columns:1fr}}.toc{{display:none}}.fig-grid,.grid-2,.grid-3,.actions,.metric-grid,.links{{grid-template-columns:1fr}} h1{{font-size:42px}}}}
</style>
</head>
<body>
<header class="hero"><div class="hero-inner">
  <div class="kicker">K-Thyro · No-Yu-data in silico dossier · 2026-05-09</div>
  <h1>CTC-EMT paper, <em>without raw CTC data?</em></h1>
  <p class="lead">가능하다. 하지만 정직한 이름은 discovery paper가 아니라 <b>public-data tissue-state / spatial / single-cell / proteomics 기반 CTC-EMT-NGS prior map</b>이다. 이 웹 dossier는 유형원 교수님 patient-level serial CTC data 없이 지금 쓸 수 있는 모든 in silico 근거와, 절대 넘으면 안 되는 claim boundary를 그림 중심으로 정리한다.</p>
  <div class="actions">
    <a href="assets/{ASSET}/no_yu_in_silico_visual_storyboard_cv2.png"><b>Visual Storyboard</b><span>CV2 contact sheet · all in silico figure panels</span></a>
    <a href="assets/{ASSET}/kthyro_no_yu_in_silico_paper_packet_2026_05_09.zip"><b>Download Packet</b><span>figures + reports + tables + source figures</span></a>
    <a href="#publish"><b>Verdict</b><span>what can be published now</span></a>
  </div>
  <div class="urlbox">전달용 웹: http://40.82.129.113:8012/papers_hub_2026_05_04/kthyro_no_yu_in_silico_ctc_emt_dossier.html</div>
  <div class="urlbox">패키지 ZIP: http://40.82.129.113:8012/papers_hub_2026_05_04/assets/{ASSET}/kthyro_no_yu_in_silico_paper_packet_2026_05_09.zip</div>
  {metric_cards([
    ("verdict", "YES", "framework/prior-map paper now"),
    ("BRAF-mutant", str(braf_n), "TCGA tumors split across states"),
    ("largest BRAF state", f"{braf_largest:.1f}%", "mutation alone is insufficient"),
    ("spatial slides", f"{z_gt2}/{slides}", "same-niche coherence z>2"),
    ("spatial spots", "57,144", "GSE250521 tissue organization"),
    ("core boundary", "NO", "no new CTC transition claim"),
  ])}
</div></header>

<div class="wrap">
<nav class="toc"><h4>Contents</h4><ol>
  <li><a href="#visual">Visual first</a></li>
  <li><a href="#publish">Publishability</a></li>
  <li><a href="#evidence">Evidence matrix</a></li>
  <li><a href="#tcga">TCGA state genetics</a></li>
  <li><a href="#spatial">Spatial organization</a></li>
  <li><a href="#scrna">Single-cell context</a></li>
  <li><a href="#protein">Proteomics</a></li>
  <li><a href="#panel">CTC-EMT-NGS panel</a></li>
  <li><a href="#gap">Dataset gap</a></li>
  <li><a href="#paper">Paper skeleton</a></li>
  <li><a href="#paths">Paths</a></li>
</ol></nav>
<main>

<section id="visual">
  <h2><span class="num">01</span>Visual First</h2>
  <p class="sub">CV2 contact sheets and figure-first overview</p>
  <div class="fig-grid wide">
    {fig_card("no_yu_in_silico_visual_storyboard_cv2.png", "CV2 storyboard: all no-Yu-data in silico panels", "Use this first for visual QC and manuscript flow.", "feature")}
    {fig_card("no_yu_decision_board_cv2.png", "CV2 decision board", "Keep this in front of the team: what is publishable now, what is not.", "feature")}
  </div>
  <div class="fig-grid">
    {fig_card("F01_no_yu_claim_ladder.png", "Claim ladder", "Defines the allowed claim envelope without patient-level CTC data.")}
    {fig_card("F02_public_data_layer_map.png", "Public data layer map", "TCGA, scRNA, spatial, proteomics, Yu literature, and future hospital cohort.")}
    {fig_card("F09_publishability_decision.png", "Publishability decision", "Framework paper yes; transition discovery and residual-risk prediction no.")}
    {fig_card("F08_ctc_emt_ngs_prior_panel.png", "CTC-EMT-NGS prior panel", "Final marker/genetic map to test later with serial blood.")}
  </div>
</section>

<section id="publish">
  <h2><span class="num">02</span>Publishability</h2>
  <p class="sub">Ruthless decision: what can be written today</p>
  <div class="grid grid-2">
    <div class="box good"><h3>Go now</h3><p><b>Public-data framework/prior-map paper.</b> The manuscript can argue that public tissue, spatial, single-cell, and proteomic layers define a conservative CTC-EMT-NGS interpretation framework for a future thyroidectomy perturbation cohort.</p></div>
    <div class="box warn"><h3>Do not claim</h3><p>No new postoperative CTC-EMT transition, no recurrence prediction, no clinical diagnostic, and no proof that public tissue states are CTC origin.</p></div>
  </div>
  {df_table(pub, max_rows=10)}
</section>

<section id="evidence">
  <h2><span class="num">03</span>Evidence Matrix</h2>
  <p class="sub">All in silico layers with allowed and forbidden claims</p>
  {df_table(ev, max_rows=12)}
</section>

<section id="tcga">
  <h2><span class="num">04</span>TCGA State Genetics</h2>
  <p class="sub">Mutation alone is not enough; this supports genetic-state mapping, not CTC biology</p>
  <div class="fig-grid">
    {fig_card("F03_tcga_braf_state_split.png", "BRAF-mutant PTC state split", "BRAF-mutant tumors distribute across multiple tissue-state labels.")}
    {fig_card("F04_tcga_driver_variance_boundary.png", "Driver-only variance boundary", "Driver categories explain only part of state-axis variance.")}
    {fig_card("tcga_braf_only_vulnerability_axis_heatmap.png", "Source: BRAF-only vulnerability heatmap", "Existing public-pilot figure.")}
    {fig_card("tcga_clinical_stage_axis_associations.png", "Clinical-risk exploratory associations", "Exploratory only; not outcome prediction.")}
  </div>
  <h3>BRAF label distribution</h3>{df_table(braf, max_rows=10)}
  <h3>Driver-only variance model</h3>{df_table(driver_only, ["axis_label","model","r_squared_pct","unexplained_pct","claim_boundary"], max_rows=12)}
  <h3>Clinical axis associations, strongest rows</h3>{df_table(clinical, max_rows=10)}
  <h3>Survival exploratory table</h3>{df_table(survival, max_rows=10)}
</section>

<section id="spatial">
  <h2><span class="num">05</span>Spatial Organization</h2>
  <p class="sub">Spatial data nominate ROI and tissue context; they do not prove CTC shedding</p>
  <div class="fig-grid">
    {fig_card("F05_spatial_tissue_state_organization.png", "Same-niche coherence", f"{z_gt2}/{slides} slides have z>2; median z={z_med:.2f}.")}
    {fig_card("spatial_representative_maps_dark.png", "Representative spatial maps", "Spatial tissue-state visual layer from public pilot.")}
    {fig_card("spatial_niche_adjacency_enrichment_heatmap.png", "Niche adjacency enrichment", "Spot-neighborhood proxy; ROI nomination only.")}
    {fig_card("spatial_condition_and_adjacency_trends.png", "Condition and adjacency trends", "Exploratory slide-level trend; n=4 per condition.")}
    {fig_card("spatial_hotspot_colocalization_heatmap.png", "Hotspot colocalization", "Co-localization support for tissue-state organization.")}
    {fig_card("spatial_interface_distance_z_by_condition.png", "Interface distance", "Spatial interface proxy; not biological replication.")}
  </div>
  <h3>Spatial coherence table</h3>{df_table(niche, ["sample_id","condition","label_set","observed","null_mean","z","empirical_p"], max_rows=18)}
  <h3>Adjacency enrichment</h3>{df_table(adjacency, max_rows=14)}
  <h3>Hotspot colocalization</h3>{df_table(hotspot, max_rows=12)}
  <h3>Interface distance tests</h3>{df_table(interface, max_rows=12)}
</section>

<section id="scrna">
  <h2><span class="num">06</span>Single-Cell Context</h2>
  <p class="sub">Marker-cell context for panel rationale; not CTC data</p>
  <div class="fig-grid">
    {fig_card("F06_scrna_marker_context.png", "scRNA marker context", "Top cell type attribution for marker modules.")}
    {fig_card("scrna_module_celltype_attribution_heatmap.png", "Source: scRNA attribution heatmap", "Existing public-pilot figure.")}
  </div>
  {df_table(scrna, max_rows=14)}
</section>

<section id="protein">
  <h2><span class="num">07</span>Proteomics</h2>
  <p class="sub">Protein directionality supports tissue-state priors; not spatial protein localization</p>
  <div class="fig-grid">
    {fig_card("F07_proteomics_dediff_direction.png", "ATC vs PTC protein directionality", "RAI loss and myeloid/TGFB gain are the strongest safe messages.")}
    {fig_card("bulk_proteomics_kthyro_contrast_heatmap.png", "Source: proteomics contrast heatmap", "Bulk proteomics proxy; not spatial.")}
    {fig_card("bulk_proteomics_kthyro_axis_boxplots.png", "Source: proteomics axis boxplots", "PTC/PDTC/ATC protein module scores.")}
  </div>
  <h3>ATC vs PTC contrasts</h3>{df_table(prot_atc, ["module","n_group_a","n_group_b","cohens_d_group_a_minus_b","mannwhitney_p","fdr_q","claim_boundary"], max_rows=12)}
  <h3>Dedifferentiation trends</h3>{df_table(prot_trend, max_rows=12)}
</section>

<section id="panel">
  <h2><span class="num">08</span>CTC-EMT-NGS Prior Panel</h2>
  <p class="sub">Use as marker rationale now; test with serial blood later</p>
  <div class="fig-grid">
    {fig_card("F08_ctc_emt_ngs_prior_panel.png", "Final marker/genetic prior panel", "Tissue-informed marker and clone-anchor modules.")}
  </div>
  {df_table(markers, max_rows=10)}
</section>

<section id="gap">
  <h2><span class="num">09</span>Dataset Gap</h2>
  <p class="sub">Why the prospective hospital cohort remains essential</p>
  {df_table(gap, max_rows=12)}
  <div class="box warn"><h3>Critical boundary</h3><p>Without Prof. Yu's patient-level serial CTC table, we cannot test postoperative E/EM/M transition. The in silico paper should explicitly state this limitation.</p></div>
</section>

<section id="paper">
  <h2><span class="num">10</span>Paper Skeleton</h2>
  <p class="sub">Manuscript-ready outline and claim defense</p>
  <div class="grid grid-2">
    <div>{markdown_excerpt(NOYU / "reports/MANUSCRIPT_SKELETON_NO_YU_DATA_KR.md")}</div>
    <div>{markdown_excerpt(NOYU / "reports/CLAIM_BOUNDARY_NO_YU_DATA.md")}</div>
  </div>
</section>

<section id="paths">
  <h2><span class="num">11</span>Paths</h2>
  <p class="sub">All generated artifacts and source files</p>
  <div class="links">
    {link_card("ZIP packet", "kthyro_no_yu_in_silico_paper_packet_2026_05_09.zip", "all figures/tables/reports")}
    {link_card("Storyboard", "no_yu_in_silico_visual_storyboard_cv2.png", "CV2 figure contact sheet")}
    {link_card("Decision board", "no_yu_decision_board_cv2.png", "publishability board")}
    {link_card("Feasibility report", "NO_YU_DATA_IN_SILICO_FEASIBILITY_KR.md", "claim-safe Korean report")}
    {link_card("Manuscript skeleton", "MANUSCRIPT_SKELETON_NO_YU_DATA_KR.md", "paper outline")}
    {link_card("Evidence table", "in_silico_evidence_matrix.tsv", "allowed/forbidden claims")}
    {link_card("Marker table", "ctc_emt_marker_prior_table.tsv", "CTC-EMT-NGS prior")}
    {link_card("Gap map", "serial_ctc_ngs_dataset_gap_map.tsv", "missing dataset map")}
  </div>
  <table class="t"><tbody>
    <tr><th>Local output dir</th><td>{esc(str(NOYU))}</td></tr>
    <tr><th>HTML source</th><td>{esc(str(HTML))}</td></tr>
    <tr><th>Live HTML</th><td>{esc(str(LIVE_HTML))}</td></tr>
    <tr><th>Live assets</th><td>{esc(str(LIVE_ASSET))}</td></tr>
  </tbody></table>
</section>

</main></div>
<footer>No-Yu-data In Silico CTC-EMT Dossier · framework/prior-map paper only · no new CTC transition claim</footer>
</body></html>"""


def main() -> None:
    copy_assets()
    text = build_html()
    HTML.write_text(text, encoding="utf-8")
    LOCAL_HTML.write_text(text, encoding="utf-8")
    LIVE.mkdir(parents=True, exist_ok=True)
    shutil.copy2(HTML, LIVE_HTML)
    print(f"Wrote {HTML}")
    print(f"Wrote {LOCAL_HTML}")
    print(f"Deployed {LIVE_HTML}")
    print(f"Assets {LIVE_ASSET}")


if __name__ == "__main__":
    main()
