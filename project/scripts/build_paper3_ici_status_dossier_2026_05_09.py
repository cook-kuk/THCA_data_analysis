#!/usr/bin/env python3
"""Build the 2026-05-09 Paper 3 / ICI status dossier.

This is a synthesis artifact only. It reads already-existing result tables and
copies already-rendered figures into the papers hub. It does not run a new
analysis.
"""

from __future__ import annotations

import csv
import html
import json
import math
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HUB = ROOT / "project" / "papers_hub_2026_05_04"
LIVE = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_DIR = HUB / "assets" / "paper3_ici_2026_05_09"
LIVE_ASSET_DIR = LIVE / "assets" / "paper3_ici_2026_05_09"
REPORT = ROOT / "project" / "reports" / "paper3_ici" / "ICI_STATUS_DOSSIER_2026_05_09.md"
HTML = HUB / "paper3_ici_status_dossier_2026_05_09.html"
LIVE_HTML = LIVE / "paper3_ici_status_dossier_2026_05_09.html"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def fnum(value: object, digits: int = 2) -> str:
    if value is None or value == "":
        return ""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(x):
        return ""
    if abs(x) >= 100:
        return f"{x:.0f}"
    if abs(x) >= 10:
        return f"{x:.1f}"
    return f"{x:.{digits}f}"


def fp(value: object) -> str:
    if value is None or value == "":
        return ""
    try:
        x = float(value)
    except (TypeError, ValueError):
        return str(value)
    if math.isnan(x):
        return ""
    if x < 1e-3:
        return f"{x:.1e}"
    return f"{x:.3f}"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def esc(value: object) -> str:
    return html.escape(str(value))


def table(headers: list[str], rows: list[list[object]], cls: str = "") -> str:
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{esc(c)}</td>" for c in row) + "</tr>")
    return f'<table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table>'


def copy_assets() -> dict[str, str]:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)

    figures = {
        "atc_meta": ROOT / "project/results/paper3_ici_track_b_lite/figures_png/fig1_atc_vs_other_meta.png",
        "pca": ROOT / "project/results/paper3_ici_track_b_lite/figures_png/fig2_pca_scatter_by_cohort.png",
        "dial_lite": ROOT / "project/results/paper3_ici_track_b_lite/figures_png/fig4_dial_lite_forest.png",
        "rai_axis": ROOT / "project/results/paper3_ici_track_b_lite/figures_png/fig7_rai_after_vs_before.png",
        "phasec_forest": ROOT / "project/results/paper11_pancancer/phase_C_ICI/results/figures/forest_logOR_response.png",
        "phasec_scores": ROOT / "project/results/paper11_pancancer/phase_C_ICI/results/figures/score_by_response_each_cohort.png",
        "hla1_forest": ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla/figures/F1_forest_HLA1_response.png",
        "hla2_forest": ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla/figures/F2_forest_HLA2_response.png",
        "hla_heatmap": ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla/figures/F4_OR_heatmap_per_cohort.png",
        "hla_survival": ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla/figures/F7_survival_forest.png",
        "p3_public_heatmap": ROOT / "project/results/p3_p9_full_execution/paper3/fig_p3_immune_heatmap.png",
    }

    urls: dict[str, str] = {}
    for key, src in figures.items():
        if not src.exists():
            continue
        dst_name = f"{key}_{src.name}"
        dst = ASSET_DIR / dst_name
        live_dst = LIVE_ASSET_DIR / dst_name
        shutil.copy2(src, dst)
        shutil.copy2(src, live_dst)
        urls[key] = f"assets/paper3_ici_2026_05_09/{dst_name}"
    return urls


def build() -> None:
    urls = copy_assets()

    track_b_dir = ROOT / "project/results/paper3_ici_track_b_lite"
    phase_c_dir = ROOT / "project/results/paper11_pancancer/phase_C_ICI/results/tables"
    hla18_dir = ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla/tables"
    p2b_dir = ROOT / "project/results/p2_braf_nature_sprint_2026_05_09"

    atc_rows = read_tsv(track_b_dir / "atc_vs_other_meta_pooled.tsv")
    rai_rows = read_tsv(track_b_dir / "gse151179_rai_after_vs_before.tsv")
    phase4 = read_json(track_b_dir / "phase4_summary.json")
    phase_c = read_tsv(phase_c_dir / "final_phase_C_ICI_summary_table.tsv")
    phase_c_surv = read_tsv(phase_c_dir / "survival_cox_per_cohort.tsv")
    hla_meta = read_tsv(hla18_dir / "T3_meta_OR_response.tsv")
    hla_strata = read_tsv(hla18_dir / "T5_within_DM1_strata_HLA1_response.tsv")
    hla_extended = read_tsv(hla18_dir / "T7_composite_DM1_inflam_extended.tsv")
    r5_adapt = read_json(p2b_dir / "r5_ici_cohorts" / "r5_adaptive_meta.json")
    r5_head = read_json(p2b_dir / "r5_ici_cohorts" / "r5_headline.json")
    r10_head = read_json(p2b_dir / "r10_ici_regimen" / "r10_headline.json")
    h27_head = read_json(p2b_dir / "h27_ici_panels" / "h27_headline.json")
    track18 = read_json(ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla/track18_summary.json")

    atc_table_rows = [
        [r["module"], fnum(r["d_pooled_n_weighted"], 2), f'{r["sign_consistency"]}/{r["of_total"]}']
        for r in sorted(atc_rows, key=lambda x: abs(float(x["d_pooled_n_weighted"])), reverse=True)
    ]
    rai_table_rows = [
        [r["module"], fnum(r["cohens_d_after_minus_before"], 2), fp(r["p"]), r["sign"]]
        for r in sorted(rai_rows, key=lambda x: abs(float(x["cohens_d_after_minus_before"])), reverse=True)
    ]
    phase_c_rows = [
        [
            r["score"],
            r["go_no_go"],
            fnum(r["pooled_OR_per_z"], 2),
            fp(r["pooled_p"]),
            f'{r["n_cohorts_d_positive"]}/{r["k_cohorts_with_response"]}',
            fnum(r["median_HR_OS"], 2),
        ]
        for r in phase_c
    ]
    hla_meta_rows = [
        [r["score"], fnum(r["OR"], 2), f'{fnum(r["lo"], 2)}-{fnum(r["hi"], 2)}', fp(r["p"]), f'{r["n_pos"]}/{r["n_total"]}']
        for r in hla_meta
    ]
    hla_strata_rows = [
        [r["cohort"], r["dm1_strata"], r["n"], r["n_R"], fnum(r["OR_hla1"], 2), fp(r["p_hla1"])]
        for r in hla_strata
        if r.get("OR_hla1")
    ]
    hla_extended_rows = [
        [r["score"], fnum(r["pooled_OR"], 2), f'{fnum(r["lo"], 2)}-{fnum(r["hi"], 2)}', fp(r["p"]), f'{r["sign_positive"]}/{r["k_cohorts"]}']
        for r in hla_extended
    ]
    surv_rows = [
        [r["cohort"], r["score"], fnum(r["HR_per_z"], 2), f'{fnum(r["HR_lo"], 2)}-{fnum(r["HR_hi"], 2)}', fp(r["p"])]
        for r in phase_c_surv
        if r["score"] in {"DM1_inflam_composite", "lineage_portable_DM1", "IFNG_T_cell_inflamed", "HLA_class_I", "checkpoint_exhaustion"}
    ]

    srcs = {
        "Track B-lite report": ROOT / "project/reports/paper3_ici/paper3_ici_track_B_lite_analysis_report.md",
        "Track A frozen bundle": ROOT / "project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md",
        "ATC meta TSV": track_b_dir / "atc_vs_other_meta_pooled.tsv",
        "RAI axis TSV": track_b_dir / "gse151179_rai_after_vs_before.tsv",
        "Phase C summary": phase_c_dir / "final_phase_C_ICI_summary_table.tsv",
        "Track18 HLA report": ROOT / "project/results/hla_deepdive_2026_05_08/track18_ici_hla/track18_report.md",
        "R5 report": p2b_dir / "r5_ici_cohorts/R5_REPORT.md",
        "R10 report": p2b_dir / "r10_ici_regimen/R10_REPORT.md",
        "H27 report": p2b_dir / "h27_ici_panels/H27_REPORT.md",
    }

    md = f"""# ICI status dossier — 2026-05-09

## Verdict

Paper 3 should stay framed as **ICI vulnerability / ICI-readiness / immunogenomic prioritization**, not thyroid ICI response prediction. The best current evidence is:

- Thyroid-side biology: ATC is inflamed-but-myeloid-suppressed and dedifferentiated (`myeloid_suppressive d={fnum(next(r for r in atc_rows if r['module']=='myeloid_suppressive')['d_pooled_n_weighted'], 2)}`, `thyroid_differentiation d={fnum(next(r for r in atc_rows if r['module']=='thyroid_differentiation')['d_pooled_n_weighted'], 2)}`).
- RAI-side bridge: GSE151179 post-RAI vs pre-RAI shows thyroid-differentiation loss (`d={fnum(next(r for r in rai_rows if r['module']=='thyroid_differentiation')['cohens_d_after_minus_before'], 2)}`, `p={fp(next(r for r in rai_rows if r['module']=='thyroid_differentiation')['p'])}`) with HLA-II / myeloid upward trends.
- External ICI coherence: four pre-treatment ICI cohorts, `n=421`, show IFNG, HLA-I, checkpoint, and cytolytic modules associated with response.
- Load-bearing response component: HLA-I expression module, not HLA-II, carries the strongest HLA class signal (`HLA-I OR={fnum(track18['phase_C_v2_reconfirm']['HLA_I']['OR'], 2)}, p={fp(track18['phase_C_v2_reconfirm']['HLA_I']['p'])}`; `HLA-II OR={fnum(track18['new_HLA_II']['OR_meta'], 2)}, p={fp(track18['new_HLA_II']['p_meta'])}`).
- BRAF-cPTC candidate layer: DM1 BRAF-cPTC has TIS-Ayers `d={fnum(h27_head['TIS_Ayers']['d_DM1_minus_DM2'], 2)}`, FDR `={fp(h27_head['TIS_Ayers']['BH_FDR'])}`; 44/110 are ICI-likely by DM1 + TIS + CD8 + HLA-I.

## Claim boundary

- Allowed: ICI vulnerability, ICI-readiness, immunogenomic prioritization, response-biology coherence in non-thyroid ICI cohorts.
- Forbidden: thyroid ICI response predictor, treatment recommendation, clinical-grade selection, causal immunotherapy mechanism.

## Current decision

Track A remains frozen. True Track B requires Paper 1 bioRxiv, Paper 2 closure, and explicit `Paper 3 Track B 시작`.

## Source paths

""" + "\n".join(f"- {name}: `{rel(path)}`" for name, path in srcs.items()) + "\n"

    REPORT.write_text(md)

    stat_cards = [
        ("ATC myeloid", f"d={fnum(next(r for r in atc_rows if r['module']=='myeloid_suppressive')['d_pooled_n_weighted'], 2)}", "4/4 sign; Track B-lite"),
        ("ATC thyroid diff", f"d={fnum(next(r for r in atc_rows if r['module']=='thyroid_differentiation')['d_pooled_n_weighted'], 2)}", "lineage collapse"),
        ("RAI post-vs-pre", f"d={fnum(next(r for r in rai_rows if r['module']=='thyroid_differentiation')['cohens_d_after_minus_before'], 2)}", f"p={fp(next(r for r in rai_rows if r['module']=='thyroid_differentiation')['p'])}"),
        ("ICI cohorts", "n=421", "4 pre-treatment cohorts"),
        ("IFNG response", f"OR={fnum(track18['phase_C_v2_reconfirm']['IFNG']['OR'], 2)}", f"p={fp(track18['phase_C_v2_reconfirm']['IFNG']['p'])}"),
        ("HLA-I response", f"OR={fnum(track18['phase_C_v2_reconfirm']['HLA_I']['OR'], 2)}", f"p={fp(track18['phase_C_v2_reconfirm']['HLA_I']['p'])}"),
        ("HLA-II response", f"OR={fnum(track18['new_HLA_II']['OR_meta'], 2)}", f"p={fp(track18['new_HLA_II']['p_meta'])}"),
        ("DM1-adaptive", f"OR={fnum(r5_adapt['mh_OR'], 2)}", f"p={fp(r5_adapt['p'])}; 5/5 sign"),
    ]
    cards_html = "\n".join(
        f'<div class="stat"><b>{esc(k)}</b><strong>{esc(v)}</strong><span>{esc(sub)}</span></div>'
        for k, v, sub in stat_cards
    )

    gallery_items = [
        ("Track B-lite ATC module meta", "atc_meta", "ATC shows myeloid/checkpoint/HLA activation plus thyroid-differentiation loss."),
        ("DIAL-lite caution", "dial_lite", "Hugo + Riaz minimal audit flips 4/7 modules; this protects the claim boundary."),
        ("RAI bridge", "rai_axis", "GSE151179 post-RAI samples show the strongest shift in thyroid-differentiation loss."),
        ("Phase C v2 response forest", "phasec_forest", "IFNG, HLA-I, checkpoint, cytolytic clear pooled response association."),
        ("HLA-I forest", "hla1_forest", "HLA-I expression module is response-associated in 4/4 cohorts."),
        ("HLA-II forest", "hla2_forest", "HLA-II is null at meta level, separating class-I from class-II biology."),
        ("HLA response heatmap", "hla_heatmap", "Per-cohort OR pattern shows the class-I load-bearing axis."),
        ("ICI survival forest", "hla_survival", "HLA-I and IFNG/checkpoint are OS-protective in the two OS-bearing cohorts."),
    ]
    gallery_html = "\n".join(
        f'''<figure><img src="{urls[key]}" alt="{esc(title)}"><figcaption><b>{esc(title)}</b><br>{esc(caption)}</figcaption></figure>'''
        for title, key, caption in gallery_items
        if key in urls
    )

    source_rows = [[name, rel(path)] for name, path in srcs.items()]

    html_doc = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Paper 3 ICI Status Dossier — 2026-05-09</title>
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=JetBrains+Mono:wght@400;500;600;700&family=Noto+Sans+KR:wght@400;500;700;800&display=swap" rel="stylesheet" />
<style>
:root {{
  --bg:#101516; --panel:#172022; --panel2:#1d282b; --ink:#edf1ed; --muted:#aab6ae;
  --line:#354245; --gold:#d7b46a; --green:#7bc88f; --red:#e1776f; --blue:#7fb4d8;
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--ink); font-family:'Noto Sans KR',system-ui,sans-serif; line-height:1.62; }}
a {{ color:var(--gold); text-decoration:none; }}
.wrap {{ max-width:1480px; margin:0 auto; padding:28px; }}
.hero {{ min-height:78vh; display:grid; align-content:center; padding:52px 0 36px; border-bottom:1px solid var(--line); }}
.kicker {{ font-family:'JetBrains Mono',monospace; color:var(--gold); font-size:12px; letter-spacing:.16em; text-transform:uppercase; }}
h1 {{ font-family:'Cormorant Garamond',serif; font-size:clamp(54px,8vw,118px); line-height:.92; margin:16px 0 18px; max-width:1050px; }}
.lead {{ max-width:980px; color:#d8ded9; font-size:20px; }}
.crumbs {{ display:flex; gap:10px; flex-wrap:wrap; margin-top:26px; font-family:'JetBrains Mono',monospace; font-size:12px; color:var(--muted); }}
.crumbs span {{ border:1px solid var(--line); padding:7px 10px; border-radius:4px; background:rgba(255,255,255,.03); }}
.stats {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin-top:30px; }}
.stat {{ border:1px solid var(--line); background:linear-gradient(180deg,rgba(255,255,255,.04),rgba(255,255,255,.015)); padding:15px 16px; border-radius:6px; min-height:112px; }}
.stat b {{ display:block; color:var(--muted); font-family:'JetBrains Mono',monospace; font-size:12px; }}
.stat strong {{ display:block; color:var(--gold); font-size:30px; margin:8px 0 3px; font-family:'Cormorant Garamond',serif; }}
.stat span {{ color:#c5cec7; font-size:13px; }}
.layout {{ display:grid; grid-template-columns:260px minmax(0,1fr); gap:28px; align-items:start; }}
.toc {{ position:sticky; top:16px; border:1px solid var(--line); border-radius:6px; padding:16px; background:#121a1c; max-height:calc(100vh - 32px); overflow:auto; }}
.toc b {{ display:block; font-family:'JetBrains Mono',monospace; color:var(--gold); margin-bottom:10px; }}
.toc a {{ display:block; color:#c8d0ca; padding:8px 0; border-top:1px solid rgba(255,255,255,.06); font-size:13px; }}
main section {{ border-top:1px solid var(--line); padding:34px 0; }}
h2 {{ font-family:'Cormorant Garamond',serif; font-size:38px; margin:0 0 12px; }}
h3 {{ margin:22px 0 8px; color:var(--gold); font-size:18px; }}
.grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
.box {{ border:1px solid var(--line); background:var(--panel); border-radius:6px; padding:18px; }}
.box.good {{ border-color:rgba(123,200,143,.55); }}
.box.warn {{ border-color:rgba(215,180,106,.65); }}
.box.bad {{ border-color:rgba(225,119,111,.55); }}
.box b.label {{ display:block; font-family:'JetBrains Mono',monospace; color:var(--gold); font-size:12px; letter-spacing:.08em; text-transform:uppercase; margin-bottom:8px; }}
table {{ width:100%; border-collapse:collapse; margin:14px 0 22px; font-size:13px; }}
th,td {{ border:1px solid var(--line); padding:9px 10px; vertical-align:top; }}
th {{ color:var(--gold); background:#131c1f; font-family:'JetBrains Mono',monospace; text-align:left; }}
td {{ background:rgba(255,255,255,.02); }}
.decision td:last-child {{ color:var(--gold); font-weight:800; }}
.gallery {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:16px; }}
figure {{ margin:0; border:1px solid var(--line); background:#0d1213; border-radius:6px; overflow:hidden; }}
figure img {{ display:block; width:100%; background:#fff; }}
figcaption {{ padding:12px 14px; color:#c7d0c8; font-size:13px; }}
code {{ font-family:'JetBrains Mono',monospace; color:#f2d58a; }}
.small {{ color:var(--muted); font-size:13px; }}
@media (max-width: 980px) {{
  .wrap {{ padding:18px; }}
  .layout {{ grid-template-columns:1fr; }}
  .toc {{ position:static; max-height:none; }}
  .stats {{ grid-template-columns:repeat(2,minmax(0,1fr)); }}
  .grid2,.gallery {{ grid-template-columns:1fr; }}
  h1 {{ font-size:56px; }}
}}
@media (max-width: 560px) {{
  .stats {{ grid-template-columns:1fr; }}
  h1 {{ font-size:44px; }}
}}
</style>
</head>
<body>
<div class="wrap">
<header class="hero">
  <div class="kicker">Paper 3 · ICI vulnerability · status dossier · 2026-05-09</div>
  <h1>ICI는 가능성은 강하지만 predictor는 아니다</h1>
  <p class="lead">현재 ICI story의 중심은 <b>thyroid ICI response prediction</b>이 아니라, DM1/ATC/RAI-refractory 축에서 보이는 <b>immunogenomic vulnerability</b>다. 외부 ICI cohort에서는 IFNG, HLA-I, checkpoint, cytolytic 축이 response biology와 맞지만, thyroid-treated cohort가 없으므로 claim guard가 논문의 생명줄이다.</p>
  <div class="crumbs">
    <span>Track A frozen</span><span>Track B-lite read-only synthesis</span><span>Phase C v2 n=421</span><span>HLA gene-expression module, not allele genotype</span><span>No thyroid ICI response predictor claim</span>
  </div>
  <div class="stats">{cards_html}</div>
</header>

<div class="layout">
<aside class="toc">
  <b>읽는 순서</b>
  <a href="#verdict">1. Verdict</a>
  <a href="#thyroid">2. Thyroid-side evidence</a>
  <a href="#external">3. External ICI coherence</a>
  <a href="#hla">4. HLA-I load-bearing axis</a>
  <a href="#braf">5. BRAF-cPTC candidate layer</a>
  <a href="#figures">6. Figure gallery</a>
  <a href="#decision">7. Decision matrix</a>
  <a href="#sources">8. Sources</a>
</aside>
<main>

<section id="verdict">
  <h2>1. Verdict</h2>
  <div class="grid2">
    <div class="box good"><b class="label">What is solid</b>ATC and post-RAI thyroid tumors move toward dedifferentiation plus immune activation. External ICI cohorts independently validate that the T-cell-inflamed components of the axis are response-associated.</div>
    <div class="box warn"><b class="label">What is not solid</b>There is no public thyroid ICI-treated raw RNA-seq cohort in this analysis. Therefore Paper 3 cannot claim response prediction or treatment selection.</div>
  </div>
  <h3>Allowed vs forbidden wording</h3>
  {table(["Allowed", "Forbidden"], [
      ["ICI vulnerability / ICI-readiness / immunogenomic prioritization", "thyroid ICI response predictor"],
      ["external coherence in urothelial + melanoma ICI cohorts", "clinical-grade treatment recommendation"],
      ["HLA-I expression module associates with response biology", "HLA allele genotype association"],
      ["candidate subgroup for future thyroid ICI validation", "validated thyroid ICI responder group"],
  ], "decision")}
</section>

<section id="thyroid">
  <h2>2. Thyroid-side evidence</h2>
  <p>Track B-lite gives the thyroid-specific biological prior: advanced/dedifferentiated thyroid tumors are not simply immune-cold. They are lineage-collapsed, checkpoint-positive, HLA-positive, and heavily myeloid-skewed, while TLS is weak or paradoxical.</p>
  {table(["Module", "ATC vs non-ATC d", "Sign"], atc_table_rows)}
  <h3>RAI-refractory bridge</h3>
  <p>GSE151179 provides a partial bridge from differentiated PTC toward the ATC-like axis. The strongest signal is thyroid-differentiation loss; HLA-II and myeloid move in the same direction but remain modest and underpowered.</p>
  {table(["Module", "post-RAI minus pre-RAI d", "p", "Direction"], rai_table_rows)}
  <p class="small">Phase 4 pooled matrix: n={esc(phase4["pooled_v3_n_samples"])} samples across {esc(len(phase4["pooled_v3_cohorts"]))} cohorts. Sampling timepoint is a proxy, not verified clinical ICI eligibility.</p>
</section>

<section id="external">
  <h2>3. External ICI coherence</h2>
  <p>Phase C v2 moved the older Hugo+Riaz-only result to four pre-treatment cohorts: IMvigor210, GSE176307, Riaz GSE91061, and MGH GSE115821. The composite DM1 scores are 4/4 directionally positive but not conventionally significant; the component T-cell modules carry the statistically clean signal.</p>
  {table(["Score", "Verdict", "OR per +1z", "p", "Sign+", "Median HR OS"], phase_c_rows)}
  <h3>OS-bearing cohorts</h3>
  {table(["Cohort", "Score", "HR per +1z", "95% CI", "p"], surv_rows)}
</section>

<section id="hla">
  <h2>4. HLA-I load-bearing axis</h2>
  <p>Track18 is the most useful 5/8 update for ICI: it separates HLA-I from HLA-II. HLA-I is response-associated and OS-protective; HLA-II is null at meta level. This means the response signal is class-I/CD8/MHC-I aligned, not a generic HLA-high artifact.</p>
  {table(["Score", "OR", "95% CI", "p", "Sign+"], hla_meta_rows)}
  <h3>DM1 × HLA-I interaction and composite revision</h3>
  <div class="grid2">
    <div class="box warn"><b class="label">Interaction</b>DM1 × HLA-I is positive in 3/3 evaluable cohorts, pooled OR={fnum(track18['interaction_pooled']['OR_int'], 2)}, p={fp(track18['interaction_pooled']['p_int'])}. Directionally useful, not significant.</div>
    <div class="box good"><b class="label">Composite v3 hint</b>Adding HLA-I to the DM1-inflam composite improves OR from 1.20 to {fnum(next(r for r in hla_extended if r['score']=='DM1_inflam_plus_HLA2')['pooled_OR'], 2)} and p from 0.13 to {fp(next(r for r in hla_extended if r['score']=='DM1_inflam_plus_HLA2')['p'])}.</div>
  </div>
  {table(["Composite", "OR", "95% CI", "p", "Sign+"], hla_extended_rows)}
  <h3>Within DM1 strata</h3>
  {table(["Cohort", "DM1 stratum", "n", "R", "HLA-I OR", "p"], hla_strata_rows)}
</section>

<section id="braf">
  <h2>5. BRAF-cPTC candidate layer</h2>
  <p>The 5/9 BRAF sprint adds a thyroid-internal candidate layer, not response validation. DM1 BRAF-cPTC is TIS-high, antigen-presentation-high, and checkpoint-rich; across extra-thyroid ICI cohorts, the HT13 × HLA-I adaptive phenotype performs better than HT13 alone.</p>
  {table(["Layer", "Metric", "Read"], [
      ["H27 BRAF-cPTC", f"TIS-Ayers d={fnum(h27_head['TIS_Ayers']['d_DM1_minus_DM2'], 2)}, FDR={fp(h27_head['TIS_Ayers']['BH_FDR'])}", "DM1 is an ICI-candidate phenotype within BRAF-cPTC"],
      ["H27 likely fraction", f"{h27_head['ici_likely_responder_n']}/110 total; {100*h27_head['ici_likely_responder_frac_DM1']:.1f}% of DM1", "candidate enrichment only"],
      ["R5 adaptive", f"OR={fnum(r5_adapt['mh_OR'], 2)} [{fnum(r5_adapt['lo'], 2)}, {fnum(r5_adapt['hi'], 2)}], p={fp(r5_adapt['p'])}", "HT13 × HLA-I, 5/5 sign"],
      ["R5 HT13 alone", f"OR={fnum(r5_head['pooled']['HT13']['pooled_OR'], 2)}, p={fp(r5_head['pooled']['HT13']['p'])}", "directional but not significant"],
      ["R10 PD-L1 stratum", f"DM1-adaptive OR={fnum(r10_head['PDL1__DM1_adaptive_HT13xHLA1']['OR'], 2)}, p={fp(r10_head['PDL1__DM1_adaptive_HT13xHLA1']['p'])}", "confounded by urothelial tumor type"],
      ["R10 PD-L1 IHC", "IC2+ beats transcript panels in IMvigor210", "HT13 does not add independent response info beyond PD-L1 IHC"],
  ])}
</section>

<section id="figures">
  <h2>6. Figure gallery</h2>
  <div class="gallery">{gallery_html}</div>
</section>

<section id="decision">
  <h2>7. Decision matrix</h2>
  {table(["Question", "Current answer", "Disposition"], [
      ["Can Paper 3 be revived as a response-predictor paper?", "No thyroid ICI-treated raw RNA-seq cohort is available.", "NO-GO"],
      ["Can Paper 3 be framed as vulnerability/readiness?", "Yes. Thyroid-side biology plus external ICI coherence supports a prioritized subgroup hypothesis.", "GO"],
      ["What is the strongest molecular gate?", "HLA-I expression module with IFNG/checkpoint/cytolytic support.", "KEEP"],
      ["What should be demoted?", "HLA-II as response predictor, HLA-G ICI response, HT13-alone response marker.", "DEMOTE"],
      ["What is the best near-term extension?", "Gide FASTQ alignment or additional processed UC/melanoma cohort to power DM1/HLA-I interaction.", "OPTIONAL"],
      ["What unlocks clinical-grade language?", "Thyroid ICI-treated cohort with raw RNA-seq and outcome labels.", "BLOCKED"],
      ["What is Track B status?", "Track A frozen; Track B still requires explicit unlock.", "FROZEN"],
      ["What belongs in Paper 1/2?", "Only brief references. ICI storyline stays Paper 3/Paper 11 territory.", "BOUNDARY"],
  ], "decision")}
</section>

<section id="sources">
  <h2>8. Sources and paths</h2>
  {table(["Artifact", "Path"], source_rows)}
  <p class="small">Generated by <code>{rel(Path(__file__))}</code>. HTML deployed to <code>{rel(HTML)}</code> and <code>{LIVE_HTML}</code>. Markdown summary at <code>{rel(REPORT)}</code>.</p>
</section>

</main>
</div>
</div>
</body>
</html>
"""

    HTML.write_text(html_doc)
    shutil.copy2(HTML, LIVE_HTML)


if __name__ == "__main__":
    build()
