#!/usr/bin/env python3
from __future__ import annotations

import html
import math
from pathlib import Path

import pandas as pd

from v17_common import FIG, RPT, TAB

PAGES_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/reports/html/pages")


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(TAB / name, sep="\t")


def metric_cards(cards: list[tuple[str, str]]) -> str:
    return "".join(
        f"<div class='tile'><b>{html.escape(v)}</b><span>{html.escape(k)}</span></div>"
        for k, v in cards
    )


def figure_card(src: str, title_kr: str, caption_kr: str, height: int = 340) -> str:
    return f"""
    <div class="fig-card">
      <div class="fig-card__head">
        <h3>{html.escape(title_kr)}</h3>
        <p>{html.escape(caption_kr)}</p>
      </div>
      <iframe loading="lazy" src="{src}" title="{html.escape(title_kr)}" style="height:{height}px"></iframe>
      <div class="fig-card__links"><a href="{src}" target="_blank" rel="noopener">새 창으로 열기</a></div>
    </div>
    """


def simple_table(df: pd.DataFrame, table_id: str) -> str:
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = "".join(f"<td>{html.escape(str(v))}</td>" for v in r.tolist())
        rows.append(f"<tr>{cells}</tr>")
    body = "".join(rows)
    return f"<table id='{table_id}' class='display compact stripe'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def layout_css() -> str:
    return """
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Noto+Sans+KR:wght@400;500;700&display=swap" rel="stylesheet">
    <link rel="stylesheet" href="../html/assets/css/main.css">
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/datatables.net-dt@1.13.8/css/jquery.dataTables.min.css">
    <style>
      :root{--bg:#07132b;--fg:#f8fafc;--accent:#f5a623}
      body.dark{font-family:"Inter","Noto Sans KR",-apple-system,BlinkMacSystemFont,sans-serif}
      body.dark h1,body.dark h2,body.dark h3,body.dark h4,body.dark p,body.dark li,body.dark td,body.dark th,body.dark a,body.dark span,body.dark div{font-family:"Inter","Noto Sans KR",-apple-system,BlinkMacSystemFont,sans-serif}
      body.dark .mono, body.dark code{font-family:"JetBrains Mono",ui-monospace,monospace}
      .v17-wrap{padding-bottom:72px}
      .v17-hero{padding:36px 32px;border-radius:24px;border:1px solid rgba(255,255,255,.08);background:
        radial-gradient(circle at 15% 18%, rgba(245,166,35,.16), transparent 28%),
        radial-gradient(circle at 85% 14%, rgba(20,184,166,.10), transparent 22%),
        linear-gradient(180deg, rgba(15,23,42,.96), rgba(2,6,23,.92))}
      .v17-hero h1{font-size:clamp(36px,4.8vw,64px);font-weight:300;letter-spacing:-.04em;line-height:1.04;margin:10px 0 10px}
      .v17-hero .subtitle{font-size:18px;color:#cbd5e1;margin:0 0 16px}
      .v17-hero .meta{display:flex;gap:12px;flex-wrap:wrap;color:#94a3b8;font-size:13px}
      .hero-summary{margin-top:18px;padding-left:18px}
      .hero-summary li{margin:6px 0;color:#dbe4f0}
      .metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px;margin-top:22px}
      .tile{padding:16px 18px;border-radius:16px;background:rgba(255,255,255,.04);border:1px solid rgba(255,255,255,.08)}
      .tile b{display:block;font-size:30px;color:#f8fafc}
      .tile span{display:block;margin-top:8px;font-size:11px;letter-spacing:.12em;text-transform:uppercase;color:#94a3b8}
      .hero-actions{display:flex;gap:10px;flex-wrap:wrap;margin-top:18px}
      .hero-actions a{display:inline-flex;align-items:center;padding:10px 14px;border-radius:999px;border:1px solid rgba(255,255,255,.12);background:rgba(255,255,255,.04)}
      .hero-actions a.primary{background:#f5a623;color:#111;border-color:#f5a623;font-weight:700}
      .v17-section{margin-top:22px}
      .v17-panel{padding:24px;border-radius:22px;background:rgba(15,23,42,.72);border:1px solid rgba(255,255,255,.08)}
      .v17-panel h2{margin:0 0 10px;font-size:30px}
      .v17-panel p{color:#cbd5e1}
      .two-col{display:grid;grid-template-columns:minmax(0,1.35fr) minmax(280px,.65fr);gap:18px}
      .card-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:16px}
      .fig-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(400px,1fr));gap:16px}
      .fig-card{border:1px solid rgba(255,255,255,.08);border-radius:18px;overflow:hidden;background:rgba(2,6,23,.42)}
      .fig-card__head{padding:14px 16px 10px;border-bottom:1px solid rgba(255,255,255,.06)}
      .fig-card__head h3{margin:0 0 4px;font-size:18px}
      .fig-card__head p{margin:0;font-size:14px;color:#94a3b8}
      .fig-card iframe{display:block;width:100%;border:0;background:#030712}
      .fig-card__links{padding:10px 16px 16px}
      .fig-card__links a{color:#7dd3fc}
      .note-box,.warn-box{padding:16px 18px;border-radius:18px}
      .note-box{background:rgba(14,165,233,.08);border:1px solid rgba(14,165,233,.2)}
      .warn-box{background:rgba(245,166,35,.08);border:1px solid rgba(245,166,35,.22)}
      .html-table-wrap{overflow:auto}
      table.dataTable{width:100%!important}
      table.dataTable thead th{background:#0f172a;color:#f8fafc;border-bottom:1px solid rgba(255,255,255,.15)}
      table.dataTable tbody td{color:#cbd5e1;background:#071120}
      .badge{display:inline-block;padding:4px 8px;border-radius:999px;background:rgba(245,166,35,.12);color:#f5a623;font-size:11px;font-weight:700;letter-spacing:.08em;text-transform:uppercase}
      .list-tight li{margin:4px 0}
      .footer-links{display:flex;gap:10px;flex-wrap:wrap;margin-top:12px}
      @media (max-width:980px){.two-col{grid-template-columns:1fr}.fig-grid{grid-template-columns:1fr}}
    </style>
    """


def subnav_html(active: str) -> str:
    items = [
        ("index.html", "Report"),
        ("methods.html", "Methods"),
        ("tables.html", "Tables"),
        ("figures.html", "Figures"),
    ]
    links = "".join(
        f"<a class='btn sm{' primary' if label.lower()==active else ' ghost'}' href='{href}'>{label}</a>"
        for href, label in items
    )
    return f"<div class='hero-actions' style='margin-top:12px'>{links}</div>"


def nav_html(active: str = "v17") -> str:
    links = [
        ("../html/index.html", "홈"),
        ("../html/pages/01_overview.html", "개요"),
        ("../html/pages/13_reports.html", "리포트"),
        ("../html/pages/v12_biological_evidence.html", "v12 Biology"),
        ("../html/pages/v13_drug_discovery.html", "v13 Drug Discovery"),
        ("../html/pages/v15_neurips_theory.html", "v15 NeurIPS"),
        ("../v17/index.html", "v17 New"),
    ]
    nav = "".join(
        f"<a class='nav-link{' active' if label=='v17 New' and active=='v17' else ''}' href='{href}' role='menuitem'>{label}</a>"
        for href, label in links
    )
    return f"""
    <nav class="topnav" aria-label="주요 메뉴">
      <div class="topnav-inner">
        <a class="brand" href="../html/index.html" aria-label="홈으로"><span class="brand-dot" aria-hidden="true"></span><span>THYROID DASH</span></a>
        <div class="nav-links" role="menubar">{nav}</div>
        <div class="nav-actions"><a class="btn sm ghost" href="../v17/index.html">v17 Hub</a></div>
      </div>
    </nav>
    """


def wrapper(title: str, body: str, toc: list[tuple[str, str]]) -> str:
    sidebar = "".join(f"<a href='#{sid}'>{label}</a>" for sid, label in toc)
    quick = "".join(f"<a href='#{sid}'>{i+1}. {label}</a>" for i, (sid, label) in enumerate(toc))
    return f"""<!doctype html>
<html lang="ko">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="dark light">
  <meta name="theme-color" content="#07132b">
  <title>{html.escape(title)}</title>
  <meta name="description" content="{html.escape(title)}">
  <link rel="icon" href="../html/assets/img/favicon.svg" type="image/svg+xml">
  {layout_css()}
</head>
<body class="dark">
  {nav_html()}
  <div class="layout v17-wrap">
    <aside class="sidebar glass" aria-label="섹션 목차">
      <h3>Section TOC</h3>
      {sidebar}
    </aside>
    <main class="content" id="main-content">
      {body}
    </main>
    <aside class="toc glass" aria-label="퀵 점프">
      <h3>Quick Jump</h3>
      {quick}
    </aside>
  </div>
  <script src="https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/datatables.net@1.13.8/js/jquery.dataTables.min.js"></script>
  <script>
  document.querySelectorAll('table.display').forEach((el)=>{{ new DataTable(el, {{pageLength:10, order:[]}}); }});
  </script>
</body>
</html>"""


def build_index() -> str:
    dark_summary = read_tsv("dark_matter_cluster_clinical.tsv")
    driver_summary = read_tsv("driver_landscape_v17_summary.tsv")
    dial = read_tsv("dial_audit_v17.tsv")
    sm = read_tsv("sample_master_v17_full.tsv")
    best_k = 2
    dia_auc = f"{dial['DIA_AUC'].min():.3f}-{dial['DIA_AUC'].max():.3f}"
    tcga = sm[(sm["dataset"] == "TCGA-THCA") & (sm["normal_vs_tumor"] != "normal")].copy()
    mismatch = tcga[
        ((tcga["driver_anchor_v17"] == "BRAF") & (tcga["molecular_subtype"] != "BRAF_like"))
        | ((tcga["driver_anchor_v17"] == "RAS") & (tcga["molecular_subtype"] != "RAS_like"))
    ]
    mismatch_rate = 100 * len(mismatch) / max(1, len(tcga[tcga["driver_anchor_v17"].isin(["BRAF", "RAS"])]))
    dark_core = dark_summary[["v17_dark_cluster", "n_samples", "tds16_mean", "rai_mean", "age_mean", "BRAF_like", "RAS_like", "unknown"]].copy()
    dark_core.columns = ["cluster", "n", "TDS16 mean", "RAI mean", "age mean", "BRAF_like", "RAS_like", "unknown"]
    fig_names = sorted([p.name for p in FIG.glob("*.html")])
    table_names = sorted([p.name for p in TAB.glob("*.tsv")])
    fig_catalog = "".join(
        f"<div class='tile'><b style='font-size:16px'>{html.escape(n.replace('.html',''))}</b><span>figure</span><div class='footer-links'><a href='../../results/v17/figs/{n}'>열기</a></div></div>"
        for n in fig_names
    )
    table_catalog = "".join(
        f"<div class='tile'><b style='font-size:16px'>{html.escape(n)}</b><span>table</span><div class='footer-links'><a href='../../results/v17/tables/{n}'>TSV</a></div></div>"
        for n in table_names
    )
    body = f"""
    <section class="v17-hero">
      <div class="subtle mono">v17 / HTML REPORT</div>
      <h1>v17 — BRAF/RAS를 넘어선 갑상선암 분자 분류</h1>
      <p class="subtitle">Cross-cohort robust subtyping of papillary thyroid carcinoma</p>
      <div class="meta"><span>국승호</span><span>·</span><span>2026-04-26</span><span>·</span><span>THYROID DASH deliverable</span></div>
      {subnav_html('report')}
      <ul class="hero-summary">
        <li>TCGA에서 BRAF/RAS-negative tumor 178명은 하나의 residual bucket이 아니라 DM1/DM2 두 군으로 안정적으로 갈렸다.</li>
        <li>Consensus stability는 0.9788로 높았고, external transfer proxy에서도 DIA-AUC 0.969-1.000을 유지했다.</li>
        <li>ComBat 후 identifiability는 1.000 근처에서 0.157-0.258까지 급락해, biology transfer와 batch fingerprint 제거가 동시에 일어남을 보였다.</li>
        <li>Trajectory에서는 TPO, TG, DIO1/2, PAX8 같은 thyroid differentiation axis가 pseudotime 증가와 함께 무너졌다.</li>
        <li>다만 TERT promoter와 raw fusion은 아직 비관측/대리변수 단계라서, 다음 단계는 Korean cohort + wet validation이다.</li>
      </ul>
      <div class="metric-grid">
        {metric_cards([("Dark Matter n", "178"), ("cluster K", str(best_k)), ("DIAL AUC", dia_auc), ("cohort", "5")])}
      </div>
      <div class="hero-actions">
        <a class="primary" href="methods.html">방법론 상세</a>
        <a href="tables.html">18개 테이블 보기</a>
        <a href="figures.html">Figure 갤러리</a>
      </div>
    </section>

    <section id="problem" class="v17-section v17-panel">
      <span class="badge">Section 1</span>
      <h2>왜 BRAF/RAS 이분법만으로 부족한가</h2>
      <div class="two-col">
        <div>
          <p>v14까지의 통합 파이프라인은 갑상선암을 BRAF-like와 RAS-like 축 위에서 안정적으로 정리했고, 5-cohort harmonization과 compact panel, DIAL audit, ComBat-seq robustness까지 이미 구축했다. 그런데 TCGA-THCA tumor 513명 안에서도 BRAF와 RAS로 깔끔하게 anchoring 되지 않는 표본이 178명 남았다. 이 집단은 이전 프레임에서는 설명의 끝이었지만, v17에서는 오히려 새로운 시작점이 된다.</p>
          <p>즉 v17의 질문은 단순하다. “BRAF/RAS로 설명되지 않는 표본은 진짜 noise인가, 아니면 아직 이름이 없던 subtype space인가?” 이번 결과는 후자 쪽으로 강하게 기운다.</p>
          <p>현재 TCGA tumor 기준 v14 molecular subtype과 v17 driver anchor의 불일치율은 <b>{mismatch_rate:.1f}%</b> 이다. 이는 surrogate subtype과 driver taxonomy가 완전히 같은 이야기를 하지 않는다는 점을 보여준다.</p>
        </div>
        <div class="warn-box">
          <h3>v14 → v17 한 줄 전환</h3>
          <p>v14가 “BRAF_like/RAS_like가 외부 코호트에서도 유지되는가”를 물었다면, v17은 “그 밖에 남는 표본도 구조를 갖는가”를 묻는다.</p>
          <p><b>BRS surrogate misclassification rate:</b> {mismatch_rate:.1f}%</p>
        </div>
      </div>
      <div class="fig-grid" style="margin-top:16px">
        {figure_card("../../results/v17/figs/driver_landscape_pie.html", "Driver-anchor 분포 도넛", "TCGA 및 통합 sample master 기준 v17 driver taxonomy의 상대적 분포를 보여준다.")}
        {figure_card("../../results/v17/figs/v14_vs_v17_confusion.html", "v14 subtype vs v17 driver confusion", "BRAF_like/RAS_like 라벨과 확장 driver taxonomy가 어디서 어긋나는지 요약한다.")}
      </div>
    </section>

    <section id="dark-matter" class="v17-section v17-panel">
      <span class="badge">Section 2</span>
      <h2>Dark Matter 코호트</h2>
      <p>TCGA-THCA tumor 중 `driver_anchor != BRAF, RAS` 조건으로 추출한 178명은 이번 분석의 중심 cohort다. TierA expression을 standardization 후 PCA로 축약하고, 반복 seed consensus clustering으로 내부 구조를 찾았다. 결과는 `DM1=109`, `DM2=69` 이고, 합의 안정도는 `0.9788` 이다.</p>
      <p>DM1은 상대적으로 lower differentiation, lower RAI, younger age, 그리고 v14 기준 BRAF_like 잔존 성분이 더 많다. DM2는 TDS/RAI 평균이 더 높고, RAS_like 비율이 상대적으로 높다. 즉 Dark Matter가 완전히 무질서한 드라이버-음성 군이 아니라, BRAF-like 잔향과 RAS-like 잔향이 다른 두 군으로 갈라지는 해석이 가능하다.</p>
      <p>아래 figure 네 장이 이 섹션의 핵심이다. UMAP은 공간 분리, consensus heatmap은 구조 안정도, marker heatmap은 분자 차이, RAI boxplot은 생물학적 해석 가능성을 보여준다.</p>
      <div class="fig-grid">
        {figure_card("../../results/v17/figs/dark_matter_umap.html", "Dark Matter UMAP", "DM1과 DM2가 저차원 embedding에서 어떻게 분리되는지 보여준다.")}
        {figure_card("../../results/v17/figs/dark_matter_consensus_heatmap.html", "Consensus heatmap", "30개 seed 반복에서 sample pair가 얼마나 일관되게 함께 묶였는지 나타낸다.")}
        {figure_card("../../results/v17/figs/dark_matter_marker_heatmap.html", "Marker heatmap", "각 cluster의 top marker가 one-vs-rest 비교에서 어떻게 갈리는지 요약한다.")}
        {figure_card("../../results/v17/figs/dark_matter_rai_boxplot.html", "RAI boxplot", "iodine uptake program이 cluster별로 어떻게 달라지는지 보여준다.")}
      </div>
      <div class="html-table-wrap" style="margin-top:18px">
        {simple_table(dark_core.round(3), "dark-core-table")}
      </div>
    </section>

    <section id="drivers" class="v17-section v17-panel">
      <span class="badge">Section 3</span>
      <h2>Driver Landscape 확장</h2>
      <p>v17 driver landscape는 로컬 TCGA mutation summary와 기존 fusion proxy를 합쳐 만든 빠른 재현판이다. 전체 1,509 sample 기준으로 `unknown=1072`, `BRAF=343`, `RAS=74`가 상위 그룹이며, 소수의 `PAX8PPARG`, `DICER1/EIF1AX/PPM1D`, `RET`, `ALK`, `NTRK`, `TP53` 축이 뒤를 잇는다.</p>
      <div class="fig-grid">
        {figure_card("../../results/v17/figs/driver_landscape_oncoplot.html", "Driver oncoplot", "TCGA tumor에서 sample-by-driver membership을 heatmap 형태로 정리한 그림.")}
        {figure_card("../../results/v17/figs/driver_landscape_pie.html", "Driver pie", "v17 taxonomy의 상대 빈도 구성.")}
        {figure_card("../../results/v17/figs/dark_matter_fusion_overlay.html", "Dark Matter overlay", "Dark Matter cluster 안에서 v17 driver 항목이 어떻게 섞이는지 보여준다.")}
      </div>
      <div class="warn-box" style="margin-top:16px">
        <h3>Caveat — fusion proxy 재사용</h3>
        <p>이번 실행은 raw TCGA fusion flat file 재호출이 아니라 기존 `v3_fusion_anchor_tcga.tsv` proxy를 재사용했다. 따라서 sensitivity 한계가 있고, rare fusion class는 under-call 되었을 가능성이 있다.</p>
      </div>
    </section>

    <section id="dial" class="v17-section v17-panel">
      <span class="badge">Section 4</span>
      <h2>Cross-Cohort 검증 (DIAL)</h2>
      <p>TCGA에서 학습한 DM1/DM2 signature를 external cohort로 transfer했을 때, separation은 거의 완벽하게 유지됐다. `GSE27155`, `GSE76039` 에서 DIA-AUC는 0.969-1.000 범위를 보였다. 반면 ComBat 후 cohort identifiability는 1.000 근처에서 0.2577, 0.1573으로 무너졌다.</p>
      <p><b>한 줄 요약:</b> biology-preserving transfer는 강하지만, batch correction은 그 과정에서 cohort fingerprint를 공격적으로 지운다.</p>
      <div class="fig-grid">
        {figure_card("../../results/v17/figs/dial_forest_plot.html", "DIAL forest plot", "cohort별, cluster별 direction-invariant AUC를 요약한다.")}
        {figure_card("../../results/v17/figs/dial_identifiability_heatmap.html", "Identifiability heatmap", "ComBat 전후에 TCGA-vs-external 구분 가능성이 얼마나 달라지는지 보여준다.")}
        {figure_card("../../results/v17/figs/dial_combat_comparison.html", "ComBat comparison", "identifiability 분포를 correction 전후로 비교한다.")}
      </div>
    </section>

    <section id="trajectory" class="v17-section v17-panel">
      <span class="badge">Section 5</span>
      <h2>탈분화 Trajectory</h2>
      <p>Trajectory는 공통 TierA gene 50개, 총 645 tumor sample에서 계산했다. pseudotime 증가와 함께 `TPO`, `TG`, `DIO1/2`, `PAX8` 같은 thyroid differentiation module은 줄고, `DUSP5`, `MET`, `LOX` 같은 축은 상대적으로 올라갔다. 이는 일부 Dark Matter cluster가 단순 driver-negative 군이 아니라 탈분화 쪽으로 기운 진행축 안에 놓일 수 있음을 시사한다.</p>
      <p>현재 구현은 full 5-cohort complete state가 아니라 usable gene-overlap subset 기반이다. 그럼에도 PTC에서 더 뒤쪽 tail로 가는 표본이 어떤 표현형을 띠는지 보는 bonus figure로는 충분히 강하다.</p>
      <div class="fig-grid">
        {figure_card("../../results/v17/figs/trajectory_umap_3d.html", "Trajectory UMAP 3D", "pseudotime으로 색칠한 global tumor manifold.")}
        {figure_card("../../results/v17/figs/trajectory_tds_gradient.html", "TDS gradient", "pseudotime에 따라 differentiation score가 어떻게 변하는지 보여준다.")}
        {figure_card("../../results/v17/figs/trajectory_high_risk_cluster.html", "High-risk tail", "pseudotime 상위 tail이 manifold 어디에 놓이는지 보여준다.")}
      </div>
      <div class="note-box" style="margin-top:16px">
        <p><b>PTC→PDTC→ATC 질문에 대한 현재 답:</b> 이번 로컬 run은 histology full projection보다는 progression-like axis 재현이 중심이었고, 그 축에서 DM1 일부가 상대적으로 뒤쪽 tail과 더 가까워 보인다. definitive claim은 Korean cohort와 higher-risk histology 확장이 필요하다.</p>
      </div>
    </section>

    <section id="limits" class="v17-section v17-panel">
      <span class="badge">Section 6</span>
      <h2>한계와 다음 단계</h2>
      <ul class="list-tight">
        <li><b>TERT promoter 비관측:</b> 현재 로컬 GDC WXS MAF는 canonical promoter hotspot을 사실상 포착하지 못했다.</li>
        <li><b>Fusion proxy 사용:</b> raw fusion caller가 아니라 기존 proxy table을 썼기 때문에 민감도 한계가 있다.</li>
        <li><b>Wet validation 부재:</b> 현재 결과는 computational subtype proposal이다.</li>
        <li><b>Korean cohort access 대기:</b> 외부 검증과 임상 일반화를 위해 가장 중요한 다음 단계다.</li>
      </ul>
    </section>

    <section id="venue" class="v17-section v17-panel">
      <span class="badge">Section 7</span>
      <h2>Paper venue 진단</h2>
      <div class="card-grid">
        <div class="warn-box"><h3>현 상태로 가능</h3><p>Bioinformatics, npj Precision Oncology, Genome Medicine</p></div>
        <div class="note-box"><h3>도전 가능</h3><p>Nature Communications — Korean cohort 추가 시</p></div>
        <div class="warn-box"><h3>매우 어려움</h3><p>Nature Cancer — wet validation + large Korean cohort + 12개월 추가 필요</p></div>
      </div>
      <p style="margin-top:16px"><b>한 줄 진단:</b> “v17 그대로 → Bioinformatics 70%, Nature Cancer 5%”.</p>
    </section>

    <section id="footer" class="v17-section v17-panel">
      <h2>재현성</h2>
      <p>Python 3.12 / `project/.venv` / local TCGA MAF cache / Plotly interactive HTML. 이번 페이지는 UTF-8 한국어 HTML로 직접 렌더링되며, figure와 table은 모두 같은 서버 루트에서 정적 제공된다.</p>
      <div class="footer-links">
        <a href="../html/index.html">v14 dashboard</a>
        <a href="../html/pages/v15_neurips_theory.html">v15 NeurIPS draft</a>
        <a href="tables.html">18 tables</a>
        <a href="figures.html">32 figures</a>
      </div>
    </section>
    <section id="all-figures" class="v17-section v17-panel">
      <h2>전체 figure catalog</h2>
      <p>이번 run에서 생성된 모든 interactive figure를 한 번에 모아 둔 섹션이다. 핵심 23개뿐 아니라 supplemental figure도 함께 포함한다.</p>
      <div class="card-grid">{fig_catalog}</div>
    </section>
    <section id="all-tables" class="v17-section v17-panel">
      <h2>전체 table catalog</h2>
      <p>이번 run의 TSV 산출물 18개를 한 페이지에서 접근할 수 있도록 모았다.</p>
      <div class="card-grid">{table_catalog}</div>
    </section>
    """
    return wrapper("v17 — BRAF/RAS를 넘어선 갑상선암 분자 분류", body, [
        ("problem", "문제 제기"),
        ("dark-matter", "Dark Matter"),
        ("drivers", "Driver landscape"),
        ("dial", "DIAL"),
        ("trajectory", "Trajectory"),
        ("limits", "한계"),
        ("venue", "Venue"),
        ("all-figures", "전체 figure"),
        ("all-tables", "전체 table"),
    ])


def build_methods() -> str:
    body = """
    <section class="v17-hero">
      <div class="subtle mono">v17 / METHODS</div>
      <h1>v17 방법론 상세</h1>
      <p class="subtitle">작동 원리, 알고리즘, 수학, 통계의 한국어 설명</p>
      {subnav_html('methods')}
    </section>
    <section id="m1" class="v17-section v17-panel"><h2>1. Dark Matter consensus clustering</h2><p>TCGA-THCA tumor 중 driver_anchor가 BRAF/RAS가 아닌 178명을 대상으로 TierA expression을 표준화한 뒤 PCA 10축으로 줄였다. 표준화는 gene별 평균 0, 분산 1을 맞춰 특정 고분산 유전자가 거리 계산을 지배하지 못하게 하려는 목적이다. PCA는 상관된 gene 집합을 저차 latent axis로 압축해 seed 반복 시 군집 안정도를 높인다.</p><p>그 다음 K=2..10, seed 30회 반복 clustering을 수행하고, sample pair가 같은 cluster에 속한 비율로 consensus matrix를 만들었다. cluster stability는 같은 cluster 내부 pair들의 평균 consensus로 계산했다. 이번 실행에서 K=2, stability=0.9788이 최대였다.</p></section>
    <section id="m2" class="v17-section v17-panel"><h2>2. Marker와 FDR</h2><p>각 cluster는 one-vs-rest 방식으로 비교했다. gene별 차이는 Welch형 t-test를 사용했고, p-value는 Benjamini-Hochberg FDR로 조정했다. Welch형 통계를 쓴 이유는 cluster 간 분산이 같다는 가정을 강하게 두지 않기 위해서다. BH-FDR는 다중검정에서 false discovery proportion을 제어하기 위한 표준 방법이다.</p></section>
    <section id="m3" class="v17-section v17-panel"><h2>3. Signature score</h2><p>RAI score와 TDS16 score는 gene set 평균으로 계산했다. 누락 gene은 0으로 넣지 않고 실제로 존재하는 gene만 평균했다. 그래야 cohort/platform 간 gene coverage 차이가 인위적인 페널티로 번지지 않는다.</p></section>
    <section id="m4" class="v17-section v17-panel"><h2>4. DIAL audit</h2><p>External cohort에는 DM1/DM2의 정답 라벨이 없어서 pseudo-label transfer를 사용했다. 먼저 TCGA에서 cluster-vs-rest logistic regression classifier를 적합하고, 외부 표본은 positive/negative centroid까지의 제곱거리 비교로 임시 라벨을 받는다. 그 pseudo-label에 대해 ROC AUC를 계산하고, 방향 뒤집힘 영향을 없애기 위해 DIA-AUC=max(AUC,1-AUC)를 사용했다.</p><p>Identifiability는 TCGA vs external batch label을 예측하는 AUC다. 값이 1이면 cohort가 거의 완벽히 구분되고, 0.5 근처면 cohort가 섞인 상태다. 이번 실행에서는 transfer AUC는 높았지만, ComBat 후 identifiability는 급감했다.</p></section>
    <section id="m5" class="v17-section v17-panel"><h2>5. ComBat과 trajectory</h2><p>Trajectory는 cohort 간 공통으로 존재하는 TierA gene 교집합 50개만 사용해 NaN을 제거했다. 이후 ComBat normalization을 적용했고, Scanpy에서 PCA -> neighbors -> UMAP -> diffusion map -> DPT 순으로 진행했다. root는 TDS16 score가 가장 높은 sample로 두었다. 즉 가장 thyroid-differentiated 상태를 출발점으로 가정했다.</p><p>Dynamic gene은 pseudotime과의 Spearman 순위상관으로 정리했다. Spearman을 쓴 이유는 선형 관계보다 monotone trend가 더 중요하기 때문이다.</p></section>
    """
    return wrapper("v17 방법론 상세", body, [("m1", "Clustering"), ("m2", "Marker/FDR"), ("m3", "Signature"), ("m4", "DIAL"), ("m5", "Trajectory")])


def build_tables() -> str:
    files = sorted([p.name for p in TAB.glob("*.tsv")])
    table_cards = []
    for i, name in enumerate(files):
        df = read_tsv(name).head(200).copy()
        table_cards.append(f"<section id='t{i}' class='v17-section v17-panel'><h2>{html.escape(name)}</h2><p>성능과 가독성을 위해 처음 200행만 직접 렌더링합니다. 원본 TSV는 아래 링크에서 전체 다운로드 가능합니다.</p><div class='footer-links'><a href='../../results/v17/tables/{name}'>원본 TSV 열기</a></div><div class='html-table-wrap'>{simple_table(df, f'tbl_{i}')}</div></section>")
    body = f"""
    <section class="v17-hero">
      <div class="subtle mono">v17 / TABLES</div>
      <h1>v17 인터랙티브 테이블 뷰</h1>
      <p class="subtitle">18개 TSV를 DataTables 검색/정렬 테이블로 제공</p>
      {subnav_html('tables')}
      <div class="metric-grid">{metric_cards([("TSV", str(len(files))), ("render limit", "200 rows"), ("search", "enabled"), ("sort", "enabled")])}</div>
    </section>
    {''.join(table_cards)}
    """
    toc = [(f"t{i}", name) for i, name in enumerate(files)]
    return wrapper("v17 tables", body, toc)


def build_figures() -> str:
    figs = sorted([p.name for p in FIG.glob("*.html")])
    cards = []
    for i, name in enumerate(figs):
        cards.append(figure_card(f"../../results/v17/figs/{name}", name.replace(".html", ""), "v17 실행 중 생성된 interactive figure.", 260))
    body = f"""
    <section class="v17-hero">
      <div class="subtle mono">v17 / FIGURES</div>
      <h1>v17 figure gallery</h1>
      <p class="subtitle">23 core figure에 supplemental figure를 더해 총 {len(figs)}개 interactive panel 제공</p>
      {subnav_html('figures')}
      <div class="metric-grid">{metric_cards([("interactive figure", str(len(figs))), ("core", "23"), ("supplemental", str(max(0, len(figs)-23))), ("format", "Plotly HTML")])}</div>
    </section>
    <section id="f0" class="v17-section"><div class="fig-grid">{''.join(cards)}</div></section>
    """
    return wrapper("v17 figures", body, [("f0", "전체 figure")])


def build_wrapper_page() -> str:
    return """<!doctype html>
<html lang="ko"><head><meta charset="UTF-8"><meta http-equiv="refresh" content="0; url=../../v17/index.html"><title>v17 redirect</title></head><body>Redirecting to <a href="../../v17/index.html">v17 report</a>.</body></html>"""


def main() -> None:
    RPT.mkdir(parents=True, exist_ok=True)
    (RPT / "index.html").write_text(build_index(), encoding="utf-8")
    (RPT / "methods.html").write_text(build_methods(), encoding="utf-8")
    (RPT / "tables.html").write_text(build_tables(), encoding="utf-8")
    (RPT / "figures.html").write_text(build_figures(), encoding="utf-8")
    (PAGES_DIR / "v17_subtyping.html").write_text(build_wrapper_page(), encoding="utf-8")


if __name__ == "__main__":
    main()
