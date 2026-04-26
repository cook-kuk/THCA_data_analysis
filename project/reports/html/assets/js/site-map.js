/* THYRAI global site map overlay — full-page navigator */
(function () {
  "use strict";

  if (window.ThyraiSiteMap) return;

  var MAP = [
    {
      title: "핵심 시작점",
      desc: "처음 보는 사람을 위한 시작 동선",
      items: [
        { label: "홈", href: "../index.html", note: "전체 연구와 기능 안내의 시작점" },
        { label: "Pipeline", href: "pipeline.html", note: "데이터 → 바이오마커 → 타깃 → 화합물 전체 흐름" },
        { label: "Research Timeline", href: "research_timeline.html", note: "v2부터 v5.1까지 버전 변화와 결론 수정 기록" },
        { label: "Reasoning Story", href: "research_reasoning.html", note: "왜 결론이 바뀌었는지 서사형 설명" },
        { label: "v5.1 Audit", href: "v4a_dial_cross_cancer.html", note: "THCA 특이적 label flip 재평가의 현재 핵심 결과" },
        { label: "v12 Biological Evidence", href: "v12_biological_evidence.html", note: "문헌 근거와 생물학적 연결" },
        { label: "Publications", href: "publications.html", note: "보고서와 공개 문서 허브" },
        { label: "Papers (DIAL)", href: "papers_archive.html", note: "Bioinformatics(OUP) · AAAI 2027 manuscript 아카이브" },
        { label: "Papers Viewer", href: "papers_viewer.html", note: "두 논문 PDF·figures·proof 인터랙티브 뷰어" },
        { label: "Reading Guide", href: "reading_guide.html", note: "5분 / 15분 / 1시간 코스 진입 가이드" }
      ]
    },
    {
      title: "주요 기능",
      desc: "이 사이트가 실제로 제공하는 기능 단위",
      items: [
        { label: "데이터 검증", href: "02_datasets.html", note: "코호트, 플랫폼, 품질, 다운로드 상태 검토" },
        { label: "샘플 탐색", href: "03_sample_master.html", note: "샘플 단위 메타데이터 검색과 confidence 확인" },
        { label: "EDA", href: "05_eda.html", note: "PCA·UMAP·tSNE로 구조와 batch effect 탐색" },
        { label: "점수 시각화", href: "06_scores.html", note: "TDS·BRS proxy·dediff 점수 비교" },
        { label: "ML baseline", href: "07_ml_baseline.html", note: "ROC, PR, calibration, confusion matrix" },
        { label: "패널 비교", href: "08_panel_comparison.html", note: "16/67/112 유전자 패널의 성능·비용 균형 비교" },
        { label: "SHAP 해석", href: "09_shap.html", note: "모델이 어떤 feature를 쓰는지 설명" },
        { label: "유전자 탐색기", href: "10_gene_explorer.html", note: "유전자별 발현·coverage·공발현 구조 확인" },
        { label: "사업성 비교", href: "12_business.html", note: "가격·workflow friction·수익 시뮬레이션" },
        { label: "한계와 실패 조건", href: "14_caveats.html", note: "과대해석을 막기 위한 리스크와 caveat" }
      ]
    },
    {
      title: "v2 상세 탐색기",
      desc: "공개 데이터 기반 THCA 상세 분석 페이지",
      items: [
        { label: "01 Overview", href: "01_overview.html", note: "연구 구조와 샘플 흐름" },
        { label: "02 Datasets", href: "02_datasets.html", note: "데이터셋 탐색기" },
        { label: "03 Sample Master", href: "03_sample_master.html", note: "샘플 단위 표와 QC" },
        { label: "04 Gene Panels", href: "04_gene_panels.html", note: "TDS16 / BRS proxy / TierA67 비교" },
        { label: "05 EDA", href: "05_eda.html", note: "탐색적 분석과 embedding" },
        { label: "06 Scores", href: "06_scores.html", note: "score distribution과 subtype 축" },
        { label: "07 ML Baseline", href: "07_ml_baseline.html", note: "기본 분류 모델 평가" },
        { label: "08 Panel Comparison", href: "08_panel_comparison.html", note: "패널 크기별 성능 비교" },
        { label: "09 SHAP", href: "09_shap.html", note: "모델 설명 가능성" },
        { label: "10 Gene Explorer", href: "10_gene_explorer.html", note: "유전자별 탐색" },
        { label: "11 Cohort Compare", href: "11_cohort_compare.html", note: "TCGA vs GEO 비교" },
        { label: "12 Business", href: "12_business.html", note: "사업/포지셔닝 시뮬레이션" },
        { label: "13 Reports", href: "13_reports.html", note: "markdown 리포트 렌더링" },
        { label: "14 Caveats", href: "14_caveats.html", note: "한계와 실패 조건" }
      ]
    },
    {
      title: "고급 분석 / 확장 페이지",
      desc: "패널·생존·메타파워·약물·양자 등 확장 분석",
      items: [
        { label: "15 Drug Discovery", href: "15_drug_discovery.html", note: "타깃-화합물 우선순위와 약물 후보" },
        { label: "16 Quantum", href: "16_quantum.html", note: "양자 알고리즘 비교 벤치마크" },
        { label: "17 Biomarker Insights", href: "17_biomarker_insights.html", note: "마커 중심 해석" },
        { label: "18 Pathway / Immune / Methylation", href: "18_pathway_immune_meth.html", note: "경로·면역·메틸레이션 통합" },
        { label: "19a CMap Network", href: "19a_cmap_network.html", note: "약물 리버설과 PPI 네트워크" },
        { label: "20 Panel Size", href: "20_panel_size.html", note: "패널 크기 sweep k=4..60 + 3개 feature selection" },
        { label: "20a Meta Power", href: "20a_meta_power.html", note: "메타분석과 통계적 파워" },
        { label: "21 Three-class Classifier", href: "21_three_class.html", note: "BRAF_like / RAS_like / dediff 3-class + SHAP" },
        { label: "21a Survival", href: "21a_survival.html", note: "생존·예후 모델" },
        { label: "23 Bethesda Sim", href: "23_bethesda_sim.html", note: "Bethesda III/IV prevalence simulation" },
        { label: "24 Multimodal", href: "24_multimodal.html", note: "methylation-only multimodal (fusion·SCNA 미보유)" },
        { label: "25 External PRJEB11591", href: "25_external_prjeb11591.html", note: "external retry 상태 점검" },
        { label: "26 Fusion Anchor", href: "26_fusion_anchor.html", note: "fusion-aware anchor" },
        { label: "27 Pooled GPL570", href: "27_pooled_gpl570.html", note: "pooled validation" },
        { label: "28 Cross-cohort Summary", href: "28_cross_cohort_summary.html", note: "코호트 간 요약 비교" }
      ]
    },
    {
      title: "감사 / 강건성 / 메서드",
      desc: "누수, batch, 재현성, 보정 실패를 솔직하게 공개하는 페이지",
      items: [
        { label: "v3 Reviewer Guide", href: "v3_reviewer_guide.html", note: "리뷰어용 동선 가이드" },
        { label: "19 Honesty Audit", href: "19_honesty_audit.html", note: "v3 honesty audit" },
        { label: "29 ComBat Rescue", href: "29_combat_rescue.html", note: "batch correction rescue 시도" },
        { label: "30 Bethesda Clinical", href: "30_bethesda_clinical.html", note: "임상적 operating point 해석" },
        { label: "31 PRJEB11591 Status", href: "31_prjeb11591_status.html", note: "외부 코호트 재시도 현황" },
        { label: "32 Honest Audit", href: "32_honest_audit.html", note: "v4 honest audit" },
        { label: "33 v4 Synthesis", href: "33_v4_synthesis.html", note: "v4 요약 정리" },
        { label: "34 v5 Non-linear Correction", href: "34_v5_nonlinear_correction.html", note: "비선형 batch correction ablation" },
        { label: "35 v5 Adversarial Method", href: "35_v5_adversarial_method.html", note: "DANN de-confounding" },
        { label: "36 v5 Cross-cancer", href: "36_v5_cross_cancer.html", note: "cross-cancer generalization 점검" },
        { label: "37 v5 Robust Targets", href: "37_v5_robust_targets.html", note: "batch-robust target shortlist" },
        { label: "38 v5 Synthesis", href: "38_v5_synthesis.html", note: "v5 전체 synthesis" },
        { label: "39 v5.1 Real-Data LODO", href: "39_v5p1_real_data.html", note: "v5.1 real-data cross-cancer DIAL (LODO)" },
        { label: "v8 Statgen Supplement", href: "v8_statgen_supplement.html", note: "methodological robustness supplement" },
        { label: "v10 AAAI Technical", href: "v10_aaai_technical.html", note: "technical supplementary framing" }
      ]
    },
    {
      title: "플랫폼 / 제품화 뷰",
      desc: "기술을 제품 관점으로 재구성한 THYRAI 레이어",
      items: [
        { label: "OMEGA-THCA", href: "platform_omega.html", note: "biomarker discovery layer" },
        { label: "CHEM-THCA", href: "platform_chem.html", note: "drug discovery layer" },
        { label: "Chem Showcase", href: "platform_chem_showcase.html", note: "chemistry showcase" },
        { label: "PANEL-THCA", href: "platform_panel.html", note: "panel design layer" },
        { label: "CLINIC-THCA", href: "platform_clinic.html", note: "clinical framing layer" },
        { label: "Investor View", href: "view_investor.html", note: "투자자용 요약" },
        { label: "Researcher View", href: "view_researcher.html", note: "연구자용 요약" },
        { label: "v6 scRNA Foundation", href: "v6_scrna_foundation.html", note: "scGPT·Geneformer·GEARS·CellRank 기반 단일세포 해석" },
        { label: "v7 Drug Repurposing", href: "v7_drug_repurposing.html", note: "drug repurposing narrative" },
        { label: "v9 Interactive DIAL", href: "v9_interactive_dial.html", note: "interactive DIAL dashboard" },
        { label: "v11 Novel Pathways", href: "v11_novel_pathways.html", note: "novel pathway circuits" },
        { label: "v13 Drug Discovery", href: "v13_drug_discovery.html", note: "v13 통합 drug discovery (8 타깃 × 10 모달리티, TROP2 ADC)" }
      ]
    },
    {
      title: "Legacy 페이지 (이전 버전 / a 버전으로 대체됨)",
      desc: "신규 a 버전 또는 후속 페이지로 대체된 이전 분석. 변화 이력 추적용으로 보존",
      items: [
        { label: "19 CMap Network (legacy)", href: "19_cmap_network.html", note: "19a로 대체됨 — 초기 CMap·PPI" },
        { label: "20 Meta Power (legacy)", href: "20_meta_power.html", note: "20a로 대체됨 — 초기 메타분석/파워" },
        { label: "21 Survival (legacy)", href: "21_survival.html", note: "21a로 대체됨 — KM/Cox 초기 버전" },
        { label: "22 Survival (v3)", href: "22_survival.html", note: "TCGA OS event 적음, 해석 주의 — v3 버전" }
      ]
    },
    {
      title: "아카이브 / 조직",
      desc: "과거 분석, 연구 기록, 팀/외부 커뮤니케이션",
      items: [
        { label: "Archive Index", href: "v4a_archive_index.html", note: "보관된 세부 페이지 허브" },
        { label: "v4a Biomarkers", href: "v4a_biomarkers.html", note: "archive biomarkers" },
        { label: "v4a Data Integration", href: "v4a_data_integration.html", note: "archive integration" },
        { label: "v4a Drug Discovery", href: "v4a_drug_discovery.html", note: "archive drug discovery" },
        { label: "v4a ML Performance", href: "v4a_ml_performance.html", note: "archive ML performance" },
        { label: "v4a Panels", href: "v4a_panels.html", note: "archive panels" },
        { label: "v4a Pathway / Immune", href: "v4a_pathway_immune.html", note: "archive pathway" },
        { label: "v4a Quantum", href: "v4a_quantum.html", note: "archive quantum" },
        { label: "v4a Limitations", href: "v4a_limitations.html", note: "archive limitations" },
        { label: "Press", href: "press.html", note: "press / public communication" },
        { label: "Team", href: "team.html", note: "팀 소개" },
        { label: "Glossary", href: "99_glossary.html", note: "용어집" }
      ]
    }
  ];

  function isPageContext() {
    return /\/pages\/[^/]+$/.test(window.location.pathname);
  }

  function normalizeHref(href) {
    if (!href) return href;
    if (/^(?:[a-z]+:|#|\/\/|\/)/i.test(href)) return href;
    if (href === "../index.html") return isPageContext() ? "../index.html" : "index.html";
    return isPageContext() ? href : "pages/" + href;
  }

  function createEl(tag, cls, html) {
    var el = document.createElement(tag);
    if (cls) el.className = cls;
    if (html != null) el.innerHTML = html;
    return el;
  }

  function buildCard(group) {
    var card = createEl("section", "site-map__group");
    var head = createEl("div", "site-map__group-head");
    head.innerHTML = "<div><h3>" + group.title + "</h3><p>" + (group.desc || "") + "</p></div><span>" + group.items.length + " pages</span>";
    card.appendChild(head);

    var list = createEl("div", "site-map__list");
    group.items.forEach(function (item) {
      var a = createEl("a", "site-map__item");
      a.href = normalizeHref(item.href);
      a.setAttribute("data-label", (item.label + " " + (item.note || "") + " " + group.title + " " + (group.desc || "")).toLowerCase());
      a.innerHTML = "<strong>" + item.label + "</strong><span>" + (item.note || "") + "</span>";
      list.appendChild(a);
    });
    card.appendChild(list);
    return card;
  }

  function init() {
    if (!document.body) return;
    if (document.querySelector(".site-map-fab")) return;

    var fab = createEl("button", "site-map-fab", "전체 페이지 맵");
    fab.type = "button";
    fab.setAttribute("aria-expanded", "false");
    fab.setAttribute("aria-controls", "site-map-overlay");

    var dock = createEl("button", "topnav__map", "전체 메뉴맵");
    dock.type = "button";
    dock.setAttribute("aria-expanded", "false");

    var overlay = createEl("div", "site-map-overlay");
    overlay.id = "site-map-overlay";
    overlay.setAttribute("aria-hidden", "true");

    var panel = createEl("div", "site-map-panel");
    panel.innerHTML =
      '<div class="site-map__top">' +
      '<div><span class="site-map__eyebrow">THYRAI navigation</span><h2>전체 페이지 맵</h2><p>현재 연구 핵심, 상세 분석, 감사 로그, 플랫폼 페이지, 아카이브까지 전부 검색하고 바로 이동할 수 있습니다.</p><div class="site-map__legend"><span>데이터셋 탐색</span><span>EDA</span><span>ML</span><span>패널 비교</span><span>약물 발굴</span><span>감사/재현성</span></div></div>' +
      '<div class="site-map__actions"><input id="site-map-search" class="site-map__search" type="search" placeholder="페이지 이름 / 기능 / 키워드 검색, 예: EDA / biomarker / audit"><button type="button" class="site-map__close" aria-label="Close">닫기</button></div>' +
      "</div>";

    var grid = createEl("div", "site-map__grid");
    MAP.forEach(function (group) { grid.appendChild(buildCard(group)); });
    panel.appendChild(grid);
    overlay.appendChild(panel);

    document.body.appendChild(fab);
    document.body.appendChild(overlay);

    var navInner = document.querySelector(".topnav__inner") || document.querySelector(".topnav-inner");
    var navCta = document.querySelector(".topnav__cta") || document.querySelector(".nav-actions");
    if (navInner) {
      if (navCta && navCta.parentNode === navInner) navInner.insertBefore(dock, navCta);
      else navInner.appendChild(dock);
    }

    var search = panel.querySelector("#site-map-search");
    var closeBtn = panel.querySelector(".site-map__close");
    var groups = Array.prototype.slice.call(panel.querySelectorAll(".site-map__group"));

    function applyFilter() {
      var q = (search.value || "").trim().toLowerCase();
      groups.forEach(function (group) {
        var visible = 0;
        Array.prototype.slice.call(group.querySelectorAll(".site-map__item")).forEach(function (item) {
          var match = !q || item.getAttribute("data-label").indexOf(q) !== -1;
          item.hidden = !match;
          if (match) visible += 1;
        });
        group.hidden = visible === 0;
      });
    }

    function open() {
      overlay.classList.add("is-open");
      overlay.setAttribute("aria-hidden", "false");
      fab.setAttribute("aria-expanded", "true");
      dock.setAttribute("aria-expanded", "true");
      document.body.classList.add("site-map-open");
      window.setTimeout(function () { search.focus(); }, 30);
    }

    function close() {
      overlay.classList.remove("is-open");
      overlay.setAttribute("aria-hidden", "true");
      fab.setAttribute("aria-expanded", "false");
      dock.setAttribute("aria-expanded", "false");
      document.body.classList.remove("site-map-open");
    }

    fab.addEventListener("click", function () {
      overlay.classList.contains("is-open") ? close() : open();
    });
    dock.addEventListener("click", function () {
      overlay.classList.contains("is-open") ? close() : open();
    });
    closeBtn.addEventListener("click", close);
    overlay.addEventListener("click", function (e) {
      if (e.target === overlay) close();
    });
    search.addEventListener("input", applyFilter);
    panel.addEventListener("click", function (e) {
      var a = e.target.closest(".site-map__item");
      if (a) close();
    });
    document.addEventListener("keydown", function (e) {
      if ((e.metaKey || e.ctrlKey) && e.shiftKey && e.key.toLowerCase() === "m") {
        e.preventDefault();
        overlay.classList.contains("is-open") ? close() : open();
      }
      if (e.key === "Escape" && overlay.classList.contains("is-open")) close();
    });
  }

  window.ThyraiSiteMap = { init: init };

  if (document.readyState !== "loading") init();
  else document.addEventListener("DOMContentLoaded", init);
})();
