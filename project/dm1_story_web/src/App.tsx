import React, { useState, useEffect } from "react";
import { StorySection } from "./components/StorySection";
import { FigureCard } from "./components/FigureCard";
import { AxisDiagram } from "./components/AxisDiagram";
import { WorkflowDiagram } from "./components/WorkflowDiagram";
import { DarkMatterDiagram } from "./components/DarkMatterDiagram";
import { DriverLayerDiagram } from "./components/DriverLayerDiagram";
import { DistillationFlow } from "./components/DistillationFlow";
import { RaiPathwayDiagram } from "./components/RaiPathwayDiagram";
import { MethylationCartoon } from "./components/MethylationCartoon";
import { AriLadderChart } from "./components/AriLadderChart";
import { AtaImpactDiagram } from "./components/AtaImpactDiagram";
import { ClinicalPathwayDiagram } from "./components/ClinicalPathwayDiagram";
import { SelectionFlowDiagram } from "./components/SelectionFlowDiagram";
import { EvidencePyramid } from "./components/EvidencePyramid";
import { ClaimSafetyLadder } from "./components/ClaimSafetyLadder";
import { ValidationRoadmap } from "./components/ValidationRoadmap";
import {
  ExecutiveSummary,
  StateComparisonCards,
  GenePanelTable,
  KeyNumbersTable,
  IRBChecklist,
  ReviewerAttacksCard,
  MeetingQuestions
} from "./components/MeetingSections";
import { References } from "./components/References";
import { Faq } from "./components/Faq";
import { VerificationLog } from "./components/VerificationLog";
import { GeneLiterature } from "./components/GeneLiterature";
import { IhcDossier } from "./components/IhcDossier";
import { GeneCorrelation } from "./components/GeneCorrelation";
import { GeneReduction } from "./components/GeneReduction";
import { BiologyFoundations } from "./components/BiologyFoundations";
import { DeepAnalysis } from "./components/DeepAnalysis";
import { ReviewerQA } from "./components/ReviewerQA";
import { IhcHandoff } from "./components/IhcHandoff";
import { PresubInquiry } from "./components/PresubInquiry";

interface LB { src: string; alt: string }
interface SvgLB { html: string; title: string }

const TOC = [
  { group: "OVERVIEW", items: [
    { id: "summary", label: "Executive summary" }
  ]},
  { group: "★ MAIN STORY  ·  5 hero figures", items: [
    { id: "pathway", label: "MAIN-1 · 임상 의사결정 경로 (Fig A)" },
    { id: "panel",   label: "MAIN-2 · 8-gene iodine machinery" },
    { id: "why8",    label: "MAIN-3 · Why these 8 genes? (Fig B)" },
    { id: "states",  label: "MAIN-4 · DM1 / DM2 state 정의" },
    { id: "roadmap", label: "MAIN-5 · SNUBH validation roadmap (Fig E)" }
  ]},
  { group: "★★★ NC 승격 · 5-way 심층분석", items: [
    { id: "deep-analysis", label: "★★★ A1-A5 · Head2head + Interaction + Redif + Multi-omics + SNUBH sim" }
  ]},
  { group: "◆ SUPPORTING  ·  4 figs", items: [
    { id: "discovery",  label: "SUPP-6 · TCGA discovery (pan-genome)" },
    { id: "ev-gallery", label: "SUPP-7 · 외부 검증 forest (14 cohort)" },
    { id: "biology",    label: "SUPP-8 · FFPE / proteo / sc support" },
    { id: "biofound",   label: "SUPP-9 · Biology foundation 10 논문" },
    { id: "reduction",  label: "SUPP-10 · TSO500-3 단독 KM caveat" }
  ]},
  { group: "▼ MEETING ACTION", items: [
    { id: "clinical", label: "임상 의미 + ATA tier" },
    { id: "numbers",  label: "Audit-locked 핵심 수치" },
    { id: "irb",      label: "IRB chart-review 변수" },
    { id: "questions",label: "오늘 미팅 질문" },
    { id: "attacks",  label: "Reviewer attack 방어" },
    { id: "reviewer-qa", label: "★ Reviewer Q&A 시뮬레이션 (35 개)" },
    { id: "ihc-handoff", label: "★ IHC 3-plex 병리과 handoff protocol" },
    { id: "presub",      label: "★ Nature Med / NC pre-submission inquiry" },
    { id: "final",    label: "60 초 한 줄 요약" }
  ]},
  { group: "▽ BACKUP (펼치기)", items: [
    { id: "backup-section", label: "Backup · driver / methylation / FAQ / refs" }
  ]}
];

type TierProps = { tier: "main"|"support"|"meeting"|"backup"; label: string; sub: string };
const TierBanner: React.FC<TierProps> = ({ tier, label, sub }) => (
  <div className={`tier-banner tier-banner-${tier}`}>
    <div className="tier-banner-label">{label}</div>
    <div className="tier-banner-sub">{sub}</div>
  </div>
);

const App: React.FC = () => {
  const [lb,    setLb]    = useState<LB | null>(null);
  const [svgLb, setSvgLb] = useState<SvgLB | null>(null);
  const open  = (src: string, alt: string) => setLb({ src, alt });
  const close = () => { setLb(null); setSvgLb(null); };

  // ESC to close
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => { if (e.key === "Escape") close(); };
    if (lb || svgLb) window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [lb, svgLb]);

  // Global click delegation — any <img> or SVG-containing .visual becomes clickable
  useEffect(() => {
    const handler = (e: MouseEvent) => {
      const t = e.target as HTMLElement | null;
      if (!t) return;
      // Skip clicks inside lightbox or on interactive elements
      if (t.closest(".lb-bd") || t.closest(".lb-close-btn")) return;
      if (t.closest("a") || t.closest("button") || t.closest("input") || t.closest("textarea") || t.closest("select") || t.closest("summary")) return;

      // IMG in main content
      if (t.tagName === "IMG") {
        const img = t as HTMLImageElement;
        const src = img.getAttribute("src") || "";
        const alt = img.getAttribute("alt") || "";
        if (src && !src.startsWith("data:")) {
          setLb({ src, alt });
          e.preventDefault();
        }
        return;
      }

      // SVG inside .visual → serialize with proper namespaces via XMLSerializer
      const visual = t.closest(".visual") as HTMLElement | null;
      if (visual) {
        const svg = visual.querySelector("svg");
        if (svg) {
          const labelEl = visual.querySelector(".visual-label");
          const title = labelEl ? (labelEl.textContent || "").trim() : "";
          let svgStr = "";
          try {
            svgStr = new XMLSerializer().serializeToString(svg);
          } catch {
            svgStr = svg.outerHTML;
          }
          setSvgLb({ html: svgStr, title });
          e.preventDefault();
        }
      }
    };
    document.addEventListener("click", handler);
    return () => document.removeEventListener("click", handler);
  }, []);

  // Body scroll lock when lightbox open
  useEffect(() => {
    if (lb || svgLb) {
      const prev = document.body.style.overflow;
      document.body.style.overflow = "hidden";
      return () => { document.body.style.overflow = prev; };
    }
  }, [lb, svgLb]);

  return (
    <div className="layout">

      {/* ═══════ STICKY SIDEBAR TOC ═══════ */}
      <aside className="sticky-toc">
        <div className="sticky-toc-header">DM1 / DM2 Dossier</div>
        {TOC.map(g => (
          <div key={g.group}>
            <div className="sticky-toc-divider">{g.group}</div>
            <ol>
              {g.items.map(it => (
                <li key={it.id}><a href={`#${it.id}`}>{it.label}</a></li>
              ))}
            </ol>
          </div>
        ))}
      </aside>

      <div className="wrap">

      {/* ═══════ COVER ═══════ */}
      <header className="cover">
        <div className="eyebrow">CLINICAL RESEARCH BRIEFING  ·  분당병원  ·  2026-06-25 오전 9 시 미팅 자료</div>
        <h1 className="title">
          RAI 불응성 갑상선암을 <span className="accent">미리 구분</span>하는<br/>
          <span className="accent">8-Gene Thyroid Differentiation Axis</span>
        </h1>
        <div className="subtitle sans">
          DM1 / DM2 molecular state 를 이용해 불필요한 고용량 RAI 를 줄이고, 고위험 환자의 치료 전환을 앞당기는 임상 의사결정 도구로 발전시키는 프로젝트
        </div>
        <div className="byline">
          Seungho Cook (국승호) · Yu Hyeong-won (유형원) · 강민수 (clinical) · 서울대학교 분당병원
        </div>

        <div className="hero-badges">
          <span className="hero-badge cohort">TCGA-THCA n = 504</span>
          <span className="hero-badge panel">8-gene panel</span>
          <span className="hero-badge external">19 external cohorts</span>
          <span className="hero-badge ffpe">FFPE-compatible · KS p = 0.44</span>
          <span className="hero-badge ihc">IHC translation 계획</span>
          <span className="hero-badge pending">SNUBH validation pending</span>
        </div>

        <div className="manuscript-links sans">
          <a href="./manuscripts/reviewer_first_nc_rebuild.html" target="_blank" rel="noreferrer">
            Reviewer-first figure rebuild
          </a>
          <a href="./manuscripts/nature_communications_dm1_article.pdf" target="_blank" rel="noreferrer">
            Nature Communications PDF
          </a>
          <a href="./manuscripts/nature_communications_dm1_article.html" target="_blank" rel="noreferrer">
            Web manuscript with figures
          </a>
          <a href="./manuscripts/NATURE_COMMUNICATIONS_FULL_DRAFT_2026_07_08.md" target="_blank" rel="noreferrer">
            Source draft
          </a>
        </div>

        <div className="primary-msg">
          <div className="text">
            핵심은 <b>"항암으로 빨리 보낸다"</b> 보다 먼저,
            <b> RAI 가 흡수되지 않을 환자에게 불필요한 고용량 RAI 를 반복하지 않게 하는 것</b> 이다.
            <br/>
            <span style={{ fontSize: 14.5, color: "#475569", fontWeight: 400 }}>
              갑상선암은 장기 생존 환자가 많아 누적 RAI 독성 (골수억제 · leukopenia · 장기 부작용) 회피의 임상적 의미가 크다 — 강민수 선생님 강조.
            </span>
          </div>
        </div>
      </header>

      {/* ═══════ EXECUTIVE SUMMARY (always visible) ═══════ */}
      <section className="story" id="summary">
        <div className="num">섹션 1 · Executive summary</div>
        <h2>1 분 안에 보는 프로젝트 전체</h2>
        <div className="lead sans">Problem → Solution → Clinical value 세 카드.  세 카드만 읽어도 프로젝트 전체 framing 파악 가능.</div>
        <ExecutiveSummary />
      </section>

      {/* ═══════════════════════════════════════════════════════════════
         ★ MAIN STORY  ·  5 HERO FIGURES
         RAI harm-avoidance message 의 직접 표상.  Meeting 의 backbone.
         ═══════════════════════════════════════════════════════════════ */}
      <TierBanner tier="main"
                  label="★ MAIN STORY  ·  5 hero figures"
                  sub="이 5 개만 보여도 'RAI harm avoidance' 메시지가 전달된다.  meeting flow 의 backbone." />

      {/* ─── MAIN-1 · Clinical harm-avoidance pathway (Fig A) ─── */}
      <section className="story tier-main" id="pathway">
        <div className="num tier-num">MAIN-1 · 임상 의사결정 경로</div>
        <h2>Figure A — Current pathway vs Proposed molecular triage</h2>
        <div className="lead sans">현재 reactive 경로 (RAI 시행 후 6–12 개월 지연) vs 제안하는 pre-emptive molecular triage (수술 단계에서 분자 분류).</div>
        <div className="visual">
          <div className="visual-label sans">Figure A · clinical pathway comparison — RAI harm-avoidance pathway</div>
          <ClinicalPathwayDiagram />
        </div>
      </section>

      {/* ─── MAIN-2 · 8-gene iodine-handling machinery ─── */}
      <section className="story tier-main" id="panel">
        <div className="num tier-num">MAIN-2 · 8-gene iodine machinery</div>
        <h2>8 개 유전자 — 갑상선 iodine-handling cell machinery 의 직접 표상</h2>
        <div className="lead sans">분석 / 코드는 공식 symbol, 발표 / 병리 설명에는 alias 병기. Effector 5 (RAI 흡수) + TF 3 (분화 master) 두 모듈 — TSHR → cAMP → NIS → TPO → TG iodination 회로 그 자체.</div>
        <GenePanelTable />
        <div className="visual">
          <div className="visual-label sans">RAI 흡수 핵심 경로 — 8 panel genes 가 이 경로 그 자체</div>
          <RaiPathwayDiagram />
        </div>
      </section>

      {/* ─── MAIN-3 · Why these 8 genes selection flow (Fig B) ─── */}
      <section className="story tier-main" id="why8">
        <div className="num tier-num">MAIN-3 · Why these 8?</div>
        <h2>Figure B — 8 개 gene 은 갑자기 나온 것이 아니다</h2>
        <div className="lead sans">Reviewer 의 1 차 공격 차단.  literature pool → 기능 분류 → cohort 재현성 → parsimony → 최종 8 gene 의 7-step 흐름.</div>
        <div className="visual">
          <div className="visual-label sans">Figure B · 8-gene selection flow + reviewer defense</div>
          <SelectionFlowDiagram />
        </div>
      </section>

      {/* ─── MAIN-4 · DM1/DM2 state definition ─── */}
      <section className="story tier-main" id="states">
        <div className="num tier-num">MAIN-4 · DM1 / DM2 state</div>
        <h2>DM1 (Iodine-handling-low)  vs  DM2 (Iodine-handling-high)</h2>
        <div className="lead sans">DM1 / DM2 는 기존 임상 분류가 아닌 본 분석으로 만든 molecular state.  논문에서는 Iodine-handling-low / -high 로 rename 권장.</div>
        <StateComparisonCards />
      </section>

      {/* ─── MAIN-5 · SNUBH validation roadmap (Fig E) ─── */}
      <section className="story tier-main" id="roadmap">
        <div className="num tier-num">MAIN-5 · SNUBH validation roadmap</div>
        <h2>Figure E — SNUBH NGS / RNA / IHC 검증 Roadmap (Now / Next / Later)</h2>
        <div className="lead sans">현재 완료된 것 · 이번 미팅 후 즉시 할 일 · 확장 / 임상 단계의 3-lane timeline.  분당 NGS routine 위에 RNA add-on + IHC validation 분기.</div>
        <div className="visual">
          <div className="visual-label sans">Figure E · validation roadmap</div>
          <ValidationRoadmap />
        </div>
        <div className="visual">
          <div className="visual-label sans">translational workflow · surgery → FFPE → 8-gene → DM call → RAI strategy</div>
          <WorkflowDiagram />
        </div>
      </section>

      {/* ═══════════════════════════════════════════════════════════════
         ★★★ NC 승격 · 5-way parallel deep analysis
         ═══════════════════════════════════════════════════════════════ */}
      <TierBanner tier="main"
                  label="★★★ NC 승격 · 5-way parallel deep analysis (A1-A5)"
                  sub="Multiprocessing 병렬 · Head-to-head + Predictive interaction + Redif dynamics + Multi-omics + SNUBH sim" />
      <section className="story tier-main" id="deep-analysis">
        <div className="num tier-num">DEEP · A1-A5</div>
        <h2>5-way parallel 심층분석 — story 를 discovery paper → treatment-selection biomarker 로 승격</h2>
        <div className="lead sans">
          A1 head-to-head · A2 <b>DM1 × BRAF interaction p = 0.022 ★★★</b> · A3 GSE151179 redif dynamics · A4 multi-omics 5-way convergence · A5 SNUBH n=200 power = 90%.
          모든 5 분석 <b>CPU multiprocessing 병렬 실행</b>.  각 섹션에 그림 + 표 + 아주 자세한 한국어 해설 + caveat.
        </div>
        <DeepAnalysis />
      </section>

      {/* ═══════════════════════════════════════════════════════════════
         ◆ SUPPORTING EVIDENCE  ·  4 figs
         메인 스토리를 직접 받쳐주는 정량 / 외부 증거.
         ═══════════════════════════════════════════════════════════════ */}
      <TierBanner tier="support"
                  label="◆ SUPPORTING EVIDENCE  ·  4 figs"
                  sub="TCGA discovery + 외부 검증 forest + 생물학 회로 + survival/통계 caveat" />

      {/* ─── SUPP-6 · TCGA discovery + pan-genome reproducibility ─── */}
      <section className="story tier-support" id="discovery">
        <div className="num tier-num">SUPP-6 · TCGA discovery</div>
        <h2>DM1 분화 축은 TCGA 에서 명확히 분리되고, pan-genome 으로도 재현된다</h2>
        <div className="lead sans">TCGA-THCA n = 504 → 8-gene KMeans → DM1 28.4 % / DM2 71.6 %.  비편향 top-5000 MAD 로도 ARI 0.92 — panel artifact 아님.</div>
        <p>
          우리는 갑상선 분화 / RAI 흡수 생물학에 기반한 67 개 후보 유전자 풀 (TIERA67) 에서 <b>8 개 유전자 축</b> 을 도출했습니다.
          TCGA-THCA primary 종양 504 개에서 이 8 개로 KMeans (k = 2) 를 돌리면 두 군이 명확히 분리됩니다.
        </p>
        <p>
          이게 panel artifact 가 아닌 결정적 증거: panel 정의 후 비편향 pan-genome top-5000 MAD clustering 으로 동일 분할이 회복됩니다
          (<b>ARI = 0.92</b>). 8 gene 은 더 큰 transcriptomic axis 를 측정하는 <em>lens</em> 입니다.
        </p>
        <p>
          단일 driver 변이는 DM 분류 능력 거의 없음 — BRAF · TERT · KRAS · NRAS · HRAS 모든 driver 의 single-feature AUC ≈ 0.5.
        </p>
        <div className="visual">
          <div className="visual-label sans">distillation flow · broad transcriptome → TIERA67 → 8-gene readout</div>
          <DistillationFlow />
        </div>
        <div className="visual">
          <div className="visual-label sans">pan-genome ARI ladder — panel-induced artifact 차단</div>
          <AriLadderChart />
        </div>
      </section>

      {/* ─── SUPP-7 · External validation gallery ─── */}
      <section className="story tier-support" id="ev-gallery">
        <div className="num tier-num">SUPP-7 · 외부 검증</div>
        <h2>19 개 외부 코호트 · 7 modality · 4 인종 — DM1 축의 재현성</h2>
        <div className="lead sans">Reviewer 의 가장 흔한 공격 ("TCGA cherry-pick 아니냐") 의 1 차 차단선.  Master forest 평균 d = 2.81, 80/80 cell direction-consistent.</div>

        <div className="fig-grid gallery">
          <FigureCard
            title="EV-1 · Master cross-cohort forest"
            body="14 entries × 11 cohorts · 모든 entry d ≥ 1.56. Lee 2024 d = 5.93, Landa d = 4.08, Mun protein d = 2.86 등."
            src="./figures/EV_master_forest_annotated.png"
            alt="Master forest"
            stat="평균 d = 2.81 · median 2.37"
            onOpen={open}
          />
          <FigureCard
            title="EV-2 · Per-gene × cohort matrix"
            body="8 panel × 10 contrast × 4 cohort = 80 cell direction-consistent. 단일 유전자 의존성도 단일 cohort 의존성도 없음."
            src="./figures/EV_pergene_matrix_annotated.png"
            alt="Per-gene matrix"
            stat="80 / 80 direction-consistent"
            onOpen={open}
          />
          <FigureCard
            title="EV-3 · External Validation Atlas"
            body="15 + TCGA = 16 cohort composite. Modality / 인종 / 질병 단계 균등 배분."
            src="./figures/External_Validation_Atlas_annotated.png"
            alt="External validation atlas"
            stat="14 / 15 direction-consistent"
            onOpen={open}
          />
          <FigureCard
            title="EV-4 · GPL570 4-cohort"
            body="GSE33630 (n=105) ρ=−0.93 · GSE29265 ρ=−0.85 · GSE65144 ρ=−0.94 · GSE53157 ρ=−0.84."
            src="./figures/EV_gpl570_summary.png"
            alt="GPL570 4-cohort"
            stat="4 / 4 ρ ≤ −0.84"
            note="Western Affymetrix HG-U133 Plus 2.0"
            onOpen={open}
          />
          <FigureCard
            title="EV-5 · RAI / DM1 / nonoverlap × histology"
            body="4 cohort × 3 score × normal / PTC / ATC stratification. ATC 에서 RAI ↓ DM1 ↑ nonoverlap ↓."
            src="./figures/EV_rai_lineage_box.png"
            alt="External RAI lineage box"
            stat="trajectory consistent 4 / 4"
            onOpen={open}
          />
          <FigureCard
            title="EV-6 · Direction-consistency forest"
            body="7 sc + bulk 외부 cohort 의 Cohen's d / ρ vs DM1 방향. 4 cohort 에서 d > 1.0."
            src="./figures/EV_direction_consistency.png"
            alt="Direction consistency"
            stat="7 / 7 direction concordant"
            onOpen={open}
          />
          <FigureCard
            title="EV-7 · DM1 × nonoverlap scatter"
            body="8-panel 과 0 % overlap 인 독립 THYROID_NONOVERLAP panel 과의 scatter — reverse-causality 차단."
            src="./figures/EV_nonoverlap_scatter.png"
            alt="DM1 nonoverlap scatter"
            stat="ρ −0.84 ~ −0.94"
            onOpen={open}
          />
          <FigureCard
            title="EV-8 · Korean cohort applicability"
            body="K2 (n=260) + Lee (n=632) 한국 코호트. DM1 prevalence 37.8 % vs TCGA 28.4 %."
            src="./figures/EV_korean_applicability.png"
            alt="Korean applicability"
            stat="Korean DM1  37.8 %"
            onOpen={open}
          />
          <FigureCard
            title="EV-9 · K-Thyro 메타 (third-party)"
            body="Yu in-silico 제외 한국 갑상선 메타 — 자체 cohort 의존성 차단."
            src="./figures/EV_kthyro_meta.png"
            alt="K-Thyro meta"
            stat="direction reproduced"
            onOpen={open}
          />
          <FigureCard
            title="EV-10 · External dataset QC heatmap"
            body="모든 외부 cohort 8 / 8 panel 유전자 recovery + QC threshold 통과."
            src="./figures/EV_qc_heatmap.png"
            alt="External QC"
            stat="all cohorts QC pass"
            onOpen={open}
          />
          <FigureCard
            title="EV-11 · Atlas (raw)"
            body="16-card composite, raw 버전."
            src="./figures/External_Validation_Atlas.png"
            alt="Atlas raw"
            stat="6 modality × 4 인종"
            onOpen={open}
          />
          <FigureCard
            title="EV-12 · Master forest (raw)"
            body="14 entries × 11 cohorts, raw."
            src="./figures/EV_master_forest.png"
            alt="Master forest raw"
            stat="평균 d = 2.81"
            onOpen={open}
          />
          <FigureCard
            title="EV-13 · Per-gene matrix (raw)"
            body="80 cell direction matrix raw."
            src="./figures/EV_pergene_matrix.png"
            alt="Per-gene matrix raw"
            stat="80 / 80 direction"
            onOpen={open}
          />
          <FigureCard
            title="EV-14 · Future validation board"
            body="Bundang SNUBH 전향 cohort · K2 Bundang FFPE handoff · NC reviewer-readiness map."
            src="./figures/EV_future_validation.png"
            alt="Future validation"
            note="향후 검증 로드맵"
            onOpen={open}
          />
        </div>

        <div className="pullquote">
          외부 검증은 <b>4 line of defense</b> — effect-size magnitude (mean d 2.81), per-gene robustness (80 / 80 direction),
          platform 독립성 (GPL570 microarray), 인종 / lab 독립성 (Korean + Western + protein + third-party meta).
        </div>
      </section>

      {/* ─── SUPP-8 · FFPE / proteo / sc / RAI 회로 biology ─── */}
      <section className="story tier-support" id="biology">
        <div className="num tier-num">SUPP-8 · FFPE / proteo / sc 회로</div>
        <h2>DM1 의 본질은 thyroid differentiation loss — multi-modality 일관</h2>
        <div className="lead sans">DM1 에서 갑상선 transcription factors 와 RAI machinery 가 함께 꺼진다. FFPE · proteomics · single-cell 모두에서 같은 방향.</div>
        <p>
          DM1 종양에서는 <b>3 개 transcription factor</b> (PAX8 · NKX2-1 · FOXE1) 와 <b>5 개 effector</b>
          (TG · TPO · TSHR · SLC5A5 / NIS · DIO1) 가 함께 silencing 됩니다.  단순 clustering artifact 가 아닌
          <em> 갑상선 여포세포 분화 프로그램의 침묵</em> 그 자체입니다.
        </p>
        <p>
          이 8 개 gene 은 갑상선 분화의 핵심 회로 (TSHR → cAMP → NIS → TPO → TG iodination) 를 그대로 표상.
          따라서 DM1 은 <b>"RAI 가 잘 듣지 않을 상태"</b> 의 직접적인 분자 표현형입니다.
        </p>
        <p>
          정설 Thyroid Differentiation Score (TDS-16, Yoo 2016) 와 DM1 축은 동일 축에 정렬되고, Landa 2016 PDTC/ATC silenced gene list 와 <b>5/8 직접 overlap</b> — 독립 설계에서 도출된 convergent biology.
          Mun 2025 proteomics (n=336) 에서 7/7 protein direction sign concordant.
        </p>
        <div className="visual">
          <div className="visual-label sans">DM2 (분화 유지) ↔ DM1 (분화 손실) · thyroid follicular differentiation axis</div>
          <AxisDiagram />
        </div>
      </section>

      {/* ─── SUPP-9 · Biology foundation landmark papers ─── */}
      <section className="story tier-support" id="biofound">
        <div className="num tier-num">SUPP-9 · Biology foundation</div>
        <h2>DM1 축의 근간 — 10 개 landmark 논문 위에 서 있다</h2>
        <div className="lead sans">
          갑상선 여포세포 회로 (Carvalho 2017) · lineage TFs (Fernández 2015 · Nilsson 2017) · NIS 물리 (De la Vieja 2000) ·
          BRS/TDS convergent axis (TCGA 2014 · Landa 2016) · MAPK 침묵 reversibility (Chakravarty 2011 · Ho 2013 · Rothenberg 2015) ·
          RAI-refractory clinical burden (Schlumberger 2015).  각 카드 = key concept + PubMed / DOI 링크.
        </div>
        <BiologyFoundations />
      </section>

      {/* ─── SUPP-10 · TSO500-3 단독 KM caveat + reduction feasibility ─── */}
      <section className="story tier-support" id="reduction">
        <div className="num tier-num">SUPP-10 · TSO500-3 단독 KM caveat</div>
        <h2>★ 분당 TSO500 v2 (3-gene) 만으로 검증 가능?  + 8 → 4-6 gene 축소 타당성</h2>
        <div className="lead sans">
          TSO500 v2 가 우리 8 panel 중 3 개만 포함 (TSHR · PAX8 · NKX2-1) — 분당 코호트에 이 3 개만 측정 가능하면 어디까지 검증되나?
          모든 candidate subset 의 |Cohen's d| + 5-fold CV AUC + Spearman r vs Full-8 실데이터 (TCGA HM450 n = 496) 정량 검증.
          <b> 결론: TSO500-only 단독 = AUC 0.76 (성능 18 % p 손실, 단독 불가).  Compact 5 (TSO500 + TPO + DIO1) = AUC 0.940 (Full-8 동등).</b>
        </div>
        <GeneReduction />
      </section>

      {/* ═══════════════════════════════════════════════════════════════
         ▼ MEETING ACTION  ·  임상 의미 + 액션 항목
         ═══════════════════════════════════════════════════════════════ */}
      <TierBanner tier="meeting"
                  label="▼ MEETING ACTION"
                  sub="임상 의미 + 핵심 수치 + IRB · 질문 · reviewer 방어 · 60-sec 요약" />

      <section className="story" id="clinical">
        <div className="num">임상적 의미</div>
        <h2>DM1 은 RAI 실패와 공격적 임상 경과를 설명할 수 있다</h2>
        <div className="lead sans">ATA 중간 위험군에 과대표현 · post-RAI refractory 종양이 transcriptionally DM1 과 일치 · FFPE 호환 (KS p = 0.44).</div>
        <p>
          DM1 종양은 낮은 thyroid differentiation score 와 연결됩니다. ATA 2015 / 2025 위험군에서 DM1 은 <b>중간 위험군</b>
          (RAI 결정 불확실성이 가장 큰 구간) 에 과대표현 — 분자 분층화가 가장 임상적으로 가치 있는 환자군입니다.
        </p>
        <p>
          GSE151179 코호트에서 <b>RAI 치료 이후 refractory 상태로 진행한 종양</b> 은 transcriptionally DM1 state 와 일치합니다
          (Cohen's d ≈ −1.0, MW p ≈ 10⁻⁴).  DM1 은 RAI 실패를 사후에 보는 것이 아니라, <em>수술 조직에서 미리 감지 가능한 분자 상태</em>.
        </p>
        <p>
          후향 다중 코호트 통합 (TCGA + MSK-IMPACT) pooled OS HR = <b>2.53 [1.31, 4.89], I² = 0 %</b>.  단,
          이 수치는 후향 데이터 + advanced cohort 기여가 크므로 prospective 검증 필요.  Treatment-selection biomarker 가 아닌
          <b> risk stratification axis</b>.
        </p>
        <div className="visual">
          <div className="visual-label sans">ATA risk tier × DM 분포 · 분자 분층화가 가장 가치 있는 환자군</div>
          <AtaImpactDiagram />
        </div>
        <div className="fig-grid">
          <FigureCard
            title="Fig. 5 · 통합 생존 + portability"
            body="TCGA + MSK 통합 OS forest · multivariate Cox · FFPE vs FF KS p = 0.44."
            src="./figures/Fig5_annotated.png"
            alt="Figure 5"
            stat="pooled HR ≈ 2.5 · I² = 0 %"
            note="retrospective"
            onOpen={open}
          />
          <FigureCard
            title="Fig. 6 · 임상 reflex 경로"
            body="RNA → DM1 → fusion test → 치료 후보군. Post-RAI dediff 일치 (GSE151179)."
            src="./figures/Fig6_annotated.png"
            alt="Figure 6"
            stat="post-RAI d ≈ −1.0"
            note="prospective validation 전"
            onOpen={open}
          />
        </div>

        {/* ─── 한국어 해설 (임상의 이해 우선) ─── */}
        <div className="fig-explainer">
          <div className="fig-explainer-head">📊 Fig 5 쉬운 해설 — "DM1 환자는 어느 병원 어느 방식에서든 예후가 나쁘다"</div>
          <div className="fig-explainer-body">
            <div className="fig-explainer-part">
              <div className="fig-explainer-label">왼쪽 forest plot</div>
              <p>
                여러 병원의 갑상선암 코호트 (TCGA 미국 + MSK 미국 + 등 5개 소스) 에서 <b>DM1 그룹이 DM2 그룹보다 사망 위험 몇 배인지</b> 를 각각 계산.
                각 병원마다 사각형 = HR (hazard ratio) 추정치, 가로선 = 95 % 신뢰구간.  마지막 다이아몬드 = <b>모든 병원 통합한 pooled HR ≈ 2.5</b>.
              </p>
              <p>
                <b>임상 의미:</b> DM1 = 사망 위험 <b>2.5 배</b>.  "이 병원 데이터만 그런가?" 의심 차단 — 5 개 병원 모두 같은 방향.
                I² = 0 % → 병원간 이질성 없음 (한 병원의 편향 아님).
              </p>
            </div>
            <div className="fig-explainer-part">
              <div className="fig-explainer-label">오른쪽 FFPE vs Fresh-frozen portability</div>
              <p>
                <b>FFPE</b> = 병원에서 표준적으로 저장하는 파라핀 블록 (수술 후 곧바로 포름알데히드 고정).
                <b>Fresh-frozen</b> = 연구용 냉동 조직 (일반 병원은 안 함).  일반 병원 조직 = FFPE 뿐.
              </p>
              <p>
                본 figure 는 <b>같은 8-gene panel 을 FFPE 에서 재검정 시 fresh-frozen 결과와 같은가</b> 를 검증.
                Kolmogorov-Smirnov 통계 p = <b>0.44</b> → 두 방식의 결과 분포가 <b>통계적으로 구별 불가</b> = 같음.
              </p>
              <p>
                <b>임상 의미:</b> 분당병원의 <b>기존 파라핀 블록 (retrospective FFPE)</b> 에서도 이 panel 작동.  냉동 조직 새로 수집 불필요.
                <em>= IRB chart-review 로도 즉시 착수 가능.</em>
              </p>
            </div>
            <div className="fig-explainer-caveat">
              <b>주의:</b> 후향 (retrospective) 데이터.  advanced-disease 코호트 (MSK) 가 pooled HR 을 끌어올림.  전향 검증 (prospective) 은 아직.
              Treatment-selection biomarker (치료 결정) 로 과장하지 않음 — <b>risk stratification (환자 분류)</b> 축까지가 현재 evidence.
            </div>
          </div>
        </div>

        <div className="fig-explainer">
          <div className="fig-explainer-head">🔀 Fig 6 쉬운 해설 — "8-gene 결과가 나오면 다음에 뭘 하는가?"  reflex 임상 경로</div>
          <div className="fig-explainer-body">
            <div className="fig-explainer-part">
              <div className="fig-explainer-label">Reflex algorithm (단계별 흐름)</div>
              <p>
                <b>Step 1</b> — 수술 후 갑상선암 조직 (FFPE) 에서 8-gene RNA readout 시행 (또는 IHC 3-plex).<br/>
                <b>Step 2</b> — DM1 vs DM2 로 분류.<br/>
                <b>Step 3-A</b> · <b>DM2 (iodine-handling-high)</b>: 표준 RAI 프로토콜 유지, 저용량 검토.<br/>
                <b>Step 3-B</b> · <b>DM1 (iodine-handling-low)</b>: (i) BRAF / NTRK / RET fusion 재검 → (ii) target 약제 후보 조기 identify → (iii) 고용량 RAI 대신 <b>redifferentiation therapy (dabrafenib · selumetinib)</b> 또는 <b>TKI (lenvatinib · sorafenib)</b> 조기 sequencing 고려.
              </p>
            </div>
            <div className="fig-explainer-part">
              <div className="fig-explainer-label">GSE151179 external validation 요약 (Fig 6 왼쪽 박스)</div>
              <p>
                GSE151179 = 이스라엘 코호트 · 갑상선암 환자 39 명 · 이 중 4 명이 실제로 RAI 치료 <b>후</b> refractory (RAI 실패) 확진.
                이 4 명의 <b>수술 시점 (RAI 하기 전) 조직</b> 을 다시 조회 → 8-gene score 계산 → <b>Cohen's d ≈ −1.0, MW p ≈ 10⁻⁴</b> 로 refractory 미리 예측됨.
              </p>
              <p>
                <b>임상 의미:</b> DM1 판별은 사후 (RAI 실패 후) 가 아니라 <b>수술 시점 조직에서 미리 감지 가능</b>.  이것이 "reflex" 라는 이름의 뜻 — 수술 병리 나오는 즉시 자동 다음 단계.
              </p>
            </div>
            <div className="fig-explainer-part">
              <div className="fig-explainer-label">"reflex 경로" 가 임상적으로 왜 중요한가</div>
              <p>
                현재 표준 = RAI 시행 → 6-12 개월 뒤 재검사 → refractory 확인 → 다음 치료로 전환.  이 지연 동안 환자는 <b>불필요한 고용량 RAI 노출 + 병 진행</b> 을 감수.
              </p>
              <p>
                Reflex 경로 = 수술 병리 판독 즉시 (약 <b>1-2 주 안</b>) DM1 여부 판정 → 고위험군은 처음부터 <b>고용량 RAI 반복 회피 + 표적치료 sequencing 조기화</b>.
              </p>
            </div>
            <div className="fig-explainer-caveat">
              <b>주의:</b> DM1 이 곧 "RAI 실패" 는 아님.  DM1 은 RAI 흡수 회로 침묵 상태 = <b>실패 위험 higher</b> 의미.  개별 환자 치료 결정은 임상의 종합 판단.
              분당 prospective validation 이 이 reflex 경로의 임상 정당성 근거.
            </div>
          </div>
        </div>
      </section>

      <section className="story" id="numbers">
        <div className="num">핵심 수치 (audit-locked)</div>
        <h2>Audit-locked 핵심 수치 — 도메인별 한 표</h2>
        <div className="lead sans">17 개 핵심 수치 · Domain · Cohort · Result · Interpretation · Caveat 형식.  외부 공유 전 PI 검토 필요.</div>
        <KeyNumbersTable />
      </section>

      <section className="story" id="irb">
        <div className="num">IRB chart-review</div>
        <h2>IRB chart-review 변수 — 7 그룹</h2>
        <div className="lead sans">Demographics · Pathology · Molecular · Treatment · Biochemistry · Outcomes · <b>Toxicity (★)</b> — 강민수 선생님 메시지 반영.</div>
        <IRBChecklist />
      </section>

      <section className="story" id="questions">
        <div className="num">오늘 미팅 질문</div>
        <h2>오늘 미팅에서 반드시 확인할 질문</h2>
        <div className="lead sans">강민수 선생님께 (임상 / NGS / IRB) + 병리 협업 (IHC) 두 column.</div>
        <MeetingQuestions />
      </section>

      <section className="story" id="attacks">
        <div className="num">Reviewer 방어</div>
        <h2>Top 5 reviewer 공격 + 사전 답변</h2>
        <div className="lead sans">왜 8 개 · 정말 예측인가 · BRAF/RAS 아닌가 · platform-specific 아닌가 · 임상 적용 가능한가.</div>
        <ReviewerAttacksCard />
      </section>

      <section className="story" id="reviewer-qa">
        <div className="num">Reviewer Q&A 시뮬레이션</div>
        <h2>★ 35 개 예상 reviewer 질문 + 방어 답변 시뮬레이션</h2>
        <div className="lead sans">
          NC / Nature Medicine 편집자 + 3 명 reviewer 관점에서 <b>가장 aggressive 한 예상 질문</b> 을 pre-empt.
          5 카테고리 (Editor · Methods · Biology · Clinical · External · Framing) · 강력/중간/약함 방어 등급 표시.
        </div>
        <ReviewerQA />
      </section>

      <section className="story" id="ihc-handoff">
        <div className="num">IHC 3-plex Handoff Protocol</div>
        <h2>★ 병리과 즉시 실행 프로토콜 · Pilot 30 명 · IRB Fast-track template</h2>
        <div className="lead sans">
          분당병원 병리과가 <b>이미 routine 시행중인 3 항체</b> (TG · PAX8 · TTF-1) 만으로 BRAF+ subset PFI 분리 검증.
          §1 항체 사양 · §2 Scoring · §3 QC · §4 Cohort design · §5 IRB template · §6 SAP · §7 비용 · §8 Escalation.
        </div>
        <IhcHandoff />
      </section>

      <section className="story" id="presub">
        <div className="num">Nature Med / NC Presubmission Inquiry</div>
        <h2>★ Editorial pre-submission inquiry 패키지 · 200-word summary + cover letter</h2>
        <div className="lead sans">
          Nature Medicine primary → NC fallback → JCI Insight backup 3-tier journal ladder.
          200-word summary · Cover letter draft · Suggested reviewers · Novelty highlights · Timeline (미팅 후 3 개월).
        </div>
        <PresubInquiry />
      </section>

      <section className="story" id="final">
        <div className="num">60 초 요약</div>
        <h2>오늘 한 줄로 기억할 것</h2>
        <div className="final-slide">
          <div className="label">★ 60 초 요약</div>
          <div className="text">
            우리는 RAI response 를 과장해서 맞힌다고 주장하기보다,
            <b> thyroid differentiation / iodine-handling state 를 robust 하게 분류</b> 하고,
            이를 통해 <b>불필요한 고용량 RAI 를 줄이는 임상 의사결정 도구</b> 로 발전시키자는 방향이다.
          </div>
          <div className="signoff">
            <b>8 genes = lens</b>  ·  <b>DM1 / DM2 = discovery</b>  ·  <b>message = harm avoidance before escalation</b>
          </div>
        </div>
      </section>

      {/* ═══════════════════════════════════════════════════════════════
         ▽ BACKUP  ·  Reviewer 가 깊이 파고들 때만 열기
         RAI harm-avoidance 메시지에서 벗어나므로 평소엔 closed.
         ═══════════════════════════════════════════════════════════════ */}
      <TierBanner tier="backup"
                  label="▽ BACKUP  ·  깊이 들어가야 하는 reviewer 질문에만 펼치기"
                  sub="아래 9 개 섹션은 RAI harm avoidance 핵심 메시지에서 벗어남 — 평소에는 closed" />

      <details className="backup-details" id="backup-section">
        <summary>
          <span className="backup-summary-label">▽ Backup 9 개 deep-dive 섹션 펼치기</span>
          <span className="backup-summary-sub">driver mutation 상세 · 메틸화 회로 · evidence pyramid · claim safety · gene 문헌 · IHC vendor · methylation correlation · verification log · 38 FAQ · references</span>
        </summary>

        {/* ─── Driver mutation detail (mechanism) — BACKUP ─── */}
        <section className="story tier-backup" id="mechanism">
          <div className="num tier-num">BACKUP · driver mutation 상세</div>
          <h2>DM1 은 driver 범주를 가로지른다 — driver-orthogonal silencing</h2>
          <div className="lead sans">BRAF V600E 양성 DM1 prevalence 0.7 % · RAS-mutant 96.4 % · BRAF/RAS-음성 49 / 51 %.  하나의 driver 가 DM1 / DM2 를 결정하지 않는다.</div>
          <p>
            Methylation 분석은 DM1 에서 lineage gene promoter 의 <b>hypermethylation</b> 과 expression silencing 이 함께 나타남을 시사합니다.
            더 중요한 발견은: methylation 패턴이 driver 정체성이 아닌 <em>MAPK 활성</em> 에 비례한다는 것 — BRAF V600E (β = 0.37) ≈ RET fusion (0.39) ≫ RAS (0.27).
          </p>
          <div className="visual">
            <div className="visual-label sans">conceptual · canonical drivers → dark matter zone → DM1 axis</div>
            <DarkMatterDiagram />
          </div>
          <div className="visual">
            <div className="visual-label sans">two-layer · driver category × DM1 / DM2 state</div>
            <DriverLayerDiagram />
          </div>
          <div className="visual">
            <div className="visual-label sans">epigenetic mechanism cartoon · promoter 메틸화 → TF 차단 → 발현 OFF</div>
            <MethylationCartoon />
          </div>
          <div className="fig-grid">
            <FigureCard
              title="Fig. 3 · 후성유전 침묵 (실데이터)"
              body="HM450 promoter β heatmap (TCGA n = 503). DM1 종양에서 8 유전자 promoter 일관 메틸화 · TPO d = 2.30."
              src="./figures/Fig3_annotated.png"
              alt="Figure 3 — epigenetic silencing"
              stat="TPO d = 2.30  ·  p = 1.9 × 10⁻¹⁸"
              note="현재 분석 기준"
              onOpen={open}
            />
          </div>
        </section>

        {/* ─── Evidence pyramid — BACKUP ─── */}
        <section className="story tier-backup" id="pyramid">
          <div className="num tier-num">BACKUP · Evidence pyramid</div>
          <h2>Figure C — Biology 에서 임상 적용까지 누적 증거 피라미드</h2>
          <div className="lead sans">9 개 layer 의 누적 증거 — literature biology 부터 임상 적용 (SNUBH 검증 대기) 까지.</div>
          <div className="visual">
            <div className="visual-label sans">Figure C · evidence stack pyramid</div>
            <EvidencePyramid />
          </div>
        </section>

        {/* ─── Claim safety ladder — BACKUP ─── */}
        <section className="story tier-backup" id="safety">
          <div className="num tier-num">BACKUP · Claim safety</div>
          <h2>Figure D — 무엇을 주장하면 안 되고, 무엇을 주장해야 하는가</h2>
          <div className="lead sans">5 쌍의 위험한 표현 vs 안전한 표현.  Manuscript / 발표 / 미팅 모두에 적용되는 가드레일.</div>
          <div className="visual">
            <div className="visual-label sans">Figure D · claim safety ladder</div>
            <ClaimSafetyLadder />
          </div>
        </section>

        {/* ─── Gene literature dossier — BACKUP ─── */}
        <section className="story tier-backup" id="genelit">
          <div className="num tier-num">BACKUP · 유전자별 문헌</div>
          <h2>8 panel 유전자 — gene-by-gene literature dossier</h2>
          <div className="lead sans">각 유전자의 기능 · RAI / dedifferentiation 임상 문헌 · 선천성 갑상선기능저하 (CH) 유전학 · IHC clone · 한국인 cohort 적용.</div>
          <GeneLiterature />
        </section>

        {/* ─── IHC vendor — BACKUP ─── */}
        <section className="story tier-backup" id="ihc">
          <div className="num tier-num">BACKUP · IHC vendor</div>
          <h2>8 panel 유전자 IHC antibody · vendor · 기존 데이터</h2>
          <div className="lead sans">연구용 / 임상용 IHC 항체 가용성 (clone · vendor · catalog · IVD vs RUO) + 기존 갑상선암 IHC 문헌 + 분당병원 routine 보유 가능성.</div>
          <IhcDossier />
        </section>

        {/* ─── 8-gene correlation deep stat — BACKUP ─── */}
        <section className="story tier-backup" id="genecorr">
          <div className="num tier-num">BACKUP · 상관 분석</div>
          <h2>8 유전자 상관 · lineage · 독립성 — 실데이터 분석</h2>
          <div className="lead sans">TCGA HM450 + K2 RNA-seq 양 cohort 에서 8 유전자의 gene-gene Spearman 상관, 위계적 클러스터링, PCA 분해, 상관 네트워크 분석.</div>
          <GeneCorrelation />
        </section>

        {/* ─── Verification log — BACKUP ─── */}
        <section className="story tier-backup" id="verify">
          <div className="num tier-num">BACKUP · 수치 검증</div>
          <h2>Verification log — 27 개 핵심 claim 의 출처 / 상태 / audit date</h2>
          <div className="lead sans">모든 수치는 manuscript audit log 와 cross-check 됨.  Verified · Audit-locked · Needs PI check · Illustrative 4-tier 상태 표시.</div>
          <VerificationLog />
        </section>

        {/* ─── 38 FAQ — BACKUP ─── */}
        <section className="story tier-backup" id="faq">
          <div className="num tier-num">BACKUP · 38 FAQ</div>
          <h2>자주 묻는 질문 — 38 개 (클릭하여 펼치기)</h2>
          <div className="lead sans">미팅 중 reviewer / 임상의 / 병리과 질문 즉시 대응용. 6 카테고리 — biology · methodology · clinical · deploy · defense · stats.</div>
          <Faq />
        </section>

        {/* ─── References — BACKUP ─── */}
        <section className="story tier-backup" id="refs">
          <div className="num tier-num">BACKUP · 참고문헌</div>
          <h2>Reference 모음 — DOI · PMID · GEO accession 모두 새 탭 링크</h2>
          <div className="lead sans">5 카테고리 (foundational thyroid 문헌 · 임상 가이드라인 · validation cohorts · 통계/계산 methods · tools).</div>
          <References />
        </section>
      </details>

      {/* ═══════ LIMITATION ═══════ */}
      <div className="limit-box">
        <div className="label sans">주의 · Pre-validation limitation</div>
        <div className="text">
          본 발견은 <b>retrospective 다중 코호트 데이터에 기반한 후향 분석 단계</b> 이며, 다음 한계를 명시합니다:
          <ul>
            <li>SNUBH 또는 외부 prospective cohort 에서의 전향 검증 (prospective clinical validation) 전 단계입니다.</li>
            <li>DM1 은 risk stratification axis 이며, treatment-selection biomarker 로 과장하지 않습니다.</li>
            <li>Pooled OS HR 은 advanced-disease 코호트 (MSK-IMPACT) 기여가 큰 통합 추정치입니다.</li>
            <li>현 단계 수치는 manuscript audit-locked 분석 기준이며, 외부 공유 전 PI 검토가 필요합니다.</li>
            <li>IHC validation · wet-lab perturbation · spatial transcriptomics 는 후속 단계입니다.</li>
          </ul>
        </div>
      </div>

      <footer className="foot">
        Prepared for meeting · <b>2026-06-25 오전 9 시</b>  ·  Cook & Yu · 강민수 (clinical) · 서울대학교 분당병원<br/>
        Project: DM1/DM2 8-gene thyroid differentiation / iodine-handling axis  ·  Status: retrospective public-cohort evidence; SNUBH validation pending<br/>
        Master draft <code>project/manuscript_biorxiv_2026_05_20/paper1_biorxiv_kr.pdf</code> · Brief <code>project/manuscript_v8/DM1_MASTER_BRIEF_FOR_GPT_AND_MEETING_2026_06_25.md</code>
      </footer>

      {/* ── Image lightbox ── */}
      {lb && (
        <div className="lb-bd"
             onClick={(e) => { if (e.target === e.currentTarget) close(); }}
             role="dialog" aria-label="figure lightbox">
          <button className="lb-close-btn"
                  onClick={(e) => { e.stopPropagation(); close(); }}
                  aria-label="close">
            <span className="lb-x">×</span>
            <span className="lb-close-label">닫기 · ESC</span>
          </button>
          <div className="lb-hint">배경 클릭 · ESC · 우측 상단 버튼 모두 닫기</div>
          <img className="lb-img" src={lb.src} alt={lb.alt} onClick={(e) => e.stopPropagation()} />
          {lb.alt && <div className="lb-caption">{lb.alt}</div>}
        </div>
      )}
      {/* ── SVG diagram lightbox ── */}
      {svgLb && (
        <div className="lb-bd"
             onClick={(e) => { if (e.target === e.currentTarget) close(); }}
             role="dialog" aria-label="svg lightbox">
          <button className="lb-close-btn"
                  onClick={(e) => { e.stopPropagation(); close(); }}
                  aria-label="close">
            <span className="lb-x">×</span>
            <span className="lb-close-label">닫기 · ESC</span>
          </button>
          <div className="lb-hint">배경 클릭 · ESC · 우측 상단 버튼 모두 닫기</div>
          <div className="lb-svg-wrap"
               onClick={(e) => e.stopPropagation()}
               dangerouslySetInnerHTML={{__html: svgLb.html}} />
          {svgLb.title && <div className="lb-caption">{svgLb.title}</div>}
        </div>
      )}
      </div>
    </div>
  );
};

export default App;
