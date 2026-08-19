import React from "react";

/* ================================================================
   FAQ — extensive Korean Q&A, native <details> for click-to-expand.
   Categories: biology · methodology · clinical · validation · IRB · framing.
   ================================================================ */

type Q = { q: string; a: React.ReactNode };
type Group = { title: string; tag: string; color: string; qs: Q[] };

const GROUPS: Group[] = [
  {
    title: "A. Biology & 8-gene panel",
    tag: "BIOLOGY",
    color: "#1E40AF",
    qs: [
      {
        q: "왜 하필 이 8 개 유전자인가?",
        a: (
          <>
            <p>
              한 줄: <b>biology-prioritized, literature-guided, functionally grouped, cross-cohort validated, parsimony-driven.</b>
            </p>
            <ol>
              <li>문헌 기반 RAI uptake / thyroid differentiation 후보 풀 수집 (TIERA67 67 gene).</li>
              <li>5 기능 범주로 분류 (uptake · organification · TSH signaling · lineage TF · metabolism).</li>
              <li>TCGA + Lee + K2 + Landa + Mun + GPL570 4-cohort 등에서 재현성 확인.</li>
              <li>Redundancy / collinearity 검토 후 parsimony 적용 → 8 gene.</li>
              <li>16-gene TDS (Yoo 2016) 와 ΔAUC NS, random 8-gene 대비 cluster recovery 압도, 4-gene subset 대비 robust.</li>
            </ol>
            <p>
              결정적 reverse-causality 차단: 8-gene 은 Yoo 2016 PLOS Genet 의 RAI biology prior 만으로 선정되었고 outcome-driven selection 이
              들어가지 않음. Landa 2016 ATC silenced gene list 와 5 / 8 overlap 이지만 panel 정의 후 알게 된 사실 — convergent biology.
            </p>
          </>
        )
      },
      {
        q: "Effector 5 vs TF 3 — 왜 두 모듈로 나누는가?",
        a: (
          <>
            <p>
              <b>Effector 5</b> (TG · TPO · TSHR · SLC5A5 / NIS · DIO1) 는 실제 RAI 흡수 / 호르몬 합성 도구이고,
              <b> TF 3</b> (PAX8 · NKX2-1 / TTF-1 · FOXE1 / TTF-2) 는 갑상선 lineage master regulator 다.
            </p>
            <p>
              TF 3 가 꺼지면 effector 5 가 줄줄이 꺼지는 위계가 발달생물학적으로 알려져 있음. DM1 에서는 두 모듈이 동시 silencing — 즉
              단일 effector 손실이 아닌 lineage program 전체의 침묵.  TCGA 분석에서 TF 3 와 effector 5 의 ARI 가 cohort 별로 다르게 dominant
              (TCGA = effector-dominant ARI 0.876, K2 = TF-dominant ARI 0.880) — 두 모듈이 독립적으로 작동.
            </p>
          </>
        )
      },
      {
        q: "DM1 / DM2 는 기존 임상 분류와 무엇이 다른가?",
        a: (
          <>
            <p>
              DM1 / DM2 는 <b>본 분석으로 만든 molecular state classification</b> 이며 기존 임상 분류 (PTC / FTC / PDTC / ATC, BRAF-like / RAS-like) 와 다름.
            </p>
            <ul>
              <li>조직형 분류는 형태 기반, DM1 / DM2 는 8-gene 발현 기반.</li>
              <li>BRAF-like / RAS-like 는 driver 변이 기반, DM1 / DM2 는 driver-orthogonal.</li>
              <li>DM1 prevalence: BRAF V600E 0.7 % · RAS+ 96 % · BRAF/RAS-neg 49 % — driver 와 cross-cut.</li>
              <li>논문 rename 권장: <b>DM1 = Iodine-handling-low</b>, <b>DM2 = Iodine-handling-high</b>.</li>
            </ul>
          </>
        )
      },
      {
        q: "MAPK 활성과 DM1 의 관계는?",
        a: (
          <>
            <p>
              Per-driver-class promoter β: <b>BRAF V600E ≈ 0.37 · RET fusion ≈ 0.39 · RAS-mutant ≈ 0.27</b>.
              BRAF 와 RET fusion 은 모두 MAPK 활성이 매우 높은 driver 이고 (BRS / TDS metric 으로 알려진 것),
              RAS 는 그보다 약함. 따라서 methylation 은 driver 정체성이 아닌 <b>MAPK 활성에 비례</b>.
            </p>
            <p>
              이게 핵심: DM1 은 driver 종류와 무관하게 MAPK 활성이 충분히 강하거나, 또는 면역 / HT route 같은 다른 경로로
              분화 silencing 에 도달하는 모든 종양을 포착한다 — 두 axis (MAPK + HT) 모두 같은 endpoint (lineage silencing) 로 수렴.
            </p>
          </>
        )
      },
      {
        q: "단일세포 데이터에서 DM1 신호는 실제 thyrocyte intrinsic 인가?",
        a: (
          <>
            <p>
              Yes. Lu 2023 (GSE193581, n = 14,624) 에서 <b>KRT8 ∩ KRT19 ∩ EPCAM 필터링</b> 후에도 DM gradient 가 명확.
              KRT8 / KRT19 / EPCAM 모두 epithelial / thyrocyte marker — 따라서 stromal 또는 immune 세포가 신호를 만드는 것이 아님.
              Pu 2021 paired 6 환자에서도 per-patient r = 0.798–0.886 (Bonferroni p &lt; 10⁻¹⁰) 로 일관.
            </p>
          </>
        )
      }
    ]
  },
  {
    title: "B. Methodology & validation",
    tag: "METHODS",
    color: "#0E7490",
    qs: [
      {
        q: "Pan-genome ARI 0.92 가 panel artifact 차단의 직접 증거인가?",
        a: (
          <>
            <p>
              Yes. 비편향 pan-genome top-5000 MAD KMeans 결과가 8-gene 기반 partition 과 ARI = 0.92 — 즉 우리가 panel 을 알려주지 않고
              전체 transcriptome 만 보여도 같은 두 군이 나옴.  ARI ladder:
            </p>
            <ul>
              <li>Driver_anchor only → ARI = −0.007 (chance 미만)</li>
              <li>8-gene panel → 0.49 (compact lens)</li>
              <li>TIERA67 67-gene → 0.90</li>
              <li>pan-genome top-1000 → 0.91</li>
              <li><b>pan-genome top-5000 → 0.92</b></li>
            </ul>
            <p>
              즉 8-gene 은 panel 이 신호를 만든 게 아니라, 더 큰 transcriptomic axis 의 compact lens.
            </p>
          </>
        )
      },
      {
        q: "Within-cohort unsupervised 가 transfer label 과 다를 때 어떻게 해석하나?",
        a: (
          <>
            <p>
              K2 는 mini-index calibration mismatch 가 있음 — TCGA-trained transfer label 의 hard DM1 호출이 K2 에서 under-call (5.4 %).
              Within-K2 unsupervised k = 2 KMeans 는 <b>Cohen's d = 1.94</b>, DM1 prevalence 60.8 %.
              TCGA-trained <i>continuous</i> p(DM1) 은 within-K2 unsupervised label 을 정확히 순위 (AUC = 0.883).
            </p>
            <p>
              결론: 축은 실재하고 K2 에서 그대로 회복됨. Hard threshold 가 cross-ethnic calibration 에 민감 — 이건 <b>calibration drift</b> 이지
              axis failure 가 아님. Within-cohort unsupervised 가 ground-truth 같은 역할.
            </p>
          </>
        )
      },
      {
        q: "Master forest 14 entries 모두 같은 방향인가? 검정 방법은?",
        a: (
          <>
            <p>
              Yes. 14 entries × 11 distinct cohorts × 4 modality 모두 d ≥ 1.56.  최저: Pu 2021 single-cell d = 1.56.
              최고: Lee 2024 within-cohort KMeans (all-cohort) d = 5.93.  평균 = 2.81, median = 2.37.
            </p>
            <p>
              검정: 각 cohort 내 k = 2 KMeans 또는 known group contrast (ATC vs PTC, PDTC vs PTC) 의 Cohen's d (pooled SD).
              Pooled meta = Fisher-z (Spearman ρ entries) 또는 random-effects (Cohen's d entries) — DerSimonian-Laird + Han-Eskin RE2.
            </p>
          </>
        )
      },
      {
        q: "80 / 80 cell direction-consistent matrix 는 어떻게 만들었나?",
        a: (
          <>
            <p>
              8 panel 유전자 × 10 contrast × 4 cohort = 80 cell 의 Cohen's d 매트릭스.  각 cell 의 effect direction (DM1-low vs DM1-high 또는
              ATC vs PTC 등) 을 확인 → 80 / 80 모두 동일 방향.  TPO 가 최대 d = +3.53 (Lee 2024 ATC vs Normal).
              가장 약한 SLC5A5 (가장 noisy 한 effector) 도 모든 contrast 에서 같은 방향 유지.
            </p>
          </>
        )
      },
      {
        q: "HM450 promoter β 산출은? 어느 probe 를 썼나?",
        a: (
          <>
            <p>
              Illumina HM450 array · TCGA-THCA n = 503 paired samples.  Per-gene 의 TSS-proximal CpG probe 평균 β.
              각 panel 유전자 마다 TSS ± 200 bp 내의 모든 island / shore CpG probe 평균.  DIO1 / SLC5A5 / TPO 등은 single-probe-per-gene
              limitation 있음 (Limitations §3.4).  Mean β 의 cohort 평균을 DM1 / DM2 별로 계산.
            </p>
          </>
        )
      },
      {
        q: "Single-feature AUC ≈ 0.5 는 무엇을 의미하나?",
        a: (
          <>
            <p>
              BRAF · TERT · KRAS · NRAS · HRAS 각 driver 변이를 단일 feature 로 DM1 vs DM2 를 예측 시 AUC 가 chance level (≈ 0.5).
              즉 어떤 single driver 도 DM 분류를 설명하지 못함.  BRAF V600E 가 가장 높지만 AUC = 0.602 — 거의 우연 수준.
              이게 driver-orthogonality 의 가장 깨끗한 quantitative 증거.
            </p>
          </>
        )
      }
    ]
  },
  {
    title: "C. Clinical translation",
    tag: "CLINICAL",
    color: "#B45309",
    qs: [
      {
        q: "이 8-gene panel 로 RAI response 를 직접 예측할 수 있나?",
        a: (
          <>
            <p>
              <b>아니오 — 그렇게 주장하면 안 됨.</b>  공개 cohort 의 RAI response label 이 일관되지 않음.  안전한 표현:
            </p>
            <p style={{ background: "#ECFDF5", padding: "10px 14px", borderLeft: "3px solid #047857", borderRadius: "6px" }}>
              "이 8-gene panel 은 RAI 반응성을 직접 확정 예측한다기보다, <b>RAI 섭취 및 갑상선 분화 상태와 관련된 분자적 상태</b> 를 분류하는
              도구로 보는 것이 타당하다."
            </p>
            <p>
              가장 직접적인 RAI response 증거는 GSE151179 post-RAI refractory 종양이 transcriptionally DM1 state 와 일치한다는 점 (d ≈ −1.0).
              그러나 이건 association 이지 prospective prediction 검증이 아님.
            </p>
          </>
        )
      },
      {
        q: "왜 \"escalation 빠르게\" 보다 \"불필요한 RAI 회피\" 가 더 강한 메시지인가?",
        a: (
          <>
            <p>
              갑상선암 환자는 <b>장기 생존자</b> 가 대부분이다.  RAI 치료의 가장 큰 임상적 부담은:
            </p>
            <ul>
              <li>골수억제 · leukopenia · 백혈구 감소</li>
              <li>장기 (decades 단위) 누적 RAI 독성</li>
              <li>침샘 / 누관 손상</li>
              <li>이차암 (특히 백혈병 · 위암) 위험</li>
            </ul>
            <p>
              "고용량 RAI 를 안 들을 환자에게 반복하지 않는다" 는 <b>harm avoidance</b> — 즉 즉시 환자에게 직접 이익.  "더 빨리 systemic"
              은 <b>escalation</b> — aggressive 환자에 한정되고 비교군 정의가 어려움. 강민수 선생님 강조 포인트는 전자가 더 임상적으로
              설득력 있다는 것.
            </p>
          </>
        )
      },
      {
        q: "ATA 중간 위험군에서 DM1 의 임상 가치는?",
        a: (
          <>
            <p>
              ATA 2015 / 2025 의 위험군 분류는 low / intermediate / high.  RAI 결정의 <b>가장 큰 불확실성 zone 은 intermediate</b>:
              low 는 안 해도 되고, high 는 무조건 해야 하지만, intermediate 는 임상의의 판단에 따라 갈림.
            </p>
            <p>
              우리 TCGA 데이터에서 DM1 종양은 intermediate tier 에 과대표현 (DM1 92 / 140 = 66 % vs DM2 195 / 360 = 54 %).
              즉 DM1 / DM2 분층화의 임상 가치가 가장 높은 환자군이 바로 intermediate-risk.
              현재 ATA 기준이 형태학 / 분자 driver 만 사용하지 분화 axis 는 사용하지 않으므로 add-on value 있음.
            </p>
          </>
        )
      },
      {
        q: "Treatment-selection biomarker 라고 부르면 왜 위험한가?",
        a: (
          <>
            <p>
              "Treatment-selection biomarker" 는 regulatory 의미가 강함 — <b>특정 치료의 적응증 / 비적응증을 결정</b> 한다는 의미.
              그러려면 prospective randomized trial 에서 검증되어야 함 (예: CDx 승인).  우리는 그런 데이터 없음.
            </p>
            <p>
              안전한 표현: <b>"risk stratification axis"</b> — ATA 위험군 분류를 보완하는 분자 정보로,
              치료 결정의 <i>입력</i> 이지 치료 자체를 결정하지 않음.  Reviewer 가 over-claim 으로 공격할 가장 큰 위험.
            </p>
          </>
        )
      },
      {
        q: "FFPE 호환성 KS p = 0.44 의 의미?",
        a: (
          <>
            <p>
              FFPE (formalin-fixed paraffin-embedded) 는 임상 일상 표본 형식. Fresh-frozen (FF) 은 연구 표본 형식.  많은 RNA-based assay
              가 FFPE 에서 잘 안 됨 (RNA degradation).
            </p>
            <p>
              우리 8-gene score 의 FFPE (Lee 2024 GSE213647 n = 632) 분포 vs FF (TCGA n = 504) 분포의 Kolmogorov-Smirnov test p = 0.44.
              즉 통계적으로 구분 안 됨 → FFPE 표본으로도 검사 가능.  임상 deploy 의 가장 큰 기술적 장벽 제거.
            </p>
          </>
        )
      },
      {
        q: "Selpercatinib 후보군 \"1000 명 당 48 명\" 계산은?",
        a: (
          <>
            <p>
              TCGA RET-fusion 빈도 ≈ 6.5 % × DM1 capture 81.8 % × 적응증 일치 ≈ 0.9 ≈ <b>4.8 %</b> = 1000 PTC 당 48 명.
              이건 illustrative population estimate 이지 prospective trial 결과 아님.
            </p>
            <p>
              임상 의미: 전체 PTC 환자에게 RET fusion 검사를 routine 으로 하면 yield 가 낮지만 (≈ 6.5 %),
              DM1 으로 enrichment 한 후 검사하면 yield 가 ≈ 81.8 % — 즉 reflex 검사 알고리즘으로 비용 / 효율 큼.
            </p>
          </>
        )
      }
    ]
  },
  {
    title: "D. IHC, NGS panel, 임상 검증",
    tag: "DEPLOY",
    color: "#047857",
    qs: [
      {
        q: "분당병원 NGS panel 에 8 개 gene 이 있나? 시나리오별 대응은?",
        a: (
          <>
            <p>
              ★ 오늘 미팅에서 강민수 선생님께 확인할 첫 질문.  3 시나리오:
            </p>
            <ol>
              <li><b>모두 포함</b> → 직접 NGS-based validation 가능, retrospective 30 명 분석 즉시 시작.</li>
              <li><b>일부만 포함</b> → 4 / 5 / 6-gene reduced model 빌드. 16-gene 대비 ΔAUC 평가.
                  포함 후보 우선순위: TG &gt; TPO &gt; TSHR &gt; PAX8 &gt; NKX2-1 &gt; DIO1 &gt; SLC5A5 &gt; FOXE1.</li>
              <li><b>거의 안 포함</b> → 별도 RNA panel 또는 IHC pivot. 진단법 / IHC 기반 특허 검토.</li>
            </ol>
          </>
        )
      },
      {
        q: "IHC 로 8 개 marker 모두 staining 가능한가?",
        a: (
          <>
            <p>
              우선순위 ★★★ (임상 routine 가능):
            </p>
            <ul>
              <li><b>PAX8</b> — 갑상선 / 신장 / 난소 진단 routine marker, antibody 안정</li>
              <li><b>NKX2-1 / TTF-1</b> — 갑상선 / 폐 진단 routine, antibody 안정</li>
              <li><b>TG</b> — 갑상선 진단 routine, antibody 안정</li>
              <li><b>TPO</b> — 임상 적용 가능, antibody 안정</li>
              <li><b>NIS (SLC5A5)</b> — hMab 17.3 등 안정, prognostic 보고 있음</li>
            </ul>
            <p>
              ★★ (가능하나 표준화 부족):
            </p>
            <ul>
              <li><b>TSHR</b> — membrane staining, 표준화 부족</li>
              <li><b>FOXE1 / TTF-2</b> — antibody 가능하나 임상 적용 적음</li>
              <li><b>DIO1</b> — antibody 가능하나 표준화 부족</li>
            </ul>
          </>
        )
      },
      {
        q: "IHC scoring 은 어떻게 하나?",
        a: (
          <>
            <p>
              병리과와 협업해 결정해야 하나, 권장 옵션:
            </p>
            <ul>
              <li><b>H-score</b> (0 – 300): intensity × positive cell %. 가장 reproducible.</li>
              <li><b>Intensity</b> 0 / 1 / 2 / 3: 간단하나 cut-off 정의 필요.</li>
              <li><b>Positive cell %</b>: 일부 marker 에 적합.</li>
            </ul>
            <p>
              Tumor cell 만 score 할지 stromal 포함할지도 결정.  Whole-slide image AI 로 자동화하면 inter-observer variability 감소.
            </p>
          </>
        )
      },
      {
        q: "30 명 retrospective validation 의 가장 현실적인 endpoint?",
        a: (
          <>
            <p>
              우선순위 (chart retrieval 가능성 기준):
            </p>
            <ol>
              <li><b>Structural incomplete response</b> (ATA 2015 기준) — 가장 쉽게 정의 가능</li>
              <li><b>Cumulative RAI dose</b> — chart 에 정확히 기록되어 있음</li>
              <li><b>Stimulated Tg ≥ 10</b> at 6 mo / 12 mo</li>
              <li><b>Post-therapy whole-body scan</b> non-avid lesion 유무</li>
              <li><b>Recurrence / persistent disease</b> 1 년 / 3 년</li>
              <li><b>★ Marrow toxicity</b> — WBC / neutrophil / platelet 변화 (강민수 메시지)</li>
            </ol>
            <p>
              n = 30 은 univariate 검증만 가능. Multivariate 는 n ≥ 100 필요 — 후속 확장.
            </p>
          </>
        )
      },
      {
        q: "IRB 통과 시 가장 까다로운 부분은?",
        a: (
          <>
            <p>
              Retrospective chart review + NGS data + 익명화 → 일반적으로 expedited review 가능.  까다로운 부분:
            </p>
            <ul>
              <li>NGS data 와 chart 의 <b>de-identification + linkage</b> 방법론.</li>
              <li>RAI dose / 부작용 chart retrospective retrieval 가능성 confirm.</li>
              <li>외부 cohort 와의 데이터 공유 (있다면) → DTA 필요.</li>
              <li>IHC validation 단계에서 FFPE block 재사용 동의 면제 가능성.</li>
            </ul>
          </>
        )
      }
    ]
  },
  {
    title: "E. Reviewer 방어 + framing",
    tag: "DEFENSE",
    color: "#B91C1C",
    qs: [
      {
        q: "Pooled HR 2.53 이 MSK advanced cohort 에 의해 driven 되었다고 공격받으면?",
        a: (
          <>
            <p>
              <b>인정하고 분리 보고.</b>  TCGA primary cohort HR = 2.30 [0.77, 6.88] (n = 504, events = 16) 단독으론 underpowered.
              MSK HR = 2.67 [1.17, 6.10] (n = 117, events = 38) 이 통계적 유의성 견인.
            </p>
            <p>
              방어: I² = 0 % 로 두 cohort 가 완벽 일치 → MSK 가 outlier 아니고 같은 axis.  Primary PTC 단독 검정력은
              prospective cohort 에서만 가능 — Limitations §3.4 에 명시.  Random-effects pooled HR 이 fixed-effect 와 거의 동일.
            </p>
          </>
        )
      },
      {
        q: "\"이거 BRAF / RAS 일 뿐 아닌가?\" 라고 공격받으면?",
        a: (
          <>
            <p>
              세 가지 layer 로 답변:
            </p>
            <ol>
              <li><b>Single-feature AUC ≈ 0.5</b>: BRAF · TERT · KRAS · NRAS · HRAS 모든 driver 가 chance level.</li>
              <li><b>BRAF / RAS-negative dark matter 49 / 51 % split</b>: driver 가 없는 환자도 DM1 / DM2 가 절반씩 — 즉 axis 가 driver 외부에서도 작동.</li>
              <li><b>BRAF V600E ≈ RET fusion methylation</b>: driver 정체성이 아닌 MAPK 활성에 비례 — driver 종류가 아닌 활성 정도가 중요.</li>
            </ol>
          </>
        )
      },
      {
        q: "Reverse-causality (panel selection bias) 공격은?",
        a: (
          <>
            <p>
              가장 미묘한 공격: "8-gene 이 Landa silenced list 와 5/8 overlap 이니, 결국 advanced thyroid cancer signature 를 재발견한 것 아닌가?"
            </p>
            <p>
              <b>3-layer reverse-causality lock:</b>
            </p>
            <ol>
              <li><b>Panel origin</b>: 8-gene 은 Yoo 2016 PLOS Genet 의 RAI biology prior 만으로 선정.
                  Landa 2016 list 와 cross-reference 한 것은 panel 정의 <i>후</i>.</li>
              <li><b>Pan-genome convergence</b>: 8-gene 없이 top-5000 MAD clustering 만 해도 동일 partition (ARI 0.92).</li>
              <li><b>Independent design</b>: Yoo 2016 (RAI uptake biology) 와 Landa 2016 (advanced disease silenced list)
                  은 서로 다른 데이터 / 다른 목적으로 도출된 list — 5/8 overlap 은 convergent biology 의 직접 증거.</li>
            </ol>
          </>
        )
      },
      {
        q: "\"단일세포 figure 가 main 에 없는데 stromal 혼동 아닌가\" 공격은?",
        a: (
          <>
            <p>
              Extended Data 에 Lu 2023 single-cell figure 한 panel 배치.  KRT8 ∩ KRT19 ∩ EPCAM 필터 후에도 DM gradient 가 thyrocyte
              자체에서 명확 — stromal / immune confound 가 아닌 thyrocyte-intrinsic 신호.
            </p>
            <p>
              추가 layer: nu-SVR deconvolution 으로 8 cell type fraction 을 보정한 후에도 DM1 effect 가 유지됨 (Effect retention ratio 47 %).
              ED9 에 4-method deconvolution composite 보유.
            </p>
          </>
        )
      },
      {
        q: "Cross-cohort heterogeneity (I² = 72.9 %) 가 너무 높지 않나?",
        a: (
          <>
            <p>
              MAPK × Panel forest 의 I² = 72.9 % — 일반적으로 높은 수치지만 본 연구는 <b>예측된 heterogeneity</b>:
            </p>
            <ul>
              <li>GSE286332 (Korean PTC + HT, n = 18): HT-route 에서 MAPK 가 driver 가 아니므로 attenuated (ρ = −0.04 NS) — <b>two-axis model 의 사전 예측 일치</b>.</li>
              <li>GSE76039 (Landa PDTC + ATC, n = 37): 이미 dedifferentiated 끝에 있어 saturation (ρ = +0.125).</li>
              <li>TCGA + Lee primary PTC 2 cohort 가 전체 신호 견인 (n = 1,204, 둘 다 p &lt; 10⁻¹²).</li>
            </ul>
            <p>
              따라서 I² 은 noise 가 아니라 <b>모델이 예측한 axis-mixing</b> — 오히려 정합성 증거.
            </p>
          </>
        )
      },
      {
        q: "전향 검증 데이터 없이 NC 같은 top journal 에 갈 수 있나?",
        a: (
          <>
            <p>
              현실적 평가:
            </p>
            <ul>
              <li><b>Nature Cancer reach</b> (≈ 30–40 %): 외부 검증 19 cohort + protein cross-modality + driver-orthogonal mechanism 으로 가능.
                  Perturbation 이 main 에 없는 게 약점.</li>
              <li><b>Cell Reports Medicine (안정)</b> (≈ 70–80 %): 같은 evidence 로 무난히 reach.</li>
              <li><b>JCI Insight / Genome Medicine</b> (≈ 80–90 %): 안전권.</li>
              <li><b>npj Precision Oncology fallback</b>: 95 %+.</li>
            </ul>
            <p>
              분당 prospective cohort (현재 0 %) 추가 시 NC 확률 50 – 65 % 로 상승. 이게 venue ladder 의 가장 큰 변수.
            </p>
          </>
        )
      },
      {
        q: "DM1 / DM2 이름을 바꿔야 한다고? 어떤 이름이 좋은가?",
        a: (
          <>
            <p>
              "DM1 / DM2" 는 내부 코드명 느낌이 강함.  논문 발표 시 권장:
            </p>
            <table style={{ width: "100%", fontSize: 13, borderCollapse: "collapse", marginTop: 8 }}>
              <thead>
                <tr style={{ background: "#F1F5F9" }}>
                  <th style={{ padding: 8, textAlign: "left" }}>후보</th>
                  <th style={{ padding: 8, textAlign: "left" }}>의미</th>
                  <th style={{ padding: 8 }}>안전도</th>
                </tr>
              </thead>
              <tbody>
                <tr><td style={{ padding: 8 }}>RAI-sensitive-like / RAI-refractory-like</td><td style={{ padding: 8 }}>직관적</td><td style={{ padding: 8, textAlign: "center" }}>★★</td></tr>
                <tr><td style={{ padding: 8 }}>Iodine-avid / Iodine-poor</td><td style={{ padding: 8 }}>임상 표현</td><td style={{ padding: 8, textAlign: "center" }}>★★</td></tr>
                <tr><td style={{ padding: 8 }}>Thyroid-differentiated / De-differentiated</td><td style={{ padding: 8 }}>분자 중립</td><td style={{ padding: 8, textAlign: "center" }}>★★★</td></tr>
                <tr style={{ background: "#ECFDF5" }}><td style={{ padding: 8 }}><b>Iodine-handling-high / -low</b></td><td style={{ padding: 8 }}>가장 안전</td><td style={{ padding: 8, textAlign: "center" }}><b>★★★★</b></td></tr>
              </tbody>
            </table>
            <p style={{ marginTop: 10 }}>
              <b>권장: Iodine-handling-low (DM1) / Iodine-handling-high (DM2)</b>.  본문에서 "previously termed DM1 / DM2" 로 cross-ref.
            </p>
          </>
        )
      }
    ]
  },
  {
    title: "F. Numbers & limitations",
    tag: "STATS",
    color: "#6B21A8",
    qs: [
      {
        q: "TCGA HR 2.30 [0.77, 6.88] 의 CI 가 1 을 가로지른다 — significant 아닌 것 아닌가?",
        a: (
          <>
            <p>
              <b>맞음.</b>  TCGA primary PTC 코호트는 events = 16 으로 underpowered.  단독으론 통계적으로 NS.
              MSK-IMPACT HR 2.67 [1.17, 6.10] (events = 38) 와 pooled (DerSimonian-Laird random-effects) 했을 때 비로소 통합 HR 2.53 [1.31, 4.89] 가 유의.
              I² = 0 % 가 두 cohort 의 정합성을 입증.
            </p>
            <p>
              논문 caption 에 명시: "primary-PTC component is underpowered alone; pooled estimate is driven by MSK-IMPACT advanced cohort."
            </p>
          </>
        )
      },
      {
        q: "TPO d = 2.30 (p = 1.9 × 10⁻¹⁸) 가 너무 좋아 보이는데 multiple testing 보정했나?",
        a: (
          <>
            <p>
              Yes.  HM450 분석에서 8 panel + 2 control 유전자 (DIO2, SLC26A4) = 10 gene 의 per-gene Mann-Whitney U.
              Bonferroni 보정 (multiplier 10) 적용해도 TPO p 는 여전히 2 × 10⁻¹⁷ 수준 — 그 정도로 큰 effect.
              FDR (Benjamini-Hochberg) 으론 q &lt; 10⁻¹⁷.
            </p>
          </>
        )
      },
      {
        q: "GSE151179 post-RAI 의 sample size 가 작다 (n = 17 post)?",
        a: (
          <>
            <p>
              Yes.  Pre n = 35, post n = 17.  effect size d ≈ −1.0 은 매우 강하지만 sample 작아 95 % CI 가 넓음.
              이 결과는 association evidence 로 사용 (transcriptomically aligned), prospective response prediction 아님.
              Limitations 에 명시.  Mun 2025 protein cohort (n = 336) 가 sample size 보강.
            </p>
          </>
        )
      },
      {
        q: "Permutation null p = 0.0002 의 의미?",
        a: (
          <>
            <p>
              우리 8-gene panel 의 DM1 vs DM2 classification AUC 를 label 을 5,000 번 random shuffle 한 null distribution 과 비교했을 때
              empirical p = 0.0002 (즉 5,000 permutation 중 1 회만 우리만큼 좋음).  Panel 자체 (column) 는 고정.
              → "AUC 가 행운인가?" 공격 차단.
            </p>
            <p>
              Thyroid-matched bio-prior null 은 더 엄격: 상위 사분위 thyroid 발현 유전자 풀 (n = 12,921) 에서 8 개 random panel 1,000 회 시행 →
              empirical p = 0.007 (상위 0.7 %).  즉 "thyroid 관련 어떤 panel 이나 작동하지 않나?" 공격에도 차단.
            </p>
          </>
        )
      },
      {
        q: "8 개 중 1 개씩 빼면 (knock-out n−1) 어떻게 되나?",
        a: (
          <>
            <p>
              n − 1 = 7 gene LOO sensitivity:  최저 Cohen's d = 1.54 (DIO1 제거 시).  나머지 7 개 모두 d &gt; 1.5 — 단일 유전자
              dominant 아님.  중복 redundancy 가 적절 (overcompressed 아님).
            </p>
          </>
        )
      },
      {
        q: "Cluster stability bootstrap 97.4 % 의 의미?",
        a: (
          <>
            <p>
              TCGA n = 504 에서 500 회 bootstrap 재clustering.  97.4 % 의 종양이 ≥ 90 % bootstrap iteration 에서 동일 cluster 에 할당.
              즉 partition 이 single-run artifact 아닌 매우 stable 한 두 군 구조.
            </p>
          </>
        )
      }
    ]
  }
];

interface FaqProps { onlyOpenFirst?: boolean; }

export const Faq: React.FC<FaqProps> = () => (
  <div className="faq">
    {GROUPS.map((g, gi) => (
      <div key={gi} className="faq-group">
        <div className="faq-group-head" style={{ borderColor: g.color }}>
          <span className="faq-tag" style={{ background: g.color }}>{g.tag}</span>
          <h4 style={{ color: g.color }}>{g.title}</h4>
          <span className="faq-count">{g.qs.length}</span>
        </div>
        {g.qs.map((qa, i) => (
          <details key={i} className="faq-item">
            <summary>
              <span className="faq-q-mark" style={{ color: g.color }}>Q{i+1}.</span>
              <span className="faq-q-text">{qa.q}</span>
              <span className="faq-arrow">▾</span>
            </summary>
            <div className="faq-a">{qa.a}</div>
          </details>
        ))}
      </div>
    ))}
    <div className="faq-foot">
      총 <b>{GROUPS.reduce((s, g) => s + g.qs.length, 0)}</b> 개 FAQ · 클릭해서 펼치기.  미팅 중 reviewer / 임상의 질문 대응용.
    </div>
  </div>
);
