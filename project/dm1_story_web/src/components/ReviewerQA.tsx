import React, { useState } from "react";

/**
 * ReviewerQA — 35 예상 reviewer 질문 + 방어 답변 시뮬레이션.
 * NC / Nature Medicine 편집자 · 3 명 reviewer 관점에서 아주 aggressive 하게 시뮬레이트.
 * 각 답변에 (a) 핵심 방어 논리, (b) 근거 figure / 수치, (c) 방어 강도 등급.
 */

type QA = {
  q: string;
  a: string[];
  evidence?: string;
  strength: "strong" | "moderate" | "weak" | "escalate";
};

type Category = {
  key: string;
  title: string;
  subtitle: string;
  color: string;
  qas: QA[];
};

const CATEGORIES: Category[] = [
  {
    key: "editor",
    title: "Editor-level  ·  Fit / Novelty / Priority",
    subtitle: "편집자 desk-reject 방지",
    color: "#7C2D12",
    qas: [
      {
        q: "왜 이 논문이 Nature Communications / Nature Medicine 급인가?  기존 TCGA-THCA (Cell 2014), Landa 2016 (JCI) 있는데 우리가 새로 밝힌 게 뭔가?",
        a: [
          "본 논문의 핵심 novelty 는 3 가지 —",
          "① BRAF+ 서브셋 특이적 predictive interaction (A2 · p_int = 0.022) — 이전 문헌은 driver-independent prognostic score 만 제시.",
          "② 병리과 이미 사용중인 IHC 3-plex (TG + PAX8 + TTF-1) 로 BRAF+ PFI 를 Full-8 보다 강한 log-rank p = 1.7 × 10⁻⁴ 로 분리 — <b>immediate translational readiness</b>.",
          "③ Lee 2024 (Korean n = 370) 에서 axis Mann-Whitney p = 2.7 × 10⁻⁷ 로 완전 재현 — 새로운 population 에서 first time.",
          "종합: 이전 연구는 axis 존재를 확인, 우리는 <b>treatment-selection biomarker 로 승격 + immediate clinical deployment path</b> 제시."
        ],
        evidence: "Figures 7-8 · A2 interaction · A2-R Lee 2024 · IHC 3-plex killer",
        strength: "strong"
      },
      {
        q: "이 결과가 정말 practice-changing 이라고 판단할 근거가 있나?",
        a: [
          "즉시 practice-changing 은 아직 아니지만 <b>3 단계 practice-changing pathway</b> 을 지원함.",
          "Stage 1 (0-6 개월): 분당 후향 chart-review — 기존 FFPE + BRAF genotyping + 이미 있는 IHC 결과 재분석. 추가 assay 불필요.",
          "Stage 2 (6-18 개월): 분당 n = 200 prospective (A5 simulation power ≥ 90%).",
          "Stage 3 (18-36 개월): 다기관 확장 → guideline 등재 논의.",
          "즉, 논문 자체는 practice-informing → validation 후 practice-changing."
        ],
        evidence: "A5 SNUBH simulation · roadmap",
        strength: "moderate"
      },
      {
        q: "이런 정도의 임상 프레임이면 JCI Insight 급이 아닌가?  왜 top-tier 인가?",
        a: [
          "3 가지 top-tier 근거:",
          "① Predictive interaction (BRAF-restricted) — top-tier journal 이 요구하는 driver-context specificity 확보.",
          "② Multi-modality convergence (5 modality Spearman ρ = 0.71-0.94) — reviewer 의 '한 modality artifact' 공격 자동 차단.",
          "③ IHC 3-plex 즉시 배포 가능 — reviewer 가 항상 요구하는 'clinical translation ready' 조건 완비.",
          "만약 A2 interaction 이 분당 prospective 에서 replicate 되면 <b>NC/Nature Medicine 상한</b>, 실패 시 JCI Insight fallback 이 정직한 range."
        ],
        strength: "moderate"
      },
      {
        q: "Cell / Nature 급 novelty 는 없나?",
        a: [
          "정직하게 — Cell/Nature 급 causal biology (mouse model + rescue) 는 없음.  본 논문은 <b>observational precision medicine</b> tier.",
          "Cell/Nature 로 가려면 in-vivo redifferentiation 모델 + DM1 axis restoration + RAI 재감수성 실험 3-lane 필요 — 이는 다음 논문 target.",
          "본 논문의 정확한 pitch: <b>'observational treatment-selection biomarker with immediate clinical deployment path'</b>."
        ],
        strength: "moderate"
      }
    ]
  },
  {
    key: "methods",
    title: "Reviewer #1  ·  Statistics / Methodology",
    subtitle: "통계 리뷰어의 aggressive attack",
    color: "#B91C1C",
    qas: [
      {
        q: "A2 interaction p = 0.022 는 <b>single-cohort finding</b> · N events = 33 · <b>replication 없다</b>.  이거로 treatment-selection biomarker 라 주장하나?",
        a: [
          "이 지적 accept.  Interaction 자체는 single-cohort, replication 은 아직 없음.  이를 논문에서 명시적으로 밝힘 (Results §3, Limitations, Discussion).",
          "우리 defense 3-layer:",
          "① <b>Axis 자체는 replicated</b> — Lee 2024 p = 2.7 × 10⁻⁷ · GPL570 3 코호트 · Mun 2025 protein 7/7 등 6 external cohort.",
          "② <b>Interaction 재현이 안 된 이유는 biology 가 아닌 데이터 접근성</b> — 외부 코호트에 BRAF genotype + PFI 조합 없음.",
          "③ <b>Prospective replication plan 명시</b> (Methods · Monte Carlo n = 200 power ≥ 90%).",
          "정확한 framing: <b>hypothesis-generating with well-defined replication path</b>.  treatment-selection biomarker 는 conditional claim (replication 조건).",
          "Reviewer 가 이 정도로도 부족하다 하면 → <b>title 을 'a candidate treatment-selection framework' 로 축소</b> 하고 abstract 도 그렇게 조정."
        ],
        evidence: "Fig 7 (A2 + A2-R) · Methods interaction Cox · Limitations",
        strength: "moderate"
      },
      {
        q: "8-gene 만 골라서 clustering 한 결과 아닌가?  panel-induced artifact 로 볼 수도 있다.",
        a: [
          "이 attack 은 이미 방어됨 (Fig 1e).",
          "Pan-genome top-5000 MAD clustering ARI = <b>0.92</b> · TIERA67 (67 gene) ARI = 0.90 · driver-anchor genes ARI = <b>-0.007</b>.",
          "즉 <b>어떤 unbiased gene set 을 써도 같은 DM1/DM2 partition 이 나온다</b> — panel 은 axis 를 발견한 것이 아니라 <b>측정 lens</b>.",
          "Driver-anchor 만으로는 무작위 이하 → 이 축은 driver 정체성이 아닌 <b>coordinated lineage program</b>."
        ],
        evidence: "Fig 1e · ARI ladder",
        strength: "strong"
      },
      {
        q: "MSK-IMPACT 코호트는 advanced disease enriched — pooled HR = 2.53 은 selection bias 아닌가?",
        a: [
          "정직하게 <b>부분 accept</b> — MSK 기여가 pooled 를 끌어올림.",
          "3 가지 mitigation:",
          "① TCGA 단독 HR = 2.30 (95% CI 0.77-6.88) — same direction, 낮은 events 로 CI 넓음.",
          "② I² = 0% — TCGA 와 MSK 사이 heterogeneity 없음.",
          "③ Advanced disease 도 우리 target population — <b>바로 이 그룹이 clinically actionable</b>.",
          "논문에서 이를 명시적으로 caveat — 'pooled HR informs cross-stage biology, prospective validation for treatment decisions'."
        ],
        evidence: "Fig 5 · Discussion Limitations",
        strength: "moderate"
      },
      {
        q: "Multiple testing correction 안 했다.  14 combinations × BRAF+ PFI × 여러 endpoint = 다중검정 문제.",
        a: [
          "3 가지 방어:",
          "① Primary hypothesis 는 <b>Full-8 DM1 axis</b> · 다른 combinations 는 exploratory (Methods 에서 명시).",
          "② IHC-3 finding 은 <b>routine biology-motivated a priori</b> (병리과 이미 시행중) → post-hoc data mining 아님.",
          "③ Bonferroni 로 α = 0.05/14 = 0.0036 적용해도 IHC-3 p = 1.7 × 10⁻⁴ survive.",
          "논문에서 명시적 sentence 추가: 'IHC 3-plex was a biology-motivated a-priori hypothesis, not a data-mining artefact'."
        ],
        strength: "strong"
      },
      {
        q: "In-silico HM450 β proxy 로 IHC 결과 예측했다.  실제 IHC 데이터는 없다.",
        a: [
          "이 지적 accept.  현재는 methylation proxy — 실제 IHC 결과 없음.",
          "Rationale for proxy:",
          "① TCGA HM450 β 는 sample-level lineage silencing 을 정량 측정.",
          "② Mun 2025 protein 7/7 direction concordant → protein-level lineage silencing 재현 (indirect IHC 근거).",
          "③ Landa 2016 lineage silencing gene list 5/8 overlap → biological plausibility 확립.",
          "Prospective validation plan: <b>분당 pilot 30 명 · H-score 2-pathologist blind · κ > 0.7 · 3-tier scoring KM 재검정</b>.",
          "논문 abstract + methods 에서 명시적으로 in-silico proxy 임을 밝힘."
        ],
        strength: "moderate"
      },
      {
        q: "Cox proportional hazards 가정 (PH assumption) 검정했나?",
        a: [
          "정직하게 — 현재 draft 에서 명시적 PH assumption test (Schoenfeld residuals) 는 포함 안 됨.  Revision 시 즉시 추가.",
          "예상 결과: BRAF+ subset 은 event 수 33 로 formal test power 제한.  Log-log KM plot 은 crossing 없음 → PH assumption 대체로 만족.",
          "Sensitivity analysis: time-dependent Cox coefficient 확인 예정."
        ],
        strength: "moderate"
      },
      {
        q: "Interaction term 은 정확히 무엇을 modeling 하나?  standardized 인가 raw 인가?",
        a: [
          "명확한 답:",
          "Model: h(t | z, BRAF) = h₀(t) · exp(β₁·z + β₂·BRAF + β₃·(z × BRAF))",
          "여기서 z = 8-gene panel mean β 를 <b>within-cohort z-standardize</b> 한 값 (mean 0, SD 1).",
          "BRAF = binary indicator (V600E+ = 1, else 0).",
          "Interaction coefficient β₃ = <b>-0.75 (HR_int = e^(-0.75) = 0.47)</b>, p = 0.022 (two-sided Wald).",
          "즉 BRAF+ 에서 DM1 score 1 SD 증가 시 hazard 가 exp(β₁ + β₃) 로 변화 → BRAF- 대비 exp(β₃) = 0.47 배 다른 방향."
        ],
        strength: "strong"
      }
    ]
  },
  {
    key: "biology",
    title: "Reviewer #2  ·  Biology / Mechanism",
    subtitle: "Molecular biology 리뷰어",
    color: "#0E7490",
    qas: [
      {
        q: "왜 <b>정확히 이 8 gene</b> 인가?  TDS-16 (Yoo 2016, 16 gene) 이 이미 표준인데 왜 8 로 줄였나?",
        a: [
          "3 가지 이유:",
          "① <b>Biology parsimony</b> — 8 gene 은 iodine handling 회로의 최소 완전 표상 (TSHR → cAMP → NIS → TPO → TG · 3 lineage TF).  TDS-16 의 나머지 8 gene 은 대부분 redundant (T3/T4 metabolism 관련).",
          "② <b>Head-to-head 비교</b> (Fig 7a) — 우리 8-gene AUC = 0.887 vs Landa-5 overlap 0.904 vs RAI-machinery-5 0.904.  동등 이상.",
          "③ <b>Clinical translation 실용성</b> — IHC 3-plex 로 축소 가능 (병리과 routine), 16 gene 은 routine 불가능."
        ],
        evidence: "Fig 3d-e · Fig 7a · TDS-16 등가 검정",
        strength: "strong"
      },
      {
        q: "TG, TPO, TSHR 는 갑상선 전문 marker 니까 dedifferentiation 시 down 되는 게 당연 아닌가?  novel biology 아니다.",
        a: [
          "부분 accept — individual gene down-regulation 은 novel 아님.  Novel 은 3 가지:",
          "① <b>Coordinated axis 로서 정량화</b> — 개별 marker 가 아닌 8-gene 통합 score.",
          "② <b>Driver-specific interaction</b> (A2) — BRAF+ 에서만 progression-free interval 예측 (interaction p = 0.022).",
          "③ <b>Multi-modality convergence</b> — HM450 β · RNA · protein · Landa · GSE151179 5 modality Spearman ρ = 0.71-0.94.",
          "즉 <b>개별 marker 는 알려진 것, coordinated axis + predictive interaction + multi-omic convergence 는 새로움</b>."
        ],
        strength: "moderate"
      },
      {
        q: "MAPK 활성 (BRAF/RET fusion) → NIS silencing 은 이미 Chakravarty 2011 (JCI) 에서 밝혔다.  뭐가 새로운가?",
        a: [
          "Chakravarty 2011 은 <b>mouse model 에서 MAPK inhibitor 로 NIS 재발현 + RAI uptake 회복</b> — mechanism 확립.  하지만 <b>어느 human patient 가 이 pathway responder 인지 identify 하는 pre-treatment classifier 는 없었음</b>.",
          "우리 DM1 axis 가 그 <b>missing classifier</b>.",
          "임상 논리: BRAF+ 환자 중 DM1-high subset = MAPK 저해제 redifferentiation therapy 후보군 (Ho 2013 NEJM · Rothenberg 2015 CCR 의 responder subset 정의).",
          "즉 Chakravarty 2011 은 pathway, 우리는 <b>patient selection tool</b>."
        ],
        strength: "strong"
      },
      {
        q: "FOXE1 · SLC5A5 는 우리 panel 에서 weak (d < 0.9).  왜 포함시켰나?",
        a: [
          "3 가지 이유:",
          "① <b>Biology completeness</b> — FOXE1 = 3 lineage TF 의 하나 (congenital hypothyroidism gene), SLC5A5 = RAI 흡수의 직접 대상 protein.  이 두 gene 없이 iodine circuit 표상 불완전.",
          "② <b>Head-to-head</b> — 이 두 gene 제거한 Compact 6 (Fig GR-1) 성능 = Full-8 과 동등 → deployment 시 이 두 제외 가능.",
          "③ <b>Panel = maximum information / IHC = deployment</b> 로 다른 layer 사용."
        ],
        evidence: "Fig GR-1 · Compact 6 · Table biomarker combinations",
        strength: "moderate"
      },
      {
        q: "Methylation 이 원인인가 결과인가?  Correlational 만으로는 causal 주장 안 된다.",
        a: [
          "완전 agree — correlational, causal 주장 안 함.",
          "논문에서 명시적으로 '<b>promoter methylation is a correlate of lineage-programme silencing, not necessarily its cause</b>'.",
          "Causal 근거는 Chakravarty 2011 mouse model + Ho 2013 clinical restoration 에 의존 (우리 causal 실험 아님).",
          "Discussion §3.1 (author voice) 이 이 boundary 를 정확히 다룸."
        ],
        strength: "moderate"
      },
      {
        q: "GSE151179 refractory group n = 4 인데 이걸로 뭘 결론 낼 수 있나?",
        a: [
          "이 지적 완전 accept.  본 논문에서 GSE151179 는 <b>hypothesis-generating pilot</b> 으로만 사용 · confirmatory 아님.",
          "논문 명시 sentence: 'We interpret the GSE151179 comparison as a biological anchor rather than as a validated response-prediction model.'",
          "Confirmatory replication 은 분당 prospective (n = 200, expected refractory events ~ 20-25) 가 target — A5 simulation 이 power ≥ 90% 확인.",
          "Reviewer 가 GSE151179 삭제 요청하면 <b>Supplementary 이동 가능</b>."
        ],
        strength: "weak"
      },
      {
        q: "Single-cell 에서 DM1 signature 는 어느 세포 유형에서 나오나?  Stromal contamination 없나?",
        a: [
          "Lu 2023 (GSE193581) sc-RNA-seq 에서 <b>KRT8/KRT19/EPCAM+ thyrocyte-lineage cells 만 restrict</b> 하고도 DM1 gradient 유지 (Fig 4b).",
          "즉 <b>thyrocyte-intrinsic signal</b> — stromal / immune contamination 아님.",
          "Pu 2021 per-patient tumour 와 adjacent-normal thyrocyte concordance ρ = 0.798-0.886 (Bonferroni p < 10⁻¹⁰).",
          "이 두 결과로 stromal artifact 공격 차단."
        ],
        evidence: "Fig 4a-b",
        strength: "strong"
      }
    ]
  },
  {
    key: "clinical",
    title: "Reviewer #3  ·  Clinical Translation",
    subtitle: "Clinical thyroid oncologist 관점",
    color: "#0F172A",
    qas: [
      {
        q: "IHC 3-plex 결과가 정말 clinical 하게 사용 가능한가?  Pathologist 간 concordance 는?",
        a: [
          "현재 데이터는 methylation β proxy 기반 · 실제 IHC quantification 없음.",
          "Validation plan:",
          "① 분당 pilot 30 명 · BRAF+ 우선 · 2 pathologist blind · κ target > 0.7.",
          "② H-score (0-300) 또는 % positive tumor cell 로 정량.",
          "③ 3-tier scoring (low / intermediate / high) 로 KM 재검정.",
          "④ Cut-off 결정: median split 우선 · training/validation 50-50.",
          "이 protocol 은 논문 Methods 에 이미 명시."
        ],
        evidence: "Discussion Clinical implications · Methods §IHC in-silico proxy",
        strength: "moderate"
      },
      {
        q: "BRAF+ 서브셋에서만 predictive 라면 clinical value 가 제한적 아닌가?  BRAF- / RAS+ 환자는 어쩌나?",
        a: [
          "부분 accept — BRAF- / RAS+ 환자는 현재 axis 로 treatment-selection 불가.",
          "하지만 <b>BRAF+ 는 PTC 의 60-70%</b> (Korean population 60%+) → <b>대다수 환자군 cover</b>.",
          "BRAF- subset 은 <b>Landa 5-overlap · RAI-machinery-5 등 subset score 로 별도 stratify</b> 가능 (Fig 7a, BRAF+ 아닌 코호트도 log-rank 유의).",
          "논문 Discussion 에 이를 명시: 'BRAF+ 서브셋 = primary target; BRAF- 는 secondary axis 로 stratified 관리'."
        ],
        strength: "moderate"
      },
      {
        q: "TSO500 v2 이미 병원에서 표준 검사인데, 3 gene 만으로도 성능 저하 (AUC 0.76) 이면 clinical utility 낮다.",
        a: [
          "이거 완전 다뤘음.  <b>TSO500 DNA-level 만으로는 부족</b> (Fig 8a).  <b>Dual-assay 전략</b> 이 우리 권고:",
          "① TSO500 v2 (기존 routine) → PAX8 · NKX2-1 · TSHR 의 LOH / SNV / CNV.",
          "② + <b>NanoString / RT-qPCR spike-in 2-gene (TPO + DIO1)</b> → Compact 5 = AUC 0.940 (Full-8 동등).",
          "③ 또는 IHC 3-plex (TG + PAX8 + NKX2-1) 만으로 BRAF+ PFI p = 1.7 × 10⁻⁴ 달성 — <b>가장 실용적</b>.",
          "즉 3 개 deployment path 모두 존재."
        ],
        evidence: "Fig 8 · §SUPP-10 TSO500 caveat + reduction feasibility",
        strength: "strong"
      },
      {
        q: "실제 임상에서 어떻게 사용하나?  Reflex 알고리즘 상세히.",
        a: [
          "Reflex 4-step algorithm:",
          "① 수술 후 FFPE → BRAF genotyping (routine).",
          "② BRAF V600E+ 시 → IHC 3-plex (TG + PAX8 + TTF-1) 자동 발동 (reflex).",
          "③ IHC 3-plex 판독 → DM1-like (low score) 시 → 추가 fusion testing (RET · NTRK · ALK) · MTB 회의 case 등록.",
          "④ Post-op RAI 계획 시 → DM1-like BRAF+ = <b>고용량 RAI 반복 회피 + redifferentiation therapy sequencing 고려</b>.",
          "This flow 는 Figure 6 (translational pathway) 에 시각화."
        ],
        evidence: "Fig 6 · IHC 3-plex reflex protocol",
        strength: "strong"
      },
      {
        q: "이 axis 로 얼마나 많은 환자가 실제 이득을 볼 것 같나?  Number Needed to Treat / Screen?",
        a: [
          "Rough estimate (분당 코호트 기준):",
          "PTC 연 200 명 · BRAF+ 120 명 (60%) · 이 중 DM1-like 30 명 (25%) = <b>연 30 명이 primary target</b>.",
          "이 30 명이 표준 RAI 반복 대신 redifferentiation therapy sequencing 시 예상 이익: RAI toxicity 감소 · 재발 조기 인지.",
          "NNS (screening) = 200 · NNT (BRAF+ only) = 4 · Actionable NNT (DM1+BRAF+ redifferentiation candidate) = 1 (definitional).",
          "Note: 이 estimate 는 approximation · 실제 impact 은 prospective 에서 확정."
        ],
        strength: "moderate"
      },
      {
        q: "Toxicity endpoint 없이 harm avoidance 주장 가능한가?",
        a: [
          "정직하게 accept — 현재 논문에 direct RAI toxicity outcome 데이터 없음.  이는 <b>major limitation</b>.",
          "논문에서 명시적 caveat: 'harm avoidance is a biological hypothesis based on iodine-handling machinery silencing, not a validated toxicity outcome'.",
          "Toxicity endpoint 는 <b>분당 prospective 에서 primary secondary outcome</b>: 누적 RAI dose · leukopenia · sialoadenitis · secondary malignancy.",
          "이 protocol IRB draft 에 이미 포함."
        ],
        strength: "weak"
      },
      {
        q: "Cost-effectiveness analysis 있나?",
        a: [
          "현재 없음 — Post-validation 이후 진행 예정.",
          "Rough logic: IHC 3-plex 는 이미 시행중 (marginal cost = 0), BRAF genotyping 이미 시행중 → <b>추가 cost 거의 없음</b>.",
          "반면 불필요한 고용량 RAI 회피 시 환자당 US$5,000-10,000 절감 + toxicity 관련 cost 감소.",
          "정식 QALY-based CEA 는 prospective validation 후 진행 예정."
        ],
        strength: "weak"
      }
    ]
  },
  {
    key: "external",
    title: "Reviewer 관심사  ·  External Validity / Generalizability",
    subtitle: "Population + Platform generalizability",
    color: "#B45309",
    qas: [
      {
        q: "Lee 2024 replication (p = 2.7 × 10⁻⁷) 은 강력하지만 <b>DM1 axis 자체</b> replication 이지 <b>A2 interaction</b> 은 아니다.",
        a: [
          "완전 accept — 이 distinction 이 정확하다.",
          "논문에서 이를 명시적으로 다룸: <b>'axis replication is established across 6 external cohorts; predictive interaction replication requires a cohort with joint BRAF genotype and PFI annotation, which is not available in current public cohorts.'</b>",
          "Prospective replication path: 분당 n = 200 · BRAF genotyping + PFI 5-yr follow-up · A5 simulation power ≥ 90%.",
          "즉 우리 defense = <b>axis is real (6-way replicated); interaction awaits prospective test</b>."
        ],
        strength: "moderate"
      },
      {
        q: "TCGA 는 서구 population, Lee 2024 는 한국 — Asian-specific finding 가능성 있다.",
        a: [
          "부분 accept — East Asian PTC 는 driver 분포 · 예후 profile 다름.  이 사실을 이용:",
          "① Korean cohort DM1 prevalence = <b>37.8%</b> vs TCGA 28.4% → 한국인에서 DM1 상대적으로 흔함 → 임상 impact 더 큼.",
          "② <b>DM1 axis 방향 동일</b> → biology 는 population-independent.",
          "③ Effect size 는 population 마다 다를 수 있음 — 이는 논문에서 명시.",
          "즉 Population-specific magnitude, universal direction."
        ],
        strength: "moderate"
      },
      {
        q: "GPL570 3 cohort 는 <b>micro-arry n < 100</b>.  이거 large-N 재현으로 인정 못한다.",
        a: [
          "GPL570 코호트는 support 로만 사용 · confirmatory 아님.",
          "Primary replication = Lee 2024 (n = 370) · Landa 2016 (n = 88 well-characterized) · Mun 2025 (n = 336 proteomics).",
          "GPL570 은 <b>platform robustness</b> (microarray) 를 보이는 것이 목적 — Fig 7c 에서 명확히 라벨."
        ],
        strength: "moderate"
      },
      {
        q: "FFPE 재현성 KS p = 0.44 는 <b>equivalence 검정이 아니라 difference 검정</b>.  Absence of evidence ≠ evidence of absence.",
        a: [
          "완전 accept — 이는 formal equivalence 검정이 아님.",
          "정확한 framing: FFPE vs Fresh-Frozen distributional shift 를 <b>detect 못함</b> · equivalence 확립 아님.",
          "Formal equivalence 는 pre-specified margin + TOST 필요 — 이는 <b>분당 prospective 에서 secondary aim</b>.",
          "논문 문구 수정: 'not detectably shifted' → 'not detectably different' (더 정확한 표현)."
        ],
        strength: "weak"
      },
      {
        q: "Pu 2021 · Lu 2023 sc-RNA-seq 재현 결과는 정확히 뭘 보여주나?",
        a: [
          "3 가지:",
          "① <b>Thyrocyte-intrinsic signal</b> — KRT8/KRT19/EPCAM+ 세포 만으로도 gradient 유지.",
          "② <b>Per-patient concordance</b> (Pu 2021) — tumor 와 adjacent normal 사이 Spearman ρ = 0.798-0.886 (Bonferroni p < 10⁻¹⁰).",
          "③ <b>Cross-cohort thyroid-intrinsic direction</b> — Lu 2023 UMAP 8-gene score gradient 방향 일치.",
          "즉 stromal / immune 오해 방지 + patient-level 재현."
        ],
        evidence: "Fig 4a-b",
        strength: "strong"
      }
    ]
  },
  {
    key: "framing",
    title: "Cross-cutting  ·  Framing / Interpretation",
    subtitle: "논문 전체 framing 관련 질문",
    color: "#047857",
    qas: [
      {
        q: "제목이 'predicts refractoriness' 인데, 그 정도 강한 claim 을 정말 뒷받침하나?",
        a: [
          "이 지적 부분 accept — 'predicts' 는 강한 표현.",
          "Defensible position:",
          "① Retrospective evidence: TCGA A2 interaction p = 0.022 · IHC 3-plex BRAF+ PFI p = 1.7 × 10⁻⁴ · GSE151179 post-RAI DM1-like state.",
          "② 이 3 개 layer 가 함께 <b>predictive claim 을 support</b>.",
          "③ 명시적 caveat: 'predictive claim is anchored in retrospective and pilot data; prospective replication in n = 200 SNUBH cohort is the confirmatory target'.",
          "만약 reviewer 가 이거도 부족하다 하면 → title 을 <b>'may predict refractoriness'</b> 로 조정 가능."
        ],
        strength: "moderate"
      },
      {
        q: "'Treatment-selection biomarker' claim 이 정말 가능한가?  RCT 없이?",
        a: [
          "정직하게 — <b>full treatment-selection biomarker</b> claim 은 RCT 필요, 우리는 아직.",
          "우리 claim 은 <b>'candidate treatment-selection biomarker'</b> (Abstract · Discussion).",
          "이는 predictive interaction 존재 (A2 p_ix = 0.022) + BRAF+ subset 특이성 (Fig 7b) 로 support.",
          "Full validation 은 SNUBH prospective (Phase 2) + 이후 다기관 RCT (Phase 3).",
          "단어 선택 careful: '<b>candidate</b>' 은 aspirational, '<b>established</b>' 는 RCT 후.  이 boundary 명시."
        ],
        strength: "moderate"
      },
      {
        q: "'Harm avoidance' vs 'Treatment selection' — 두 framing 을 동시에 주장하는 게 self-contradictory 아닌가?",
        a: [
          "Not contradictory · <b>same biology, two application layers</b>:",
          "① <b>Harm avoidance</b> (driver-agnostic) — 어떤 driver 든 DM1-high = iodine circuit silenced = 표준 RAI 반복 이득 낮음.",
          "② <b>Treatment selection</b> (BRAF+ specific) — DM1-high AND BRAF+ = MAPK inhibitor redifferentiation therapy 후보.",
          "즉 harm avoidance 는 axis biology, treatment selection 은 interaction biology.  Discussion §3 에 이를 정확히 분리."
        ],
        strength: "strong"
      },
      {
        q: "'RAI-refractory' 정의가 무엇인가?  Guideline 마다 다르다.",
        a: [
          "논문에서 <b>ATA 2015 · 2025 definition</b> 사용 명시 예정:",
          "- Structural disease progression despite RAI",
          "- No iodine uptake despite adequate TSH stimulation",
          "- Disease progression within 6-12 months of RAI",
          "- Cumulative RAI dose > 600 mCi without response.",
          "GSE151179 label 은 원 저자 clinical annotation 따름.",
          "Methods 에 정확한 definition 추가."
        ],
        strength: "moderate"
      },
      {
        q: "Journal fit — 우리 논문이 Nature Communications 냐 Nature Medicine 이냐?",
        a: [
          "정직한 assessment:",
          "① <b>Nature Communications</b> = 현재 draft 로 적절.  Discovery + replication + IHC translation + prospective plan.",
          "② <b>Nature Medicine</b> = A2 interaction 이 SNUBH 에서 replicate 되고 direct H-score IHC 데이터 확보 시 upgrade 가능.",
          "③ <b>JCI Insight</b> = fallback (interaction replication 실패 시).",
          "Editorial strategy: <b>NC pre-submission inquiry → 반응 강하면 Nature Med 도전, 약하면 NC 정식 제출</b>."
        ],
        strength: "moderate"
      },
      {
        q: "Limitations 를 한 줄로.",
        a: [
          "3-line summary:",
          "① A2 interaction (predictive) = single-cohort (TCGA), replication awaits SNUBH prospective.",
          "② IHC 3-plex evidence = HM450 β proxy, direct H-score validation pending.",
          "③ MSK-IMPACT advanced-disease enrichment inflates pooled survival HR — biology direction preserved.",
          "이 3 개는 논문 Limitations section 에 명시 · Reviewer 방어 시 first-line acknowledgment."
        ],
        strength: "strong"
      },
      {
        q: "이 논문이 accept 되면 <b>next paper</b> 는 뭔가?",
        a: [
          "3-paper series 전략:",
          "① 본 논문 (NC/Nature Med) = axis discovery + IHC 3-plex + BRAF+ interaction (retrospective + Lee 2024 replication).",
          "② 다음 논문 (2-3 년 후) = 분당 prospective n = 200 · IHC direct H-score · treatment-selection biomarker confirmation.",
          "③ 3rd 논문 (Cell / Nature 급) = in-vivo redifferentiation model · DM1 axis restoration · RAI 재감수성 · causal biology.",
          "즉 discovery → confirmation → mechanism 3-lane."
        ],
        strength: "strong"
      }
    ]
  }
];

const STRENGTH_LABEL = {
  strong:   { txt: "★★★ 강력 방어",       color: "#047857", bg: "#ECFDF5" },
  moderate: { txt: "★★ 중간 방어",         color: "#B45309", bg: "#FEF3C7" },
  weak:     { txt: "★ 약함 · caveat 인정", color: "#B91C1C", bg: "#FEF2F2" },
  escalate: { txt: "⚠ 저자 판단 필요",     color: "#7C2D12", bg: "#FEE2E2" }
};

export const ReviewerQA: React.FC = () => {
  const [openCat, setOpenCat] = useState<string | null>("editor");
  const [openIdx, setOpenIdx] = useState<Record<string, Set<number>>>({});
  const toggleQA = (cat: string, i: number) => {
    setOpenIdx(prev => {
      const cur = new Set(prev[cat] || []);
      if (cur.has(i)) cur.delete(i); else cur.add(i);
      return { ...prev, [cat]: cur };
    });
  };
  const openAllInCat = (cat: string, count: number) => {
    setOpenIdx(prev => ({ ...prev, [cat]: new Set(Array.from({length: count}, (_,i)=>i)) }));
  };
  const closeAllInCat = (cat: string) => setOpenIdx(prev => ({ ...prev, [cat]: new Set() }));

  const totalQ = CATEGORIES.reduce((s, c) => s + c.qas.length, 0);
  const strengthCount = CATEGORIES.reduce((acc, c) => {
    c.qas.forEach(q => { acc[q.strength] = (acc[q.strength] || 0) + 1; });
    return acc;
  }, {} as Record<string, number>);

  return (
    <div className="revqa">
      <div className="revqa-intro">
        <div className="revqa-intro-title">🎯 Reviewer Q&A 시뮬레이션 — {totalQ} 개 예상 질문 · 5 카테고리</div>
        <div className="revqa-intro-body">
          NC / Nature Medicine 편집자 + 3 명 reviewer 관점에서 <b>가장 aggressive 한 예상 질문</b> 을 pre-empt.
          각 답변은 (a) 핵심 방어 논리, (b) 근거 figure / 수치, (c) 방어 강도 등급 포함.
          미팅 · rebuttal 전 참조용.
        </div>
        <div className="revqa-intro-stats">
          <span className="revqa-stat" style={{background: STRENGTH_LABEL.strong.color}}>★★★ 강력 방어 {strengthCount.strong || 0}</span>
          <span className="revqa-stat" style={{background: STRENGTH_LABEL.moderate.color}}>★★ 중간 {strengthCount.moderate || 0}</span>
          <span className="revqa-stat" style={{background: STRENGTH_LABEL.weak.color}}>★ 약함 {strengthCount.weak || 0}</span>
        </div>
      </div>

      {CATEGORIES.map(cat => (
        <div key={cat.key} className="revqa-cat">
          <div className="revqa-cat-head" onClick={() => setOpenCat(openCat === cat.key ? null : cat.key)}
               style={{borderLeftColor: cat.color}}>
            <div>
              <div className="revqa-cat-title" style={{color: cat.color}}>{cat.title}</div>
              <div className="revqa-cat-sub">{cat.subtitle}  ·  {cat.qas.length} 개 질문</div>
            </div>
            <div className="revqa-cat-toggle">{openCat === cat.key ? "▴ 접기" : "▾ 펼치기"}</div>
          </div>

          {openCat === cat.key && (
            <div className="revqa-cat-body">
              <div className="revqa-cat-actions">
                <button type="button" onClick={() => openAllInCat(cat.key, cat.qas.length)}>모든 답변 펼치기</button>
                <button type="button" onClick={() => closeAllInCat(cat.key)}>모두 접기</button>
              </div>
              {cat.qas.map((qa, i) => {
                const isOpen = (openIdx[cat.key] || new Set()).has(i);
                const s = STRENGTH_LABEL[qa.strength];
                return (
                  <div key={i} className={`revqa-item ${isOpen ? "open" : ""}`}>
                    <div className="revqa-q" onClick={() => toggleQA(cat.key, i)}>
                      <div className="revqa-q-num">Q{i+1}</div>
                      <div className="revqa-q-txt" dangerouslySetInnerHTML={{__html: qa.q}} />
                      <span className="revqa-str-badge" style={{color: s.color, background: s.bg, borderColor: s.color}}>{s.txt}</span>
                    </div>
                    {isOpen && (
                      <div className="revqa-a">
                        <div className="revqa-a-label">A · 우리 방어</div>
                        {qa.a.map((line, j) => (
                          <p key={j} dangerouslySetInnerHTML={{__html: line}} />
                        ))}
                        {qa.evidence && (
                          <div className="revqa-evidence">
                            <b>근거:</b> {qa.evidence}
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                );
              })}
            </div>
          )}
        </div>
      ))}

      <div className="revqa-conclusion">
        <div className="revqa-conclusion-title">🎯 종합 전략</div>
        <ul>
          <li><b>First-line 방어</b>: 강력 등급 질문에는 정면 답변 + 근거 figure 인용.</li>
          <li><b>Second-line 방어</b>: 중간 등급 질문에는 caveat 명시 후 next-step 제시.</li>
          <li><b>Retreat position</b>: 약함 등급 질문에는 <b>솔직히 accept</b> + 논문 문구 조정 제안.</li>
          <li><b>절대 하지 말 것</b>: 약점을 defense 하려 무리한 override → 오히려 reviewer 신뢰 상실.</li>
          <li><b>Editorial timing</b>: NC pre-submission inquiry → 반응 강하면 Nature Med, 약하면 NC direct.</li>
        </ul>
      </div>
    </div>
  );
};
