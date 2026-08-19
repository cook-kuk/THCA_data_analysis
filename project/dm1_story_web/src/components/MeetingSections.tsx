import React from "react";

/* ================================================================
   ExecutiveSummary — Problem / Solution / Clinical Value 3 카드
   ================================================================ */
export const ExecutiveSummary: React.FC = () => (
  <div className="exec-summary">
    <div className="exec-card exec-problem">
      <div className="exec-badge">A. 문제</div>
      <h3>현재는 RAI 를 해본 뒤에야 refractory 여부를 알 수 있다</h3>
      <p>
        평균 <b>6 – 12 개월 지연</b>이 발생하고, 그 사이 불필요한 고용량 ¹³¹I 반복으로
        골수억제 · leukopenia · 장기 독성 위험이 누적된다.
        갑상선암은 장기 생존 환자가 많아 누적 RAI 독성의 임상적 의미가 크다.
      </p>
    </div>

    <div className="exec-card exec-solution">
      <div className="exec-badge">B. 해법</div>
      <h3>8 개 유전자로 iodine-handling molecular state 를 미리 분류한다</h3>
      <p>
        <b>TG · TPO · TSHR · SLC5A5/NIS · DIO1 · PAX8 · NKX2-1/TTF-1 · FOXE1/TTF-2</b>
        — 갑상선 분화 / RAI 흡수 회로의 핵심 8 유전자로 종양을
        <b> DM1 (Iodine-handling-low) · DM2 (Iodine-handling-high)</b> 두 분자 상태로 분류한다.
      </p>
    </div>

    <div className="exec-card exec-value">
      <div className="exec-badge">C. 임상 가치</div>
      <h3>불필요한 고용량 RAI 회피 + 고위험군 조기 치료 전환</h3>
      <p>
        Primary: <b>"RAI 가 안 들을 환자에게 불필요한 고용량 RAI 를 반복하지 않게 한다"</b><br/>
        Secondary: aggressive 환자에서는 RAI 를 오래 끌지 않고 더 빨리 systemic 치료로 전환.
        새로운 치료법을 주장하는 것이 아닌 임상 의사결정 보조 도구.
      </p>
    </div>
  </div>
);

/* ================================================================
   StateComparisonCards — DM1 vs DM2 큰 비교
   ================================================================ */
export const StateComparisonCards: React.FC = () => (
  <div className="state-compare">
    <div className="state-card state-dm2">
      <div className="state-header">
        <div className="state-tag">DM2</div>
        <div className="state-newname">Iodine-handling-<b>HIGH</b></div>
      </div>
      <table className="state-table">
        <tbody>
          <tr><th>내부 코드</th><td>DM2</td></tr>
          <tr><th>제안 명칭</th><td>Iodine-handling-high</td></tr>
          <tr><th>유병률 (TCGA)</th><td>71.6 %  (360 / 504)</td></tr>
          <tr><th>분자 상태</th><td>갑상선 분화 유지 · iodine-handling on</td></tr>
          <tr><th>8-gene 발현</th><td>높음 (TF₃ + effector₅ 모두 활성)</td></tr>
          <tr><th>HM450 mean β</th><td>0.253  (낮은 promoter 메틸화)</td></tr>
          <tr><th>임상 해석</th><td>RAI-responsive-like · 표준 RAI 치료 적용</td></tr>
          <tr><th>주의</th><td>institutional validation 필요</td></tr>
        </tbody>
      </table>
    </div>

    <div className="state-card state-dm1">
      <div className="state-header">
        <div className="state-tag">DM1</div>
        <div className="state-newname">Iodine-handling-<b>LOW</b></div>
      </div>
      <table className="state-table">
        <tbody>
          <tr><th>내부 코드</th><td>DM1</td></tr>
          <tr><th>제안 명칭</th><td>Iodine-handling-low</td></tr>
          <tr><th>유병률 (TCGA)</th><td>28.4 %  (140 / 504)</td></tr>
          <tr><th>분자 상태</th><td>dedifferentiated · iodine-handling off</td></tr>
          <tr><th>8-gene 발현</th><td>낮음 (TF₃ + effector₅ 동시 silencing)</td></tr>
          <tr><th>HM450 mean β</th><td>0.385  (+52 % promoter 과메틸화)</td></tr>
          <tr><th>임상 해석</th><td>RAI-refractory-like / aggressive-like</td></tr>
          <tr><th>주의</th><td>치료 선택 biomarker 아님 — risk stratification 만</td></tr>
        </tbody>
      </table>
    </div>

    <div className="state-equation">
      <span className="eq-part">8-gene expression profile</span>
      <span className="eq-arrow">→</span>
      <span className="eq-part">molecular state</span>
      <span className="eq-arrow">→</span>
      <span className="eq-part eq-final">risk stratification axis</span>
    </div>
  </div>
);

/* ================================================================
   GenePanelTable — 8 유전자 + IHC feasibility
   ================================================================ */
export const GenePanelTable: React.FC = () => {
  const rows = [
    { sym: "TG",      alias: "thyroglobulin",      group: "effector", dir: "down", d: 0.86, tier: "moderate", role: "갑상선 호르몬 precursor · iodination 기질",  ihc: "★★★  임상 IHC 표준 marker (Tg)" },
    { sym: "TPO",     alias: "thyroperoxidase",    group: "effector", dir: "down", d: 2.30, tier: "strong",   role: "Iodide oxidation · organification 핵심 효소", ihc: "★★★  antibody 안정 · 임상 적용 흔함" },
    { sym: "TSHR",    alias: "TSH receptor",       group: "effector", dir: "down", d: 1.20, tier: "strong",   role: "TSH 신호 진입점 — 분화 신호 유지",            ihc: "★★    membrane staining · 표준화 필요" },
    { sym: "SLC5A5",  alias: "NIS",                group: "effector", dir: "down", d: 0.22, tier: "mild",     role: "I⁻ 흡수 symporter — RAI uptake 직접 결정",     ihc: "★★★  hMab 17.3 등 안정 · 임상 prognostic 보고 있음" },
    { sym: "DIO1",    alias: "type-1 deiodinase",  group: "effector", dir: "down", d: 1.24, tier: "strong",   role: "T4 → T3 변환 · 분화 metabolism",              ihc: "★★    antibody 가능하나 표준화 부족" },
    { sym: "PAX8",    alias: "—",                  group: "TF",       dir: "down", d: 0.97, tier: "strong",   role: "Thyroid lineage master TF — 분화 유지",       ihc: "★★★  병리 진단 routine marker" },
    { sym: "NKX2-1",  alias: "TTF-1",              group: "TF",       dir: "down", d: 0.63, tier: "moderate", role: "Thyroid · lung 분화 TF — 발달 핵심",          ihc: "★★★  병리 진단 routine marker" },
    { sym: "FOXE1",   alias: "TTF-2",              group: "TF",       dir: "down", d: 0.84, tier: "moderate", role: "Thyroid 분화 발달 TF",                         ihc: "★★    antibody 가능 · 임상 적용 적음" }
  ];
  return (
    <div>
      <table className="gene-table">
        <thead>
          <tr>
            <th>공식 HGNC</th>
            <th>Alias / 임상명</th>
            <th>Group</th>
            <th>DM1 발현 방향</th>
            <th>HM450 β  d (DM1−DM2)</th>
            <th>RAI biology 에서의 역할</th>
            <th>IHC feasibility</th>
          </tr>
        </thead>
        <tbody>
          {rows.map(r => (
            <tr key={r.sym} className={r.group === "TF" ? "row-tf" : "row-effector"}>
              <td><code>{r.sym}</code></td>
              <td>{r.alias}</td>
              <td>
                <span className={`group-pill ${r.group}`}>
                  {r.group === "TF" ? "TF (lineage)" : "Effector (uptake)"}
                </span>
              </td>
              <td>
                <span className="dir-cell down">
                  <span className="dir-arrow">▼</span> DOWN-regulated
                </span>
              </td>
              <td className={`d-cell tier-${r.tier}`}>
                <b>+{r.d.toFixed(2)}</b>
                <span className="d-tier">{r.tier}</span>
              </td>
              <td>{r.role}</td>
              <td>{r.ihc}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <div className="dir-table-note">
        ★ <b>모든 8 panel 유전자가 DM1 에서 DOWN-regulated (silenced)</b>.  Promoter β Cohen's d 가 클수록 DM1 에서 더 강한 hypermethylation → 더 강한 mRNA silencing.
        TPO (d = 2.30) 가 가장 강한 silencing, SLC5A5 (d = 0.22 promoter NS) 가 mRNA level 에서만 silencing.
        TCGA HM450 n = 503  ·  audit_2026_04_30/round5/r5_2_per_gene_methylation_DM.tsv.
      </div>
      <div className="gene-group-summary">
        <div className="ggs effector">
          <div className="ggs-tag">Effector 5</div>
          <div className="ggs-genes">TG · TPO · TSHR · SLC5A5 · DIO1</div>
          <div className="ggs-role">Iodine handling · 호르몬 합성 · RAI uptake 도구</div>
        </div>
        <div className="ggs tf">
          <div className="ggs-tag">TF 3</div>
          <div className="ggs-genes">PAX8 · NKX2-1 · FOXE1</div>
          <div className="ggs-role">Thyroid lineage / 분화 master regulator</div>
        </div>
      </div>
      <div className="warn-line">
        ⚠ 분석 / 코드 / 테이블에서는 <b>공식 HGNC symbol</b> 로 통일.  발표 / 병리 설명에서만 alias 병기.
      </div>
    </div>
  );
};

/* ================================================================
   KeyNumbersTable — Audit-locked 핵심 수치
   ================================================================ */
export const KeyNumbersTable: React.FC = () => {
  const rows = [
    { dom: "Discovery", cohort: "TCGA-THCA n = 504", result: "DM1 28.4 % / DM2 71.6 % (KMeans k=2)", interp: "발견 코호트 partition", caveat: "primary PTC only" },
    { dom: "Discovery", cohort: "TCGA top-5000 MAD",  result: "ARI = 0.92 vs 8-gene reference",         interp: "panel artifact 아님", caveat: "—" },
    { dom: "Driver",    cohort: "TCGA",                result: "BRAF · TERT · RAS single-feature AUC ≈ 0.5", interp: "driver-orthogonal axis", caveat: "single-feature 기준" },
    { dom: "Driver",    cohort: "TCGA BRAF V600E+",   result: "DM1 prevalence 0.7 %",                  interp: "BRAF tumor 거의 모두 DM2", caveat: "—" },
    { dom: "Driver",    cohort: "TCGA RAS-mutant",     result: "DM1 prevalence 96.4 %",                  interp: "RAS-route → dedifferentiation", caveat: "—" },
    { dom: "Driver",    cohort: "TCGA BRAF/RAS−",     result: "DM1 49.1 % / DM2 50.9 %",                interp: "dark matter — 분층화 가치 max", caveat: "n=173" },
    { dom: "Mechanism", cohort: "TCGA HM450 n = 503",  result: "TPO d = 2.30 · p = 1.9 × 10⁻¹⁸",        interp: "epigenetic silencing 강함", caveat: "promoter β 기준" },
    { dom: "Mechanism", cohort: "TCGA HM450",          result: "mean β  DM1 0.385 vs DM2 0.253 (+52 %)", interp: "8-gene 동시 silencing", caveat: "—" },
    { dom: "External",  cohort: "Master forest 14 ent × 11 coh", result: "mean d = 2.81 · median 2.37 · all ≥ 1.56", interp: "cross-cohort 강건", caveat: "—" },
    { dom: "External",  cohort: "Per-gene × cohort matrix", result: "80 / 80 cell direction-consistent", interp: "단일 유전자 의존성 없음", caveat: "—" },
    { dom: "External",  cohort: "Lee 2024 FFPE n = 632", result: "within-cohort KMeans d = 5.93",         interp: "Korean PTC 재현 강력", caveat: "—" },
    { dom: "External",  cohort: "Mun 2025 protein n = 336", result: "thyroid_diff d = −1.91 · 7 / 7 sign", interp: "단백질 수준 검증", caveat: "—" },
    { dom: "External",  cohort: "Landa 2016 PDTC + ATC", result: "5 / 8 panel overlap with silenced list", interp: "convergent biology (RC lock)", caveat: "독립 설계" },
    { dom: "Clinical",  cohort: "TCGA + MSK pooled OS", result: "HR 2.53 [1.31, 4.89] · I² = 0 %",        interp: "retrospective 후향 hazard",  caveat: "MSK 기여 큼 · primary alone underpowered" },
    { dom: "Clinical",  cohort: "GSE151179 post-RAI",   result: "thyroid_diff d ≈ −1.0 · p ≈ 10⁻⁴",      interp: "post-RAI ≈ DM1 transcriptional", caveat: "pre n=35 / post n=17" },
    { dom: "Deploy",    cohort: "FFPE vs FF",           result: "Kolmogorov-Smirnov p = 0.44",            interp: "FFPE-compatible assay",      caveat: "TCGA FF vs Lee FFPE" },
    { dom: "sc",        cohort: "Lu 2023 thyrocyte n = 14,624", result: "DM gradient 유지 (KRT8 ∩ KRT19 ∩ EPCAM)", interp: "thyrocyte-intrinsic 신호", caveat: "stromal/immune confound 아님" }
  ];
  return (
    <table className="numbers-table">
      <thead>
        <tr>
          <th>Domain</th><th>Cohort / Dataset</th><th>Result</th><th>Interpretation</th><th>Caveat</th>
        </tr>
      </thead>
      <tbody>
        {rows.map((r, i) => (
          <tr key={i} className={`dom-${r.dom.toLowerCase()}`}>
            <td><span className={`dom-pill dom-${r.dom.toLowerCase()}`}>{r.dom}</span></td>
            <td>{r.cohort}</td>
            <td><b>{r.result}</b></td>
            <td>{r.interp}</td>
            <td className="caveat-cell">{r.caveat}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
};

/* ================================================================
   IRBChecklist — IRB chart-review 변수 (그룹별)
   ================================================================ */
export const IRBChecklist: React.FC = () => {
  const groups = [
    { name: "Demographics",    vars: ["나이", "성별"] },
    { name: "Pathology",       vars: ["Histology (PTC · FTC · PDTC · ATC)", "Tumor size", "Lymph node metastasis", "Distant metastasis", "TNM stage"] },
    { name: "Molecular",       vars: ["BRAF V600E", "RAS (HRAS · KRAS · NRAS)", "TERT promoter", "RET / NTRK / ALK fusion", "NGS panel 전체 결과"] },
    { name: "Treatment",       vars: ["수술일", "RAI 시행 여부", "RAI dose", "Cumulative RAI dose", "Post-therapy WBS uptake", "Systemic / targeted therapy 기록"] },
    { name: "Biochemistry",    vars: ["Stimulated Tg", "TgAb"] },
    { name: "Outcomes",        vars: ["Recurrence", "Persistent disease", "Progression", "PFS", "Structural incomplete response"] },
    { name: "Toxicity ★",      vars: ["WBC 변화", "Neutrophil 변화", "Platelet 변화", "Leukopenia", "Marrow suppression", "기타 RAI 부작용"], emph: true }
  ];
  return (
    <div className="irb-grid">
      {groups.map(g => (
        <div key={g.name} className={`irb-group ${g.emph ? "irb-emph" : ""}`}>
          <div className="irb-title">{g.name}</div>
          <ul>{g.vars.map(v => <li key={v}>{v}</li>)}</ul>
          {g.emph && <div className="irb-emph-note">★ 강민수 선생님 메시지 — 누적 RAI 독성 secondary endpoint 로 추적</div>}
        </div>
      ))}
    </div>
  );
};

/* ================================================================
   ReviewerAttacksCard — top 5 reviewer 공격 + 답변
   ================================================================ */
export const ReviewerAttacksCard: React.FC = () => {
  const attacks = [
    {
      q: "왜 하필 이 8 개 gene 인가?",
      a: "Biology-prioritized · literature-guided · functionally grouped · cross-cohort validated · parsimony-driven.",
      detail: "TIERA67 67-gene 후보 풀 (Yoo 2016 prior) → 5 기능 범주로 분류 → 19 외부 코호트 재현성 → redundancy / parsimony → 최종 8 개. 16-gene · random 8-gene · pan-genome 과 비교."
    },
    {
      q: "정말 RAI response 를 예측하나?",
      a: "Definitive response prediction 이 아닌, RAI avidity / refractory-like biology 와 관련된 molecular state classifier 로 표현.",
      detail: "공개 cohort 의 RAI response label 이 일관되지 않으므로 \"directly predicts response\" 라고 쓰지 않음. GSE151179 post-RAI 종양이 transcriptionally DM1-aligned (d ≈ −1.0) 라는 점이 가장 가까운 직접 증거."
    },
    {
      q: "BRAF / RAS 일 뿐 아닌가?",
      a: "Canonical drivers 단일 feature AUC = chance level (~0.5).  BRAF / RAS-negative dark matter 가 DM1 49 % / DM2 51 % 로 거의 절반 분할.",
      detail: "Driver 변이로는 DM 분류가 불가능. Methylation 은 MAPK output 에 비례하지만 driver 정체성에는 비례하지 않음 (BRAF ≈ RET fusion ≫ RAS)."
    },
    {
      q: "Platform-specific 아닌가?",
      a: "RNA-seq · GPL570 microarray · HM450 methylation · proteogenomic · single-cell · FFPE 모두에서 일관.",
      detail: "Lee 2024 FFPE Korean d=5.93, Mun 2025 protein d=−1.91, GPL570 4 cohort 모두 ρ ≤ −0.84, FFPE vs FF KS p=0.44. Cross-platform robust."
    },
    {
      q: "임상에서 바로 쓸 수 있나?",
      a: "아직 아니다.  SNUBH NGS 호환성 확인 → retrospective validation → IHC translation → prospective validation 4 단계 필요.",
      detail: "현 단계는 retrospective public-cohort evidence + external validation 단계.  Treatment-selection biomarker 가 아닌 risk stratification axis 로만 기술."
    }
  ];
  return (
    <div className="attacks-grid">
      {attacks.map((a, i) => (
        <div key={i} className="attack-card">
          <div className="attack-num">⚠ {i+1}</div>
          <div className="attack-q">"{a.q}"</div>
          <div className="attack-a"><b>답변:</b> {a.a}</div>
          <div className="attack-detail">{a.detail}</div>
        </div>
      ))}
    </div>
  );
};

/* ================================================================
   MeetingQuestions — 강민수 선생님 / 병리 협업 질문 2 column
   ================================================================ */
export const MeetingQuestions: React.FC = () => (
  <div className="mq-grid">
    <div className="mq-col">
      <div className="mq-head">강민수 선생님께 — 임상 데이터 / NGS</div>
      <ol>
        <li>분당서울대병원 NGS panel 에 8 개 gene 이 각각 포함되어 있나요?<br/>
            <span className="mq-sub">TG · TPO · TSHR · SLC5A5 · DIO1 · PAX8 · NKX2-1 · FOXE1</span></li>
        <li>갑상선암 환자 중 NGS 시행 환자는 몇 명 정도 있나요?</li>
        <li>그 환자들의 <b>RAI 치료 이력</b> (시행 / cumulative dose / post-therapy WBS / Tg response / recurrence) 을 retrospective 로 추적 가능한가요?</li>
        <li>"RAI 불응성" 정의를 어떤 기준으로 잡는 게 임상적으로 가장 적절한가요?<br/>
            <span className="mq-sub">예: structural incomplete response · Tg ≥ 10 · post-therapy scan non-avid · distant met 진행</span></li>
        <li>약 30 명 internal validation 에서 가장 현실적인 endpoint 는?<br/>
            <span className="mq-sub">RAI uptake · stimulated Tg · structural incomplete response · recurrence · cumulative RAI dose</span></li>
        <li><b>★ 누적 RAI 독성</b> (WBC 감소 · leukopenia · marrow suppression) 을 chart 에서 retrieve 가능한가요?</li>
      </ol>
    </div>
    <div className="mq-col">
      <div className="mq-head">병리 협업 — IHC translation</div>
      <ol>
        <li>수술 후 보관 FFPE block 에서 8 개 marker IHC 가 가능한가요?</li>
        <li>8 개 중 실제 IHC antibody 가 임상 등급으로 안정적인 marker 는?<br/>
            <span className="mq-sub">우선 후보: PAX8 · NKX2-1 (TTF-1) · TG · TPO (병리 routine)</span></li>
        <li>병리 scoring 은 H-score · intensity 0/1/2/3 · positive cell % 중 어떤 방식이 적절한가요?</li>
        <li>Cancer cell-specific staining 과 stromal staining 을 분리해서 볼 수 있나요?</li>
        <li>환자 당 몇 slide 까지 가능하고, 비용 / 시간은?</li>
        <li>장기 비전: H&E + spatial transcriptomics AI 와 IHC 의 협업 가능성?</li>
      </ol>
    </div>
  </div>
);
