import React from "react";

/**
 * IhcHandoff — 분당병원 병리과 handoff document.
 * IHC 3-plex protocol (TG + PAX8 + TTF-1) + pilot 30 명 IRB fast-track template.
 */

export const IhcHandoff: React.FC = () => (
  <div className="ihchandoff">

    {/* Cover */}
    <div className="ihh-cover">
      <div className="ihh-cover-badge">병리과 즉시 실행 가능 · IRB fast-track template</div>
      <div className="ihh-cover-title">IHC 3-plex (TG + PAX8 + TTF-1) BRAF+ Pilot Protocol</div>
      <div className="ihh-cover-sub">
        분당병원 병리과 · 이미 routine 시행중 3 antibody + BRAF genotyping 조합 → DM1 axis 재현 검증<br/>
        Target: <b>Retrospective FFPE 30 명 pilot (BRAF+ 20 · BRAF− 10)</b> · 6 개월 완료 목표
      </div>
    </div>

    {/* Section 1 · Antibody specification table */}
    <div className="ihh-block">
      <div className="ihh-block-num">§1</div>
      <div className="ihh-block-title">항체 사양 · vendor · catalog · dilution</div>
      <table className="ihh-tbl">
        <thead>
          <tr>
            <th>Marker</th><th>Clone</th><th>Vendor · Catalog</th><th>Grade</th>
            <th>Dilution</th><th>Retrieval</th><th>Detection</th><th>Positive Control</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td><b>TG</b> (Thyroglobulin)</td>
            <td>2H11 + 6E1<br/>(cocktail)</td>
            <td>Dako / Agilent<br/>A0251 (polyclonal)<br/>또는 M0781</td>
            <td>IVD</td>
            <td>1:500 - 1:1000<br/>polyclonal</td>
            <td>Citrate pH 6.0<br/>HIER 20 min</td>
            <td>DAB · Envision FLEX<br/>Hematoxylin counter</td>
            <td>Normal thyroid<br/>(cytoplasmic + colloid)</td>
          </tr>
          <tr>
            <td><b>PAX8</b></td>
            <td>MRQ-50</td>
            <td>Roche / Cell Marque<br/>760-4638</td>
            <td>IVD (RUO some markets)</td>
            <td>Ready-to-use<br/>또는 1:100</td>
            <td>EDTA pH 8.0<br/>HIER 30 min</td>
            <td>DAB · UltraView<br/>Hematoxylin counter</td>
            <td>Normal thyroid<br/>(nuclear only)</td>
          </tr>
          <tr>
            <td><b>TTF-1 / NKX2-1</b></td>
            <td>8G7G3/1</td>
            <td>Roche / Ventana<br/>790-4398<br/>또는 Dako M3575</td>
            <td>IVD</td>
            <td>Ready-to-use<br/>또는 1:200</td>
            <td>EDTA pH 8.0<br/>HIER 30 min</td>
            <td>DAB · OptiView<br/>Hematoxylin counter</td>
            <td>Normal lung Type II<br/>(nuclear)</td>
          </tr>
        </tbody>
      </table>
      <div className="ihh-note">
        <b>Note:</b> 위 3 항체는 <b>분당병원 병리과가 이미 갑상선암 diagnosis 에 routine 시행</b> — 새 reagent 도입 불필요.
        Dako A0251 / Roche 760-4638 / 790-4398 은 대부분 국내 대학병원 공통 stock.  이 protocol 은 <b>기존 SOP 재활용</b>.
      </div>
    </div>

    {/* Section 2 · Scoring template */}
    <div className="ihh-block">
      <div className="ihh-block-num">§2</div>
      <div className="ihh-block-title">Scoring 프로토콜 · H-score + % positive</div>
      <div className="ihh-2col">
        <div>
          <h4>Primary scoring · H-score (0-300)</h4>
          <p className="ihh-p">
            H-score = <b>(1 × %weak) + (2 × %moderate) + (3 × %strong)</b>
          </p>
          <ul className="ihh-ul">
            <li>Weak (1+): 광학현미경 400× 에서 배경 대비 겨우 인식</li>
            <li>Moderate (2+): 명확한 signal, subcellular 위치 인식</li>
            <li>Strong (3+): 강한 signal, 즉시 판독</li>
            <li>Sum: 0 (전혀 없음) ~ 300 (모두 강한 3+)</li>
          </ul>
        </div>
        <div>
          <h4>Secondary · % Positive tumor cells</h4>
          <p className="ihh-p">
            <b>Tumor cell 중 positive 비율 (0-100 %)</b> — H-score 와 병기 시 robust.
          </p>
          <ul className="ihh-ul">
            <li>TG: cytoplasmic + colloid 두 위치</li>
            <li>PAX8: nuclear only (cytoplasmic 은 background)</li>
            <li>TTF-1: nuclear only</li>
            <li>Tumor area ≥ 200 cells 계수 · Peripheral necrotic zone 제외</li>
          </ul>
        </div>
      </div>

      <div className="ihh-tier-table">
        <h4>3-tier binning (survival KM 용)</h4>
        <table className="ihh-tbl">
          <thead>
            <tr><th>Tier</th><th>H-score (mean of 3)</th><th>DM1 axis 해석</th><th>임상 함의</th></tr>
          </thead>
          <tbody>
            <tr className="ihh-tier-lo"><td><b>Low (DM1-like)</b></td><td>0 - 100</td><td>Iodine-handling-low · lineage silenced</td><td>고용량 RAI 반복 이득 낮음 후보</td></tr>
            <tr className="ihh-tier-md"><td><b>Intermediate</b></td><td>101 - 200</td><td>Partial silencing</td><td>보수적 관찰 · fusion testing 우선</td></tr>
            <tr className="ihh-tier-hi"><td><b>High (DM2-like)</b></td><td>201 - 300</td><td>Iodine-handling-high · differentiated</td><td>표준 RAI 프로토콜 적용</td></tr>
          </tbody>
        </table>
      </div>
    </div>

    {/* Section 3 · QC + pathologist concordance */}
    <div className="ihh-block">
      <div className="ihh-block-num">§3</div>
      <div className="ihh-block-title">Quality Control · 2-pathologist blind + concordance</div>
      <ol className="ihh-ol">
        <li><b>Blind readout</b> — 병리과 senior 2 명이 임상 정보 blind 로 판독. 각자 H-score 기록.</li>
        <li><b>Concordance target</b> — Cohen's κ ≥ 0.7 (substantial agreement) · Weighted κ 계산.</li>
        <li><b>Discordant case</b> — κ &lt; 0.7 시 3rd 판독자 tie-break, 또는 consensus review.</li>
        <li><b>Batch effect control</b> — 모든 slide 를 하나의 IHC run 에서 processing (staining batch bias 방지).</li>
        <li><b>Positive/Negative controls</b> — Normal thyroid tissue (positive) · Non-thyroid tumor (negative) 각 run.</li>
        <li><b>Digital pathology</b> — QuPath 또는 Aperio HALO 자동 quantification 병기 (semi-quantitative validation).</li>
      </ol>
    </div>

    {/* Section 4 · Pilot cohort design */}
    <div className="ihh-block">
      <div className="ihh-block-num">§4</div>
      <div className="ihh-block-title">Pilot cohort design · Retrospective 30 명 (BRAF+ 20 · BRAF− 10)</div>
      <table className="ihh-tbl">
        <thead>
          <tr><th>Criteria</th><th>포함 기준</th><th>제외 기준</th></tr>
        </thead>
        <tbody>
          <tr>
            <td>Diagnosis</td>
            <td>Classical PTC 또는 FVPTC · 수술적 완전 절제 · 병리 확진</td>
            <td>Anaplastic thyroid cancer · Medullary thyroid cancer · Follicular carcinoma</td>
          </tr>
          <tr>
            <td>Molecular</td>
            <td>BRAF V600E genotyping 완료 (routine) · TSO500 v2 결과 available</td>
            <td>Driver mutation unknown · Sequencing failure</td>
          </tr>
          <tr>
            <td>Tissue</td>
            <td>Primary tumor FFPE block 5 μm section available · 최근 5 년 이내</td>
            <td>Recurrence-only sample · Poor FFPE quality · &lt; 200 tumor cells</td>
          </tr>
          <tr>
            <td>Follow-up</td>
            <td>수술 후 ≥ 2 년 follow-up · RAI treatment / response documented</td>
            <td>수술 후 6 개월 이내 사망 (RAI outcome 불가)</td>
          </tr>
          <tr>
            <td>Endpoint</td>
            <td>Recurrence · Progression · RAI refractory (ATA 2015 criteria)</td>
            <td>—</td>
          </tr>
        </tbody>
      </table>
      <div className="ihh-note">
        <b>Sampling strategy:</b> BRAF+ 20 명 (DM1-high recurrent 10 + DM2-like non-recurrent 10) + BRAF− 10 명 (control) = <b>30 명</b>.
        Extreme sampling 으로 예비 effect size 최대화 (validation cohort 는 balanced sampling 별도 계획).
      </div>
    </div>

    {/* Section 5 · IRB fast-track template */}
    <div className="ihh-block">
      <div className="ihh-block-num">§5</div>
      <div className="ihh-block-title">IRB Fast-track Template · 8 항목</div>
      <ol className="ihh-ol">
        <li><b>연구 제목:</b> "갑상선유두암 IHC 3-plex 재분석과 BRAF V600E-mutant tumor 의 radioiodine refractoriness 예측 관련성 후향 연구"</li>
        <li><b>연구 목적:</b> 이미 routine 시행중인 TG · PAX8 · TTF-1 IHC 데이터를 후향 재분석하여 BRAF+ subset 에서 무진행생존기간 (PFI) 분리 여부 검증</li>
        <li><b>대상 및 방법:</b> 최근 5 년 분당병원 갑상선유두암 수술 30 명 · FFPE block 재사용 · IHC 재판독 (blind H-score) · Cox regression</li>
        <li><b>동의 면제 근거:</b>
          <ul>
            <li>기존 FFPE 조직 재사용 (환자 추가 procedure 없음)</li>
            <li>연구자 개인정보 접근 없음 (병리 코드 익명화)</li>
            <li>후향 관찰 · 임상 결정에 영향 없음</li>
            <li>사망 · 소재불명 환자 위주 (동의 획득 실무 불가)</li>
          </ul>
        </li>
        <li><b>예상 위험:</b> Minimal — 조직 재판독만, 환자 직접 접촉 없음</li>
        <li><b>예상 이익:</b> 미래 갑상선암 환자의 불필요한 고용량 RAI 회피 근거 확보</li>
        <li><b>데이터 관리:</b>
          <ul>
            <li>병리 코드 → 익명 ID 변환 (양방향 lookup key 는 병리과장만 보관)</li>
            <li>Excel 파일 병원 secure server 보관 · 외부 export 금지</li>
            <li>분석 종료 후 5 년 보관 → 폐기</li>
          </ul>
        </li>
        <li><b>Timeline:</b> IRB 승인 → 병리 재판독 4 주 → 분석 4 주 → 결과 보고 8 주 = <b>총 4 개월</b></li>
      </ol>
    </div>

    {/* Section 6 · Statistical analysis plan */}
    <div className="ihh-block">
      <div className="ihh-block-num">§6</div>
      <div className="ihh-block-title">Statistical Analysis Plan (SAP)</div>
      <div className="ihh-2col">
        <div>
          <h4>Primary hypothesis</h4>
          <p className="ihh-p">
            BRAF V600E+ subset (n = 20) 에서 IHC 3-plex mean H-score tertile 이 PFI 를 유의하게 stratify.
          </p>
          <h4>Primary endpoint</h4>
          <ul className="ihh-ul">
            <li>PFI (수술 → 재발 or 진행) log-rank test, tertile split</li>
            <li>Cox HR per +1 SD H-score (BRAF+ subset)</li>
            <li>Target: p &lt; 0.05, HR &gt; 1.5 or &lt; 0.67</li>
          </ul>
        </div>
        <div>
          <h4>Secondary endpoints</h4>
          <ul className="ihh-ul">
            <li>RAI-refractory 임상 판정 vs H-score correlation</li>
            <li>Individual marker (TG, PAX8, TTF-1) contribution decomposition</li>
            <li>Pathologist concordance Cohen's κ</li>
            <li>Digital pathology (QuPath) automated vs manual concordance</li>
          </ul>
          <h4>Sensitivity analyses</h4>
          <ul className="ihh-ul">
            <li>Age adjusted Cox</li>
            <li>Stage III/IV subset 분석</li>
            <li>Follow-up ≥ 3 년 restrict</li>
          </ul>
        </div>
      </div>
    </div>

    {/* Section 7 · Cost + timeline */}
    <div className="ihh-block">
      <div className="ihh-block-num">§7</div>
      <div className="ihh-block-title">비용 · 일정 · 인력</div>
      <table className="ihh-tbl">
        <thead>
          <tr><th>항목</th><th>세부</th><th>비용 (KRW)</th><th>담당</th><th>기간</th></tr>
        </thead>
        <tbody>
          <tr><td>IHC 재판독</td><td>30 slide × 3 marker = 90 slide (기존 stain 재사용)</td><td>0 (재사용)</td><td>병리 senior 2 명</td><td>4 주</td></tr>
          <tr><td>병리 판독료 (blind)</td><td>2 pathologist × 30 case × 3 marker</td><td>~ 900,000</td><td>병리과</td><td>4 주</td></tr>
          <tr><td>Digital pathology</td><td>QuPath open-source software</td><td>0</td><td>연구원 1 명</td><td>2 주</td></tr>
          <tr><td>통계 분석</td><td>Cox regression · κ · KM</td><td>0 (내부)</td><td>Cook (분석)</td><td>2 주</td></tr>
          <tr><td>IRB submission</td><td>Fast-track template</td><td>0</td><td>연구 coordinator</td><td>2 주</td></tr>
          <tr><td><b>합계</b></td><td></td><td><b>~ 900,000 KRW</b></td><td></td><td><b>총 4 개월</b></td></tr>
        </tbody>
      </table>
      <div className="ihh-note ihh-note-hi">
        <b>결정적 이점:</b> 새 assay 개발 · 새 reagent 구매 · 새 pathologist 훈련 <b>모두 불필요</b>.
        기존 SOP 로 즉시 착수 가능 → <b>비용 &lt; 1M KRW · 4 개월 completion</b>.  이후 분당 prospective (n = 200) 착수 근거 확보.
      </div>
    </div>

    {/* Section 8 · Escalation path */}
    <div className="ihh-block">
      <div className="ihh-block-num">§8</div>
      <div className="ihh-block-title">Escalation Path · Retrospective → Prospective</div>
      <div className="ihh-escalation">
        <div className="ihh-step">
          <div className="ihh-step-num">Step 1</div>
          <div className="ihh-step-body">
            <b>Pilot 30 명 (본 protocol)</b> · Retrospective FFPE · 6 개월 · 비용 &lt; 1M KRW.
            성공 기준: BRAF+ subset log-rank p &lt; 0.05 · κ ≥ 0.7.
          </div>
        </div>
        <div className="ihh-arrow">↓</div>
        <div className="ihh-step">
          <div className="ihh-step-num">Step 2</div>
          <div className="ihh-step-body">
            <b>Prospective 200 명</b> · A5 SNUBH simulation power ≥ 90% 확인 · 2 년 등록 + 5 년 follow-up · 비용 ~ 20M KRW (staining routine + follow-up 조사).
            IRB 별도 정식 신청.
          </div>
        </div>
        <div className="ihh-arrow">↓</div>
        <div className="ihh-step">
          <div className="ihh-step-num">Step 3</div>
          <div className="ihh-step-body">
            <b>다기관 확장</b> · Aju · Samsung · Asan · Severance 협업 · n ~ 1000 · 5 년.
            Guideline 등재 논의 · KFDA / ATA committee submission.
          </div>
        </div>
      </div>
    </div>

    <div className="ihh-footer">
      본 handoff document 는 <b>강민수 선생님 (분당병원 병리과) · 유형원 선생님 (외과) · Cook 국승호 (연구팀)</b> 공동 실행.
      2026-06-25 미팅 후 <b>2 주 이내 IRB Fast-track submission target</b>.  IRB 승인 후 병리 판독 즉시 시작.
    </div>
  </div>
);
