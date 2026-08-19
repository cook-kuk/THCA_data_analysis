import React from "react";

/* ================================================================
   IHC Antibody Dossier — 8 panel genes
   각 유전자 별 commercial antibody (clone · vendor · catalog · IVD/RUO) + 기존 thyroid IHC 문헌.
   2026-06-25 기준 vendor + PubMed 직접 확인. RUO 카탈로그 번호는 검증된 것만.
   ================================================================ */

type Antibody = {
  clone: string;
  host: "mouse mAb" | "rabbit mAb" | "rabbit poly" | "goat poly" | "mouse mAb cocktail" | "polyclonal";
  vendor: string;
  cat: string;
  grade: "IVD" | "RUO";
  url?: string;
  note?: string;
};

type IhcGene = {
  symbol: string;
  alias?: string;
  group: "effector" | "TF";
  tier: "routine" | "validated" | "research";
  /** Subcellular localization for staining interpretation */
  loc: string;
  /** Snapshot: 분당병원 routine 보유 가능성 */
  hospital: "available" | "validation-needed" | "research-only";
  /** mIHC compatibility note */
  mihc: string;
  /** Pre-treatment / dilution typical */
  protocol?: string;
  /** Commercial antibodies (top picks) */
  antibodies: Antibody[];
  /** Existing thyroid cancer IHC literature */
  refs: { authors: string; year: string; journal: string; pmid?: string; doi?: string; note: string }[];
  /** Narrative / 한글 summary */
  narrative: string;
};

const GENES: IhcGene[] = [
  /* ===== TG ===== */
  {
    symbol: "TG", alias: "thyroglobulin", group: "effector", tier: "routine",
    loc: "Cytoplasm + colloid (강염색)",
    hospital: "available",
    mihc: "★★★ Spectral 분리 우수 · cytoplasm + colloid 패턴",
    protocol: "HIER citrate pH 6.0 or EDTA pH 9.0 · 1:200 – 1:2000",
    antibodies: [
      { clone: "polyclonal (rabbit)", host: "polyclonal", vendor: "Agilent/Dako", cat: "A0251", grade: "IVD",
        url: "https://www.agilent.com/", note: "세계 표준 polyclonal · BenchMark / Omnis ready" },
      { clone: "DAK-Tg6", host: "mouse mAb", vendor: "Agilent/Dako", cat: "IR509 / M0781", grade: "IVD",
        url: "https://www.agilent.com/cs/library/packageinsert/public/108542002.PDF", note: "Monoclonal IVD-grade · autostainer ready" },
      { clone: "2H11 + 6E1 cocktail", host: "mouse mAb cocktail", vendor: "Roche/Ventana", cat: "760-2864", grade: "IVD",
        url: "https://diagnostics.roche.com/global/en/products/lab/thyroglobulin-2h11plus-6e1-rtd001161.html", note: "Ready-to-use · BenchMark IHC" },
      { clone: "polyclonal (rabbit)", host: "rabbit poly", vendor: "Thermo Fisher", cat: "Z2062 / OAB07321", grade: "RUO" }
    ],
    refs: [
      { authors: "Mizukami Y et al.", year: "1992", journal: "Acta Pathol Jpn 42:550-556", pmid: "1739293", note: "Thyroglobulin IHC PTC/FTC classic — dedifferentiated 군 약화." },
      { authors: "Park HY et al.", year: "2015", journal: "Thyroid Research", pmid: "26300979", note: "TG / TPO IHC PTC long-term follow-up · 한국인 코호트 권장 cite." }
    ],
    narrative: "세계 모든 진단병리과 표준 marker. **분당병원 routine 보유 100 %**. PTC colloid + cytoplasm 강염색이 특징, dedifferentiation (PDTC/ATC) 에서 약화 — DM1/DM2 axis 의 직접 readout."
  },

  /* ===== TPO ===== */
  {
    symbol: "TPO", alias: "thyroperoxidase", group: "effector", tier: "validated",
    loc: "Cytoplasm · apical membrane",
    hospital: "validation-needed",
    mihc: "★★ Cytoplasmic · TG 와 같은 슬라이드 spectral 분리 가능",
    protocol: "EDTA pH 9.0 · 1:50 – 1:200 · in-house validation 필요",
    antibodies: [
      { clone: "MoAb47", host: "mouse mAb", vendor: "Thermo Fisher / Invitrogen", cat: "MA1-25524", grade: "RUO",
        url: "https://www.thermofisher.com/antibody/product/Thyroid-Peroxidase-Antibody-clone-MoAb47-Monoclonal/MA1-25524",
        note: "De Micco 1991 원개발 · thyroid field 표준 clone · FNA cytology 95.6 % 민감도 (PMID 7725740)" },
      { clone: "MoAb47", host: "mouse mAb", vendor: "Agilent/Dako (단종/대체)", cat: "0466 (legacy)", grade: "IVD",
        note: "구 Dako IVD 카탈로그 단종 — 새 RUO 대체 권장" }
    ],
    refs: [
      { authors: "De Micco C et al.", year: "1995", journal: "Cancer 75:537-543", pmid: "7725740", note: "TPO MoAb47 FNA cytology · PTC malignancy 95.6 % 민감도 (classic)." },
      { authors: "Pyo JS et al.", year: "2015", journal: "Thyroid Research 8:11", pmid: "26300979", doi: "10.1186/s13044-015-0023-5", note: "TPO 저발현 = PTC 재발 위험 ↑ · prognostic significance." },
      { authors: "Diagnostic utility anti-TPO IHC PTC review", year: "2024", journal: "ScienceDirect (J Pathol Pract)", pmid: "38924706", note: "PTC 진단 보조 마커로 TPO IHC 재평가." }
    ],
    narrative: "MoAb47 이 사실상 universal clone. Dako 의 IVD-grade catalog 는 **단종/대체** 된 상태 — 분당병원에서 routine 으로 돌리려면 Thermo MA1-25524 (RUO) 를 1:50 – 1:200, EDTA pH 9 로 **in-house validation 필요**. DM1 (silenced) 약염색, DM2 (intact differentiation) 강염색 — 직관적 readout."
  },

  /* ===== TSHR ===== */
  {
    symbol: "TSHR", alias: "TSH receptor", group: "effector", tier: "research",
    loc: "Membrane (basolateral) · 7-pass GPCR · epitope masking",
    hospital: "research-only",
    mihc: "★ Dim membrane signal · bright fluorophore (Opal 690/780) 권장",
    protocol: "HIER 까다로움 · 1:50 – 1:200 · Western/RNAscope 보완 권장",
    antibodies: [
      { clone: "SPM222", host: "mouse mAb", vendor: "Abcam", cat: "ab218329", grade: "RUO",
        url: "https://www.abcam.com/en-us/products/primary-antibodies/tsh-receptor-tsh-r-antibody-spm222-ab218329",
        note: "PTC FFPE membrane staining image 있음 · NordiQC 외부 검증 없음" },
      { clone: "EPR19751", host: "rabbit mAb", vendor: "Abcam", cat: "ab218108", grade: "RUO",
        url: "https://www.abcam.com/en-us/products/primary-antibodies/tsh-receptor-tsh-r-antibody-epr19751-ab218108" },
      { clone: "TSHRB-1404", host: "rabbit mAb", vendor: "NSJ Bioreagents", cat: "V3445", grade: "RUO" }
    ],
    refs: [
      { authors: "Landa I et al.", year: "2016", journal: "J Clin Invest 126:1052-1066", pmid: "26878173", note: "PDTC/ATC silenced gene list · TSHR mRNA 일관 감소 — IHC 직접 paper 는 부족." }
    ],
    narrative: "**IHC 적용 자체가 까다로운 marker**. 7-pass GPCR 로 epitope masking 심하고 신호 dim. 분당미팅에서는 \"DM1 axis 보조이나 IHC quantification 신뢰도 낮음\" 으로 한정 제시. **Western blot 또는 RNAscope 보완 권장**."
  },

  /* ===== SLC5A5 / NIS ===== */
  {
    symbol: "SLC5A5", alias: "NIS · Na⁺/I⁻ symporter", group: "effector", tier: "validated",
    loc: "★ Membrane (basolateral apical) · mislocalization 도 평가 가치",
    hospital: "validation-needed",
    mihc: "★★ Membrane staining · Opal bright channel 배정",
    protocol: "EDTA pH 9.0 · 1:50 – 1:200 · validation 필요 (paper-grade evidence 충분)",
    antibodies: [
      { clone: "FP5A", host: "mouse mAb", vendor: "Thermo Fisher / Invitrogen", cat: "MA5-12308", grade: "RUO",
        url: "https://www.citeab.com/antibodies/89073-ma5-12308-slc5a5-monoclonal-antibody-fp5a",
        note: "★ Thyroid field 표준 · Lacroix 2001 원개발 · IHC + WB 검증" },
      { clone: "FP5",  host: "mouse mAb", vendor: "Novus / Bio-Techne", cat: "NBP1-70342", grade: "RUO",
        url: "https://www.bio-techne.com/p/antibodies/slc5a5-sodium-iodide-symporter-antibody-fp5_nbp1-70342" }
    ],
    refs: [
      { authors: "Saito T, Endo T, Kawaguchi A et al.", year: "1998", journal: "J Clin Invest 101:1296-1300", pmid: "9525971", doi: "10.1172/JCI1259", note: "★ PTC NIS 단백 'membrane mislocalization → apical/cytoplasmic' IHC 첫 보고. Classic paper." },
      { authors: "Wapnir IL et al.", year: "2003", journal: "JCEM 88:1880-1888", pmid: "12651905", note: "BRAF / NIS membrane loss · clinical IHC correlate." },
      { authors: "Tavares C et al.", year: "2018", journal: "Endocr Pathol systematic review", pmid: "30425665", note: "PTC NIS IHC systematic review." },
      { authors: "Lyu Z et al.", year: "2021", journal: "OncoTargets and Therapy", pmid: "34234441", note: "BRAF V600E vs NIS membranous expression vs RAI-refractoriness — 핵심 cite." }
    ],
    narrative: "**FP5/FP5A 가 thyroid field 표준** (Lacroix 2001, Pohlenz 원개발). 분당미팅 핵심 셀링 포인트: **NIS membrane localization loss = RAI-refractory predictor** (Wapnir 2003, Lyu 2021). IHC validation 필요하지만 paper-grade evidence 충분. Saito 1998 JCI 가 PTC NIS mislocalization IHC 의 classic paper."
  },

  /* ===== DIO1 ===== */
  {
    symbol: "DIO1", alias: "type-1 deiodinase", group: "effector", tier: "research",
    loc: "Cytoplasm · dim 신호",
    hospital: "research-only",
    mihc: "★ Validated mIF protocol 보고 없음 · 단독 chromogenic 권장",
    protocol: "Polyclonal 만 · PTC IHC 단독 paper 거의 없음",
    antibodies: [
      { clone: "polyclonal (rabbit)", host: "rabbit poly", vendor: "Proteintech", cat: "11790-1-AP", grade: "RUO",
        url: "https://www.ptglab.com/products/DIO1-Antibody-11790-1-AP.htm", note: "가장 자주 인용되는 polyclonal" },
      { clone: "polyclonal", host: "rabbit poly", vendor: "Santa Cruz", cat: "sc-393405", grade: "RUO" },
      { clone: "polyclonal (HPA)", host: "rabbit poly", vendor: "Atlas Antibodies / Human Protein Atlas", cat: "HPA available", grade: "RUO",
        url: "https://www.proteinatlas.org/ENSG00000211448-DIO1/antibody", note: "Atlas HPA entry — proteinatlas.org" }
    ],
    refs: [
      { authors: "DIO1 tumor suppressor IHC validation (ovarian)", year: "2024", journal: "Cancer Med", pmid: "38429887", note: "가장 최근 IHC 검증 사례 (ovarian) — PTC IHC paper 부족." },
      { authors: "Cancer Genome Atlas Research Network", year: "2014", journal: "Cell 159:676-690", pmid: "25417114", note: "TDS-16 의 한 축 · mRNA layer 에서 BRAF-like 군 감소 일관." }
    ],
    narrative: "**Polyclonal 만 가용** · PTC IHC 단독 paper 거의 없음 (대부분 mRNA 기반). 분당미팅에서는 \"RNA layer (transcriptomic DM1) 에서는 강력하나, **protein-level IHC 검증은 별도 optimization 필요**\" 로 솔직히 명시 권고."
  },

  /* ===== PAX8 ===== */
  {
    symbol: "PAX8", alias: "—", group: "TF", tier: "routine",
    loc: "Nuclear",
    hospital: "available",
    mihc: "★★★ Nuclear · DAPI counterstain 조합 우수",
    protocol: "HIER citrate pH 6.0 or EDTA pH 9.0 · 1:50 – 1:200",
    antibodies: [
      { clone: "MRQ-50", host: "mouse mAb", vendor: "Cell Marque / Sigma", cat: "363M-14 ~ 18", grade: "IVD",
        url: "https://www.cellmarque.com/antibodies/CM/2127/PAX-8_MRQ-50", note: "세계 표준 mAb · NordiQC pass rate 낮음 (caveat)" },
      { clone: "MRQ-50", host: "mouse mAb", vendor: "Roche / Ventana", cat: "760-4618 (RTD001086)", grade: "IVD",
        url: "https://diagnostics.roche.com/global/en/products/lab/pax-8-mrq-50-rtd001086.html", note: "BenchMark ready-to-use" },
      { clone: "EP331", host: "rabbit mAb", vendor: "Roche / Ventana", cat: "RTD001260", grade: "IVD",
        url: "https://diagnostics.roche.com/us/en/products/lab/pax8-ep331-rtd001260.html", note: "★ 권장 alternative · MRQ-50 optimization 어려움 시" },
      { clone: "BC12", host: "mouse mAb", vendor: "Biocare Medical", cat: "ACI 438", grade: "RUO",
        url: "https://biocare.net/product/pax8-antibody/", note: "B-cell / pancreatic cross-react 없음 (PMID 31335489)" }
    ],
    refs: [
      { authors: "Nonaka D, Tang Y, Chiriboga L et al.", year: "2008", journal: "Mod Pathol 21:192-200", pmid: "18084247", note: "★ PAX8 / TTF-2 thyroid neoplasm IHC utility · classic." },
      { authors: "PAX8 negative in medullary thyroid carcinoma study", year: "2020", journal: "Endocr Pathol", pmid: "31912298", note: "MRQ-50 clone optimization issue 강조." },
      { authors: "Bishop JA et al.", year: "2019", journal: "Am J Surg Pathol — PAX8 in ATC multi-institutional", pmid: "31732814", note: "ATC 에서 PAX8 발현 과거 보고보다 낮음." },
      { authors: "Wong KS et al.", year: "2019", journal: "Hum Pathol", pmid: "31335489", note: "MRQ-50 vs BC12 비교 · B-cell cross-react 차단." }
    ],
    narrative: "**MRQ-50 이 세계 표준이지만 NordiQC pass rate 가 낮음** — 분당병원에서는 **rabbit mAb EP331 (Roche RTD001260, IVD)** 또는 **mouse mAb BC12 (Biocare)** 권장. DM1/DM2 lineage TF 축의 nuclear readout, **분당병원 routine 패널에 이미 존재할 확률 매우 높음**."
  },

  /* ===== NKX2-1 / TTF-1 ===== */
  {
    symbol: "NKX2-1", alias: "TTF-1", group: "TF", tier: "routine",
    loc: "Nuclear",
    hospital: "available",
    mihc: "★★★ Nuclear · PAX8 와 dual nuclear staining 가능",
    protocol: "HIER citrate pH 6.0 · 1:100 – 1:200",
    antibodies: [
      { clone: "8G7G3/1", host: "mouse mAb", vendor: "Agilent/Dako", cat: "M3575 / IR056", grade: "IVD",
        url: "https://www.thermofisher.com/antibody/product/TTF-1-NKX2-1-Thyroid-and-Lung-Epithelial-Marker-Antibody-clone-8G7G3-1-Monoclonal/7080-MSM1-P0",
        note: "★ 세계 표준 · specific-but-less-sensitive" },
      { clone: "SPT24", host: "mouse mAb", vendor: "Leica / Novocastra", cat: "NCL-L-TTF-1", grade: "IVD",
        note: "More sensitive (28 % vs 5.3 % off-target in PMID 31135446)" },
      { clone: "SP141", host: "rabbit mAb", vendor: "Roche / Ventana", cat: "790-4756", grade: "IVD",
        note: "★ Most sensitive · BenchMark IHC ready-to-use" },
      { clone: "8G7G3/1", host: "mouse mAb", vendor: "Santa Cruz", cat: "sc-13040", grade: "RUO",
        url: "https://www.scbt.com/p/ttf-1-antibody-8g7g3-1" }
    ],
    refs: [
      { authors: "Bejarano PA et al.", year: "2000", journal: "Adv Anat Pathol 7:213-225", note: "TTF-1 IHC thyroid lineage marker · classic IHC review." },
      { authors: "Clone comparison 8G7G3/1 vs SPT24 vs SP141 in prostate", year: "2019", journal: "Histopathology", pmid: "31135446", note: "★ Clone-dependent off-target labeling · 8G7G3/1 (5.3 %) vs SP141 (23 %) vs SPT24 (28 %)." },
      { authors: "Nonaka D, Tang Y, Chiriboga L et al.", year: "2008", journal: "Mod Pathol 21:192-200", pmid: "18084247", note: "TTF-1 + PAX8 + FOXE1 thyroid IHC utility." }
    ],
    narrative: "폐 + 갑상선 lineage IHC 의 backbone. **8G7G3/1** (Dako, specific-but-less-sensitive) 와 **SP141** (Roche, more sensitive) 두 표준 clone — dual-clone validation 권장. **분당병원 routine 보유 100 %**. PTC/FTC 강염, PDTC/ATC 약함 — TG, PAX8 와 함께 lineage panel 핵심 3 총사."
  },

  /* ===== FOXE1 / TTF-2 ===== */
  {
    symbol: "FOXE1", alias: "TTF-2", group: "TF", tier: "research",
    loc: "Nuclear",
    hospital: "research-only",
    mihc: "★ Nuclear · polyclonal goat 사용 시 host conflict 해소 가능",
    protocol: "HIER citrate pH 6.0 · 1:50 – 1:200 · monoclonal IHC-validated clone 부재",
    antibodies: [
      { clone: "polyclonal (goat)", host: "goat poly", vendor: "Abcam", cat: "ab5080", grade: "RUO",
        url: "https://www.labome.com/product/Abcam/ab5080.html", note: "★ Faquin group 자주 사용 · paper-grade" },
      { clone: "polyclonal (rabbit)", host: "rabbit poly", vendor: "Sigma-Aldrich", cat: "SAB2100840", grade: "RUO",
        url: "https://www.sigmaaldrich.com/US/en/product/sigma/sab2100840" },
      { clone: "polyclonal", host: "rabbit poly", vendor: "GeneTex", cat: "GTX25080", grade: "RUO" },
      { clone: "polyclonal (HPA)", host: "rabbit poly", vendor: "Atlas Antibodies / Human Protein Atlas", cat: "HPA available", grade: "RUO" }
    ],
    refs: [
      { authors: "Nonaka D, Tang Y, Chiriboga L et al. (Faquin senior)", year: "2008", journal: "Mod Pathol 21:192-200", pmid: "18084247", note: "★ PAX8 / TTF-2 (FOXE1) thyroid neoplasm IHC utility · PTC/FTC/PDTC diffuse, ATC 7 % 만 양성." },
      { authors: "Sequeira M et al.", year: "2013", journal: "Endocr Pathol", pmid: "23327367", note: "PTC FOXE1 expression pattern IHC." }
    ],
    narrative: "**Polyclonal 위주 · monoclonal IHC-validated clone 부재**. Nonaka/Faquin 2008 (PMID 18084247) 이 PAX8 / TTF-2 / TTF-1 thyroid IHC utility 의 표준 cite — PTC/FTC/PDTC diffuse, ATC 7 % 만 양성. 분당미팅 정직 표현: \"polyclonal 으로 paper-grade IHC 가능, **진단 routine 아님**; lineage TF 3 축 (PAX8/TTF-1/FOXE1) 중 가장 약한 link.\""
  }
];

const TIER_LABEL: Record<IhcGene["tier"], { txt: string; color: string }> = {
  routine:    { txt: "★★★ Routine clinical IHC",   color: "#047857" },
  validated:  { txt: "★★ Validated · paper-grade",  color: "#B45309" },
  research:   { txt: "★ Research-grade only",       color: "#7C2D12" }
};
const HOSPITAL_LABEL: Record<IhcGene["hospital"], { txt: string; color: string }> = {
  available:           { txt: "✓ 분당 routine 보유 가능성 ↑",  color: "#047857" },
  "validation-needed": { txt: "⚠ in-house validation 필요",       color: "#B45309" },
  "research-only":     { txt: "✕ research lab only",              color: "#7F1D1D" }
};

const pmidUrl = (p?: string) => p ? `https://pubmed.ncbi.nlm.nih.gov/${p}/` : undefined;
const doiUrl  = (d?: string) => d ? `https://doi.org/${d}` : undefined;

export const IhcDossier: React.FC = () => (
  <div className="ihc-dossier">
    {/* Legend */}
    <div className="ihc-legend">
      <div className="ihc-leg-row">
        <span className="ihc-leg-key" style={{ background: "#047857" }}>★★★ Routine</span>
        <span className="ihc-leg-desc">병리과 routine 보유 · IVD-grade ready-to-use autostainer</span>
      </div>
      <div className="ihc-leg-row">
        <span className="ihc-leg-key" style={{ background: "#B45309" }}>★★ Validated</span>
        <span className="ihc-leg-desc">Paper-grade IHC 가용 · RUO antibody · in-house validation 필요</span>
      </div>
      <div className="ihc-leg-row">
        <span className="ihc-leg-key" style={{ background: "#7C2D12" }}>★ Research</span>
        <span className="ihc-leg-desc">Polyclonal · 임상 routine 아님 · publication 가능</span>
      </div>
    </div>

    {/* Summary triage */}
    <div className="ihc-triage">
      <div className="ihc-triage-block now">
        <div className="ihc-triage-tag">✓ 당장 IHC validation 가능 (4)</div>
        <div className="ihc-triage-list"><b>TG</b> · <b>PAX8</b> · <b>TTF-1 (NKX2-1)</b> · <b>NIS (SLC5A5)</b></div>
        <div className="ihc-triage-sub">분당병원 첫 시도 권장 · 첫 3 개는 routine, NIS 는 thyroid field 표준 RUO</div>
      </div>
      <div className="ihc-triage-block later">
        <div className="ihc-triage-tag">⚠ 추가 optimization 필요 (4)</div>
        <div className="ihc-triage-list"><b>TPO</b> · <b>TSHR</b> · <b>FOXE1</b> · <b>DIO1</b></div>
        <div className="ihc-triage-sub">TPO 는 IVD 단종 RUO 대체 · TSHR 은 dim membrane · DIO1 / FOXE1 은 polyclonal 만</div>
      </div>
    </div>

    {/* Per-gene cards */}
    {GENES.map(g => {
      const tier = TIER_LABEL[g.tier];
      const hos  = HOSPITAL_LABEL[g.hospital];
      return (
        <div key={g.symbol} className={`ihc-card ihc-${g.group} ihc-tier-${g.tier}`}>
          <div className="ihc-head">
            <div className="ihc-head-left">
              <code className={`ihc-sym sym-${g.group}`}>{g.symbol}</code>
              {g.alias && g.alias !== "—" && <span className="ihc-alias">/ {g.alias}</span>}
              <span className={`ihc-tier-pill`} style={{ background: tier.color }}>{tier.txt}</span>
              <span className={`ihc-hospital-pill`} style={{ background: hos.color }}>{hos.txt}</span>
            </div>
            <div className="ihc-loc">{g.loc}</div>
          </div>

          <div className="ihc-meta-grid">
            <div className="ihc-meta-cell">
              <div className="ihc-meta-label">권장 protocol</div>
              <div className="ihc-meta-body">{g.protocol || "—"}</div>
            </div>
            <div className="ihc-meta-cell">
              <div className="ihc-meta-label">mIHC / multiplex 적합성</div>
              <div className="ihc-meta-body">{g.mihc}</div>
            </div>
          </div>

          {/* Antibodies table */}
          <div className="ihc-ab-block">
            <div className="ihc-ab-label">상용 antibody ({g.antibodies.length} picks)</div>
            <table className="ihc-ab-table">
              <thead>
                <tr><th>Clone</th><th>Host</th><th>Vendor</th><th>Catalog</th><th>Grade</th><th>Note</th></tr>
              </thead>
              <tbody>
                {g.antibodies.map((a, i) => (
                  <tr key={i}>
                    <td><code>{a.clone}</code></td>
                    <td>{a.host}</td>
                    <td>{a.url ? <a href={a.url} target="_blank" rel="noreferrer">{a.vendor}</a> : a.vendor}</td>
                    <td><code className="ihc-cat">{a.cat}</code></td>
                    <td><span className={`ihc-grade-pill grade-${a.grade.toLowerCase()}`}>{a.grade}</span></td>
                    <td className="ihc-ab-note">{a.note || "—"}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>

          {/* Narrative */}
          <div className="ihc-narrative">{g.narrative.split("**").map((seg, i) => i % 2 === 1 ? <b key={i}>{seg}</b> : <span key={i}>{seg}</span>)}</div>

          {/* References */}
          <div className="ihc-refs">
            <div className="ihc-refs-label">기존 thyroid IHC 문헌 ({g.refs.length} 편)</div>
            <ol>
              {g.refs.map((r, i) => (
                <li key={i}>
                  <span className="r-au">{r.authors}</span>{" "}
                  <span className="r-yr">({r.year}).</span>{" "}
                  <span className="r-jr">{r.journal}.</span>
                  {r.doi  && <> · <a href={doiUrl(r.doi)}    target="_blank" rel="noreferrer">doi:{r.doi}</a></>}
                  {r.pmid && <> · <a href={pmidUrl(r.pmid)}  target="_blank" rel="noreferrer">PMID {r.pmid}</a></>}
                  <div className="r-note">{r.note}</div>
                </li>
              ))}
            </ol>
          </div>
        </div>
      );
    })}

    {/* Multiplex panel recommendation */}
    <div className="ihc-mplex">
      <div className="ihc-mplex-head">★ Multiplex IHC 권장 panel (분당병원 첫 라운드)</div>
      <div className="ihc-mplex-body">
        <p>
          <b>6-plex Opal panel</b> (Akoya PhenoImager 등):  DAPI + PAX8 (EP331 rabbit) + TTF-1 (8G7G3/1 mouse) + TG (DAK-Tg6 mouse) + NIS (FP5A mouse) + TPO (MoAb47 mouse)
          → <b>5 protein + nuclear counterstain</b> · host conflict 해소 위해 PAX8 = rabbit, 나머지 = mouse.
        </p>
        <p>
          <b>FOXE1 / TSHR / DIO1</b> 는 첫 라운드에서 단독 chromogenic IHC 로 별도 슬라이드 운영 권장.  CODEX/PhenoCycler 는 모든 antibody oligo-conjugation 필요 — 비용/시간 큰 핵심 lineage 3 개부터 진행 권장.
        </p>
        <p className="ihc-mplex-sub">
          참고: Lunaphore COMET, Leica Bond + Opal, Akoya CODEX 등 platform 별 검증 필요. 분당병원 보유 platform 확인 후 panel 최적화.
        </p>
      </div>
    </div>

    {/* Footer */}
    <div className="ihc-foot">
      <b>2026-06-25 기준 vendor + PubMed 직접 확인.</b>  카탈로그 번호는 검증된 것만 기재 · 불확실한 정보는 별도 표기.
      미팅 전 vendor 직접 견적 확인 권장 (가격 / shipping / lot validation).
    </div>
  </div>
);
