import React from "react";

/**
 * Foundational biology reference cards — 10 landmark papers on which the DM1/8-gene axis rests.
 *
 * NOTE: Original figures are © respective publishers and are NOT reproduced here.  Each card
 * summarizes the key concept in our own words and provides a DOI / PubMed link to the original.
 * Where useful, we provide our own SVG concept sketch (not derived from the paper's figure).
 */

type Foundation = {
  tier: string;
  topic: string;
  citation: string;
  journal: string;
  key_finding: string;
  fig_desc: string;
  supports: string;
  pmid?: string;
  doi?: string;
  extra?: string;
};

const FOUNDATIONS: Foundation[] = [
  {
    tier: "FOUNDATION 1",
    topic: "Thyrocyte anatomy · iodine handling circuit",
    citation: "Carvalho DP, Dupuy C.  (2017)",
    journal: "Molecular and Cellular Endocrinology  ·  458 : 6-15",
    key_finding: "Thyroid follicular cell iodine transport circuit — TSH → TSHR → cAMP → NIS (basolateral) → intracellular I⁻ → Pendrin (apical) → TPO iodinates Tg tyrosines → MIT/DIT → T3/T4.",
    fig_desc: "Figure 1 · schematic of thyrocyte polarized architecture (colloid on apical side, blood on basolateral side).",
    supports: "우리 8-gene 중 5 개 (TSHR · NIS/SLC5A5 · TPO · TG · DIO1) 가 이 회로의 정확한 노드.  Panel = 회로 그 자체.",
    pmid: "28153798",
    doi: "10.1016/j.mce.2017.01.038"
  },
  {
    tier: "FOUNDATION 2",
    topic: "Lineage TFs · PAX8 + NKX2-1 + FOXE1",
    citation: "Fernández LP, López-Márquez A, Santisteban P.  (2015)",
    journal: "Nature Reviews Endocrinology  ·  11 (1) : 29-42",
    key_finding: "PAX8, NKX2-1 (TTF-1), FOXE1 이 갑상선 여포세포의 3 대 lineage TF 로서 TG · TPO · TSHR · NIS enhancer 를 공동 activation.  이 3 TF 어느 하나만 없어도 갑상선 유전자 발현 침묵.",
    fig_desc: "Figure 2 · TF binding site architecture at TG · TPO · TSHR · NIS promoters — cooperative activation model.",
    supports: "우리 panel 의 TF 3 개 (PAX8 · NKX2-1 · FOXE1) 를 선택한 근거.  DM1 에서 3 TF 동시 침묵이 effector silencing 을 설명.",
    pmid: "25350068",
    doi: "10.1038/nrendo.2014.186"
  },
  {
    tier: "FOUNDATION 3",
    topic: "NIS (SLC5A5) — the RAI symporter",
    citation: "De la Vieja A, Dohan O, Levy O, Carrasco N.  (2000)",
    journal: "Physiological Reviews  ·  80 (3) : 1083-1105",
    key_finding: "Na⁺/I⁻ symporter (NIS / SLC5A5) — basolateral membrane symporter 가 blood I⁻ 를 Na⁺ 와 함께 세포 안으로 능동수송.  이것이 ¹³¹I RAI 치료의 물리적 기반.",
    fig_desc: "Figure 4 · NIS topology, transport stoichiometry (2 Na⁺ : 1 I⁻), TSH/cAMP regulation.",
    supports: "NIS silencing → RAI 흡수 실패의 direct 메커니즘.  DM1 에서 NIS β 상승 (침묵) 이 관찰됨.",
    pmid: "10893432",
    doi: "10.1152/physrev.2000.80.3.1083"
  },
  {
    tier: "FOUNDATION 4",
    topic: "TCGA-THCA integrated genomics · BRS / TDS",
    citation: "Cancer Genome Atlas Research Network.  (2014)",
    journal: "Cell  ·  159 (3) : 676-690",
    key_finding: "TCGA-THCA n = 496 통합 유전체 분석.  BRAF-like vs RAS-like transcriptomic 축 도입.  Thyroid Differentiation Score (TDS) 정량화.  BRAF+ 에서 TDS 감소, RAS+ 에서 유지.",
    fig_desc: "Figure 4-5 · TDS distribution across drivers · BRS scoring · genomic-transcriptomic linkage.",
    supports: "우리 DM1 축이 TDS (Thyroid Differentiation Score) 와 동일 방향으로 정렬 — convergent biology 증명.  BRS-BRAF-like end 가 DM1 과 부분 overlap.",
    pmid: "25417114",
    doi: "10.1016/j.cell.2014.09.050"
  },
  {
    tier: "FOUNDATION 5",
    topic: "PDTC / ATC dedifferentiation · silenced gene list",
    citation: "Landa I, Ibrahimpasic T, Boucai L, et al.  (2016)",
    journal: "Journal of Clinical Investigation  ·  126 (3) : 1052-1066",
    key_finding: "Poorly differentiated + anaplastic thyroid carcinoma (PDTC/ATC) n = 149 유전체 분석.  분화 손실 시 silenced 되는 lineage effector list (TG · TPO · TSHR · NIS · DIO1 등).",
    fig_desc: "Figure 3 · progressive dedifferentiation signature — WDTC → PDTC → ATC 축을 따라 lineage gene 침묵.  Table 3 · silenced gene list.",
    supports: `우리 8-gene 중 <b>5 개가 Landa 2016 silenced list 와 직접 overlap</b> — 독립 설계에서 도출된 convergent biology.  DM1 = 「less-differentiated」 축의 재현.`,
    pmid: "26878173",
    doi: "10.1172/JCI85271"
  },
  {
    tier: "FOUNDATION 6",
    topic: "MAPK-driven NIS silencing · redifferentiation possible",
    citation: "Chakravarty D, Santos E, Ryder M, et al.  (2011)",
    journal: "Journal of Clinical Investigation  ·  121 (12) : 4700-4711",
    key_finding: "BRAF V600E conditional 마우스 모델에서 MAPK 활성이 NIS 를 direct 하게 silencing.  MEK inhibitor 로 MAPK 차단 시 NIS 재발현 + ¹³¹I 흡수 회복.",
    fig_desc: "Figure 4-6 · MEK inhibitor 처리 시 NIS mRNA/protein 회복 · in vivo ¹²⁴I PET uptake 증가.",
    supports: "DM1 침묵이 <em>reversible</em> 하다는 mechanistic 근거.  MAPK 축 → lineage silencing → RAI 실패의 causal chain.",
    pmid: "22105174",
    doi: "10.1172/JCI46382"
  },
  {
    tier: "FOUNDATION 7",
    topic: "Clinical proof · Selumetinib restores RAI uptake",
    citation: "Ho AL, Grewal RK, Leboeuf R, et al.  (2013)",
    journal: "New England Journal of Medicine  ·  368 : 623-632",
    key_finding: "RAI-refractory metastatic DTC 20 명 pilot trial.  Selumetinib (MEK inhibitor) 4 주 후 ¹²⁴I PET 재검 → 12/20 환자에서 iodine uptake 증가, 8/20 에서 ¹³¹I 재치료 자격 회복.  BRAF+ 보다 NRAS+ 반응 강함.",
    fig_desc: "Figure 1 · pre vs post treatment ¹²⁴I PET/CT.  Figure 2 · MEK inhibitor 반응 vs driver mutation.",
    supports: "임상에서 dedifferentiation-driven RAI failure 가 실제 존재하며 partially reversible 임을 human 에서 증명.  DM1 stratification 의 임상 utility 근거.",
    pmid: "23406027",
    doi: "10.1056/NEJMoa1209288"
  },
  {
    tier: "FOUNDATION 8",
    topic: "BRAF+ redifferentiation with dabrafenib",
    citation: "Rothenberg SM, McFadden DG, Palmer EL, et al.  (2015)",
    journal: "Clinical Cancer Research  ·  21 (5) : 1028-1035",
    key_finding: "BRAF V600E RAI-refractory PTC 10 명 pilot.  Dabrafenib (BRAF inhibitor) 6 주 후 ¹²⁴I whole-body scan → 6/10 iodine uptake 회복 + ¹³¹I 재치료 후 tumor shrinkage.",
    fig_desc: "Figure 1 · patient-level ¹²⁴I uptake pre/post.  Figure 2 · dosimetry recovery.  Case examples.",
    supports: "BRAF-DM1 subset 에서 targeted redifferentiation 이 임상 가능한 첫 human proof.  분당 프로토콜 templating 근거.",
    pmid: "25549723",
    doi: "10.1158/1078-0432.CCR-14-2915"
  },
  {
    tier: "FOUNDATION 9",
    topic: "Thyroid organogenesis · TF hierarchy",
    citation: "Nilsson M, Fagman H.  (2017)",
    journal: "Development  ·  144 (12) : 2123-2140",
    key_finding: "발생 단계에서 갑상선 정체성 확립 — PAX8, NKX2-1, FOXE1, HHEX 4 TF 가 여포세포 lineage commitment 를 결정.  TF 발현이 없으면 갑상선 자체가 형성되지 않음 (congenital hypothyroidism).",
    fig_desc: "Figure 3-4 · TF expression timeline in developing thyroid · knock-out phenotypes.",
    supports: "우리 3 TF (PAX8 · NKX2-1 · FOXE1) 가 lineage-defining 이며 종양에서 침묵 시 identity 손실 → dedifferentiation.",
    pmid: "28634182",
    doi: "10.1242/dev.145615"
  },
  {
    tier: "FOUNDATION 10",
    topic: "RAI-refractory advanced DTC clinical outcomes",
    citation: "Schlumberger M, Tahara M, Wirth LJ, et al.  (2015)",
    journal: "New England Journal of Medicine  ·  372 : 621-630",
    key_finding: "RAI-refractory advanced DTC 392 명 phase 3 (SELECT trial).  Lenvatinib vs placebo → median PFS 18.3 vs 3.6 개월 (HR = 0.21).  이 환자군의 unmet need + magnitude 확립.",
    fig_desc: "Figure 2 · PFS Kaplan-Meier curves.  Figure 3 · tumor shrinkage waterfall.",
    supports: "RAI-refractory 환자군이 임상적으로 명확히 정의 가능 + TKI 시대에도 여전히 심각한 unmet need — DM1 이 pre-identify 하면 sequencing 조기화 가능.",
    pmid: "25671254",
    doi: "10.1056/NEJMoa1406470"
  }
];

const linkify = (f: Foundation) => {
  const links: JSX.Element[] = [];
  if (f.pmid) {
    links.push(
      <a key="pmid" href={`https://pubmed.ncbi.nlm.nih.gov/${f.pmid}/`} target="_blank" rel="noreferrer" className="biofound-link biofound-link-pmid">
        PMID {f.pmid}
      </a>
    );
  }
  if (f.doi) {
    links.push(
      <a key="doi" href={`https://doi.org/${f.doi}`} target="_blank" rel="noreferrer" className="biofound-link biofound-link-doi">
        DOI {f.doi}
      </a>
    );
  }
  return links;
};

export const BiologyFoundations: React.FC = () => (
  <div className="biofound">
    <div className="biofound-intro">
      <h3>DM1 축의 8 개 gene 은 이 10 개 landmark 논문이 세운 갑상선 분화 생물학 위에 서 있다</h3>
      <p>
        4 개 categorization : (1) <b>세포 회로 · TF hierarchy</b> — FDN 1, 2, 3, 9 ·
        (2) <b>Convergent transcriptomic axis</b> (BRS/TDS · Landa silenced list) — FDN 4, 5 ·
        (3) <b>MAPK-driven silencing reversible</b> — FDN 6, 7, 8 ·
        (4) <b>Clinical burden of RAI failure</b> — FDN 10.
      </p>
      <p className="biofound-copyright-note">
        <b>저작권 안내:</b> 원 논문 figure 는 각 출판사 저작권 · 여기서는 각 논문의 핵심 concept 요약 + PubMed / DOI 링크만 제공.
        미팅 중 clicking 즉시 원문 figure 확인 가능.
      </p>
    </div>

    <div className="biofound-grid">
      {FOUNDATIONS.map((f, i) => (
        <div key={i} className="biofound-card">
          <div className="biofound-card-header">
            <span className="biofound-tier">{f.tier}</span>
            <span className="biofound-topic">{f.topic}</span>
          </div>
          <div className="biofound-cite">
            <b>{f.citation}</b>
            <div className="biofound-journal">{f.journal}</div>
          </div>
          <div className="biofound-row">
            <div className="biofound-label">Key finding</div>
            <div className="biofound-body">{f.key_finding}</div>
          </div>
          <div className="biofound-row">
            <div className="biofound-label">Original figure</div>
            <div className="biofound-body biofound-figdesc">{f.fig_desc}</div>
          </div>
          <div className="biofound-row supports">
            <div className="biofound-label supports">→ 우리 논문 근거</div>
            <div className="biofound-body" dangerouslySetInnerHTML={{__html: f.supports}} />
          </div>
          <div className="biofound-links">
            {linkify(f)}
          </div>
        </div>
      ))}
    </div>

    <div className="biofound-conclusion">
      <div className="biofound-conclusion-head">종합 — 왜 이 10 개 논문 위에 우리 결과가 서는가</div>
      <div className="biofound-conclusion-body">
        위 10 논문은 (a) 갑상선 여포세포 회로 (b) 그 회로를 조절하는 3 lineage TF (c) 종양에서 이 회로의 침묵 = dedifferentiation
        (d) 그 침묵이 RAI failure 임상 표현형과 직결됨 을 각각 확립했다.
        우리 <b>8-gene DM1 축</b> 은 이 4 층 위에 <em>정량 axis</em> 를 얹은 것 —
        기존 biology 를 부정하지 않고, <b>수술 시점 조직에서 미리 감지 가능한 분자 표현형</b> 으로 재구성한 것.
      </div>
    </div>
  </div>
);
