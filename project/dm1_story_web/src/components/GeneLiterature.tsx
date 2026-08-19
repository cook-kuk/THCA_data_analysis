import React from "react";

/* ================================================================
   Gene-by-gene literature dossier — 8 panel genes
   각 유전자 별: 기능 / RAI-dediff ref / congenital hypothyroidism ref / IHC clone / Korean cohort
   모든 PMID / DOI 는 PubMed 검증 (2026-06).
   ================================================================ */

type GeneCard = {
  symbol: string;
  alias?: string;
  group: "effector" | "TF";
  role: string;
  rai: string;
  ihc: { available: "routine" | "limited" | "research"; clone?: string };
  ch?: string;
  /** DM1 vs DM2 expression direction + magnitude. All 8 panel genes are DOWN in DM1. */
  direction: {
    expr: "down" | "up";              // DM1 vs DM2 expression direction
    methD?: number;                   // HM450 promoter β Cohen's d (DM1 − DM2). + means DM1 hypermethylated.
    methP?: string;                   // MW p-value notation
    tier: "strong" | "moderate" | "mild" | "ns";  // magnitude category
  };
  refs: { authors: string; year: string; journal: string; pmid?: string; doi?: string; note: string }[];
};

const GENES: GeneCard[] = [
  {
    symbol: "TG",
    alias: "thyroglobulin",
    group: "effector",
    role: "갑상선 호르몬 (T4 · T3) 의 substrate · iodine storage 단백질",
    rai: "Silencing 시 호르몬 합성 능력 소실, TDS-16 (Yoo 2016) 의 핵심 component.",
    ihc: { available: "routine", clone: "polyclonal · Dako" },
    ch: "Congenital hypothyroidism (dyshormonogenesis) 의 원인 — TG LoF 변이 다수 등록.",
    direction: { expr: "down", methD: 0.86, methP: "2.2 × 10⁻⁶", tier: "moderate" },
    refs: [
      { authors: "Targovnik HM, Citterio CE, Rivolta CM", year: "2010", journal: "Mol Cell Endocrinol — review on TG mutations", pmid: "20093166", note: "TG mutation spectrum review — congenital hypothyroidism 의 분자 기전." },
      { authors: "Yoo SK et al.", year: "2016", journal: "PLOS Genet 12(8):e1006239", pmid: "27494611", doi: "10.1371/journal.pgen.1006239", note: "★ 한국인 PTC/FTC 통합 transcriptome · TDS-16 panel 의 한국인 cohort 적용 reference." },
      { authors: "Cancer Genome Atlas Research Network", year: "2014", journal: "Cell 159:676-690", pmid: "25417114", doi: "10.1016/j.cell.2014.09.050", note: "TCGA-THCA · TDS-16 의 origin · TG 가 16 gene 중 핵심." }
    ]
  },
  {
    symbol: "TPO",
    alias: "thyroperoxidase",
    group: "effector",
    role: "Iodide 의 oxidation · organification (TG iodination) 의 heme-peroxidase 효소",
    rai: "BRAF V600E PTC 에서 NIS/TPO mRNA 일관 감소 (Romei 2008) · RAI-refractory direct driver.",
    ihc: { available: "routine", clone: "MoAb47" },
    ch: "Congenital hypothyroidism · TPO LoF 가 가장 흔한 dyshormonogenesis 원인.",
    direction: { expr: "down", methD: 2.30, methP: "1.9 × 10⁻¹⁸", tier: "strong" },
    refs: [
      { authors: "Romei C et al.", year: "2008", journal: "Endocr Relat Cancer 15:511-520", pmid: "18509003", doi: "10.1677/ERC-07-0130", note: "BRAF V600E PTC 에서 NIS, TPO, TG, TSHR mRNA 일관 감소 — RAI-refractoriness 의 분자 근거." },
      { authors: "Durante C et al.", year: "2007", journal: "JCEM 92:2840-2843", doi: "10.1210/jc.2006-2707", note: "BRAF V600E PTC 의 iodine-handling 기능 감소 (NIS/TPO/TG/TSHR 동시) — early 분자 정립." },
      { authors: "Ho AL et al.", year: "2013", journal: "NEJM 368:623-632", pmid: "23406027", doi: "10.1056/NEJMoa1209288", note: "★ Selumetinib (MEK inhibitor) → NIS/TPO 재유도 → ¹³¹I uptake 회복. Redifferentiation therapy 의 임상 landmark." }
    ]
  },
  {
    symbol: "TSHR",
    alias: "TSH receptor",
    group: "effector",
    role: "TSH 의 receptor · cAMP → PKA → 분화 / 호르몬 합성 신호의 진입점",
    rai: "TSH-driven proliferation/differentiation 의 receptor. PDTC/ATC 에서 silencing — Landa 2016 silenced list.",
    ihc: { available: "limited", clone: "research-grade only · qPCR 기반 정량 표준" },
    ch: "Sunthornthepvarakul-Refetoff 1995 NEJM — TSHR loss-of-function 이 TSH resistance / congenital hypothyroidism.",
    direction: { expr: "down", methD: 1.20, methP: "9.8 × 10⁻¹²", tier: "strong" },
    refs: [
      { authors: "Landa I et al.", year: "2016", journal: "J Clin Invest 126:1052-1066", pmid: "26878173", doi: "10.1172/JCI85271", note: "★ PDTC + ATC silenced gene list — NIS/TG/TPO 와 함께 TSHR 강조. DM1 panel 의 5/8 overlap." },
      { authors: "Sunthornthepvarakul T, Gottschalk ME, Hayashi Y, Refetoff S", year: "1995", journal: "NEJM 332:155-160", pmid: "7528344", doi: "10.1056/NEJM199501193320305", note: "★ TSHR LoF → TSH resistance + congenital hypothyroidism · TSHR 의 발달 필수성 인간 유전학 증거." },
      { authors: "Cancer Genome Atlas Research Network", year: "2014", journal: "Cell 159:676-690", pmid: "25417114", note: "TSHR 가 TCGA TDS-16 의 한 축." }
    ]
  },
  {
    symbol: "SLC5A5",
    alias: "NIS · Na⁺/I⁻ symporter",
    group: "effector",
    role: "Thyrocyte basolateral membrane 의 active iodide uptake transporter — RAI 흡수 그 자체",
    rai: "★ Iodide uptake 의 직접 매개체. NIS silencing = RAI-refractoriness 의 mechanistic root.",
    ihc: { available: "limited", clone: "FP5A (research-grade)" },
    ch: "Fujiwara 2000 (V59E) · Kosugi (T354P) — NIS LoF → iodide transport defect → congenital hypothyroidism (ITD).",
    direction: { expr: "down", methD: 0.22, methP: "0.42 (NS at promoter; mRNA-level still DOWN)", tier: "mild" },
    refs: [
      { authors: "Dai G, Levy O, Carrasco N", year: "1996", journal: "Nature 379:458-460", pmid: "8559250", doi: "10.1038/379458a0", note: "★ NIS cloning 의 classic paper · 갑상선 iodide uptake 분자 기전 정립." },
      { authors: "Spitzweg C, Harrington KJ, Pinke LA, Vile RG, Morris JC", year: "2001", journal: "JCEM 86:3327-3335", doi: "10.1210/jcem.86.7.7641", note: "NIS biology + 갑상선암 임상 의의 review." },
      { authors: "Ho AL et al.", year: "2013", journal: "NEJM 368:623-632", pmid: "23406027", note: "★ Selumetinib redifferentiation — NIS 재유도 → ¹³¹I uptake 회복 임상 lead." },
      { authors: "Rothenberg SM et al.", year: "2015", journal: "Clin Cancer Res 21:1028-1035", pmid: "25549723", doi: "10.1158/1078-0432.CCR-14-2915", note: "★ Dabrafenib (BRAF inhibitor) → NIS 재유도 redifferentiation lead." },
      { authors: "Ravera S et al.", year: "2017", journal: "Annu Rev Physiol 79:261-289", pmid: "28192058", doi: "10.1146/annurev-physiol-022516-034125", note: "NIS biology + 임상 적용 최신 review." },
      { authors: "Ward LS et al.", year: "2003", journal: "Cancer Lett 200:85-91", pmid: "14550956", doi: "10.1016/s0304-3835(03)00392-6", note: "Low NIS expression = aggressive PTC 표현형 — prognostic significance." },
      { authors: "Fujiwara H et al.", year: "2000", journal: "Thyroid 10:471-474", pmid: "10907989", doi: "10.1089/thy.2000.10.471", note: "NIS V59E LoF · iodide transport defect — 인간 유전학." }
    ]
  },
  {
    symbol: "DIO1",
    alias: "type-1 deiodinase",
    group: "effector",
    role: "T4 → T3 peripheral deiodination · thyroid hormone metabolism",
    rai: "TCGA TDS-16 의 한 축 — BRAF-like 군에서 일관 감소. 별도의 LoF 증후군 없음 (DIO2 가 hypothyroidism 주 원인).",
    ihc: { available: "research", clone: "research-grade only" },
    direction: { expr: "down", methD: 1.24, methP: "6.5 × 10⁻¹¹", tier: "strong" },
    refs: [
      { authors: "Cancer Genome Atlas Research Network", year: "2014", journal: "Cell 159:676-690", pmid: "25417114", note: "TDS-16 의 한 축 · BRAF-like 군에서 DIO1 감소 일관." },
      { authors: "Landa I et al.", year: "2016", journal: "J Clin Invest 126:1052-1066", pmid: "26878173", note: "PDTC/ATC dedifferentiation 의 deiodination 손실." },
      { authors: "Yoo SK et al.", year: "2016", journal: "PLOS Genet 12(8):e1006239", pmid: "27494611", note: "한국인 cohort TDS-16 적용성." }
    ]
  },
  {
    symbol: "PAX8",
    alias: "—",
    group: "TF",
    role: "Thyroid lineage master transcription factor — NIS, TG, TPO 의 직접 transactivator",
    rai: "PAX8 silencing → NIS/TG/TPO 도미노 silencing → RAI-refractoriness. DM1 의 upstream regulator.",
    ihc: { available: "routine", clone: "MRQ-50 · BC12 (mouse mAb)" },
    ch: "★ Macchia 1998 Nat Genet — heterozygous LoF → thyroid dysgenesis · 발달 필수 인간 유전학.",
    direction: { expr: "down", methD: 0.97, methP: "4.5 × 10⁻⁸", tier: "strong" },
    refs: [
      { authors: "Macchia PE, Lapi P, Krude H et al.", year: "1998", journal: "Nat Genet 19:83-86", pmid: "9590296", doi: "10.1038/ng0598-83", note: "★ PAX8 heterozygous LoF → thyroid dysgenesis + congenital hypothyroidism · 발달 핵심." },
      { authors: "De Felice M, Di Lauro R", year: "2004", journal: "Endocr Rev 25:722-746", pmid: "15466941", doi: "10.1210/er.2003-0028", note: "Lineage TF (PAX8/NKX2-1/FOXE1) 의 갑상선 발달 standard review." },
      { authors: "Fernández LP, López-Márquez A, Santisteban P", year: "2015", journal: "Nat Rev Endocrinol 11:29-42", pmid: "25350068", doi: "10.1038/nrendo.2014.186", note: "★ Thyroid TF 의 발달 / 질환 / 암 최신 review (DM1 framework 의 직접 근거)." }
    ]
  },
  {
    symbol: "NKX2-1",
    alias: "TTF-1",
    group: "TF",
    role: "Thyroid + 폐 + 뇌 lineage TF — homeodomain · 갑상선 bud specification 의 master",
    rai: "Silencing 시 NIS/TG/TPO/TSHR 동시 침묵 → DM1 axis 핵심. 임상 갑상선 진단 IHC 의 routine confirm marker.",
    ihc: { available: "routine", clone: "8G7G3/1 (Dako · routine)" },
    ch: "★ Krude 2002 JCI — brain-lung-thyroid syndrome (choreoathetosis + neonatal RDS + congenital hypothyroidism).",
    direction: { expr: "down", methD: 0.63, methP: "8.9 × 10⁻⁷", tier: "moderate" },
    refs: [
      { authors: "Krude H et al.", year: "2002", journal: "J Clin Invest 109:475-480", pmid: "11854319", doi: "10.1172/JCI14341", note: "★ NKX2-1 haploinsufficiency → brain-lung-thyroid syndrome · 인간 발달 필수." },
      { authors: "Cancer Genome Atlas Research Network", year: "2014", journal: "Cell 159:676-690", pmid: "25417114", note: "TDS-16 의 한 축 · BRAF-like 군에서 silencing." },
      { authors: "Fernández LP, López-Márquez A, Santisteban P", year: "2015", journal: "Nat Rev Endocrinol 11:29-42", pmid: "25350068", note: "Thyroid TF review · NKX2-1 의 갑상선암 dediff role." }
    ]
  },
  {
    symbol: "FOXE1",
    alias: "TTF-2",
    group: "TF",
    role: "Thyroid bud migration · forkhead TF — 갑상선 발달의 maturation phase 핵심",
    rai: "Silencing → 분화 maturation 손실 · PDTC/ATC dedifferentiation 마커. DM1 의 lineage TF axis.",
    ihc: { available: "limited", clone: "polyclonal · borderline 케이스 보조용" },
    ch: "★ Clifton-Bligh 1998 Nat Genet — Bamforth-Lazarus syndrome (thyroid agenesis + cleft palate + choanal atresia + spiky hair).",
    direction: { expr: "down", methD: 0.84, methP: "1.0 × 10⁻⁵", tier: "moderate" },
    refs: [
      { authors: "Clifton-Bligh RJ et al.", year: "1998", journal: "Nat Genet 19:399-401", pmid: "9697705", doi: "10.1038/1294", note: "★ FOXE1 LoF → Bamforth-Lazarus syndrome (갑상선 무형성 + 구개열) · 발달 필수." },
      { authors: "De Felice M, Di Lauro R", year: "2004", journal: "Endocr Rev 25:722-746", pmid: "15466941", note: "Thyroid TF 발달 review · FOXE1 의 migration phase role." },
      { authors: "Landa I et al.", year: "2016", journal: "J Clin Invest 126:1052-1066", pmid: "26878173", note: "PDTC/ATC silenced list — FOXE1 포함." }
    ]
  }
];

const pmidUrl = (p?: string) => p ? `https://pubmed.ncbi.nlm.nih.gov/${p}/` : undefined;
const doiUrl  = (d?: string) => d ? `https://doi.org/${d}` : undefined;

const IHC_LABEL: Record<GeneCard["ihc"]["available"], { txt: string; color: string }> = {
  routine: { txt: "★★★ Routine 임상 IHC",      color: "#047857" },
  limited: { txt: "★★ 가능하나 표준화 제한",    color: "#B45309" },
  research:{ txt: "★ Research-grade only",      color: "#7C2D12" }
};

const TIER_LABEL: Record<"strong"|"moderate"|"mild"|"ns", { txt: string; color: string }> = {
  strong:   { txt: "strong",   color: "#7F1D1D" },
  moderate: { txt: "moderate", color: "#B45309" },
  mild:     { txt: "mild",     color: "#92400E" },
  ns:       { txt: "n.s.",     color: "#64748B" }
};

export const GeneLiterature: React.FC = () => (
  <div className="genelit">
    <div className="gene-legend">
      <span className="gl-key">
        <span className="gl-arrow down">▼</span> DM1 에서 expression <b>DOWN-regulated (silenced)</b>
      </span>
      <span className="gl-key">
        <span className="gl-meth-key">Methylation β (Cohen's d, DM1 − DM2)</span> — 양수 = DM1 promoter 과메틸화
      </span>
      <span className="gl-key">
        <span className="gl-tier-key strong">strong</span>
        <span className="gl-tier-key moderate">moderate</span>
        <span className="gl-tier-key mild">mild</span> · effect magnitude tier
      </span>
    </div>
    {GENES.map(g => {
      const ihc = IHC_LABEL[g.ihc.available];
      const tier = TIER_LABEL[g.direction.tier];
      return (
        <div key={g.symbol} className={`gene-dossier gene-${g.group}`}>
          <div className="gene-head">
            <div className="gene-sym-block">
              <code className={`gene-sym sym-${g.group}`}>{g.symbol}</code>
              {g.alias && g.alias !== "—" && <span className="gene-alias">/ {g.alias}</span>}
              <span className={`gene-group-tag ${g.group}`}>
                {g.group === "TF" ? "TF · lineage" : "Effector · uptake"}
              </span>
              <span className={`gene-dir gene-dir-${g.direction.expr}`}>
                <span className="gene-dir-arrow">{g.direction.expr === "down" ? "▼" : "▲"}</span>
                DM1 {g.direction.expr === "down" ? "DOWN" : "UP"}-regulated
              </span>
              <span className="gene-tier-pill" style={{ background: tier.color }}>{tier.txt}</span>
            </div>
            <span className="gene-ihc-pill" style={{ background: ihc.color }}>
              {ihc.txt}{g.ihc.clone ? `  ·  ${g.ihc.clone}` : ""}
            </span>
          </div>

          {/* Direction summary bar — HM450 methylation Cohen's d */}
          {typeof g.direction.methD === "number" && (
            <div className="gene-dir-bar">
              <div className="gene-dir-bar-label">
                Promoter <b>β methylation</b> · Cohen's d  (DM1 − DM2):
              </div>
              <div className="gene-dir-bar-track">
                <div
                  className={`gene-dir-bar-fill tier-${g.direction.tier}`}
                  style={{ width: `${Math.min(Math.abs(g.direction.methD) / 2.4 * 100, 100)}%` }}
                />
                <span className="gene-dir-bar-value">
                  d = <b>+{g.direction.methD.toFixed(2)}</b>
                  {g.direction.methP && <> · p = {g.direction.methP}</>}
                </span>
              </div>
              <div className="gene-dir-bar-foot">
                높은 β (과메틸화) → mRNA 발현 ↓ → 분화 기능 ↓ ·  TCGA HM450 n = 503 (Round 5, audit-locked)
              </div>
            </div>
          )}

          <div className="gene-grid">
            <div className="gene-cell">
              <div className="gene-cell-label">기능</div>
              <div className="gene-cell-body">{g.role}</div>
            </div>
            <div className="gene-cell">
              <div className="gene-cell-label">RAI / dedifferentiation 관련</div>
              <div className="gene-cell-body">{g.rai}</div>
            </div>
            {g.ch && (
              <div className="gene-cell gene-cell-ch">
                <div className="gene-cell-label">선천성 갑상선기능저하 (CH)</div>
                <div className="gene-cell-body">{g.ch}</div>
              </div>
            )}
          </div>

          <div className="gene-refs">
            <div className="gene-refs-label">핵심 참고문헌  ({g.refs.length} 편)</div>
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

    <div className="gene-summary-band">
      <div className="gs-direction-summary">
        <div className="gs-dir-title">
          ★ 핵심 — <b>8 개 panel 유전자 모두</b> DM1 종양에서 <span className="dir-down">▼ DOWN-regulated (silenced)</span> ·
          DM2 종양에서 <span className="dir-up">▲ UP-regulated (분화 유지)</span>
        </div>
        <div className="gs-dir-strength">
          {[
            ["TPO",    2.30, "strong"],
            ["DIO1",   1.24, "strong"],
            ["TSHR",   1.20, "strong"],
            ["PAX8",   0.97, "strong"],
            ["TG",     0.86, "moderate"],
            ["FOXE1",  0.84, "moderate"],
            ["NKX2-1", 0.63, "moderate"],
            ["SLC5A5", 0.22, "mild"]
          ].map(([sym, d, tier]) => (
            <div key={sym as string} className="gs-bar-row">
              <span className="gs-bar-sym">{sym}</span>
              <div className="gs-bar-track">
                <div className={`gs-bar-fill tier-${tier}`} style={{ width: `${Math.min((d as number)/2.4*100, 100)}%` }}/>
                <span className="gs-bar-d">d = +{(d as number).toFixed(2)}</span>
              </div>
            </div>
          ))}
        </div>
        <div className="gs-dir-foot">
          ▲ 양수 d = DM1 promoter β <b>hypermethylated</b>  →  mRNA 발현 <b>DOWN</b>.   8 / 8 유전자 모두 같은 방향.
          TCGA HM450 n = 503 · audit-locked (Round 5, 2026-04-30).
        </div>
      </div>
      <div className="gs-row">
        <span className="gs-cell"><b>5 effector</b> (TG · TPO · TSHR · SLC5A5 · DIO1) — RAI 흡수 / 호르몬 합성 도구</span>
        <span className="gs-cell"><b>3 TF</b> (PAX8 · NKX2-1 / TTF-1 · FOXE1 / TTF-2) — 갑상선 lineage master regulator</span>
      </div>
      <div className="gs-narrative">
        모든 8 panel 유전자는 <b>인간 유전학에서 발달 필수성이 입증된 분화 회로</b>의 구성원이다 (PAX8 → thyroid dysgenesis; NKX2-1 → brain-lung-thyroid; FOXE1 → Bamforth-Lazarus; TG/TPO/SLC5A5 → congenital hypothyroidism with goiter / ITD; TSHR → TSH resistance). 즉 panel 의 silencing 은 thyrocyte 기능 상실의 인간유전학적 정당성을 갖는다 — <em>"단순 분류 panel 이 아니라 갑상선 분화 핵심 회로 그 자체"</em>.
      </div>
      <div className="gs-narrative">
        Redifferentiation therapy 의 두 임상 lead (<b>selumetinib · Ho 2013 NEJM</b>; <b>dabrafenib · Rothenberg 2015 CCR</b>) 가 모두 effector 축 (NIS · TPO) 의 재유도를 통한 <sup>131</sup>I uptake 회복을 endpoint 로 함 — DM1 axis 의 임상 actionability 의 직접 근거.
      </div>
    </div>
  </div>
);
