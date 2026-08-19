import React from "react";

/**
 * "Why these 8 genes?" — 7-step selection flow + reviewer defense.
 */
export const SelectionFlowDiagram: React.FC = () => {
  const steps = [
    { id: 1, label: "문헌 기반 후보군",       sub: "RAI uptake · iodine metabolism · thyroid differentiation 관련 유전자 수집",   color: "#94A3B8" },
    { id: 2, label: "기능별 분류",             sub: "Uptake · organification · TSH signaling · lineage TF · metabolism 5 범주",  color: "#0E7490" },
    { id: 3, label: "TIERA67 67-gene 후보 풀", sub: "Yoo 2016 PLOS Genet RAI biology prior 기반",                                  color: "#0E7490" },
    { id: 4, label: "Cross-cohort 재현성",     sub: "TCGA + Lee + K2 + Landa 등에서 일관 분리 여부 확인",                          color: "#B45309" },
    { id: 5, label: "Redundancy / collinearity", sub: "유전자 중복성 · 상관성 제거",                                              color: "#B45309" },
    { id: 6, label: "Parsimony / model 축소",  sub: "Compact panel — 임상 deploy 가능한 최소 set",                                 color: "#B91C1C" },
    { id: 7, label: "최종 8-gene panel",       sub: "TF₃ ⊕ effector₅ = TG · TPO · TSHR · SLC5A5 · DIO1 · PAX8 · NKX2-1 · FOXE1",  color: "#B91C1C" }
  ];

  return (
    <svg viewBox="0 0 1200 760" className="diagram" preserveAspectRatio="xMidYMid meet">
      <defs>
        <marker id="sfArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto">
          <path d="M0,0 L10,5 L0,10 z" fill="#475569" />
        </marker>
      </defs>

      <text x="600" y="34" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="16" fontWeight="700" fill="#B45309">
        Figure B · 8 개 gene 은 갑자기 나온 것이 아니다 — biology-guided 7-step distillation
      </text>

      {/* 7 vertical step boxes */}
      {steps.map((s, i) => {
        const y = 72 + i * 76;
        return (
          <g key={s.id}>
            {/* number circle */}
            <circle cx="90" cy={y + 30} r="22" fill={s.color}/>
            <text x="90" y={y + 36} textAnchor="middle" fontFamily="Inter" fontSize="18" fontWeight="700" fill="#FFFFFF">{s.id}</text>
            {/* box */}
            <rect x="124" y={y} width="640" height="60" rx="10" fill="#FFFFFF" stroke={s.color} strokeWidth="1.6"/>
            <text x="144" y={y + 26} fontFamily="Noto Sans KR" fontSize="15" fontWeight="700" fill="#0F172A">{s.label}</text>
            <text x="144" y={y + 46} fontFamily="Noto Sans KR" fontSize="12" fill="#475569">{s.sub}</text>
            {/* arrow */}
            {i < steps.length - 1 && (
              <line x1="90" y1={y + 52} x2="90" y2={y + 76} stroke="#475569" strokeWidth="1.6" markerEnd="url(#sfArrow)"/>
            )}
          </g>
        );
      })}

      {/* Reviewer attack defense box on the right */}
      <g>
        <rect x="800" y="80" width="380" height="280" rx="14" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="1.8"/>
        <text x="990" y="110" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#B91C1C">⚠ Reviewer 가 가장 먼저 공격할 부분</text>
        <text x="820" y="142" fontFamily="Noto Sans KR" fontSize="13" fontStyle="italic" fill="#7F1D1D">"왜 하필 이 8 개 gene 인가?"</text>

        <line x1="820" y1="156" x2="1160" y2="156" stroke="#FCA5A5" strokeWidth="1"/>

        <text x="820" y="180" fontFamily="Noto Sans KR" fontSize="12.5" fontWeight="700" fill="#0F172A">답변 — 4 가지 차단선</text>
        {[
          ["①", "Biology-prioritized — RAI 흡수 / 분화 생물학 기반"],
          ["②", "Literature-guided — Yoo 2016 prior (outcome 무관)"],
          ["③", "Cohort-validated — 19 외부 코호트 direction-consistent"],
          ["④", "Parsimony-driven — 16-gene 대비 동등 + 임상 deploy 가능"]
        ].map(([n, txt], i) => (
          <g key={i}>
            <text x="828" y={205 + i*28} fontFamily="Noto Sans KR" fontSize="13" fontWeight="700" fill="#B91C1C">{n}</text>
            <text x="852" y={205 + i*28} fontFamily="Noto Sans KR" fontSize="12" fill="#1F2937">{txt}</text>
          </g>
        ))}
      </g>

      {/* Comparison panel on the right bottom */}
      <g>
        <rect x="800" y="380" width="380" height="200" rx="14" fill="#FFFBEB" stroke="#B45309" strokeWidth="1.6"/>
        <text x="990" y="408" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#7C2D12">비교 검증 (manuscript 필수)</text>
        {[
          ["16-gene Yoo 2016 TDS", "본 panel 와 동등 (ΔAUC NS)"],
          ["Random 8-gene", "현저히 낮은 cluster recovery"],
          ["4-gene subset", "임상 NGS panel 호환 대비"],
          ["Pan-genome top-5000", "ARI 0.92 (panel artifact 아님)"]
        ].map(([k, v], i) => (
          <g key={i}>
            <text x="820" y={436 + i*30} fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#7C2D12">{k}</text>
            <text x="820" y={452 + i*30} fontFamily="Noto Sans KR" fontSize="11.5" fill="#92400E">→ {v}</text>
          </g>
        ))}
      </g>

      {/* Final message */}
      <text x="600" y="660" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#0F172A">
        결론: ML 이 아무 gene 이나 고른 것이 아니다.
      </text>
      <text x="600" y="684" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#475569" fontStyle="italic">
        thyroid biology 에 기반한 후보군에서 출발 → 공개 데이터 검증 → compact 8-gene panel 로 정리
      </text>
      <text x="600" y="710" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#B45309" fontWeight="600">
        8 genes = lens · DM1/DM2 = discovery
      </text>
    </svg>
  );
};
