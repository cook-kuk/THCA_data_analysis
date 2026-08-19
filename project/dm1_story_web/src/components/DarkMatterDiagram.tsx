import React from "react";

/**
 * Canonical drivers → partially explain disease → unexplained clinical behavior → DM1 axis
 */
export const DarkMatterDiagram: React.FC = () => (
  <svg viewBox="0 0 900 360" className="diagram" preserveAspectRatio="xMidYMid meet">
    <defs>
      <marker id="dmArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#475569" />
      </marker>
    </defs>

    {/* Left: canonical drivers */}
    <g>
      <rect x="40" y="50" width="220" height="200" rx="14" fill="#EFF6FF" stroke="#1E40AF" strokeWidth="1.5" />
      <text x="150" y="80" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#1E40AF">기존 driver 분류</text>
      <text x="150" y="120" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fill="#1E3A8A">BRAF V600E</text>
      <text x="150" y="148" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fill="#1E3A8A">RAS-mutant</text>
      <text x="150" y="176" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fill="#1E3A8A">RET / NTRK / ALK fusion</text>
      <text x="150" y="220" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#475569" fontStyle="italic">대부분의 PTC 를 설명</text>
    </g>

    {/* Middle: dark matter zone */}
    <g>
      <rect x="340" y="50" width="220" height="200" rx="14" fill="#1E293B" stroke="#0F172A" strokeWidth="1.5" />
      <text x="450" y="80" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#FBBF24">설명되지 않는 영역</text>
      <text x="450" y="115" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fill="#E5E7EB">BRAF / RAS / fusion</text>
      <text x="450" y="135" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fill="#E5E7EB">으로 분류되지 않거나</text>
      <text x="450" y="155" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fill="#E5E7EB">같은 군 안에서도</text>
      <text x="450" y="175" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fill="#E5E7EB">임상 경과가 다른 환자</text>
      <text x="450" y="220" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#FCD34D" fontStyle="italic">molecular dark matter</text>

      {/* dotted star markers */}
      {[[380,200],[420,205],[460,195],[500,210],[520,200]].map((p,i) => (
        <circle key={i} cx={p[0]} cy={p[1]} r="2" fill="#FBBF24" />
      ))}
    </g>

    {/* Right: DM1 axis */}
    <g>
      <rect x="640" y="50" width="220" height="200" rx="14" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="1.5" />
      <text x="750" y="80" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#B91C1C">DM1 분화 축</text>
      <text x="750" y="118" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#7F1D1D">숨겨진 분화 상태</text>
      <text x="750" y="138" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#7F1D1D">driver-orthogonal</text>
      <text x="750" y="158" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#7F1D1D">lineage state</text>
      <text x="750" y="178" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#7F1D1D">RAI 흡수 능력 감소</text>
      <text x="750" y="220" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#B45309" fontStyle="italic">8-gene readout 으로 측정</text>
    </g>

    {/* arrows */}
    <line x1="262" y1="150" x2="338" y2="150" stroke="#475569" strokeWidth="1.6" markerEnd="url(#dmArrow)" />
    <text x="300" y="140" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#475569">설명 안 됨</text>

    <line x1="562" y1="150" x2="638" y2="150" stroke="#475569" strokeWidth="1.6" markerEnd="url(#dmArrow)" />
    <text x="600" y="140" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#475569">새로운 축</text>

    {/* bottom caption */}
    <text x="450" y="295" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="600" fill="#0F172A">
      DM1 = canonical driver 분류 너머의 숨겨진 분화 축 (molecular dark matter)
    </text>
    <text x="450" y="320" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#64748B" fontStyle="italic">
      driver mutation 으로 설명되지 않는 임상 차이를 읽는 lens
    </text>
  </svg>
);
