import React from "react";

/**
 * DM2 (differentiated, RAI machinery on)  ←————→  DM1 (dedifferentiated, RAI machinery off)
 * 자체 SVG; 가로축 + 두 끝 라벨 + 분화 표시기들.
 */
export const AxisDiagram: React.FC = () => (
  <svg viewBox="0 0 900 360" className="diagram" preserveAspectRatio="xMidYMid meet">
    <defs>
      <linearGradient id="axisGradient" x1="0" x2="1">
        <stop offset="0%" stopColor="#1E40AF" />
        <stop offset="50%" stopColor="#94A3B8" />
        <stop offset="100%" stopColor="#B91C1C" />
      </linearGradient>
      <marker id="arrowR" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#0F172A" />
      </marker>
    </defs>

    {/* axis bar */}
    <rect x="80" y="170" width="740" height="22" rx="11" fill="url(#axisGradient)" opacity="0.94" />

    {/* end labels */}
    <g>
      <rect x="40" y="100" width="200" height="58" rx="10" fill="#EFF6FF" stroke="#1E40AF" strokeWidth="1.5" />
      <text x="140" y="124" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="17" fontWeight="700" fill="#1E40AF">DM2</text>
      <text x="140" y="146" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#1E3A8A">분화된 상태 · RAI machinery on</text>

      <rect x="660" y="100" width="200" height="58" rx="10" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="1.5" />
      <text x="760" y="124" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="17" fontWeight="700" fill="#B91C1C">DM1</text>
      <text x="760" y="146" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#7F1D1D">탈분화 상태 · RAI machinery off</text>
    </g>

    {/* tickmarks for genes */}
    {[
      { x: 165, label: "TG", level: "high" },
      { x: 230, label: "TPO", level: "high" },
      { x: 305, label: "TSHR", level: "high" },
      { x: 380, label: "PAX8", level: "high" },
      { x: 455, label: "NKX2-1", level: "mid" },
      { x: 530, label: "FOXE1", level: "mid" },
      { x: 605, label: "SLC5A5", level: "low" },
      { x: 680, label: "DIO1", level: "low" }
    ].map(g => (
      <g key={g.label}>
        <line x1={g.x} y1="200" x2={g.x} y2="232" stroke="#0F172A" strokeWidth="1.2" />
        <text x={g.x} y="252" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#0F172A">{g.label}</text>
      </g>
    ))}

    {/* directional arrow */}
    <g>
      <line x1="80" y1="290" x2="820" y2="290" stroke="#0F172A" strokeWidth="1.5" markerEnd="url(#arrowR)" />
      <text x="450" y="313" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fontWeight="600" fill="#0F172A">
        분화 손실 · RAI 흡수 능력 감소
      </text>
    </g>

    {/* tag */}
    <text x="450" y="56" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14.5" fontWeight="700" fill="#B45309">
      Thyroid follicular differentiation axis (DM1 / DM2)
    </text>
  </svg>
);
