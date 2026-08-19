import React from "react";

/**
 * Two-layer diagram:
 *   Layer 1 (top): driver categories — BRAF V600E / RAS / fusion / driver-neg
 *   Layer 2 (bottom): DM1 / DM2 differentiation state, cutting across drivers
 * Diagonal lines/links show that DM1 status is not driver-determined.
 */
export const DriverLayerDiagram: React.FC = () => {
  const drivers = [
    { x: 90,  label: "BRAF V600E", dm1pct: 0.7 },
    { x: 290, label: "RAS-mutant", dm1pct: 96.4 },
    { x: 490, label: "Fusion+", dm1pct: 60 },
    { x: 690, label: "Driver-neg", dm1pct: 49.1 }
  ];
  return (
    <svg viewBox="0 0 900 380" className="diagram" preserveAspectRatio="xMidYMid meet">
      <defs>
        <marker id="dlArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
          <path d="M0,0 L10,5 L0,10 z" fill="#475569" />
        </marker>
      </defs>

      <text x="450" y="32" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14.5" fontWeight="700" fill="#B45309">
        Driver 변이 ≠ DM 분화 상태  ·  DM1 은 driver 범주를 가로지른다
      </text>

      {/* layer 1 title */}
      <text x="40" y="80" fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#475569">Layer 1 · driver 변이</text>
      {/* layer 2 title */}
      <text x="40" y="280" fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#475569">Layer 2 · 분화 상태 (DM1 / DM2)</text>

      {/* driver boxes */}
      {drivers.map(d => (
        <g key={d.label}>
          <rect x={d.x} y={90} width={150} height={56} rx={10} fill="#F1F5F9" stroke="#475569" strokeWidth={1.2} />
          <text x={d.x + 75} y={114} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#0F172A">{d.label}</text>
          <text x={d.x + 75} y={134} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#64748B">DM1 prevalence  {d.dm1pct} %</text>
        </g>
      ))}

      {/* DM2 box */}
      <rect x={150} y={290} width={240} height={56} rx={10} fill="#EFF6FF" stroke="#1E40AF" strokeWidth={1.5} />
      <text x={270} y={314} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#1E40AF">DM2 (분화 유지)</text>
      <text x={270} y={332} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#1E3A8A">전형적 well-differentiated PTC</text>

      {/* DM1 box */}
      <rect x={510} y={290} width={240} height={56} rx={10} fill="#FEF2F2" stroke="#B91C1C" strokeWidth={1.5} />
      <text x={630} y={314} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#B91C1C">DM1 (분화 손실)</text>
      <text x={630} y={332} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#7F1D1D">RAI machinery off · aggressive</text>

      {/* arrows showing driver → state crossing */}
      {drivers.map(d => {
        const sx = d.x + 75;
        const sy = 148;
        // proportional split: most BRAF → DM2, most RAS → DM1, fusion ~ split, neg ~ split
        const toDM1Frac = d.dm1pct / 100;
        return (
          <g key={"link-" + d.label}>
            {/* line to DM2 */}
            <line x1={sx} y1={sy} x2={270} y2={290} stroke="#1E40AF" strokeWidth={1 + 4*(1 - toDM1Frac)} opacity={0.55} />
            {/* line to DM1 */}
            <line x1={sx} y1={sy} x2={630} y2={290} stroke="#B91C1C" strokeWidth={1 + 4*toDM1Frac} opacity={0.55} />
          </g>
        );
      })}

      <text x="450" y="370" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fontStyle="italic" fill="#475569">
        Driver mutation 은 출발점.  DM1 은 disease state.  하나의 driver 가 DM1 / DM2 를 결정하지 않는다.
      </text>
    </svg>
  );
};
