import React from "react";

/**
 * Surgery → FFPE tissue → 8-gene readout → DM call → risk stratification → RAI strategy → prospective validation
 */
export const WorkflowDiagram: React.FC = () => {
  const steps = [
    { label: "1.  수술 조직", sub: "Primary PTC tissue", color: "#1E40AF" },
    { label: "2.  FFPE 표본", sub: "FFPE / FF — KS p = 0.44", color: "#1E40AF" },
    { label: "3.  8-gene readout", sub: "RNA → 8-유전자 점수", color: "#0E7490" },
    { label: "4.  DM1 / DM2 호출", sub: "분화 축 분류", color: "#B45309" },
    { label: "5.  위험 분층화", sub: "ATA 중간 위험군에서 의미", color: "#B45309" },
    { label: "6.  RAI 전략 보조", sub: "치료 결정의 분자 입력", color: "#B91C1C" },
    { label: "7.  전향 검증", sub: "prospective validation required", color: "#047857" }
  ];
  const w = 1200, h = 220;
  const stepW = (w - 80) / steps.length;
  return (
    <svg viewBox={`0 0 ${w} ${h}`} className="diagram" preserveAspectRatio="xMidYMid meet">
      <defs>
        <marker id="wfArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
          <path d="M0,0 L10,5 L0,10 z" fill="#94A3B8" />
        </marker>
      </defs>
      <text x={w/2} y="34" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14.5" fontWeight="700" fill="#B45309">
        DM1 임상 적용 워크플로
      </text>
      {steps.map((s, i) => {
        const x = 40 + i * stepW;
        return (
          <g key={i}>
            <rect x={x} y="60" width={stepW - 20} height="80" rx="10" fill="#FFFFFF" stroke={s.color} strokeWidth="1.6" />
            <text x={x + (stepW - 20) / 2} y="88" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fontWeight="700" fill={s.color}>{s.label}</text>
            <text x={x + (stepW - 20) / 2} y="115" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#475569">{s.sub}</text>
            {i < steps.length - 1 && (
              <line x1={x + (stepW - 20) + 2} y1="100" x2={x + stepW - 2} y2="100" stroke="#94A3B8" strokeWidth="1.4" markerEnd="url(#wfArrow)" />
            )}
          </g>
        );
      })}
      <text x={w/2} y="180" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fontStyle="italic" fill="#64748B">
        Wet-lab · drug perturbation · spatial validation can follow after prospective clinical anchoring
      </text>
    </svg>
  );
};
