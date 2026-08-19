import React from "react";

/**
 * Two-column unsafe vs safe claim table — visual.
 */
export const ClaimSafetyLadder: React.FC = () => {
  const pairs = [
    {
      bad: "이 8 개 gene 으로 RAI response 를 완벽히 예측한다",
      good: "RAI avidity / refractory-like biology 와 관련된 molecular state 를 분류한다"
    },
    {
      bad: "이건 RAI response predictor 다",
      good: "Iodine-handling / thyroid differentiation state classifier 다"
    },
    {
      bad: "Treatment-selection biomarker 로 사용 가능하다",
      good: "Risk stratification axis · ATA 중간 위험군의 RAI 결정 보조"
    },
    {
      bad: "즉시 임상 적용 가능하다",
      good: "Retrospective + external validation 기반 · prospective / IHC validation 필요"
    },
    {
      bad: "TCGA + MSK HR 2.5 = DM1 환자는 더 빨리 죽는다",
      good: "후향적 통합 HR 2.53 [1.31, 4.89] · advanced cohort 기여 큼 · primary PTC 단독 underpowered"
    }
  ];

  return (
    <svg viewBox="0 0 1200 600" className="diagram" preserveAspectRatio="xMidYMid meet">
      <text x="600" y="34" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="16" fontWeight="700" fill="#B45309">
        Figure D · Claim safety ladder — 무엇을 주장하면 안 되고, 무엇을 주장해야 하는가
      </text>

      {/* Headers */}
      <rect x="40" y="64" width="560" height="40" rx="8" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="1.5"/>
      <text x="320" y="89" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#B91C1C">⚠ 위험한 표현 (Reviewer attack 유발)</text>
      <rect x="620" y="64" width="540" height="40" rx="8" fill="#ECFDF5" stroke="#047857" strokeWidth="1.5"/>
      <text x="890" y="89" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#047857">✓ 안전한 표현 (Manuscript / 발표 권장)</text>

      {pairs.map((p, i) => {
        const y = 120 + i * 90;
        return (
          <g key={i}>
            <rect x="40" y={y} width="560" height="78" rx="8" fill="#FFFFFF" stroke="#FCA5A5" strokeWidth="1"/>
            <text x="60" y={y + 30} fontFamily="Noto Sans KR" fontSize="13" fill="#7F1D1D" fontWeight="600">✕  {p.bad}</text>

            {/* arrow */}
            <text x="610" y={y + 44} textAnchor="middle" fontFamily="Inter" fontSize="20" fontWeight="700" fill="#475569">→</text>

            <rect x="620" y={y} width="540" height="78" rx="8" fill="#FFFFFF" stroke="#6EE7B7" strokeWidth="1"/>
            <text x="640" y={y + 30} fontFamily="Noto Sans KR" fontSize="13" fill="#065F46" fontWeight="600">✓  {p.good}</text>
          </g>
        );
      })}
    </svg>
  );
};
