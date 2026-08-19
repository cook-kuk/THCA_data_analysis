import React from "react";

/**
 * Clinical impact: ATA risk tier × DM cluster mosaic with DM1 over-representation in intermediate tier.
 * The "intermediate" tier is highlighted — RAI decision uncertainty zone.
 */
export const AtaImpactDiagram: React.FC = () => {
  // illustrative counts mirroring main Fig 6B mosaic
  const dm1 = { low: 18, mid: 92, high: 30 };
  const dm2 = { low: 120, mid: 195, high: 45 };
  const tot1 = dm1.low + dm1.mid + dm1.high;
  const tot2 = dm2.low + dm2.mid + dm2.high;
  const W = 900, H = 360;
  return (
    <svg viewBox={`0 0 ${W} ${H}`} className="diagram" preserveAspectRatio="xMidYMid meet">
      <text x={W/2} y="32" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14.5" fontWeight="700" fill="#B45309">
        ATA 2015 / 2025 위험군 × DM1 / DM2 분포  —  분자 분층화가 가장 가치 있는 환자군
      </text>

      {/* Bars */}
      {[
        { y: 90,  label: "DM1  (n ≈ 140)",  data: dm1, total: tot1, accent: "#B91C1C" },
        { y: 200, label: "DM2  (n ≈ 360)",  data: dm2, total: tot2, accent: "#1E40AF" }
      ].map((row, ri) => {
        const segs: [number, string, string][] = [
          [row.data.low,  "#047857", "low"],
          [row.data.mid,  "#B45309", "intermediate"],
          [row.data.high, "#B91C1C", "high"]
        ];
        let xx = 200;
        const wTotal = 600;
        return (
          <g key={row.label}>
            <text x={190} y={row.y + 38} textAnchor="end" fontFamily="Noto Sans KR" fontSize="13" fontWeight="700" fill={row.accent}>{row.label}</text>
            {segs.map(([n, col, tierLabel], i) => {
              const w = (n / row.total) * wTotal;
              const seg = (
                <g key={i}>
                  <rect x={xx} y={row.y} width={w} height={60} fill={col} stroke="#FFFFFF" strokeWidth="1.5"/>
                  <text x={xx + w/2} y={row.y + 32} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fontWeight="700" fill="#FFFFFF">{n}</text>
                  <text x={xx + w/2} y={row.y + 50} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#FFFFFF">{tierLabel}</text>
                </g>
              );
              xx += w;
              return seg;
            })}
            <text x={210 + wTotal + 20} y={row.y + 38} fontFamily="Noto Sans KR" fontSize="11" fill="#475569">
              {((row.data.mid / row.total) * 100).toFixed(0)} % intermediate
            </text>
          </g>
        );
      })}

      {/* Spotlight on intermediate tier */}
      <g>
        <rect x={395} y={75} width={235} height={195} rx="6" fill="none" stroke="#F59E0B" strokeWidth="3" strokeDasharray="6 4"/>
        <text x={510} y={300} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fontWeight="700" fill="#B45309">
          ATA intermediate-risk zone
        </text>
        <text x={510} y={320} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#7C2D12">
          RAI 치료 결정의 불확실성이 가장 큰 구간
        </text>
        <text x={510} y={340} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fontStyle="italic" fill="#0F172A">
          DM1 환자가 이 zone 에 과대표현 → 분자 분층화의 임상 가치 가장 큼
        </text>
      </g>
    </svg>
  );
};
