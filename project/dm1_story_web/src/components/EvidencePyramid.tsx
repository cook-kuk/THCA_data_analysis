import React from "react";

/**
 * Evidence stack pyramid — bottom = broad biology, top = clinical translation.
 */
export const EvidencePyramid: React.FC = () => {
  const layers = [
    { y: 80,  w: 980, label: "Clinically translatable RAI decision support", sub: "SNUBH validation pending", color: "#B91C1C", text: "#FFFFFF", italic: true },
    { y: 130, w: 920, label: "Single-cell support — Lu 2023 thyrocyte-intrinsic gradient (n = 14,624)", color: "#B45309", text: "#FFFFFF" },
    { y: 180, w: 860, label: "Proteogenomic support — Mun 2025 (n = 336, d = −1.91, 7/7 sign-consistent)", color: "#B45309", text: "#FFFFFF" },
    { y: 230, w: 800, label: "FFPE compatibility — KS p = 0.44 (Lee FFPE n = 632 vs TCGA FF)", color: "#0E7490", text: "#FFFFFF" },
    { y: 280, w: 740, label: "External cohorts — 19 cohorts · master forest mean d = 2.81", color: "#0E7490", text: "#FFFFFF" },
    { y: 330, w: 680, label: "Epigenetic mechanism — TPO d = 2.30, p = 1.9 × 10⁻¹⁸", color: "#1E40AF", text: "#FFFFFF" },
    { y: 380, w: 620, label: "Driver-orthogonal analysis — BRAF/RAS/TERT AUC ≈ 0.5", color: "#1E40AF", text: "#FFFFFF" },
    { y: 430, w: 560, label: "TCGA discovery — n = 504, DM1 28.4 %", color: "#334155", text: "#FFFFFF" },
    { y: 480, w: 500, label: "Literature-guided biology — Yoo 2016 RAI prior · TIERA67", color: "#334155", text: "#FFFFFF" }
  ];

  return (
    <svg viewBox="0 0 1200 580" className="diagram" preserveAspectRatio="xMidYMid meet">
      <text x="600" y="38" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="16" fontWeight="700" fill="#B45309">
        Figure C · Evidence stack — biology 에서 임상 적용까지 누적 증거 피라미드
      </text>

      {layers.map((L, i) => (
        <g key={i}>
          <rect x={600 - L.w/2} y={L.y} width={L.w} height="42" rx="8"
                fill={L.color} stroke="#FFFFFF" strokeWidth="2"/>
          <text x="600" y={L.y + 27} textAnchor="middle"
                fontFamily="Noto Sans KR" fontSize="13" fontWeight="700"
                fill={L.text} fontStyle={L.italic ? "italic" : "normal"}>
            {L.label}
          </text>
          {L.sub && (
            <text x="600" y={L.y + 56} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#B45309" fontStyle="italic">
              {L.sub}
            </text>
          )}
        </g>
      ))}

      <text x="600" y="556" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#64748B" fontStyle="italic">
        아래 layer 가 더 광범위한 biology, 위 layer 로 갈수록 임상 적용에 가까워짐.  각 layer 모두 독립 증거.
      </text>
    </svg>
  );
};
