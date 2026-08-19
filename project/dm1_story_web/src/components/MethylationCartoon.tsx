import React from "react";

/**
 * DM2 vs DM1 promoter cartoon:
 *   DM2 = open promoter (low β), TFs binding, expression on
 *   DM1 = methylated promoter (high β), TFs blocked, expression off
 */
export const MethylationCartoon: React.FC = () => (
  <svg viewBox="0 0 900 360" className="diagram" preserveAspectRatio="xMidYMid meet">
    <defs>
      <marker id="mcArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#475569" />
      </marker>
    </defs>

    <text x="450" y="32" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14.5" fontWeight="700" fill="#B45309">
      DM2 vs DM1 promoter — 후성유전 침묵의 분자 cartoon
    </text>

    {/* ===== LEFT: DM2 (open promoter) ===== */}
    <g>
      <rect x="40" y="60" width="380" height="240" rx="14" fill="#EFF6FF" stroke="#1E40AF" strokeWidth="1.5"/>
      <text x="230" y="86" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#1E40AF">DM2 (분화 유지)</text>

      {/* DNA strand */}
      <line x1="70" y1="170" x2="390" y2="170" stroke="#0F172A" strokeWidth="3"/>
      <line x1="70" y1="184" x2="390" y2="184" stroke="#0F172A" strokeWidth="3"/>

      {/* CpG sites — open circles (unmethylated) */}
      {[120, 160, 200, 240, 280, 320, 360].map(x => (
        <circle key={x} cx={x} cy="177" r="6" fill="#FFFFFF" stroke="#1E40AF" strokeWidth="1.5"/>
      ))}
      <text x="230" y="220" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#1E3A8A">CpG unmethylated · β ≈ 0.25</text>

      {/* TF binding */}
      <ellipse cx="180" cy="138" rx="22" ry="14" fill="#1E40AF"/>
      <text x="180" y="142" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fontWeight="700" fill="#FFFFFF">PAX8</text>
      <ellipse cx="260" cy="138" rx="22" ry="14" fill="#1E40AF"/>
      <text x="260" y="142" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fontWeight="700" fill="#FFFFFF">FOXE1</text>
      <line x1="180" y1="152" x2="180" y2="167" stroke="#1E40AF" strokeWidth="1.5"/>
      <line x1="260" y1="152" x2="260" y2="167" stroke="#1E40AF" strokeWidth="1.5"/>

      {/* RNA produced */}
      <text x="230" y="248" textAnchor="middle" fontSize="16" fill="#047857" fontWeight="700">∿∿∿∿  mRNA  ∿∿∿∿</text>
      <text x="230" y="272" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#047857" fontWeight="600">표적 유전자 발현 ON</text>
      <text x="230" y="290" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#047857" fontStyle="italic">TG · TPO · TSHR · NIS · DIO1 활성</text>
    </g>

    {/* ===== RIGHT: DM1 (methylated promoter) ===== */}
    <g>
      <rect x="480" y="60" width="380" height="240" rx="14" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="1.5"/>
      <text x="670" y="86" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#B91C1C">DM1 (분화 손실)</text>

      {/* DNA strand */}
      <line x1="510" y1="170" x2="830" y2="170" stroke="#0F172A" strokeWidth="3"/>
      <line x1="510" y1="184" x2="830" y2="184" stroke="#0F172A" strokeWidth="3"/>

      {/* CpG sites — filled circles (methylated) */}
      {[560, 600, 640, 680, 720, 760, 800].map(x => (
        <g key={x}>
          <circle cx={x} cy="177" r="6" fill="#B91C1C" stroke="#7F1D1D" strokeWidth="1.5"/>
          <text x={x} y="200" textAnchor="middle" fontSize="9" fontWeight="700" fill="#B91C1C">CH₃</text>
        </g>
      ))}
      <text x="670" y="226" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#7F1D1D">CpG methylated · β ≈ 0.39</text>

      {/* TF blocked - "X" through them */}
      <g opacity="0.45">
        <ellipse cx="620" cy="138" rx="22" ry="14" fill="#94A3B8"/>
        <text x="620" y="142" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fontWeight="700" fill="#FFFFFF">PAX8</text>
        <ellipse cx="700" cy="138" rx="22" ry="14" fill="#94A3B8"/>
        <text x="700" y="142" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fontWeight="700" fill="#FFFFFF">FOXE1</text>
      </g>
      {/* X marks */}
      <text x="620" y="143" textAnchor="middle" fontFamily="Inter" fontSize="22" fontWeight="900" fill="#B91C1C">✕</text>
      <text x="700" y="143" textAnchor="middle" fontFamily="Inter" fontSize="22" fontWeight="900" fill="#B91C1C">✕</text>
      <text x="660" y="115" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#7F1D1D" fontStyle="italic">TF 결합 차단</text>

      {/* No RNA */}
      <text x="670" y="248" textAnchor="middle" fontSize="18" fill="#94A3B8" fontWeight="700">— · — · —</text>
      <text x="670" y="272" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#B91C1C" fontWeight="600">표적 유전자 발현 OFF</text>
      <text x="670" y="290" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#7F1D1D" fontStyle="italic">RAI 흡수 능력 ↓ · 분화 손실</text>
    </g>

    {/* Bottom caption */}
    <text x="450" y="330" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fontWeight="600" fill="#0F172A">
      DM1 종양은 promoter hypermethylation 으로 분화 TF 결합과 RAI 효소 발현이 동시 차단된다.
    </text>
    <text x="450" y="350" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#64748B" fontStyle="italic">
      TCGA HM450 n = 503  ·  TPO Cohen's d = 2.30, p = 1.9 × 10⁻¹⁸  ·  mean 8-gene β  DM1 0.385 vs DM2 0.253
    </text>
  </svg>
);
