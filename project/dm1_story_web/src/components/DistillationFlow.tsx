import React from "react";

/**
 * Funnel: Broad transcriptome ~20,000 genes → TIERA67 67-gene pool → 8-gene DM1 readout
 * with the panel members shown as a flowing column on the right.
 */
export const DistillationFlow: React.FC = () => (
  <svg viewBox="0 0 900 380" className="diagram" preserveAspectRatio="xMidYMid meet">
    <defs>
      <linearGradient id="funnelGrad" x1="0" x2="1">
        <stop offset="0%" stopColor="#94A3B8" />
        <stop offset="100%" stopColor="#B91C1C" />
      </linearGradient>
      <marker id="distArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#475569" />
      </marker>
    </defs>
    <text x="450" y="32" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14.5" fontWeight="700" fill="#B45309">
      Distillation flow · 광범위한 transcriptome → 8-유전자 DM1 readout
    </text>

    {/* Trapezoid stages from wide to narrow */}
    {/* Stage 1: full transcriptome */}
    <polygon points="80,80 540,80 480,150 140,150" fill="#F1F5F9" stroke="#94A3B8" strokeWidth="1.5"/>
    <text x="310" y="106" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#0F172A">Step 1.  Broad transcriptome</text>
    <text x="310" y="128" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#475569">≈ 20,000 protein-coding 유전자 · pan-genome 후보군</text>

    {/* Stage 2: TIERA67 candidate pool */}
    <polygon points="155,170 465,170 405,240 215,240" fill="#FEF3C7" stroke="#B45309" strokeWidth="1.5"/>
    <text x="310" y="196" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#7C2D12">Step 2.  TIERA67 후보 풀</text>
    <text x="310" y="218" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#7C2D12">67 thyroid lineage / RAI biology 유전자 · Yoo 2016 prior 기반</text>

    {/* Stage 3: 8-gene readout */}
    <polygon points="230,260 390,260 370,320 250,320" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="2"/>
    <text x="310" y="285" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#7F1D1D">Step 3.  8-gene DM1 readout</text>
    <text x="310" y="305" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#7F1D1D">측정 lens — DM1 / DM2 호출</text>

    {/* Right column: 8 panel genes */}
    <g>
      <rect x="600" y="80" width="220" height="240" rx="10" fill="#FFFFFF" stroke="#B91C1C" strokeWidth="1.5"/>
      <text x="710" y="105" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fontWeight="700" fill="#B91C1C">8 panel 유전자</text>
      {[
        ["TG", "thyroglobulin · effector"],
        ["TPO", "thyroperoxidase · effector"],
        ["TSHR", "TSH receptor · effector"],
        ["SLC5A5", "NIS · I⁻ symporter"],
        ["DIO1", "deiodinase · effector"],
        ["PAX8", "lineage TF"],
        ["NKX2-1", "lineage TF"],
        ["FOXE1", "lineage TF"]
      ].map(([g, role], i) => (
        <g key={g}>
          <text x="618" y={130 + i*22} fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#0F172A">{g}</text>
          <text x="668" y={130 + i*22} fontFamily="Noto Sans KR" fontSize="11" fill="#64748B">{role}</text>
        </g>
      ))}
    </g>

    {/* Connecting arrows */}
    <line x1="540" y1="115" x2="600" y2="115" stroke="#475569" strokeWidth="1.6" markerEnd="url(#distArrow)" />
    <text x="570" y="105" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#475569">측정</text>

    {/* Bottom caption */}
    <text x="310" y="350" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#64748B" fontStyle="italic">
      Reverse-causality lock — panel 은 Landa 2016 silenced list 와 무관하게 RAI 생물학 prior 만으로 선정.
    </text>
  </svg>
);
