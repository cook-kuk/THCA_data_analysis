import React from "react";

/**
 * Current pathway (left) vs Proposed molecular triage (right).
 * Large SVG, projector-readable.
 */
export const ClinicalPathwayDiagram: React.FC = () => (
  <svg viewBox="0 0 1200 540" className="diagram" preserveAspectRatio="xMidYMid meet">
    <defs>
      <marker id="cpArrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#475569" />
      </marker>
      <marker id="cpArrowRed" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#B91C1C" />
      </marker>
      <marker id="cpArrowSage" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#047857" />
      </marker>
    </defs>

    {/* Headers */}
    <g>
      <rect x="40" y="34" width="540" height="40" rx="10" fill="#F1F5F9" stroke="#94A3B8"/>
      <text x="310" y="60" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="15" fontWeight="700" fill="#0F172A">현재 임상 경로 — Reactive</text>
      <rect x="620" y="34" width="540" height="40" rx="10" fill="#FEF3C7" stroke="#F59E0B"/>
      <text x="890" y="60" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="15" fontWeight="700" fill="#7C2D12">제안 임상 경로 — Pre-emptive molecular triage</text>
    </g>

    {/* ===== LEFT: current pathway ===== */}
    {[
      { y: 96, label: "수술", sub: "thyroidectomy" },
      { y: 158, label: "RAI 시험 투여", sub: "고용량 ¹³¹I" },
      { y: 220, label: "6 – 12 개월 추적", sub: "post-therapy WBS · Tg · imaging" },
      { y: 282, label: "Refractory 인지 (지연)", sub: "구조적 지속 / 재발 확인 시점" },
      { y: 344, label: "고용량 RAI 반복 또는 지연된 systemic 치료", sub: "골수억제 · leukopenia · 장기 독성 위험", danger: true }
    ].map((step, i, arr) => (
      <g key={i}>
        <rect x="80" y={step.y} width="460" height="46" rx="8"
              fill={step.danger ? "#FEF2F2" : "#FFFFFF"}
              stroke={step.danger ? "#B91C1C" : "#94A3B8"} strokeWidth={step.danger ? "1.6" : "1.2"}/>
        <text x="100" y={step.y + 22} fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill={step.danger ? "#7F1D1D" : "#0F172A"}>{i+1}. {step.label}</text>
        <text x="100" y={step.y + 38} fontFamily="Noto Sans KR" fontSize="11" fill={step.danger ? "#B91C1C" : "#475569"}>{step.sub}</text>
        {i < arr.length - 1 && (
          <line x1="310" y1={step.y + 46} x2="310" y2={arr[i+1].y - 2} stroke="#94A3B8" strokeWidth="1.4" markerEnd="url(#cpArrow)"/>
        )}
      </g>
    ))}
    <text x="310" y="424" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#B91C1C" fontWeight="600" fontStyle="italic">
      문제 — 평균 6 – 12 개월 지연 · 불필요한 고용량 RAI · 누적 독성
    </text>

    {/* ===== RIGHT: proposed pathway ===== */}
    {/* 1. surgery + assay */}
    <rect x="660" y="96" width="460" height="46" rx="8" fill="#FFFFFF" stroke="#0E7490" strokeWidth="1.4"/>
    <text x="680" y="118" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#0E7490">1. 수술 + 분자 분석 동시</text>
    <text x="680" y="134" fontFamily="Noto Sans KR" fontSize="11" fill="#155E75">FFPE · NGS / RNA / IHC assay</text>
    <line x1="890" y1="142" x2="890" y2="156" stroke="#0E7490" strokeWidth="1.4" markerEnd="url(#cpArrowSage)"/>

    {/* 2. 8-gene score */}
    <rect x="660" y="158" width="460" height="46" rx="8" fill="#FFFBEB" stroke="#B45309" strokeWidth="1.5"/>
    <text x="680" y="180" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#7C2D12">2. 8-gene Iodine-handling score</text>
    <text x="680" y="196" fontFamily="Noto Sans KR" fontSize="11" fill="#92400E">TG · TPO · TSHR · SLC5A5/NIS · DIO1 · PAX8 · NKX2-1 · FOXE1</text>
    <line x1="890" y1="204" x2="890" y2="218" stroke="#B45309" strokeWidth="1.4" markerEnd="url(#cpArrow)"/>

    {/* 3. Iodine-handling-high vs low */}
    <rect x="660" y="220" width="460" height="46" rx="8" fill="#F1F5F9" stroke="#334155" strokeWidth="1.4"/>
    <text x="890" y="248" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#0F172A">3. Iodine-handling-high (DM2)  vs  -low (DM1) 분류</text>

    {/* 4-A and 4-B split */}
    <line x1="780" y1="266" x2="780" y2="290" stroke="#1E40AF" strokeWidth="1.6" markerEnd="url(#cpArrowSage)"/>
    <line x1="1000" y1="266" x2="1000" y2="290" stroke="#B91C1C" strokeWidth="1.6" markerEnd="url(#cpArrowRed)"/>

    {/* DM2 path */}
    <rect x="640" y="292" width="280" height="68" rx="10" fill="#EFF6FF" stroke="#1E40AF" strokeWidth="1.6"/>
    <text x="780" y="316" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#1E40AF">DM2 · Iodine-handling-high</text>
    <text x="780" y="334" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#1E3A8A">전형적 RAI 치료 — 표준 ATA 가이드라인 적용</text>
    <text x="780" y="350" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#1E3A8A" fontStyle="italic">RAI-responsive-like</text>

    {/* DM1 path */}
    <rect x="940" y="292" width="240" height="68" rx="10" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="1.6"/>
    <text x="1060" y="316" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13.5" fontWeight="700" fill="#B91C1C">DM1 · Iodine-handling-low</text>
    <text x="1060" y="334" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#7F1D1D">불필요한 고용량 RAI 회피</text>
    <text x="1060" y="350" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#7F1D1D" fontStyle="italic">집중 추적 · 조기 systemic 논의</text>

    {/* benefit boxes */}
    <rect x="660" y="376" width="460" height="60" rx="8" fill="#ECFDF5" stroke="#047857" strokeWidth="1.5"/>
    <text x="890" y="398" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="13" fontWeight="700" fill="#047857">임상 가치</text>
    <text x="890" y="416" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#065F46">① 불필요한 고용량 RAI 누적 회피 → 골수독성 ↓</text>
    <text x="890" y="430" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#065F46">② 고위험군 조기 치료 전환 → progression-free 기간 ↑</text>

    {/* Bottom message */}
    <text x="600" y="496" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#B45309">
      핵심: "항암으로 빨리 보낸다"보다 먼저, RAI 가 흡수되지 않을 환자에게 불필요한 고용량 RAI 를 반복하지 않게 한다
    </text>
    <text x="600" y="518" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11.5" fill="#64748B" fontStyle="italic">
      Figure A · Current pathway (reactive) vs Proposed molecular triage (pre-emptive)
    </text>
  </svg>
);
