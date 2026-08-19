import React from "react";

/**
 * Clinician-friendly "same-group agreement %" chart replacing the raw ARI ladder.
 *
 * Message (임상의 이해 우선):
 *   "환자를 어떤 유전자 세트로 나눠도 같은 두 그룹이 나온다"
 *   → 우리 8-gene 이 만든 artifact 가 아니라 종양 안의 실제 분자 축.
 *
 * Bars show "how many % of patients land in the SAME DM1/DM2 group when the split is
 * re-computed using different gene sets" — intuitive to any clinician.
 * (ARI ~0.92 ≈ ~96 % agreement · ARI ~0.49 ≈ ~74 % · ARI ~0 ≈ ~50 % chance.)
 */
export const AriLadderChart: React.FC = () => {
  const data = [
    { label: "우리 8-gene panel  (기준)",         pct: 100, color: "#B91C1C", note: "정의상 100 %",           tier: "ref" },
    { label: "TIERA67  (분화 후보 67 gene)",       pct:  95, color: "#0E7490", note: "다른 후보 pool",         tier: "yes" },
    { label: "1,000-gene 자동 선택",              pct:  96, color: "#047857", note: "biology 무지 · variance", tier: "yes" },
    { label: "5,000-gene 비편향 pan-genome",       pct:  96, color: "#047857", note: "panel 완전 배제",         tier: "yes" },
    { label: "무작위 100-gene × 100 회 평균",       pct:  73, color: "#B45309", note: "무작위 sampling 대조군",  tier: "med" },
    { label: "BRAF / RAS driver 변이만",           pct:  50, color: "#94A3B8", note: "동전 던지기 수준",         tier: "no" }
  ];

  const W = 1200, H = 620;
  const padL = 380, padR = 140, padT = 200, padB = 100;
  const plotW = W - padL - padR;
  const plotH = H - padT - padB;
  const rowH = plotH / data.length;
  const scaleX = (p: number) => padL + (p / 100) * plotW;

  return (
    <svg viewBox={`0 0 ${W} ${H}`} className="diagram" preserveAspectRatio="xMidYMid meet">
      {/* ─── Title ─── */}
      <text x={W/2} y="34" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="16" fontWeight="800" fill="#7C2D12">
        환자를 어떤 유전자 세트로 나눠도 → 같은 두 그룹이 나온다
      </text>
      <text x={W/2} y="58" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#64748B" fontStyle="italic">
        즉 DM1 / DM2 는 우리 8 panel 이 만든 분류가 아니라, 종양 안에 이미 존재하는 분자 축.
      </text>

      {/* ─── Clinical analogy box ─── */}
      <g transform="translate(60,80)">
        <rect x="0" y="0" width={W - 120} height="94" rx="10" fill="#FEF3C7" stroke="#F59E0B" strokeWidth="1.5" />
        <text x="20" y="26" fontFamily="Noto Sans KR" fontSize="13" fontWeight="800" fill="#7C2D12">
          임상 비유 · 왜 이 실험이 중요한가
        </text>
        <text x="20" y="50" fontFamily="Noto Sans KR" fontSize="12" fill="#78350F">
          환자 100 명을 (a) 3 개 lab 로 두 그룹으로 나눔  →  (b) 30 개 lab 로 다시 나눔  →  (c) 300 개 lab 로 다시 나눔.
        </text>
        <text x="20" y="70" fontFamily="Noto Sans KR" fontSize="12" fill="#78350F">
          만약 (a) (b) (c) 모두 <tspan fontWeight="700" fill="#7C2D12">거의 같은 두 그룹으로 배정</tspan> 된다면 →
          "이 분류는 lab 선택 artifact 가 아니라 환자 자체의 실제 상태다".
        </text>
        <text x="20" y="88" fontFamily="Noto Sans KR" fontSize="11.5" fill="#78350F" fontStyle="italic">
          우리도 8 개 대신 67 · 1,000 · 5,000 gene 으로 다시 나눠봤다.  결과는 아래.
        </text>
      </g>

      {/* ─── Column headers ─── */}
      <text x={padL - 12} y={padT - 10} textAnchor="end" fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#0F172A">
        어떤 유전자 세트로 나눴는가
      </text>
      <text x={padL + plotW/2} y={padT - 10} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#0F172A">
        우리 8-gene 결과와 같은 그룹에 배정된 환자 비율 (%)
      </text>

      {/* ─── Bars ─── */}
      {data.map((d, i) => {
        const y = padT + i * rowH + 8;
        const barH = rowH - 20;
        const barW = scaleX(d.pct) - padL;
        return (
          <g key={d.label}>
            {/* Row label (Y axis) */}
            <text x={padL - 12} y={y + barH/2 + 5} textAnchor="end"
                  fontFamily="Noto Sans KR" fontSize="13" fontWeight={d.tier === "ref" ? 800 : 600}
                  fill={d.tier === "ref" ? "#B91C1C" : "#0F172A"}>
              {d.label}
            </text>

            {/* Gray background bar (full range) */}
            <rect x={padL} y={y} width={plotW} height={barH} rx="4" fill="#F1F5F9" />

            {/* Colored value bar */}
            <rect x={padL} y={y} width={Math.max(barW, 2)} height={barH} rx="4" fill={d.color}
                  opacity={d.tier === "no" ? 0.6 : 0.92} />

            {/* Percent label */}
            <text x={padL + barW + 8} y={y + barH/2 + 5} fontFamily="Inter" fontSize="16" fontWeight="800" fill={d.color}>
              {d.pct} %
            </text>

            {/* Note */}
            {d.note && (
              <text x={padL - 12} y={y + barH/2 + 22} textAnchor="end"
                    fontFamily="Noto Sans KR" fontSize="10" fill="#64748B" fontStyle="italic">
                {d.note}
              </text>
            )}
          </g>
        );
      })}

      {/* ─── Reference lines ─── */}
      <line x1={scaleX(50)} y1={padT - 4} x2={scaleX(50)} y2={H - padB + 4}
            stroke="#94A3B8" strokeDasharray="4 3" strokeWidth="1" />
      <text x={scaleX(50)} y={padT - 8} textAnchor="middle" fontFamily="Inter" fontSize="10" fill="#64748B">
        50 % (우연)
      </text>
      <line x1={scaleX(90)} y1={padT - 4} x2={scaleX(90)} y2={H - padB + 4}
            stroke="#047857" strokeDasharray="4 3" strokeWidth="1" />
      <text x={scaleX(90)} y={padT - 8} textAnchor="middle" fontFamily="Inter" fontSize="10" fill="#047857">
        90 %
      </text>

      {/* X axis */}
      <line x1={padL} y1={H - padB} x2={padL + plotW} y2={H - padB} stroke="#475569" strokeWidth="1" />
      {[0, 25, 50, 75, 100].map(t => (
        <g key={t}>
          <line x1={scaleX(t)} y1={H - padB} x2={scaleX(t)} y2={H - padB + 6} stroke="#475569" strokeWidth="1" />
          <text x={scaleX(t)} y={H - padB + 20} textAnchor="middle" fontFamily="Inter" fontSize="11" fill="#475569">{t}%</text>
        </g>
      ))}

      {/* ─── Conclusion band ─── */}
      <g transform={`translate(60,${H - 70})`}>
        <rect x="0" y="0" width={W - 120} height="52" rx="10" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="1.5" />
        <text x="20" y="22" fontFamily="Noto Sans KR" fontSize="13" fontWeight="800" fill="#B91C1C">
          결론 —  이 그룹 분류는 우리 8 panel 이 "발명"한 것이 아니다.
        </text>
        <text x="20" y="42" fontFamily="Noto Sans KR" fontSize="12" fill="#7F1D1D">
          종양 자체에 이미 존재하는 축을 <tspan fontWeight="700">8 개 유전자로 "측정" 할 뿐</tspan> —
          그래서 8 개만 봐도 이 축을 신뢰할 수 있고, 임상에서 rapid readout 이 가능하다.
        </text>
      </g>
    </svg>
  );
};
