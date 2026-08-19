import React from "react";

/**
 * Now / Next / Later validation roadmap.
 */
export const ValidationRoadmap: React.FC = () => {
  const lanes = [
    {
      title: "NOW",
      subtitle: "지금까지 완료",
      color: "#047857",
      bg: "#ECFDF5",
      items: [
        "Public cohort discovery (TCGA n=504)",
        "19 외부 코호트 cross-validation",
        "8-gene symbol / alias HGNC 통일",
        "Epigenetic mechanism (HM450 n=503)",
        "Proteogenomic 검증 (Mun 2025 n=336)",
        "Manuscript draft + figure 정리",
      ]
    },
    {
      title: "NEXT",
      subtitle: "이번 미팅 후 즉시",
      color: "#B45309",
      bg: "#FFFBEB",
      items: [
        "SNUBH NGS panel 에 8 gene 포함 여부 확인",
        "갑상선암 NGS 환자 수 count",
        "RAI 치료 이력 / refractory 정의",
        "약 30 명 retrospective validation 설계",
        "IRB 초안 작성",
        "4 / 6-gene reduced model (NGS subset 대비)"
      ]
    },
    {
      title: "LATER",
      subtitle: "확장 / 임상 단계",
      color: "#B91C1C",
      bg: "#FEF2F2",
      items: [
        "FFPE block 8-marker IHC validation",
        "H-score / intensity / positive % scoring",
        "병리 협업 (cancer vs stromal staining)",
        "진단법 / IHC 기반 특허 검토",
        "Prospective validation cohort 모집",
        "Spatial transcriptomics + H&E AI (2 차)"
      ]
    }
  ];

  return (
    <svg viewBox="0 0 1200 540" className="diagram" preserveAspectRatio="xMidYMid meet">
      <text x="600" y="34" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="16" fontWeight="700" fill="#B45309">
        Figure E · Validation roadmap — Now / Next / Later
      </text>

      {lanes.map((L, i) => {
        const x = 40 + i * 390;
        return (
          <g key={i}>
            {/* lane header */}
            <rect x={x} y={64} width={370} height={64} rx={10} fill={L.color}/>
            <text x={x + 185} y={92} textAnchor="middle" fontFamily="Inter" fontSize="22" fontWeight="900" fill="#FFFFFF" letterSpacing="3">{L.title}</text>
            <text x={x + 185} y={114} textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#FFFFFF" opacity="0.9">{L.subtitle}</text>

            {/* lane body */}
            <rect x={x} y={132} width={370} height={380} rx={10} fill={L.bg} stroke={L.color} strokeWidth="1.4"/>
            {L.items.map((it, j) => (
              <g key={j}>
                <circle cx={x + 22} cy={158 + j * 56} r={5} fill={L.color}/>
                <text x={x + 40} y={163 + j * 56} fontFamily="Noto Sans KR" fontSize="12.5" fill="#0F172A">{it}</text>
              </g>
            ))}
          </g>
        );
      })}
    </svg>
  );
};
