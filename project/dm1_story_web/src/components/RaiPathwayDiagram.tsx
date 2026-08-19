import React from "react";

/**
 * Thyrocyte cell-biology cartoon of iodine-handling circuit — no overlapping labels.
 * Layout: colloid (top) → apical membrane → cell (nucleus | ER/Golgi | DIO1) → basolateral membrane → blood.
 * Adapted from Kopp 2013 · Carvalho & Dupuy 2017 · DeGroot Thyroid Manager.
 */
export const RaiPathwayDiagram: React.FC = () => (
  <svg viewBox="0 0 1500 1150" className="diagram" preserveAspectRatio="xMidYMid meet" xmlns="http://www.w3.org/2000/svg">
    <defs>
      <linearGradient id="rpColloid" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stopColor="#FEF3C7" />
        <stop offset="100%" stopColor="#FDE68A" />
      </linearGradient>
      <linearGradient id="rpCell" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stopColor="#FFF7ED" />
        <stop offset="100%" stopColor="#FFEBD5" />
      </linearGradient>
      <linearGradient id="rpBlood" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stopColor="#FEE2E2" />
        <stop offset="100%" stopColor="#FECACA" />
      </linearGradient>
      <linearGradient id="rpNuc" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0%" stopColor="#DBEAFE" />
        <stop offset="100%" stopColor="#BFDBFE" />
      </linearGradient>
      <marker id="rp-ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#475569" />
      </marker>
      <marker id="rp-ah-r" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#B91C1C" />
      </marker>
      <marker id="rp-ah-b" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#1E40AF" />
      </marker>
      <marker id="rp-ah-t" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="#0E7490" />
      </marker>
    </defs>

    {/* ── Title (dedicated top band) ── */}
    <text x="750" y="34" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="18" fontWeight="700" fill="#7C2D12">
      갑상선 여포세포 · 방사성 요오드 흡수 세포 회로
    </text>
    <text x="750" y="60" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12.5" fill="#64748B" fontStyle="italic">
      Adapted from Kopp 2013 · Carvalho &amp; Dupuy 2017 · DeGroot Thyroid Manager
    </text>

    {/* ═══════ 1. COLLOID ═══════ */}
    <rect x="40" y="90" width="1420" height="140" rx="8" fill="url(#rpColloid)" stroke="#B45309" strokeWidth="1.8" />
    <text x="65" y="118" fontFamily="Noto Sans KR" fontSize="15" fontWeight="700" fill="#7C2D12">Follicular Lumen  (Colloid)</text>
    <text x="65" y="140" fontFamily="Noto Sans KR" fontSize="11.5" fill="#7C2D12" fontStyle="italic">Iodinated Tg · MIT / DIT · T3 / T4 저장고</text>

    {/* Iodinated Tg (left) */}
    <g transform="translate(400,182)">
      <ellipse cx="0" cy="0" rx="34" ry="20" fill="#B91C1C" opacity="0.7" />
      <circle cx="-22" cy="-8" r="5" fill="#F59E0B" />
      <circle cx="22" cy="-6" r="5" fill="#F59E0B" />
      <circle cx="-10" cy="12" r="5" fill="#F59E0B" />
      <circle cx="14" cy="10" r="5" fill="#F59E0B" />
      <text x="0" y="4" textAnchor="middle" fontFamily="Inter" fontSize="12" fontWeight="800" fill="#FFFFFF">Tg-I</text>
    </g>
    <text x="400" y="145" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#7F1D1D">이오드화 thyroglobulin</text>

    {/* MIT · DIT (center-left) */}
    <g transform="translate(680,182)">
      <circle cx="-11" cy="0" r="7" fill="#F59E0B" stroke="#B45309" strokeWidth="1" />
      <circle cx="11" cy="0" r="8" fill="#F59E0B" stroke="#B45309" strokeWidth="1" />
    </g>
    <text x="680" y="145" textAnchor="middle" fontFamily="Inter" fontSize="11.5" fontWeight="700" fill="#7C2D12">MIT · DIT</text>

    {/* T3 · T4 (center-right) */}
    <g transform="translate(920,182)">
      <ellipse cx="0" cy="0" rx="17" ry="11" fill="#F59E0B" opacity="0.9" />
      <text x="0" y="4" textAnchor="middle" fontFamily="Inter" fontSize="11" fontWeight="800" fill="#FFFFFF">T3 / T4</text>
    </g>
    <text x="920" y="145" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#7C2D12">갑상선 호르몬</text>

    {/* Free I⁻ (right) */}
    <g transform="translate(1220,182)">
      <circle cx="-18" cy="0" r="6" fill="#0E7490" />
      <circle cx="0" cy="0" r="6" fill="#0E7490" />
      <circle cx="18" cy="0" r="6" fill="#0E7490" />
    </g>
    <text x="1220" y="145" textAnchor="middle" fontFamily="Inter" fontSize="13" fontWeight="800" fill="#0E7490">I⁻</text>

    {/* ═══════ 2. APICAL MEMBRANE ═══════ */}
    <rect x="40" y="245" width="1420" height="42" fill="#FCA5A5" stroke="#B91C1C" strokeWidth="1.2" />
    <text x="65" y="272" fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#7F1D1D">Apical membrane</text>

    {/* TPO — panel */}
    <g transform="translate(490,266)">
      <rect x="-55" y="-22" width="110" height="44" rx="6" fill="#FED7AA" stroke="#B45309" strokeWidth="2.5" />
      <text x="0" y="-2" textAnchor="middle" fontFamily="Inter" fontSize="16" fontWeight="800" fill="#7C2D12">TPO</text>
      <text x="0" y="14" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fill="#7C2D12">peroxidase</text>
    </g>
    <text x="490" y="316" textAnchor="middle" fontFamily="Inter" fontSize="10.5" fontWeight="700" fill="#7C2D12" fontStyle="italic">◆ panel gene</text>

    {/* Pendrin — context */}
    <g transform="translate(750,266)">
      <rect x="-38" y="-22" width="76" height="44" rx="6" fill="#F1F5F9" stroke="#94A3B8" strokeWidth="1.2" />
      <text x="0" y="-2" textAnchor="middle" fontFamily="Inter" fontSize="12" fontWeight="600" fill="#475569">Pendrin</text>
      <text x="0" y="14" textAnchor="middle" fontFamily="Inter" fontSize="9.5" fill="#64748B">SLC26A4</text>
    </g>

    {/* DUOX2 — context */}
    <g transform="translate(970,266)">
      <rect x="-38" y="-22" width="76" height="44" rx="6" fill="#F1F5F9" stroke="#94A3B8" strokeWidth="1.2" />
      <text x="0" y="-2" textAnchor="middle" fontFamily="Inter" fontSize="12" fontWeight="600" fill="#475569">DUOX2</text>
      <text x="0" y="14" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="9.5" fill="#64748B">H₂O₂</text>
    </g>

    {/* I⁻ efflux annotation (dedicated label track above membrane in colloid area) */}
    <text x="750" y="222" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#0E7490" fontStyle="italic">Pendrin: 세포내 I⁻ → 콜로이드</text>

    {/* ═══════ 3. CYTOPLASM ═══════ */}
    <rect x="40" y="287" width="1420" height="540" fill="url(#rpCell)" stroke="#B45309" strokeWidth="1.8" />
    <text x="65" y="315" fontFamily="Noto Sans KR" fontSize="15" fontWeight="700" fill="#7C2D12">Thyrocyte (갑상선 여포세포)</text>
    <text x="65" y="335" fontFamily="Noto Sans KR" fontSize="11" fill="#7C2D12" fontStyle="italic">cytoplasm  ·  8 panel gene 대부분이 여기서 발현</text>

    {/* ─── NUCLEUS (top-left region) ─── */}
    <ellipse cx="290" cy="530" rx="205" ry="140" fill="url(#rpNuc)" stroke="#1E40AF" strokeWidth="2.8" />
    <ellipse cx="290" cy="530" rx="180" ry="118" fill="none" stroke="#3B82F6" strokeWidth="1" strokeDasharray="4 4" opacity="0.55" />
    <text x="290" y="418" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14" fontWeight="800" fill="#1E40AF">Nucleus</text>
    <text x="290" y="440" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#1E3A8A" fontStyle="italic">3 lineage TFs · co-activation</text>

    {/* TF cards - larger spacing */}
    <g transform="translate(290,486)">
      <rect x="-90" y="-20" width="180" height="34" rx="5" fill="#FFFFFF" stroke="#1E40AF" strokeWidth="1.5" />
      <text x="0" y="4" textAnchor="middle" fontFamily="Inter" fontSize="17" fontWeight="800" fill="#1E40AF">PAX8</text>
    </g>
    <g transform="translate(290,528)">
      <rect x="-90" y="-20" width="180" height="34" rx="5" fill="#FFFFFF" stroke="#1E40AF" strokeWidth="1.5" />
      <text x="0" y="4" textAnchor="middle" fontFamily="Inter" fontSize="17" fontWeight="800" fill="#1E40AF">NKX2-1</text>
    </g>
    <g transform="translate(290,570)">
      <rect x="-90" y="-20" width="180" height="34" rx="5" fill="#FFFFFF" stroke="#1E40AF" strokeWidth="1.5" />
      <text x="0" y="4" textAnchor="middle" fontFamily="Inter" fontSize="17" fontWeight="800" fill="#1E40AF">FOXE1</text>
    </g>
    <text x="290" y="622" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#1E3A8A" fontStyle="italic">→ TG · TPO · TSHR · NIS enhancer 공동 activation</text>

    {/* ─── TG blob (top-center) ─── */}
    <g transform="translate(720,400)">
      <circle cx="0" cy="0" r="26" fill="#B91C1C" opacity="0.5" />
      <circle cx="0" cy="0" r="18" fill="#B91C1C" opacity="0.85" />
      <text x="0" y="5" textAnchor="middle" fontFamily="Inter" fontSize="16" fontWeight="800" fill="#FFFFFF">TG</text>
    </g>
    <text x="720" y="356" textAnchor="middle" fontFamily="Inter" fontSize="11" fontWeight="700" fill="#B91C1C" fontStyle="italic">◆ panel gene</text>

    {/* ─── ER / Golgi (middle-center) ─── */}
    <g transform="translate(620,470)">
      <path d="M0 0 Q50 -6 100 0 Q150 6 200 0" fill="none" stroke="#B45309" strokeWidth="1.6" />
      <path d="M0 18 Q50 12 100 18 Q150 24 200 18" fill="none" stroke="#B45309" strokeWidth="1.6" />
      <path d="M0 36 Q50 30 100 36 Q150 42 200 36" fill="none" stroke="#B45309" strokeWidth="1.6" />
      <path d="M0 54 Q50 48 100 54 Q150 60 200 54" fill="none" stroke="#B45309" strokeWidth="1.6" />
      <text x="100" y="82" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#B45309">ER / Golgi</text>
      <text x="100" y="98" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#B45309" fontStyle="italic">TG 합성 → 소포 → apical 분비</text>
    </g>

    {/* ─── DIO1 (right) ─── */}
    <g transform="translate(1180,560)">
      <rect x="-100" y="-30" width="200" height="60" rx="8" fill="#FED7AA" stroke="#7C2D12" strokeWidth="2.5" />
      <text x="0" y="-4" textAnchor="middle" fontFamily="Inter" fontSize="17" fontWeight="800" fill="#7C2D12">DIO1</text>
      <text x="0" y="16" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="11" fill="#7C2D12">T4 → T3 deiodinase</text>
    </g>
    <text x="1180" y="518" textAnchor="middle" fontFamily="Inter" fontSize="11" fontWeight="700" fill="#7C2D12" fontStyle="italic">◆ panel gene</text>

    {/* ─── Cytoplasmic I⁻ (small, near center-lower) ─── */}
    <g transform="translate(880,720)">
      <circle cx="-14" cy="0" r="5" fill="#0E7490" />
      <circle cx="0" cy="0" r="5" fill="#0E7490" />
      <text x="20" y="4" fontFamily="Inter" fontSize="12" fontWeight="800" fill="#0E7490">I⁻</text>
      <text x="0" y="-16" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fill="#0E7490" fontStyle="italic">세포내 I⁻</text>
    </g>

    {/* ─── ARROWS (clean, well-separated paths) ─── */}
    {/* TF nucleus → TG synthesis (up-right curve) */}
    <path d="M475 480 Q580 440 700 415" fill="none" stroke="#1E40AF" strokeWidth="1.8" strokeDasharray="6 4" markerEnd="url(#rp-ah-b)" opacity="0.75" />
    <text x="580" y="435" fontFamily="Noto Sans KR" fontSize="10" fill="#1E40AF" fontStyle="italic">→ TG</text>

    {/* TF nucleus → TPO (up to apical) */}
    <path d="M420 420 Q450 360 490 300" fill="none" stroke="#1E40AF" strokeWidth="1.8" strokeDasharray="6 4" markerEnd="url(#rp-ah-b)" opacity="0.75" />
    <text x="440" y="380" fontFamily="Noto Sans KR" fontSize="10" fill="#1E40AF" fontStyle="italic">→ TPO</text>

    {/* TF nucleus → TSHR (down-right to basolateral) */}
    <path d="M460 630 Q490 730 520 800" fill="none" stroke="#1E40AF" strokeWidth="1.8" strokeDasharray="6 4" markerEnd="url(#rp-ah-b)" opacity="0.75" />
    <text x="480" y="720" fontFamily="Noto Sans KR" fontSize="10" fill="#1E40AF" fontStyle="italic">→ TSHR</text>

    {/* TF nucleus → NIS (long down-right) */}
    <path d="M495 610 Q650 700 900 800" fill="none" stroke="#1E40AF" strokeWidth="1.8" strokeDasharray="6 4" markerEnd="url(#rp-ah-b)" opacity="0.6" />
    <text x="700" y="740" fontFamily="Noto Sans KR" fontSize="10" fill="#1E40AF" fontStyle="italic">→ NIS</text>

    {/* TG → colloid apical secretion */}
    <path d="M720 375 Q720 320 710 290" fill="none" stroke="#B91C1C" strokeWidth="2.2" markerEnd="url(#rp-ah-r)" />
    <text x="750" y="330" fontFamily="Noto Sans KR" fontSize="10.5" fill="#7F1D1D" fontStyle="italic">apical secretion</text>

    {/* T4 endocytosis from colloid to DIO1 */}
    <g transform="translate(1180,260)">
      <ellipse cx="0" cy="0" rx="14" ry="9" fill="#F59E0B" />
      <text x="0" y="4" textAnchor="middle" fontFamily="Inter" fontSize="10" fontWeight="800" fill="#FFFFFF">T4</text>
    </g>
    <path d="M1180 275 L1180 528" fill="none" stroke="#F59E0B" strokeWidth="2" strokeDasharray="6 4" markerEnd="url(#rp-ah)" />
    <text x="1200" y="410" fontFamily="Noto Sans KR" fontSize="10.5" fill="#B45309" fontStyle="italic">T4 endocytosis</text>

    {/* DIO1 → T3 exit downward */}
    <path d="M1180 595 L1180 680 L1300 780" fill="none" stroke="#F59E0B" strokeWidth="2" markerEnd="url(#rp-ah)" />
    <text x="1250" y="670" fontFamily="Noto Sans KR" fontSize="10.5" fill="#B45309" fontStyle="italic">→ T3 활성</text>

    {/* I⁻ trajectory: basolateral NIS → cytoplasm → apical Pendrin (long curved arrow) */}
    <path d="M880 795 Q880 750 850 720" fill="none" stroke="#0E7490" strokeWidth="1.6" strokeDasharray="4 3" markerEnd="url(#rp-ah-t)" opacity="0.85" />
    <path d="M870 715 Q800 500 760 288" fill="none" stroke="#0E7490" strokeWidth="1.6" strokeDasharray="4 3" markerEnd="url(#rp-ah-t)" opacity="0.7" />

    {/* ═══════ 4. BASOLATERAL MEMBRANE ═══════ */}
    <rect x="40" y="827" width="1420" height="42" fill="#FCA5A5" stroke="#B91C1C" strokeWidth="1.2" />
    <text x="65" y="854" fontFamily="Noto Sans KR" fontSize="12" fontWeight="700" fill="#7F1D1D">Basolateral membrane</text>

    {/* TSHR — panel */}
    <g transform="translate(520,848)">
      <rect x="-65" y="-22" width="130" height="44" rx="6" fill="#BFDBFE" stroke="#1E40AF" strokeWidth="2.5" />
      <text x="0" y="-2" textAnchor="middle" fontFamily="Inter" fontSize="16" fontWeight="800" fill="#1E40AF">TSHR</text>
      <text x="0" y="14" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fill="#1E40AF">TSH receptor</text>
    </g>
    <text x="520" y="898" textAnchor="middle" fontFamily="Inter" fontSize="10.5" fontWeight="700" fill="#1E40AF" fontStyle="italic">◆ panel gene</text>

    {/* Gs / cAMP / PKA annotation (own dedicated card, above basolateral membrane) */}
    <g transform="translate(520,790)">
      <rect x="-85" y="-18" width="170" height="34" rx="5" fill="#FFFFFF" stroke="#0E7490" strokeWidth="1.2" strokeDasharray="4 3" />
      <text x="0" y="0" textAnchor="middle" fontFamily="Inter" fontSize="11.5" fontWeight="700" fill="#0E7490">Gs → cAMP → PKA</text>
      <text x="0" y="12" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="9.5" fill="#0E7490">신호 증폭</text>
    </g>

    {/* NIS / SLC5A5 — panel */}
    <g transform="translate(900,848)">
      <rect x="-80" y="-22" width="160" height="44" rx="6" fill="#FED7AA" stroke="#B45309" strokeWidth="2.5" />
      <text x="0" y="-2" textAnchor="middle" fontFamily="Inter" fontSize="16" fontWeight="800" fill="#7C2D12">NIS  (SLC5A5)</text>
      <text x="0" y="14" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10" fill="#7C2D12">Na⁺ / I⁻ symporter</text>
    </g>
    <text x="900" y="898" textAnchor="middle" fontFamily="Inter" fontSize="10.5" fontWeight="700" fill="#7C2D12" fontStyle="italic">◆ panel gene</text>

    {/* ═══════ 5. BLOOD CAPILLARY ═══════ */}
    <rect x="40" y="920" width="1420" height="130" fill="url(#rpBlood)" stroke="#B91C1C" strokeWidth="1.8" />
    <text x="65" y="948" fontFamily="Noto Sans KR" fontSize="14" fontWeight="700" fill="#7F1D1D">Blood capillary</text>
    <text x="65" y="968" fontFamily="Noto Sans KR" fontSize="11" fill="#7F1D1D" fontStyle="italic">TSH (뇌하수체) + I⁻ + Na⁺ (섭취) 공급 · T3 / T4 방출</text>

    {/* TSH molecules → TSHR */}
    <g transform="translate(520,1010)">
      <ellipse cx="-14" cy="0" rx="12" ry="9" fill="#F59E0B" stroke="#B45309" strokeWidth="1.2" />
      <ellipse cx="14" cy="0" rx="12" ry="9" fill="#F59E0B" stroke="#B45309" strokeWidth="1.2" />
      <text x="0" y="26" textAnchor="middle" fontFamily="Inter" fontSize="12" fontWeight="800" fill="#7C2D12">TSH</text>
    </g>
    <path d="M520 998 L520 872" fill="none" stroke="#B45309" strokeWidth="2.4" markerEnd="url(#rp-ah)" />

    {/* Na⁺ + I⁻ → NIS */}
    <g transform="translate(900,1010)">
      <circle cx="-20" cy="0" r="6" fill="#0E7490" />
      <circle cx="0" cy="0" r="6" fill="#0E7490" />
      <circle cx="20" cy="0" r="6" fill="#0E7490" />
      <text x="0" y="26" textAnchor="middle" fontFamily="Inter" fontSize="12" fontWeight="800" fill="#0E7490">I⁻ + Na⁺</text>
    </g>
    <path d="M900 998 L900 872" fill="none" stroke="#0E7490" strokeWidth="2.4" markerEnd="url(#rp-ah)" />

    {/* T3 / T4 exit */}
    <g transform="translate(1300,1010)">
      <path d="M0 -140 L0 -20" fill="none" stroke="#B45309" strokeWidth="2.4" markerEnd="url(#rp-ah)" />
      <ellipse cx="0" cy="6" rx="16" ry="10" fill="#F59E0B" />
      <text x="0" y="10" textAnchor="middle" fontFamily="Inter" fontSize="10" fontWeight="800" fill="#FFFFFF">T3 · T4</text>
      <text x="0" y="34" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="10.5" fill="#7C2D12" fontStyle="italic">전신 순환</text>
    </g>

    {/* ═══════ 6. Clinical takeaway ═══════ */}
    <g transform="translate(40,1078)">
      <rect x="0" y="0" width="1420" height="60" rx="10" fill="#FEF2F2" stroke="#B91C1C" strokeWidth="2" />
      <text x="710" y="26" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="14.5" fontWeight="800" fill="#B91C1C">
        ¹³¹I RAI 는 이 회로를 사용한다 — NIS 흡수 · TPO 산화 · TG 통합
      </text>
      <text x="710" y="48" textAnchor="middle" fontFamily="Noto Sans KR" fontSize="12" fill="#7F1D1D">
        DM1 에서 8 panel gene 동시 침묵 → 회로 shutdown → RAI 흡수 실패 → 임상 = 불필요한 고용량 RAI 회피
      </text>
    </g>
  </svg>
);
