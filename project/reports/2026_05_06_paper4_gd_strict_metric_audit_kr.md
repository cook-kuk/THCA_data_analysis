# Paper 4 — Korean GD HLA Pan-Asian: strict-metric 정밀 감사 (KR)
**Date:** 2026-05-06 · **Author:** Seungho Cook
**Scope:** Paper 4 = Korean GD HLA backlog (4/4 gating, Pan-Asian meta path)
**Discipline:** Multi-test 보정 + cross-paper direction audit + Korean replication readiness 까지만. Final HLA association 주장 금지. Risk allele 주장 금지. Paper 1 미접촉.

---

## 0. TL;DR

- **Chu 2018 GD anchor (n = 1468 case vs 1490 control, Chinese GD)**: 6 alleles 모두 **BH-q < 0.05 ✓ AND Bonferroni p < 0.05 ✓** — 강한 anchor signal.
- 4 alleles GD-enriched (A*02:07, B*46:01, C*01:02, DPB1*05:01), 2 alleles GD-depleted/protective (DQB1*02:01, DRB1*07:01).
- Cross-paper direction (Paper 4 GD vs Paper 2 PTC S3 strict allele-vs-allele):
  - **BOTH same direction:** A*02:07 (enrich), B*46:01 (enrich), DQB1*02:01 (depletion), DRB1*07:01 (depletion).
  - **OPPOSITE direction:** **C*01:02** (GD enrich vs PTC depletion) and **DPB1*05:01** (GD enrich vs PTC ~null).
- Pan-Asian meta status: 3 alleles still **single_source_anchor_only** (A*02:07, DQB1*02:01, DRB1*07:01). 3 alleles screening_random_effects (B*46:01 k=3, C*01:02 k=2, DPB1*05:01 k=3) but Korean source 없음.
- **Korean GD case-control NGS HLA cohort = NOT YET ASSEMBLED**. 이것이 Paper 4 backlog gating 4/4 의 핵심.

본 패킷은 Chu 2018 + Pan-Asian meta까지의 strengthening. Korean replication 필요.

---

## 1. Chu 2018 anchor — strict multi-test (within 6 candidates)

`paper4_gd_chu2018_bh_fdr.tsv` 참조.

| Allele | OR | 95% CI | raw p | **BH-q** | **Bonferroni p** | Interpretation |
|---|---|---|---|---|---|---|
| A*02:07 | 2.10 | 1.70-2.59 | 2.07×10⁻¹² | **2.07×10⁻¹² ✓** | **1.24×10⁻¹¹ ✓** | GD-enriched |
| B*46:01 | 2.38 | 1.99-2.86 | 8.78×10⁻²¹ | **1.76×10⁻²⁰ ✓** | **5.27×10⁻²⁰ ✓** | GD-enriched |
| C*01:02 | 1.83 | 1.57-2.12 | 5.40×10⁻¹⁵ | **8.10×10⁻¹⁵ ✓** | **3.24×10⁻¹⁴ ✓** | GD-enriched |
| DPB1*05:01 | 1.90 | 1.69-2.14 | 1.73×10⁻²⁶ | **1.04×10⁻²⁵ ✓** | **1.04×10⁻²⁵ ✓** | ★ Asian Graves' top |
| DQB1*02:01 | 0.57 | 0.49-0.66 | 2.31×10⁻¹³ | **2.77×10⁻¹³ ✓** | **1.39×10⁻¹² ✓** | GD-depleted/protective |
| DRB1*07:01 | 0.43 | 0.36-0.51 | 2.49×10⁻²¹ | **7.47×10⁻²¹ ✓** | **1.49×10⁻²⁰ ✓** | GD-depleted/protective |

**모든 alleles BH-q AND Bonferroni 0.05 통과**. 단일 anchor (Chu 2018 China GD) 내에서는 robust.

`P4_F17_chu2018_bh_bonferroni.png` — forest + multi-test 표시.

---

## 2. Cross-paper direction — Paper 4 GD vs Paper 2 PTC S3 (strict allele)

같은 6 alleles, 서로 다른 disease (GD vs PTC), 서로 다른 데이터 (Chu 2018 GD case-control vs Korean PTC pool S3 strict allele-vs-allele).

| Allele | GD OR (Chu 2018) | GD direction | PTC S3 OR | PTC S3 direction | 일치? |
|---|---|---|---|---|---|
| A*02:07 | 2.10 (enrich strong) | enrichment | 1.26 (enrich weak) | enrichment | **BOTH enrichment** ✓ |
| B*46:01 | 2.38 (enrich strong) | enrichment | 1.05 (~null) | enrichment | **BOTH enrichment** ✓ |
| C*01:02 | 1.83 (enrich strong) | enrichment | 0.68 (depletion strict q=7e-4) | depletion | **OPPOSITE** ⚠ |
| DPB1*05:01 | 1.90 (enrich strong) | enrichment | 0.98 (~null) | depletion | **OPPOSITE** ⚠ |
| DQB1*02:01 | 0.57 (depleted) | depletion | 0.018 (zero-cell SUSPENDED) | depletion | **BOTH depletion** (PTC suspended by QC audit) |
| DRB1*07:01 | 0.43 (depleted) | depletion | 0.84 (depletion-trend) | depletion | **BOTH depletion** ✓ |

`P4_F18_gd_vs_ptc_direction.png` 참조.

### 해석 (signature direction 차원)

- **Coincident direction (4/6)**: A*02:07, B*46:01, DQB1*02:01, DRB1*07:01 — Asian Graves'와 Korean PTC가 같은 면역-유전적 background와 일관.
- **Discordant direction (2/6)**: C*01:02, DPB1*05:01 — GD에서는 enrichment지만 PTC strict 비교에서는 ~null 또는 depletion.
- **Yu 2026-05-04 분리벽 (separation wall)**: Paper 2 = HT-isolated, Paper 4 = GD. 이 cross-paper direction 분석은 Paper 2/Paper 4가 *서로 다른* 질환임을 보여주는 데이터로만 사용. Paper 2 → Paper 4 추론, 또는 그 반대 추론은 절대 금지.

---

## 3. Pan-Asian meta readiness

`paper4_gd_korean_replication_readiness.tsv` 참조.

| Allele | Sources | k | Pooled OR | p_random | Meta status |
|---|---|---|---|---|---|
| A*02:07 | Chu 2018 only | 1 | 2.10 | 2.1×10⁻¹² | **single_source_anchor_only** |
| B*46:01 | Chu 2018; Shin 2019; Chen 2011 | 3 | 2.04 | 7.6×10⁻³ | screening_random_effects |
| C*01:02 | Chu 2018; Shin 2019 | 2 | 1.85 | 4.6×10⁻¹⁶ | screening_random_effects |
| DPB1*05:01 | Chu 2018; Shin 2019; Chen 2011 | 3 | 2.17 | 4.1×10⁻⁹ | screening_random_effects |
| DQB1*02:01 | Chu 2018 only | 1 | 0.57 | 2.3×10⁻¹³ | **single_source_anchor_only** |
| DRB1*07:01 | Chu 2018 only | 1 | 0.43 | 2.5×10⁻²¹ | **single_source_anchor_only** |

**3/6 alleles single-source anchor only — heterogeneity test 불가.** Korean source는 어떤 allele에서도 없음.

`P4_F19_korean_replication_readiness.png` — 3-lane current → required → graduation.

---

## 4. Claim grade per allele (v2)

| Allele | Claim grade | Direction (anchor) | Direction (PTC strict) | Korean replication | Next step priority |
|---|---|---|---|---|---|
| A*02:07 | anchor_strong_no_korean_replication | GD enrichment (Chu 2018, BH-q ✓ Bonferroni ✓) | PTC enrichment (weak, n.s.) | NEEDED | high |
| B*46:01 | anchor_strong_partial_pan_asian_replication | GD enrichment (Chu+Shin+Chen meta) | PTC ~null | NEEDED | high |
| C*01:02 | anchor_strong_korean_direction_discordant | GD enrichment (Chu+Shin) | PTC depletion ⚠ | NEEDED + direction question | **high** |
| DPB1*05:01 | anchor_strong_partial_pan_asian + korean_discordant | GD enrichment (Chu+Shin+Chen) | PTC ~null | NEEDED | medium |
| DQB1*02:01 | anchor_strong_no_korean_replication | GD depletion (Chu 2018) | PTC depletion (zero-cell SUSPENDED by QC audit) | NEEDED + zero-cell audit applies cross-paper | medium-high |
| DRB1*07:01 | anchor_strong_no_korean_replication | GD depletion (Chu 2018) | PTC depletion-trend | NEEDED | high |

---

## 5. 가장 강하게 허용되는 한 문장

> "Asian Graves' disease anchor (Chu 2018, n = 1468 vs 1490)에서 6개 candidate HLA alleles 모두 BH-q AND Bonferroni < 0.05를 통과하지만, 6개 중 3개는 single-source anchor only 상태이고 Korean GD HLA matched-control NGS replication은 아직 존재하지 않으므로 Pan-Asian profile은 Korean cohort 검증 후 발표 가능한 상태가 아니다."

이게 본 자료에서 허용된 최강 표현. 이보다 강한 association/risk-allele 표현은 over-claim.

---

## 6. 절대 금지 표현

- "Korean GD HLA association"
- "DPB1*05:01 / B*46:01 / C*01:02 / A*02:07 = Korean GD risk allele"
- "DQB1*02:01 / DRB1*07:01 = Korean GD protective allele"
- "Paper 2 PTC C*01:02 depletion + Paper 4 GD C*01:02 enrichment = HT-vs-GD HLA distinction"
  — 두 paper는 서로 다른 데이터, 다른 design, 다른 질환. **Cross-paper 추론 금지.**
- "Causal HLA susceptibility for Korean GD"
- "Clinical Korean GD risk prediction"
- "Final Pan-Asian meta-analysis"

---

## 7. Korean replication plan (Paper 4 backlog gating 통과 조건)

1. **독립 Korean GD case-control NGS HLA typing** — case ≥ 200, control ≥ 200, matched.
2. Per-allele Fisher in Korean cohort + random-effects meta with ≥ 2 independent Asian GD cohorts.
3. **Allele-level harmonization** — 양측 모두 individual-level genotype data (carrier-vs-allele 메트릭 mismatch 없도록).
4. Heterogeneity test (Cochran's Q + I²); flag if I² > 50%.
5. Pre-registered hypothesis + analysis plan (random-effects + heterogeneity + correction).
6. **DQB1*02:01의 zero-cell QC audit 결과를 GD에도 cross-apply** — Korean GD에서 DQB1 callability + *02:01 vs *02:02 + rare *02:xxx 분포 audit 동일 적용.

졸업 후 표현: "Korean GD HLA pattern replicates the Pan-Asian profile" (NOT "Korean GD risk allele").

---

## 8. 첨부 자료

- `project/results/paper4_gd_hla/paper4_gd_chu2018_bh_fdr.tsv` — Chu 2018 + multi-test
- `project/results/paper4_gd_hla/paper4_gd_vs_paper2_ptc_direction.tsv` — cross-paper direction
- `project/results/paper4_gd_hla/paper4_gd_korean_replication_readiness.tsv` — readiness audit
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F17_chu2018_bh_bonferroni.png`
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F18_gd_vs_ptc_direction.png`
- `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F19_korean_replication_readiness.png`
- 웹: http://40.82.129.113/paper4-hla/ (extended with strict-metric section)

---

**End of Paper 4 GD strict-metric audit. No final association claim. Korean replication is the gate.**
