# 2026-05-06 KEY FINDINGS — 한 화면 인덱스
**작성:** Seungho Cook · **모드:** professor-facing · **원칙:** 모두 exploratory / strengthening level — 어떠한 final association도 주장하지 않음

---

## 0. "지금 뭐 봐야 되나" — 가장 빠른 경로

| 우선 | URL | 무엇 |
|---|---|---|
| 1 | http://40.82.129.113:8012/papers_hub_2026_05_04/index.html | 8-papers 허브, top에 2026-05-06 strengthening banner + Paper 2/4 카드에 새 버튼 |
| 2 | http://40.82.129.113/paper2-hla/ | **Paper 2 Pillar I — HLA forest v2 (DQB1*02:01 QC 감사 + claim_grade v2)** |
| 3 | http://40.82.129.113/paper2-ht/ | **Paper 2 Pillar II — HT biology (NEW)** |
| 4 | http://40.82.129.113/paper4-hla/#strict-metric | **Paper 4 GD strict-metric strengthening (anchor at #strict-metric)** |
| 5 | http://40.82.129.113/gpl570_validation_pack.html | Paper 1 GPL570 external validation (이전 세션) |

웹은 port 80(`/var/www/papers/`)과 port 8012(serve_secure, project root) 둘 다 200 OK.

---

## 1. Paper 1 — DM1 dark matter (driver-orthogonal differentiation axis)

**상태:** SHIPPED. GPL570 external validation은 마감.

**Key numbers:**
- 4 GPL570 cohorts (GSE76039 anchor + GSE33630 + GSE29265 + GSE65144) — 총 216 thyroid samples.
- DM1_like vs THYROID_NONOVERLAP Spearman ρ ∈ [−0.846, −0.941], anchor −0.925.
- 30/30 lineage axis × dataset × ATC contrast cells direction-consistent.
- TACSTD2: 1/5 cells consistent → supplement only (정직한 negative).
- GSE33630 ATC vs Normal: THYROID_NONOVERLAP **d = −6.53**, DM1_like d = +5.43.

**Paper 1 자료:**
- 분석 스크립트: `project/notebooks_or_scripts/p_gpl570_validation_first_pass.py`
- 출력: `project/results/p_gpl570_validation/`
- 보고서: `project/reports/2026_05_04_gpl570_external_validation_report.md`
- Self-contained 웹: `project/papers_hub_2026_05_04/gpl570_validation_pack.html`

---

## 2. Paper 2 Pillar I — HLA forest (v2) ⚠ DQB1*02:01 SUSPENDED

**페이지:** http://40.82.129.113/paper2-hla/

**v2 core statement (페이지 상단 빨간 박스):**
> S0 carrier-vs-allele forest는 6/6 significant 신호를 보였지만, **strict S3 allele-vs-allele harmonization 후에는 C*01:02, DQB1*02:01 두 개 depletion 후보만 남음. Enrichment 후보는 0.**

### 2.1 v2 claim grade table (per allele)

| Allele | Claim grade (v2) | S0 q | S3 direct q | S0 dir | S3 dir | Priority |
|---|---|---|---|---|---|---|
| A*02:07 | **S0_only_exploratory** | 1.1×10⁻⁵ | 0.43 (n.s.) | enrich | enrich (OR=1.26) | medium |
| B*46:01 | **S0_only_exploratory** | 1.1×10⁻⁵ | 0.80 (n.s.) | enrich | enrich (OR=1.05) | medium |
| C*01:02 | **survives_strict_metric_but_direction_flips** | 4.3×10⁻⁴ | **7.4×10⁻⁴ ✓** | enrich | **depletion** (OR=0.68) | **high** |
| DPB1*05:01 | **source_sensitive_artifact_likely** | 9.9×10⁻¹⁴ | 0.80 (n.s.) | enrich | depletion (~null) | medium |
| DQB1*02:01 | **survives_strict + zero_cell_caution** ⚠ SUSPENDED by QC audit | 2.5×10⁻⁶ | **5.3×10⁻⁸** | depletion | depletion (zero-cell) | **high** |
| DRB1*07:01 | **S0_only_exploratory** | 4.3×10⁻⁴ | 0.43 (n.s.) | enrich | depletion (OR=0.84) | medium |

### 2.2 DQB1*02:01 zero-cell QC audit (2026-05-06)

**Verdict: 기술적 의심 — arcasHLA subtype-resolution artifact 가능성 높음. 일시 SUSPEND 상태.**

| 검사 | 결과 | 통과? |
|---|---|---|
| DQB1 locus callability ≥ 95% | **72.2%** (631/874) — A/B/C는 ≥99.9% | ❌ |
| DQB1*02 family rate (Korean baseline ~13–17%) | **13.3%** | ✅ |
| 표준 DQB1*02:01이 DQB1*02 family 내 dominant | **0 obs** (DQB1*02:02=53, *02:276=14, *02:253=10, *02:225Q=4, *02:199=3, *02:214N=1, …) | ❌ |
| Rare-subtype tail < 5% | **11.1%** (subtype ≥ 50) | ❌ |
| 3 subcohorts 일관성 | K2 0/163, Lee2024 0/464, GSE286332 0/4 | systematic |

핵심: DQB1*02 family는 Korean baseline rate대로 나타나는데, 표준 *02:01 자리만 비어있고 mass가 rare subtype tail로 분산됨. probability of zero observed *02:01 if true AF=2.1% is ≈ 10⁻¹¹ — biologically near-impossible without typing artifact.

**Files:**
- `project/notebooks_or_scripts/paper2_hla_professor_strengthening.py`
- `project/notebooks_or_scripts/paper2_dqb1_zero_cell_qc_audit.py`
- `project/results/p2_pillar1_forest_v2/paper2_hla_claim_grade_v2.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_hla_metric_sensitivity_strengthened.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_hla_allele_vs_allele_direct.tsv`
- `project/results/p2_pillar1_forest_v2/paper2_dqb1_zero_cell_qc_audit.tsv`
- `project/papers_hub_2026_05_04/assets/paper2_hla/F16, F17, F18, F19_validation_roadmap, F20_DQB1_zero_cell_qc_matrix, F21_DQB1_locus_callability_barplot`
- 보고서: `2026_05_06_paper2_hla_professor_packet_kr.md` · `2026_05_06_paper2_hla_one_page_summary_kr.md` · `2026_05_06_paper2_hla_claim_boundary.md` · `2026_05_06_paper2_hla_next_validation_plan.md` · `2026_05_06_paper2_dqb1_zero_cell_qc_report.md`

---

## 3. Paper 2 Pillar II — HT biology (NEW page) ★

**페이지:** http://40.82.129.113/paper2-ht/

**가장 강하게 허용되는 한 문장:**
> TCGA n=500 PTCs에서 DM2 tumors는 DM1 대비 HT-like signature, HLA-II, IFN-response, TLS 등 면역 axis가 일관되게 elevated이며 HT-call rate가 4–5× 높다 — Hashimoto-overlap PTC가 DM1과 분리된 immune-active 분자 background에 위치할 가능성을 시사하는 **signature-level exploratory 신호**이며, case-control association이나 causal HT→PTC progression이 아니다.

### 3.1 TCGA HT signature × DM (n=500)

| Signature | Cohen d (DM2 − DM1) | MW p |
|---|---|---|
| HLA_II | **+1.41** | 8×10⁻³³ |
| HLA_I | **+1.42** | 6×10⁻³⁴ |
| sig_score (HT signature 종합) | +0.92 | 2×10⁻²⁷ |
| Stromal | +1.05 | 5×10⁻²³ |
| T_cell | +0.88 | 9×10⁻²³ |
| TLS | +0.72 | 8×10⁻²³ |
| IFN_resp | +0.73 | 3×10⁻¹⁸ |
| B_cell (단일 marker) | +0.18 | 0.57 (n.s.) |

### 3.2 TCGA HT-call × DM crosstabs (6 calling rules)

DM2 HT+ rate가 DM1 대비 **4–5×** (top30: 37.5% vs 10.7%, OR=0.20, p=6×10⁻¹⁰).
6/6 calling rules 동일 방향.

### 3.3 Lu 2023 single-cell GSE193581 (n=6 PTC)

PTC05 outlier: 67% hashi-like cells + median HLA-II = +1.46. 나머지 5개 샘플 hashi-like rate < 10%로 명확히 분리.

### 3.4 GSE286332 PTC vs PTC+HT (n=9 vs 9)
모든 7개 alleles Fisher p ≥ 0.20 — **underpowered**. Allele-level 결론 불가, signature-level cohort로만 사용.

**Files:**
- 분석 스크립트: `project/notebooks_or_scripts/paper2_ht_biology_strengthening.py`
- 출력 TSV: `project/results/paper2_ht_biology/`
- 그림: `project/papers_hub_2026_05_04/assets/paper2_ht/HT_F1, HT_F2, HT_F3, HT_F4`
- KR 패킷: `project/reports/2026_05_06_paper2_ht_biology_packet_kr.md`

---

## 4. Paper 4 — Korean GD HLA (strict-metric strengthening) ★

**페이지:** http://40.82.129.113/paper4-hla/#strict-metric (앵커 직접)

**가장 강하게 허용되는 한 문장:**
> Asian Graves' anchor (Chu 2018, n=1468 vs 1490)에서 6개 후보 HLA alleles 모두 BH-q AND Bonferroni < 0.05를 통과하지만, **6개 중 3개는 single-source anchor only이고 Korean GD HLA matched-control NGS replication은 아직 존재하지 않으므로 Pan-Asian profile은 Korean cohort 검증 후 발표 가능한 상태가 아니다.**

### 4.1 Chu 2018 anchor — multi-test 통과 (6/6 ✓)

| Allele | OR | raw p | BH-q | Bonferroni p | 방향 |
|---|---|---|---|---|---|
| A*02:07 | 2.10 | 2.07×10⁻¹² | 2.07×10⁻¹² ✓ | 1.24×10⁻¹¹ ✓ | GD-enriched |
| B*46:01 | 2.38 | 8.78×10⁻²¹ | 1.76×10⁻²⁰ ✓ | 5.27×10⁻²⁰ ✓ | GD-enriched |
| C*01:02 | 1.83 | 5.40×10⁻¹⁵ | 8.10×10⁻¹⁵ ✓ | 3.24×10⁻¹⁴ ✓ | GD-enriched |
| DPB1*05:01 | 1.90 | 1.73×10⁻²⁶ | 1.04×10⁻²⁵ ✓ | 1.04×10⁻²⁵ ✓ | ★ Asian Graves' top |
| DQB1*02:01 | 0.57 | 2.31×10⁻¹³ | 2.77×10⁻¹³ ✓ | 1.39×10⁻¹² ✓ | GD-depleted |
| DRB1*07:01 | 0.43 | 2.49×10⁻²¹ | 7.47×10⁻²¹ ✓ | 1.49×10⁻²⁰ ✓ | GD-depleted |

### 4.2 Cross-paper direction (GD vs PTC S3, 같은 6 alleles)

| Allele | GD direction | PTC S3 direction | 일치? |
|---|---|---|---|
| A*02:07 | enrichment | enrichment (n.s.) | **BOTH enrichment** ✓ |
| B*46:01 | enrichment | enrichment (n.s.) | **BOTH enrichment** ✓ |
| C*01:02 | enrichment | depletion (q=7e-4) | **OPPOSITE** ⚠ |
| DPB1*05:01 | enrichment | depletion (~null) | **OPPOSITE** ⚠ |
| DQB1*02:01 | depletion | depletion (zero-cell SUSPENDED) | **BOTH depletion** |
| DRB1*07:01 | depletion | depletion-trend (n.s.) | **BOTH depletion** ✓ |

**Yu 2026-05-04 분리벽 enforce:** Paper 2/Paper 4는 서로 다른 질환·design — cross-paper 인과 추론 절대 금지.

### 4.3 Pan-Asian meta + Korean replication readiness

3/6 alleles **single_source_anchor_only** (A*02:07, DQB1*02:01, DRB1*07:01) — Chu 2018만 있음.
3/6 alleles screening_random_effects with Chu+Shin±Chen (B*46:01 k=3, C*01:02 k=2, DPB1*05:01 k=3) — 그러나 Korean source 0.

**Korean GD case-control NGS HLA cohort = NOT YET ASSEMBLED → 4/4 backlog gating 미통과.**

**Files:**
- 분석 스크립트: `project/notebooks_or_scripts/paper4_gd_strict_metric_strengthening.py`
- 출력 TSV: `project/results/paper4_gd_hla/paper4_gd_chu2018_bh_fdr.tsv`, `paper4_gd_vs_paper2_ptc_direction.tsv`, `paper4_gd_korean_replication_readiness.tsv`
- 그림: `project/papers_hub_2026_05_04/assets/paper4_hla/P4_F17, P4_F18, P4_F19`
- KR 감사: `project/reports/2026_05_06_paper4_gd_strict_metric_audit_kr.md`

---

## 5. 절대 금지 (전체)

- "Final case-control association"
- "Causal HLA susceptibility" / "Korean PTC HLA risk allele" / "Korean GD risk allele"
- "Clinical risk prediction" / "Patient selection"
- "Carrier 빈도와 allele 빈도 동등 비교 가능"
- "DPB1*05:01 / B*46:01 / A*02:07 / C*01:02 / DRB1*07:01 = Korean PTC OR Korean GD risk allele"
- "Paper 2 PTC C*01:02 depletion + Paper 4 GD C*01:02 enrichment = HT-vs-GD HLA distinction" (cross-paper inference 금지)
- "HT가 PTC를 일으킨다" / "HT causes PTC dedifferentiation"
- "HT signature → 8-gene-low → Paper 9 perturbation chain" (quantitative connection 금지)

---

## 6. 진척 및 다음 단계 사이드맵

| Track | 단계 | 다음 액션 |
|---|---|---|
| Paper 1 | shipped | manuscript 본문 (voice-protected, user keyboard) |
| Paper 2 Pillar I | strengthening v2 lock | matched-control NGS + DQB1 orthogonal NGS-typer |
| Paper 2 Pillar II HT | NEW first pass | independent PTC vs PTC+HT (n ≥ 50 vs 50) + replicated signature |
| Paper 3 ICI | bundle frozen | Track A 유지, Track B = Paper 1 bioRxiv 후 |
| Paper 4 GD | backlog 4/4 | **Korean GD case-control NGS HLA pilot** |
| Paper 9 perturbation | first pass + sprint1+2 | additional 5-hour sprint (예정) |

---

## 7. 본 세션의 수확 — 한 줄

> Paper 2 HLA 패킷에 metric-mismatch strengthening, DQB1 QC suspension, claim_grade v2가 들어가 over-claim 가능성이 봉쇄되었고, 별도로 Paper 2 Pillar II (HT biology, NEW)와 Paper 4 GD strict-metric strengthening이 각각 separate 페이지로 발행되어, 8-papers 허브 (포트 8012/80 둘 다)에서 한 클릭으로 도달 가능합니다. **모두 exploratory, final association 없음, Paper 1 미접촉.**

---

**End of key-findings index. 모든 commit은 보류 상태 (cleanup plan G3 voice-protected 제외; 새 추가분은 G4 확장 또는 별도 그룹으로 묶을 수 있음 — 지시 시 진행).**
