# Paper 2 Pillar II — Hashimoto-overlap PTC (HT biology) 정밀 패킷 (KR)
**Date:** 2026-05-06 · **Author:** Seungho Cook
**Scope:** Paper 2 Pillar II = Hashimoto-overlap PTC molecular biology (HT-isolated, Yu 2026-05-04 결정)
**Pillar I (HLA forest)와는 분리된 트랙입니다.** 이 패킷은 signature/biology 측면.
**Discipline:** signature-level association만. Causal HT → PTC progression 주장 금지. Paper 1 미접촉. Paper 9 quantitative connection 없음.

---

## 0. TL;DR

- **TCGA n=500** 에서 DM2 (8-gene low) tumors는 DM1 대비 HLA-II / HLA-I / T-cell / IFN-response / TLS / HT signature score 모두 elevated. **HLA-II Cohen d = +1.41, p = 8×10⁻³³.**
- TCGA HT-call rate: DM2 22.8% vs DM1 5.7% (OR = 0.205, p = 2.1×10⁻⁶) — **DM2가 DM1보다 HT-like rate 4×**.
- **6 calling rules (GMM / Otsu / top10 / top20 / top30 / resid_otsu) 모두 동일 방향**, OR 0.20-0.32, p ≤ 7×10⁻³.
- Lu 2023 single-cell GSE193581 (n=6 PTC) — per-PTC median DM ↔ HLA-II coordinated; PTC05 outlier 67% hashi-like + highest HLA-II.
- GSE286332 PTC vs PTC+HT (n=9 vs 9) — allele-level Fisher 모두 **underpowered** (p ≥ 0.20). Allele-level 결론 불가, signature-level만 사용.
- **이 결과들은 case-control 인과 association이 아니며, signature-level (TCGA snapshot) hypothesis-generating 수준**입니다.

---

## 1. 핵심 caveat

본 패킷은 **HT pathology label이 직접 측정된 case/control 비교가 아닙니다**. 모든 분석은:
- TCGA에서 "HT signature score 기반으로 HT-like 호출" → DM2 vs DM1 enrichment, OR
- Lu 2023 single-cell sample 6개 (n이 작음), OR
- GSE286332에서 라벨된 PTC vs PTC+HT (n=9 vs 9, 너무 작음)

각각 *signature-level proxy* 또는 *underpowered cohort*입니다. **Causal HT → PTC progression 주장 금지.**

---

## 2. TCGA HT signature × DM groups (n = 500)

`paper2_ht_tcga_signature_dm_groups.tsv` 참조.

| Signature | mean DM2 | mean DM1 | Cohen d (DM2 vs DM1) | MW p |
|---|---|---|---|---|
| sig_score (HT signature composite) | +0.232 | −0.597 | **+0.92** | 1.6×10⁻²⁷ |
| **HLA_I** | +0.306 | −0.787 | **+1.42** | 6×10⁻³⁴ |
| **HLA_II** | +0.302 | −0.777 | **+1.41** | 8×10⁻³³ |
| Stromal | +0.201 | −0.517 | +1.05 | 5×10⁻²³ |
| T_cell | +0.209 | −0.538 | +0.88 | 9×10⁻²³ |
| TLS | +0.150 | −0.385 | +0.72 | 8×10⁻²³ |
| IFN_resp | +0.168 | −0.433 | +0.73 | 3×10⁻¹⁸ |
| B_cell | +0.042 | −0.109 | +0.18 | 0.57 (n.s.) |

**HLA-II / HLA-I / IFN / TLS / T-cell 모두 DM2 elevated**, B-cell 단일 marker만 약함 (TLS와 다른 결과 — TLS는 multi-marker composite).

## 3. TCGA HT-call × DM crosstab (6 calling rules)

| Calling rule | DM1 HT+ % | DM2 HT+ % | OR (HT+ in DM1 vs DM2) | Fisher p |
|---|---|---|---|---|
| hashi_GMM | 5.71 | 22.78 | 0.205 | 2.1×10⁻⁶ |
| hashi_otsu | 7.14 | 24.44 | 0.238 | 4.5×10⁻⁶ |
| hashi_top10 | 4.29 | 12.22 | 0.322 | 7.4×10⁻³ |
| hashi_top20 | 7.86 | 24.72 | 0.260 | 1.0×10⁻⁵ |
| hashi_top30 | 10.71 | 37.50 | 0.200 | 6.4×10⁻¹⁰ |
| hashi_resid_otsu | 23.57 | 51.67 | 0.289 | 8.0×10⁻⁹ |

**6/6 calling rules 모두 같은 방향 (DM1 protective from HT-like).** Top-30 + resid_otsu의 OR ≈ 0.20-0.29와 가장 strong한 effect는 sample size + threshold-stringency 효과 둘 다 반영.

→ **DM2 (8-gene low / dark matter) tumors는 HT-like signature를 4–5× 더 자주 가진다** *(at signature level on TCGA snapshot)*.

## 4. Lu 2023 single-cell GSE193581 (n = 6 PTC)

`paper2_ht_lu2023_per_sample.tsv` 참조.

| Sample | n_cells | median_DM | pct_hashi_like | median_HLA_II |
|---|---|---|---|---|
| PTC06 | 417 | 0.378 | 3.8% | 0.034 |
| PTC07 | 2727 | 0.362 | 3.1% | 0.036 |
| PTC01 | 3565 | 0.279 | 1.5% | −0.217 |
| **PTC05** | **629** | **0.197** | **66.9%** | **+1.455** |
| PTC03 | 86 | 0.194 | 1.2% | −0.378 |
| PTC04 | 1197 | 0.113 | 9.7% | −0.033 |

**PTC05 outlier**: 67% hashi-like cells + median HLA-II = +1.46 (z-score). 다른 샘플은 hashi-like rate < 10%로 명확히 분리됨. Single-cell level에서도 hashi-like cells과 HLA-II expression이 coordinated.
**Caveat**: n=6 PTC samples은 generalization 근거로 부족; signature-direction 일치만 시사.

## 5. GSE286332 PTC vs PTC+HT (n=9 vs 9)

7 alleles tested allele-level Fisher — **모두 p ≥ 0.20** (underpowered). Allele-level 결론 불가.
**Useful only as**: 하시모토-overlap PTC label이 있는 RNA-seq 코호트로서 signature-level transfer 검증에 사용 (v17 D4-P2의 기반).

## 6. v17 D4-P2 signature transfer summary (참고)

GSE286332 PTC+HT signature (143 up + 46 dn genes) → TCGA에서:
- bimodality skew 1.25, kurtosis 1.628, bimodality coefficient 0.55 → 신호가 GMM/Otsu로 분리 가능
- DM2 enrichment 최대 5× (top30 OR=0.20 = 1/0.20 ≈ 5×)
- 자세한 내용: `project/results/d4p2_tcga_hashimoto_signature/D4P2_summary.json`

---

## 7. Evidence grade per track

| 트랙 | n | 등급 | 사용 가능 표현 | 사용 금지 표현 |
|---|---|---|---|---|
| TCGA HT signature × DM | 500 | **strong_exploratory** | "DM2 PTCs는 TCGA에서 HT-like signature를 더 자주 보이고 HLA-II 등 면역 axis가 elevated" | "HT가 DM2 PTC를 일으킨다" / "HT가 PTC dedifferentiation cause" |
| Lu 2023 single-cell | 6 | supportive_low_n | "single-cell pseudo-bulk에서 DM ↔ HLA-II coordinated 방향 일치" | n=6 기반 generalization |
| GSE286332 allele Fisher | 18 | **underpowered_bridge** | "signature-level cohort로만 사용 가능" | allele-level association |
| v17 D4-P2 signature transfer | 500 | strong_exploratory | "GSE286332 HT signature가 TCGA에서 bimodal하게 분포하고 DM2-enriched" | causal HT→PTC progression |

---

## 8. 가장 강하게 허용되는 한 문장

> "TCGA n=500 PTCs에서 DM2 tumors는 DM1 대비 HT-like signature, HLA-II, IFN-response, TLS 등 면역 axis가 일관되게 elevated이며 HT-call rate가 4-5× 높다 — 이는 Hashimoto-overlap PTC가 DM1과 분리된 immune-active 분자 background에 위치할 가능성을 시사하는 **signature-level exploratory 신호**이며, 이를 case-control association이나 causal HT→PTC progression으로 해석하지 않는다."

## 9. 절대 금지 표현

- "HT가 DM2 PTC를 일으킨다" / "HT causes PTC dedifferentiation"
- "Hashimoto-overlap PTC는 8-gene-low subtype이다" *(signature-level overlap이지 동치 아님)*
- "HT signature가 PTC risk biomarker"
- "Final case-control association"
- "Clinical risk prediction" / "Patient selection"
- "HT signature → 8-gene-low → Paper 9 perturbation" 같은 quantitative chain claim

---

## 10. 다음 검증 단계 (HT validation roadmap)

`HT_F4_validation_roadmap.png` 참조.

1. **독립 PTC vs PTC+HT 큰 코호트** (n ≥ 50 vs 50) with centrally reviewed pathology — TCGA + GSE286332에 의존하지 않는 데이터.
2. **HT signature 비-TCGA bulk-RNA replication** — Korean PTC bulk RNA-seq에서 HT signature 동일 방향 확인.
3. **Single-cell HT-like cluster 재현** — ≥ 2 cohorts에서 TLS / B-cell / IGHV markers 함께 검증.
4. **Pre-registered analysis plan** — HT-vs-non-HT signature comparison plan을 가설 등록 후 실행.
5. **Causal language 금지** — "Hashimoto-overlap PTC is molecularly distinct" graduation 표현까지가 한계, "HT causes PTC" 절대 사용 금지.

## 11. 졸업 기준 (graduation criteria)

- DM2 ↔ HT-signature 방향이 ≥ 2 independent cohorts에서 stable (signature level).
- Effect size > pre-specified MID in matched-cohort 비교.
- 같은 환자에 직접 HT pathology + signature 둘 다 측정 (label transfer 아닌 데이터).
- 졸업 후 표현: "Hashimoto-overlap PTC is molecularly distinct" (NOT "HT causes PTC").

---

## 12. 첨부 자료

- `project/results/paper2_ht_biology/paper2_ht_tcga_signature_dm_groups.tsv` — DM2 vs DM1 signature audit
- `project/results/paper2_ht_biology/paper2_ht_lu2023_per_sample.tsv` — Lu 2023 per-PTC summary
- `project/results/paper2_ht_biology/paper2_ht_evidence_grade.tsv` — track-level grade
- `project/papers_hub_2026_05_04/assets/paper2_ht/HT_F1_tcga_signature_forest.png` — Cohen d forest
- `project/papers_hub_2026_05_04/assets/paper2_ht/HT_F2_tcga_hashi_dm_crosstab.png` — 3-panel crosstab
- `project/papers_hub_2026_05_04/assets/paper2_ht/HT_F3_lu2023_single_cell_scatter.png` — Lu 2023 scatter
- `project/papers_hub_2026_05_04/assets/paper2_ht/HT_F4_validation_roadmap.png` — validation roadmap
- 웹: http://40.82.129.113/paper2-ht/

---

**End of HT biology packet. Signature-level exploratory only — case-control association 결과 아님.**
