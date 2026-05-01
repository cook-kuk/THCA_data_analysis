---
title: "오후 미팅 브리프 — 2026-04-29 (유형원 교수님)"
purpose: 오전 미팅에서 제기된 3개 의문에 대한 forensic 결과를 1-pager로 정리. 미팅 들어가기 전 5분 안에 읽도록.
triage: A=blocker check, B=stat sanity check, D=cohort identity confusion
---

# 핵심 메시지 3줄

1. **8-gene 셀렉션은 RAI bias 의심받을 수 있지만 "의도적 design choice"였다.** Driver 유전자(BRAF/TERT)를 candidate pool에서 명시적으로 제외하는 코드(`TIERA67_CLEAN_CATEGORIES`, `rerun_v2.py:198`)가 있고, manuscript title 자체가 "**RAI-responsiveness biomarker**"로 framing되어 있음. **Paper는 차단되지 않음. Methods 한 문장만 reframe.**

2. **TERT-only triple-negative-ish가 worst인 결과는 group definition artifact다.** 4-way category에서 TERT+ (n=36)이 BRAF/RAS와 mutually-exclusive하게 정의되어 있는데, 실제로는 36명 중 25명(69%)이 BRAF+. "TERT-only로 진짜 trip-neg-otherwise"인 환자는 ~5-10명 추정 → small N에 의한 outlier-driven HR. **Firth correction + 8-cell breakdown으로 재계산 필요.**

3. **🚨 K2 ≠ 분당 SNUH.** 미팅에서 보신 Korean validation set은 **K2 = PRJEB11591 = SNU-GMI Yoo 2016 공개 cohort** (n=260). 진짜 분당 SNUH cohort는 outreach email 단계 (2026-04-27 작성), **데이터 수령 안 됨**. 미팅에서 이 점 명확히 짚어야 misunderstanding 안 생김.

---

# A. 8-gene 포렌식 결과 (BLOCKER 해제)

**Selection mechanism (실제):**
- Pool: 67-gene 큐레이션 (TIERA67) — 7 카테고리 (TDS_core, MAPK_output, Driver_anchor, Aggressive_marker, Dediff_invasion, Immune_stromal, Thyroid_lineage_extra)
- **Driver_anchor 제외** (`{k: v for k, v in TIERA67_CATEGORIES.items() if k != "Driver_anchor"}`) → 55-gene clean pool
- RandomForest feature importance로 ranking
- Top-8: NIS/SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1

**왜 BRAF/TERT가 안 나왔나:** Driver_anchor 카테고리에 있어서 candidate pool에서 명시적으로 빠짐. **의도적**. label leakage 방지 (DM1/DM2 cluster가 BRAF-like/RAS-like subtype과 confound되는 거 막으려고).

**Manuscript v6 disclosure 상태:**
- Title: "An **8-gene RAI-responsiveness biomarker** outperforms BRAF V600E status..."
- Abstract: "DM1/DM2 axis statistically distinct from BRAF/RAS dichotomy (ρ=0.49)"
- 즉 RAI framing은 이미 paper의 핵심 narrative

**Leak-free 재검증 이미 했음 (R1-B):**
- 30-gene MAPK+immune+EMT panel (8-gene과 zero overlap)로 cluster 재학습
- 8-gene이 새 cluster를 AUC 0.925로 예측 (ΔAUC +0.130 vs BRAF)
- 즉 8-gene이 측정하는 "differentiation axis"는 panel 외부에서도 reproducible

**미팅 발언 권고:**
> "8-gene이 RAI 쪽으로 몰린 건 의도된 design choice였습니다. Driver gene을 candidate pool에서 명시적으로 빼는 코드가 있고요, 그렇게 한 이유는 BRAF-like/RAS-like label과 confound 막기 위해서였습니다. 그래서 결과적으로 differentiation axis 위에 sit하는 panel이 나왔고, 이건 BRAF mutation status랑 orthogonal한 axis (Spearman ρ=0.49)로 paper에서 framing되어 있습니다. Methods 한 문장만 더 명확하게 쓰면 reviewer-proof 합니다."

상세: `audit_report_8gene.md` (smoking-gun line numbers + Q&A)

---

# B. TERT-only paradox 해소

**TCGA-THCA n=504 4-group breakdown (existing data, `FINAL_crosstab_4group.tsv`):**

| Group | N | OS events | Event rate |
|-------|---|-----------|------------|
| BRAF only (TERT-) | 250 | 3 | 1.2% |
| RAS only (TERT-) | 48 | 0 | 0.0% |
| TERT+ (any BRAF/RAS) | 36 | 6 | **16.7%** |
| Triple-neg | 170 | 7 | 4.1% |

**숨은 문제:** "TERT+" 36명 안에서:
- BRAF+ & TERT+: ~25명 (69%) ← 진짜로 worst
- RAS+ & TERT+: 소수
- 트리플 negative-otherwise & TERT+: ~5-10명 ← 미팅에서 "이상하게 worst"라고 본 그룹

**원인:** 4-group 분류가 TERT를 첫 우선순위로 두고 있어서, TERT+ 안의 hetero를 가려버림. 8-cell (2^3) breakdown 안 하면 small-N artifact를 진짜 biology처럼 오해.

**미팅 발언 권고:**
> "TERT-only triple-negative-ish가 worst로 나온 건 group definition 때문이고, 그 cell의 N이 5-10명대로 작아서 HR이 outlier-driven인 것 같습니다. 8-cell (BRAF±, RAS±, TERT±)로 다시 breakdown하면 literature 일관 (BRAF+TERT+ worst, HR 8.51 from Xing 2014) 회복될 것 같고, Firth correction 적용해서 supplementary fig로 자세히 정리하겠습니다."

상세 데이터: `project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv`

---

# D. 🚨 K2 ≠ 분당 SNUH — 미팅에서 반드시 짚을 것

**현재 분석에 들어간 Korean validation set의 정체:**
- 이름: K2 / KOREAN_K2
- 실체: **PRJEB11591 = Yoo et al. 2016 (SNU-GMI) RNA-seq cohort, 공개 데이터 (EBI ENA)**
- N: 260 samples (cPTC 77, FVPTC 48, FTC 30, FA 23, normal 81, unknown 1)
- Modality: paired-end RNA-seq (HiSeq 2000), kallisto pseudoalignment 8-gene mini-index
- 파일: `project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv`

**진짜 분당 SNUH cohort 상태:**
- Outreach email v2 작성: 2026-04-27 (`email_bundang_KR_v2.md`, n≈100 RNA-seq 제안)
- **데이터 수령 안 됨, 협업 합의 안 됨**
- v12_snubh_shortlist.py 는 prospective gene candidate list (validation 후보), 실제 data import 아님

**미팅 발언 권고:**
> "교수님께서 보신 Korean validation 그래프는 사실 K2 cohort, 즉 Yoo Seunggeun 박사님의 2016 SNU-GMI 공개 데이터입니다. 분당 SNUH는 4월 27일에 outreach email v2 보낼 준비 단계고, 아직 데이터 수령 전입니다. 분당 cohort 협업이 시작되면 그때 진짜 'Bundang validation'이 됩니다. 지금 plot은 '공개 Korean cohort에서의 외부 validation'으로 caption 정정하겠습니다."

이거 안 짚으면 manuscript reviewer 또는 collaborators가 분당 데이터로 오해할 위험 있음.

---

# 우선순위 (오후 미팅 후)

| 시간 | 작업 |
|------|------|
| 오늘 밤 | A의 manuscript edit suggestion 적용 / B의 8-cell breakdown + Firth Cox 재돌림 |
| 내일 | C robustness benchmark (overnight job) / F sc wrap-up figure 5 / H 262 FFPE QC |
| 다음 주 | G trajectory / E MSK caveat / J Wang 정체확인 |
| 별도 트랙 | I NRG1 separate paper / Graves pivot (별도 prompt 파일) |

---

# 미팅 들어가기 직전 체크

- [ ] audit_report_8gene.md 한 번 훑어 (5분 안에 됨)
- [ ] "K2 = SNU-GMI ≠ 분당" 입에 익히기
- [ ] 4-way matrix는 "small-N artifact 의심, Firth로 재검증 중" 한 줄
- [ ] Graves pivot 결정 미팅에서 다룰지 말지 사전 정리 (별도 prompt 파일 참조)

— end —
