# 유 교수님께 — v17 npj 투고 전 마지막 리뷰 (BOOST sprint 후)

날짜: 2026-04-27
분량: A4 1–2 페이지

## 1. 한 페이지 요약

- **헤드라인 (변하지 않음):** 8-gene RAI panel이 BRAF V600E 단독 baseline 보다 ΔAUC = **+0.132** 만큼 우수하다 (CV AUC 0.954 vs 0.822, n=513 TCGA-THCA).
- **이번 BOOST sprint 에서 추가된 3 개 robustness 레이어 (reviewer 가 반드시 묻는 것):**
  1. **Multi-cohort meta-analysis** (4개 외부 cohort, n=290) → pooled AUC = **0.980 (95% CI 0.869–0.997, I² = 0%)**, 4 cohort 모두 AUC > 0.96.
  2. **Decision Curve Analysis** (Vickers–Elkin 2006) → 8-gene 전략이 threshold 0.05–0.95 의 모든 91 지점에서 dominant. ("어떤 임계값을 쓰든 8-gene 이 우월하다" — 가장 강한 형태.)
  3. **9-strata subgroup forest** → 9 개 중 **7 개**에서 AUC ≥ 0.85; Stage III/IV (가장 위험군이자 임상 결정이 가장 어려운 군) 에서 **AUC 0.996**.
- **인터랙티브 RAI calculator** (HTML 1 페이지, 서버 없음) — 비전산 reviewer 도 직접 8-gene log2(TPM+1) 값을 넣어보고 logit P(DM2) score 를 실시간으로 받을 수 있음.
- **TERT (R8) 는 honest 하게 정정** — univariate Cox HR 6.31 → multivariate (stage + age + sex 보정) HR 1.88, p = 0.29. "stage-and-age-independent prognostic marker" 가 아니라 "Stage III/IV 식별을 위한 molecular handle" 로 reframe. 1000-iter bootstrap median p = 2.6×10⁻⁵, leave-one-out worst-case p < 10⁻³ 은 그대로 유지.

## 2. 진짜 기여와 솔직한 약점

### 진짜 기여 (npj scope 와 정확히 맞음)
1. 임상에서 바로 쓸 수 있는 양적 biomarker (BRAF V600E 단독을 능가하는 것을 정량화).
2. Decision Curve 가 모든 임계값에서 dominant — 인구단위 임상 유용성이 검증됨.
3. 4 cohort meta-analysis 의 I² = 0% — 갑상선 transcriptomic biomarker 에서 매우 드문 tight 한 결과.
4. DM1/DM2 axis 가 BRAF/RAS dichotomy 와 statistically distinct (Spearman ρ = 0.49) 하면서 immune hot/cold landscape 에 mapping (Cohen's d = +1.683).
5. 두 개의 ongoing thyroid TROP2-ADC trial (NCT06235216, NCT07521670) 에 correlative biomarker overlay 로 즉시 제안 가능.

### 솔직한 약점 (limitations 로 모두 명시함, 13 개)
1. Wet-lab validation 없음.
2. Korean cohort 아직 통합 안 됨 (현재 SNUH / 분당서울대 outreach 진행 중) — meta-analytic I² = 0% 가 강화는 하지만 대체하지는 못함.
3. DM1/DM2 cluster 자체는 OS independent 가 아님 (TCGA-THCA 의 exceptional prognosis 한계, ~16 OS events).
4. PRISM screen n=5 vs n=5 underpowered.
5. scRNA n=7 patients only.
6. 4 meta cohort 중 GSE76039 는 trajectory analysis 에 사용된 cohort — meta 결과는 strict held-out 가 아니라 transferability check.
7. DCA dominance 는 population-level — individual-decision calibration (intercept α=0.04, slope β=0.91) 은 다음 단계.
8. TERT calls 가 BAM re-call 이 아니라 cBioPortal mirror.
9. TERT subset 6 events 만 있음 (그래서 bootstrap + LOO 다 보고함).
10. Univariate HR 6.31 → multivariate HR 1.88 의 stage adjustment caveat.
11. Pan-cancer transfer 는 signature overlap 만 보았음.
12. 8-gene panel 이 RAI uptake score 의 함수이므로 직접 self-prediction 은 informative 하지 않음 — DM1/DM2 가 proxy target.
13. 외부 cohort label 은 unsupervised stratification proxy (직접적 RAI response truth 가 아님).

## 3. 교수님께 여쭤보고 싶은 5 가지 결정 사항

| # | 결정 항목 | 현 상태 / 후보 | 교수님 의견 부탁드립니다 |
|---|-----------|--------------|------------------------|
| 1 | **저자 표기 (Yu Kyungho 풀네임)** | manuscript / cover letter 의 placeholder `[Yu Kyungho]` — 한글 + 영문 표기 결정 필요 | 사용하실 영문 spelling (예: Kyungho Yu vs Kyung-Ho Yu) |
| 2 | **소속 affiliation** | placeholder `[Department TBD], [Affiliation TBD], Seoul, Republic of Korea` | 현 소속 (분당서울대, 서울대 의대 등) 정확히 |
| 3 | **이메일** | placeholder `[email TBD]` | corresponding 으로 노출하실 이메일 |
| 4 | **분당서울대 cohort 협업 timing** | 현재 outreach DRAFT (`outreach/email_*_DRAFT.md`) 작성만 됨, 미발송 | 투고 전 발송할지, revision round 에서 발송할지 — 그리고 첫 contact 받을 분 |
| 5 | **인용 추가 (한국 논문)** | 현재 ref [12] Yoo SK et al. PLoS Genet 2017 만 한국 group | 추가로 인용해야 할 한국 그룹 논문 (특히 갑상선 분자 분류 / TERT / RAI uptake) |

## 4. 본 reviewer 입장에서 가장 위험한 attack 5개와 우리의 답변

| Attack | Reviewer 가 묻는 이유 | 우리 답변 |
|--------|---------------------|---------|
| "TCGA-THCA 한 cohort 에서만 train" | overfitting 의심 | 4 cohort meta I² = 0%, pooled AUC 0.980 |
| "AUC 0.954 가 인상적이긴 한데 임상적으로 의미 있나?" | discrimination ≠ utility | DCA 가 모든 91 thresholds 에서 dominant — 인구단위 utility 입증됨 |
| "FVPTC 같은 minority subgroup 에서는 모델이 무너지는 것 아닌가?" | 형평성 / robustness | 9-strata forest, 7/9 ≥ 0.85; FVPTC 0.823 (CI overlaps 0.85), 가장 위험한 Stage III/IV 는 0.996 |
| "TERT HR 6.31 은 stage 보정하면 사라지지 않나?" | 정확하게 — 그래서 우리가 먼저 보고함 | 보정 후 HR 1.88, p = 0.29. "Stage III/IV molecular handle" 로 reframe. Bootstrap + LOO 는 그대로 유지 |
| "비전산 reviewer 가 어떻게 모델을 직접 확인하나?" | 재현성 / 투명성 | 인터랙티브 calculator (HTML 1 페이지) 동봉 — 서버 없이 8 gene log2(TPM+1) 입력해서 score 직접 확인 |

## 5. 교수님께서 검토하실 파일 (총 4 개)

1. `submission/npj/manuscript_v4.pdf` — 본문 (BOOST sprint 반영된 최종)
2. `submission/npj/cover_letter_v3.pdf` — 커버 레터
3. `reports/v17_boost/index.html` — 한국어 검토 대시보드 (이 문서의 인터랙티브 버전)
4. `reports/v17_boost/rai_calculator.html` — 인터랙티브 RAI score calculator

(서버 주소: 40.82.129.113:8765 — 위 4 개 파일이 모두 그 서버에서 직접 열림)

## 6. 다음 액션

- 교수님 의견 수렴 → manuscript / cover letter 의 placeholder 5 개 fill in.
- reviewer_bundle.zip 를 v4 로 regenerate.
- npj Editorial Manager 에 업로드 (SUBMIT_INSTRUCTIONS.md 참조).
- 분당서울대 cohort outreach: revision round 에 보낼 수 있도록 draft 보강.
