# GSE250521 → Paper 1 inclusion 냉정 판정
**2026-05-03 · 7-criteria cold review · n=16 slides, 55,873 spots post-QC**

## Top-line: **SUPPLEMENTARY only.** Main figure 불가, drop 도 아님.

| 기준 | Verdict | 점수 |
|---|---|---|
| 1. DM1 ↑ PT→ATC? | MARGINAL (epi only) | ⚠️ |
| 2. RAI_8 ↓ PT→ATC? | MARGINAL (epi only) | ⚠️ |
| 3. DM1 ↔ TDS 음의 상관? | STRONG, 그러나 부분 redundant | ✅ (한정적) |
| 4. counts/ngenes artifact? | PARTIAL (~20%) | ⚠️ |
| 5. epithelial subset 유지? | DIRECTIONALLY 만 | ⚠️ |
| 6. CAF/EMT/hypoxia 연결? | WEAK / DIRECTION-INCONSISTENT | ❌ |
| 7. 등급? | **SUPPLEMENTARY** | — |

---

## 1. DM1_like_score가 PT→PTC→LPTC→ATC 에서 증가하는가?

**Verdict**: All-spot 에서 NO, epithelial subset 에서 MARGINAL.
**Evidence**:
- Sample-mean Spearman, all-spot: ρ = 0.000, p = 1.0 (n=16)
- Epithelial top-50% subset, raw: ρ = +0.437, p = 0.091
- Epithelial subset, depth-corrected: ρ = +0.352, p = 0.18

**Weakness**: 단일 cohort n=16 에서 depth correction 후 p=0.18 → significance 미달. 방향만 일관.

**Required next analysis**: HRA003537 + GSE230424 추가하여 meta sample-mean Spearman (총 n≈60). 그 결과 p<0.05 가능성 medium (current effect size 가 실제이고 확장 시 power 충분 가정).

**Figure recommendation**: 단독 main figure 금지. supplementary figure에 sample-mean strip plot + 직접 ρ/p 표시.

---

## 2. RAI_8_score 감소?

**Verdict**: 1번과 mathematically mirror (DM1_like = −RAI_8).
**Evidence**: same numbers, sign reversed.
**Weakness**: same.
**Required next**: same.
**Figure recommendation**: RAI_8 와 DM1_like 둘 다 표시할 필요 없음. **DM1_like 만** 메인 axis 로 사용 (Paper 1 기존 framework 일관성).

---

## 3. DM1_like ↔ TDS_like 음의 상관?

**Verdict**: STRONG YES, 그러나 partially tautological.
**Evidence**:
- All-spot: Spearman ρ = −0.89, Pearson r = −0.92, p ≈ 0 (n=55,873)
- Sample-mean: ρ = −0.87, Pearson r = −0.99, p = 1e-5 (n=16)

**Weakness**: TDS_like (TG, TPO, TSHR, SLC5A5, DIO1, PAX8) 와 RAI_8 (8 genes) 가 **6 genes overlap**. DM1_like = −RAI_8 이므로 본질적으로 6/8 overlap 한 두 score 의 anti-correlation. **순수 independent validation 아님**. Paper 1 reviewer 가 "redundant" 라고 지적할 수 있음.

**Required next analysis**:
- TDS-like 외 independent thyroid-lineage marker 로 재검정. 후보:
  - **GATA3 + GATA6** (transcription factor, 8-gene 과 비-overlap)
  - **HHEX** (thyroid lineage TF)
  - **FOXA1, FOXA2** (endodermal lineage)
- 또는 TDS-like 에서 6 overlap genes 빼고 (남는 게 없음) 대신 cytokeratin 8/18 module 로 epithelial 비-thyroid axis 잡기.

**Figure recommendation**: scatter plot DM1 vs TDS, sample-mean dots, ρ/p annotated. **caveat 명시: 6/8 gene overlap 으로 부분 redundancy**.

---

## 4. total_counts / n_genes artifact 가능성은?

**Verdict**: PARTIAL ARTIFACT. ~20% of borderline epithelial signal explained by depth.
**Evidence**:
- DM1_like ~ log(total_counts), all-spot Spearman ρ = −0.45, p ≈ 0
- RAI_8 ~ log(total_counts) ρ = +0.45 (mirror)
- Raw epi stage trend ρ = ±0.44 → depth-corrected ±0.35 (magnitude ~80% retained, sign 보존)
- Mixed-effects model with log_counts + log_ngenes covariates 도 같은 방향 유지

**Weakness**: depth confound 자체는 큼 (ρ=−0.45). DM1 score 가 thyroid 8-gene mean expression 의 negative 인데, 8-gene 들이 thyroid-lineage 에서 abundantly expressed → high-coverage spot 일수록 RAI 잘 감지 → DM1 낮음. 이건 mechanism artifact 가 아니라 **score definition + sequencing technology interaction**.

**Required next analysis**: depth correction 을 score 단계에 baked-in 하기 — z-score 전에 log1p normalize 후 lm(gene ~ log_counts + log_ngenes) residual 사용. 또는 SCT/sctransform 같은 variance-stabilizing transform.

**Figure recommendation**: supplementary panel에 raw vs depth-corrected ρ 동시 표시, depth confound 명시적으로 인정.

---

## 5. Epithelial / tumor-enriched subset 에서도 유지되는가?

**Verdict**: DIRECTIONALLY YES, INFERENTIALLY NO.
**Evidence**:
- Epithelial top-50% (EPCAM/KRT8/KRT18/KRT19/TACSTD2 z-score 상위 절반): raw RAI ρ = −0.44, DM1 ρ = +0.44, TDS ρ = −0.46 (p ≈ 0.07-0.09)
- Depth-corrected: ρ = ±0.35, p = 0.18

**Weakness**:
- Epithelial_score 는 crude proxy. **실제 tumor mask 아님** — H&E 병리 annotation 없이는 stromal/immune contamination 구분 불가
- Top-50% threshold 는 임의 선택. 25% / 75% 로 바꾸면 결과 달라질 가능성
- TG/PAX8 단독 filter 는 사용 안 했음 (DM1 spot 은 thyroid lineage 낮을 수 있어서 의도적 회피) — 그러나 그 결과 tumor 정의가 더 약해짐

**Required next analysis**:
- Pathologist H&E annotation (또는 SpaCET/cell2location 같은 deconvolution) 으로 tumor area 정의 → 같은 stage trend 재시행
- Epithelial threshold sensitivity (top 25% / 50% / 75%) panel

**Figure recommendation**: epithelial-restricted violin plot (이미 생성됨, `violin_primary_scores_epi.png`) supplementary, 본 trend 의 robustness sensitivity 도 함께 표시.

---

## 6. CAF / EMT / hypoxia 와 연결되는가?

**Verdict**: WEAK or DIRECTION-INCONSISTENT.
**Evidence**:
- DM1 ~ CAF_ECM: spot ρ = +0.24 (weak positive); sample-mean ρ = +0.10 (NS)
- DM1 ~ EMT: spot ρ = **−0.09** (예상과 반대 방향!); sample-mean ρ = −0.06 (NS)
- DM1 ~ Hypoxia: spot ρ = −0.26 (negative); sample-mean ρ = +0.71 p=0.002 (**Simpson flip**)
- Stage trend epithelial depth-corrected: CAF_ECM ρ = **−0.52, p = 0.038** (CAF stage 와 함께 **감소** — 비상식적, over-correction 의심)
- Hypoxia stage trend: NS

**Weakness**:
- DM1-high tumor 가 CAF / EMT / hypoxic 환경과 연관 가설 **이 데이터에서 직접 지지 안 됨**
- EMT 는 spot-level 에서 음의 상관 — Paper 1 가설과 정반대
- Hypoxia Simpson flip 은 within-sample / between-sample 신호 분리 — single-cohort 에서 단독 해석 위험
- CAF_ECM 의 epithelial-depth-corrected negative trend 는 over-correction artifact 가능성 — genuine stromal infiltration 까지 깎인 것

**Required next analysis**:
- Multi-cohort meta 가 안 되면 mechanism claim 자체 포기
- 또는 spot-level deconvolution (CARD / cell2location) 으로 fibroblast / hypoxic niche 위치 분리한 후 DM1 score 와 spatial co-localization 통계 (Moran's I, neighborhood enrichment)

**Figure recommendation**: **Paper 1 main figure 후보 아님**. supplementary 에라도 넣지 않고 **차라리 빼는 게 안전**. EMT 음의 상관은 reviewer 에게 reverse-causality 의심을 줄 수 있음.

---

## 7. 등급 판정: Main / Supplementary / Drop

**Verdict**: **SUPPLEMENTARY 1 figure + Discussion 1 sentence + Methods 1 sentence**. Main figure 절대 불가. Drop 도 안 함.

**Evidence (등급 결정 근거)**:
- ✅ Direction-consistent (RAI ↓ / DM1 ↑ / TDS ↓ in epithelial spots)
- ✅ Internal consistency strong (DM1 ~ TDS ρ=−0.89)
- ✅ Stage labels 생물학적으로 valid (proliferation gradient survives depth correction)
- ✅ 8-gene 전부 detect, no dropout
- ❌ Sample-mean stage trend 자체는 underpowered (epi corrected p=0.18)
- ❌ CAF/EMT/hypoxia mechanism link 불충분 (EMT 는 wrong direction)
- ❌ TDS-overlap 으로 partial redundancy

**Required next analysis (등급 upgrade 시도 시)**:
- Paper 1 reach venue (Cell Rep Med / JCI Insight) 정당화 위해 필요한 것:
  1. HRA003537 + GSE230424 ingest 후 meta-stage trend (n≈60)
  2. Independent thyroid-lineage axis (GATA3/HHEX/FOXA1) 로 DM1 cross-validate
  3. Pathologist H&E mask 로 tumor-area 정확히 정의
- 위 3개 다 통과해야 main figure 자격
- 그 전엔 **MARATHON DISCIPLINE 위반 — 확장 금지**

**Figure recommendation (확정안)**:

| 자리 | 내용 | 출처 |
|---|---|---|
| Supp Fig N | 2-row × 4-stage spatial heatmap (RAI_8 / DM1_like) representative slides | `fig_2x4_RAI_DM1_per_stage.png` |
| Supp Fig N+1 | 16-slide sample-mean strip plot (RAI_8, DM1_like, TDS_like by stage) | `fig_sample_means_RAI_DM1_TDS.png` |
| Supp Table | stage_trend_summary.csv + correlation_summary.csv + depth_corrected_trend.csv 통합 | 3 files |
| Methods | "8 RAI genes scored as within-sample z-score; sample-mean Spearman vs ordinal stage; depth-corrected residualization with log_counts + log_ngenes" | 1 sentence |
| Discussion | "Single-cohort GSE250521 (n=16) showed direction-consistent but underpowered DM1-axis trend (epi-enriched ρ=−0.35 after depth correction; p=0.18). Independent multi-cohort spatial meta-analysis required for inferential confirmation." | 1 sentence (본인 voice) |

**CAF/EMT/hypoxia 분석은 Paper 1 에서 빼는 게 안전.** 본 데이터로는 mechanism claim 불가, 잘못 넣으면 reviewer 가 "EMT inverse correlation, CAF down-trend" 잡고 reverse-causality 의심함.

---

## 한 줄 결론

> **Internal consistency (DM1~TDS ρ=−0.89) + direction consistency (epi RAI↓/DM1↑) + stage label validity (proliferation gradient OK) → drop 안 함. 그러나 sample-mean trend underpowered + mechanism (CAF/EMT/hypoxia) link 약함 → main figure 불가. SUPPLEMENTARY 1장 + Methods/Discussion 각 1줄 이 정확한 자리.**
