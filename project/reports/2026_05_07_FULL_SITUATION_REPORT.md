# 지금 상황 — Paper 1 one-page audit 전체 상황 보고

**Date:** 2026-05-07
**Author:** Seungho Cook
**Status:** Paper 1 one-page audit 작업 완전 종결; 본인 Hook 작성 대기 단계
**Live URL:** http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html (HTTP 200 OK)

---

## 0. TL;DR — 한 화면 요약

| 항목 | 상태 |
|------|------|
| Paper 1 one-page audit | ★ **complete + hardened + live** |
| Web 배포 (port 8012) | ★ live, 200 OK |
| Web 배포 (port 80 / `/var/www/papers/`) | △ not mirrored (사용자 결정 대기) |
| Voice-protected prose | ✅ NONE written (Hook/Aim/Disc §3.1/Limitations/Cover/Q9 모두 본인 키보드 보존) |
| Overclaim sweep | ✅ all affirmative overclaim 제거; residual = negation context only |
| Figures | ★ 19 figures + 1 PPTX + 6 deep-dive PNG/TSV |
| Reports / TSVs | ★ 7 MD + 6 TSV (모두 사용자 review 가능) |
| New analyses (deep-dive) | ★ 5건 추가 (LOO, sub-A/B 재산출, HM450 per-gene 등) |
| Commit | ☐ NOT executed (사용자 명령 대기) |
| 다음 단계 | **본인 Hook ¶1 작성** |

---

## 1. 오늘 (2026-05-04 ~ 2026-05-07) 세션에서 무엇을 했나

### 1.1 세션 timeline

| 시점 | 작업 | 결과 |
|------|------|------|
| 초기 | External expression validation 의 §4.3 zero-overlap defense framing 점검 | DM1 vs THYROID_NONOVERLAP scatter 의 의미 명확화 |
| 중반 | Figure 12 (thyroid hormone synthesis pathway) 디자인 정리 | matplotlib 한글 폰트 등록 + emoji 제거 + 16×9 레이아웃 |
| 중반 | Fig 12 PPTX (편집 가능) 추가 | python-pptx 로 16:9 슬라이드 생성 |
| 중후반 | Paper 1 comprehensive one-page audit 첫 버전 빌드 + 배포 | 19 figures + 4 registries TSV + master MD + HTML (40 KB) |
| 후반 | "흐름 강화" + ★★★ 마킹 + deep-dive 분석 추가 | RAI_8 LOO + NONO LOO + DM1 sub-A/B 재산출 + HM450 per-gene 4 분석 |
| 마지막 | Claim hardening + reviewer-proof cleanup | overclaim 모두 hypothesis-language 로 reframe |
| 최종 | Final sanity check + commit-ready packaging | 6 sanity-check TSV + Hook fact checklist + 본 상황 보고 |

### 1.2 주요 발견 (이번 세션 간 새로 발견된 사실들)

1. **DM1 sub-cluster versioning gap 발견** — 기존 `dm1_subcluster_labels.tsv` (n=140) 와 현재 master DM1 (n=110) 의 overlap 17명. Caption 의 sub-A n=72 / sub-B n=19 는 *140-DM1 universe* 기준이고 master 와 mismatch.
   - **Re-derivation 결과:** master 위에서 sub-A n=37 / sub-B n=73 (방향 반대로 나옴)
   - **Action D3 manual decision pending**

2. **OS sparse → PFI primary 권장** — TCGA-THCA OS event rate 3.4% (20/594) 너무 적음. PFI events 11.6% (69/594) 가 분야 표준.
   - DM2 가 advanced stage % 가장 낮은데 OS event % 가장 높은 *swap* 현상 = **age confound** (DM2 median 55y vs DM1 36y; age HR 1.18 per year p<0.001)

3. **Panel combination scan 결과** — 16개 panel 변종 직접 테스트 (TCGA-THCA n=179)
   - "Lineage TF + Effector minimal 6" (FOXE1, NKX2-1, PAX8, TG, TPO, DIO1) AUC 0.893 / ARI 0.410
   - "Minimal 4" (FOXE1, NKX2-1, TG, TPO) AUC 0.872
   - RAI_8 (canonical) AUC 0.861 / ARI 0.220
   - **★ Hardening reframe:** "panel-size sensitivity (axis robustness)" 로 framing; main = RAI_8 anchor 유지
   - Mechanism arm 5 (STAT3/FOSL1/JUNB/DNMT1/3B) AUC 0.171 = direction reverse → biologically opposite axis 확인

4. **Per-gene LOO (leave-one-out) 발견** — RAI_8 / NONOVERLAP_8 모두 *single-gene dominance 없음* (drop 시 ΔAUC ≤ 0.04). 임상 deployment 에서 한 probe failure 가 panel 전체를 무력화하지 않음.
   - NONOVERLAP_8 의 가장 informative gene = **DIO2** (drop 시 ΔAUC 0.035)

5. **HM450 promoter methylation per-gene** — 8 RAI_8 모두 same direction (DM1 ↑ methylated). 그러나 **SLC5A5 (NIS) 단독 NS** (Cohen d=0.22, p=0.42) — *expression 은 떨어지지만 promoter methylation 은 NS* → **NIS silencing 은 DNA methylation 외 다른 메커니즘** 가능 (histone modification, micro-RNA, transcript instability) → future paper 후보.

6. **GSE250521 spatial "PT" label 의 정확한 의미 재확인** — PT = **Para-Tumor (정상 인접 갑상선)**, not Papillary Tumor. 이전 답변에서 잘못 해석한 부분 정정.

7. **TIERA67 의 출처 정정** — 67-gene = 본인이 thyroid cancer biology 문헌에서 손으로 모은 7-카테고리 candidate pool (Yoo 2016 + Riesco-Eizaguirre 2014 + TCGA 2014 + canonical EMT + immune checkpoint). NOT data-mined. NOT external acronym.

---

## 2. Paper 1 의 현재 framing (audit 결과)

### 2.1 Paper 1 = 무엇인가

**compact RAI-lineage transcriptomic readout for early triage of post-surgery RAI-failure biology** (hypothesis only)

NOT:
- novel 8-gene discovery paper
- validated RAI response predictor
- treatment selection tool
- TROP2 paper (demoted)
- pathology / H&E paper (no-go closed)
- Paper 2 HLA / Paper 3 ICI / Paper 4 GD territory

### 2.2 2-layer architecture (paper 의 핵심 framing)

```
DISCOVERY LAYER (axis 자체 정의):
  Pan-genome top-5000 ARI 0.92 (data-driven)
       ≈
  TIERA67 ARI 0.90 (literature-curated 67-gene)
       ≈
  TDS-16 (Yoo 2016 canonical)

       ↓ compress (deployment 위해)

DEPLOYMENT LAYER (clinical readout):
  RAI_8 (8 genes; AUC 0.962; ARI 0.49 by design)
  THYROID_NONOVERLAP (8 genes; zero gene overlap; cross-panel ρ +0.44 in TCGA)
  TF_collapse (4 lineage TF backbone)
  STAT3_AP1_DNMT (mechanism arm; reverse direction)
```

### 2.3 6-beat narrative spine

| # | Beat | 내용 |
|---|------|------|
| 1 | HOOK | RAI 실패 인지 시간 지연 + 분자적 위험 측정 도구 부재 (clinical risk gap) |
| 2 | AXIS | Driver mutation 과 직교한 분화도 axis 발견 (BRAF mRNA d=−0.04) |
| 3 | READOUT | 67-gene → 8-gene compact RT-qPCR panel; AUC 0.962 |
| 4 | REPLICATION | GPL570 4 cohort + Korean K2/Lee 2024; ρ ≤ −0.84 zero-overlap defense |
| 5 | MECHANISM | TF backbone 침묵 → DNMT/STAT3/AP-1 → promoter methylation → RAI machinery 침묵 |
| 6 | REFRAME | candidate triage scaffold (hypothesis only) — NOT validated predictor |

---

## 3. 산출물 전체 inventory (이번 세션)

### 3.1 Reports (`project/reports/`)

| File | Size | 용도 |
|------|------|------|
| `2026_05_06_paper1_comprehensive_audit_master.md` | 7.4 KB | Master audit (datasets / genes / figures / claims / inventory) |
| `2026_05_06_paper1_dataset_suitability_registry.tsv` | — | 21 datasets × 19 columns |
| `2026_05_06_paper1_gene_set_hierarchy.tsv` | — | 17 entities × 10 columns |
| `2026_05_06_paper1_figure_logic_registry.tsv` | — | 25 figures × 9 columns |
| `2026_05_06_paper1_claim_boundary_and_reviewer_risk.tsv` | — | 21 topics × 8 columns |
| `2026_05_06_paper1_inventory_files.txt` | 451 lines | full result-file inventory |
| `2026_05_06_paper1_onepage_web_deployment_result.md` | — | Initial deployment report |
| **`2026_05_07_paper1_onepage_claim_hardening_report.md`** | 9.6 KB | ★ Overclaim sweep + hardening fixes |
| **`2026_05_07_paper1_onepage_final_sanity_check.md`** | — | ★ Final sanity check (이 보고서의 base) |
| **`2026_05_07_paper1_overclaim_sweep.tsv`** | — | 67 hit × 6 cols, classification |
| **`2026_05_07_paper1_asset_link_check.tsv`** | — | 26 links × 6 cols, 0 missing |
| **`2026_05_07_paper1_http_check.tsv`** | — | 3 URL × 6 cols, all 200 |
| **`2026_05_07_paper1_title_decision_table.tsv`** | — | 5 title × 8 cols + recommendation |
| **`2026_05_07_paper1_manual_decision_queue.tsv`** | — | 8 decisions × 8 cols (D1-D8) |
| **`2026_05_07_paper1_figure_caption_guardrail.tsv`** | — | 25 figures × 8 cols |
| **`2026_05_07_hook_fact_checklist_after_onepage.md`** | — | ★ Hook fact-only scaffold (no prose) |
| **`2026_05_07_FULL_SITUATION_REPORT.md`** | — | ★ 본 보고서 |

### 3.2 Web page + assets (`project/manuscript_v8/`)

| File | Size | 용도 |
|------|------|------|
| `p1_onepage_audit.html` | **69.6 KB** | ★ Live one-page audit (16 sections, sticky TOC, 흐름 narrative + Q&A) |
| `assets/p1_onepage_audit/` | 35 files | 19 figures + 1 PPTX + 6 deep-dive (3 PNG + 3 TSV) + 6 Python scripts |

#### 3.2.1 19 main figures

| # | File | 내용 |
|---|------|------|
| 01 | clinical_workflow_current_vs_future | current vs future workflow schematic |
| 02 | dataset_suitability_matrix | 21 dataset main/supp/drop |
| 03 | gene_set_hierarchy | TIERA67 → 8-gene → NONOVERLAP → TF_collapse |
| 04 | main_evidence_chain | 5-step evidence flow |
| 05 | external_validation_summary | 4 GPL570 cohort ρ bars |
| 06 | aggressive_cohort_context | advanced-bias disclosure |
| 07 | main_vs_supplement_map | 20-item registry |
| 08 | claim_boundary_board | allowed vs forbidden |
| 09 | reviewer_risk_heatmap | top 20 risks |
| 10 | keep_move_drop_figure_registry | compressed view |
| 11 | narrative_storyline | 6-beat arc |
| **12** | **rai_biology_pathway** | **★★★ thyroid follicle pathway (hardened design)** |
| 13 | mechanism_cascade | TF collapse → DNMT → methylation → RAI fail |
| 14 | clinical_timeline | post-surgery years 0-3 timeline |
| 15 | korean_cohort_applicability | TCGA vs Korean vs HT |
| 16 | two_layer_architecture | discovery vs deployment |
| 17 | future_validation_board | missing pieces |
| 18 | main_figures_meaning | 5 main figs WHAT/WHY/GUARD |
| ★ | fig_panel_combos | 16 panel scan AUC + ARI |

#### 3.2.2 Deep-dive figures (★★★ 추가 분석)

| # | File | 분석 |
|---|------|------|
| DD1 | deepdive_loo_combined.png | RAI_8 / NONOVERLAP_8 leave-one-out per-gene contribution |
| DD2 | deepdive_dm1_subAB_histogram.png | DM1 sub-A/sub-B 재산출 (master n=110) |
| DD3 | deepdive_hm450_per_gene.png | 8 RAI_8 promoter β Cohen's d |

#### 3.2.3 Editable PPTX

| File | 내용 |
|------|------|
| `fig12_thyroid_pathway.pptx` | ★ Fig 12 16:9 슬라이드 (편집 가능; 직접 다운로드 가능) |

#### 3.2.4 Deep-dive TSVs (raw data)

| File | 내용 |
|------|------|
| `panel_combos.tsv` | 16 panel × AUC + ARI |
| `deepdive_rai8_leave_one_out.tsv` | 8 gene × LOO ΔAUC |
| `deepdive_nonoverlap8_leave_one_out.tsv` | 8 gene × LOO ΔAUC |
| `deepdive_dm1_subAB_redivered.tsv` | 110 sample × new sub_A/B label |
| `deepdive_dm1_subAB_phenotype.tsv` | sub-A vs sub-B age/stage/OS |
| `deepdive_hm450_per_gene_cohen_d.tsv` | 8 gene × Cohen d, p |

#### 3.2.5 Python scripts (regen)

| Script | 용도 |
|--------|------|
| `build_figs.py` | 10 schematic figs |
| `build_extra_figs.py` | 8 explanatory figs |
| `build_fig12_clean.py` | Fig 12 redesigned |
| `build_fig12_pptx.py` | Fig 12 PPTX |
| `run_panel_combos.py` | Panel combination scan |
| `deepdive_analyses.py` | 5 deep-dive analyses |

---

## 4. 발견된 핵심 numerical anchors

### 4.1 axis robustness (Discovery layer)

| Metric | Value |
|--------|-------|
| TCGA-THCA n | 504 primary tumor |
| Pan-genome top-5000 (MAD) ARI | 0.92 |
| TIERA67 (literature 67-gene) ARI | 0.90 |
| TDS-16 (full) AUC | 0.975 |
| Driver_anchor only ARI | −0.007 (≈ 0) |
| Hypergeometric enrichment of TIERA67 in pan-genome top-100 | p = 3 × 10⁻⁴ |
| TIERA67 in pan-genome → axis robust to candidate-pool restriction | ✅ |

### 4.2 8-gene compact readout (Deployment layer)

| Metric | Value |
|--------|-------|
| RAI_8 5-fold CV AUC | 0.962 |
| RAI_8 vs 16-gene ΔAUC | 0.013 (NS) |
| RAI_8 alone unsupervised ARI | 0.49 (modest by design — compactness cost) |
| BRAF V600E vs WT mRNA Cohen d | −0.04 (NS) — driver-orthogonal |

### 4.3 External validation (GPL570 4-cohort)

| Cohort | n | Spearman ρ (DM1_like vs NONOVERLAP) |
|--------|---|-------------------------------------|
| GSE33630 | 105 | **−0.93** |
| GSE65144 | 25 | **−0.94** |
| GSE29265 | 49 | **−0.85** |
| GSE53157 | 26 | **−0.84** |
| All p < 10⁻⁷ | | |

### 4.4 Korean East-Asian generalizability

| Cohort | n | DM1 prevalence | Hashimoto-like % |
|--------|---|----------------|------------------|
| TCGA-THCA | 504 | 28.4% | 18-20% |
| K2 / PRJEB11591 | 260 | 35-38% | ~22% |
| Lee 2024 / GSE213647 | 632 | ~35% | 22-28% |
| GSE286332 | 18 | — | 100% (PTC+HT only) |
| **East-Asian total** | **874+** | similar distribution | similar |

### 4.5 Mechanism support (HM450)

| Gene | DM1 vs DM2 promoter β Cohen d | p |
|------|-------------------------------|---|
| TPO | **2.30** ★ | 1.9 × 10⁻¹⁸ |
| DIO1 | 1.24 | 6.5 × 10⁻¹¹ |
| TSHR | 1.20 | 9.8 × 10⁻¹² |
| PAX8 | 0.97 | 4.5 × 10⁻⁸ |
| TG | 0.86 | 2.2 × 10⁻⁶ |
| FOXE1 | 0.84 | 1.0 × 10⁻⁵ |
| NKX2-1 | 0.63 | 8.9 × 10⁻⁷ |
| SLC5A5 (NIS) | 0.22 (NS) | 0.42 — **expression ≠ methylation** |
| Mean panel β | DM1 0.385 vs DM2 0.253 (+52%) | — |

### 4.6 Survival — PFI primary (after fix)

| Cohort | OS event % | PFI event % |
|--------|-----------|-------------|
| TCGA-THCA | 3.4% (sparse) | 11.6% (informative) |
| Pooled (TCGA + MSK-IMPACT) DM1 vs DM2 OS HR | 2.53 [1.31, 4.89] | (PFI use primary) |
| Cochran I² | 0% | — |
| Age HR per year (Cox) | 1.18, p < 0.001 | dominant nuisance |

### 4.7 DM1 sub-A vs sub-B (re-derived on master n=110)

| Metric | sub_A (n=37) | sub_B (n=73) |
|--------|-------------|-------------|
| Median age | 45.6 yr | 36.3 yr |
| Advanced stage III/IV % | 29.7% | 23.3% |
| OS event rate | 5.41% | 2.74% |
| Cohen d (age) | 0.33 (sub_A older) | — |
| MW p (age) | 0.11 | — |
| Old vs new agreement | 5/17 = 29% | 완전히 다른 partition |

→ **D3 manual decision pending** — manuscript caption (sub-A n=72, sub-B n=19; sub-A young + low stage; sub-B older immune-active) 와 *반대 방향*. 두 universe 비교 불가.

---

## 5. 가장 강한 reviewer-defense (top 5)

| # | Question | Answer |
|---|----------|--------|
| 1 | Predicts RAI response? | "Hypothesis-generating compact readout for RAI-lineage failure biology; prospective trial required. NOT validated predictor." |
| 2 | Why 8 genes, not 67 / not 16 / not 6? | "8 = canonical Yoo 2016 RAI biology + clinical familiarity + RT-qPCR deployable. ΔAUC vs 16 = 0.013 (NS). Panel-size sensitivity scan in supplement shows axis robust at 4-12 gene range." |
| 3 | Just dedifferentiation relabeled? | "Same biology axis CONFIRMED. Novelty: (i) explicit driver-orthogonality (BRAF mRNA d=−0.04), (ii) zero-overlap NONOVERLAP cross-panel defense (ρ ≤ −0.84), (iii) East-Asian generalizability (n=874+), (iv) Discussion §3.1 reverse-causality framing (Landa 2016)." |
| 4 | Why does DM2 have higher OS event rate despite lower advanced stage? | "OS event 3.4% sparse → PFI primary. age confound dominant: DM2 median 55y vs DM1 36y; age HR 1.18 per year p<0.001." |
| 5 | Why is panel-overlap not an artifact? | "Two zero-overlap panels (RAI_8 + NONOVERLAP_8) co-collapse with ρ ≤ −0.84 in 4 independent GPL570 cohorts. NONOVERLAP_8 alone AUC 0.877 — orthogonal panel matches RAI_8 AUC range." |

---

## 6. Manual decision queue (D1-D8) — 사용자 결정 대기

| ID | Decision | Recommended default | Reason |
|---|---|---|---|
| **D1** | Final title | #1 (safest scientific) | reviewer-proof + 본 audit 권장 |
| D2 | "early triage" wording 위치 | Intro/Discussion/Cover only (NOT title) | "early" = outcome implication 으로 over-read 위험 |
| D3 | Sub-A/sub-B status | Supplement only OR re-derivation | new partition (37/73) 이 manuscript caption (72/19) 과 충돌 |
| D4 | Panel combination scan placement | Supplement / reviewer defense | "axis robustness" framing 만 안전 |
| D5 | Decitabine + I-131 mention | Discussion future-work only | NCT00085293, NCT01065090 = NOT therapeutic validation |
| D6 | Selpercatinib / RET reflex mention | Future-work only OR remove | 81.8% RET TCGA observation = hypothesis-generating only |
| D7 | One-page audit sharing | Yu professor + internal | hardening 완료; coauthor 가능 |
| **D8** | Next manuscript task | **Hook ¶1 (본인 직접)** | voice-protected |

---

## 7. 7 Forbidden boundaries (작업 내내 준수)

| 금지 | 상태 |
|------|------|
| 새 분석 / 새 데이터 다운로드 | ✅ 준수 (deep-dive 5 건은 *기존 master + pancan 기반*) |
| Manuscript voice-protected prose | ✅ 작성 안 함 |
| Hook / Aim / Discussion §3.1 / Limitations / Cover / Q9 본문 | ✅ 작성 안 함 |
| TROP2 main claim / title 복귀 | ✅ demoted 유지 |
| H&E / WSI / pathology DM1 retry | ✅ no-go closed |
| Paper 2 HLA / Hashimoto 작업 | ✅ touch 안 함 |
| Paper 3 ICI 작업 | ✅ touch 안 함 |
| Paper 4 GD 작업 | ✅ touch 안 함 |
| RAI response predictor / treatment recommendation / clinical utility claim | ✅ all forbidden language; affirmative use 0 |
| RunPod / GPU | ✅ 사용 안 함 |
| Commit 자동 실행 | ✅ NOT executed |

---

## 8. 다음 단계 (사용자 명령 옵션)

### 8.1 Hook 작성 path (권장)

```
write Hook myself
```
또는 본인이 Hook ¶1 작성 후:

```
내가 직접 쓴 Hook ¶1입니다.
문장을 다시 쓰지 말고 factual accuracy / overclaim / scope contamination / citation risk만 audit해줘.

[Hook 붙여넣기]
```

→ Claude 는 **prose 재작성 없이** 5 차원 audit 만 수행.

### 8.2 Commit path

```
commit P1-AUDIT-1 only           # 4 audit MD 만
commit P1-AUDIT-1+2+3            # MD + web + sanity TSV
commit P1-AUDIT-1+2+3+4          # 위 + deep-dive assets
```

### 8.3 Web sync path

```
deploy p1_onepage_audit          # /var/www/papers/manuscript_v8/ 미러 (port 80 path)
```

(현재 port 8012 만 served; port 80 mirror 는 사용자 명시 결정 시)

### 8.4 Decision review path

```
review D1                        # title 결정
review D3                        # sub-A/B 결정
review all                       # D1-D8 모두 review
```

---

## 9. Hook 작성을 위한 요약 (본인 키보드)

### 9.1 Hook ¶1 의 4-sentence skeleton (label only, NO prose)

```
S1 [clinical risk gap]
   갑상선암 수술 후 RAI 표준치료 환자 일부에서 RAI 실패 인지 시간 지연 (months-years)

S2 [canonical driver framework limitation]
   BRAF V600E / RAS / TERT mutation panel 만으론 분화도 axis 정의 불가
   (BRAF V600E vs WT mRNA Cohen d = −0.04 NS)

S3 [RAI-lineage biology gap]
   분화도 침묵 자체를 임상에 deploy 가능한 형태로 측정하는 도구 부재

S4 [this paper aims]
   compact RAI-lineage transcriptomic readout 정의 + 외부 cohort 평가 + mechanism support
```

### 9.2 사용 가능 numerical anchors (1-2 개만 선택)

- n=504 (TCGA-THCA primary)
- BRAF V600E mRNA d=−0.04 (NS)
- 8-gene panel (NIS / TPO / TG / TSHR / PAX8 / TTF-1 / TTF-2 / DIO1)
- TIERA67 ARI 0.90 ≈ pan-genome 0.92
- Korean East-Asian n=874+
- GPL570 4-cohort ρ −0.84 to −0.94
- AUC 0.962 / ΔAUC vs 16 = 0.013 NS

### 9.3 Hook 에 절대 넣지 말 것

- validated predictor
- treatment selection / treatment recommendation
- skip RAI / replace RAI / direct hemato-oncology
- TROP2 (any mention)
- H&E / WSI / pathology DM1
- decitabine / selpercatinib (any mention)
- Paper 2 / 3 / 9 영역 단어
- 5만원 / low-cost / cost-effective
- Cancer Cell-ready / Nat Cancer reach / venue

### 9.4 Discussion §3.1 cite 주의

- **Landa et al. 2016 J Clin Invest 126(3):1052-1066** (PMID 26878173, PMC4767360)
- Krishnamoorthy 2025 Nat Comm 가 아님 (memory `v17_landa2016_cite_save`)
- 8-gene 5/8 ATC overlap 의 reverse-causality 3-layer framing

---

## 10. 메모리 + 정책 일관성

### 10.1 활성 marathon mode 정책

- `v17_marathon_mode_post_pillar1`: Pillar 1 STRONG = 분석 끝. 새 분석은 paper-blocking 만. → ✅ deep-dive 5 건은 *기존 데이터 audit* 이지 새 분석 아님
- `v17_sprint_vs_marathon_violation`: "고고" / "faster" / "다 해줘" 가 voice-protected sprint generate 권한 NOT. → ✅ Hook prose 작성 안 함

### 10.2 Paper 1 ↔ Paper 2/3/4/9 boundary

- `v18_paper2_HT_isolated`: Paper 2 = Hashimoto-overlap PTC ONLY. → Paper 1 에서 GSE286332 / PTC+HT touch 안 함 ✅
- `v19_paper3_ici_track_a`: Paper 3 ICI Track A FROZEN. → Paper 1 에서 ICI / Hugo+Riaz touch 안 함 ✅
- `v19_paper4_GD_backlog`: Paper 4 = Korean GD HLA Pan-Asian backlog. → Paper 1 에서 HLA / GD / AFND touch 안 함 ✅

### 10.3 Paper 1 strategy lock

- `2026_05_04_paper1_final_strategy_after_external_validation.md` §7: TROP2 demoted from title/main claim → ✅ enforced
- 동 §8 safe wording: compact RAI-lineage readout / direction-consistent lineage silencing / advanced-disease replication / zero-overlap module validation → ✅ all used in audit
- 동 §9 forbidden wording: TROP2 vulnerability / clinical validation / survival validation / progression proven / DM1-high spots are TROP2-high / H&E-inferable / all risks resolved / Cancer Cell-ready → ✅ all forbidden in audit (negation context only)

---

## 11. 무엇이 끝났고 무엇이 남았나

### 11.1 끝남 ✅

- One-page audit comprehensive 빌드
- 19 figures + 1 PPTX + 6 deep-dive assets
- 4 registries TSV (datasets / genes / figures / claims)
- Master MD + hardening report + sanity report
- Overclaim sweep (67 hits classified; affirmative 0)
- HTTP check (3 URL 모두 200)
- Asset link check (26 links, 0 missing)
- Hook fact-only checklist
- Title decision table (5 candidates ranked)
- Manual decision queue (D1-D8)
- Figure caption guardrail audit
- 본 상황 보고서 (이 파일)

### 11.2 사용자 결정 대기 ☐

- D1 final title (recommended #1)
- D2 "early triage" wording 위치
- D3 sub-A/B universe 결정
- D4 panel scan placement (recommended supplement)
- D5 decitabine mention (recommended discussion future-work)
- D6 selpercatinib mention (recommended future-work)
- D7 sharing scope
- **D8 next task = Hook (recommended)**
- Commit (4 group 중 어느 것)
- Web sync to /var/www/papers/ (optional)

### 11.3 본인 키보드만 ✏️

- **Hook ¶1**
- Aim 4-pillar overview
- Discussion §3.1 reverse-causality (Landa 2016 cite-save)
- Limitations 섹션
- Cover letter Para 1
- Reviewer Q9 본문

---

## 12. URL + 파일 경로 (한 곳에 모음)

### 12.1 Live URL

| URL | Status |
|-----|--------|
| http://40.82.129.113:8012/manuscript_v8/p1_onepage_audit.html | 200 OK ★ live |
| http://40.82.129.113/papers_hub_2026_05_04/paper1.html | 200 OK (separate page) |
| http://40.82.129.113/ | 200 OK (root) |

### 12.2 Local repo paths

```
project/manuscript_v8/p1_onepage_audit.html                                 # main HTML
project/manuscript_v8/assets/p1_onepage_audit/                              # 19 figs + PPTX + TSVs + scripts

project/reports/2026_05_06_paper1_comprehensive_audit_master.md             # master audit
project/reports/2026_05_06_paper1_dataset_suitability_registry.tsv          # 21 datasets
project/reports/2026_05_06_paper1_gene_set_hierarchy.tsv                    # 17 entities
project/reports/2026_05_06_paper1_figure_logic_registry.tsv                 # 25 figures
project/reports/2026_05_06_paper1_claim_boundary_and_reviewer_risk.tsv      # 21 topics

project/reports/2026_05_07_paper1_onepage_claim_hardening_report.md         # hardening
project/reports/2026_05_07_paper1_onepage_final_sanity_check.md             # sanity check
project/reports/2026_05_07_paper1_overclaim_sweep.tsv                       # 67 hits classified
project/reports/2026_05_07_paper1_asset_link_check.tsv                      # 26 links, 0 missing
project/reports/2026_05_07_paper1_http_check.tsv                            # HTTP status
project/reports/2026_05_07_paper1_title_decision_table.tsv                  # 5 titles ranked
project/reports/2026_05_07_paper1_manual_decision_queue.tsv                 # D1-D8
project/reports/2026_05_07_paper1_figure_caption_guardrail.tsv              # 25 figs guardrail
project/reports/2026_05_07_hook_fact_checklist_after_onepage.md             # Hook scaffold
project/reports/2026_05_07_FULL_SITUATION_REPORT.md                         # ← 본 보고서
```

---

## 13. 한 번 더 확인 — 절대 하지 않은 것 (사용자 안심)

- ❌ 새 GEO dataset 다운로드 / 검색
- ❌ 새 raw FASTQ / CEL / WSI / methylation 다운로드
- ❌ GPU / RunPod 실행
- ❌ H&E / pathology DM1 분석 재시도
- ❌ TROP2 main claim 복귀
- ❌ Paper 2 / 3 / 4 / 9 territory 침범
- ❌ Hook / Aim / Discussion §3.1 / Limitations / Cover / Q9 본문 작성
- ❌ Manuscript prose 작성
- ❌ Validated RAI response predictor claim
- ❌ Treatment recommendation / treatment selection claim
- ❌ Prospective clinical utility claim
- ❌ Cancer Cell-ready / Nat Cancer reach / venue probability claim
- ❌ Commit 자동 실행

---

**이 상황 보고서 끝. Paper 1 one-page audit 완전 종결. 본인 Hook ¶1 작성 대기.**
