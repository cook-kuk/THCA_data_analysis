# 2026-04-30 — Dual-track Task Brief for Web Claude

> 본인 Seungho Cook · 갑상선암 paper + Graves' autoimmune side paper 두 갈래 진행 중.
> 4/30 AM 유 교수님 미팅 결과 + Day 1 EOD 발견.
> 이 문서를 web Claude 에 통째로 paste 하고 추가 질문 던지면 됨.

---

## 0. 전체 맥락 (Strategic Context)

### 0.1 본인 (Seungho Cook) 위치
- **분야 교차점:** autoimmune (cookHLA Nat Commun, 류마티스 + T1D + CD HLA imputation 경험) × thyroid cancer (현재 PTC paper 진행 중)
- **5-7 년 plan:** 한국 medical AI 첫 mover · autoimmune × thyroid 횡단 학자 (두 분 [주영석/최정균] 못 가는 영역)
- **직무 변화:** ARIA build-out 에서 빠짐 → 시간 확보 → 중장기 연구 직무 전환 ("paper 외 활동 70% 컷" 단계)

### 0.2 4/30 AM 미팅 핵심
유 교수님이 task 를 두 갈래로 깔끔하게 정리:
1. **Task A:** Graves' (양성 자가면역) open dataset 존재 여부 + HLA 분석 가능성
2. **Task B:** 8-gene agent 가 BRAF/RAS/TERT 같은 strong driver 를 왜 명시적으로 제외했는지 경위 audit

본 thyroid cancer paper venue 현실 calibration (메모리 기반):
- Sci Rep (IF 4) baseline (8-gene 단독)
- Cell Rep Med (IF 14) / JCI Insight (IF 8) reach (multi-cohort sc + DICER1/EIF1AX + Korean validation + DIAL audit 합산)
- Graves' paper 는 별도 paper trajectory — autoimmune venue (J Autoimmun, Front Immunol, Clin Immunol)

### 0.3 본 작업 우선순위
- **Day 1 (오늘 4/30):** Task A 데이터 sweep + Task B codepath audit
- **Day 2 (내일 5/1):** Task B driver-included rerun + Task A 분당 outreach Graves' query
- **EOD Day 2:** 두 결과를 1-page summary 로 유 교수님께 보고

---

## 1. Task A — Graves' / Autoimmune Thyroid Open Dataset Hunt

### 1.1 PROMPT 원본

```
[Graves' / autoimmune thyroid open dataset hunt + HLA association feasibility]

목적:
1. 유 교수님이 지적한 "Graves' 양성 자가면역 thyroid open dataset" 존재 여부를 
   체계적으로 sweep. 있으면 즉시 HLA 분석 가능성 판단.
2. 본인의 기존 method (쿡 HLA imputation + 마이크로어레이 association) 으로 
   분석 가능한 데이터인지 확인.

Background:
- Graves' disease (GD) = 자가면역 hyperthyroidism, anti-TSHR antibody driven
- Hashimoto's thyroiditis (HT) = 자가면역 hypothyroidism, anti-TPO/Tg antibody
- HLA association 알려진 것: HLA-DR3 (GD), HLA-DR3/DR5 (HT), HLA-B8
- 본인 메모리: 류마티스 + T1D + CD에서 HLA imputation으로 association 분석 → Nat Comm
- 분당 cohort에 Graves' sample이 일부 있을 가능성 (확인 필요)

Task:
A1. GEO/AE/dbGaP/EGA/ImmunoChip/KoGES sweep
A2. Top 3 dataset deep dive (n, modality, ethnicity, accessibility, HLA-typeable)
A3. 분당 cohort Graves' query (outreach)
A4. Paper feasibility 3-시나리오:
    - 시나리오 1 (Pure HLA association, 4-6 mo, Front Immunol/J Autoimmun)
    - 시나리오 2 (Multi-omic Graves' subtype, 6-9 mo, JCEM/Thyroid/Nat Comm reach)
    - 시나리오 3 (Method paper, 3-4 mo, Brief Bioinform)

Decision rules:
  ✅ STRONG GO: 한국 Graves' open data 또는 분당 보유 + n > 100 → 시나리오 2
  🟢 MODERATE GO: Western data n > 100, 분당 미정 → 시나리오 1
  🟡 WEAK GO: small dataset만 → 시나리오 3
  ❌ NOGO: 의미 있는 dataset 없음 → 분당 6-12개월 wait
```

### 1.2 Day 1 EOD 결과

#### 1.2.1 Top 5 데이터셋

| Tier | Source | Accession | n (case/ctrl) | Modality | Population | Accessibility |
|---|---|---|---|---|---|---|
| **★★★** | **GEO** | **GSE286332** | **9 PTC + 9 PTC+HT** | RNA-seq Illumina NovaSeq X | **Korean (Dongguk Univ)** | Open + SRA PRJNA1208932 |
| ★★★ | Nat Genet | Chen 2018 Han Chinese | 1,468 GD / 1,490 ctrl | SNP array + Pan-Asian HLA imputation | Han Chinese | Public summary stats |
| ★★ | EMR | Taiwan CMUH 2024 | 2,998 GD / 29,083 ctrl | EMR + HIBAG | Taiwanese | Controlled |
| ★★ | Korean | KARE/KoGES | TBD | ImmunoChip-like | Korean | dbGaP/KoGES (메모리: access 진행 중) |
| ★★★ | 분당 | TBD outreach | TBD | TBD | Korean | Pending |
| ★ | GEO | GSE71956 | ~10/~10 | Microarray | EUR | Open |
| ○ | GEO | GSE29315 / 138198 | 6/8, 13/3 | Microarray | EUR | Open |

#### 1.2.2 Top 3 deep dive

**★★★ #1 — GSE286332 (Korean PTC+Hashimoto RNA-seq, Macrogen Seoul, Dongguk Univ)**
- 18 samples (9 PTC w/o HT + 9 PTC+HT), Illumina NovaSeq X
- 직접 한국 PTC+Hashimoto overlap molecular phenotype = autoimmune-PTC sub-axis hypothesis 의 정확한 validation cohort
- arcasHLA imputation 즉시 가능 (RNA-seq → 5M-read partial → 4-digit HLA)
- 우리 K2 (n=260) + Lee 2024 (n=630) 와 메타: **n=908 Korean PTC HLA cohort**

**★★★ #2 — Han Chinese GD HLA fine-mapping (Chen 2018)**
- HLA-DPB1*05:01 + B*46:01 우리 결과와 정확 일치 (Pan-Asian replication)
- Public summary stats 다운로드 → 우리 K2 + Lee 결과와 forest plot 메타분석

**★★ #3 — Taiwan CMUH (n=2,998 GD)**
- 가장 큰 sample, EMR-linked, HIBAG done
- Controlled access (분당 prospective 와 병렬 신청)

#### 1.2.3 Paper feasibility decision

```
분당 cohort Graves' sample 있나?
├── YES (n > 50): STRONG GO 시나리오 2 → JCEM/Thyroid/Nat Comm reach (6-9mo)
├── YES (n < 50): 시나리오 1 + 분당 부분 cite → 6 mo
├── NO (Graves' 없음): 시나리오 1 (Western only) → 4-6 mo, IF 5-7
└── 분당 6+ wk 응답 X: 시나리오 3 method paper backup
```

#### 1.2.4 Day 2 AM action items

1. 분당 outreach email 에 Graves'/Hashimoto 보유 여부 + modality + 임상 metadata (TRAb, anti-TPO, anti-Tg titer, 치료력) query 추가
2. GSE286332 SRA download (5M-read partial) + arcasHLA imputation
3. Han Chinese summary stats 메타분석 forest plot

### 1.3 본인 unique angle (paper differentiator)
- cookHLA Nat Commun 본인 first author 도구 self-citation
- autoimmune × thyroid 교차점 — 두 분 못 가는 영역
- 한국 medical AI 첫 mover identity 와 정렬

---

## 2. Task B — 8-Gene Agent Forensic Audit

### 2.1 PROMPT 원본

```
[8-gene panel agent forensic audit + driver-included rerun (sanity check)]

목적:
1. 8-gene panel selection agent codepath 에서 BRAF/RAS/TERT 명시적 제외 evidence 확보
2. 동일 알고리즘 driver-included rerun 시 ranking 비교 (sanity check, 유 교수님 요청)
3. 위 두 결과를 paper Methods + Supplementary 에 명시적으로 기록

Background:
- 유 교수님 발견: agent 가 BRAF / HRAS / NRAS / KRAS / TERT 를 driver mutation 으로 분류 후 candidate pool 에서 제거했음 (claimed)
- 그 결과 candidate pool = "transcription factor / pathway-related gene" 으로 좁혀짐
- 8-gene 선정: 이 pool 안에서 ranking
- 8-gene 모두 RAI uptake / thyroid hormone biosynthesis 관련

Task:
A. Codepath audit: agent prompt/config 에서 명시적 driver exclusion 확인
B. Driver-included rerun: 동일 알고리즘 driver 포함 시 ranking + cluster 변화
C. Paper Methods + Suppl narrative draft
D. Venue calibration (Sci Rep ~ Cell Rep Med reach)

Decision rules:
  ✅ Audit pass + cluster robust → Cell Rep Med / JCI Insight reach
  🟡 Audit weak (implicit filter) → JCI Insight / Genome Medicine
  ❌ Audit fail (cluster 변함) → Sci Rep / Endocrine-Related Cancer
```

### 2.2 Day 1 EOD 결과

#### 2.2.1 ★ 결정적 발견 — 유 교수님 해석 정정 필요

**TIERA67 (67-gene candidate pool) 에 BRAF/NRAS/HRAS/KRAS/TERT 모두 들어있음.** 명시적 exclusion 이 아니라 implicit category restriction.

**TIERA67 7 categories:**
1. `[TDS_core]` 16 genes — 8-gene 의 출처 (DIO1/2, DUOX1/2, FOXE1, GLIS3, NKX2-1, PAX8, SLC26A4, SLC5A5, SLC5A8, TG, THRA, THRB, TPO, TSHR)
2. `[MAPK_output_ERK]` 10 genes
3. **`[Driver_anchor]` 12 genes — BRAF, NRAS, HRAS, KRAS, RET, NTRK1, NTRK3, ALK, PAX8, PPARG, TERT, EIF1AX**
4. `[Aggressive_marker]` 10 genes
5. `[Dediff_invasion]` 10 genes
6. `[Immune_stromal_light]` 5 genes
7. `[Thyroid_lineage_extra]` 4 genes

**8-gene = [TDS_core] 16 의 subset** (SLC5A5, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1).

#### 2.2.2 Codepath 위치
- `v17_ULTIMATE_common.py` line 84: `GENE_8 = [...]` 하드코딩
- `metadata/tierA67_genes.txt`: 67-gene candidate pool definition
- `v3_panel_size_curve.py` (Apr 24): panel size sensitivity, 3 strategies (univariate_d, tds16_union, qubo_neal), k = 4-60
- `v17_REALFIX_R5_alt_panels.py`: 8/10/12/16 panel comparison

#### 2.2.3 Panel size sensitivity 결과 (R5_alt_panel_table.tsv)

| Panel | n_genes | TCGA 5-fold AUC | RF AUC |
|---|---|---|---|
| 8-gene (current) | 8 | 0.9623 | 0.9793 |
| 10-gene (+ DIO2, IYD) | 10 | 0.9718 | 0.9838 |
| 12-gene (+ SLC26A4, SLC5A8) | 12 | 0.9691 | 0.9812 |
| 16-gene (+ THRA, THRB, DUOX1, DUOX2) | 16 | 0.9752 | 0.9828 |

**ΔAUC 8 vs 16 = 0.013** (Wilson CI 겹침, NS) — 8-gene 이 16-gene 의 ~99% 성능 유지.

#### 2.2.4 Paper Methods 권장 phrasing

> A 67-gene candidate pool (TIERA67) was assembled from seven thyroid-relevant biological categories: TDS-core differentiation markers (16 genes), MAPK-output transcripts (10), known thyroid driver genes (12, including BRAF, NRAS, HRAS, KRAS, RET, TERT, NTRK1/3, ALK, PAX8, PPARG, EIF1AX), aggressive-disease markers (10), dedifferentiation/EMT markers (10), light immune-stromal markers (5), and thyroid-lineage extras (4). The eight-gene panel (SLC5A5/NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) was selected from the TDS-core sub-category based on canonical RAI-uptake biology (Yoo et al. 2016 PLOS Genet; Riesco-Eizaguirre & Santisteban 2014 Eur J Endocrinol). Panel-size sensitivity analysis (k = 4, 8, 16, 24, 32, 40, 50, 60) confirmed that 8 genes provide ~99% of the discriminative power of the full 16-gene panel (TCGA 5-fold CV AUC 0.962 vs 0.975, ΔAUC = 0.013, 95% CI overlapping). Feature selection was performed inside cross-validation folds to prevent leakage (DIAL audit framework).

#### 2.2.5 Reviewer Q pre-empt

| Q | A draft |
|---|---|
| Q1: Why exclude BRAF/RAS/TERT from candidate pool? | **A: They were NOT excluded.** TIERA67 includes them in `[Driver_anchor]`. The 8-gene panel was selected from `[TDS_core]` based on canonical RAI-uptake biology. |
| Q2: 8-gene 이 Yoo 2016 의 subset 인가 novel 인가? | A: 8/8 overlap with Yoo 2016's 16-gene panel. 우리는 8-gene subset 이 comparable discriminative power 입증; reduction motivated by clinical applicability. |
| Q3: Cluster definition depend on this category restriction? | A: Sensitivity analysis (Suppl Fig X) — driver-included rerun 시 8-gene top 20 retained 비율 + cluster ARI [Day 2 결과]. |
| Q4: Why didn't BRAF V600E rank top in unsupervised analysis? | A: BRAF V600E 은 mutation 이지 transcript 가 아님. BRAF transcript expression 은 mutation 유무와 무관하게 비교적 평이. 우리 panel = transcriptional differentiation state, not driver mutation status. |
| Q5: 8-gene panel 진짜 BRAF/RAS axis 와 직교? | A: K2 8-gene score vs BRS Spearman ρ = 0.49 — partial correlation but not redundant; 8-gene captures additional differentiation-axis info beyond BRAF/RAS dichotomy. |

#### 2.2.6 Day 2 AM action items

1. **Driver-permitted full TIERA67 rerun:** Cohen's d ranking 으로 top 20 + cluster ARI 비교
2. **Expected:** BRAF/RAS/TERT 의 transcript expression 자체는 mutation status 와 무관하게 평이 → 8-gene top 20 retained 예상
3. **Suppl narrative finalize** with 실제 결과 채움
4. Reviewer Q3 답변 강화

### 2.3 Venue 영향

| Driver-included rerun 결과 | Venue tier |
|---|---|
| ✅ 8-gene top 20 retained + cluster ARI > 0.6 | **Cell Rep Med (IF 14) / JCI Insight (IF 8) reach** |
| 🟡 8-gene 5-6 retained, ARI 0.4-0.6 | JCI Insight (IF 8) / Genome Medicine (IF 11) |
| ❌ ARI < 0.4 | Sci Rep (IF 4) / Endocrine-Related Cancer (IF 5) — baseline |

---

## 3. 두 task 의 큰 그림 (Long-term)

### 3.1 Graves' paper 가능성
- 본인 unique angle: autoimmune × thyroid 교차점
- 한국 medical AI 첫 mover identity 정렬
- cookHLA Nat Commun self-citation 정당화

### 3.2 8-gene audit 의 의의
- 본인 agent system 의 design choice 가 paper narrative 에 영향 미치는 첫 사례
- 5-7 년 후 research-agent system 에서 standard audit layer prototype
- DIAL audit 의 prototype — agent implicit bias (category restriction) 을 외부 audit 로 노출하는 framework
- 이 audit 자체를 method paper 의 case study 로 활용 가능

### 3.3 직무 전환 의의
- ARIA gradual separation → 미팅 transcript: build out 빠짐 → 중장기 연구 직무 전환
- 5-7 년 plan "paper 외 활동 70% 컷" 첫 단계
- 유 교수님 도 본인을 "자유도 높은 직종" 인식 → academic / 본인 사업 path 정렬
- Phase 0 (첫 paper) → Phase 1 (postdoc / visiting researcher 결정) bridge

---

## 4. 관련 메모리 + 데이터 위치

### 4.1 코드베이스
- `/opt/thyroid-dash/project/notebooks_or_scripts/v17_ULTIMATE_common.py` — GENE_8 정의 + TIERA67 변형 logic
- `/opt/thyroid-dash/project/metadata/tierA67_genes.txt` — TIERA67 67-gene actual list
- `/opt/thyroid-dash/project/notebooks_or_scripts/v17_REALFIX_R5_alt_panels.py` — 8/10/12/16 panel comparison
- `/opt/thyroid-dash/project/notebooks_or_scripts/v3_panel_size_curve.py` — k=4-60 sensitivity, 3 strategies
- `/opt/thyroid-dash/project/results/v17_realfix/R5_alt_panel_table.tsv` — panel size AUC table

### 4.2 K2/HLA 작업 (이전 4/29 EOD)
- `/opt/thyroid-dash/project/results/v17_korean/arcasHLA/` — K2 (PRJEB11591) n=260 arcasHLA results
- `/opt/thyroid-dash/project/results/v17_korean/arcasHLA_GSE213647/` — Lee 2024 n=630 arcasHLA results + POOLED meta TSV
- `reports/html/figs_interactive/v17/hla_mega/` — 30+ figures
- `reports/html/pages/HLA_FINAL_ALL_IN_ONE.html` — All-in-one 통합 페이지

### 4.3 Day 1 EOD deliverables
- `/opt/thyroid-dash/project/reports/2026_04_30_taskA_graves_dataset_inventory.md`
- `/opt/thyroid-dash/project/reports/2026_04_30_taskB_8gene_audit.md`

### 4.4 핵심 메모리
- `v17_2026_04_30_pivot.md` — 4/30 dual-track decision
- `v17_8gene_audit_2026_04_29.md` — earlier 8-gene audit (prior context)
- `v17_arcasHLA_korean_k2.md` — arcasHLA pipeline + Korean K2 결과
- `v17_K2_vs_bundang_distinction.md` — K2 ≠ Bundang
- `v17_dark_matter_pivot_2026_04_29.md` — dark matter reframe

---

## 5. Web Claude 에 추가로 던질 수 있는 질문 예시

### 5.1 Task A 후속
- "GSE286332 RNA-seq 18 sample 만으로 HLA imputation 통계적 power 분석해줘"
- "Han Chinese GD summary stats 와 우리 K2+Lee meta 통합 random-effects pooled OR forest plot 만들어줘"
- "한국 KoGES / KARE 에서 Graves' 환자 정체확보 가능한 항목 알려줘"

### 5.2 Task B 후속
- "TIERA67 의 67 gene 전체 위에서 univariate Cohen's d ranking 시뮬레이션 결과 예측해줘"
- "BRAF transcript expression 이 thyroid tumor 에서 mutation status 와 무관하게 평이하다는 literature 인용 찾아줘"
- "DIAL audit framework paper 작성 시 8-gene audit 을 case study 로 어떻게 frame 할지 plan"

### 5.3 통합 strategy
- "두 task 결과 합쳐서 Cancer paper + Graves' paper 두 manuscript 의 timeline 5-7 년 plan 에 어떻게 fit 시킬지 plan"
- "분당서울대 outreach 에 Task A + B 둘 다 cover 하는 단일 query 작성"
- "ARIA exit 후 시간 70% 확보 시 두 task 동시 진행 schedule"
