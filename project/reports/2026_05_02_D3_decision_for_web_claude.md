# 2026-05-02 — Day 3 Decision Brief for Web Claude

> 본인 Seungho Cook · 4/30 Yu professor 미팅 dual-track 후속.
> 4/30 Day 1 (Task A/B) + 5/1 Day 2 (P1-P5) + 5/2 Day 3 일부 (P4 audit + P5 mediation) 완료.
> Day 3 잔여 7-prompt 중 Prompt 1 (SRA + arcasHLA + Han Chinese forest meta) 이 Azure 비용 발생 → 본인이 결정 못 정함.
> 이 문서를 web Claude 에 paste 하고 결정 받기 위함.

---

## 0. 전체 맥락 (4/30 미팅 → 5/2 현재)

### 0.1 본인 위치
- **분야 교차점:** autoimmune (cookHLA Nat Commun 본인 first-author, 류마티스+T1D+CD HLA imputation 경험) × thyroid cancer (현 PTC paper 진행 중)
- **5-7 년 plan:** 한국 medical AI 첫 mover, autoimmune × thyroid 횡단 학자 (두 분 [주영석/최정균] 못 가는 영역)
- **직무 변화:** ARIA build-out 빠짐 → 시간 70% 확보 → 중장기 연구 직무 전환 ("paper 외 활동 70% 컷" 1단계)

### 0.2 Cancer paper venue (post-Day 2)
- 4-Pillar 기반 **Cell Rep Med (IF 14) / JCI Insight (IF 8) reach** — 분당 outreach 없이도 가능
- 분당 Graves' 데이터 추가시 Nat Commun reach
- Phase 2 Graves' paper 별도 trajectory (autoimmune venue)

### 0.3 4-Pillar 현황 (5/1 master 기준)
1. **Korean Pan-Asian HLA cohort (n=890):** K2 (PRJEB11591 n=260) + Lee 2024 (GSE213647 n=630) arcasHLA RNA-seq imputation 완료. DPB1*05:01 56% replicates Kim 2014.
2. **GSE286332 PTC vs PTC+HT 분자 dissection:** 10,380 DEGs, 8-gene RAI d=−1.60, HLA-II d=+3.65, IFN-γ FDR=2e-4, KEGG Type I diabetes FDR=0.
3. **Driver mRNA neutrality (P1):** BRAF d=−0.04 (NS); driver AUC < 0.61; driver-only cluster ARI=0.
4. **Pan-genome robustness (P4):** TIERA67 ARI=0.90 ≈ pan-genome top-5000 ARI=0.92; 8-gene alone ARI=0.49 (honest); driver-only ARI=0.

---

## 1. Day 3 status (5/2 현재 시점)

7개 prompt 중 2개 완료 (오전 quick wins).

### 1.1 ✅ Prompt 4 — Wang 2024 supplementary audit (★실은 Chen XF 2024)

| Item | Finding |
|---|---|
| 정확한 인용 | Chen XF et al. 2024 Endocr Connect 13(11):e240301 (PMID 39235852, PMC11562686). 책임저자에 yulongwang@fudan.edu.cn 있어서 "Wang 2024" 라고 메모리 태그됨. |
| Patient-level supplementary | ❌ 없음 (단일 PDF 208 KB, aggregate only) |
| Follow-up data (OS/RFS/RAI) | ❌ **없음.** 논문 명시: *"insufficient follow-up time; have not recorded recurrence or death cases yet."* |
| HR for prognostic claims | ❌ 없음 (descriptive stat + p-value only) |
| DICER1/EIF1AX | ❌ DICER1 미언급; EIF1AX 2 benign 만 |
| mol_subtype / BRS | ❌ 미사용 (pathology + mutation status only) |

**Decision: Scenario C 확정 — mutation landscape comparator only.** Discussion 한 줄 인용. Phase 2 author contact 가능 (3-6 mo wait).

K2 의 DICER1/EIF1AX 4/4 (4.4%) finding 은 **고유 — Asian alternative driver landscape claim 유지**.

File: `project/reports/2026_05_02_D3P4_wang2024_audit.md`

---

### 1.2 ✅ Prompt 5 — GSE286332 P(DM1) gradient + mediation

5/1 P3 에서 PTC+HT 18/18 모두 DM2 분류된 paradox 해소 시도. 결과 STRONG.

#### Spearman ρ (P_DM1 vs covariates, n=18)
| | ρ | p |
|---|---|---|
| **8-gene RAI** | **+0.84** | 1.4e-5 |
| HLA-II | −0.81 | 5.4e-5 |
| HLA-I | −0.75 | 4e-4 |
| immune | −0.70 | 1.4e-3 |

#### OLS decomposition (z-standardized)
- R² = **0.756** (R² adj 0.704, F=14.5, p=1e-4)
- HLA-II β=−0.130 p=0.039 ★
- g8_RAI β=+0.054 p=0.065
- immune β=+0.097 p=0.120

#### Single-predictor R²
- **HLA-II alone: 0.663** (87.7% of full)
- g8_RAI alone: 0.654
- immune alone: 0.570

#### Baron-Kenny mediation (5,000-iter bootstrap)

| Mediator | a | b | indirect | %mediated | Boot 95% CI | p_emp |
|---|---|---|---|---|---|---|
| **HLA-II** | +1.67 | −0.11 | **−0.18** | **140%** | [−0.31, −0.03] | **0.023** ★ |
| **g8_RAI** | −1.00 | +0.08 | −0.08 | 63% | [−0.15, −0.03] | **0.002** ★★ |
| immune | +1.61 | −0.07 | −0.11 | 88% | [−0.26, +0.04] | 0.120 NS |
| HLA-I | +1.47 | −0.08 | −0.11 | 87% | [−0.20, −0.006] | 0.042 ★ |

#### Within-PTC severity gradient (n=9)
- ρ(g8_RAI, P_DM1) = **+0.73** (p=0.025) — **pre-clinical Hashimoto-like spectrum 신호**
- ρ(HLA-II, P_DM1) = −0.53 (p=0.14, NS at n=9)

**Decision: STRONG.** "PTC+HT → HLA-II infiltration → P_DM1 ↓" mechanistic claim 확정. 5/1 P3 의 18/18 DM2 paradox 정량 해소 — DM_call 이 binary 일 뿐 underlying P_DM1 spectrum 은 causal 하게 structured. 8-gene RAI 가 parallel partial mediator (63%, p=0.002).

Files:
- `project/reports/2026_05_02_D3P5_pdm1_gradient.md`
- `project/results/d3p5_pdm1_gradient/{D3P5_summary.json, mediation_results.json, decomposition.json, ...}`

---

## 2. ★ 결정 필요 — Prompt 1 (Azure 비용 발생)

### 2.1 Task 정의

```
GSE286332 arcasHLA partial RNA-seq imputation + Han Chinese GD HLA forest meta

목적:
- GSE286332 18 sample 의 4-digit HLA allele 추출 → K2(260) + Lee(630) + GSE286332(18) = n=908 → 
  PTC: 917, PTC+HT: 9 별도 layer (또는 sub-analysis)
- Chen 2018 Han Chinese GD HLA fine-mapping summary stats 다운로드
- Pan-Asian forest plot meta: Korean PTC+HT vs Chen Chinese GD vs Korean baseline
- Pillar 1 → "Korean Pan-Asian HLA cohort n=908 + Pan-Asian Graves' replication" 확장
```

### 2.2 차단성

**5/1 master 13.1 immediate task** 로 명시되었지만 실행 보류 중. Pillar 1 (HLA arm) 을 n=890 → n=908 + Pan-Asian Graves' replication 으로 확장하는 핵심 분석.

만약 PTC+HT 9 sample 이 Han Chinese GD 와 비슷한 risk allele profile 보이면:
- "**PTC+HT 는 Graves'-adjacent HLA background**" 입증 → paper Discussion 에 강력한 mechanism layer
- Pillar 1 을 단순 Korean PTC HLA 에서 Pan-Asian autoimmune-overlap HLA 로 확장 가능

만약 다르면:
- "PTC+HT 는 Hashimoto-specific HLA, Graves' 와 별개 axis" → 여전히 유효한 finding (다만 narrative 다름)

### 2.3 인프라 옵션 비교

| 옵션 | 도구 위치 | 비용 | 시간 | Blast radius | 컨펌 필요 |
|---|---|---|---|---|---|
| **(A) Azure burst VM** | Standard_D32s_v5 Korea Central | ~$5–12 (6–12시간) | 6–12hr background | Azure Subscription 사용 (이전 K2/Lee burst 과 동일 패턴, 메모리 검증됨) | 🔴 **사용자 컨펌 필수** |
| **(B) Main VM 로컬 설치** | sra-tools + arcasHLA + IMGT/HLA reference | $0 | 12–24hr serial (main VM CPU/disk 점유) | 다른 분석 작업 영향 (main VM 단일) | 🟡 가능하지만 main VM 점유 |
| **(C) Prompt 1 보류** | — | $0 | 0 | 없음 | 🟢 즉시 진행 가능 |

#### (A) Azure burst 세부
- 이전 K2 (n=260) + Lee (n=630) arcasHLA 처리 패턴 검증됨 (메모리 `v17_arcasHLA_korean_k2.md` — $4.80 burst pattern)
- 18 sample × 5M reads paired-end FASTQ download → arcasHLA genotype → genes.json/genotype.json 결과 받기
- ENA HTTPS 차단 이슈 → main VM relay pattern 으로 우회 (이전에 검증)
- 이전 burst 인스턴스 계속 운영 안 한 상태 (auto-shutdown). 다시 spin up 필요.

#### (B) Main VM 로컬
- sra-tools + arcasHLA + IMGTHLA git clone + bowtie2 index 빌드 = 5-10 GB 디스크 + 빌드 시간 + 18 sample × paired FASTQ download (~5-10 GB)
- 직렬 처리: 18 sample × ~30-60 min/sample = 9-18 hr CPU
- main VM 의 다른 ongoing 분석에 영향
- 1회성 task 라 재사용성 낮음 (K2/Lee 는 이미 burst 에서 처리됨)

#### (C) 보류
- D3-D8 schedule 에서 Prompt 1 background 였음. 보류해도 남은 6 prompt 중 paper-changing 효과 큰 것은 **Prompt 2 (TCGA Hashimoto-like cluster identification)** — 5/1 P3 의 PTC+HT 18/18 DM2 paradox 일반화 검증, local TCGA 데이터만 필요
- Prompt 3 (K2 calibration fix), 6 (BCR repertoire), 7 (DM1 sub-cluster) 도 모두 local data
- Prompt 1 은 paper enabler 지만 Pillar 1 의 n=890 → n=908 확장은 marginal (18 sample 추가). Han Chinese forest meta 가 더 큰 가치.

### 2.4 본인 판단 후보

#### 후보 1: (A) Azure burst, 즉시 spin
- Pros: 5/1 master 13.1 명시 task 즉시 closure. Han Chinese forest meta 는 Pillar 1 의 Pan-Asian replication 핵심.
- Cons: ~$5-12 비용. 본인이 ARIA exit 직후라 budget 모니터링 중.
- Risk: 이전 burst pattern 검증되었지만 1회성 spin up + auto-shutdown 관리 + 결과 retrieve 다시 손이 가야 함.

#### 후보 2: (C) 보류, Prompt 2 먼저
- Pros: Local TCGA 데이터만으로 paper-changing 효과 큼. 비용 $0. PTC+HT axis 일반화 입증되면 Cancer paper 의 Discussion + Methods 가 한층 더 단단해짐.
- Cons: Prompt 1 의 Han Chinese forest meta 는 별도 trajectory 로 미뤄짐. Cancer paper 제출까지 별도 time slot 필요.

#### 후보 3: (B) Main VM 로컬, 백그라운드
- Pros: 비용 $0. Prompt 1 closure.
- Cons: 9-18hr 동안 main VM 점유. Prompt 2/3/6/7 병렬 진행 어려움. 1회성 install 재사용성 낮음.

#### 후보 4: 하이브리드 — Prompt 2 즉시 시작 (오늘) + Prompt 1 Azure burst 다음 주 (paper draft 직전)
- Pros: 두 마리 토끼. Prompt 2 의 paper-changing 효과 즉시 검증. Prompt 1 은 Pillar 1 expansion 으로 manuscript draft 시작 직전 closure.
- Cons: Azure burst 는 결국 spin 함 (비용 발생). 시간 분산.

---

## 3. Web Claude 에게 던질 수 있는 추가 질문

### 3.1 결정 자체
- "후보 1/2/3/4 중 5-7년 plan + Cancer paper 제출 timeline + ARIA exit 직후 budget 상황 종합해서 어떤 게 합리적?"
- "Prompt 1 의 Han Chinese forest meta 가 Cell Rep Med / JCI Insight 수준 venue 결정에 정말 필수적인가? 아니면 Pillar 1 의 K2+Lee n=890 만으로도 paper 의 HLA arm 충분한가?"
- "분당 Bundang outreach 응답 (시나리오 A vs B vs C) 에 따라 Prompt 1 의 우선순위가 어떻게 바뀌나?"

### 3.2 Prompt 1 분석 자체
- "GSE286332 18 sample 의 arcasHLA imputation power 분석. 4-digit allele call rate 80% 이상 가능한지 (5M reads partial 기준)."
- "Chen 2018 Han Chinese GD HLA fine-mapping summary stats 의 access 방법 (PMC6161647 supplementary 또는 dbGaP)."
- "Pan-Asian forest meta 에서 phenotype heterogeneity 처리 — GD vs PTC+HT vs PTC 를 random-effects 로 묶을 수 있나, 아니면 sub-stratify?"

### 3.3 D5 schedule 의 Prompt 2 (paper-changing 후보)
- "TCGA-THCA clinical metadata 에 lymphocytic_thyroiditis flag 있는지. GDC API 에서 retrievable 한 항목."
- "GSE286332 PTC+HT signature (top 100 DEGs) 를 TCGA 에 transfer 했을 때 single-sample GSEA (ssGSEA) 가 적절한지, 아니면 단순 Z-mean 더 robust 한지."
- "TCGA Hashimoto-like 환자가 DM2 enriched 라면 그게 진짜 mechanism layer 인지, 아니면 confounder (예: stromal infiltration 일반) 인지 검증 방법."

### 3.4 통합 strategy
- "5-7년 plan 에서 Cancer paper Phase 0 (지금) + Bundang Graves' paper Phase 1 + DIAL audit method paper Phase 2 의 timeline. Prompt 1 closure 가 어느 phase 에 가장 큰 leverage?"
- "ARIA exit 후 시간 70% 확보 — 그 시간을 D3-D8 7-prompt closure 에 쓸지, Cancer paper draft 시작 (5/1 master 13.2 D3-D4 작업) 에 쓸지."
- "본인 사업 / 공기업 / academic path 결정 timeline 과 paper 제출 timeline 의 관계."

---

## 4. 관련 메모리 + 데이터 위치

### 4.1 5/1 deliverables (이미 완료)
- `project/reports/2026_05_01_FULL_RESULTS_dual_track.md` (893 lines, 47KB) — Master document P1-P5 + Task A/B 통합
- `project/reports/2026_05_01_day2_yu_1pager.md` — Yu professor 1-pager
- `project/reports/2026_04_30_P3_GSE286332_critical_result.md` — P3 critical brief
- `project/reports/2026_04_30_taskA_graves_dataset_inventory.md`
- `project/reports/2026_04_30_taskB_8gene_audit.md`
- `project/reports/2026_04_30_dual_track_for_web_claude.md` — 이전 web Claude brief

### 4.2 5/2 deliverables (오늘)
- `project/reports/2026_05_02_D3P4_wang2024_audit.md` — Chen XF 2024 audit
- `project/reports/2026_05_02_D3P5_pdm1_gradient.md` — P_DM1 gradient + mediation

### 4.3 results/ 디렉토리 (P1-P5 + D3 일부)
```
project/results/p1_driver_mrna_audit/      — P1 driver mRNA × mutation
project/results/p2_power_planB/             — P2 power table + Plan B
project/results/p3_gse286332/               — P3 DEG + GSEA + 8-gene + DM + HLA
project/results/p4_pangenome_vs_tiera67/    — P4 ARI + pan-genome ranking
project/results/p5_8gene_vs_hla_autocorr/   — P5 residualization + 2-cohort meta
project/results/d3p5_pdm1_gradient/         — D3-P5 mediation
```

### 4.4 K2/HLA 작업 (이전 4/29)
```
project/results/v17_korean/arcasHLA/         — K2 PRJEB11591 n=260 arcasHLA results
project/results/v17_korean/arcasHLA_GSE213647/  — Lee 2024 n=630 arcasHLA results
reports/html/figs_interactive/v17/hla_mega/  — 30+ HLA figures
reports/html/pages/HLA_FINAL_ALL_IN_ONE.html — All-in-one HLA page
```

### 4.5 핵심 메모리
- `v17_2026_04_30_pivot.md` — 4/30 dual-track decision
- `v17_gse286332_strong_go.md` — P3 STRONG GO (5/1 추가)
- `v17_8gene_pangenome_robustness.md` — P4 robustness (5/1 추가)
- `v17_arcasHLA_korean_k2.md` — Azure burst $4.80 pattern 검증
- `v17_korean_k2_calibration.md` — K2 raw-TPM ≠ centered profile (Prompt 3 의 root cause)
- `v17_K2_vs_bundang_distinction.md` — K2 ≠ Bundang
- `v17_dark_matter_pivot_2026_04_29.md` — dark matter reframe

---

## 5. 한 줄 요약

**Prompt 1 (Azure burst $5-12) vs Prompt 2 (TCGA local, paper-changing) — 어디 먼저 가야 5-7년 plan + Cancer paper 제출 + ARIA exit budget 다 만족?** Web Claude 가 판단해서 후보 1/2/3/4 중 추천 + 그 후 D3-D8 schedule 재구성 부탁.
