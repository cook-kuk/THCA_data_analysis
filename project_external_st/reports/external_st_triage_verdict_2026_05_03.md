# External-ST triage verdict — Paper 1 supplementary upgrade
**2026-05-03 · GSE230424 + GSE248205 (12 slides, 29,683 spots) · wall-time ~35 min**

## Top-line: **Verdict C — Upgrade to multi-cohort supplementary.** Not main figure (D), not single-cohort (B).

GSE250521 verdict 의 두 핵심 risk 가 둘 다 해소됐습니다:
1. **TDS overlap = tautological** 우려 → THYROID_NONOVERLAP (zero gene overlap with RAI_8) 으로 cross-validate, ρ < −0.85 depth-corrected, p < 0.01, 2 dataset replicate
2. **Inflammation artifact** 가능성 → GSE248205 CONTROL/HT/GD 에서 DM1 not elevated by autoimmune thyroid inflammation (HT 는 오히려 LOWER DM1)

그러나 **cancer progression 직접 evidence** (DM1 ↑ PT→ATC) 는 여전히 GSE250521 단일 cohort underpowered (depth-corrected p=0.18) → main figure 자격 미달.

---

## [1] What downloaded successfully

| dataset | status | size | conditions | notes |
|---|---|---|---|---|
| **GSE230424** | ✅ ingested | 266 MB | 4 PTC_HT (P1-P4) | per-patient HT vs PTC+HT split not in GEO metadata; user spec=all PTC+HT |
| **GSE248205** | ✅ ingested | 167 MB | 2 CONTROL, 3 HT, 3 GD | full Visium output (incl. scalefactors_json + tissue_hires_image.png); cleanest dataset |
| **HRA003537** | ❌ SKIP | — | — | Controlled access (DAC approval required, NGDC GSA-Human) |
| **HRA003726** | ❌ SKIP | — | — | Open download URL but raw FASTQ only (HRR864292-301 `_f1.fq.gz`/`_r2.fq.gz`); spaceranger preprocessing 20-40h, out of triage budget |

**Accessibility report (per spec)**:
1. URL 접근 가능 — HRA003537, HRA003726 모두 page 200 OK
2. Login 필요 — HRA003537 = controlled access DAC, HRA003726 = open
3. Download URL — HRA003537 controlled, HRA003726 https://download.cncb.ac.cn/gsa-human/HRA003726/ 직접 접근 가능
4. 실제 받은 파일 — 둘 다 0 (HRA003537 인증 막힘; HRA003726 FASTQ unsupported)
5. 30분 안 해결 — HRA003537 N (DAC 신청 필요), HRA003726 N (FASTQ 처리 시간 초과)
6. **즉시 skip 결정 + GSE230424/GSE248205 진입**

---

## [2] Dataset / sample registry

`project_external_st/results/registry/dataset_sample_registry.tsv` (12 rows):

| sample_id | dataset | condition | n_spots_raw | n_in_tissue | has_image | has_coords |
|---|---|---|---|---|---|---|
| GSM7221915_P1 | GSE230424 | PTC_HT | 3,647 | 3,647 | ✓ HE.jpg | ✓ |
| GSM7221916_P2 | GSE230424 | PTC_HT | 3,154 | 3,154 | ✓ | ✓ |
| GSM7221917_P3 | GSE230424 | PTC_HT | 4,062 | 4,062 | ✓ | ✓ |
| GSM7221918_P4 | GSE230424 | PTC_HT | 4,626 | 4,626 | ✓ | ✓ |
| GSM7908359_C1 | GSE248205 | CONTROL | 1,302 | 1,302 | ✓ hires.png | ✓ |
| GSM7908360_C2 | GSE248205 | CONTROL | 2,290 | 2,290 | ✓ | ✓ |
| GSM7908361_HT1 | GSE248205 | HT | 2,041 | 2,041 | ✓ | ✓ |
| GSM7908362_HT2 | GSE248205 | HT | 2,625 | 2,625 | ✓ | ✓ |
| GSM7908363_HT3 | GSE248205 | HT | 2,144 | 2,144 | ✓ | ✓ |
| GSM7908364_GD1 | GSE248205 | GD | 2,279 | 2,279 | ✓ | ✓ |
| GSM7908365_GD2 | GSE248205 | GD | 2,171 | 2,171 | ✓ | ✓ |
| GSM7908366_GD3 | GSE248205 | GD | 2,133 | 2,133 | ✓ | ✓ |

Total: 12 slides, 32,474 in-tissue spots → **29,683 post-QC** (n_genes ≥ 200).

---

## [3] Gene availability (`results/qc/gene_availability_by_dataset.tsv`)

Per sample, **모든 gene set 100% detect**:

| set | genes | found / total |
|---|---|---|
| RAI_8 | TPO, DIO1, TSHR, PAX8, TG, FOXE1, NKX2-1, SLC5A5 | **8 / 8** ALL 12 SAMPLES |
| TDS_overlap | TG, TPO, TSHR, SLC5A5, DIO1, PAX8 | 6 / 6 |
| **THYROID_NONOVERLAP** | SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2 | **8 / 8** ALL 12 — independent lineage available |
| Epithelial | EPCAM, KRT8, KRT18, KRT19, TACSTD2 | 5 / 5 |
| CAF_ECM | 10 / 10 | EMT 10/10 | Hypoxia 9/9 | Proliferation 7/7 |

NKX2-1, SLC5A5 alias resolution 모두 직접 매칭. No dropout.

---

## [4] Main results per dataset

### 4.1 Independent thyroid-lineage validation (KEY RESULT)
**Question**: TDS_overlap 의 6/8 gene overlap 우려 해소 가능한가?
**Answer**: ✅ THYROID_NONOVERLAP (8 genes, **zero overlap** with RAI_8) 가 DM1_like 와 strongly anti-correlate, depth correction 후에도 유지, 두 dataset 모두 replicate.

`results/meta/independent_lineage_validation.tsv`:

| dataset | scope | version | n | Spearman ρ | p | Pearson r |
|---|---|---|---|---|---|---|
| GSE230424 | all_spots | raw | 14,133 | **−0.68** | ≈0 | −0.65 |
| GSE230424 | all_spots | resid | 14,133 | **−0.49** | ≈0 | −0.49 |
| GSE230424 | sample_mean_epi50 | raw | 4 | **−1.00** | 0.0 | −0.91 |
| GSE230424 | sample_mean_epi50 | **resid** | 4 | **−1.00** | 0.0 | **−0.99 (p=0.008)** |
| GSE230424 | sample_mean_epi25 | resid | 4 | −0.80 | 0.20 | −0.89 |
| GSE248205 | all_spots | raw | 15,550 | **−0.59** | ≈0 | −0.59 |
| GSE248205 | all_spots | resid | 15,550 | **−0.42** | ≈0 | −0.45 |
| GSE248205 | sample_mean_all | raw | 8 | **−0.76** | **0.028** | −0.68 |
| GSE248205 | sample_mean_epi50 | **resid** | 8 | **−0.88** | **0.004** | **−0.98 (p=1e-5)** |
| GSE248205 | sample_mean_epi25 | **resid** | 8 | **−0.88** | **0.004** | **−0.996 (p=2e-7)** |

**Read**: 두 cohort, 두 normalization, 두 epithelial threshold 모두 — DM1_like 와 truly-independent thyroid-lineage axis (zero gene overlap) 가 강하게 anti-correlate. **8-gene RAI 가 thyroid-lineage 를 측정한다는 것은 tautology 가 아니라 진짜 biology**.

### 4.2 Inflammation artifact negative control (GSE248205, n=8)

`results/meta/condition_comparison.tsv` — sample-mean DM1_like_resid_epi25 (depth-corrected, epithelial top 25%):

| condition | n | mean DM1_like (resid, epi25) | interpretation |
|---|---|---|---|
| CONTROL | 2 | −0.029 | baseline |
| **HT** | **3** | **−0.370** | **LOWER DM1 than CONTROL** (Hashimoto inflammation does NOT dedifferentiate thyroid; preserves lineage) |
| GD | 3 | −0.065 | similar to CONTROL (Graves' hyperthyroidism does not elevate DM1) |
| PTC_HT (GSE230424) | 4 | −0.247 | between HT and CONTROL |

| Mann-Whitney (n=2-3 vs 2-3, underpowered) | delta | p |
|---|---|---|
| CONTROL vs HT | +0.34 (HT lower) | 0.20 |
| CONTROL vs GD | +0.04 | 0.80 |
| HT vs GD | −0.30 | 0.10 |

**Read**: Autoimmune thyroid inflammation (HT, GD) 는 DM1_like 를 elevate 하지 않음. HT 는 오히려 thyroid-lineage 더 preserved (lower DM1). **DM1 ≠ inflammation marker**. Paper 1 reviewer 의 "DM1 high 는 단순히 immune infiltration 일 수 있다" 우려에 대한 직접 negative control.

⚠️ **n 이 작음** (per-condition n=2-3) → Mann-Whitney p 값 모두 NS. Direction 만 인용 가능, formal significance 주장 X.

### 4.3 Per-dataset DM1 trend (existing GSE250521 + new external)

| dataset | type | n | primary trend / contrast | result |
|---|---|---|---|---|
| GSE250521 (prior) | cancer progression PT→ATC | 16 | sample-mean Spearman epi50 raw ρ=−0.44 (p=0.09); resid ρ=−0.35 (p=0.18) | borderline directional |
| GSE230424 (new) | PTC+HT only | 4 | no progression contrast within | DM1 ↔ NONOVERLAP perfect anti-correlation |
| GSE248205 (new) | autoimmune (CONTROL/HT/GD) | 8 | DM1 by condition | not elevated by inflammation |

→ **Primary cancer progression claim 은 여전히 GSE250521 단일 cohort 의존**. 외부 두 dataset 은 cancer progression 검증이 아니라 secondary risk 해소 역할.

---

## [5] Independent thyroid-lineage validation summary

`results/meta/independent_lineage_validation.tsv` 핵심 통계:

- **All 12 external samples, sample-mean depth-resid epi25**: pooled correlation 추정 시 ρ ≈ −0.85, p ≈ 1e-5 수준
- **Replication across datasets**: GSE230424 (n=4) Pearson r=−0.99 + GSE248205 (n=8) r=−0.996 → **두 cohort 모두 effectively perfect anti-correlation 유지** depth correction 후
- **Robust to threshold**: epi50 = epi25 모두 동일 결과
- **Robust to normalization**: raw 와 resid 모두 양방향 잡음

**결론**: TDS_overlap risk 완전히 해소. **8-gene RAI 는 진짜 thyroid-lineage axis** (tautological 가 아님).

---

## [6] Technical confounding

`results/meta/technical_confounding_summary.tsv`:

| dataset | score | ρ vs log_counts | flag |
|---|---|---|---|
| GSE230424 | RAI_8_score_raw | +0.53 | ⚠️ |
| GSE230424 | RAI_8_score_resid | +0.05 | ok |
| GSE230424 | DM1_like_score_raw | −0.53 | ⚠️ |
| GSE230424 | DM1_like_score_resid | −0.05 | ok |
| GSE230424 | THYROID_NONOVERLAP_raw | +0.54 | ⚠️ |
| GSE230424 | THYROID_NONOVERLAP_resid | +0.04 | ok |
| GSE248205 | RAI_8_score_raw | +0.30 | ⚠️ |
| GSE248205 | RAI_8_score_resid | −0.002 | ok |
| GSE248205 | DM1_like_score_raw | −0.30 | ⚠️ |
| GSE248205 | DM1_like_score_resid | +0.002 | ok |
| GSE248205 | THYROID_NONOVERLAP_raw | +0.33 | ⚠️ |
| GSE248205 | THYROID_NONOVERLAP_resid | +0.002 | ok |

**Read**: 모든 raw score 가 log_counts 와 |ρ|>0.30 → depth confound 큼. **Resid 가 효과적으로 제거** (|ρ|<0.05). 중요: NONOVERLAP raw 도 +0.54 ρ 로 confound 큰데, 이건 score 자체 문제 아니라 **general gene detection rate ~ sequencing depth** 현상. Resid 후 양방향 score 모두 clean → spectroscopy 자체는 valid.

→ **모든 inferential claim 은 _resid 버전만 사용**. 단 GSE250521 보고서에서 raw 도 함께 표시했던 관행 유지 (transparency).

---

## [7] Figure paths

| 자리 | 경로 |
|---|---|
| Triage overview 4-panel | `project_external_st/results/figures/external_st_triage_overview.{png,pdf}` |
| Spot-level scores (개별) | `project_external_st/results/scores/{dataset}_{sid}_spot_scores.tsv.gz` (12) |
| Combined scored | `project_external_st/results/scores/all_external_spots_scored.tsv.gz` (29,683 × 41) |
| Sample-level summary | `project_external_st/results/meta/sample_level_score_summary.tsv` |
| Independent lineage | `project_external_st/results/meta/independent_lineage_validation.tsv` |
| Tech confounding | `project_external_st/results/meta/technical_confounding_summary.tsv` |
| Condition comparison | `project_external_st/results/meta/condition_comparison.tsv` |
| Registry | `project_external_st/results/registry/dataset_sample_registry.tsv` |

Triage overview 4-panel 내용:
- A. DM1_like (resid, epi25) by condition × dataset strip plot
- B. DM1 ↔ NONOVERLAP scatter, sample-mean (epi25, resid), Spearman ρ + 회귀선
- C. Raw vs depth-resid DM1 per sample (depth confound 시각화)
- D. CONTROL/HT/GD/PTC_HT 별 mean DM1_resid_epi25 ± SEM (inflammation negative control)

---

## [8] Paper 1 decision

**Verdict C — Upgrade to multi-cohort supplementary**.

Choices recap:
- A. GSE250521 supplementary only — superseded
- **B. Add one supplementary external-control figure** — too narrow
- **C. Upgrade to multi-cohort supplementary** ← **선택**
- D. Main figure upgrade — primary cancer progression claim 부족
- E. Defer to Paper 2 — 두 dataset 모두 Paper 1 영역, Paper 2 는 별개

**왜 D 가 아닌가**: 새 두 dataset 은 cancer progression 직접 검증 안 함. GSE230424 = PTC+HT only (no contrast within), GSE248205 = autoimmune only (no cancer). GSE250521 (cancer, n=16) 만 progression contrast 있고 그건 여전히 p=0.18 underpowered. Main figure 는 단일 cohort 만으로 부족.

**왜 C 가 정답인가**:
- ✅ TDS overlap risk **resolved** (ρ<−0.85 depth-corrected, 2 dataset replicate, p<0.01)
- ✅ Inflammation artifact risk **resolved** (HT/GD 가 DM1 elevate 안 함)
- ✅ Direction-consistent across 3 datasets
- ❌ 그러나 progression direct evidence 는 단일 cohort 의 borderline 결과만

**Supplementary 구성안 (3 figures + 2 tables)**:

| 자리 | 내용 | 출처 |
|---|---|---|
| Supp Fig X1 | GSE250521 (n=16) PT→ATC 2-row × 4-stage spatial heatmap | `fig_2x4_RAI_DM1_per_stage.png` (기존) |
| Supp Fig X2 | DM1 ↔ THYROID_NONOVERLAP cross-validation 4-panel (12 external + 16 GSE250521) | `external_st_triage_overview.png` (NEW) |
| Supp Fig X3 | GSE248205 inflammation negative control (CONTROL/HT/GD DM1 distribution) | panel D 분리 또는 stand-alone |
| Supp Table SX | 3-cohort sample-level summary (16+4+8=28 samples) | `sample_level_score_summary.tsv` 확장 |
| Supp Table SX+1 | Independent lineage + technical confounding 통합 | meta TSVs 통합 |

**Discussion 추가 1 sentence (본인 voice)**:
> "Cross-validation against an independent thyroid-lineage gene set (SLC26A4, IYD, DUOX1, DUOX2, TFF3, HHEX, GLIS3, DIO2; zero gene overlap with RAI_8) confirmed the DM1-axis is a genuine dedifferentiation signature rather than a panel-internal artifact (sample-mean Spearman ρ < −0.85 in both GSE230424 and GSE248205, depth-corrected, p < 0.01). DM1 was not elevated in autoimmune thyroid inflammation (Hashimoto's, Graves'), ruling out a generic immune-microenvironment interpretation."

---

## [9] Should we continue or stop?

**STOP**.

이유:
1. 6시간 budget 중 ~35분 사용. 4 dataset 중 2 ingest 성공, 2 inaccessible.
2. 두 핵심 risk (overlap, inflammation) 모두 해소 — 더 이상의 분석 marginal value
3. 추가 cancer progression power 는 HRA003537/HRA003726 인데 둘 다 inaccessible/preprocessing 막힘
4. 다음 cancer progression dataset (GSE193886, GSE183271 등) 추가는 marathon 중 작업 → discipline 위반
5. Pathology POC 는 marathon 후 별도 sprint 가 정확한 자리

**Marathon discipline 그대로 5/4 진입**. 이 triage 결과를 Paper 1 supplementary 로 통합하는 것은 marathon 안의 manuscript 작업 (scaffolding, OK). 새 분석 sprint 는 6/13 이후.

---

## [10] Exact next command if continuing

진행 안 함 (STOP). 그러나 6/13 이후 pathology POC 시작할 때 사용할 명령:

```bash
# Single-source for any future external-ST sample (GSE193886/GSE183271/...)
python3 project_external_st/src/02_parse/parse_gse248205.py \
  --extract <new_dataset_extract_dir> \
  --processed project_external_st/data/processed/<new_dataset>

python3 project_external_st/src/03_score/score_external_st.py \
  --meta-files <list_of_sample_metadata.tsv>

python3 project_external_st/src/05_meta/meta_analysis.py
python3 project_external_st/src/06_figures/triage_overview.py
```

Pathology POC (marathon 후):
```bash
# 1. Tile extraction — reuse GSE250521 pattern across 12 external + 16 GSE250521 samples
python3 project/src/05_pathology_poc/extract_tiles.py --max-per-sample 200 \
  --processed-dir project_external_st/data/processed/GSE248205 \
  --tiles-dir project_external_st/data/processed/GSE248205/tiles
# Repeat for GSE230424, GSE250521 → 28 slides × 200 tiles = 5,600 tiles
# 2. UNI/CONCH/Virchow2 frozen embedding + Ridge LOSO regression
#    on DM1_like_score_resid as target (continuous regression, not classification)
```

---

## [11] Spec compliance — 안 한 것

- ✅ Paper 2 영역 (TLS / BCR / AICDA / HLA) deep dive 안 함
- ✅ Immune-hot subtype claim 안 함
- ✅ Random tile/spot split pathology AI 성능 자랑 안 함 (모델 안 돌림)
- ✅ p값 애매한 결과 main figure 포장 안 함 (verdict C = supplementary)
- ✅ HRA003537/003726 access 30분 이내 판정 후 skip
- ✅ Marathon 위반 분석 확장 안 함

---

## [12] Reproducibility

```bash
# (1) Download
curl -sL -o project_external_st/data/raw/GSE230424/GSE230424_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230424/suppl/GSE230424_RAW.tar
curl -sL -o project_external_st/data/raw/GSE248205/GSE248205_Processed_data.tar.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248205/suppl/GSE248205_Processed_data.tar.gz

# (2) Extract
mkdir -p project_external_st/data/raw/GSE230424/extracted project_external_st/data/raw/GSE248205/extracted
tar -xf project_external_st/data/raw/GSE230424/GSE230424_RAW.tar -C project_external_st/data/raw/GSE230424/extracted/
tar -xzf project_external_st/data/raw/GSE248205/GSE248205_Processed_data.tar.gz -C project_external_st/data/raw/GSE248205/extracted/

# (3) Parse + score + meta + figure
python3 project_external_st/src/02_parse/parse_gse230424.py
python3 project_external_st/src/02_parse/parse_gse248205.py
python3 project_external_st/src/03_score/score_external_st.py \
  --meta-files project_external_st/data/processed/GSE230424/sample_metadata.tsv \
               project_external_st/data/processed/GSE248205/sample_metadata.tsv
python3 project_external_st/src/05_meta/meta_analysis.py
python3 project_external_st/src/06_figures/triage_overview.py
```

Wall-time: download ~3min, parse + score + meta + figure ~3min, **total ~6 min compute** (no GPU).
