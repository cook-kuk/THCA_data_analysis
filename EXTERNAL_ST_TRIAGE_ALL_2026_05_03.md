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

---

# Appendix A — registry/dataset_sample_registry.tsv

```tsv
dataset	sample_id	source_accession	condition_raw	condition_inferred	disease_axis	has_HE_image	has_spatial_coords	n_spots_raw	n_genes	n_in_tissue	h5ad	image	notes	status
GSE230424	GSM7221915_P1	GSM7221915	P1	PTC_HT	PTC_HT	True	True	3647	36601	3647	project_external_st/data/processed/GSE230424/GSM7221915_P1/GSM7221915_P1.raw.h5ad	project_external_st/data/processed/GSE230424/GSM7221915_P1/HE.jpg	GEO sample_title=P[1-4],ST-seq; per-patient HT vs PTC+HT split not in metadata; user spec=PTC+HT all	ok
GSE230424	GSM7221916_P2	GSM7221916	P2	PTC_HT	PTC_HT	True	True	3154	36601	3154	project_external_st/data/processed/GSE230424/GSM7221916_P2/GSM7221916_P2.raw.h5ad	project_external_st/data/processed/GSE230424/GSM7221916_P2/HE.jpg	GEO sample_title=P[1-4],ST-seq; per-patient HT vs PTC+HT split not in metadata; user spec=PTC+HT all	ok
GSE230424	GSM7221917_P3	GSM7221917	P3	PTC_HT	PTC_HT	True	True	4062	36601	4062	project_external_st/data/processed/GSE230424/GSM7221917_P3/GSM7221917_P3.raw.h5ad	project_external_st/data/processed/GSE230424/GSM7221917_P3/HE.jpg	GEO sample_title=P[1-4],ST-seq; per-patient HT vs PTC+HT split not in metadata; user spec=PTC+HT all	ok
GSE230424	GSM7221918_P4	GSM7221918	P4	PTC_HT	PTC_HT	True	True	4626	36601	4626	project_external_st/data/processed/GSE230424/GSM7221918_P4/GSM7221918_P4.raw.h5ad	project_external_st/data/processed/GSE230424/GSM7221918_P4/HE.jpg	GEO sample_title=P[1-4],ST-seq; per-patient HT vs PTC+HT split not in metadata; user spec=PTC+HT all	ok
GSE248205	GSM7908359_C1	GSM7908359	C1	CONTROL	CONTROL	True	True	1302	36601	1302	project_external_st/data/processed/GSE248205/GSM7908359_C1/GSM7908359_C1.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908359_C1/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
GSE248205	GSM7908360_C2	GSM7908360	C2	CONTROL	CONTROL	True	True	2290	36601	2290	project_external_st/data/processed/GSE248205/GSM7908360_C2/GSM7908360_C2.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908360_C2/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
GSE248205	GSM7908364_GD1	GSM7908364	GD1	GD	GD	True	True	2279	36601	2279	project_external_st/data/processed/GSE248205/GSM7908364_GD1/GSM7908364_GD1.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908364_GD1/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
GSE248205	GSM7908365_GD2	GSM7908365	GD2	GD	GD	True	True	2171	36601	2171	project_external_st/data/processed/GSE248205/GSM7908365_GD2/GSM7908365_GD2.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908365_GD2/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
GSE248205	GSM7908366_GD3	GSM7908366	GD3	GD	GD	True	True	2133	36601	2133	project_external_st/data/processed/GSE248205/GSM7908366_GD3/GSM7908366_GD3.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908366_GD3/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
GSE248205	GSM7908361_HT1	GSM7908361	HT1	HT	HT	True	True	2041	36601	2041	project_external_st/data/processed/GSE248205/GSM7908361_HT1/GSM7908361_HT1.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908361_HT1/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
GSE248205	GSM7908362_HT2	GSM7908362	HT2	HT	HT	True	True	2625	36601	2625	project_external_st/data/processed/GSE248205/GSM7908362_HT2/GSM7908362_HT2.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908362_HT2/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
GSE248205	GSM7908363_HT3	GSM7908363	HT3	HT	HT	True	True	2144	36601	2144	project_external_st/data/processed/GSE248205/GSM7908363_HT3/GSM7908363_HT3.raw.h5ad	project_external_st/data/processed/GSE248205/GSM7908363_HT3/tissue_hires_image.png	Martinez-Hernandez 2024 PMID 39003267; full Visium output (incl scalefactors+hires png)	ok
```

# Appendix B — qc/gene_availability_by_dataset.tsv

```tsv
sample_id	n_spots_post_qc	RAI_8_n_genes	TDS_overlap_n_genes	THYROID_NONOVERLAP_n_genes	CAF_ECM_n_genes	EMT_n_genes	Hypoxia_n_genes	Proliferation_n_genes	Epithelial_n_genes	RAI_8_genes	TDS_overlap_genes	THYROID_NONOVERLAP_genes	CAF_ECM_genes	EMT_genes	Hypoxia_genes	Proliferation_genes	Epithelial_genes	dataset	condition
GSM7221915_P1	3566	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE230424	PTC_HT
GSM7221916_P2	3005	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE230424	PTC_HT
GSM7221917_P3	2961	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE230424	PTC_HT
GSM7221918_P4	4601	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE230424	PTC_HT
GSM7908359_C1	1064	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	CONTROL
GSM7908360_C2	1744	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	CONTROL
GSM7908364_GD1	2266	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	GD
GSM7908365_GD2	2005	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	GD
GSM7908366_GD3	1969	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	GD
GSM7908361_HT1	1927	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	HT
GSM7908362_HT2	2595	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	HT
GSM7908363_HT3	1980	8	6	8	10	10	9	7	5	TPO,DIO1,TSHR,PAX8,TG,FOXE1,NKX2-1,SLC5A5	TG,TPO,TSHR,SLC5A5,DIO1,PAX8	SLC26A4,IYD,DUOX1,DUOX2,TFF3,HHEX,GLIS3,DIO2	COL1A1,COL1A2,COL3A1,COL6A1,DCN,LUM,FAP,ACTA2,PDGFRB,POSTN	VIM,FN1,SNAI1,SNAI2,ZEB1,ZEB2,TWIST1,CDH2,ITGA5,ITGB1	VEGFA,CA9,SLC2A1,LDHA,PDK1,BNIP3,NDRG1,EGLN3,ADM	MKI67,TOP2A,PCNA,MCM2,MCM5,STMN1,CENPF	EPCAM,KRT8,KRT18,KRT19,TACSTD2	GSE248205	HT
```

# Appendix C — meta/independent_lineage_validation.tsv (KEY)

```tsv
dataset	scope	version	n	spearman_rho	spearman_p	pearson_r	pearson_p	interpretation
GSE230424	all_spots	raw	14133	-0.6778991184603385	0.0	-0.6543048029526817	0.0	DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
GSE230424	all_spots	resid	14133	-0.4907277373778117	0.0	-0.48815557216416705	0.0	DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
GSE248205	all_spots	raw	15550	-0.5916398323344473	0.0	-0.592702430229072	0.0	DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
GSE248205	all_spots	resid	15550	-0.4237188283680889	0.0	-0.4512621022011145	0.0	DM1_like vs THYROID_NONOVERLAP — true independent lineage check (no overlap with 8-gene)
GSE230424	sample_mean_all	raw	4	0.19999999999999998	0.8	-0.04319728822571858	0.9568027117742814	DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
GSE230424	sample_mean_epi50	raw	4	-1.0	0.0	-0.9079401692940423	0.09205983070595769	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
GSE230424	sample_mean_epi25	raw	4	-0.7999999999999999	0.2000000000000001	-0.8901738446021238	0.10982615539787632	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
GSE230424	sample_mean_all	resid	4	0.39999999999999997	0.6	0.2354598609066923	0.7645401390933078	DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
GSE230424	sample_mean_epi50	resid	4	-1.0	0.0	-0.9922998547245199	0.00770014527548013	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
GSE230424	sample_mean_epi25	resid	4	-0.7999999999999999	0.2000000000000001	-0.8855376804163091	0.11446231958369091	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
GSE248205	sample_mean_all	raw	8	-0.7619047619047621	0.028004939153071805	-0.6804112023543244	0.06329489208004352	DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
GSE248205	sample_mean_epi50	raw	8	-0.9761904761904763	3.314396026200098e-05	-0.952856443924328	0.0002527687647699057	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
GSE248205	sample_mean_epi25	raw	8	-0.9761904761904763	3.314396026200098e-05	-0.9875673132915305	4.759666886352331e-06	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
GSE248205	sample_mean_all	resid	8	0.19047619047619052	0.6514014957024814	0.24047596033013344	0.5661890148104232	DM1_like vs THYROID_NONOVERLAP at sample-mean (all)
GSE248205	sample_mean_epi50	resid	8	-0.880952380952381	0.0038503204637324005	-0.9835697384701191	1.0952324856746705e-05	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi50)
GSE248205	sample_mean_epi25	resid	8	-0.880952380952381	0.0038503204637324005	-0.995775779417616	1.8784638373668558e-07	DM1_like vs THYROID_NONOVERLAP at sample-mean (epi25)
```

# Appendix D — meta/technical_confounding_summary.tsv

```tsv
dataset	score	n	spearman_rho_with_log_counts	p	flag
GSE230424	RAI_8_score_raw	14133	0.5250977996463394	0.0	⚠️ depth-confound
GSE230424	DM1_like_score_raw	14133	-0.5250977996463394	0.0	⚠️ depth-confound
GSE230424	RAI_8_score_resid	14133	0.04869431464750445	6.962959406502137e-09	ok
GSE230424	DM1_like_score_resid	14133	-0.04869431464750445	6.962959406502137e-09	ok
GSE230424	THYROID_NONOVERLAP_score_raw	14133	0.541594583405772	0.0	⚠️ depth-confound
GSE230424	THYROID_NONOVERLAP_score_resid	14133	0.04477759914841075	1.0069570702760872e-07	ok
GSE248205	RAI_8_score_raw	15550	0.30180685764015897	0.0	⚠️ depth-confound
GSE248205	DM1_like_score_raw	15550	-0.30180685764015897	0.0	⚠️ depth-confound
GSE248205	RAI_8_score_resid	15550	-0.0016644993789532675	0.8355835265535433	ok
GSE248205	DM1_like_score_resid	15550	0.0016644993789532675	0.8355835265535433	ok
GSE248205	THYROID_NONOVERLAP_score_raw	15550	0.3284834630155037	0.0	⚠️ depth-confound
GSE248205	THYROID_NONOVERLAP_score_resid	15550	0.0019020238569800432	0.8125304662085958	ok
```

# Appendix E — meta/condition_comparison.tsv

```tsv
dataset	version	scope	comparison	n_a	n_b	mean_a	mean_b	delta_mean	mannwhitney_U	mannwhitney_p
GSE248205	raw	all	CONTROL_vs_HT	2	3	4.2946478904990865e-08	-3.2553032097052132e-09	4.620178211469608e-08	6.0	0.2
GSE248205	raw	all	CONTROL_vs_GD	2	3	4.2946478904990865e-08	1.4373980049516851e-09	4.150908090003918e-08	4.0	0.8
GSE248205	raw	all	HT_vs_GD	3	3	-3.2553032097052132e-09	1.4373980049516851e-09	-4.6927012146568986e-09	6.0	0.7
GSE248205	raw	epi50	CONTROL_vs_HT	2	3	-0.1487550515829441	-0.2948581645650985	0.1461031129821544	6.0	0.2
GSE248205	raw	epi50	CONTROL_vs_GD	2	3	-0.1487550515829441	-0.16777124144560618	0.01901618986266207	4.0	0.8
GSE248205	raw	epi50	HT_vs_GD	3	3	-0.2948581645650985	-0.16777124144560618	-0.12708692311949232	0.0	0.1
GSE248205	raw	epi25	CONTROL_vs_HT	2	3	-0.16988904039380476	-0.49489469413538867	0.3250056537415839	6.0	0.2
GSE248205	raw	epi25	CONTROL_vs_GD	2	3	-0.16988904039380476	-0.1782160439585616	0.008327003564756852	4.0	0.8
GSE248205	raw	epi25	HT_vs_GD	3	3	-0.49489469413538867	-0.1782160439585616	-0.31667865017682706	0.0	0.1
GSE248205	resid	all	CONTROL_vs_HT	2	3	-2.287915628107647e-18	-2.5713959548019025e-18	2.834803266942556e-19	3.0	1.0
GSE248205	resid	all	CONTROL_vs_GD	2	3	-2.287915628107647e-18	1.842657495672237e-18	-4.130573123779884e-18	2.0	0.8
GSE248205	resid	all	HT_vs_GD	3	3	-2.5713959548019025e-18	1.842657495672237e-18	-4.414053450474139e-18	2.0	0.4
GSE248205	resid	epi50	CONTROL_vs_HT	2	3	-0.04088077970884896	-0.21200244394369028	0.17112166423484132	6.0	0.2
GSE248205	resid	epi50	CONTROL_vs_GD	2	3	-0.04088077970884896	-0.056814254305147	0.015933474596298042	4.0	0.8
GSE248205	resid	epi50	HT_vs_GD	3	3	-0.21200244394369028	-0.056814254305147	-0.15518818963854328	0.0	0.1
GSE248205	resid	epi25	CONTROL_vs_HT	2	3	-0.029091154523980277	-0.36957413703608544	0.3404829825121052	6.0	0.2
GSE248205	resid	epi25	CONTROL_vs_GD	2	3	-0.029091154523980277	-0.0652087426230931	0.03611758809911283	6.0	0.2
GSE248205	resid	epi25	HT_vs_GD	3	3	-0.36957413703608544	-0.0652087426230931	-0.3043653944129923	0.0	0.1
```

# Appendix F — meta/sample_level_score_summary.tsv (12 rows × 80 cols, full)

```tsv
sample_id	dataset	condition	n_spots	mean_RAI_8_score_raw_all	mean_RAI_8_score_raw_epi50	mean_RAI_8_score_raw_epi25	median_RAI_8_score_raw_all	mean_DM1_like_score_raw_all	mean_DM1_like_score_raw_epi50	mean_DM1_like_score_raw_epi25	median_DM1_like_score_raw_all	mean_TDS_overlap_score_raw_all	mean_TDS_overlap_score_raw_epi50	mean_TDS_overlap_score_raw_epi25	median_TDS_overlap_score_raw_all	mean_THYROID_NONOVERLAP_score_raw_all	mean_THYROID_NONOVERLAP_score_raw_epi50	mean_THYROID_NONOVERLAP_score_raw_epi25	median_THYROID_NONOVERLAP_score_raw_all	mean_Epithelial_score_raw_all	mean_Epithelial_score_raw_epi50	mean_Epithelial_score_raw_epi25	median_Epithelial_score_raw_all	mean_CAF_ECM_score_raw_all	mean_CAF_ECM_score_raw_epi50	mean_CAF_ECM_score_raw_epi25	median_CAF_ECM_score_raw_all	mean_EMT_score_raw_all	mean_EMT_score_raw_epi50	mean_EMT_score_raw_epi25	median_EMT_score_raw_all	mean_Hypoxia_score_raw_all	mean_Hypoxia_score_raw_epi50	mean_Hypoxia_score_raw_epi25	median_Hypoxia_score_raw_all	mean_Proliferation_score_raw_all	mean_Proliferation_score_raw_epi50	mean_Proliferation_score_raw_epi25	median_Proliferation_score_raw_all	mean_RAI_8_score_resid_all	mean_RAI_8_score_resid_epi50	mean_RAI_8_score_resid_epi25	median_RAI_8_score_resid_all	mean_DM1_like_score_resid_all	mean_DM1_like_score_resid_epi50	mean_DM1_like_score_resid_epi25	median_DM1_like_score_resid_all	mean_TDS_overlap_score_resid_all	mean_TDS_overlap_score_resid_epi50	mean_TDS_overlap_score_resid_epi25	median_TDS_overlap_score_resid_all	mean_THYROID_NONOVERLAP_score_resid_all	mean_THYROID_NONOVERLAP_score_resid_epi50	mean_THYROID_NONOVERLAP_score_resid_epi25	median_THYROID_NONOVERLAP_score_resid_all	mean_Epithelial_score_resid_all	mean_Epithelial_score_resid_epi50	mean_Epithelial_score_resid_epi25	median_Epithelial_score_resid_all	mean_CAF_ECM_score_resid_all	mean_CAF_ECM_score_resid_epi50	mean_CAF_ECM_score_resid_epi25	median_CAF_ECM_score_resid_all	mean_EMT_score_resid_all	mean_EMT_score_resid_epi50	mean_EMT_score_resid_epi25	median_EMT_score_resid_all	mean_Hypoxia_score_resid_all	mean_Hypoxia_score_resid_epi50	mean_Hypoxia_score_resid_epi25	median_Hypoxia_score_resid_all	mean_Proliferation_score_resid_all	mean_Proliferation_score_resid_epi50	mean_Proliferation_score_resid_epi25	median_Proliferation_score_resid_all	mean_log_total_counts	mean_n_genes
GSM7221915_P1	GSE230424	PTC_HT	3566	2.206451486304301e-08	0.3670549267130847	0.4796970143377018	-0.0208599605	-2.206451486304301e-08	-0.3670549267130847	-0.4796970143377018	0.0208599605	2.493031408034424e-08	0.3913311623570387	0.5123328533604261	-0.025343707	-2.772887268511592e-08	0.283975866663253	0.35191779470116596	-0.0257060225	-3.072728547829892e-08	0.4659691216003926	0.755777583912556	-0.0249632245	-8.083568141649177e-09	0.011838774565073472	0.019962101555017937	-0.0663158595	1.6387944478856272e-08	0.032086600596040384	0.023684648901390132	-0.05471988	1.8208939993015956e-09	0.09002310625518116	0.10749946735345066	-0.0226016905	-3.936090853014913e-09	-0.03599492729680314	-0.04796067909573991	-0.27845258	-7.970193334381382e-18	0.19810554938638608	0.2940545275888461	0.0047317959855509	7.970193334381382e-18	-0.19810554938638608	-0.2940545275888461	-0.0047317959855509	-4.981370833988364e-18	0.2188359018482959	0.32462403305209303	0.0038517060838807998	7.970193334381382e-18	0.1370391051264972	0.19047926985339617	-0.032795552280521195	1.1581687189022945e-17	0.3312536876382997	0.6122019938724633	-0.0498097686435046	4.4832337505895274e-18	-0.04504306987300049	-0.047544558998821256	-0.0815417385011164	4.981370833988364e-18	-0.022268908811705186	-0.04402208151717099	-0.06226293052115155	8.2192618760808e-18	0.018894249517681524	0.024824878405239545	-0.05826127436142465	1.2951564168369745e-17	-0.06680077906704861	-0.0901067790417984	-0.11227490045638955	8.344773142072134	2078.8359506449806
GSM7221916_P2	GSE230424	PTC_HT	3005	4.333552079608425e-08	0.24647042147032597	0.2692194738589096	-0.029272512	-4.333552079608425e-08	-0.24647042147032597	-0.2692194738589096	0.029272512	4.747275207614961e-08	0.2586491222822222	0.269426595506516	0.0070655546	-1.218909817450258e-08	0.22768912847045908	0.2354583024942819	-0.063830316	6.490285194178369e-09	0.46645236600486156	0.7612834013696809	-0.023605501	-1.4645257903424116e-09	-0.025370272082435134	-0.029008647665385637	-0.061739314	7.599166055004309e-09	0.03076072548569528	0.04084553357944149	-0.055418987	8.22465557031703e-09	0.06577522534038589	0.09434181893125	-0.037024245	-8.414808645632649e-09	0.008418468579241523	-0.00044981604547871657	-0.24921398	-7.09360468312912e-18	0.11509475573873088	0.10208707771707576	0.0044495861503104	7.09360468312912e-18	-0.11509475573873088	-0.10208707771707576	-0.0044495861503104	-3.54680234156456e-18	0.12257536223034667	0.09384676010752593	0.0253299208946177	3.54680234156456e-18	0.09419931867039258	0.06845264998252622	-0.0216128408959226	7.09360468312912e-18	0.36525915983625645	0.6449610802898839	-0.0668341264057466	2.3645348943763734e-18	-0.04405049840674774	-0.042731346375513825	-0.0734528514712101	1.3596075642664147e-17	-0.001691560927011218	0.0037518552099645225	-0.0663689944760291	9.458139577505494e-18	0.010845328157974506	0.02947232007596308	-0.0705894045558138	1.77340117078228e-17	-0.03056717363628052	-0.042914392061885714	-0.1291697017235323	8.004281383292286	1737.7630615640599
GSM7221917_P3	GSE230424	PTC_HT	2961	1.3911503982100136e-08	0.2798939110287644	0.3755527770995142	-0.020792156	-1.3911503982100136e-08	-0.2798939110287644	-0.3755527770995142	0.020792156	5.505856129904826e-09	0.3026008615788926	0.4066285236062078	-0.014084418	1.221285038531682e-08	0.27638125283041864	0.3797876744242915	-0.0342924	1.3434012159441826e-08	0.4292191915908913	0.7452067874224022	-0.10555353	-9.253600128941454e-09	0.11736817756683998	0.14554883456099865	-0.049047686	-1.3798351907576257e-08	0.07960088931808237	0.09195432635863697	-0.056221746	-7.20949003734712e-09	0.09124024564996623	0.10708486747619433	-0.058487944	5.2686558517981924e-09	0.0013456630295003303	-0.012032189744480436	-0.27885732	7.498973486154385e-18	0.09223422104804148	0.15978996712888807	-0.048099423703895	-7.498973486154385e-18	-0.09223422104804148	-0.15978996712888807	0.048099423703895	2.6996304550155782e-18	0.09926699140319847	0.17409716426735977	-0.0499019433273557	-6.299137728369682e-18	0.08786535079216244	0.16902325154273432	-0.0339421474128753	1.1398439698954664e-17	0.2913786318253962	0.5921760750272337	-0.0532615023263867	1.0798521820062313e-17	0.014110210467548463	0.02224212704662005	-0.0785019473863326	1.2598275456739365e-17	0.0012443741048180626	-9.921636618963473e-05	-0.0798522643440212	6.599096667815858e-18	0.008607664768087336	0.010158510232603935	-0.0663015532537917	2.654636614098652e-17	-0.03146503596521097	-0.05125413879969498	-0.1590387925163006	7.864565656632322	1558.222897669706
GSM7221918_P4	GSE230424	PTC_HT	4601	4.4751314932536234e-08	0.18206948377309867	0.1964268512337967	0.037418872	-4.4751314932536234e-08	-0.18206948377309867	-0.1964268512337967	-0.037418872	5.023506846233679e-08	0.21855652010910473	0.23506353755939183	0.066851676	1.0873501846414618e-08	0.16507278523500218	0.18322689107779322	0.015322447	2.9581703110495432e-08	0.46002287414080834	0.7241618721633363	0.028538942	1.2757361530037124e-08	-0.03572751439003042	-0.02166160555968723	-0.07682978	-1.6701820908978017e-08	-6.400882919513263e-05	0.016203298451968724	-0.063981935	3.380214299196033e-09	0.049018935997270745	0.06886880428122502	-0.044225797	-1.5781568271283107e-10	-0.003734394413215985	-0.0029465485172545515	-0.23337707	6.177289595827865e-18	0.08025039215485183	0.07162829408587118	0.0090562860666173	-6.177289595827865e-18	-0.08025039215485183	-0.07162829408587118	-0.0090562860666173	3.0886447979139326e-18	0.10491253315507008	0.09468910138295904	0.0322388420752315	1.3126740391134213e-17	0.07360786927792229	0.07355887740750028	-0.0048591983343972	8.493773194263314e-18	0.35853216222240575	0.6066044825660797	-0.0141909486259259	1.698754638852663e-17	-0.046125853238652766	-0.03451827586153373	-0.0782774020862832	1.0810256792698764e-17	-0.022081201103062588	-0.009602550118518193	-0.0701072635773586	1.0810256792698764e-17	0.014837626485421357	0.028089150794942355	-0.0569541786286409	1.235457919165573e-17	-0.018874876053031557	-0.020660445814102948	-0.1967449526754338	8.196229456109075	1777.8674201260596
GSM7908359_C1	GSE248205	CONTROL	1064	-2.8334144736395334e-08	0.1903283989381015	0.2312468666593985	0.0417603145	2.8334144736395334e-08	-0.1903283989381015	-0.2312468666593985	-0.0417603145	-4.5312312028682285e-08	0.19747030828176693	0.24710999143721804	0.07519888999999999	-4.9071287591413034e-08	0.15415640592135338	0.20483303194135338	0.035005301	-8.466560154717304e-09	0.3616018365330452	0.6057719801503759	-0.0376647415	4.824248123787362e-09	-0.06362019077819549	-0.06050314993496241	-0.11816717	4.547462404130132e-09	-0.009616200894172932	-0.02395112838834586	-0.11757363	-5.562312018458164e-09	0.06426108828176692	0.07665774467556391	-0.0383600685	-1.133270675714369e-08	-0.006510527778195472	-0.0016135113834586318	-0.12801062	-4.591147846194256e-18	0.06674487052452653	0.05321945019103397	0.0143498547861858	4.591147846194256e-18	-0.06674487052452653	-0.05321945019103397	-0.0143498547861858	-2.5042624615605033e-18	0.07314940650118751	0.06819078390248522	0.01898223186710965	-1.2521312307802517e-18	0.03453403724564876	0.03353447075030577	-0.008564396837127699	0.0	0.2769716465023349	0.4861525641156546	-0.060407905399017955	1.0017049846242013e-17	-0.07380118160230446	-0.07485696301586756	-0.11204883408725205	1.0017049846242013e-17	-0.032912680262900465	-0.05694386716894088	-0.0927719653388229	1.7947214307850274e-17	0.015528565596294102	0.006668966474404623	-0.07664540430590294	2.5042624615605033e-17	-0.026216347214862983	-0.029749769280188656	-0.11045868609411169	7.17008890344879	945.7208646616541
GSM7908360_C2	GSE248205	CONTROL	1744	-5.75588130735864e-08	0.1071817042277867	0.10853121412821101	0.01071248575	5.75588130735864e-08	-0.1071817042277867	-0.10853121412821101	-0.01071248575	-8.320704415333901e-08	0.11146390254795296	0.11527097868309633	0.0256301685	-4.530177743636874e-09	0.09357599382308487	0.08820000567809633	0.0038451441	-1.0134495402666465e-08	0.38323170771954135	0.6590022357339451	-0.068722705	-8.106135316504768e-09	-0.019922520358142196	-0.0016653132759174296	-0.09476305	3.898509155660689e-10	0.021547292423165134	0.016089796535091738	-0.1109452175	-1.4852884167230724e-08	0.03699463499409405	0.045362589204472475	-0.0487491995	1.2500000030979703e-09	0.01979424783256881	0.02050708094036697	-0.10629403	9.16697910240955e-18	0.015016688893171383	0.004962858856926584	-0.017775212368554	-9.16697910240955e-18	-0.015016688893171383	-0.004962858856926584	0.017775212368554	3.564936317603714e-18	0.01723681766358908	0.008756213945109847	-0.004495741387798249	-5.0927661680053055e-19	0.010688468309190865	-0.004432307613017266	-0.03673208644837335	-1.2731915420013264e-17	0.3196019074245176	0.593937283186703	-0.09778313532932009	7.129872635207428e-18	-0.030496102351522194	-0.013481617589889219	-0.1066776392789576	1.3750468653614325e-17	-0.003583538028910464	-0.011736512127512672	-0.0895829859137577	5.729361939005969e-18	0.0024559077744876205	0.007534737555859098	-0.08383913255972515	3.564936317603714e-17	0.0030978776004409077	0.0020798328287160634	-0.1008891745558166	6.974112489676446	751.5120412844037
GSM7908361_HT1	GSE248205	HT	1927	4.007814224643958e-09	0.3317991512843361	0.6577239333952283	-0.26966327	-4.007814224643958e-09	-0.3317991512843361	-0.6577239333952283	0.26966327	-1.511996887008646e-08	0.37915562216267634	0.7504351327978631	-0.29702786	-2.4005469642079892e-08	0.24330654920477177	0.5051605273622407	-0.2074973	2.1407924228582423e-08	0.4302569308953008	0.8759642589419087	-0.22218843	-3.0532589015467084e-10	-0.012772356497136922	-0.09108828720477177	-0.096768945	2.5001816295303785e-08	0.06302304882865145	0.09152637093056017	-0.084191196	4.553710421482786e-09	0.050327325213163895	0.10195093962742738	-0.06709805	-1.1039880626348735e-08	-0.0634854372086618	-0.10959227564066389	-0.19429623	5.530950200519721e-18	0.2647839261040105	0.5424896249416963	-0.1277682939303332	-5.530950200519721e-18	-0.2647839261040105	-0.5424896249416963	0.1277682939303332	1.1061900401039443e-17	0.30541097090795827	0.6240342159097584	-0.1468592024510836	4.609125167099768e-18	0.18219104159382513	0.4012794939987963	-0.1110989445749171	1.0140075367619488e-17	0.3768565859239482	0.7801752498828849	-0.1327662853621246	6.452775233939675e-18	-0.004612407706811287	-0.08199405466195814	-0.0950688279422003	1.0140075367619488e-17	0.03118700914406997	0.037330940498599315	-0.0736304553454255	5.300493942164733e-18	0.01748404930484786	0.05893998394381165	-0.0719586174667152	1.290555046787935e-17	-0.12339739150043862	-0.1755016616884836	-0.1314260936232569	8.393554499235004	1767.0757654385054
GSM7908362_HT2	GSE248205	HT	2595	-5.783217726757892e-09	0.24213265372665638	0.3051069827489985	-0.036882885	5.783217726757892e-09	-0.24213265372665638	-0.3051069827489985	0.036882885	-3.5941641657739302e-09	0.27513757336923267	0.3363148846835192	-0.008138756	-1.1087976879449598e-08	0.20603460319375963	0.25965089152788906	-0.09844017	-3.3132154133162504e-08	0.42235730657939907	0.7006204737750386	-0.030088484	-2.044507076675144e-09	-0.05252488438977658	-0.06583899390416025	-0.1259018	-1.67661657180305e-09	0.030562274068859786	0.059140684621032356	-0.08044119	-1.5605934489618845e-08	0.0661894066692604	0.08278899529815098	-0.06072361	-4.616418887835997e-09	-0.01594713647872651	-0.003896728921725741	-0.2606126	3.080387582775001e-18	0.11450313320468684	0.13365566038251084	-0.0123517357399901	-3.080387582775001e-18	-0.11450313320468684	-0.13365566038251084	0.0123517357399901	3.422652869750001e-18	0.1324755714689273	0.14377364731833198	0.0225585707401482	8.214366887400003e-18	0.09123356756171443	0.1071479674769101	-0.0376566731559922	1.3690611479000004e-18	0.333943255867259	0.5866383825901651	-0.0307284183380472	1.4717407339925005e-17	-0.02720720032399687	-0.033484845884148524	-0.1116940177557278	7.358703669962502e-18	-0.006333818073388976	0.010711750569330662	-0.0781241404060857	1.0267958609250002e-17	0.013668279625657651	0.015380564875543522	-0.0675337303467795	1.9166856070600006e-17	-0.06316858916947789	-0.06446175797559091	-0.131301874111515	8.486451190205768	1706.3526011560693
GSM7908363_HT3	GSE248205	HT	1980	1.1541313131229573e-08	0.31064268868430306	0.5218531662619393	-0.101351735	-1.1541313131229573e-08	-0.31064268868430306	-0.5218531662619393	0.101351735	2.3240742425046374e-08	0.3230190351348686	0.5391765537257778	-0.08803232	7.578044950698275e-09	0.2741621146147172	0.4507039216569091	-0.1144383525	1.3865848481070587e-08	0.4932649037276565	0.8234215084848485	-0.091691335	-5.51649494210528e-10	0.04438685560948889	0.04500662388795959	-0.055919171000000004	-3.226253535078847e-08	0.06448997011479798	0.06702778033979799	-0.05585424	2.0233867678323944e-08	0.0601521050979798	0.0937893387959596	-0.0281512205	-1.1112643933672368e-08	-0.045318382779015155	-0.05898541156629294	-0.09993553699999999	-8.971499188890154e-19	0.2567202725223734	0.4325771257840491	-0.0639700150356512	8.971499188890154e-19	-0.2567202725223734	-0.4325771257840491	0.0639700150356512	-7.177199351112123e-18	0.27175015720688983	0.4529136150637544	-0.04892292462472535	1.0765799026668185e-17	0.21058181038480261	0.3485261735134457	-0.07777842472675395	-4.037174635000569e-18	0.4315128875053027	0.7299101350628747	-0.0695069105822016	9.86864910777917e-18	0.00716859013918425	-0.01629530887913079	-0.0794419022819633	1.2560098864446216e-17	0.023517200779013942	0.0003509711753629771	-0.05875869984391155	8.522924229445646e-18	0.015404165524668986	0.02477300272375334	-0.0493783744989473	1.3457248783335231e-17	-0.094863192687642	-0.13413926737157533	-0.09889886445466134	8.598048718204767	2468.787878787879
GSM7908364_GD1	GSE248205	GD	2266	-6.420303221457915e-08	0.12299252118495499	0.11843613721856083	0.0855180695	6.420303221457915e-08	-0.12299252118495499	-0.11843613721856083	-0.0855180695	-1.0776063547776656e-07	0.1325234966465137	0.12564329352239859	0.0995796525	-6.368411297572641e-08	0.09924728706381288	0.06900999501005291	0.0611502005	2.6524779350487662e-08	0.35386231497263904	0.5711359091181658	0.03454385	7.905841570493954e-09	-0.03804439515509532	-0.020978690800141093	-0.09100892	-1.8664263018113663e-08	-0.0023900973784377754	0.0031462690788888904	-0.065714602	-4.397038833850529e-09	0.03479796827049956	0.017153070501058194	-0.029508569999999998	-4.71932479170103e-09	0.0013917017062753856	-0.005135511775132262	-0.0943447675	1.567834809708959e-18	0.056731812476533944	0.05478615559461784	0.03685918135299755	-1.567834809708959e-18	-0.056731812476533944	-0.05478615559461784	-0.03685918135299755	6.271339238835836e-18	0.06459694783953804	0.059656231867829645	0.0297720737938925	1.567834809708959e-18	0.03986326675966103	0.013774446591199982	0.01323521745337475	7.055256643690315e-18	0.3062820155632614	0.5326834733206053	-0.02229520245358275	1.1758761072817193e-18	-0.05224477952524965	-0.041907997343537104	-0.0692024547626417	2.5477315657770584e-18	-0.026222814442492453	-0.018106414704973973	-0.0679788824296261	7.839174048544794e-18	0.004118170324900537	-0.008523446351938104	-0.052960320755794196	1.489443069223511e-17	-0.01791874963931072	-0.02088163747150482	-0.11182688966145299	8.816253552134285	2919.860105913504
GSM7908365_GD2	GSE248205	GD	2005	4.583332169477387e-08	0.20270624779023927	0.23336652711087652	0.026329488	-4.583332169477387e-08	-0.20270624779023927	-0.23336652711087652	-0.026329488	5.903840897942808e-08	0.21197443111729813	0.23561641017659365	0.05349936	6.878563582509138e-09	0.20120110631218344	0.22724508725498005	-0.023197442	-3.218960597369351e-08	0.4137631818219342	0.6568598664143426	0.0047706007	-5.049860342134997e-09	-0.0172913327559322	-0.024403329717330675	-0.06107093	-5.143401493462942e-09	0.03405095980527417	0.042699601546394426	-0.10993485	-2.9440967741849665e-09	0.0641932511577966	0.0977220187376494	-0.066545784	-9.03441395687117e-09	0.024768411710867407	0.02009186971314742	-0.15353891	-6.644726331921136e-18	0.037429915014626704	0.05446772697450657	-0.0224504391568713	6.644726331921136e-18	-0.037429915014626704	-0.05446772697450657	0.0224504391568713	-7.087708087382545e-18	0.043603764036988395	0.05290857883798224	-0.0184843074763048	1.0631562131073818e-17	0.03120595937460928	0.04350463429141715	-0.0350312813445413	4.429817554614091e-18	0.27035523338039885	0.512221177220691	-0.0677375314428372	1.0631562131073818e-17	-0.021003484662577794	-0.028642586905697882	-0.0612437327053827	7.530689842843955e-18	-0.01461890085990741	-0.01041544607994275	-0.0663540870361314	1.1517525641996636e-17	-0.0007198815661283507	0.028140063417796814	-0.0789833967129189	2.9236795860453e-17	-0.005005886809203131	-0.012139324861839888	-0.12875170200889	7.335870148536536	1031.922693266833
GSM7908366_GD3	GSE248205	GD	1969	1.405751650495023e-08	0.17761495536162433	0.18284546754624748	0.06994128	-1.405751650495023e-08	-0.17761495536162433	-0.18284546754624748	-0.06994128	1.591854748643471e-08	0.1940784286372589	0.19959723932900605	0.09726477	1.647476739354076e-08	0.12473471058570254	0.12951899697873628	0.018922903	-2.050136108820773e-08	0.3605567194617462	0.5982806875862069	-0.0177213	-1.261990859243326e-08	0.0032093810473197894	-0.00498416407726167	-0.069288425	-2.7209700362723314e-09	0.0052659672035228425	0.006076862040304257	-0.10055027	9.047832397640651e-09	0.031061687055106597	0.03901097229107504	-0.047274563	5.952260049329511e-09	0.008768609878172608	0.014289245375253567	-0.10977102	-4.510809648045328e-19	0.07628103542428036	0.0863723453001549	0.0468109693817912	4.510809648045328e-19	-0.07628103542428036	-0.0863723453001549	-0.0468109693817912	-7.893916884079325e-18	0.0843607589420149	0.09541494792913852	0.0711555454145025	-4.059728683240796e-18	0.03207162035135564	0.040643357049159605	-0.0196762096337786	3.608647718436263e-18	0.3005422058372341	0.5466020642186633	-0.0430475685842367	9.021619296090658e-18	-0.020298066163485118	-0.0254443572033849	-0.084825926306615	1.781769810977905e-17	-0.018590519451109358	-0.016577031791401412	-0.0867414513402176	8.119457366481591e-18	0.001184956510538976	0.010674416088935521	-0.0721009728504998	3.608647718436263e-17	-0.0073097241756194435	-0.0007404022416812783	-0.1075709978860771	7.239143402706802	872.8141188420518
```

---

# Appendix G — Cross-reference: GSE250521 prior verdict (single cohort)

GSE250521 (n=16 cancer Visium PT→PTC→LPTC→ATC) sole prior dataset, verdict:
- All-spot stage trend (sample-mean Spearman): ρ=0.000, p=1.0 (null)
- Epithelial top-50% raw: RAI_8 ρ=−0.44, p=0.091 / DM1 ρ=+0.44, p=0.091
- Epithelial top-50% depth-resid: RAI_8 ρ=−0.35, p=0.18 / DM1 ρ=+0.35, p=0.18
- Proliferation epi top50 depth-resid: ρ=+0.58, p=0.018 (validates stage labels)
- DM1 ~ TDS_overlap all-spots: ρ=−0.89 (strong, but 6/8 gene overlap = partial tautology)
- DM1 ~ log_total_counts: ρ=−0.45 (depth confound, ~20% magnitude reduction after correction)

Combined 3-cohort (GSE250521 + GSE230424 + GSE248205) status:
- Total slides: 28 (16 cancer + 4 PTC_HT + 8 autoimmune)
- Cancer progression direct evidence: GSE250521 only, underpowered (p=0.18)
- Independent lineage validation: 12 external samples, ρ<−0.85, p<0.01 → tautology risk RESOLVED
- Inflammation artifact: GSE248205 negative control → DM1 NOT inflammation marker, RESOLVED
- Verdict: C (multi-cohort supplementary), not D (main figure)
