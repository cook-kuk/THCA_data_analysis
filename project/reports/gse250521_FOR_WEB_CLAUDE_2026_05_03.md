# GSE250521 Paper 1 inclusion 판정 — Web Claude 2nd opinion 요청
**작성: 2026-05-03 (Seungho Cook + Claude Code Opus 4.7)**

## 0. 요청

Web Claude (claude.ai) 야, 아래 7-criteria cold review 의 verdict (Supplementary only, main figure 불가) 에 동의하는지 독립적으로 판단해줘. 특히:
1. **DM1_like ~ TDS_like ρ=−0.89 가 6/8 gene overlap 으로 partial tautological** 이라는 우려가 정당한가
2. **EMT 음의 상관 (ρ=−0.09 spot, −0.06 sample)** 을 reverse-causality risk 로 간주해 Paper 1 에서 제외하는 게 맞는가
3. **단일 cohort n=16, depth-corrected p=0.18** 인 이 결과를 supplementary 로 포함하는 vs drop 의 trade-off 가 어느 쪽이 옳은가
4. Paper 1 main figure 자격 upgrade 하려면 (HRA003537 / GSE230424 meta) marathon mode (5/4-6/13 manuscript writing, 6/13 bioRxiv target) 을 위반해서까지 할 가치가 있는가

답변 형식:
- 동의/부분동의/반대 + 이유
- 추가 분석 권장 (있다면)
- Marathon discipline 깨고 확장할지 vs 현 상태 supplementary 진행할지

---

## 1. Paper 1 context (Web Claude 가 알아야 할 framework)

**Paper 1**: BRAF/RAS-negative PTC dark matter 에서 8-gene driver-excluded sub-stratifier (DM1/DM2 subtypes). bioRxiv target **2026-06-13**, marathon writing 5/4-6/13. Base venue Sci Rep, reach Cell Rep Med / JCI Insight.

**8 RAI genes**: TPO, DIO1, TSHR, PAX8, TG, FOXE1, NKX2-1, SLC5A5

**Score definitions**:
- `RAI_8_score` = within-sample mean z-score of 8 RAI genes (radioactive iodine uptake / thyroid lineage preservation)
- `DM1_like_score` = `−RAI_8_score` (mathematically mirrored; high DM1 = dedifferentiated)
- `TDS_like_score` = within-sample mean z-score of 6 thyroid-differentiation genes: TG, TPO, TSHR, SLC5A5, DIO1, PAX8 — **6/8 overlap with RAI_8**
- CAF_ECM (10 genes), EMT (10), Hypoxia (9), Proliferation (7), Epithelial (5)

**Hypothesis (a priori)**:
1. RAI_8 decreases PT→PTC→LPTC→ATC
2. DM1_like increases PT→PTC→LPTC→ATC
3. DM1_like and TDS_like negative correlation
4. DM1_like correlates with CAF/EMT/hypoxia/ECM (mechanism for dedifferentiation)

**Cross-paper isolation rules** (web Claude judgment 시 준수):
- Paper 1 = cancer dedifferentiation only (이 분석)
- Paper 2 = Hashimoto-overlap PTC mechanism (TLS / BCR / AICDA) — touch X
- Paper 3 = Korean Graves' HLA — touch X

---

## 2. Data — GSE250521 (Lu et al. thyroid Visium ST)

| stage | label | n_slides | n_spots_post_QC |
|---|---|---|---|
| PT (para-tumor) | N-1..4 | 4 | 14,740 |
| PTC | PTC-1..4 | 4 | 16,223 |
| LPTC | LPTC-1..4 | 4 | 13,634 |
| ATC | ATC-1..4 | 4 | 11,276 |
| **total** | | **16** | **55,873** |

QC: in_tissue + n_genes ≥ 200, within-sample log1p z-score normalize. paired scRNA-seq (9 samples, GSM7980876-7980884) **excluded** from this analysis.

**Gene availability** (per sample, all 16 slides identical): 8/8 RAI · 6/6 TDS · 10/10 CAF_ECM · 10/10 EMT · 9/9 Hypoxia · 7/7 Proliferation · 5/5 Epithelial. **No dropout**.

---

## 3. Stage trend results (PT=0, PTC=1, LPTC=2, ATC=3)

### 3.1 All spots, sample-mean Spearman (n=16 slides)

| score | ρ | p |
|---|---|---|
| RAI_8_score | 0.000 | 1.000 |
| DM1_like_score | 0.000 | 1.000 |
| TDS_like_score | −0.109 | 0.687 |
| CAF_ECM_score | −0.218 | 0.417 |
| EMT_score | +0.291 | 0.274 |
| Hypoxia_score | +0.085 | 0.755 |
| Proliferation_score | −0.243 | 0.365 |

**Read**: All-spot stage trend at sample-mean level is **null across the board**.

### 3.2 Epithelial top-50% subset, sample-mean Spearman (n=16)

Epithelial mask = top 50% of `Epithelial_score` (mean z of EPCAM/KRT8/KRT18/KRT19/TACSTD2) **within each sample**. NOT a tumor mask — just epithelial-enriched.

| score | ρ (raw) | p (raw) | ρ (depth-corrected) | p (depth-corrected) |
|---|---|---|---|---|
| RAI_8_score | **−0.437** | **0.091** | −0.352 | 0.182 |
| DM1_like_score | **+0.437** | **0.091** | +0.352 | 0.182 |
| TDS_like_score | **−0.461** | **0.072** | −0.327 | 0.216 |
| CAF_ECM_score | −0.206 | 0.444 | **−0.521** | **0.038** |
| EMT_score | +0.279 | 0.296 | −0.352 | 0.182 |
| Hypoxia_score | +0.243 | 0.365 | −0.049 | 0.858 |
| Proliferation_score | **+0.800** | **2.0e-4** | **+0.582** | **0.018** |

Depth correction = within-sample residualization of score against `log(total_counts) + log(n_genes_by_counts)` before sample-mean.

**Read**:
- Raw epithelial RAI_8/DM1/TDS direction-consistent (RAI ↓, DM1 ↑, TDS ↓) at p ≈ 0.07-0.09
- After depth correction, magnitude reduces ~20%, p inflates to 0.18 (NS)
- **Proliferation gradient is the only RAI-relevant biology that survives depth correction** (p=0.018) — confirms stage labels are biologically real
- CAF_ECM goes **down with stage** in epithelial subset after depth correction (p=0.038) — likely over-correction artifact, not interpreted

### 3.3 Mixed-effects model (mixedlm_stage_ord)

All p > 0.99 across all scores and subsets. **Misspecified** — stage is sample-level, so `(1|sample_id)` random intercept absorbs all stage variance. Reported in CSV for transparency, **not used for inference**.

---

## 4. Correlation matrix (DM1_like as reference)

### 4.1 All spots (n=55,873)

| pair | Spearman ρ | p | Pearson r |
|---|---|---|---|
| DM1_like ~ TDS_like | **−0.888** | ≈0 | −0.919 |
| DM1_like ~ log(total_counts) | **−0.448** | ≈0 | −0.399 |
| DM1_like ~ CAF_ECM | +0.245 | ≈0 | +0.293 |
| DM1_like ~ EMT | **−0.090** | 5e-100 | −0.032 |
| DM1_like ~ Hypoxia | −0.258 | ≈0 | −0.187 |
| DM1_like ~ Proliferation | −0.182 | ≈0 | −0.091 |

### 4.2 Sample-mean (n=16)

| pair | Spearman ρ | p | Pearson r |
|---|---|---|---|
| DM1_like ~ TDS_like | **−0.874** | 1.0e-5 | −0.987 |
| DM1_like ~ CAF_ECM | +0.100 | 0.71 | +0.205 |
| DM1_like ~ EMT | −0.056 | 0.84 | −0.017 |
| DM1_like ~ Hypoxia | **+0.706** | 0.002 | +0.534 |
| DM1_like ~ Proliferation | +0.003 | 0.99 | +0.168 |

**Simpson-flip on Hypoxia**: spot-level ρ=−0.26 (negative) vs sample-mean ρ=+0.71 (positive). Aggregate level alignment with hypoxia is real, but spot-level overwhelmed by within-sample technical variance.

---

## 5. 7-Criteria cold verdict

| # | 기준 | Verdict |
|---|---|---|
| 1 | DM1 ↑ PT→ATC? | ⚠️ MARGINAL (epi only, depth-corrected p=0.18) |
| 2 | RAI_8 ↓ PT→ATC? | ⚠️ MARGINAL (1번과 mathematically mirror) |
| 3 | DM1 ↔ TDS 음의 상관? | ✅ STRONG (ρ=−0.89) **but 6/8 gene overlap → partial tautology** |
| 4 | counts/ngenes artifact? | ⚠️ PARTIAL (~20% magnitude reduction after depth correction) |
| 5 | Epithelial subset 유지? | ⚠️ DIRECTIONALLY YES, INFERENTIALLY NO |
| 6 | CAF/EMT/hypoxia 연결? | ❌ WEAK / **EMT spot-level ρ=−0.09 wrong direction** |
| 7 | 등급? | **SUPPLEMENTARY 1 figure + Methods 1줄 + Discussion 1줄** |

### 5.1 핵심 risk 두 개

**(a) TDS_like 와 6/8 gene overlap**
- TDS_like = mean_z(TG, TPO, TSHR, SLC5A5, DIO1, PAX8)
- RAI_8 = mean_z(TPO, DIO1, TSHR, PAX8, TG, FOXE1, NKX2-1, SLC5A5)
- 공통 6 genes: TG, TPO, TSHR, SLC5A5, DIO1, PAX8
- DM1_like = −RAI_8 → DM1_like ~ TDS_like 의 ρ=−0.89 는 부분적으로 동어반복
- **Reviewer 가 지적할 가능성 높음**

**(b) EMT spot-level 음의 상관**
- DM1 high = dedifferentiated tumor 가설 → EMT 양의 상관 기대
- 실제: spot-level ρ=−0.09 (방향 반대), sample-mean ρ=−0.06 (NS)
- Paper 1 mechanism claim 직접 충돌 → reviewer 가 reverse-causality 의심 가능

### 5.2 Verdict reasoning

- ✅ Direction-consistent (RAI ↓ / DM1 ↑ / TDS ↓ in epithelial spots): Paper 1 framework 와 호환
- ✅ Internal consistency (DM1 ~ TDS ρ=−0.89): score 가 thyroid lineage 측정한다는 sanity check
- ✅ Stage labels valid (proliferation gradient survives depth correction p=0.018): cohort 자체 신뢰
- ✅ 8-gene 전부 detect, no dropout: methods note 가치
- ❌ Sample-mean stage trend underpowered (epi corrected p=0.18)
- ❌ Mechanism (CAF/EMT/hypoxia) 직접 지지 안 됨
- ❌ TDS overlap 으로 partial redundancy

→ **Drop 안 함, main figure 도 안 됨, supplementary 가 정확한 자리.**

### 5.3 Paper 1 활용 권장 (확정안)

| 자리 | 내용 | 출처 파일 |
|---|---|---|
| Supp Fig N | 2-row × 4-stage spatial heatmap (RAI_8 / DM1_like) | `project/results/figures_for_advisor/fig_2x4_RAI_DM1_per_stage.png` |
| Supp Fig N+1 | 16-slide sample-mean strip plot (RAI_8 / DM1_like / TDS_like) | `project/results/figures_for_advisor/fig_sample_means_RAI_DM1_TDS.png` |
| Supp Table | stage_trend_summary + correlation_summary + depth_corrected | 3 CSV |
| Methods 1줄 | "8 RAI genes scored as within-sample z-score; sample-mean Spearman vs ordinal stage; depth correction by within-sample residualization on log_counts + log_ngenes" | new |
| Discussion 1줄 (본인 voice) | "Single-cohort GSE250521 (n=16) Visium spatial profiling showed direction-consistent but underpowered DM1-axis trend (epithelial-enriched ρ=−0.35 after depth correction, p=0.18); independent multi-cohort spatial meta-analysis required for inferential confirmation." | new |

**CAF/EMT/hypoxia 분석은 Paper 1 에서 완전히 제외** (mechanism claim 위험 / EMT wrong direction).

### 5.4 등급 upgrade 시도 시 필요한 분석 (marathon discipline 위반)

1. HRA003537 (Chinese Visium thyroid ~30 slides) ingest
2. GSE230424 (additional thyroid Visium) ingest
3. Combined n≈60 으로 meta sample-mean Spearman (random-effects)
4. Independent thyroid-lineage axis (GATA3 / HHEX / FOXA1) 로 DM1 cross-validate (TDS overlap 우회)
5. Pathologist H&E annotation 으로 tumor mask (현재 epithelial top-50% 대체)

비용: 1 주 ingest + 1 주 통합 + 1 주 reanalysis = 3 주 → marathon (5/4-6/13, 6 weeks) 의 절반. bioRxiv 6/13 target 슬립 risk.

---

## 6. 사용자 marathon 상황

- 5/4 marathon Day 1 시작 (오늘 5/3 = Day 0)
- 5/4-6/13 = 6 weeks manuscript writing only
- Voice-protected sections (Hook / Aim / Disc 3.1 / Limitations / Cover Letter Para 1 / Q9): 본인 키보드
- Default = scaffolding / infra / analysis 만 Claude 가능
- 새 분석은 **paper-blocking only** (memory v17_marathon_mode_post_pillar1)

본 분석 자체는 marathon 시작 전에 끝남 (총 wall-time 5분). 등급 upgrade 위해 marathon 중간에 multi-cohort 확장하면 marathon discipline 위반.

---

## 7. 의사결정 옵션 (web Claude 가 추천 줘)

**옵션 A — 현 상태 supplementary 진행 (본인 + Claude Code 추천)**
- 1 supp fig + 1 supp table + Methods 1줄 + Discussion 1줄
- CAF/EMT/hypoxia 제외
- Marathon 그대로 진행, bioRxiv 6/13
- Risk: Paper 1 reach venue 정당화에 spatial 기여 약함

**옵션 B — 확장 후 main figure 시도 (marathon 위반)**
- HRA003537 + GSE230424 ingest, n≈60 meta
- 결과 따라 main figure 또는 supplementary
- Marathon 3 주 소진, bioRxiv 6/13 → 7/4 슬립
- Risk: 합쳐도 p<0.05 안 나올 수 있음

**옵션 C — Drop**
- 본 분석 commit 후 Paper 1 에서 spatial 빼고 marathon
- Risk: direction-consistent 결과 + 5분 wall-time 분석 버리는 것 (의미 있는 internal consistency 결과 손실)

---

## 8. 부록 — 산출 파일 inventory

```
project/data/raw/GSE250521/                       — 1.27 GB tar + 96 extracted gz (gitignored)
project/data/processed/GSE250521/                 — 16 raw + scored .h5ad, tiles 320 MB (gitignored)
project/data/processed/GSE250521/sample_metadata.tsv    — 16 sample stage map
project/results/00_qc/scoring_summary.tsv               — gene set availability per sample
project/results/00_qc/{sid}.spot_scores.tsv.gz   — 16 per-sample spot tables
project/results/01_spatial_score/all_spots_scored.tsv.gz       — 55,873 × 22 combined
project/results/01_spatial_score/per_sample/*.png             — 16 per-slide RAI/DM1/TDS heatmaps
project/results/02_stage_trend/stage_trend_summary.csv         — primary numeric output
project/results/02_stage_trend/correlation_summary.csv         — DM1 vs everything
project/results/02_stage_trend/depth_corrected_trend.csv       — depth-residualized re-test
project/results/02_stage_trend/violin_*.png                   — all-spot + epithelial violin
project/results/03_pathology_poc/tile_metadata.tsv.gz         — 3,200 tiles, 800/stage
project/results/figures_for_advisor/fig_2x4_RAI_DM1_per_stage.{png,pdf}
project/results/figures_for_advisor/fig_sample_means_RAI_DM1_TDS.{png,pdf}
project/reports/gse250521_dm1_spatial_validation_report.md          — 1-page advisor brief
project/reports/gse250521_meaningful_or_not_2026_05_03.md           — meaningful/negative split
project/reports/gse250521_paper1_inclusion_verdict_2026_05_03.md    — 7-criteria cold review
project/reports/gse250521_FOR_WEB_CLAUDE_2026_05_03.md              — this file
project/src/00_download/download_gse250521.sh
project/src/01_parse_visium/parse_gse250521.py
project/src/02_qc_score/score_8gene_dm1.py
project/src/03_spatial_plots/plot_spatial_scores.py
project/src/04_stage_trend/stage_trend_analysis.py
project/src/04_stage_trend/depth_corrected_trend.py
project/src/05_pathology_poc/extract_tiles.py
project/src/05_pathology_poc/README.md   — foundation-model embedding skeleton
```

전체 wall-time end-to-end: **~5 분** (no GPU). Reproducible with the 7 commands in `project/reports/gse250521_dm1_spatial_validation_report.md` §11.
