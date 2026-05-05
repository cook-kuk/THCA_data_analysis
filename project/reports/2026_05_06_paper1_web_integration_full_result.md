# Paper 1 web integration full result

**Date:** 2026-05-06  
**Live URL:** http://40.82.129.113/papers_hub_2026_05_04/paper1.html  
**Index URL:** http://40.82.129.113/  
**Scope:** Paper 1 전체 흐름 웹 통합. 논문 본문 prose가 아니라 비전공자도 이해 가능한 해설형 report / figure guide / table guide / claim boundary.

---

## 1. Integrated Decision

Paper 1의 최종 전략 정체성은 다음으로 고정했다.

> **A thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations**

핵심 결정:

- GPL570 external expression validation은 **lineage axis에 대해 HIT**.
- GSE76039 + GPL570 pack은 external expression validation 근거로 통합.
- 8-gene은 “새 discovery panel”이 아니라 **compact RAI-lineage readout**.
- TROP2/TACSTD2는 GPL570 bulk external replication이 heterogeneous/negative이므로 **title/main claim에서 demoted**.
- Paper 1의 본체는 **driver-orthogonal thyroid-lineage differentiation axis**.
- 외부 dataset 추가, 새 분석, raw alignment, RunPod/GPU, WSI/H&E 재시도 없음.

---

## 2. Files Changed / Created

### Web files

```text
project/papers_hub_2026_05_04/paper1.html
project/three_papers_index.html
project/papers_hub_2026_05_04/assets/paper1/
```

### Strategy memo

```text
project/reports/2026_05_04_paper1_final_strategy_after_external_validation.md
```

### Public server sync

The live server root is:

```text
/var/www/papers/
```

Synced files:

```text
/var/www/papers/index.html
/var/www/papers/three_papers_index.html
/var/www/papers/papers_hub_2026_05_04/paper1.html
/var/www/papers/papers_hub_2026_05_04/lightbox.js
/var/www/papers/papers_hub_2026_05_04/assets/paper1/
```

---

## 3. Current Live Page Structure

The live Paper 1 page is organized for non-specialist reading:

1. **한 줄 결론**
2. **비전공자 배경**
3. **전체 논리 흐름**
4. **핵심 결과와 그림 해설**
5. **외부 검증: GSE76039 + GPL570 pack**
6. **표로 보는 claim boundary**
7. **Figure / table plan**
8. **TROP2 demotion**
9. **다음 행동**

Design intent:

- 먼저 “무슨 논문인지”를 쉽게 설명.
- 그 다음 “왜 driver mutation만으로 부족한지” 설명.
- 8-gene이 무엇을 하는지 과장 없이 설명.
- 외부 validation 그림들을 Results flow 안에 배치.
- TROP2 overclaim을 명시적으로 차단.

---

## 4. Plain-Language Summary

갑상선은 원래 요오드를 모아 갑상선호르몬을 만드는 기관이다. 정상 갑상선 세포는 요오드 흡수와 갑상선 기능에 필요한 유전자를 켜고 있다. 암이 더 공격적인 상태로 변하면 이 “갑상선다운 기능”을 잃는다. 이 기능 상실을 expression data에서 하나의 축으로 읽은 것이 Paper 1의 핵심이다.

쉽게 말하면:

- **BRAF/RAS mutation**은 자동차의 엔진 종류에 가깝다.
- **Thyroid-lineage differentiation axis**는 지금 차가 정상 주행 중인지, 고장 모드로 가는지 보여주는 계기판에 가깝다.

Paper 1은 “새로운 8-gene panel을 발견했다”가 아니라, 갑상선암에서 **driver mutation과 별개로 작동하는 lineage-silencing expression state**가 있고, 이 축이 TCGA와 외부 GEO expression datasets에서 반복적으로 보인다는 주장으로 정리했다.

---

## 5. Main Logic Flow

1. **Driver mutation만으로는 부족하다.**  
   BRAF/RAS/TERT transcript signal과 driver-only clustering은 DM1/DM2 axis를 잘 설명하지 못한다.

2. **더 넓은 transcriptome에서도 같은 축이 나온다.**  
   TIERA67과 pan-genome top-MAD gene clustering이 높은 ARI로 같은 구조를 재현한다.

3. **8-gene은 discovery claim이 아니라 compact readout이다.**  
   RAI/thyroid-lineage biology를 간단히 측정하는 계기판 역할이다.

4. **Zero-overlap validation이 핵심 방어다.**  
   8-gene과 겹치지 않는 THYROID_NONOVERLAP gene set도 같은 방향으로 움직인다.

5. **External expression cohort에서도 재현된다.**  
   GSE76039와 GPL570 pack에서 advanced disease가 lineage-silenced direction으로 간다.

6. **TROP2는 내려놓는다.**  
   TCGA bulk에서는 신호가 있지만 GPL570 bulk external replication은 heterogeneous/negative라 제목과 main claim에서 제외한다.

---

## 6. Figure Assets Integrated

Assets copied into:

```text
project/papers_hub_2026_05_04/assets/paper1/
```

Current asset count: **20 PNG files**.

Key figures:

```text
F3_driver_neutrality.png
F4_pangenome_robustness.png
s_tcga_4panel.png
SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.png
m_master_regulator_volcano.png
q1_gsea_hallmark.png
q5_thyroid_tf_network.png
external_rai_lineage_boxplots.png
external_dm1_nonoverlap_scatter_grid.png
external_direction_consistency_forest.png
external_dataset_qc_heatmap.png
gse76039_dm1_lineage_boxplot.png
gse76039_dm1_vs_nonoverlap_scatter.png
gse76039_mechanism_heatmap.png
n1_dark_matter.png
q3_trop2_spatial.png
o_multivariate_cox.png
s_tcga_km_dm1.png
fig_2x4_RAI_DM1_per_stage.png
fig_sample_means_RAI_DM1_TDS.png
```

---

## 7. Figure Explanations

### 7.1 Driver orthogonality

**Figure:** `F3_driver_neutrality.png`

Purpose:

- Shows that BRAF/RAS/TERT driver information alone does not define the DM1/DM2 expression axis.
- Supports the claim that Paper 1 is not another BRAF/RAS dichotomy paper.
- The 4.1 section was widened because the original figure’s yellow annotation looked visually crowded at normal page width.

UI update:

- The 4.1 figure now uses a wider display class: `driver-wide`.
- All figures can be clicked to zoom via `lightbox.js`.

### 7.2 Pan-genome robustness

**Figure:** `F4_pangenome_robustness.png`

Purpose:

- Shows the axis is not an artifact of using only 8 genes.
- TIERA67 and pan-genome top-MAD gene sets recover the same axis.
- Driver_anchor alone performs near random.

Plain-language meaning:

- The signal is broad transcriptome biology.
- The 8-gene readout is useful because it is compact and interpretable, not because it alone created the axis.

### 7.3 TCGA zero-overlap anchor

**Figure:** `s_tcga_4panel.png`

Purpose:

- Shows TCGA bulk DM1_like vs THYROID_NONOVERLAP anti-correlation.
- THYROID_NONOVERLAP shares zero genes with RAI_8.

Plain-language meaning:

- If two independent gene sets tell the same story, the signal is less likely to be a panel artifact.

### 7.4 Spatial non-overlap validation

**Figure:** `SuppFig_X2_DM1_thyroid_nonoverlap_crossvalidation.png`

Purpose:

- Shows the same axis in spatial transcriptomics slides.
- Supports the idea that the lineage axis is not only a bulk-expression artifact.

Boundary:

- This is supportive spatial validation.
- It is not H&E inference and not WSI prediction.

### 7.5 Mechanism support

Figures:

```text
m_master_regulator_volcano.png
q1_gsea_hallmark.png
q5_thyroid_tf_network.png
```

Purpose:

- Supports lineage transcription-factor collapse.
- Shows STAT3/AP1/DNMT directionality.
- Shows pathway-level consistency.

Boundary:

- This is supportive mechanism evidence.
- It does not prove causal TF-to-DNMT-to-methylation biology.
- Use “supports” or “consistent with,” not “demonstrates.”

### 7.6 External GPL570 validation

Figures:

```text
external_rai_lineage_boxplots.png
external_dm1_nonoverlap_scatter_grid.png
external_direction_consistency_forest.png
external_dataset_qc_heatmap.png
```

Purpose:

- Shows lineage scores by histology group.
- Shows DM1_like vs THYROID_NONOVERLAP zero-overlap anti-correlation.
- Shows effect-size direction consistency.
- Shows panel gene coverage.

Plain-language meaning:

- The axis reappears in independent thyroid expression datasets.
- The result is external expression validation, not clinical validation.

### 7.7 GSE76039 anchor

Figures:

```text
gse76039_dm1_lineage_boxplot.png
gse76039_dm1_vs_nonoverlap_scatter.png
gse76039_mechanism_heatmap.png
```

Purpose:

- Anchors the advanced-disease PDTC/ATC signal.
- Supports advanced-disease replication.

Boundary:

- No survival claim.
- No mutation claim.
- No PTC-to-PDTC-to-ATC progression proof.

### 7.8 Honest negatives

Figures:

```text
n1_dark_matter.png
q3_trop2_spatial.png
```

Purpose:

- Shows dark-matter subset outcome claim is not headline-grade.
- Shows TROP2 spot-level colocalization is negative.

Boundary:

- Do not claim “DM1-high spots are TROP2-high.”
- Do not claim TROP2 external validation.

---

## 8. External Validation Summary

### Cohorts

| Dataset | Role | Samples | Status |
|---|---|---:|---|
| GSE76039 | Advanced-disease anchor | 37 | Used |
| GSE33630 | Main external replication | 105 | Used |
| GSE65144 | Main supporting | 25 | Used |
| GSE29265 | Supplement | 49 | Used |
| GSE53157 | Sensitivity | 26 used | Used, underpowered PDTC |
| GSE126698 | Tier 1 candidate | 28 metadata | Dropped |
| GSE53072 / GSE60542 / GSE120177 | Hold | NA | Not downloaded |

### Main result

The external expression validation supports:

- RAI_8 silencing in advanced disease.
- DM1_like increase when RAI_8 decreases.
- THYROID_NONOVERLAP decrease with DM1_like increase.
- TDS_like and TF_collapse direction consistency.
- STAT3/AP1/DNMT supportive directionality.

It does **not** support:

- TROP2/TACSTD2 as a robust bulk-microarray external replication claim.
- Survival validation.
- Clinical validation.
- Fusion validation.
- Progression proof.

---

## 9. Claim Boundary Table

| Claim | Status | Plain-language meaning | Boundary |
|---|---|---|---|
| Driver-orthogonal thyroid-lineage differentiation axis | MAIN | Driver mutation과 별개로 갑상선 기능 상실 expression axis가 있다. | Title/main claim OK |
| 8-gene compact RAI-lineage readout | MAIN TOOL | 축을 간단히 측정하는 계기판 | Not discovery premise |
| Zero-overlap module validation | MAIN DEFENSE | 겹치지 않는 gene set도 같은 축을 읽는다 | Strong anti-artifact argument |
| GPL570 external expression validation | MAIN/SUPP | 외부 thyroid expression data에서도 같은 방향 | Expression validation only |
| STAT3/AP1/DNMT mechanism arm | SUPPORTIVE | 가능한 mechanism arm이 같은 방향으로 움직인다 | Not causal proof |
| TACSTD2/TROP2 | DEMOTED | TCGA signal은 있으나 GPL570 bulk replication 불안정 | No title/main claim |
| H&E/WSI inference | DROP | 이미지에서 DM1 예측은 현재 제외 | No retry |

---

## 10. Wording Rules

### Safe wording

- external expression validation
- direction-consistent lineage silencing
- advanced-disease replication
- zero-overlap module validation
- driver-orthogonal thyroid-lineage differentiation axis
- compact RAI-lineage readout
- supportive mechanism arm

### Forbidden wording

- TROP2-targetable as headline
- TROP2 vulnerability in title/main claim
- progression proven
- clinical validation
- survival validation
- fusion validation
- DM1-high spots are TROP2-high
- H&E-inferable
- all risks resolved
- Cancer Cell-ready / Nat Cancer reach / venue probability claims

---

## 11. UI Updates Completed

### 11.1 Live page update

Page:

```text
http://40.82.129.113/papers_hub_2026_05_04/paper1.html
```

The page was rewritten as a compact but complete explanatory report:

- clearer title
- non-specialist background
- result flow
- figure-by-figure explanations
- table-by-table claim boundary
- TROP2 demotion section
- next action section

### 11.2 Main index update

Page:

```text
http://40.82.129.113/
```

The Paper 1 card now uses:

```text
A thyroid-lineage differentiation axis stratifies thyroid cancer beyond canonical driver mutations
```

and links to:

```text
papers_hub_2026_05_04/paper1.html
```

### 11.3 Image click-to-zoom

All figure images now use:

```text
lightbox.js
```

Behavior:

- Click any figure image to zoom.
- Click background or press `Esc` to close.
- Captions are shown in the lightbox.

### 11.4 4.1 Driver orthogonality visual fix

The 4.1 driver orthogonality figure was widened using:

```text
figure.driver-wide
```

Reason:

- The original inline width made the figure annotation visually crowded.
- The source PNG is high-resolution (`3000 x 2200`), so browser-side larger display is appropriate.

---

## 12. Verification

HTTP checks passed:

```text
http://40.82.129.113/                                                200
http://40.82.129.113/papers_hub_2026_05_04/paper1.html               200
http://40.82.129.113/papers_hub_2026_05_04/lightbox.js               200
http://40.82.129.113/papers_hub_2026_05_04/assets/paper1/F3_driver_neutrality.png 200
http://40.82.129.113/papers_hub_2026_05_04/assets/paper1/external_rai_lineage_boxplots.png 200
http://40.82.129.113/papers_hub_2026_05_04/assets/paper1/external_dm1_nonoverlap_scatter_grid.png 200
```

Overclaim sweep on active web files:

```text
No active hits for:
TROP2-targetable
TROP2 vulnerability in title
tumor-population TROP2 vulnerability
TROP2 elevated in advanced disease confirmed
progression proven
clinical validation
survival validation
fusion validation
Cancer Cell-ready
Nat Cancer reach
all risks resolved
H&E-inferable
```

Gitignore check:

```text
project/results/p_external_expression_validation/raw/
```

is present in `.gitignore`.

---

## 13. Commit Proposal

Do not commit automatically.

### W1 web integration

```text
project/papers_hub_2026_05_04/paper1.html
project/three_papers_index.html
project/papers_hub_2026_05_04/assets/paper1/
```

Suggested message:

```text
web: integrate Paper 1 lineage-axis external validation explainer
```

### W2 strategy memo

```text
project/reports/2026_05_04_paper1_final_strategy_after_external_validation.md
project/reports/2026_05_06_paper1_web_integration_full_result.md
```

Suggested message:

```text
docs: record Paper 1 final strategy after external validation
```

### W3 external validation canonical result

```text
project/reports/2026_05_04_external_expression_FULL_RESULTS.md
project/results/p_external_expression_validation/scripts/run_external_validation.py
project/results/p_external_expression_validation/scripts/make_figures.py
```

Suggested message:

```text
docs+scripts: add external expression validation canonical report
```

### W4 processed outputs and figures

```text
project/results/p_external_expression_validation/sample_metadata.tsv
project/results/p_external_expression_validation/external_gene_coverage.tsv
project/results/p_external_expression_validation/external_score_tests.tsv
project/results/p_external_expression_validation/external_direction_consistency.tsv
project/results/p_external_expression_validation/external_spearman.tsv
project/results/p_external_expression_validation/external_sample_scores.tsv.gz
project/results/p_external_expression_validation/*_expression_gene_log.tsv.gz
project/results/p_external_expression_validation/*.png
```

Suggested message:

```text
data: add Paper 1 external expression validation outputs
```

### W5 gitignore

```text
.gitignore
```

Suggested message:

```text
chore: ignore external expression raw GEO downloads
```

Exclude:

```text
project/results/p_external_expression_validation/raw/
Paper 2/3/4 files
voice-protected prose
WSI / tiles / embeddings / RunPod runtime
```

---

## 14. Final Status

External validation integrated. TROP2 demoted. No more datasets. Return to voice-hook.
