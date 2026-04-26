# SNUBH prospective validation shortlist

_Build: 2026-04-25 · 20 candidates from v12 top-20 novel deep dive,
ranked by composite SNUBH validation priority score._

## Scoring

`snubh_priority_score = 0.30·min(|d|/2.5,1) + 0.20·rep_rate + 0.20·novelty + 0.20·druggable + 0.10·min(BRAF_pp/5,1)`

- `|d|` = absolute Cohen's d in TCGA THCA (BRAF vs RAS).
- `rep_rate` = fraction of external cohorts (GSE27155, GSE126698) replicating direction.
- `novelty` = 1.0 (PubMed-novel <5 thyroid papers), 0.5 (emerging 5-24), 0 (established ≥25).
- `druggable` = 1.0 if in the v7 8-target druggable shortlist, else 0.
- `BRAF_pp` = PubMed papers co-mentioning the gene with BRAF (any cancer).

## Top tier (score ≥ 0.65) — first wave for SNUBH validation

| Rank | Gene | Score | log2FC | Cohen's d | FDR | Replicated | Novelty | Drug? | Rationale |
|---:|---|---:|---:|---:|---|---:|---|---|---|
| 1 | **PLEKHA6** | 0.906 | +2.35 | +2.22 | 3.0e-54 | 2 | novel | yes | v7-druggable · PubMed-novel (<5 thyroid) · 2 BRAF-context paper(s) · large effect (|d|=2.22) · fully replicated (2/2 cohorts) |
| 2 | **PTPRE** | 0.891 | +2.38 | +2.26 | 3.8e-45 | 2 | novel | yes | v7-druggable · PubMed-novel (<5 thyroid) · 1 BRAF-context paper(s) · large effect (|d|=2.26) · fully replicated (2/2 cohorts) |
| 3 | **CYP1B1** | 0.846 | +3.47 | +2.22 | 9.0e-53 | 2 | emerging | yes | v7-druggable · emerging (5-24 thyroid) · 4 BRAF-context paper(s) · large effect (|d|=2.22) · fully replicated (2/2 cohorts) |
| 4 | **B3GNT3** | 0.842 | +3.78 | +2.02 | 8.0e-52 | 2 | novel | yes | v7-druggable · PubMed-novel (<5 thyroid) · large effect (|d|=2.02) · fully replicated (2/2 cohorts) |
| 5 | **GABRB2** | 0.809 | +3.70 | +2.24 | 2.9e-47 | 2 | emerging | yes | v7-druggable · emerging (5-24 thyroid) · 2 BRAF-context paper(s) · large effect (|d|=2.24) · fully replicated (2/2 cohorts) |
| 6 | **PDLIM4** | 0.710 | +3.44 | +2.25 | 3.1e-40 | 2 | novel | no | PubMed-novel (<5 thyroid) · 2 BRAF-context paper(s) · large effect (|d|=2.25) · fully replicated (2/2 cohorts) |
| 7 | **CREB5** | 0.704 | +2.09 | +2.03 | 3.6e-49 | 2 | novel | no | PubMed-novel (<5 thyroid) · 3 BRAF-context paper(s) · large effect (|d|=2.03) · fully replicated (2/2 cohorts) |
| 8 | **LDLR** | 0.702 | +2.24 | +2.02 | 2.8e-53 | 2 | established | yes | v7-druggable · 3 BRAF-context paper(s) · large effect (|d|=2.02) · fully replicated (2/2 cohorts) |
| 9 | **CST6** | 0.665 | +3.87 | +2.04 | 1.2e-41 | 2 | novel | no | PubMed-novel (<5 thyroid) · 1 BRAF-context paper(s) · large effect (|d|=2.04) · fully replicated (2/2 cohorts) |
| 10 | **LY6E** | 0.659 | +2.13 | +2.16 | 5.8e-43 | 2 | novel | no | PubMed-novel (<5 thyroid) · large effect (|d|=2.16) · fully replicated (2/2 cohorts) |

## Second tier (0.50 ≤ score < 0.65)

| Rank | Gene | Score | Novelty | Drug? | Rationale |
|---:|---|---:|---|---|---|
| 11 | BNC1 | 0.635 | novel | no | PubMed-novel (<5 thyroid) · 2 BRAF-context paper(s) · large effect (|d|=1.63) · fully replicated (2/2 cohorts) |
| 12 | KCNQ3 | 0.635 | novel | no | PubMed-novel (<5 thyroid) · large effect (|d|=1.96) · fully replicated (2/2 cohorts) |
| 13 | SPOCK2 | 0.634 | novel | no | PubMed-novel (<5 thyroid) · 1 BRAF-context paper(s) · large effect (|d|=1.79) · fully replicated (2/2 cohorts) |
| 14 | DCSTAMP | 0.609 | emerging | no | emerging (5-24 thyroid) · 2 BRAF-context paper(s) · large effect (|d|=2.24) · fully replicated (2/2 cohorts) |
| 15 | ST6GALNAC5 | 0.605 | novel | no | PubMed-novel (<5 thyroid) · large effect (|d|=1.71) · fully replicated (2/2 cohorts) |
| 16 | ITGA3 | 0.597 | emerging | no | emerging (5-24 thyroid) · 2 BRAF-context paper(s) · large effect (|d|=2.14) · fully replicated (2/2 cohorts) |
| 17 | KCNN4 | 0.585 | emerging | no | emerging (5-24 thyroid) · large effect (|d|=2.37) · fully replicated (2/2 cohorts) |
| 18 | TAGLN2 | 0.562 | emerging | no | emerging (5-24 thyroid) · large effect (|d|=2.18) · fully replicated (2/2 cohorts) |
| 19 | BID | 0.555 | established | no | 44 BRAF-context paper(s) · large effect (|d|=2.13) · fully replicated (2/2 cohorts) |

## Third tier (score < 0.50)

| Rank | Gene | Score | Notes |
|---:|---|---:|---|
| 20 | SYT12 | 0.497 | emerging (5-24 thyroid) · large effect (|d|=2.47) · partially replicated (1/2 cohorts) |

## Top-5 detailed synopses for clinical-collaborator briefing

### 1. PLEKHA6  (SNUBH score 0.906)

- **Druggable**: yes (v7 chemistry support; target literature: 0 thyroid / 4 cancer / 2 BRAF papers; druggability flag = 1).
- **TCGA effect**: log2FC = +2.35, Cohen's d = +2.22, FDR = 3.0e-54; replicated in 2 of 2 external cohorts.
- **PubMed depth**: 0 thyroid · 2 BRAF-context — classification: novel.
- **Clinical**: 0 thyroid trial(s), 0 any-cancer trial(s); translatability tier = medium.
- **Recommended SNUBH validation modality**: targeted RT-qPCR + IHC pilot (n≈30 BRAF V600E vs 30 RAS-like) to confirm subtype-specific over-expression before any therapeutic-modality discussion.

### 2. PTPRE  (SNUBH score 0.891)

- **Druggable**: yes (v7 chemistry support; target literature: 4 thyroid / 14 cancer / 1 BRAF papers; druggability flag = 1).
- **TCGA effect**: log2FC = +2.38, Cohen's d = +2.26, FDR = 3.8e-45; replicated in 2 of 2 external cohorts.
- **PubMed depth**: 4 thyroid · 1 BRAF-context — classification: novel.
- **Clinical**: 0 thyroid trial(s), 0 any-cancer trial(s); translatability tier = medium.
- **Recommended SNUBH validation modality**: targeted RT-qPCR + IHC pilot (n≈30 BRAF V600E vs 30 RAS-like) to confirm subtype-specific over-expression before any therapeutic-modality discussion.

### 3. CYP1B1  (SNUBH score 0.846)

- **Druggable**: yes (v7 chemistry support; target literature: 11 thyroid / 450 cancer / 4 BRAF papers; druggability flag = 1).
- **TCGA effect**: log2FC = +3.47, Cohen's d = +2.22, FDR = 9.0e-53; replicated in 2 of 2 external cohorts.
- **PubMed depth**: 16 thyroid · 4 BRAF-context — classification: emerging.
- **Clinical**: 0 thyroid trial(s), 0 any-cancer trial(s); translatability tier = medium.
- **Recommended SNUBH validation modality**: targeted RT-qPCR + IHC pilot (n≈30 BRAF V600E vs 30 RAS-like) to confirm subtype-specific over-expression before any therapeutic-modality discussion.

### 4. B3GNT3  (SNUBH score 0.842)

- **Druggable**: yes (v7 chemistry support; target literature: 0 thyroid / 42 cancer / 0 BRAF papers; druggability flag = 1).
- **TCGA effect**: log2FC = +3.78, Cohen's d = +2.02, FDR = 8.0e-52; replicated in 2 of 2 external cohorts.
- **PubMed depth**: 0 thyroid · 0 BRAF-context — classification: novel.
- **Clinical**: 0 thyroid trial(s), 0 any-cancer trial(s); translatability tier = medium.
- **Recommended SNUBH validation modality**: targeted RT-qPCR + IHC pilot (n≈30 BRAF V600E vs 30 RAS-like) to confirm subtype-specific over-expression before any therapeutic-modality discussion.

### 5. GABRB2  (SNUBH score 0.809)

- **Druggable**: yes (v7 chemistry support; target literature: 9 thyroid / 18 cancer / 2 BRAF papers; druggability flag = 1).
- **TCGA effect**: log2FC = +3.70, Cohen's d = +2.24, FDR = 2.9e-47; replicated in 2 of 2 external cohorts.
- **PubMed depth**: 9 thyroid · 2 BRAF-context — classification: emerging.
- **Clinical**: 0 thyroid trial(s), 0 any-cancer trial(s); translatability tier = medium.
- **Recommended SNUBH validation modality**: targeted RT-qPCR + IHC pilot (n≈30 BRAF V600E vs 30 RAS-like) to confirm subtype-specific over-expression before any therapeutic-modality discussion.

