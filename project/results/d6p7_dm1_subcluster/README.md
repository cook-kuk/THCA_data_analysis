---
dir: results/d6p7_dm1_subcluster
pillar: Paper 2 Pillar V — DM1 sub-B = NBNR cluster (Korean K2 bridge)
status: PARTIAL
generator: notebooks_or_scripts/v17_D6P7_dm1_subcluster.py
---

# d6p7_dm1_subcluster — Paper 2 Pillar V DM1 internal sub-cluster

TCGA DM1 (n=140) 내부 KMeans k=2 on TIERA67 expression → **sub-A n=84 (60%, 69% RAS+)** + **sub-B n=56 (40%, 96% mutation-negative)**. NBNR (No-BRAF-No-RAS) sub-cluster identification.

## Files

| File | Description | Used by |
|---|---|---|
| `D6P7_summary.json` | Top-level: cluster sizes, mutation breakdown, Hashimoto-like rate per sub | Pillar V anchor |
| `dm1_subcluster_labels.tsv` | Per-sample DM1 sub-A / sub-B label + mutation status | Brief **Fig 14 ★★★** |
| `dm1_subBvA_deg.tsv` | sub-B vs sub-A DEG list (PyDESeq2) | Pillar V signature for `../d8c_dm1_subB_x_K2_NBNR/` transfer |
| `subcluster_clinical.tsv` | TCGA clinical + sub-cluster join | Pillar V Hashimoto-like 12.5% (sub-B) vs 3.6% (sub-A) |
| `subcluster_score_profile.tsv` | TIERA67 expression z-score profile per sub-cluster | Score-based transfer |
| `subcluster_scores.tsv` | Per-sample sub-A vs sub-B signature score | Korean transfer input |

## Brief usage cross-reference

| Brief Fig | Source |
|---|---|
| **Fig 14 ★★★** — DM1 sub-A vs sub-B mutation landscape | `dm1_subcluster_labels.tsv` (96% mutation-negative paper-defining) |
| Fig 15 — Korean GSE213647 sub-B-like rate (companion) | `../d8c_dm1_subB_x_K2_NBNR/korean_subB_score.tsv` |

## Key claims

- DM1 (BRAF/RAS-positive cluster) 내부에 **96% mutation-negative sub-cluster** 존재 (unsupervised KMeans 가 mutation status 정보 없이 발견)
- Sub-A: 84명 중 RAS+ 51 (69%), BRAF+ 1 — RAS+ FVPTC core
- Sub-B: 56명 중 RAS+ 2, BRAF+ 1, mutation-negative 47/49 (96%) — **NBNR cluster**
- Sub-B Hashimoto-like 12.5% vs sub-A 3.6%, Fisher OR 0.26, p=0.09 (trend)
- 유 교수님 K2 NBNR ETE-aggressive cohort phenotype 의 **TCGA-side molecular equivalent**

## Cross-paper boundary

- Mutation landscape (sub-A vs sub-B) → Paper 1 (Cancer paper, DM1 dark matter)
- Hashimoto convergence (sub-B Hashimoto-like 12.5%) → **Paper 2 Pillar V** (이 brief)
- → Discussion Q04 advisor 결정 의제: 본인 default = sub-B mutation = Paper 1, sub-B Hashimoto = Paper 2
