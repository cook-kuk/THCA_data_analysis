---
dir: results/d8c_dm1_subB_x_K2_NBNR
pillar: Paper 2 Pillar V — DM1 sub-B Korean transfer (GSE213647)
status: PARTIAL
generator: notebooks_or_scripts/v17_D8C_dm1_subB_x_K2_NBNR.py
---

# d8c_dm1_subB_x_K2_NBNR — Paper 2 Pillar V Korean transfer

TCGA DM1 sub-B signature (`../d6p7_dm1_subcluster/dm1_subBvA_deg.tsv`) → Korean Lee 2024 GSE213647 n=632 transfer. **K2 NBNR phenotype 의 generalizable axis** 검정.

## Files

| File | Description | Used by |
|---|---|---|
| `D8C_summary.json` | Top-level: per-cutoff sub-B-like rate, Hashimoto convergence | Pillar V Korean arm |
| `korean_subB_score.tsv` | Per-sample (n=632) sub-B signature score + GMM/Otsu binary call | Brief **Fig 15** |
| `korean_subB_x_hashimoto.tsv` | Korean sub-B-like × Hashimoto-like cross-tab | Pillar V replication |

## Brief usage cross-reference

| Brief Fig | Source |
|---|---|
| Fig 15 — Korean GSE213647 sub-B-like rate (47-53%) | `korean_subB_score.tsv` |

## Key claims

- Korean GSE213647 n=632 sub-B-like rate: **GMM 47.2%, Otsu 52.5%** — 한국 일반 PTC 의 약 50% 가 sub-B-like (NBNR-like)
- TCGA reference 40% (sub-B / DM1 = 56/140) — Korean rate higher
- 유 교수님 K2 NBNR ETE-aggressive cohort phenotype 의 generalizable axis

## Companion

- TCGA discovery arm: `../d6p7_dm1_subcluster/` (sub-A vs sub-B mutation landscape)
- Korean GSE213647 raw: GEO accession GSE213647 (processed via `notebooks_or_scripts/v17_ULTIMATE_U1C_gse213647.py`)
- K2 raw: ENA PRJEB11591 (Yoo 2016 SNU-GMI)
