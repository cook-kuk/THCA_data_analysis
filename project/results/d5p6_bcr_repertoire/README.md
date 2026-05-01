---
dir: results/d5p6_bcr_repertoire
pillar: Paper 2 Pillar IV — Antigen-driven specificity (BCR clonal + TLS + AICDA)
status: STRONG
generator: notebooks_or_scripts/v17_D5P6_bcr_repertoire.py
---

# d5p6_bcr_repertoire — Paper 2 Pillar IV antigen-driven evidence

GSE286332 raw FASTQ → IGHV/IGKV/IGLV V-gene 동정 + Shannon entropy + clonality (1 − normalized entropy) + Cabrita 2020 12-gene TLS signature + AICDA expression.

## Files

| File | Description | Used by |
|---|---|---|
| `D5P6_summary.json` | Top-level: n_IG genes per sample, group_comparison Cohen d + p, AICDA, TLS, spearman_clo_g8_tls | Brief Fig 13 (spearman matrix data) |
| `group_comparison.tsv` | TLS / IGHV / IGKV clonality / entropy — Cohen's d (PTC+HT vs PTC) | Brief **Fig 12 ★★★** |
| `per_sample_diversity.tsv` | Per-sample IGHV/IGKV/IGLV entropy + clonality | Pillar IV input |
| `tls_score_per_sample.tsv` | Per-sample Cabrita 2020 12-gene TLS signature | Pillar IV input |

## Brief usage cross-reference

| Brief Fig | Source |
|---|---|
| **Fig 12 ★★★** — TLS d=+1.96, IGHV clonality d=+2.09 | `group_comparison.tsv` |
| Fig 13 — Spearman matrix (8-gene × IGHV × TLS × HLA-II) | `D5P6_summary.json:spearman_clo_g8_tls` (key embeds matrix) |

## Key claims

- TLS Cabrita 2020 Cohen's d = **+1.96** (massive, in PTC+HT)
- IGHV clonality d = **+2.09** (p=1.1e−3) — antigen-specific clonal expansion direct signal
- IGKV clonality d = +1.41 (p=0.034)
- AICDA up — somatic hypermutation activation
- 8-gene RAI × IGHV clonality ρ=−0.67 (p=0.002), TLS ρ=−0.79 (p=1e−4)
- IGHV clonality × TLS ρ=+0.82 (p=3e−5) — TLS = clonal expansion spatial substrate

## Mechanism interpretation

Generic immune infiltration 만으로는 IGHV total ↑ 만 보임. **Clonality** (normalized entropy 의 inverse) 와 AICDA 와 TLS 가 함께 상승해야 antigen-driven response. 본 dir 의 4-way evidence (clonality + AICDA + TLS + IGHV total) 가 단순 면역 침윤 가설을 직접 기각.
