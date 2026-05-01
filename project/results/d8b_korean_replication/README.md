---
dir: results/d8b_korean_replication
pillar: Paper 2 Pillar III — Korean replication arm (GSE213647 Lee 2024)
status: STRONG
generator: notebooks_or_scripts/v17_D8B_korean_replication.py
---

# d8b_korean_replication — Paper 2 Pillar III Korean arm

GSE286332 PTC+HT-vs-PTC top DEG signature → Korean Lee 2024 GSE213647 n=632 transfer. **Cross-cohort generalization confirmation** (TCGA + Korean 두 독립 코호트).

## Files

| File | Description | Used by |
|---|---|---|
| `D8B_summary.json` | Top-level: per-cutoff Hashimoto-like prevalence Korean arm | Brief Fig 8 |
| `korean_GSE213647_hashimoto_scores.tsv` | Per-sample (n=632) signature score + binary calls | Pillar III Korean input |

## Brief usage cross-reference

| Brief Fig | Source |
|---|---|
| Fig 8 — Hashimoto-like prevalence × 4 thresholds (Korean arm) | `D8B_summary.json` |

## Key claims

- Korean GSE213647 Hashimoto-like prevalence: **GMM 22.8%, Otsu 28.2%, top 20% 20.1%, top 30% 30.1%**
- TCGA 18-30% 와 일관 — 두 독립 한국인 코호트에서도 generalize
- Bimodality coefficient Korean=0.471 — bimodal (axis 분리 가능)

## Companion

- TCGA discovery arm: `../d4p2_tcga_hashimoto_signature/`
- Korean GSE213647 raw: GEO accession GSE213647, processed via `notebooks_or_scripts/v17_ULTIMATE_U1C_gse213647.py`
