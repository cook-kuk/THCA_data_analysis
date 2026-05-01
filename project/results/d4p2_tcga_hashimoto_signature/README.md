---
dir: results/d4p2_tcga_hashimoto_signature
pillar: Paper 2 Pillar III — Cross-cohort generalization (TCGA arm)
status: STRONG
generator: notebooks_or_scripts/v17_D4P2_tcga_hashimoto_signature.py
---

# d4p2_tcga_hashimoto_signature — Paper 2 Pillar III TCGA arm

GSE286332 PTC+HT-vs-PTC top DEG (padj<0.01, |LFC|>1) signature transferred to TCGA-THCA n=500. **DM1/DM2 × Hashimoto-like cross-tab paper-defining anchor**.

## Files

| File | Description | Used by |
|---|---|---|
| `D4P2_summary.json` | Top-level: prevalence × 4 thresholds, crosstabs (DM1/DM2 × Hashimoto-like for each cutoff), residualized OR, HLA-II saturation | Brief Fig 8, 9, 10, 11 |
| `tcga_signature_scores.tsv` | Per-sample (n=500) z-mean signature score (143 up + 46 down) | Pillar III input |
| `tcga_with_clinical_mutations.tsv` | TCGA-THCA clinical + mutation joined for crosstab | DM1/DM2 cross-tab + residualization |

## Brief usage cross-reference

| Brief Fig | Source key |
|---|---|
| Fig 8 — Hashimoto-like prevalence × 4 thresholds | `D4P2_summary.json:prevalence_4cutoffs` |
| **Fig 9 ★★★** — TCGA DM1/DM2 × Hashimoto-like (OR=0.20, p=6.4e−10) | `D4P2_summary.json:crosstabs.hashi_top30` |
| Fig 10 — Stromal+immune residualized OR | `D4P2_summary.json:crosstabs.hashi_resid_otsu` |
| Fig 11 — HLA-II saturation analysis (DM1/DM2 × Hashimoto subset) | `D4P2_summary.json:hla_residualization_*` |

## Key claims

- 4 cutoff methods (GMM, Otsu, top 20%, top 30%) → Hashimoto-like prevalence 18-30% TCGA-side (robust)
- DM1 환자 중 Hashimoto-like 10.7% vs DM2 37.5% → **Fisher OR 0.20, p = 6.4×10⁻¹⁰**
- Stromal+immune residualized: OR 0.29, p=8e−9 (specific axis, not generic immune confounder)
- HLA-II saturation: full TCGA d=−1.41, excl Hashimoto+ d=−1.60 (stronger), within Hashimoto+ d=+0.14 NS → two pathways

## Companion

- Korean replication arm: `../d8b_korean_replication/` (GSE213647 n=632)
- GSE286332 source signature: `../p3_gse286332/deg_ptcht_vs_ptc.tsv` (top 143 up + 46 down)
