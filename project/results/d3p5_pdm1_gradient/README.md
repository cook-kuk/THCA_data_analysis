---
dir: results/d3p5_pdm1_gradient
pillar: Paper 2 § 3.6 Mediation (causal arrow)
status: STRONG
generator: notebooks_or_scripts/v17_D3P5_pdm1_gradient.py
---

# d3p5_pdm1_gradient — Paper 2 § 3.6 Mediation

Baron-Kenny 5,000-iter bootstrap mediation analysis: PTC+HT → mediator (HLA-II / 8-gene RAI / HLA-I / generic immune) → P(DM1) ↓. **HLA-II 가 dominant mediator (140%)**.

## Files

| File | Description | Used by |
|---|---|---|
| `D3P5_summary.json` | Top-level: mediation per mediator + within-PTC severity gradient | § 3.6 anchor |
| `mediation_results.json` | Per-mediator a/b/indirect/% mediated/bootstrap CI/p_emp | Brief **Fig 16 ★★★** + Table 1 |
| `decomposition.json` | Indirect / direct / total effect decomposition | § 3.6 supplementary |
| `merged_18sample.tsv` | GSE286332 n=18 with all mediator columns + outcome | Mediation input |
| `spearman_pdm1_vs_covariates.tsv` | ρ(P_DM1, each covariate) | Mediator selection |
| `within_ptc_severity.json` | PTC-only n=9 ρ(g8_RAI, P_DM1) = +0.73 (p=0.025) | Brief **Fig 17** |

## Brief usage cross-reference

| Brief Fig / Table | Source |
|---|---|
| **Fig 16 ★★★** — Mediation % (HLA-II 140%, generic immune 88% NS) | `mediation_results.json` |
| Fig 17 — Within-PTC severity gradient (pre-clinical Hashimoto-like spectrum) | `within_ptc_severity.json` |
| Table 1 — Mediation summary 4 mediators | `mediation_results.json` |

## Key claims

- **HLA-II 140% mediated** (p_emp=0.023, 95% CI [−0.31, −0.03]) — over-mediation, dominant
- 8-gene RAI 63% partial (p=0.002)
- HLA-I 87% (p=0.042)
- Generic immune 88% **NS** (p=0.120) — confounder, not mediator
- → **HLA-II 가 PTC+HT → P(DM1) ↓ 경로의 dominant mediator**, generic immune-infiltration 가설 직접 기각
- Within-PTC 9명 ρ(g8_RAI, P_DM1) = +0.73 (p=0.025) — pre-clinical Hashimoto-like spectrum 신호

## Paradox 해소

GSE286332 18/18 환자가 모두 DM2 호출 (categorical) → P_DM1 spectrum (continuous 0.005–0.304) 이 PTC+HT 와 negatively correlated 로 mediation analysis 가 categorical label 의 함정을 우회.
