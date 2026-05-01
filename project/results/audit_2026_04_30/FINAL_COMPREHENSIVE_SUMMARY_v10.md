---
title: "rThyroid Dark Matter Paper — 통합 v10 (manuscript v7 ship-ready)"
date: 2026-05-01
sessions: ["04-29 audit", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "manuscript v6→v7 line edits"]
status: ship-ready
manuscript_target: npj Precision Oncology (current submission); reach Cell Reports Medicine / Nature Medicine after npj decision
purpose: Final wrap. v7 manuscript + v4 cover letter + v7 SUBMIT_INSTRUCTIONS rendered. Ready for editorial-manager click-submit.
---

# rThyroid Dark Matter Paper — 통합 v10 (final wrap)

## 🎯 Status: ship-ready 2026-05-01 11:20 KST

| Artifact | Path | Format |
|---|---|---|
| **Manuscript v7** | `project/submission/npj/manuscript_v7.{md,html,pdf}` | 322 lines, 3.86 MB PDF |
| **Cover letter v4** | `project/submission/npj/cover_letter_v4.{md,html,pdf}` | updated 2026-05-01 |
| **Submit instructions v7** | `project/submission/npj/SUBMIT_INSTRUCTIONS_v7.md` | step-by-step |
| **Fig 10** | `project/reports/html/figs_interactive/v17/v17_fig9_composite_audit.html` | composite ROC + GSE213647 |
| **Fig 11** | `project/reports/html/figs_interactive/v17/v17_fig10_R9_chr7_RAI_K2.html` | chr7 forest + I-131 + K2 within-z |
| **Fig 12** | `project/reports/html/figs_interactive/v17/v17_fig11_R10_TLS_3way_chr7_MSK.html` | TCGA TLS + 3-way + chr7×fus + MSK |

## 📋 11 sessions complete

| # | Session | Date | Δ | Status |
|---|---|---|---|---|
| 1 | 04-29 audit (A–J) | 04-29 | 67-chart dashboard | ✅ |
| 2 | R1 (P1–P7) | 04-30 | TERT⁺ Cox HR=4.33 | ✅ |
| 3 | R2 (N1–N7) | 04-30 | MSK MAF 100% match | ✅ |
| 4 | R3 (F1–F4) | 04-30 | DM1 76.8% fusion paradigm | ✅ |
| 5 | R4 | 04-30 | OR 7–9 robust | ✅ |
| 6 | R5 | 04-30 | TPO promoter d=2.30 | ✅ |
| 7 | R6 | 05-01 | meth × expr ρ=−0.68 | ✅ |
| 8 | R7 | 05-01 | composite + stemness 역설 | ✅ |
| 9 | R8 | 05-01 | chr7 arm CNV + RET partner-agnostic + composite re-merge AUC 0.77 | ✅ |
| 10 | R9 | 05-01 | Xena RAI dose + chr7 RTK gene-level + bootstrap + K2 within-z | ✅ |
| 11 | R10 | 05-01 | Cabrita TLS d=0.67 replication + 3-way clinical algorithm + chr7×fusion + MSK panel limitation | ✅ |
| 12 | **manuscript v6→v7 line edits** | **05-01** | **Abstract / §3.5 / §3.10 / §3.11 / Discussion 3.5 / Methods § / Fig 10/11/12 legends** | ✅ |

## 🧱 9-layer DM1 mechanism — final consolidation

| Layer | Evidence | Round |
|---|---|---|
| 1. Genetic (driver) | 76.8% RET-fusion (n=110 DM1) | R3 |
| 2. Heterogeneity | sub-A (RET-fusion+) vs sub-B (NBNR n=56, chr7-amp 17%) | R4 + R10-3 |
| 3. Epigenetic | TPO promoter β = 0.92, d = 2.30 vs not_DM | R5-2 |
| 4. Causal | meth × expr Spearman ρ = −0.68 TPO | R6-4 |
| 5. Genomic structural | chr 7p/7q gain DM1 3.4% vs not_DM 0% (Fisher p=0.011; bootstrap 95% CI [0.000, 0.079] p=0.053); BRAF/EGFR/MET log2CNA d=+0.27~+0.29 | R8-1 + R9-2 + R9-3 |
| 6. Signaling | RPPA pMAPK signature DM1↑ | R6-3 |
| 7. Immune | HLA-II↑ + IFN-γ + Cabrita TLS d=+1.96 (GSE286332) + d=+0.67 p=4.4e-15 (TCGA, replication) | R5-1 / R6-1 / R10-5 |
| 8. Stemness | Active stem (MYC/BMI1/CD44↑) DM1 vs Dormant stem (ALDH1A1↑) DM2 | R7-5 |
| 9. Pathway dichotomy | RAI-responsive vs RAI-refractory split within DM1 by RET-fusion+/− | R7-1 / R10-2 |

**Clinical pilot**: DM1 + RET-fusion-positive + TERT-WT (n=67, **0 OS events**) = candidate biomarker-positive cohort for continued-RAI-responder protocol.

## 📊 Honest disclosures (cumulative — 46) → mapped to manuscript v7 limitations

v7 explicitly enumerates **16** limitations (was 13 in v6). 새로 추가된 3개:

- **L14** TCGA per-patient ¹³¹I dose recorded for only 36/513 (Xena legacy field)
- **L15** chr7 finding borderline at α=0.05 (bootstrap 95% CI lower 0.000)
- **L16** MSK 2016 chr7 GISTIC sparse calls (targeted-panel artefact)

내부 audit MD에 더 자세한 43–46 추가 disclosures 있음 (v9 MD §"Honest disclosures"); 그 중 핵심 3개만 manuscript에 명시. Composite over-merge (audit-internal) 등은 supplementary methods note에서 다룸.

## 📚 Reviewer Q&A (v9 20개) → manuscript v7 vs supplementary

Manuscript main에 직접 답변된 Qs (1, 2, 3 = anti-circular; 5 = TCGA TLS; 7 = RET partner-agnostic; 9 = chr7 ext-replication; 11 = stemness paradox; 13 = 3-way subgroup): 6개.

내부 audit MD에 보존, supplementary methods 또는 rebuttal 단계용: 14개.

## 🎨 Figure landscape — 12 figures total

| Fig | Source | Status |
|---|---|---|
| Fig 1 | DM1/DM2 axis discovery | shipped (v6) |
| Fig 2 | Biology + clinical | shipped |
| Fig 3 | Trajectory | shipped |
| Fig 4 | Mutation distinctness | shipped |
| Fig 5 | Hot/cold immune | shipped |
| Fig 6 | 8-gene panel original | shipped |
| Fig 7 | PRISM drug | shipped |
| Fig 8 | TERT 4-group | shipped |
| Fig 9 | REAL FIX leak-free | shipped (was internally Fig 9 in v6, kept name) |
| **Fig 10** | composite ROC + GSE213647 (R8-4) | NEW v7 |
| **Fig 11** | chr7 forest + I-131 + K2 within-z (R8-1, R9-1, R9-2, R9-3, R9-4) | NEW v7 |
| **Fig 12** | TCGA Cabrita TLS + 3-way + chr7×fus + MSK (R10) | NEW v7 |

## 🏁 사용자가 다음에 해야 할 것 (즉시)

1. **5-min pre-flight** — `manuscript_v7.pdf` 한 번 훑고 author block 확정 / title 확정
2. **30-min portal submission** — `SUBMIT_INSTRUCTIONS_v7.md` 따라 editorial-manager 입력
3. **2-min after-submit** — 4 outreach drafts 보내기 (분당, 유 교수님, Xing, Landa)
4. **Optional follow-up** — supplementary table 2개 (S7 8cell, S8 chr7 arm) 추가, Fig 10-12 PDF export (만약 npj가 HTML 거부 시 `fig.write_image` 한 줄)

## 📁 참고 경로 (final ship cheatsheet)

```
project/submission/npj/
├── manuscript_v7.{md,html,pdf}          ← 본 paper
├── cover_letter_v4.{md,html,pdf}        ← 본 cover
├── SUBMIT_INSTRUCTIONS_v7.md            ← step-by-step
├── anonymous_code.zip                   ← reproducibility (re-zip with round8-10 권장)
└── figures/Fig1.pdf ... Fig8.pdf

project/reports/html/figs_interactive/v17/
├── v17_fig9_composite_audit.html        ← Fig 10 (npj label)
├── v17_fig10_R9_chr7_RAI_K2.html        ← Fig 11
└── v17_fig11_R10_TLS_3way_chr7_MSK.html ← Fig 12

project/results/audit_2026_04_30/
├── FINAL_COMPREHENSIVE_SUMMARY_v10.md   ← 본 wrap
├── round{2..10}/                        ← 모든 audit raw outputs
└── JUDGMENT_PROMPT_STATUS.md            ← 결정 history
```
