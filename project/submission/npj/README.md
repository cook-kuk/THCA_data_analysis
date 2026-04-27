# v17 npj submission package — reviewer orientation

This folder is the complete submission bundle for "An 8-gene RAI-responsiveness biomarker outperforms BRAF V600E status in papillary thyroid carcinoma" (target: _npj Precision Oncology_, Article).

## ⚡ v4 BOOST sprint update (2026-04-27 12:00)

Three reviewer-anticipated robustness layers added to the v4 manuscript:

1. **4-cohort meta-analysis** — pooled AUC 0.980 (95% CI 0.869–0.997, **I² = 0%**) across GSE27155 / GSE29265 / GSE33630 / GSE76039 (n = 290).
2. **Decision Curve Analysis** — 8-gene LogReg dominant at **all 91 thresholds** (0.05–0.95), the strongest possible Vickers–Elkin outcome.
3. **9-strata subgroup forest** — 7/9 strata AUC ≥ 0.85; Stage III/IV peaks at **AUC 0.996**.

Plus a browser-based **interactive RAI calculator** (no server, no install) at `../reports/v17_boost/rai_calculator.html`.

## Where to start

| File | What it is |
|------|-----------|
| **`manuscript_v4.pdf`** | **Read this first.** Latest manuscript with BOOST validation cascade. MD/HTML/DOCX siblings. |
| `manuscript_v1.pdf` | Original v1 manuscript (kept for diff). |
| **`cover_letter_v3.pdf`** | Cover letter v3 — adds BOOST layers + interactive tool URL. |
| `cover_letter.pdf` | Original v1 cover letter (kept for diff). |
| `figures/master_panel.html` | Browser-friendly index of all 26 figures (8 main + 15 supp + 3 BOOST 6E/6F/6G + interactive calculator card). |
| `figures/figure6E_dca.{png,pdf}` | DCA curve, 8-gene vs BRAF vs treat-all/none. |
| `figures/figure6F_multi_cohort_forest.{png,pdf}` | 4-cohort meta forest. |
| `figures/figure6G_subgroup_forest.{png,pdf}` | 9-strata subgroup forest. |
| `../reports/v17_boost/rai_calculator.html` | Interactive 8-gene RAI score calculator. |
| `../reports/v17_boost/index.html` | Korean review dashboard (한 페이지 요약 for 유 교수님). |
| `tables/Supplementary_Tables.xlsx` | All 7 supplementary tables. |
| `v17p35_REVIEWER_DEFENSE.md` | Internal: 20 anticipated reviewer attacks with rebuttals. |
| `SUBMIT_INSTRUCTIONS.md` | Step-by-step npj Editorial Manager submission walkthrough. |
| `SUBMISSION_CHECKLIST.md` | Pre-flight check. |
| `anonymous_code.zip` | Code archive with author identifiers stripped. |
| `reviewer_bundle.zip` | Single-shot upload bundle (will be regenerated to v4). |

## Headline numbers (v4)

- 8-gene panel CV AUC **0.954** vs BRAF V600E baseline **0.822** (ΔAUC = **+0.132**)
- 4-cohort meta-analysis pooled AUC **0.980** (95% CI 0.869–0.997, I² = **0%**)
- Decision Curve Analysis: 8-gene dominant at **all 91 thresholds** (0.05–0.95)
- 9-strata subgroup forest: **7/9 strata AUC ≥ 0.85**, Stage III/IV peak AUC **0.996**
- Hot/Cold composite Cohen's d = **+1.683** (Mann-Whitney p = 4.0×10⁻¹⁸)
- TERT 4-group survival logrank p = **3.78×10⁻⁵** (n = 504); univariate Cox HR 6.31 → **multivariate HR 1.88, p = 0.29 (honest)**
- External GSE76039 validation AUC **0.974** (correct-direction)
- 4 / 4 BRAF V600E CCLE thyroid cell lines correctly DM1

## Layout

```
submission/npj/
├── manuscript_v4.{md,html,pdf,docx}      ← latest, BOOST-validated
├── manuscript_v1.{md,html,pdf,docx}      ← original (for diff)
├── cover_letter_v3.{md,html,pdf,docx}    ← latest, mentions BOOST + calculator
├── cover_letter.{md,html,pdf,docx}       ← original (for diff)
├── figures/
│   ├── Fig1..8 + SuppFig1..15 (existing 23)
│   ├── figure6E_dca.{png,pdf}            ← NEW
│   ├── figure6F_multi_cohort_forest.{png,pdf}  ← NEW
│   ├── figure6G_subgroup_forest.{png,pdf}      ← NEW
│   └── master_panel.html                 ← updated, 26 figs
├── tables/
├── reviewer_bundle.zip                   ← regenerated to v4
├── anonymous_code.zip
└── ../reports/v17_boost/
    ├── rai_calculator.html               ← interactive calculator
    └── index.html                        ← Korean review dashboard for PI
```

## Reproducibility

Pipeline scripts at `../../notebooks_or_scripts/v17p35_*.py`, `../../notebooks_or_scripts/v17_FINAL_*.py`, and `../../notebooks_or_scripts/v17_BOOST_*.py`. Random seeds (`random_state=42`, `np.random.seed(42)`) throughout.
