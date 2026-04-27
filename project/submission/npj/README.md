# v17 npj submission package — reviewer orientation

This folder is the complete submission bundle for **"An 8-gene RAI-responsiveness biomarker outperforms BRAF V600E status in papillary thyroid carcinoma"** (target: *npj Precision Oncology*, Article).

Last updated: 2026-04-27 (v17 KOREAN sprint K2 + Fig 8 Cox-label fix).

## Where to start

If you are a **reviewer** or **PI reading this for the first time**, open the bundled web dashboard:

| Entry point | What it gives you |
|------|-----------|
| **`http://40.82.129.113:8765/`** (or open `index.html` locally) | Landing page → links to every artifact below. |
| `START_HERE.html` | One-page overview for new readers. |
| `easy_explainer.html` | Plain-Korean explanation of DM1/DM2 + the panel (no jargon). |
| `reviewer_faq.html` | 20+ anticipated reviewer questions with rebuttals. |
| `analyses_done.html` / `data_models_performance.html` | Granular log of what was run + every metric. |
| `glossary.html` | Defined terms reference. |
| `korean_dashboard_v3.html` | v17 KOREAN sprint — PRJEB11591 K2 pilot results, Fig K2 + 3D, external-data request detail tables, 4 outreach drafts. |

## Active manuscript

| File | Notes |
|------|-------|
| **`manuscript_v6.{md,html,pdf,docx}`** | **Current submission text.** ~12.9 K words. Korean cohort initial integration (PRJEB11591, n = 9 pilot) included; Fig 8 Cox HRs synced to lifelines computation (univariate joint 4-group HR = 4.33 → multivariate stage·age·sex HR = 0.95). |
| `manuscript_v6_ULTIMATE.md` | Same content variant in long form. |
| **`cover_letter_v6_ULTIMATE.{md,html,pdf,docx}`** | Latest cover letter. |
| `manuscript_v{1,4,5}.*` + `cover_letter_v{3,4,5}.*` | Earlier versions kept for diff. |

## Headline numbers

- **8-gene panel CV AUC 0.962** (95% CI 0.940–0.979) vs BRAF V600E baseline 0.849 → **ΔAUC = +0.113** (95% CI excludes 0).
- **4-cohort meta-analysis pooled AUC 0.980** (I² = 0%, all four AUC > 0.96), all 4 leak-free panel variants tested.
- **Decision Curve Analysis**: 8-gene strategy dominant across 0.05–0.95 threshold range.
- **Subgroup forest** (9 strata): 7/9 AUC ≥ 0.85; Stage III/IV peaks at AUC 0.996.
- **Hot/Cold composite Cohen's d = +1.683** (Mann-Whitney p = 4.0 × 10⁻¹⁸).
- **TERT 4-group survival** logrank p = 3.78 × 10⁻⁵ (n = 504); joint 4-group Cox univariate HR = 4.33 → **multivariate (stage + age + sex) HR = 0.95, p = 0.95** — honestly reframed as advanced-stage molecular handle, not stage-independent prognostic.
- **External GSE76039** validation AUC = 0.935 (leak-free re-trained model).
- **Korean cohort initial integration (PRJEB11591, Yoo 2016, n = 9 pilot)**: kallisto single-end + 8-gene mini-index + scale-invariant within-sample-centered LogReg → 9/9 DM2 (mean p = 0.898). Full 262-run + ground-truth-labeled AUC committed for revision round.

## Figures

42 PNG / 41 PDF / 41 HTML interactive figures in `figures/`. Highlights:

- `Fig1.{png,pdf,html}` — DM1/DM2 axis discovery (panel b: BRS-surrogate × DM1/DM2 confusion, BRAF-like × DM1 = 365/93%, RAS-like × DM2 = 78/70%).
- `Fig2`–`Fig7` — DM1/DM2 phenotype, panel CV, immune integration, GSE76039 transfer, DCA, subgroup forest, alt-panels, BHT-101 anchoring.
- `Fig5.{png,pdf}` — Hot/Cold immune (panel a: GSEA NES bars properly populated; legend grouped by cell types vs DM clusters).
- `Fig8.{png,pdf}` — TERT 4-group survival; **panel b shows univariate AND multivariate Cox HR side-by-side per group**, marker shape distinguishes models.
- `Korean_K2.{png,pdf}` — within-sample-centered 8-gene heatmap (9 Korean + TCGA centroids) + per-sample p(DM2) bars.
- `Korean_K2_3D.html` — Plotly Scatter3d, rotatable: TCGA DM1/DM2 + 9 Korean projected onto same panel-PCA space.
- `SuppFig1`–`SuppFig15` — robustness, calibration, scRNA, alt-panels.
- `master_panel.html` — browse all figures in one page.

## Data + tables

- `tables/Supplementary_Tables.xlsx` — all supplementary tables.
- `results/v17_realfix/` — leak-free R1/R2/R3/R4 outputs (cluster labels, AUCs, predictions).
- `results/v17_korean/` — K2 outputs: `K2_8gene_tpm_matrix_v4.tsv`, `K2_korean_predictions_v4.tsv`, `K2_korean_summary_v4.json`, `K1A_prjeb11591_runs.tsv` (262 runs).
- `results/v17_ultimate/` — U1A 4-variant table, U2D K-metrics, U3B Yang/TCGA/MSK comparison, U3C supplementary table S20.

## PDF + HTML viewer exports

`pdf_exports/` (47 MB, listed in `.gitignore`) — every `.md` / `.json` / `.tsv` file in the package has a paired `.pdf` (printable) and `.html` (browser viewer with toolbar: 📥 download original, 📕 download PDF, search filter for tables ≥ 10 rows). Regenerate via:

```bash
cd /opt/thyroid-dash/project && .venv/bin/python notebooks_or_scripts/v17_md_to_pdf.py
```

Korean text is rendered via NanumGothic (in `~/.local/share/fonts/`); the converter chain is `pandoc → html5 → weasyprint`.

## Outreach

Active outreach drafts (NOT sent — author dispatches manually):

| File | Target | Ask |
|------|--------|-----|
| `outreach/K3A_EGA_DAR_application.md` | EGA DAC for `EGAD00001004845` (Yoo 2019 advanced Korean) | Data Access Request |
| `outreach/K3B_email_park_yj.md` | Park YJ Group (SNU) | Han 2023 BL/RL framework alignment + cohort discussion |
| `outreach/K3C_email_bundang_v3.md` | Yu Hyeong Won (분당SNUH) | Contemporary Korean cohort (n = 50–100) |
| `outreach/K3D_email_yoo_sk.md` | Yoo SK (PRJEB11591 first author) | PRJEB11591 sample-level metadata + DAR endorsement |

## Reproducibility

- Random seed: `random_state=42`, `np.random.seed(42)` throughout.
- Figure-generation scripts: `../../notebooks_or_scripts/v17p35_SYNTH1_Fig*.py`, `../../notebooks_or_scripts/v17_FINAL_F1B_fig8.py`, `../../notebooks_or_scripts/v17_KOREAN_K2_*.py`.
- Korean K2 pipeline: `v17_KOREAN_K2_v4_single_end.py` (kallisto quant) → `v17_KOREAN_K2_v4_fix.py` (within-sample-centered LogReg) → `v17_KOREAN_K2_figure.py` (matplotlib) + `v17_KOREAN_K2_3d_umap.py` (Plotly 3D).
- PDF/HTML viewer pipeline: `v17_md_to_pdf.py`.

## Submission instructions

`SUBMIT_INSTRUCTIONS.md` + `SUBMIT_INSTRUCTIONS_scenarioA.md` + `SUBMISSION_CHECKLIST.md` — step-by-step npj Editorial Manager walkthrough + pre-flight check.

## Bundles

- `anonymous_code.zip` — code archive with author identifiers stripped.
- `reviewer_bundle.zip` — single-shot upload bundle (regenerate before submit).
