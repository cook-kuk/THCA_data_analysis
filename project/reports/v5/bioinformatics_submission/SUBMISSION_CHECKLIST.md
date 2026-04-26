# Bioinformatics Submission Checklist (v2 — all weaknesses addressed)

Target venue: *Bioinformatics* (Oxford University Press) — Original Paper.
Submission portal: <https://academic.oup.com/bioinformatics>

## Summary of fixes vs v1

| # | Weakness flagged | Fix |
|---|---|---|
| 1 | Word count short (~2,950) | Expanded to **5,597 words (texcount)**; body 5,398 |
| 2 | Fig 6 used synthetic PCA | Re-rendered from `data_processed/v5_cross_cancer/THCA/X_combined.npz` (real $n=392$ PCA, top-3000 variance, within-BRAF histogram) |
| 3 | Fig 7 used random noise for 45/50 pathways | Re-rendered from `v8_pathway_perpath_dial.tsv` (all 50 Hallmarks, real DIAL) |
| 4 | Fig 5 used constant error bar | Re-rendered from `v8_metasoft_forest_data.tsv` with real Hanley–McNeil 95% CIs |
| 5 | Fig 3 unlabelled as schematic | Explicitly labelled as schematic in panel titles and figure super-title |
| 6 | Fig 4, Fig 8 used hardcoded fallback data | Re-rendered from `v8_combatseq_vs_combat.tsv` and `v8_quantum_comparison.tsv` (real) |
| 7 | Lemma proof hand-waved sign-flip | Rewritten as **Lemma: first-moment annihilation** (provable); Corollary states residual predictive content; a "what the lemma does not prove" paragraph explicitly disclaims the probabilistic flip claim |
| 8 | `han2012mvalue` venue wrong | Kept venue as PLoS Genetics 2012 (it is correct); added clarifying note |
| 9 | TCGA/Agrawal duplicated | Removed `agrawal2014thca` entry |
| 10 | `pu2021thyroid` hallucinated | Removed; retargeted to `giordano2005thyroid` and `landa2016thyroid` |
| 11 | Abstract 206 words | Trimmed to **164 words** |
| 12 | SNUBH n=15 prospective claim unverified | Softened to "independent prospective clinical cohort whose acquisition is in progress" — no specific cohort asserted |
| 13 | LGG specificity weakness not flagged | New Limitations paragraph explicitly notes LGG is within-TCGA TSS split; conservative specificity denominator (15/15 cross-platform non-THCA) reported |
| 14 | Supplementary Table 1 duplicated main Table 1 | Replaced with subtype-scheme + exclusion paragraph |
| 15 | hyperref `[H]` duplicate-destination warnings | Replaced all `[H]` with `[htbp]`; removed `float` package; warnings gone |
| 16 | Authorship placeholder | Still placeholder — requires user to set real names/ORCIDs before submission |
| 17 | `bioinformatics.cls` not used | Still `article` class — swap `.cls` before submission (no content changes required) |
| 18 | Dashboard URL not verified | URL listed verbatim; user to verify reachability before submission |

## Pre-submission verification

- [x] **Abstract** 164 words (structured: Motivation / Results / Conclusion / Availability / Contact)
- [x] **Main text** 5,398 body words, 5,597 total — within 5,500–6,500 target (texcount -inc -sum)
- [x] **Main PDF** 18 pages; compiles clean under pdflatex + bibtex (no errors, no duplicate-label warnings)
- [x] **Figures** 8 (PDF vector + PNG 300 DPI); Fig 3 explicitly labelled schematic; all others driven by real TSVs or the harmonised NPZ matrix
- [x] **Table 1** (cohort availability) via `\input{figures/table1.tex}`, booktabs formatted
- [x] **Supplementary** 9 pages (S1–S9); compiles clean
- [x] **References** 58 BibTeX entries, covering batch correction (ComBat, SVA, ComBat-seq, Harmony, scVI, MNN, FastRNA), Han-lab statgen (METASOFT, m-value, BUHMBOX, FastRNA, Sul LMM), DE/statistics (DESeq2, pydeseq2, limma, edgeR, GSVA, ssGSEA, Hallmark, GSEA), quantum ML (Havlicek, Schuld, Huang, Biamonte), thyroid biology (TCGA 2014, Landa 2016, Giordano 2005, Xing, Nikiforov, Fagin, Kondo), subtype classification (PAM50, CMS, GBM, Pan-Cancer), ML auditing (Gebru, Mitchell, D'Amour), AUC/meta statistics (Hanley–McNeil, DerSimonian–Laird, Higgins I², VanderWeele E-value), tools (SciPy, NumPy, matplotlib, scikit-learn, XGBoost, statsmodels, UMAP)
- [x] **Cover letter** drafted with editor salutation, contributions, significance, recommended reviewers, COI, data availability
- [x] **Code repository** URL: dashboard at <http://40.82.129.113:8012/>
- [x] **Data availability statement** included in cover letter and Conclusion
- [x] **Ethics statement** — public data only; no primary patient consent required
- [x] **Conflict of interest** declared (none)
- [ ] **ORCID IDs** — to be filled by authors at submission time
- [ ] **Corresponding author** — to be confirmed (user email `kukshomr@gmail.com` listed)

## User action items

1. **Confirm authorship order.** Draft currently lists "THCA-Dash Project Contributors" as a placeholder.
2. **Finalize recommended reviewers.** Draft suggests Leek, Johnson, Han, Kendziorski, Sun.
3. **Verify dashboard URL reachability.** <http://40.82.129.113:8012/> must be up through review.
4. **Swap in OUP `bioinformatics.cls`.** No content change required; only `\documentclass` line.
5. **Optional bioRxiv pre-print** before submission.
6. **Decide public GitHub release timing** (pre-submission recommended for reviewer transparency).

## Files produced

```
reports/v5/bioinformatics_submission/
├── abstract.tex                               (164 words)
├── main.tex             →  main.pdf           (18 pages, 5,597 texcount words)
├── references.bib                             (58 entries)
├── cover_letter.md
├── SUBMISSION_CHECKLIST.md
├── figures/
│   ├── table1.tex                             (cohort availability, booktabs)
│   ├── fig1_heatmap.{pdf,png}                 (real v5p1_dial_all_cancers.tsv)
│   ├── fig2_scatter.{pdf,png}                 (real)
│   ├── fig3_flip_geom.{pdf,png}               (schematic, labelled as such)
│   ├── fig4_combatseq.{pdf,png}               (real v8_combatseq_vs_combat.tsv)
│   ├── fig5_forest.{pdf,png}                  (real forest TSV, Hanley-McNeil CI)
│   ├── fig6_mechanism.{pdf,png}               (real PCA + real KS + real FastRNA)
│   ├── fig7_pathway.{pdf,png}                 (real 50-pathway DIAL)
│   └── fig8_quantum.{pdf,png}                 (real v8_quantum_comparison.tsv)
└── supplementary/
    └── supplementary.tex → supplementary.pdf  (9 pages, S1–S9)
```

After user sign-off, submit via the OUP Online Editorial System at
<https://academic.oup.com/bioinformatics>.
