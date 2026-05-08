# Next analysis plan — proteogenomic_v1

*Marathon-mode rules apply: any of these that count as new sprint analysis must wait for explicit user greenlight.*

## 1. Recurrence-risk metabolic subtype replication (MTBLS3339)
- Use harmonized sample index → metabolite long matrix.
- Test thyroid-hormone pathway metabolites (T3/T4/MIT/DIT/iodotyrosine + tyrosine, iodide proxies if measured) for recurrence-risk-stratified diff.
- Cross-reference with Wang 2024 RNA panel and our 8-gene panel: do the metabolites correlate with the same dediff axis?
- Output: `metabolomics_recurrence_risk_diff.tsv` + a `metabolite_x_8gene_corr.tsv`.

## 2. Advanced DTC CC1/CC2/CC3 proteomic characterization (CellRepMed2026)
- Map our DIAL-lite 7-module gene list onto the protein matrix.
- Compute module mean per sample, group by CC1/CC2/CC3.
- Compare CC3 (immunogenic) module profile to Mun 2025 ATC profile from `paper3_mun2025_dediff_layer/module_ATC_vs_PTC.tsv`. If CC3 ≈ ATC, the advanced-DTC immunogenic phenotype is on the same dediff axis we already characterized — strong replication evidence for Paper 3 reviewer-reserve.

## 3. Immune / metabolic / lineage axis comparison
- Combine #1 + #2: metabolomic phenotype (MTBLS) ↔ proteomic phenotype (CRM).
- Boundary: same direction at the *cohort-level summary*, NOT sample-level merge.
- This is the cleanest reviewer-reserve framing: 'three independent cohorts (Wang RNA+protein, Mun protein+phospho, CRM advanced DTC protein, MTBLS metabolomics) point to the same dediff axis.'

## 4. External public-data validation hooks
- TCGA-THCA HM450 methylation (n=503): already on disk, integrated into Paper 1 reviewer-reserve memo.
- TCGA-THCA miRNA-seq from GDC: deferred (sprint-class).
- Lu 2023 ATC scRNA per-cell-type 7-module: deferred (sprint-class).

## Boundary reminders
- HLA validation stays on a separate track (`v17_arcasHLA_korean_k2`).
- Track B (Paper 3) remains FROZEN until Paper 1 bioRxiv + Paper 2 A/B/C + explicit unfreeze (`v19_paper3_ici_track_a`).
- Voice-protected manuscript sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) are author-keyboard only.
