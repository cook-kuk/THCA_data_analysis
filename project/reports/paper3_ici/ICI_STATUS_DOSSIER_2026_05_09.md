# ICI status dossier — 2026-05-09

## Verdict

Paper 3 should stay framed as **ICI vulnerability / ICI-readiness / immunogenomic prioritization**, not thyroid ICI response prediction. The best current evidence is:

- Thyroid-side biology: ATC is inflamed-but-myeloid-suppressed and dedifferentiated (`myeloid_suppressive d=2.53`, `thyroid_differentiation d=-2.51`).
- RAI-side bridge: GSE151179 post-RAI vs pre-RAI shows thyroid-differentiation loss (`d=-1.01`, `p=3.7e-04`) with HLA-II / myeloid upward trends.
- External ICI coherence: four pre-treatment ICI cohorts, `n=421`, show IFNG, HLA-I, checkpoint, and cytolytic modules associated with response.
- Load-bearing response component: HLA-I expression module, not HLA-II, carries the strongest HLA class signal (`HLA-I OR=1.35, p=0.014`; `HLA-II OR=1.02, p=0.839`).
- BRAF-cPTC candidate layer: DM1 BRAF-cPTC has TIS-Ayers `d=1.35`, FDR `=1.3e-08`; 44/110 are ICI-likely by DM1 + TIS + CD8 + HLA-I.

## Claim boundary

- Allowed: ICI vulnerability, ICI-readiness, immunogenomic prioritization, response-biology coherence in non-thyroid ICI cohorts.
- Forbidden: thyroid ICI response predictor, treatment recommendation, clinical-grade selection, causal immunotherapy mechanism.

## Current decision

Track A remains frozen. True Track B requires Paper 1 bioRxiv, Paper 2 closure, and explicit `Paper 3 Track B 시작`.

## Source paths

- Track B-lite report: `project/reports/paper3_ici/paper3_ici_track_B_lite_analysis_report.md`
- Track A frozen bundle: `project/reports/paper3_ici/PAPER3_ICI_TRACK_A_BUNDLE.md`
- ATC meta TSV: `project/results/paper3_ici_track_b_lite/atc_vs_other_meta_pooled.tsv`
- RAI axis TSV: `project/results/paper3_ici_track_b_lite/gse151179_rai_after_vs_before.tsv`
- Phase C summary: `project/results/paper11_pancancer/phase_C_ICI/results/tables/final_phase_C_ICI_summary_table.tsv`
- Track18 HLA report: `project/results/hla_deepdive_2026_05_08/track18_ici_hla/track18_report.md`
- R5 report: `project/results/p2_braf_nature_sprint_2026_05_09/r5_ici_cohorts/R5_REPORT.md`
- R10 report: `project/results/p2_braf_nature_sprint_2026_05_09/r10_ici_regimen/R10_REPORT.md`
- H27 report: `project/results/p2_braf_nature_sprint_2026_05_09/h27_ici_panels/H27_REPORT.md`
