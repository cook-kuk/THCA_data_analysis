# PAAD + THCA cancer-vaccine combo and neoantigen-patient triage - 2026-05-09

## Bottom line

This is a decision dossier scaffold, not clinical advice.

The priority split is clear:

- **PAAD**: vaccine-first is defensible only in the **resected/MRD low-burden** setting. The strongest current anchors are personalized mRNA neoantigen vaccine + atezolizumab + mFOLFIRINOX, and off-the-shelf lymph-node-targeted mKRAS vaccine for MRD-positive KRAS disease.
- **THCA**: vaccine-first is not defensible for routine PTC. The usable strategy is **patient-selection-first**: ATC/PDTC/progressive RAI-refractory DTC, intact antigen presentation, actionable backbone therapy, and immune-active HT/TLS/IFN/HLA-II context when possible.
- **Neoantigen patient selection** must require WES/RNA/HLA typing plus HLA LOH/B2M/presentation checks. The local neoantigen model is useful for ranking and abstention, not for definitive immunogenicity proof.

## Files

- `evidence_matrix.tsv` - external and local evidence anchors.
- `patient_triage_matrix.tsv` - practical GO / COMBO-GO / ABSTAIN logic by disease.
- `combination_priority_matrix.tsv` - prioritized combination strategies.
- Hub page: `project/papers_hub_2026_05_04/paad_thca_vaccine_combo_triage_2026_05_09.html`.
- Figure assets mirrored into `project/papers_hub_2026_05_04/manuscript_artifacts/` so local `:8012` pages can render the same PAAD/THCA figures as the `/var/www` deployment.

## 2026-05-09 follow-up fix

The first deployed version referenced `manuscript_artifacts/*.png` files that existed under `/var/www/papers/...` but not under the repo hub directory served by `:8012`, so `03 Disease Logic` showed broken images. Fixed by mirroring the 78 MB `manuscript_artifacts` bundle into `project/papers_hub_2026_05_04/manuscript_artifacts/` and expanding `03 Disease Logic` to a 6-panel figure block:

- PAAD KRAS distribution.
- PAAD top vaccine candidates.
- PAAD HLA coverage.
- THCA TMB landscape.
- THCA Korean HLA gap.
- THCA paper pipeline.

## Highest-value combinations

| Rank | Segment | Combination | Read |
|---:|---|---|---|
| 1 | PAAD resected/personalized | autogene cevumeran + atezolizumab + mFOLFIRINOX | Best clinical anchor; responder RFS signal at 3.2-year follow-up. |
| 2 | PAAD MRD+ KRAS-public | ELI-002 7P/2P | Best scalable shared-neoantigen strategy; MRD-selected. |
| 3 | PAAD resected KRAS broad | mKRAS-VAX + nivolumab + ipilimumab | Strong immunogenicity, but checkpoint toxicity and small-n DFS caution. |
| 4 | THCA BRAF V600E ATC | dabrafenib + trametinib +/- pembrolizumab | Backbone first; vaccine is only a trial add-on concept. |
| 5 | THCA RAI-refractory DTC/PDTC | lenvatinib + pembrolizumab +/- vaccine if antigen-positive | Stronger near-term combo backbone than vaccine alone. |
| 6 | THCA HT/TLS immune-active high-risk PTC | vaccine + PD-1 only in high-risk/progressive trial setting | Biology is plausible; thyroid-specific vaccine evidence is missing. |

## Source anchors

- PDAC personalized mRNA vaccine long follow-up: https://www.nature.com/articles/s41586-024-08508-4
- ELI-002 2P final AMPLIFY-201: https://www.nature.com/articles/s41591-025-03876-4
- mKRAS-VAX + nivolumab/ipilimumab: https://www.nature.com/articles/s41467-026-68324-4
- Autogene cevumeran phase 2 registry: https://clinicaltrials.gov/study/NCT05968326
- ELI-002 7P registry: https://clinicaltrials.gov/study/NCT05726864
- FDA dabrafenib/trametinib ATC approval: https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-dabrafenib-plus-trametinib-anaplastic-thyroid-cancer-braf-v600e-mutation
- Lenvatinib + pembrolizumab RR-DTC phase 2: https://pmc.ncbi.nlm.nih.gov/articles/PMC11883846/
- NCI NEO-COMBAT XL trial: https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2025-01421&r=1
- Local method blueprint: `project/results/p_neo_bayesian_2026_05_09/NEOANTIGEN_METHOD_BLUEPRINT_2026_05_09.md`

## Caveats to preserve

- PDAC vaccine data are promising but early-phase and enriched for low-burden/resected/MRD windows.
- THCA has very low TMB; thyroid vaccine translation should not be overclaimed from PAAD.
- HT/TLS thyroid tumors may be immune-active, but checkpoint-autoimmunity risk is not a side issue.
- HLA presentation failure is an absolute vaccine blocker.
- Public neoantigen predictors remain leakage-sensitive; report uncertainty and abstain on OOD/high-risk cases.
