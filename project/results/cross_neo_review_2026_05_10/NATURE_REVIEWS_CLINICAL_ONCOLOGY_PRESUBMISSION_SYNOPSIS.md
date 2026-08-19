# Pre-Submission Synopsis: Nature Reviews Clinical Oncology

## Proposed Article Type
Review or Perspective

## Proposed Title
From neoantigen discovery to vaccine-ready prioritization

## Alternative Titles
1. Vaccine-ready neoantigen prioritization in the era of personalized cancer vaccines
2. The next bottleneck in personalized neoantigen vaccines: ranking manufacturable targets
3. From binding prediction to vaccine candidate selection in personalized cancer immunotherapy

## One-Sentence Pitch
Personalized neoantigen vaccines are becoming clinically credible, but the computational bottleneck has shifted from finding HLA binders to ranking a small number of manufacturable vaccine candidates under class imbalance, incomplete negative labels, source shift and benchmark contamination.

## Rationale
Recent clinical studies in melanoma, pancreatic cancer and renal cell carcinoma have renewed confidence that personalized neoantigen vaccines can induce durable antitumor T cell responses. At the same time, computational methods for neoantigen prediction have rapidly expanded from binding and presentation predictors to immunogenicity, TCR-propensity, structure-aware, protein language model and patient-context models. However, the clinical utility of these tools is increasingly limited by evaluation problems: public benchmark overlap, source/study shift, HLA shortcuts, incomplete negative labels and AUROC-centric reporting despite the practical need to select only a top-ranked handful of vaccine candidates.

This Review would synthesize the clinical and computational landscape around a central translational question: how should neoantigen candidates be prioritized for vaccine manufacture and immune testing?

## Proposed Scope
- Clinical motivation from recent personalized neoantigen vaccine studies.
- Why MHC binding and presentation prediction are necessary but insufficient.
- Predictor families: NetMHCpan, MHCflurry, MixMHCpred, PRIME, BigMHC, DeepImmuno, TransPHLA, T-SCAPE, IMPROVE and emerging PLM/structure/TCR-aware models.
- Benchmarking failure modes: exact peptide-HLA overlap, near-peptide overlap, source-protein overlap, source/study shift, HLA allele shortcuts and public-corpus reuse.
- Practical metrics: AUPRC, top-k precision, enrichment over prevalence, recall@k, calibration, OOD detection and abstention.
- A proposed evaluation contract for vaccine-ready neoantigen prioritizers.

## Why Now
The field has reached a point where clinical vaccine pipelines are credible and computational predictors are abundant, but the connection between predictor performance and manufacturable vaccine decisions remains underdeveloped. A timely Review can set evaluation standards before the next wave of personalized vaccine trials and AI-based prioritization systems.

## Distinctiveness
This would not be a generic review of neoantigen vaccine platforms or a catalogue of algorithms. It would focus on the translational decision layer: how candidate ranking should be benchmarked and reported when only a small number of peptides can be manufactured, assayed or administered.

## Proposed Display Items
1. Neoantigen vaccine prioritization funnel from variants to manufacturable peptides.
2. Predictor taxonomy mapped to biological layers: binding, presentation, TCR recognition, tumor context and product triage.
3. Benchmark contamination and source-shift failure modes.
4. Metrics for vaccine-ready ranking: AUPRC, top-k precision, enrichment and abstention.
5. Proposed evaluation contract for clinical/product-grade neoantigen prioritization.

## Proposed Boxes
- Box 1: Why AUROC can mislead in low-prevalence neoantigen ranking.
- Box 2: Public benchmark overlap and what a clean split should exclude.
- Box 3: Practical reporting checklist for neoantigen prioritization studies.

## Author Positioning
The manuscript will be written as a balanced field review. Any ongoing CROSS-Neo work would be treated only as motivation for future benchmarking standards, not as unpublished evidence or a dominant focus.

## Expected Length
4,500-6,000 words, 4-6 display items, 100-150 references.

## Editorial Fit
Nature Reviews Clinical Oncology covers immunotherapy, experimental therapies, medical oncology, statistics, genetics and translational cancer research. The proposed article sits at the interface of personalized vaccine clinical translation and computational decision-making.
