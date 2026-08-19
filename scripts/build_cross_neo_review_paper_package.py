#!/usr/bin/env python3
"""Create a review-first manuscript package for neoantigen prioritization.

This package is scaffolding plus evidence tables. It does not touch the
voice-protected THCA manuscript sections.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "project/results/cross_neo_review_2026_05_10"


REFERENCES = [
    {
        "short": "BigMHC",
        "year": 2023,
        "type": "predictor",
        "claim": "Pan-allelic deep models for MHC-I presentation and transfer-learned neoepitope immunogenicity.",
        "why_it_matters": "Strong public competitor and example of presentation-to-immunogenicity transfer learning.",
        "caveat": "Uses public datasets including TESLA, IEDB, NEPdb and prior predictor resources; clean benchmark overlap must be audited.",
        "url": "https://www.nature.com/articles/s42256-023-00694-6",
    },
    {
        "short": "PRIME",
        "year": 2021,
        "type": "predictor",
        "claim": "Combines antigen presentation and TCR-recognition propensity for CD8 epitope immunogenicity.",
        "why_it_matters": "Canonical immunogenicity predictor beyond binding affinity.",
        "caveat": "Binding and TCR-recognition components can still inherit public-corpus and allele-motif biases.",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7897774/",
    },
    {
        "short": "IMPROVE",
        "year": 2024,
        "type": "predictor_dataset",
        "claim": "Feature model from broad-scale validation of 17,520 neoepitope candidates with 467 recognized across 70 patients.",
        "why_it_matters": "One of the strongest current examples of clinically grounded, patient-level neoepitope screening data.",
        "caveat": "Includes patient and tumor-context features that may not be available in all deployment settings.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/38633261/",
    },
    {
        "short": "Muller_Immunity",
        "year": 2023,
        "type": "benchmark_dataset",
        "claim": "Harmonized WES/RNA-seq neoantigen datasets and ML ranking improved immunogenicity prediction.",
        "why_it_matters": "Shows source harmonization and feature engineering can materially change rankings.",
        "caveat": "Even harmonized datasets remain sparse relative to HLA/peptide/source space.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/37816353/",
    },
    {
        "short": "TESLA",
        "year": 2020,
        "type": "benchmark_consortium",
        "claim": "Consortium benchmark identified key parameters for tumor epitope immunogenicity.",
        "why_it_matters": "Central public benchmark used by many downstream predictors.",
        "caveat": "Because it is central, overlap-sensitive validation must be treated carefully.",
        "url": "https://www.cell.com/cell/fulltext/S0092-8674(20)31119-9",
    },
    {
        "short": "Rojas_PDACH_mRNA",
        "year": 2023,
        "type": "clinical_vaccine",
        "claim": "Personalized RNA neoantigen vaccine induced T cell responses in pancreatic cancer; responders had delayed recurrence.",
        "why_it_matters": "Clinical motivation for better neoantigen selection in low-response settings.",
        "caveat": "Small early-phase study; response correlation is not the same as broad clinical validation.",
        "url": "https://www.nature.com/articles/s41586-023-06063-y",
    },
    {
        "short": "PDAC_extended_followup",
        "year": 2025,
        "type": "clinical_vaccine",
        "claim": "Extended follow-up reported persistent neoantigen-specific CD8 T cells and continued RFS association.",
        "why_it_matters": "Supports durable immunologic relevance of correctly selected vaccine targets.",
        "caveat": "Still a specialized clinical setting and not a direct predictor benchmark.",
        "url": "https://www.nature.com/articles/s41586-024-08508-4",
    },
    {
        "short": "RCC_neoantigen_vaccine",
        "year": 2025,
        "type": "clinical_vaccine",
        "claim": "Phase I personalized neoantigen vaccine in renal cell carcinoma generated antitumor immunity.",
        "why_it_matters": "Shows feasibility even in lower mutation-burden tumors.",
        "caveat": "Nine-patient phase I study; not definitive efficacy evidence.",
        "url": "https://www.nature.com/articles/s41586-024-08507-5",
    },
    {
        "short": "mRNA_4157_KEYNOTE_942",
        "year": 2024,
        "type": "clinical_vaccine",
        "claim": "Individualized mRNA neoantigen therapy plus pembrolizumab improved melanoma disease-control signals in phase 2b.",
        "why_it_matters": "Commercial and clinical proof that ranking matters at product scale.",
        "caveat": "Vaccine design pipeline is proprietary and not a transparent predictor benchmark.",
        "url": "https://doi.org/10.1016/S0140-6736(23)02268-7",
    },
    {
        "short": "Performance_ceiling_commentary",
        "year": 2023,
        "type": "perspective",
        "claim": "Current dataset biases and limited high-quality data constrain neoantigen prediction generalizability.",
        "why_it_matters": "Frames the review's central problem: benchmark quality can dominate apparent method quality.",
        "caveat": "Perspective article; should be paired with quantitative benchmark evidence.",
        "url": "https://www.nature.com/articles/s43018-023-00675-z",
    },
    {
        "short": "Personalized_vaccine_review_2025",
        "year": 2025,
        "type": "review",
        "claim": "Recent review of personalized neoantigen vaccine platforms and implementation challenges.",
        "why_it_matters": "Useful clinical and translational framing.",
        "caveat": "Broad review; this manuscript should differentiate by focusing on prioritization and benchmarking.",
        "url": "https://www.frontiersin.org/journals/oncology-reviews/articles/10.3389/or.2025.1541326/full",
    },
]


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def write_reference_matrix() -> None:
    pd.DataFrame(REFERENCES).to_csv(OUT / "review_reference_evidence_matrix.tsv", sep="\t", index=False)


def write_strategy() -> None:
    text = """# Review-First Strategy

## Working Title
From neoantigen discovery to vaccine-ready prioritization: why benchmark contamination, source shift, and top-k utility now matter more than AUROC

## Core Thesis
Personalized neoantigen vaccines are clinically credible again, but the computational bottleneck has shifted. The field no longer only needs better MHC binding prediction. It needs contamination-controlled, source-aware, patient-context-aware ranking systems that optimize the handful of candidates actually manufactured or assayed.

## Why Review First
1. It lets us define the problem before presenting CROSS-Neo.
2. It makes public-overlap and benchmark leakage a field-wide issue, not a defensive caveat about our model.
3. It positions top-k precision, AUPRC, OOD behavior, and abstention as practical clinical/product metrics.
4. It creates the conceptual bridge for the original CROSS-Neo paper: a contamination-controlled pan-allele prioritizer rather than another binding predictor.

## Target Journal Shape
- Fast review / perspective: Nature Reviews Clinical Oncology, Nature Reviews Immunology, Cancer Discovery review, Trends in Cancer, JITC review, Frontiers if speed matters.
- If business timing dominates: preprint + white paper first, then invited-style review.

## Do Not Claim
- Do not claim external validation for CROSS-Neo.
- Do not claim quantum advantage.
- Do not claim public predictors are invalid; say they are powerful but benchmark-overlap-sensitive.
- Do not claim neoantigen vaccines already have universal clinical efficacy.
"""
    (OUT / "REVIEW_FIRST_STRATEGY.md").write_text(text)


def write_manuscript_skeleton() -> None:
    text = """# From Neoantigen Discovery to Vaccine-Ready Prioritization

## Abstract Draft
Personalized neoantigen vaccines have re-emerged as a clinically plausible immunotherapy strategy, supported by durable T cell responses in pancreatic cancer, renal cell carcinoma, and melanoma studies. Yet the translational bottleneck is no longer simply whether candidate peptides bind HLA. In clinical and product settings, only a small number of candidates can be manufactured, assayed, or administered, making prioritization a top-k ranking problem under severe class imbalance, incomplete negative labels, source shift, and public benchmark overlap. This review synthesizes recent progress in neoantigen immunogenicity prediction, including binding-centric predictors, transfer-learned immunogenicity models, feature-based patient-context models, and structure/TCR-aware approaches. We argue that the next generation of neoantigen prioritizers should be evaluated by contamination-controlled splits, AUPRC, top-k precision, enrichment over prevalence, calibration, and OOD/abstention behavior rather than AUROC alone. We conclude with a practical evaluation contract for vaccine-ready computational prioritization.

## 1. Clinical Re-Emergence of Personalized Neoantigen Vaccines
- Melanoma: individualized mRNA therapy plus pembrolizumab.
- Pancreatic cancer: autogene cevumeran and durable CD8 T cell responses.
- Renal cell carcinoma: feasibility and immunogenicity in a lower mutation-burden tumor.
- Practical implication: better candidate ranking matters because manufacturing slots are scarce.

## 2. Why Binding Prediction Was Necessary but Is Not Sufficient
- MHC binding and presentation are required filters.
- Immunogenicity also depends on mutant-WT contrast, TCR-facing features, clonality, expression, tumor context, immune state, and tolerance.
- Binding-centric success can inflate apparent utility when evaluated on overlap-sensitive public corpora.

## 3. Current Predictor Families
- Binding/presentation predictors: NetMHCpan, MHCflurry, MixMHCpred, BigMHC-EL.
- Immunogenicity predictors: PRIME, DeepImmuno, BigMHC-IM, TransPHLA, T-SCAPE, IMPROVE.
- Dataset-harmonization and feature models: Muller et al., IMPROVE.
- Structure/TCR-aware and emerging PLM approaches.

## 4. Benchmarking Failure Modes
- Exact peptide-HLA overlap.
- Near-peptide and source-protein overlap.
- HLA allele and supertype shortcuts.
- Source/study/assay shift.
- Positive-label source bias and ambiguous negatives.
- AUROC-only reporting under low prevalence.

## 5. Metrics That Match Vaccine Decisions
- AUPRC over AUROC for rare positives.
- top5/top10 precision and recall@k.
- enrichment over prevalence.
- Brier/ECE calibration.
- abstention coverage versus precision.
- per-HLA, per-study, source-heldout, and time-heldout reporting.

## 6. Toward Vaccine-Ready Prioritization
- Mutant-WT counterfactual modeling.
- Retrieval evidence with explicit contamination flags.
- Structure/geometry and TCR-facing features.
- PU-aware treatment of negatives.
- Pan-allele generalization with HLA as context, not a memorized boundary.
- OOD-aware abstention and wet-lab triage.

## 7. Practical Evaluation Contract
- Train-only retrieval indices.
- No public predictor scores as clean training features if used as comparators.
- Exact/near/source/HLA/study/time split families.
- Separate clean manuscript claims from product-assisted triage claims.

## 8. Bridge to Original Work
- CROSS-Neo should be introduced later as an implementation of this review's contract.
- The review creates the rationale: not another MHC-binding predictor, but a contamination-controlled prioritization stack.

## Figure Plan
See `review_figure_plan.tsv`.

## Table Plan
See `review_table_plan.tsv`.
"""
    (OUT / "review_manuscript_skeleton.md").write_text(text)


def write_table_and_figure_plans() -> None:
    figures = [
        {
            "figure": "Fig 1",
            "title": "The neoantigen vaccine prioritization funnel",
            "message": "Candidate space collapses from thousands of variants to a handful of manufacturable targets.",
            "data_needed": "schematic; no new data",
        },
        {
            "figure": "Fig 2",
            "title": "Predictor families and what biological step they model",
            "message": "Binding, presentation, immunogenicity, TCR recognition, tumor context and product triage are distinct layers.",
            "data_needed": "reference matrix",
        },
        {
            "figure": "Fig 3",
            "title": "Benchmark contamination and source-shift failure modes",
            "message": "Exact overlap, near overlap, HLA shortcut, and study shift can inflate apparent performance.",
            "data_needed": "CROSS-Neo overlap audit plus public examples",
        },
        {
            "figure": "Fig 4",
            "title": "Metrics aligned to vaccine decisions",
            "message": "Top-k precision, AUPRC, enrichment and abstention are more deployment-relevant than AUROC alone.",
            "data_needed": "strong competitor gauntlet and simulated prevalence curves",
        },
        {
            "figure": "Fig 5",
            "title": "Evaluation contract for vaccine-ready rankers",
            "message": "A checklist for split safety, public-comparator separation, OOD and calibration.",
            "data_needed": "contract schematic",
        },
    ]
    tables = [
        {
            "table": "Table 1",
            "title": "Clinical neoantigen vaccine signals and unresolved prioritization needs",
            "rows": "melanoma, PDAC, RCC, other selected trials",
        },
        {
            "table": "Table 2",
            "title": "Computational predictor families",
            "rows": "NetMHCpan, MHCflurry, MixMHCpred, PRIME, BigMHC, DeepImmuno, IMPROVE, TransPHLA, T-SCAPE",
        },
        {
            "table": "Table 3",
            "title": "Benchmark failure modes and required controls",
            "rows": "exact overlap, near overlap, HLA shortcut, source shift, label ambiguity, public predictor leakage",
        },
        {
            "table": "Table 4",
            "title": "Recommended metrics for neoantigen rankers",
            "rows": "AUPRC, top-k precision, enrichment, Brier, ECE, abstention coverage, per-HLA and per-study metrics",
        },
    ]
    pd.DataFrame(figures).to_csv(OUT / "review_figure_plan.tsv", sep="\t", index=False)
    pd.DataFrame(tables).to_csv(OUT / "review_table_plan.tsv", sep="\t", index=False)


def write_original_bridge() -> None:
    text = """# Bridge From Review Paper to CROSS-Neo Original Paper

## Review Paper Role
Define the field problem:
- public predictors are powerful but overlap-sensitive;
- vaccine deployment is a top-k ranking problem;
- negative labels are often incomplete;
- source/HLA/study shifts can dominate performance;
- evaluation must be contamination-controlled and OOD-aware.

## Original CROSS-Neo Paper Role
Demonstrate one implementation:
- pan-allele mutant-WT counterfactual prioritizer;
- hard-decoy and PU-aware pressure for ambiguous negatives;
- explicit public-comparator separation;
- retrieval evidence marked as safe/unsafe;
- source-aware and OOD-aware triage;
- business/product mode separated from clean scientific claims.

## Paper Pairing
1. Review title: From neoantigen discovery to vaccine-ready prioritization.
2. Original title: CROSS-Neo: contamination-controlled pan-allele neoantigen prioritization for vaccine candidate triage.

## Clean Claim Boundary for Original Paper
- Internal locked/split-safe validation only unless an independent prospective set is obtained.
- No quantum advantage claim.
- Public predictor scores are comparators or product-assist signals, not clean model features.
- Main outcome should be AUPRC/top-k/enrichment and failure-mode audit.
"""
    (OUT / "CROSS_NEO_ORIGINAL_PAPER_BRIDGE.md").write_text(text)


def write_review_dossier() -> None:
    text = """# Review Paper Build Dossier

## One-Sentence Thesis
Neoantigen vaccine success now depends less on finding more candidate binders and more on ranking the few candidates worth manufacturing under contamination, source shift and incomplete-label uncertainty.

## Strongest Angle
This review is not a generic vaccine-platform review. It is a computational-translational review focused on vaccine-ready prioritization: how to evaluate, benchmark and deploy neoantigen rankers.

## What Makes It Timely
- Recent clinical vaccine signals in melanoma, pancreatic cancer and renal cell carcinoma.
- New stronger predictors including BigMHC and IMPROVE.
- Recognition that dataset overlap, low prevalence and source shift can create misleading benchmark wins.
- Product pressure: the clinically relevant output is a short candidate list, not a calibrated probability for every peptide.

## Immediate Writing Order
1. Lock title, abstract, and figure plan.
2. Fill Table 2 predictor taxonomy from the reference evidence matrix.
3. Draft Section 4 benchmarking failure modes using CROSS-Neo audits as internal examples, without overclaiming.
4. Draft Section 5 metrics contract.
5. Only then introduce CROSS-Neo as the follow-on original work.
"""
    (OUT / "REVIEW_BUILD_DOSSIER.md").write_text(text)


def main() -> None:
    ensure_dirs()
    write_reference_matrix()
    write_strategy()
    write_manuscript_skeleton()
    write_table_and_figure_plans()
    write_original_bridge()
    write_review_dossier()
    summary = {
        "out_dir": str(OUT),
        "review_first": True,
        "original_follow_on": "CROSS-Neo contamination-controlled pan-allele prioritizer",
        "n_reference_seeds": len(REFERENCES),
    }
    (OUT / "review_package_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[review-package] wrote {OUT}")


if __name__ == "__main__":
    main()
