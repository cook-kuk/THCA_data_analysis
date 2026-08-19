#!/usr/bin/env python3
"""Build a higher-impact framing package for the CROSS-Neo review.

The goal is to elevate the manuscript from "review of predictors" to a
field-defining Perspective/Grand Challenge on vaccine-ready AI standards.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "project/results/cross_neo_review_2026_05_10"
FIG = OUT / "figures"


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def add_box(ax, x, y, w, h, text, color, fs=9.5):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        facecolor=color,
        edgecolor="none",
        alpha=0.96,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color="white", fontsize=fs, weight="bold")


def add_arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, linewidth=2.0, color="#5f6770"))


def make_flagship_figure() -> Path:
    path = FIG / "fig0_high_impact_grand_challenge_overview.png"
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.text(
        0.04,
        0.94,
        "Grand challenge: vaccine-ready AI for personalized cancer immunotherapy",
        fontsize=17,
        weight="bold",
        color="#16202a",
    )
    ax.text(
        0.04,
        0.89,
        "The missing standard is not another predictor; it is trustworthy top-k prioritization for manufacturable neoantigen vaccines.",
        fontsize=11,
        color="#4b5563",
    )

    # Three-column structure
    add_box(ax, 0.05, 0.66, 0.24, 0.12, "Clinical reality\nvaccines are credible again", "#2f6f9f", 10)
    ax.text(
        0.17,
        0.56,
        "melanoma\npancreatic cancer\nrenal cell carcinoma\nmanufacturing slots",
        ha="center",
        va="center",
        fontsize=9.3,
        color="#16202a",
    )

    add_box(ax, 0.38, 0.66, 0.24, 0.12, "Evaluation gap\nleaderboards can mislead", "#b36b2c", 10)
    ax.text(
        0.50,
        0.55,
        "public overlap\nsource shift\nHLA shortcuts\nambiguous negatives\nAUROC-only reporting",
        ha="center",
        va="center",
        fontsize=9.3,
        color="#16202a",
    )

    add_box(ax, 0.71, 0.66, 0.24, 0.12, "Field standard\nvaccine-ready ranking", "#4b8f73", 10)
    ax.text(
        0.83,
        0.55,
        "AUPRC\ntop-k precision\nenrichment\nOOD abstention\nprovenance",
        ha="center",
        va="center",
        fontsize=9.3,
        color="#16202a",
    )

    add_arrow(ax, (0.30, 0.72), (0.37, 0.72))
    add_arrow(ax, (0.63, 0.72), (0.70, 0.72))

    add_box(ax, 0.10, 0.25, 0.20, 0.13, "Review/Perspective\nsets the rules", "#5b5b8a", 10)
    add_box(ax, 0.40, 0.25, 0.20, 0.13, "CROSS-Neo original\nimplements them", "#326b7c", 10)
    add_box(ax, 0.70, 0.25, 0.20, 0.13, "Product demo\ntriages candidates", "#7d6f2f", 10)
    add_arrow(ax, (0.31, 0.315), (0.39, 0.315))
    add_arrow(ax, (0.61, 0.315), (0.69, 0.315))

    ax.text(
        0.08,
        0.12,
        "High-impact claim: personalized cancer vaccines need a reporting and benchmarking standard for the AI layer that selects what gets manufactured.",
        fontsize=11.5,
        weight="bold",
        color="#16202a",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def write_high_impact_strategy(fig0: Path) -> None:
    text = f"""# High-Impact Repositioning

## New Category

Do not pitch this as a neoantigen prediction review.

Pitch it as:

> A Grand Challenge for vaccine-ready AI in personalized cancer immunotherapy.

## High-Impact Thesis

Personalized neoantigen vaccines are now clinically plausible, but the field lacks a trustworthy standard for the AI layer that chooses which candidates are manufactured, assayed or administered. The next bottleneck is not peptide discovery; it is contamination-controlled, source-aware, top-k prioritization under incomplete-label uncertainty.

## Best Title

**Vaccine-ready AI for personalized cancer immunotherapy**

## Strong Alternative Titles

1. **The missing benchmark standard for personalized neoantigen vaccines**
2. **From neoantigen prediction to vaccine-ready AI**
3. **Personalized cancer vaccines need top-k AI, not AUROC leaderboards**
4. **A reporting standard for computational neoantigen prioritization**

## High-Impact Target Route

1. **Nature Medicine Perspective or Comment**
   - Best if framed as clinical AI standard for personalized vaccine decision-making.
   - Use a sharp, accessible, policy/standards angle.

2. **Nature Cancer Perspective**
   - Best fit for an opinionated cancer-field benchmarking standard.
   - Strong title: "Neoantigen prediction needs vaccine-ready benchmarks, not another AUROC leaderboard."

3. **Nature Reviews Immunology / Nature Reviews Clinical Oncology**
   - Best for authoritative review.
   - Must be broad and not dominated by CROSS-Neo.

4. **Nature Biomedical Engineering Perspective**
   - Best if framed as engineering/validation/deployment of computational systems for therapy selection.

5. **Cancer Cell / Cancer Discovery**
   - Strong if the story is made provocative and figure-led.

## What Makes It Big

- It connects clinical vaccine momentum with AI evaluation failure.
- It proposes a field standard rather than another model.
- It gives reviewers a vocabulary for public overlap, source shift, PU labels and top-k ranking.
- It creates space for CROSS-Neo as the implementation paper.

## Flagship Figure

`{fig0}`

## Claim Boundary

- Review/Perspective: standards and field synthesis.
- CROSS-Neo original: internal locked implementation.
- Product demo: retrospective business triage.
- No external validation claim.
- No quantum advantage claim.
"""
    (OUT / "HIGH_IMPACT_REPOSITIONING.md").write_text(text)


def write_target_matrix() -> None:
    rows = [
        {
            "rank": 1,
            "journal": "Nature Medicine",
            "format": "Perspective or Comment",
            "impact_logic": "clinical AI standard for personalized vaccine decisions",
            "recommended_title": "Vaccine-ready AI for personalized cancer immunotherapy",
            "why_it_could_hit": "timely clinical vaccine momentum plus AI validation/decision-standard problem",
            "risk": "must be broad medical interest and not too algorithmic",
            "action": "send presubmission abstract/synopsis",
        },
        {
            "rank": 2,
            "journal": "Nature Cancer",
            "format": "Perspective",
            "impact_logic": "cancer-field benchmark reset",
            "recommended_title": "Neoantigen prediction needs vaccine-ready benchmarks, not another AUROC leaderboard",
            "why_it_could_hit": "clear viewpoint, timely, provocative but balanced",
            "risk": "must avoid self-promotion and unpublished-data dependence",
            "action": "send sharper synopsis with Fig 0 and checklist",
        },
        {
            "rank": 3,
            "journal": "Nature Reviews Immunology",
            "format": "Review or Perspective proposal",
            "impact_logic": "immunology field standard for cancer vaccine antigen selection",
            "recommended_title": "From neoantigen prediction to vaccine-ready immunogenicity",
            "why_it_could_hit": "updates classic personalized cancer vaccine review with AI-era prioritization",
            "risk": "commissioned-review barrier; needs broad immunology beyond algorithms",
            "action": "pitch only after synopsis is clinically broadened",
        },
        {
            "rank": 4,
            "journal": "Nature Reviews Clinical Oncology",
            "format": "Review or Perspective proposal",
            "impact_logic": "translational oncology review",
            "recommended_title": "From neoantigen discovery to vaccine-ready prioritization",
            "why_it_could_hit": "clinical-trial relevance and practical candidate-selection bottleneck",
            "risk": "commissioned-review barrier",
            "action": "send 1-page clinical synopsis",
        },
        {
            "rank": 5,
            "journal": "Nature Biomedical Engineering",
            "format": "Perspective",
            "impact_logic": "validation/deployment standard for computational therapy-design systems",
            "recommended_title": "Benchmarking vaccine-ready AI systems for neoantigen prioritization",
            "why_it_could_hit": "computational systems for therapy optimization fit the journal scope",
            "risk": "may require stronger engineering/deployment angle",
            "action": "prepare 1000-word presubmission synopsis",
        },
        {
            "rank": 6,
            "journal": "Cancer Cell",
            "format": "Review/Perspective pitch",
            "impact_logic": "big cancer-immunology framework",
            "recommended_title": "The new bottleneck in personalized cancer vaccines",
            "why_it_could_hit": "clear cancer-field story with clinical and computational stakes",
            "risk": "hard invitation barrier and need senior field coauthor",
            "action": "use after senior immunology/oncology coauthor alignment",
        },
        {
            "rank": 7,
            "journal": "Cancer Discovery",
            "format": "Perspective/Review",
            "impact_logic": "oncology translational standard",
            "recommended_title": "Toward vaccine-ready neoantigen prioritization",
            "why_it_could_hit": "strong practical relevance to trial design",
            "risk": "less naturally computational unless clinical angle dominates",
            "action": "secondary prestige route",
        },
        {
            "rank": 8,
            "journal": "JITC",
            "format": "Review",
            "impact_logic": "relevant and realistic immunotherapy audience",
            "recommended_title": "A practical evaluation contract for neoantigen vaccine prioritization",
            "why_it_could_hit": "high fit and faster path",
            "risk": "lower category-setting prestige",
            "action": "backup direct route",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "high_impact_target_matrix.tsv", sep="\t", index=False)


def write_neoprior_checklist() -> None:
    rows = [
        {
            "item": "NEO-PRIOR-01",
            "requirement": "State whether the study is vaccine triage, immunogenicity prediction, binding/presentation prediction, or benchmark comparison.",
            "why": "Prevents treating all peptide-HLA scores as the same biological endpoint.",
        },
        {
            "item": "NEO-PRIOR-02",
            "requirement": "Report exact peptide-HLA, exact peptide, near-peptide and source-protein overlap against train and public corpora.",
            "why": "Separates true generalization from retrieval/memorization-sensitive gains.",
        },
        {
            "item": "NEO-PRIOR-03",
            "requirement": "Use train-only retrieval, preprocessing, scaling, calibration, feature selection and threshold selection.",
            "why": "Prevents test-fold leakage.",
        },
        {
            "item": "NEO-PRIOR-04",
            "requirement": "Include exact, near, source/study, HLA allele, HLA supertype and time/assay holdouts when feasible.",
            "why": "Random splits are insufficient for personalized vaccine deployment.",
        },
        {
            "item": "NEO-PRIOR-05",
            "requirement": "Report AUPRC, top5/top10 precision, enrichment over prevalence and recall@k as primary metrics.",
            "why": "Vaccine manufacture is a top-k decision under low prevalence.",
        },
        {
            "item": "NEO-PRIOR-06",
            "requirement": "Report AUROC only as a secondary metric.",
            "why": "AUROC can look acceptable while top-k utility collapses.",
        },
        {
            "item": "NEO-PRIOR-07",
            "requirement": "Treat negative labels as potentially ambiguous unless assay design proves true non-immunogenicity.",
            "why": "Many neoantigen negatives are untested, assay-limited or context-dependent.",
        },
        {
            "item": "NEO-PRIOR-08",
            "requirement": "Report calibration, Brier/ECE and OOD/abstention coverage.",
            "why": "Deployment needs confidence-aware candidate triage.",
        },
        {
            "item": "NEO-PRIOR-09",
            "requirement": "Separate public predictor scores used as comparators from model training features.",
            "why": "Clean claims require avoiding benchmark-comparator leakage.",
        },
        {
            "item": "NEO-PRIOR-10",
            "requirement": "Publish a data provenance table mapping each dataset to its role and claim boundary.",
            "why": "Makes clinical, benchmark and product claims auditable.",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "NEO_PRIOR_reporting_checklist.tsv", sep="\t", index=False)


def write_high_impact_synopsis() -> None:
    text = """# Nature Medicine / Nature Cancer High-Impact Synopsis

## Proposed Title
Vaccine-ready AI for personalized cancer immunotherapy

## Standfirst
Personalized neoantigen vaccines have re-entered clinical oncology, but the AI layer that selects candidate antigens is still evaluated like a prediction leaderboard rather than a therapeutic decision system.

## 150-Word Abstract
Personalized neoantigen vaccines are becoming clinically credible across melanoma, pancreatic cancer and renal cell carcinoma, renewing the need to select patient-specific vaccine targets. Yet the computational layer that prioritizes neoantigens is often evaluated with retrospective AUROC leaderboards, random splits and public benchmark corpora that may not reflect the clinical decision: choosing a small number of manufacturable candidates for a new patient. We argue that neoantigen prediction should be reframed as vaccine-ready AI: contamination-controlled, source-aware top-k ranking under incomplete-label uncertainty. This Perspective proposes NEO-PRIOR, a practical reporting and evaluation checklist for neoantigen prioritization, covering public-corpus overlap, exact and near holdouts, HLA and study shift, AUPRC, top-k precision, enrichment, calibration and OOD-aware abstention. Establishing such standards is essential before AI-driven vaccine candidate selection can be interpreted as clinically meaningful rather than benchmark-specific.

## Killer Sentence
The next advance in personalized cancer vaccines may come less from predicting more HLA binders than from trusting the shortlist that decides what gets manufactured.

## Core Claims
1. Personalized vaccine trials have made neoantigen selection clinically consequential.
2. The clinically relevant output is a short ranked list, not a genome-wide AUROC.
3. Benchmark contamination and source shift can make weak deployment models look strong.
4. Public predictors are valuable, but their scores must be separated from clean model claims.
5. NEO-PRIOR can become the minimum reporting standard for vaccine-ready neoantigen AI.

## Display Items
1. Grand challenge overview: clinical reality -> evaluation gap -> field standard.
2. NEO-PRIOR reporting checklist.
3. Benchmark failure modes and split controls.
4. Vaccine-ready metric panel: AUPRC, top-k, enrichment, calibration, OOD.

## Why It Is High Impact
This is not a methods review. It is a standards paper for a clinically emerging AI decision layer.
"""
    (OUT / "HIGH_IMPACT_NATURE_MEDICINE_NATURE_CANCER_SYNOPSIS.md").write_text(text)


def write_email() -> None:
    text = """# High-Impact Editor Email

## Nature Medicine

Subject: Presubmission enquiry: Perspective on vaccine-ready AI for personalized cancer immunotherapy

Dear Editors,

I am writing to ask whether Nature Medicine would consider a Perspective provisionally titled "Vaccine-ready AI for personalized cancer immunotherapy."

Personalized neoantigen vaccines have re-entered clinical oncology, with encouraging signals in melanoma, pancreatic cancer and renal cell carcinoma. However, the computational layer that selects patient-specific vaccine targets is still often evaluated like a retrospective prediction leaderboard rather than a therapeutic decision system. In practice, the relevant question is whether an AI system can prioritize the small number of neoantigens worth manufacturing or testing for a new patient without relying on public-corpus overlap, HLA shortcuts or source-specific label bias.

We propose a Perspective arguing that neoantigen prediction should be reframed as vaccine-ready AI: contamination-controlled, source-aware top-k ranking under incomplete-label uncertainty. The article would introduce NEO-PRIOR, a practical reporting checklist covering exact and near overlap, HLA and source holdouts, AUPRC, top-k precision, enrichment over prevalence, calibration and OOD-aware abstention.

The manuscript would not present unpublished original results. It would synthesize recent clinical and computational developments and propose a field standard for evaluating AI systems that select personalized cancer vaccine candidates.

I attach a short synopsis and proposed figure plan and would be grateful for your advice on whether this topic would be of interest.

Sincerely,

Seungho Cook

## Nature Cancer

Subject: Presubmission enquiry: Perspective on vaccine-ready benchmarks for neoantigen AI

Dear Editors,

I am writing to ask whether Nature Cancer would consider a Perspective provisionally titled "Neoantigen prediction needs vaccine-ready benchmarks, not another AUROC leaderboard."

The central argument is that personalized cancer vaccines have made neoantigen selection clinically consequential, but benchmark practices have not kept pace with deployment. A model's clinically relevant output is a short candidate list for manufacture or immune testing, not a genome-wide AUROC under overlap-sensitive retrospective splits.

We propose to introduce NEO-PRIOR, a reporting checklist for vaccine-ready neoantigen AI that covers public-corpus overlap, exact and near holdouts, source and HLA shift, incomplete negative labels, AUPRC, top-k precision, enrichment, calibration and OOD-aware abstention.

The Perspective would be balanced: public predictors are powerful and useful, but evaluation standards now determine whether their performance translates to vaccine decision-making.

Sincerely,

Seungho Cook
"""
    (OUT / "HIGH_IMPACT_EDITOR_EMAILS.md").write_text(text)


def main() -> None:
    ensure_dirs()
    fig0 = make_flagship_figure()
    write_high_impact_strategy(fig0)
    write_target_matrix()
    write_neoprior_checklist()
    write_high_impact_synopsis()
    write_email()
    print(f"[high-impact-review] wrote {OUT}")


if __name__ == "__main__":
    main()
