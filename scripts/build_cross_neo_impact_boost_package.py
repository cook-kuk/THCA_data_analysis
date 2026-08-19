#!/usr/bin/env python3
"""Create an impact-boost package for the CROSS-Neo two-paper strategy.

This upgrades the story from two manuscripts to a field-standard platform:
NEO-PRIOR consensus checklist + benchmark challenge + implementation paper.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "project/results/cross_neo_two_paper_strategy_2026_05_10"
OUT = ROOT / "impact_boost"
FIG = OUT / "figures"
REVIEW = REPO / "project/results/cross_neo_review_2026_05_10"
METHODS = REPO / "project/results/cross_neo_methods_benchmark_2026_05_10"
PRODUCT = REPO / "project/results/cross_neo_product_demo_2026_05_10"


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def add_box(ax, x, y, w, h, text, color, fs=10):
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


def add_arrow(ax, p1, p2):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=15, color="#5d6570", linewidth=2))


def fig_three_layer_model() -> Path:
    path = FIG / "fig_impact_three_layer_model.png"
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Impact upgrade: from papers to field infrastructure", fontsize=18, weight="bold", color="#17202a")
    ax.text(0.04, 0.89, "The high-impact version is a standard, a benchmark, and an implementation.", fontsize=11, color="#4b5563")

    layers = [
        ("1. Standard\nNEO-PRIOR", "Minimum reporting\nand evaluation checklist", "#5b5b8a"),
        ("2. Benchmark\nNeoBench-Vax", "Public score matrix,\nlocked splits, model cards", "#2f6f9f"),
        ("3. Implementation\nCROSS-Neo", "Pan-allele prioritizer\nand abstention system", "#4b8f73"),
    ]
    xs = [0.06, 0.37, 0.68]
    for (title, detail, color), x in zip(layers, xs):
        add_box(ax, x, 0.60, 0.25, 0.15, title, color, 12)
        ax.text(x + 0.125, 0.46, detail, ha="center", va="center", fontsize=10.5, color="#17202a")
    add_arrow(ax, (0.32, 0.675), (0.36, 0.675))
    add_arrow(ax, (0.63, 0.675), (0.67, 0.675))

    ax.text(0.08, 0.25, "Paper 1", fontsize=12, weight="bold", color="#17202a")
    ax.text(0.08, 0.20, "Nature Medicine / Nature Cancer\nPerspective: define NEO-PRIOR", fontsize=10, color="#17202a")
    ax.text(0.40, 0.25, "Public Asset", fontsize=12, weight="bold", color="#17202a")
    ax.text(0.40, 0.20, "benchmark registry + score matrix\n+ provenance/model cards", fontsize=10, color="#17202a")
    ax.text(0.70, 0.25, "Paper 2", fontsize=12, weight="bold", color="#17202a")
    ax.text(0.70, 0.20, "Nature Biomedical Engineering /\nNMI: implement and stress-test", fontsize=10, color="#17202a")
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def fig_stakeholder_map() -> Path:
    path = FIG / "fig_editorial_stakeholder_map.png"
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Who cares if NEO-PRIOR exists?", fontsize=18, weight="bold", color="#17202a")

    center = (0.40, 0.48, 0.20, 0.14)
    add_box(ax, *center, "NEO-PRIOR\nvaccine-ready AI\nstandard", "#5b5b8a", 11)
    stakeholders = [
        ("Clinical trialists\nselect endpoints\nand candidates", 0.07, 0.70, "#2f6f9f"),
        ("Vaccine companies\nprioritize manufacture\nslots", 0.38, 0.75, "#4b8f73"),
        ("AI reviewers\navoid leakage\nleaderboards", 0.70, 0.70, "#b36b2c"),
        ("Regulators\nneed provenance\nand calibration", 0.08, 0.24, "#9a3e3e"),
        ("Immunologists\nseparate binding\nfrom recognition", 0.39, 0.14, "#326b7c"),
        ("Patients\nneed reliable\nshortlists", 0.71, 0.24, "#7d6f2f"),
    ]
    for text, x, y, color in stakeholders:
        add_box(ax, x, y, 0.20, 0.12, text, color, 9)
        add_arrow(ax, (x + 0.10, y + 0.06), (0.50, 0.55))
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def write_impact_memo(fig1: Path, fig2: Path) -> None:
    text = f"""# Impact Boost Memo

## Upgrade

The story should not be "two papers about CROSS-Neo."

The upgraded story is:

> We define the reporting standard for vaccine-ready AI in personalized cancer immunotherapy, release a benchmark infrastructure, and then show CROSS-Neo as the first implementation.

## New Three-Layer Architecture

1. **NEO-PRIOR Standard**
   - Field-level checklist and reporting standard.
   - Paper 1 target: Nature Medicine / Nature Cancer.

2. **NeoBench-Vax Benchmark**
   - Public score matrix, locked split registry, model cards, data provenance, overlap audit.
   - This makes the work look like infrastructure, not just an opinion.

3. **CROSS-Neo Implementation**
   - Benchmark/validation system that operationalizes NEO-PRIOR.
   - Paper 2 target: Nature Biomedical Engineering / Nature Machine Intelligence.

## Why This Raises Impact

- Standards papers are bigger than model papers.
- Benchmark infrastructure is more reusable than one dataset result.
- CROSS-Neo becomes the first proof-of-concept rather than the only reason the work exists.
- Source-heldout failures become part of honest deployment standards.
- The product story becomes credible because it is backed by a field-standard framework.

## New Names

- **NEO-PRIOR**: Neoantigen Prioritization Reporting and Evaluation Standard.
- **NeoBench-Vax**: benchmark registry and public score matrix for vaccine-ready neoantigen prioritization.
- **CROSS-Neo**: contamination-controlled pan-allele prioritization implementation.

## High-Impact Figures

- Three-layer model: `{fig1}`
- Stakeholder map: `{fig2}`

## One-Line Editorial Hook

Personalized cancer vaccines now need a benchmark standard for the AI layer that decides what gets manufactured.
"""
    (OUT / "IMPACT_BOOST_MEMO.md").write_text(text)


def write_neoprior_statement() -> None:
    text = """# NEO-PRIOR Consensus-Style Statement

## Full Name
NEO-PRIOR: Neoantigen Prioritization Reporting and Evaluation Standard

## Purpose
NEO-PRIOR defines minimum reporting requirements for computational systems that rank neoantigen candidates for personalized vaccine manufacture, immune testing or clinical prioritization.

## Scope
NEO-PRIOR applies to:
- binding/presentation predictors used in vaccine pipelines;
- immunogenicity predictors;
- pan-allele neoantigen rankers;
- public benchmark papers;
- clinical vaccine candidate-selection pipelines;
- product/demo triage systems.

## Minimum Items

1. State the task: binding, presentation, immunogenicity, vaccine triage or benchmark comparison.
2. Report exact peptide-HLA and exact peptide overlap.
3. Report near-peptide and source-protein-window overlap.
4. Disclose public predictor scores used as features, comparators or product assists.
5. Use train-only retrieval, preprocessing, scaling, calibration and threshold selection.
6. Include exact, near, source/study, HLA allele/supertype and time/assay holdouts when feasible.
7. Treat negative labels as ambiguous unless true non-immunogenicity is experimentally established.
8. Report AUPRC, top-k precision, enrichment over prevalence and recall@k as primary metrics.
9. Report AUROC as secondary.
10. Report calibration, Brier/ECE and OOD/abstention coverage.
11. Publish data provenance and claim-boundary tables.
12. Separate clean scientific validation from product/demo retrospective triage.

## Why It Matters
Neoantigen vaccines depend on a short candidate list. The wrong benchmark can reward memorization, HLA shortcuts or source-specific bias, while the right benchmark can reveal when a system should abstain.

## Intended Use
NEO-PRIOR is a proposed reporting checklist, not a regulatory standard. It is designed to make neoantigen AI papers auditable, comparable and clinically interpretable.
"""
    (OUT / "NEO_PRIOR_CONSENSUS_STYLE_STATEMENT.md").write_text(text)


def write_neobench_plan() -> None:
    text = """# NeoBench-Vax Benchmark Challenge Plan

## Purpose
NeoBench-Vax turns the two-paper package into reusable infrastructure.

## Assets
- canonical benchmark registry;
- public predictor score matrix;
- locked split manifests;
- exact/near/source/HLA/study/time holdouts;
- public overlap audit;
- model cards;
- data provenance cards;
- top-k and AUPRC reporting templates;
- failure-mode table.

## Challenge Tracks

### Track A: Clean Comparator Track
No public predictor scores as training features. Public predictors are comparators only.

### Track B: Product-Assisted Track
Public predictor scores allowed, but clearly labeled as assisted triage, not clean scientific validation.

### Track C: Source-Stress Track
Models are ranked by source-heldout top-k enrichment and abstention behavior.

### Track D: OOD/Abstention Track
Models are rewarded for abstaining when source/HLA/length/retrieval evidence is unsafe.

## Leaderboard Metrics
1. AUPRC
2. top5 precision
3. top10 precision
4. enrichment@10
5. recall@10
6. Brier/ECE
7. abstention coverage versus precision
8. source-heldout collapse penalty

## Why This Helps the Papers
Paper 1 proposes the reporting standard. NeoBench-Vax makes the standard concrete. Paper 2 shows CROSS-Neo as the first system evaluated through the benchmark.
"""
    (OUT / "NEOBENCH_VAX_CHALLENGE_PLAN.md").write_text(text)


def write_editor_hook_matrix() -> None:
    rows = [
        {
            "editor_type": "clinical oncology",
            "hook": "AI now decides which neoantigens are manufactured for personalized vaccines.",
            "asset": "NEO-PRIOR standard + clinical vaccine framing",
        },
        {
            "editor_type": "cancer biology/immunology",
            "hook": "Binding prediction and immunogenicity recognition must be benchmarked as distinct biological layers.",
            "asset": "predictor taxonomy and benchmark failure figure",
        },
        {
            "editor_type": "AI/ML",
            "hook": "Overlap-sensitive biomedical leaderboards can invert under top-k and source-heldout evaluation.",
            "asset": "NeoBench-Vax locked splits and competitor score matrix",
        },
        {
            "editor_type": "biomedical engineering",
            "hook": "Personalized vaccine design needs validation and deployment criteria for candidate-selection systems.",
            "asset": "CROSS-Neo validation-system paper",
        },
        {
            "editor_type": "industry/product",
            "hook": "The framework converts model scores into auditable candidate queues with abstention boundaries.",
            "asset": "customer-safe queue and product model card",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "editor_hook_matrix.tsv", sep="\t", index=False)


def write_coauthor_brief() -> None:
    text = """# Senior Coauthor / Advisor Recruitment Brief

## Why We Need One
Nature Medicine, Nature Cancer, Nature Reviews and Cancer Cell routes become more plausible if a senior immuno-oncology or cancer-vaccine expert is attached. The expert does not need to own CROSS-Neo; they strengthen the clinical and immunology framing of NEO-PRIOR.

## Ideal Profiles

1. Personalized neoantigen vaccine trialist
   - melanoma, pancreatic cancer, renal cell carcinoma or mRNA vaccine trial experience.

2. Tumor immunologist
   - T cell recognition, antigen processing, HLA biology, cancer vaccine expertise.

3. Clinical computational oncologist
   - experience translating ML predictions into treatment or trial decisions.

4. Benchmark/reproducibility expert
   - strengthens NEO-PRIOR and NeoBench-Vax credibility.

## Ask
We are preparing a Perspective proposing NEO-PRIOR, a reporting standard for vaccine-ready neoantigen AI, plus a companion benchmark implementation. We would value senior input on clinical framing, reporting checklist, and claim boundaries.

## What To Send
- HIGH_IMPACT_NATURE_MEDICINE_NATURE_CANCER_SYNOPSIS.md
- fig0_high_impact_grand_challenge_overview.png
- NEO_PRIOR_CONSENSUS_STYLE_STATEMENT.md
- TWO_PAPER_COVER_PITCH.md
"""
    (OUT / "SENIOR_COAUTHOR_RECRUITMENT_BRIEF.md").write_text(text)


def write_title_bank() -> None:
    rows = [
        {"paper": "Paper 1", "rank": 1, "title": "Vaccine-ready AI for personalized cancer immunotherapy", "tone": "Nature Medicine"},
        {"paper": "Paper 1", "rank": 2, "title": "Neoantigen prediction needs vaccine-ready benchmarks, not another AUROC leaderboard", "tone": "Nature Cancer"},
        {"paper": "Paper 1", "rank": 3, "title": "The missing benchmark standard for personalized neoantigen vaccines", "tone": "Nature Reviews"},
        {"paper": "Paper 1", "rank": 4, "title": "From neoantigen prediction to vaccine-ready AI", "tone": "broad review"},
        {"paper": "Paper 2", "rank": 1, "title": "Benchmarking vaccine-ready AI systems for neoantigen prioritization", "tone": "Nature Biomedical Engineering"},
        {"paper": "Paper 2", "rank": 2, "title": "Contamination-controlled benchmarking of neoantigen prioritization models", "tone": "Nature Machine Intelligence"},
        {"paper": "Paper 2", "rank": 3, "title": "CROSS-Neo: pan-allele neoantigen prioritization under leakage, source shift and label uncertainty", "tone": "Patterns/JITC"},
        {"paper": "Paper 2", "rank": 4, "title": "A benchmark and model-card framework for personalized neoantigen vaccine candidate selection", "tone": "Patterns"},
    ]
    pd.DataFrame(rows).to_csv(OUT / "high_impact_title_bank.tsv", sep="\t", index=False)


def main() -> None:
    ensure_dirs()
    fig1 = fig_three_layer_model()
    fig2 = fig_stakeholder_map()
    write_impact_memo(fig1, fig2)
    write_neoprior_statement()
    write_neobench_plan()
    write_editor_hook_matrix()
    write_coauthor_brief()
    write_title_bank()
    print(f"[impact-boost] wrote {OUT}")


if __name__ == "__main__":
    main()
