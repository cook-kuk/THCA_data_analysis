#!/usr/bin/env python3
"""Build a unified two-paper impact dossier.

Paper 1: high-impact standards Perspective.
Paper 2: benchmark/validation implementation paper.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "project/results/cross_neo_two_paper_strategy_2026_05_10"
FIG = OUT / "figures"
REVIEW = REPO / "project/results/cross_neo_review_2026_05_10"
METHODS = REPO / "project/results/cross_neo_methods_benchmark_2026_05_10"


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def add_box(ax, x, y, w, h, text, color, fs=10):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        facecolor=color,
        edgecolor="none",
        alpha=0.96,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color="white", fontsize=fs, weight="bold")


def add_arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=16, linewidth=2.2, color="#5f6770"))


def make_impact_map() -> Path:
    path = FIG / "two_paper_impact_map.png"
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Two-paper impact architecture", fontsize=18, weight="bold", color="#17202a")
    ax.text(0.04, 0.895, "One paper defines the field standard; the second implements and stress-tests it.", fontsize=11, color="#4b5563")

    add_box(ax, 0.06, 0.62, 0.36, 0.17, "Paper 1\nStandards Perspective\nNEO-PRIOR", "#5b5b8a", 12)
    add_box(ax, 0.58, 0.62, 0.36, 0.17, "Paper 2\nBenchmark / Validation Article\nCROSS-Neo", "#2f6f9f", 12)
    add_arrow(ax, (0.43, 0.705), (0.57, 0.705))

    ax.text(
        0.08,
        0.46,
        "Claim:\nPersonalized cancer vaccines need\nvaccine-ready AI standards, not\nAUROC leaderboards.",
        fontsize=10.5,
        color="#17202a",
        va="top",
    )
    ax.text(
        0.60,
        0.46,
        "Claim:\nCROSS-Neo is a contamination-\ncontrolled implementation and\nstress-test of that standard.",
        fontsize=10.5,
        color="#17202a",
        va="top",
    )

    add_box(ax, 0.08, 0.22, 0.30, 0.11, "Nature Medicine /\nNature Cancer", "#9a3e3e", 10)
    add_box(ax, 0.60, 0.22, 0.30, 0.11, "Nature Biomedical Engineering /\nNature Machine Intelligence", "#4b8f73", 9.5)

    ax.text(
        0.08,
        0.10,
        "Result: category leadership before model publication.",
        fontsize=11.5,
        weight="bold",
        color="#17202a",
    )
    ax.text(
        0.60,
        0.10,
        "Result: method becomes the first implementation of the standard.",
        fontsize=11.5,
        weight="bold",
        color="#17202a",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def write_master_strategy(fig_path: Path) -> None:
    text = f"""# CROSS-Neo Two-Paper High-Impact Strategy

## Decision

Yes: this should be **two papers**, not one.

The impact goes up because the papers do different jobs:

1. **Paper 1 defines the field standard.**
2. **Paper 2 implements the standard.**

Trying to put both into one paper weakens both: the Review becomes self-promotional, and the Methods paper becomes underpowered as a universal model claim. Separating them makes the story cleaner and higher-impact.

## Paper 1: Standards Perspective

### Role
Set the field agenda.

### Best Title
**Vaccine-ready AI for personalized cancer immunotherapy**

### Stronger Provocative Title
**Neoantigen prediction needs vaccine-ready benchmarks, not another AUROC leaderboard**

### Core Contribution
NEO-PRIOR: a minimum reporting and evaluation standard for neoantigen prioritization.

### Target Order
1. Nature Medicine Perspective / Comment
2. Nature Cancer Perspective
3. Nature Reviews Immunology
4. Nature Reviews Clinical Oncology
5. Nature Biomedical Engineering Perspective
6. Cancer Cell / Cancer Discovery

### What It Must Not Do
- Do not center CROSS-Neo.
- Do not present unpublished model results as evidence.
- Do not claim public predictors are bad.
- Do not claim external validation.

### Winning Sentence
The next advance in personalized cancer vaccines may come less from predicting more HLA binders than from trusting the shortlist that decides what gets manufactured.

## Paper 2: Benchmark / Validation Implementation

### Role
Show that the standard can be operationalized.

### Best Title
**Benchmarking vaccine-ready AI systems for neoantigen prioritization**

### Alternative Title
**CROSS-Neo: contamination-controlled pan-allele prioritization for personalized neoantigen vaccine triage**

### Core Contribution
A leakage-controlled benchmark and prioritization stack with explicit public-comparator separation, top-k primary metrics, source stress tests and OOD/abstention boundaries.

### Target Order
1. Nature Biomedical Engineering Article
2. Nature Machine Intelligence Analysis
3. Cell Reports Medicine
4. Patterns
5. JITC
6. npj Precision Oncology

### What It Must Not Do
- Do not claim external validation.
- Do not claim quantum advantage.
- Do not hide source-heldout failures.
- Do not present product-demo testset-aware scores as clean generalization.

## Why This Raises Impact

| Problem | One-paper approach | Two-paper approach |
|---|---|---|
| Review seems self-serving | high risk | avoided |
| Method seems underpowered | high risk | framed as implementation/stress-test |
| Public predictor criticism | defensive | becomes field-standard discussion |
| Small-n issue | central weakness | handled through benchmark/abstention framing |
| Business value | product claim only | category leadership + demo |

## Launch Order

1. Send Paper 1 presubmission to Nature Medicine and Nature Cancer.
2. Use the response to sharpen the standards framing.
3. In parallel finalize Paper 2 benchmark package.
4. Submit Paper 2 only after Paper 1 pitch is in motion, so CROSS-Neo can be framed as implementing NEO-PRIOR.

## Shared Claim Boundary

- Review/Perspective: field standard, no new original model claims.
- Methods/Benchmark: internal locked benchmark and implementation, no external validation.
- Product/demo: retrospective triage artifact, not manuscript validation.

## Impact Map

`{fig_path}`
"""
    (OUT / "TWO_PAPER_HIGH_IMPACT_STRATEGY.md").write_text(text)


def write_submission_matrix() -> None:
    rows = [
        {
            "paper": "Paper 1",
            "role": "Standards Perspective",
            "primary_title": "Vaccine-ready AI for personalized cancer immunotherapy",
            "primary_target": "Nature Medicine",
            "secondary_target": "Nature Cancer",
            "fallback": "Nature Reviews Immunology / Nature Reviews Clinical Oncology / Trends in Cancer",
            "main_asset": str(REVIEW / "HIGH_IMPACT_NATURE_MEDICINE_NATURE_CANCER_SYNOPSIS.md"),
            "central_deliverable": "NEO-PRIOR reporting standard",
            "claim_boundary": "field synthesis and standards; no original CROSS-Neo claim",
        },
        {
            "paper": "Paper 2",
            "role": "Benchmark / Validation Article",
            "primary_title": "Benchmarking vaccine-ready AI systems for neoantigen prioritization",
            "primary_target": "Nature Biomedical Engineering",
            "secondary_target": "Nature Machine Intelligence",
            "fallback": "Cell Reports Medicine / Patterns / JITC",
            "main_asset": str(METHODS / "NATURE_BIOMEDICAL_ENGINEERING_SYNOPSIS.md"),
            "central_deliverable": "CROSS-Neo benchmarked implementation",
            "claim_boundary": "internal locked benchmark; no external validation or quantum advantage",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "two_paper_submission_matrix.tsv", sep="\t", index=False)


def write_cover_pitch() -> None:
    text = """# Two-Paper Cover Pitch

## Category Pitch

We are not proposing another neoantigen predictor. We are proposing a field standard for vaccine-ready AI in personalized cancer immunotherapy, and a companion implementation that stress-tests whether such a standard can be operationalized.

## Paper 1 Pitch

Personalized neoantigen vaccines have made antigen selection clinically consequential, but the AI layer that selects vaccine targets is still evaluated like a retrospective prediction leaderboard. We propose NEO-PRIOR, a reporting standard for vaccine-ready neoantigen AI that makes top-k utility, contamination control, source shift, ambiguous negatives, calibration and OOD abstention central to evaluation.

## Paper 2 Pitch

CROSS-Neo is the companion implementation: a contamination-controlled pan-allele benchmark and prioritization stack evaluated under locked split families, public-comparator separation and source-heldout stress. It is deliberately framed as internal locked validation and failure-mode analysis, not as external clinical validation.

## Why Editors Should Care

- Personalized cancer vaccines are clinically active again.
- AI candidate selection is becoming a therapeutic decision layer.
- Current benchmark practices are not enough for that decision layer.
- The proposed standard is useful even beyond CROSS-Neo.
- The companion implementation shows feasibility and exposes real deployment limits.
"""
    (OUT / "TWO_PAPER_COVER_PITCH.md").write_text(text)


def write_action_checklist() -> None:
    rows = [
        {"order": 1, "action": "Send Paper 1 Nature Medicine presubmission synopsis", "file": str(REVIEW / "HIGH_IMPACT_EDITOR_EMAILS.md"), "status": "ready"},
        {"order": 2, "action": "Send Paper 1 Nature Cancer presubmission synopsis", "file": str(REVIEW / "HIGH_IMPACT_NATURE_MEDICINE_NATURE_CANCER_SYNOPSIS.md"), "status": "ready"},
        {"order": 3, "action": "Identify senior immuno-oncology coauthor/advisor", "file": "not local", "status": "needed"},
        {"order": 4, "action": "Polish Fig 0 as graphical abstract", "file": str(REVIEW / "figures/fig0_high_impact_grand_challenge_overview.png"), "status": "draft"},
        {"order": 5, "action": "Finalize Paper 2 NBE presubmission", "file": str(METHODS / "NATURE_BIOMEDICAL_ENGINEERING_SYNOPSIS.md"), "status": "ready"},
        {"order": 6, "action": "Strengthen Paper 2 with additional competitor manifest if time allows", "file": str(METHODS / "METHODS_BENCHMARK_MANUSCRIPT_SKELETON.md"), "status": "optional"},
        {"order": 7, "action": "Keep product demo separate from clean manuscript claims", "file": str(REPO / "project/results/cross_neo_product_demo_2026_05_10"), "status": "ready"},
    ]
    pd.DataFrame(rows).to_csv(OUT / "two_paper_action_checklist.tsv", sep="\t", index=False)


def main() -> None:
    ensure_dirs()
    fig_path = make_impact_map()
    write_master_strategy(fig_path)
    write_submission_matrix()
    write_cover_pitch()
    write_action_checklist()
    print(f"[two-paper-impact] wrote {OUT}")


if __name__ == "__main__":
    main()
