#!/usr/bin/env python3
"""Build a max-impact strategy package for the CROSS-Neo two-paper program.

The package deliberately separates:
1) field-level standard setting (NEO-PRIOR),
2) reusable benchmark infrastructure (NeoBench-Vax), and
3) the CROSS-Neo implementation.

It avoids claims of external validation, quantum advantage, or clean public
predictor superiority. The goal is to prepare editorial, reviewer, clinical,
and business-facing scaffolding without touching protected manuscript prose.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from textwrap import dedent


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/cross_neo_two_paper_strategy_2026_05_10/impact_max"
FIG = OUT / "figures"
TRANSFER = ROOT / "project/results/transfer_packages"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dedent(text).strip() + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def draw_box(ax, xy, wh, title, body, fc="#f8fafc", ec="#334155", title_color="#0f172a"):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    rect = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.3,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(rect)
    ax.text(x + 0.03, y + h - 0.08, title, fontsize=12.5, weight="bold", color=title_color, va="top")
    ax.text(x + 0.03, y + h - 0.16, body, fontsize=9.5, color="#334155", va="top", linespacing=1.25)


def save_editor_one_screen() -> None:
    import matplotlib.pyplot as plt

    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(13.5, 7.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("#ffffff")
    ax.text(
        0.05,
        0.95,
        "Max-Impact Story: Vaccine-Ready AI Needs a Standard, a Benchmark, and an Implementation",
        fontsize=16,
        weight="bold",
        color="#0f172a",
        va="top",
    )
    ax.text(
        0.05,
        0.90,
        "Strict claim boundary: internal/locked evaluation only; no external validation or quantum-advantage claim.",
        fontsize=10.5,
        color="#475569",
        va="top",
    )
    draw_box(
        ax,
        (0.05, 0.60),
        (0.27, 0.23),
        "1. NEO-PRIOR",
        "Reporting standard for\nneoantigen prioritization AI:\noverlap, leakage, top-k,\ncalibration, OOD, provenance.",
        fc="#ecfeff",
        ec="#0891b2",
    )
    draw_box(
        ax,
        (0.365, 0.60),
        (0.27, 0.23),
        "2. NeoBench-Vax",
        "Benchmark registry:\nlocked splits, public-score matrix,\nmodel cards, failure modes,\nsource-stress and abstention tracks.",
        fc="#f0fdf4",
        ec="#16a34a",
    )
    draw_box(
        ax,
        (0.68, 0.60),
        (0.27, 0.23),
        "3. CROSS-Neo",
        "First implementation:\ncontamination-controlled,\nOOD-aware, source-stress tested,\nbusiness triage ready.",
        fc="#fff7ed",
        ec="#ea580c",
    )
    ax.annotate("", xy=(0.36, 0.715), xytext=(0.32, 0.715), arrowprops=dict(arrowstyle="->", lw=2, color="#64748b"))
    ax.annotate("", xy=(0.675, 0.715), xytext=(0.635, 0.715), arrowprops=dict(arrowstyle="->", lw=2, color="#64748b"))
    draw_box(
        ax,
        (0.05, 0.26),
        (0.42, 0.22),
        "Why Editors Care",
        "The bottleneck is no longer only predicting binders.\nThe clinical decision is a short manufacturing list.\nBad benchmarks reward memorization, leakage and source bias.",
        fc="#f8fafc",
        ec="#475569",
    )
    draw_box(
        ax,
        (0.53, 0.26),
        (0.42, 0.22),
        "Why This Is Defensible",
        "Public predictors stay comparators, not clean features.\nAUPRC/top-k lead; AUROC is secondary.\nFailures are reported as deployment boundaries, not hidden.",
        fc="#f8fafc",
        ec="#475569",
    )
    ax.text(
        0.5,
        0.10,
        "Editorial one-liner: Personalized cancer vaccines need auditable AI for the shortlist that decides what gets manufactured.",
        ha="center",
        fontsize=12,
        weight="bold",
        color="#0f172a",
    )
    fig.savefig(FIG / "fig_editor_one_screen_story.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def save_ecosystem() -> None:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.axis("off")
    fig.patch.set_facecolor("#ffffff")
    circles = [
        (0, 0.45, 0.43, "#ecfeff", "#0891b2", "NEO-PRIOR\nstandard"),
        (-0.48, -0.25, 0.43, "#f0fdf4", "#16a34a", "NeoBench-Vax\nbenchmark"),
        (0.48, -0.25, 0.43, "#fff7ed", "#ea580c", "CROSS-Neo\nimplementation"),
    ]
    for x, y, r, fc, ec, label in circles:
        ax.add_patch(patches.Circle((x, y), r, facecolor=fc, edgecolor=ec, linewidth=2.2, alpha=0.88))
        ax.text(x, y, label, ha="center", va="center", fontsize=15, weight="bold", color="#0f172a")
    ax.text(0, -0.02, "vaccine-ready AI\nfor clinical prioritization", ha="center", va="center", fontsize=12, weight="bold", color="#334155")
    annotations = [
        (-0.92, 0.88, "Perspective / standard\nfield-level contribution"),
        (0.86, 0.88, "Reviewer defense\nclaim-boundary table"),
        (-1.0, -0.86, "Reusable infrastructure\nleaderboard + model cards"),
        (0.78, -0.88, "Business demo\ntriage + abstention"),
    ]
    for x, y, txt in annotations:
        ax.text(x, y, txt, ha="center", fontsize=10.5, color="#334155")
    ax.text(0, 1.08, "Impact-Max Ecosystem", ha="center", fontsize=18, weight="bold", color="#0f172a")
    ax.text(0, 0.99, "A standards-first program, not a single-model performance claim", ha="center", fontsize=11, color="#475569")
    fig.savefig(FIG / "fig_neoprior_neobench_crossneo_ecosystem.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def save_value_chain() -> None:
    import matplotlib.pyplot as plt

    steps = [
        ("Tumor\nvariants", "candidate universe"),
        ("HLA + peptide\ncontext", "presentation constraints"),
        ("Prioritization\nAI", "shortlist risk"),
        ("Wet-lab\nselection", "assay budget"),
        ("Manufacture", "vaccine slots"),
        ("Clinical\nmonitoring", "immune readout"),
    ]
    colors = ["#e0f2fe", "#ecfeff", "#fef9c3", "#ffedd5", "#fce7f3", "#ede9fe"]
    fig, ax = plt.subplots(figsize=(14, 5.8))
    ax.set_xlim(0, len(steps))
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.1, 0.94, "Clinical and Business Value Chain", fontsize=16, weight="bold", color="#0f172a", va="top")
    ax.text(0.1, 0.88, "The commercial product is decision support for scarce assay/manufacturing slots, not a claim of universal biology.", fontsize=10.5, color="#475569", va="top")
    for i, ((title, body), color) in enumerate(zip(steps, colors)):
        draw_box(ax, (i + 0.08, 0.36), (0.78, 0.30), title, body, fc=color, ec="#64748b")
        if i < len(steps) - 1:
            ax.annotate("", xy=(i + 1.02, 0.51), xytext=(i + 0.88, 0.51), arrowprops=dict(arrowstyle="->", lw=1.8, color="#64748b"))
    ax.text(3.0, 0.18, "NEO-PRIOR asks: can the model rank the shortlist, explain risk, and abstain when evidence is unsafe?", ha="center", fontsize=12, weight="bold", color="#0f172a")
    fig.savefig(FIG / "fig_clinical_translation_value_chain.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def save_risk_map() -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    labels = [
        "external validation",
        "quantum advantage",
        "public predictor leakage",
        "source/study shortcut",
        "near-peptide memorization",
        "negative-label ambiguity",
        "top-k instability",
        "calibration/OOD",
    ]
    risk = np.array([0.85, 0.75, 0.70, 0.82, 0.76, 0.72, 0.65, 0.60])
    control = np.array([0.55, 0.72, 0.90, 0.62, 0.78, 0.68, 0.58, 0.73])
    fig, ax = plt.subplots(figsize=(12, 6.8))
    y = np.arange(len(labels))
    ax.barh(y + 0.18, risk, height=0.32, color="#fb7185", label="uncontrolled claim risk")
    ax.barh(y - 0.18, control, height=0.32, color="#38bdf8", label="NEO-PRIOR control strength")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("relative score")
    ax.set_title("Internal Claim-Boundary Risk Map", fontsize=15, weight="bold")
    ax.legend(loc="lower right")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.savefig(FIG / "fig_claim_boundary_risk_map.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    TRANSFER.mkdir(parents=True, exist_ok=True)

    write_text(
        OUT / "MAX_IMPACT_STRATEGY.md",
        """
        # CROSS-Neo Max-Impact Strategy

        ## The Upgrade
        The strongest framing is not "we built a better neoantigen predictor."

        The strongest framing is:

        > Personalized cancer vaccines need auditable AI for the shortlist that decides what gets manufactured.

        This makes the program bigger than a model. It becomes a field-standard package with three linked assets:

        1. **NEO-PRIOR**: a reporting and evaluation standard.
        2. **NeoBench-Vax**: a public benchmark and leaderboard infrastructure.
        3. **CROSS-Neo**: the first implementation tested under that standard.

        ## Paper Architecture

        ### Paper 1: Field-Standard / Perspective
        **Working title:** Vaccine-ready AI for personalized cancer immunotherapy

        **Target story:** cancer vaccine AI is at risk of optimizing the wrong benchmark. The field needs a minimum reporting standard before model improvements can be clinically trusted.

        **Target journals:** Nature Medicine, Nature Cancer, Nature Reviews Clinical Oncology, Nature Reviews Immunology, Cancer Cell.

        **Do not center:** unpublished CROSS-Neo performance.

        ### Paper 2: Benchmark / Implementation
        **Working title:** Benchmarking vaccine-ready AI systems for neoantigen prioritization

        **Target story:** a leakage-safe benchmark, claim-boundary framework, and first implementation show how neoantigen AI should be evaluated before clinical or product use.

        **Target journals:** Nature Biomedical Engineering, Nature Machine Intelligence, Cell Reports Medicine, Patterns, JITC, npj Precision Oncology.

        **Do not claim:** external validation, quantum advantage, or universal superiority over public predictors.

        ## Impact Lever
        A model paper competes with every other model paper. A standard and benchmark paper defines the rules by which model papers are judged.

        ## Positioning Against Existing Predictors
        Public predictors such as MHCflurry, NetMHCpan, BigMHC, PRIME and MixMHCpred are not enemies. They become benchmark comparators or product-assisted inputs under explicit disclosure.

        The claim is not "we replace binding predictors." The claim is:

        > Binding and presentation predictors are necessary but insufficient for vaccine-ready prioritization unless overlap, source bias, top-k precision, OOD risk and negative-label ambiguity are reported.

        ## Business Translation
        The business product is not a magic immunogenicity oracle. It is an auditable triage layer that:

        - concentrates likely positives into a smaller assay/manufacturing queue;
        - flags when a sample is retrieval-supported versus clean-no-reference;
        - abstains under unsafe HLA/source/peptide-distribution shift;
        - separates clean scientific evaluation from product-assisted triage.

        ## Immediate Next Moves
        1. Recruit one senior cancer-vaccine or immuno-oncology advisor.
        2. Convert NEO-PRIOR into a two-page checklist.
        3. Convert NeoBench-Vax into a public challenge/leaderboard mockup.
        4. Make one editorial figure that shows standard + benchmark + implementation.
        5. Keep CROSS-Neo performance claims conservative until public-overlap and source-heldout limitations are resolved.
        """,
    )

    checklist_rows = [
        ("Task", "Define whether the model predicts binding, presentation, immunogenicity, vaccine triage, or benchmark ranking.", "required", "Prevents task drift."),
        ("Data provenance", "Report source dataset, study, patient, assay, date/year, HLA typing method, and candidate-generation route when available.", "required", "Needed for source/study shift."),
        ("Public-score disclosure", "List all public predictor scores used as features, comparators, filters, or product assists.", "required", "Prevents hidden comparator leakage."),
        ("Exact overlap", "Report exact peptide and exact peptide-HLA overlap between train, test, public corpora, and comparator training data where available.", "required", "Controls retrieval leakage."),
        ("Near overlap", "Report near-peptide, mutation-pair, source-window and TCR motif overlap where available.", "required", "Controls analog memorization."),
        ("Train-only pipeline", "Build retrieval, scaling, imputation, calibration, feature selection and thresholding inside train folds only.", "required", "Core leakage control."),
        ("Split contract", "Include exact, near, source/study, HLA allele/supertype, time/assay and retrieval-clean splits when feasible.", "required", "Tests robustness."),
        ("Primary metrics", "Use AUPRC, top-k precision, enrichment over prevalence and recall@k as primary metrics.", "required", "Matches vaccine shortlist use."),
        ("Secondary metrics", "Report AUROC as secondary, with prevalence baseline.", "required", "Avoids AUROC-only claims."),
        ("Calibration", "Report Brier, ECE and calibration plots.", "required", "Supports decision risk."),
        ("OOD/abstention", "Report OOD flags, abstention coverage, and precision after abstention.", "required", "Supports safe non-decision."),
        ("Negative labels", "Describe whether negatives are true non-immunogenic, untested, assay-limited, or unlabeled.", "required", "Avoids false certainty."),
        ("Source bias", "Report per-source prevalence, feature shift, score shift and source-heldout performance.", "required", "Controls dataset shortcuts."),
        ("HLA robustness", "Report per-HLA and HLA-heldout/supertype-heldout performance.", "required", "Controls allele shortcuts."),
        ("Top-k stability", "Report bootstrap or repeated split uncertainty for top5/top10.", "recommended", "Important for small n."),
        ("Model cards", "Publish model card with intended use, forbidden use, inputs, leakage controls and known failure modes.", "recommended", "Clinical/product transparency."),
        ("Claim boundary", "Separate internal validation, locked retrospective tests, source-heldout stress tests, and external prospective validation.", "required", "Prevents overclaim."),
        ("Assay budget", "Report performance at realistic assay/manufacturing queue sizes.", "recommended", "Business relevance."),
        ("Comparator fairness", "Do not train on public predictor scores when those predictors are clean comparators.", "required", "Fair benchmark."),
        ("Reproducibility", "Provide split manifests, row identifiers, seeds and scripts.", "recommended", "Auditability."),
    ]
    write_tsv(
        OUT / "NEO_PRIOR_FULL_CHECKLIST.tsv",
        [
            {"domain": d, "item": i, "level": l, "why_it_matters": w}
            for d, i, l, w in checklist_rows
        ],
        ["domain", "item", "level", "why_it_matters"],
    )

    write_text(
        OUT / "NEOBENCH_VAX_LEADERBOARD_GOVERNANCE.md",
        """
        # NeoBench-Vax Leaderboard and Governance

        ## Purpose
        NeoBench-Vax is a benchmark infrastructure for vaccine-ready neoantigen prioritization. It is designed to make models auditable under NEO-PRIOR rather than to reward a single headline AUROC.

        ## Tracks

        | Track | Name | Public predictor scores as features | Intended claim |
        |---|---|---:|---|
        | A | Clean Comparator | No | Scientific comparison under strict leakage control |
        | B | Product-Assisted Triage | Yes, disclosed | Retrospective decision-support utility only |
        | C | Source-Stress | No or disclosed by subtrack | Robustness under source/study/prevalence shift |
        | D | OOD/Abstention | No or disclosed by subtrack | Reliability and safe non-decision behavior |

        ## Ranking Metrics
        Primary leaderboard rank should be multi-objective:

        1. AUPRC over prevalence baseline.
        2. top10 precision and enrichment@10.
        3. top5 precision for manufacturing-scarce use.
        4. retrieval-clean-only AUPRC and top-k.
        5. source-heldout collapse penalty.
        6. calibration and abstention benefit.

        ## Model Card Requirements
        Every submission should state:

        - allowed inputs;
        - forbidden inputs;
        - public predictor use;
        - split contract;
        - training data provenance;
        - known overlap risks;
        - source/HLA/time generalization boundary;
        - whether negative labels are hard negatives or unlabeled.

        ## Governance Principle
        A model can lead a product-assisted track and still be ineligible for clean scientific claims. This separation is central to the benchmark.
        """,
    )

    write_text(
        OUT / "CLINICAL_BUSINESS_TRANSLATION_MEMO.md",
        """
        # Clinical and Business Translation Memo

        ## Product Wedge
        CROSS-Neo should be packaged as an auditable vaccine-candidate triage system, not as a standalone biological truth engine.

        ## Buyer/User Pain
        Personalized vaccine teams face too many candidates and too few assay/manufacturing slots. They need ranked shortlists, uncertainty flags, and a way to justify which candidates move forward.

        ## Practical Product Modules
        1. Candidate queue: top-k ranking with enrichment and prevalence baseline.
        2. Evidence tags: clean-no-reference, near-retrieval-supported, exact-retrieval-supported, OOD.
        3. Abstention mode: flag unsafe source/HLA/length/structure/retrieval cases.
        4. Comparator panel: optional public predictor-assisted view, clearly labeled.
        5. Audit report: NEO-PRIOR checklist, data provenance, overlap and claim-boundary table.

        ## Demo Claim
        Acceptable:

        > In internal retrospective testing, the product-assisted queue concentrated positives in the top candidate list better than any single available score in this dataset.

        Not acceptable:

        > Externally validated, universally superior, or quantum-advantaged cancer vaccine predictor.

        ## Business-Safe Differentiation
        The differentiator is not a secret score. It is the audit layer: leakage controls, source-risk labeling, top-k decision metrics, and abstention.
        """,
    )

    write_text(
        OUT / "EDITORIAL_TRIANGLE_PITCH.md",
        """
        # Editorial Triangle Pitch

        ## One-Sentence Pitch
        Personalized cancer vaccines need auditable AI for the shortlist that decides what gets manufactured.

        ## Why Now
        Neoantigen vaccine trials are increasingly real, but computational prioritization remains hard to compare because papers mix different tasks, overlaps, public predictors, negative-label assumptions and evaluation splits.

        ## What We Contribute
        - NEO-PRIOR: a reporting and evaluation standard.
        - NeoBench-Vax: a benchmark/challenge infrastructure.
        - CROSS-Neo: a first contamination-controlled implementation and stress-tested example.

        ## Why It Is High Impact
        The work changes the unit of contribution from a model to a standard. It gives editors a broad clinical AI governance story while giving method reviewers concrete benchmark artifacts.

        ## Editorial Guardrail
        The paper should be explicit that current results are internal/locked retrospective tests. Source-heldout collapse is reported as a deployment boundary and a reason for standardization, not hidden.
        """,
    )

    objection_rows = [
        ("This is just another MHC binding predictor.", "No. The proposed scope is vaccine-ready prioritization: overlap control, source stress, top-k ranking, calibration, OOD and auditability. Binding predictors are comparators or product-assisted inputs."),
        ("Small n makes the model claims weak.", "The main contribution is the reporting standard and benchmark contract. CROSS-Neo is a first implementation, with small-n uncertainty explicitly reported."),
        ("Public predictor scores may leak.", "Clean track forbids them as features. Product-assisted track permits them only with disclosure and no clean scientific claim."),
        ("Source-heldout performance collapses.", "This is a central finding and supports the NEO-PRIOR requirement for source-stress reporting and abstention."),
        ("Quantum branch is overclaimed.", "No quantum advantage is claimed. QK is a fixed fallback branch evaluated as an implementation component."),
        ("AUROC looks fine but top-k is unstable.", "AUPRC, top-k precision and enrichment are primary. AUROC is explicitly secondary."),
        ("Negatives are not true negatives.", "NEO-PRIOR requires negative-label ambiguity reporting and supports PU-aware evaluation."),
        ("How is this clinically useful without prospective validation?", "It is not presented as externally validated. The practical utility is retrospective triage/audit infrastructure and a prospective-ready evaluation framework."),
    ]
    write_tsv(
        OUT / "reviewer_objection_response_matrix.tsv",
        [{"objection": o, "response": r} for o, r in objection_rows],
        ["objection", "response"],
    )

    gap_rows = [
        ("Field confusion", "Different papers optimize different tasks.", "NEO-PRIOR task declaration", "Paper 1"),
        ("Leakage risk", "Exact/near public overlap can inflate performance.", "Overlap audit + train-only pipeline", "Paper 1/2"),
        ("Wrong metric", "AUROC can hide poor shortlist utility.", "AUPRC/top-k/enrichment primary", "Paper 1/2"),
        ("Source shortcut", "Models learn study/prevalence/assay artifacts.", "Source-stress track", "Paper 2"),
        ("Clinical uncertainty", "Models do not know when not to decide.", "OOD/abstention track", "Paper 2/product"),
        ("Comparator ambiguity", "Public predictors used as both features and baselines.", "Clean vs product-assisted tracks", "Paper 1/2/product"),
        ("Small-n instability", "Neoantigen immunogenicity datasets remain limited.", "Bootstrap, repeated splits, claim boundary", "Paper 2"),
        ("Business trust", "Customers need audit trail, not only score.", "Candidate queue + NEO-PRIOR report", "Product"),
    ]
    write_tsv(
        OUT / "high_impact_gap_to_asset_matrix.tsv",
        [{"gap": g, "why_it_matters": w, "asset": a, "paper_or_product": p} for g, w, a, p in gap_rows],
        ["gap", "why_it_matters", "asset", "paper_or_product"],
    )

    stakeholder_rows = [
        ("Nature Medicine / Nature Cancer editor", "Clinical AI governance story", "NEO-PRIOR + vaccine-ready AI framing"),
        ("Nature Biomedical Engineering editor", "Deployable biomedical AI system and validation infrastructure", "NeoBench-Vax + CROSS-Neo implementation"),
        ("Cancer vaccine trialist", "Justifiable candidate shortlist", "top-k queue, abstention, source-risk labels"),
        ("Computational immunology reviewer", "Fair comparison and leakage control", "split manifests, overlap audit, public-score disclosure"),
        ("Biotech partner", "Assay/manufacturing prioritization product", "customer-safe queue and audit report"),
        ("Skeptical ML reviewer", "Avoid overfit and target leakage", "nested/fold-safe methods, claim boundary, source-heldout stress"),
    ]
    write_tsv(
        OUT / "stakeholder_value_matrix.tsv",
        [{"stakeholder": s, "cares_about": c, "asset_to_show": a} for s, c, a in stakeholder_rows],
        ["stakeholder", "cares_about", "asset_to_show"],
    )

    title_rows = [
        ("standards", "Vaccine-ready AI for personalized cancer immunotherapy"),
        ("standards", "A reporting standard for neoantigen AI in personalized cancer vaccines"),
        ("standards", "From neoantigen prediction to vaccine-ready prioritization"),
        ("benchmark", "Benchmarking vaccine-ready AI systems for neoantigen prioritization"),
        ("benchmark", "Leakage-safe evaluation of AI systems for personalized neoantigen vaccines"),
        ("benchmark", "NeoBench-Vax: a benchmark framework for vaccine-ready neoantigen prioritization"),
        ("product", "An auditable triage layer for personalized cancer vaccine candidate selection"),
        ("product", "Decision support for neoantigen vaccine shortlisting under uncertainty"),
    ]
    write_tsv(
        OUT / "max_impact_title_bank.tsv",
        [{"track": t, "title": title} for t, title in title_rows],
        ["track", "title"],
    )

    write_text(
        OUT / "FIGURE_TABLE_BLUEPRINT.md",
        """
        # Figure and Table Blueprint

        ## Paper 1: Standard / Perspective
        - Figure 1: Vaccine-ready AI value chain from variants to manufactured vaccine.
        - Figure 2: Where benchmark leakage enters neoantigen AI.
        - Figure 3: NEO-PRIOR checklist as a clinical AI reporting standard.
        - Figure 4: Clean scientific validation versus product-assisted triage.
        - Table 1: NEO-PRIOR minimum reporting items.
        - Table 2: Failure modes and required controls.

        ## Paper 2: Benchmark / Implementation
        - Figure 1: NeoBench-Vax benchmark architecture.
        - Figure 2: Locked split contract and overlap audit.
        - Figure 3: v0/v1/v2 performance under AUPRC/top-k/source stress.
        - Figure 4: Source-heldout collapse and abstention boundary.
        - Figure 5: CROSS-Neo implementation and claim-boundary decision tree.
        - Table 1: Dataset provenance and split manifest.
        - Table 2: Model comparison under clean/product/source-stress tracks.
        - Table 3: Reviewer-facing failure-mode matrix.

        ## Product Demo
        - Figure 1: Candidate triage dashboard.
        - Figure 2: Evidence tag distribution.
        - Figure 3: Top-k enrichment versus comparator panel.
        - Table 1: Customer-safe candidate queue.
        - Table 2: NEO-PRIOR audit report.
        """,
    )

    write_text(
        OUT / "CLAIM_BOUNDARY_AND_DEFENSE.md",
        """
        # Claim Boundary and Defense

        ## Claims We Can Make
        - We propose NEO-PRIOR as a reporting and evaluation standard.
        - We define NeoBench-Vax as a benchmark infrastructure for vaccine-ready neoantigen prioritization.
        - CROSS-Neo is a first implementation evaluated under internal/locked split-safe tests.
        - AUPRC, top-k precision and enrichment are the primary decision metrics.
        - Source-heldout failure is a deployment warning and motivates abstention.

        ## Claims We Should Not Make
        - External clinical validation.
        - Universal superiority over all public predictors.
        - Quantum advantage.
        - Public predictors are clean baselines if their scores were used as features.
        - Negatives are definitive biological non-immunogenic cases unless experimentally established.

        ## How To Say It Safely
        Use:

        > internal locked retrospective evaluation

        Avoid:

        > externally validated

        Use:

        > product-assisted triage

        Avoid:

        > clean superiority claim

        Use:

        > fixed quantum-kernel fallback branch

        Avoid:

        > quantum advantage
        """,
    )

    write_text(
        OUT / "SENIOR_ADVISOR_EMAIL_SHORTLIST_TEMPLATE.md",
        """
        # Senior Advisor Email Template

        Subject: Draft standard for vaccine-ready neoantigen AI evaluation

        Dear [Name],

        We are preparing a Perspective proposing NEO-PRIOR, a reporting and evaluation standard for AI systems that rank neoantigen candidates for personalized cancer vaccine manufacture and immune testing.

        The motivation is that current computational papers often mix binding, presentation, immunogenicity and triage tasks, while public-data overlap, near-peptide retrieval, source/study shift, negative-label ambiguity and top-k decision utility are not reported consistently.

        We are also preparing a companion benchmark implementation, NeoBench-Vax/CROSS-Neo, as a concrete example of how the standard can be operationalized. We are not claiming external validation or universal model superiority; the goal is to make the evaluation contract clinically interpretable and auditable.

        We would value your feedback on the clinical framing and whether the checklist captures what vaccine teams would need before trusting an AI-generated shortlist.

        Best,
        Seungho
        """,
    )

    save_editor_one_screen()
    save_ecosystem()
    save_value_chain()
    save_risk_map()

    summary = {
        "package": "cross_neo_max_impact",
        "output_dir": str(OUT),
        "figures": [
            str(FIG / "fig_editor_one_screen_story.png"),
            str(FIG / "fig_neoprior_neobench_crossneo_ecosystem.png"),
            str(FIG / "fig_clinical_translation_value_chain.png"),
            str(FIG / "fig_claim_boundary_risk_map.png"),
        ],
        "claim_boundary": "internal/locked only; no external validation; no quantum advantage",
        "primary_story": "NEO-PRIOR standard + NeoBench-Vax benchmark + CROSS-Neo implementation",
    }
    (OUT / "impact_max_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
