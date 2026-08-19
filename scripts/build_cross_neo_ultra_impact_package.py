#!/usr/bin/env python3
"""Build the ultra-impact editorial launch package for CROSS-Neo.

This package is intentionally different from the model/results folders. It is
an editor-facing strategy layer: journal fit, presubmission pitch, public
benchmark launch plan, consortium plan, and defensible claim boundaries.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from textwrap import dedent


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
BASE = ROOT / "project/results/cross_neo_two_paper_strategy_2026_05_10"
OUT = BASE / "impact_ultra"
FIG = OUT / "figures"


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


def box(ax, xy, wh, title, body, fc, ec="#334155", fs_title=12, fs_body=9.2):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    patch = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.015,rounding_size=0.020",
        linewidth=1.25,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(x + 0.025, y + h - 0.055, title, fontsize=fs_title, fontweight="bold", color="#0f172a", va="top")
    ax.text(x + 0.025, y + h - 0.125, body, fontsize=fs_body, color="#334155", va="top", linespacing=1.22)


def fig_submission_routes() -> None:
    import matplotlib.pyplot as plt

    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")
    ax.text(0.04, 0.96, "Ultra-Impact Submission Routes", fontsize=17, fontweight="bold", color="#0f172a", va="top")
    ax.text(0.04, 0.915, "Two papers plus public infrastructure. Claims stay internal/locked until prospective validation exists.", fontsize=10.5, color="#475569")
    box(ax, (0.05, 0.66), (0.27, 0.20), "Route A\nNature Medicine", "Perspective/Comment angle:\nvaccine-ready AI governance,\nclinical shortlist risk,\nNEO-PRIOR standard.", "#ecfeff", "#0891b2")
    box(ax, (0.365, 0.66), (0.27, 0.20), "Route B\nNature Cancer", "Cancer immunotherapy angle:\nneoantigen vaccine triage,\nsource bias, auditability,\nclinical translation boundary.", "#fdf2f8", "#db2777")
    box(ax, (0.68, 0.66), (0.27, 0.20), "Route C\nNature Biomedical Engineering", "Article/Perspective angle:\ncomputational health system,\nbenchmark infrastructure,\nvalidation and deployment.", "#f0fdf4", "#16a34a")
    box(ax, (0.20, 0.35), (0.27, 0.19), "Route D\nNature Machine Intelligence", "Analysis/Comment angle:\nAI reliability, benchmark design,\nuncertainty, public-score leakage,\nnegative-label ambiguity.", "#eef2ff", "#4f46e5")
    box(ax, (0.55, 0.35), (0.27, 0.19), "Route E\nCancer Cell / Cancer Discovery", "If senior clinical coauthor joins:\nstrong cancer-vaccine narrative,\nreview/opinion plus benchmark roadmap.", "#fff7ed", "#ea580c")
    ax.annotate("", xy=(0.50, 0.61), xytext=(0.50, 0.56), arrowprops=dict(arrowstyle="-|>", lw=2, color="#64748b"))
    box(ax, (0.28, 0.08), (0.44, 0.16), "Core Asset That Travels Across All Routes", "NEO-PRIOR + NeoBench-Vax + CROSS-Neo\nstandard, public benchmark, first implementation.\nNo external-validation or quantum-advantage claim.", "#f8fafc", "#475569", fs_title=12.5, fs_body=10)
    fig.savefig(FIG / "fig_ultra_submission_routes.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def fig_public_infrastructure() -> None:
    import matplotlib.pyplot as plt

    modules = [
        ("Registry", "datasets\nIDs\nprovenance"),
        ("Splits", "exact\nnear\nsource\nHLA\ntime"),
        ("Score Matrix", "public predictors\nclean track\nproduct track"),
        ("Audits", "overlap\nretrieval\nclaim boundary"),
        ("Model Cards", "inputs\nforbidden uses\nfailure modes"),
        ("Leaderboard", "AUPRC\ntop-k\nabstention\ncollapse penalty"),
    ]
    fig, ax = plt.subplots(figsize=(14, 6.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.04, 0.94, "NeoBench-Vax as Public Infrastructure", fontsize=17, fontweight="bold", color="#0f172a", va="top")
    ax.text(0.04, 0.89, "This is the asset that converts a model project into a reusable field resource.", fontsize=10.5, color="#475569")
    x0 = 0.05
    for i, (title, body) in enumerate(modules):
        x = x0 + i * 0.155
        box(ax, (x, 0.42), (0.13, 0.27), title, body, "#f8fafc", "#334155", fs_title=10.5, fs_body=8.6)
        if i < len(modules) - 1:
            ax.annotate("", xy=(x + 0.15, 0.555), xytext=(x + 0.132, 0.555), arrowprops=dict(arrowstyle="->", lw=1.7, color="#64748b"))
    box(ax, (0.20, 0.13), (0.60, 0.16), "Launch Promise", "Every submitted model is evaluated under the same provenance, overlap, split, metric and claim-boundary contract.", "#ecfeff", "#0891b2", fs_title=12, fs_body=10)
    fig.savefig(FIG / "fig_ultra_neobench_public_infrastructure.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def fig_checklist_radar() -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    labels = ["Task", "Provenance", "Overlap", "Splits", "Top-k", "Calibration", "OOD", "Negatives", "Source", "Claim"]
    values = np.array([0.95, 0.85, 0.95, 0.90, 0.90, 0.78, 0.84, 0.82, 0.88, 0.95])
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    values_closed = np.r_[values, values[0]]
    angles_closed = np.r_[angles, angles[0]]
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, polar=True)
    ax.plot(angles_closed, values_closed, color="#0891b2", lw=2.5)
    ax.fill(angles_closed, values_closed, color="#67e8f9", alpha=0.30)
    ax.set_xticks(angles)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_yticklabels([])
    ax.set_ylim(0, 1)
    ax.set_title("NEO-PRIOR Coverage Map", fontsize=16, fontweight="bold", pad=25)
    fig.savefig(FIG / "fig_ultra_neoprior_coverage_map.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def fig_evidence_stack() -> None:
    import matplotlib.pyplot as plt

    layers = [
        ("Clinical Need", "shortlists decide assay/manufacturing slots", "#fee2e2", "#dc2626"),
        ("Evaluation Gap", "overlap, source bias, negatives, top-k instability", "#ffedd5", "#ea580c"),
        ("NEO-PRIOR", "minimum reporting and evaluation standard", "#ecfeff", "#0891b2"),
        ("NeoBench-Vax", "public benchmark, score matrix and model cards", "#f0fdf4", "#16a34a"),
        ("CROSS-Neo", "first contamination-controlled implementation", "#eef2ff", "#4f46e5"),
        ("Product", "auditable triage and abstention report", "#fdf2f8", "#db2777"),
    ]
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.06, 0.95, "Evidence Stack for a High-Impact Story", fontsize=17, fontweight="bold", color="#0f172a", va="top")
    for i, (title, body, fc, ec) in enumerate(layers):
        y = 0.78 - i * 0.12
        box(ax, (0.12 + i * 0.015, y), (0.74 - i * 0.03, 0.085), title, body, fc, ec, fs_title=11.5, fs_body=9)
    ax.text(0.5, 0.06, "Impact rises because each layer is reusable beyond one model result.", ha="center", fontsize=11.5, fontweight="bold", color="#0f172a")
    fig.savefig(FIG / "fig_ultra_evidence_stack.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    write_text(
        OUT / "ULTRA_IMPACT_EDITORIAL_DOSSIER.md",
        """
        # CROSS-Neo Ultra-Impact Editorial Dossier

        ## Executive Upgrade
        The program should be pitched as a three-part field intervention:

        1. **NEO-PRIOR**: a reporting and evaluation standard for vaccine-ready neoantigen AI.
        2. **NeoBench-Vax**: public benchmark infrastructure with split manifests, overlap audits, model cards and leaderboard tracks.
        3. **CROSS-Neo**: the first implementation and stress-tested example under the standard.

        ## Strongest Editorial Thesis
        Personalized cancer vaccines now depend on computational shortlists, but the field lacks a consistent way to judge whether a shortlist is trustworthy.

        ## Why This Is Bigger Than a Model
        A model paper asks whether one score is better. This package asks what evidence any model must provide before its shortlist can influence assay, manufacturing or clinical prioritization.

        ## Paper 1
        **Title:** Vaccine-ready AI for personalized cancer immunotherapy

        **Contribution:** NEO-PRIOR standard plus the clinical and benchmark rationale.

        **Best target:** Nature Medicine or Nature Cancer Perspective/Comment route.

        **Core display items:**
        - clinical value chain;
        - failure modes in current evaluation;
        - NEO-PRIOR checklist;
        - clean validation versus product-assisted triage.

        ## Paper 2
        **Title:** Benchmarking vaccine-ready AI systems for neoantigen prioritization

        **Contribution:** NeoBench-Vax benchmark and CROSS-Neo implementation.

        **Best target:** Nature Biomedical Engineering Article or Perspective route; Nature Machine Intelligence as an AI-benchmark alternative.

        **Core display items:**
        - benchmark architecture;
        - locked split contract;
        - model comparison with AUPRC/top-k/source stress;
        - claim-boundary and abstention decision tree.

        ## What Must Stay Conservative
        - No external validation claim.
        - No quantum advantage claim.
        - No universal superiority claim over public predictors.
        - Public predictor scores must be separated into clean-comparator versus product-assisted tracks.
        - Source-heldout collapse is a finding and boundary, not a detail to hide.

        ## One-Sentence Pitch
        Personalized cancer vaccines need auditable AI for the shortlist that decides what gets manufactured.
        """,
    )

    journal_rows = [
        {
            "journal": "Nature Medicine",
            "best_use": "Paper 1 Perspective or Comment style pitch",
            "official_fit_basis": "Translational and clinical research with originality, timeliness and impact on human health; publishes Reviews, Perspectives and other contextual content.",
            "why_it_fits": "The story is a clinical AI governance standard for personalized cancer vaccine shortlists.",
            "main_risk": "Without senior clinical vaccine voice, may be seen as too methodological.",
            "source_url": "https://www.nature.com/nm/aims",
        },
        {
            "journal": "Nature Cancer",
            "best_use": "Paper 1 cancer-focused Perspective or Comment",
            "official_fit_basis": "Covers major cancer advances across life, physical, applied and social sciences, including immunotherapy, therapy delivery, clinical work and methodology.",
            "why_it_fits": "The framing is cancer immunotherapy and neoantigen vaccine candidate selection.",
            "main_risk": "Needs strong cancer-vaccine clinical framing, not only ML benchmarking.",
            "source_url": "https://www.nature.com/natcancer/aims",
        },
        {
            "journal": "Nature Biomedical Engineering",
            "best_use": "Paper 2 Article; alternative Perspective",
            "official_fit_basis": "Covers experimental and computational systems, methods, technologies and processes that may improve health or healthcare; includes AI, computational medicine, precision medicine and vaccines.",
            "why_it_fits": "NeoBench-Vax plus CROSS-Neo is a computational health technology and validation/deployment infrastructure.",
            "main_risk": "Article route needs robust system, clear validation package and reproducibility artifacts.",
            "source_url": "https://www.nature.com/natbiomedeng/aims",
        },
        {
            "journal": "Nature Machine Intelligence",
            "best_use": "Paper 2 Analysis/Comment alternative",
            "official_fit_basis": "Publishes AI/ML research and discusses impacts across science, healthcare, society and industry.",
            "why_it_fits": "Benchmark design, uncertainty, leakage, negative-label ambiguity and claim boundaries are AI evaluation issues.",
            "main_risk": "Needs broader ML contribution beyond domain-specific benchmark.",
            "source_url": "https://www.nature.com/natmachintell/aims",
        },
        {
            "journal": "Nature Biomedical Engineering Perspective",
            "best_use": "If Article data package remains too small",
            "official_fit_basis": "Perspectives are forward-looking, balanced, up to about 2500 words and up to 4 display items per current content-type page.",
            "why_it_fits": "A compact standard plus benchmark proposal can fit better than an overclaimed original Article.",
            "main_risk": "Perspective should not be dominated by unpublished own results.",
            "source_url": "https://www.nature.com/natbiomedeng/content",
        },
    ]
    write_tsv(
        OUT / "JOURNAL_FIT_OFFICIAL_SCOPE_MATRIX.tsv",
        journal_rows,
        ["journal", "best_use", "official_fit_basis", "why_it_fits", "main_risk", "source_url"],
    )

    write_text(
        OUT / "PRESUBMISSION_INQUIRY_PACK.md",
        """
        # Presubmission Inquiry Pack

        ## Universal Subject Line
        Presubmission inquiry: reporting standard for vaccine-ready neoantigen AI

        ## Nature Medicine / Nature Cancer Version
        Dear Editors,

        We are preparing a Perspective proposing NEO-PRIOR, a reporting and evaluation standard for AI systems that rank neoantigen candidates for personalized cancer vaccine manufacture, immune testing and clinical prioritization.

        The motivation is that candidate selection for personalized vaccines now depends on computational shortlists, yet published evaluations often mix binding, presentation, immunogenicity and triage tasks, while public-data overlap, near-peptide memorization, source/study shift, negative-label ambiguity, calibration and top-k utility are not reported consistently.

        We propose a minimum reporting checklist and a clean separation between scientific validation and product-assisted triage. A companion benchmark framework, NeoBench-Vax, illustrates how models can be evaluated using locked splits, overlap audits, model cards and abstention-aware metrics. We do not claim external validation or universal model superiority; instead, we argue that the field needs an auditable evaluation contract before candidate-ranking AI can be interpreted clinically.

        We would be grateful for your view on whether this topic would be suitable as a Perspective or Comment.

        ## Nature Biomedical Engineering Version
        Dear Editors,

        We are preparing a manuscript on NeoBench-Vax, a benchmark and validation infrastructure for vaccine-ready neoantigen prioritization systems, together with CROSS-Neo as a first contamination-controlled implementation.

        The work treats neoantigen prioritization as a biomedical decision-support system rather than a single-score predictor. The benchmark includes locked exact, near-peptide, source/study, HLA and retrieval-clean splits; public-overlap audits; model-card requirements; AUPRC/top-k/enrichment as primary metrics; calibration and OOD/abstention reporting; and a clean separation between public predictors as comparators versus product-assisted features.

        Our claim boundary is conservative: internal locked retrospective evaluation only, no external validation claim, and no quantum-advantage claim. The central contribution is a reproducible evaluation framework for computational systems that may influence assay and manufacturing shortlists in personalized cancer vaccine workflows.

        We would welcome guidance on whether this is better suited as an Article or Perspective.

        ## Nature Machine Intelligence Version
        Dear Editors,

        We are preparing an analysis of AI evaluation failure modes in neoantigen prioritization for personalized cancer vaccines. The manuscript focuses on leakage-sensitive benchmarks, public predictor reuse, ambiguous negative labels, source/study shift, uncertainty, abstention and top-k decision metrics.

        We introduce NEO-PRIOR as a reporting checklist and NeoBench-Vax as a benchmark framework, with CROSS-Neo as a first implementation. The broader AI contribution is an evaluation contract for small-n biomedical ranking problems where the practical objective is a short, high-precision candidate list rather than a calibrated population-level classifier.

        We would appreciate your advice on fit for an Analysis, Comment or related format.
        """,
    )

    write_text(
        OUT / "NEO_PRIOR_TWO_PAGE_CHECKLIST.md",
        """
        # NEO-PRIOR Two-Page Checklist

        ## Page 1: Minimum Reporting Items
        1. Define the prediction task: binding, presentation, immunogenicity, triage or benchmark ranking.
        2. Report data provenance: study, patient, assay, HLA, candidate-generation route and date/year when available.
        3. Disclose public predictor use as feature, comparator, filter or product assist.
        4. Report exact peptide and exact peptide-HLA overlap.
        5. Report near-peptide, source-window, mutation-pair and TCR motif overlap when available.
        6. Keep retrieval, scaling, imputation, calibration, feature selection and thresholding inside train folds.
        7. Use exact, near, source/study, HLA allele/supertype, time/assay and retrieval-clean splits where feasible.
        8. Treat negatives as ambiguous unless experimentally established as true non-immunogenic cases.
        9. Report AUPRC, top-k precision, enrichment over prevalence and recall@k as primary metrics.
        10. Report AUROC as secondary, not as the primary evidence of shortlist utility.

        ## Page 2: Deployment and Claim Boundary
        11. Report Brier/ECE, calibration curves and score reliability.
        12. Report OOD flags and abstention coverage versus precision.
        13. Report per-source and source-heldout results.
        14. Report per-HLA and HLA-heldout/supertype-heldout results.
        15. Publish split manifests, row IDs, seeds and scripts where possible.
        16. Publish model cards with allowed inputs, forbidden inputs and known failures.
        17. Separate clean scientific validation from product-assisted triage.
        18. State whether public predictor scores are forbidden, allowed or disclosed.
        19. State whether the evidence is internal, locked retrospective, external retrospective or prospective.
        20. Include a one-paragraph claim boundary in the abstract or discussion.
        """,
    )

    write_text(
        OUT / "NEOBENCH_VAX_PUBLIC_LAUNCH_SPEC.md",
        """
        # NeoBench-Vax Public Launch Spec

        ## Minimum Public Website Sections
        1. Benchmark overview.
        2. Dataset registry and provenance cards.
        3. Locked split downloads.
        4. Public overlap audit reports.
        5. Clean comparator leaderboard.
        6. Product-assisted leaderboard.
        7. Source-stress and OOD/abstention leaderboard.
        8. Model-card template.
        9. Submission instructions.
        10. Claim-boundary policy.

        ## Leaderboard Columns
        - model name;
        - allowed inputs;
        - public predictor score use;
        - split;
        - n;
        - prevalence;
        - AUPRC;
        - AUPRC / prevalence;
        - top5 precision;
        - top10 precision;
        - enrichment@10;
        - recall@10;
        - Brier;
        - ECE;
        - abstention coverage;
        - source-heldout collapse penalty;
        - claim eligibility.

        ## Governance Rule
        A model may be useful in the product-assisted track while being ineligible for clean scientific superiority claims. This rule is the core trust mechanism.

        ## First Public Release
        Release as a static website and GitHub-style benchmark repository before or alongside Paper 2 submission. The first version can include internal locked examples, missing-data manifests and public-overlap status, provided the claim boundary is explicit.
        """,
    )

    write_text(
        OUT / "CONSORTIUM_AND_ADVISORY_BOARD_PLAN.md",
        """
        # Consortium and Advisory Board Plan

        ## Goal
        Convert NEO-PRIOR from an author proposal into a community-facing standard.

        ## Minimum Advisor Mix
        - one personalized cancer vaccine trialist;
        - one tumor immunologist with T cell/HLA expertise;
        - one computational immunology benchmark expert;
        - one clinical AI or medical device evaluation expert;
        - one biotech/product translation advisor.

        ## What To Ask For
        1. Does the checklist include the evidence a vaccine team would need before trusting a shortlist?
        2. Which reporting items are mandatory versus recommended?
        3. Which public datasets/comparators must be included?
        4. Which failure modes should trigger abstention rather than ranking?
        5. What prospective validation would be required after this retrospective framework?

        ## Why This Raises Impact
        A standard with expert input can be pitched as a field resource. A model-only manuscript cannot.
        """,
    )

    risk_rows = [
        {"risk": "Seen as only a small-n model paper", "impact": "high", "mitigation": "Lead with NEO-PRIOR standard and NeoBench-Vax infrastructure; keep CROSS-Neo as first implementation."},
        {"risk": "External validation overclaim", "impact": "high", "mitigation": "Use internal/locked retrospective wording everywhere."},
        {"risk": "Quantum advantage skepticism", "impact": "medium", "mitigation": "Do not claim quantum advantage; call it fixed QK fallback branch."},
        {"risk": "Public predictor leakage concern", "impact": "high", "mitigation": "Separate clean comparator and product-assisted tracks."},
        {"risk": "Source-heldout collapse", "impact": "high", "mitigation": "Make it a benchmark failure mode and abstention requirement."},
        {"risk": "Nature Medicine demands stronger clinical voice", "impact": "high", "mitigation": "Recruit senior vaccine/immuno-oncology advisor before presubmission."},
        {"risk": "NBE asks for more system validation", "impact": "medium", "mitigation": "Emphasize reproducible infrastructure, model cards, locked split registry and demo pipeline."},
        {"risk": "NMI wants broader AI novelty", "impact": "medium", "mitigation": "Frame as small-n biomedical ranking evaluation under ambiguous negatives and leakage-sensitive public benchmarks."},
    ]
    write_tsv(OUT / "ULTRA_IMPACT_RISK_REGISTER.tsv", risk_rows, ["risk", "impact", "mitigation"])

    write_text(
        OUT / "HIGH_IMPACT_SUBMISSION_DECISION_TREE.md",
        """
        # High-Impact Submission Decision Tree

        ## If a senior cancer-vaccine coauthor joins quickly
        First submit Paper 1 as a Nature Medicine or Nature Cancer Perspective/Comment presubmission.

        ## If the clinical coauthor is delayed
        Submit the standard as Nature Biomedical Engineering Perspective or Nature Machine Intelligence Analysis while continuing advisor recruitment.

        ## If Paper 2 data package improves
        Aim Nature Biomedical Engineering Article with NeoBench-Vax and CROSS-Neo as the system.

        ## If source-heldout remains weak
        Do not hide it. Position it as evidence for source-stress reporting and abstention.

        ## If public-overlap audit remains incomplete
        Do not claim clean superiority against public tools. Use missing-data manifest and product-assisted labels.

        ## If demo pressure is immediate
        Use the product-assisted track only. Keep scientific claims separate.
        """,
    )

    write_text(
        OUT / "ONE_PAGE_GRAND_CHALLENGE_MANIFESTO.md",
        """
        # Grand Challenge Manifesto

        Personalized cancer vaccines require a short list of candidates. That list is increasingly shaped by computational models, but the field lacks a shared standard for deciding whether a ranked list is trustworthy.

        The central issue is not whether one more model can improve AUROC. The central issue is whether any model can show that its shortlist is not driven by public-data overlap, near-peptide memorization, HLA shortcuts, source-specific assay bias, or ambiguous negative labels.

        NEO-PRIOR defines the minimum evidence that vaccine-ready neoantigen AI should report. NeoBench-Vax makes the standard operational through locked splits, public-overlap audits, model cards and leaderboard tracks. CROSS-Neo is the first implementation, evaluated conservatively under this framework.

        The clinical promise is not a magic oracle. The practical promise is auditable prioritization: fewer candidates, clearer evidence tags, better top-k enrichment, and explicit abstention when the evidence is unsafe.
        """,
    )

    fig_submission_routes()
    fig_public_infrastructure()
    fig_checklist_radar()
    fig_evidence_stack()

    summary = {
        "package": "cross_neo_ultra_impact",
        "output_dir": str(OUT),
        "core_pitch": "Personalized cancer vaccines need auditable AI for the shortlist that decides what gets manufactured.",
        "new_assets": [
            "official journal-fit matrix",
            "presubmission inquiry pack",
            "NEO-PRIOR two-page checklist",
            "NeoBench-Vax public launch spec",
            "consortium/advisory board plan",
            "ultra-impact risk register",
            "submission decision tree",
            "four editorial figures",
        ],
        "claim_boundary": "internal/locked retrospective only unless future prospective evidence is added",
    }
    (OUT / "ultra_impact_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
