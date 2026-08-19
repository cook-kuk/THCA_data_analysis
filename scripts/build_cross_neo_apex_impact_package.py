#!/usr/bin/env python3
"""Create the apex launch package for NEO-PRIOR / NeoBench-Vax / CROSS-Neo."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from textwrap import dedent


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
BASE = ROOT / "project/results/cross_neo_two_paper_strategy_2026_05_10"
OUT = BASE / "impact_apex"
FIG = OUT / "figures"


def text(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dedent(body).strip() + "\n", encoding="utf-8")


def tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def box(ax, x, y, w, h, title, body, fc, ec="#475569"):
    import matplotlib.patches as patches

    ax.add_patch(
        patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.014,rounding_size=0.02",
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.15,
        )
    )
    ax.text(x + 0.02, y + h - 0.045, title, fontsize=11.5, weight="bold", color="#0f172a", va="top")
    ax.text(x + 0.02, y + h - 0.105, body, fontsize=8.9, color="#334155", va="top", linespacing=1.18)


def make_figures() -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    FIG.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(15, 8.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.05, 0.96, "Apex Launch Story", fontsize=19, weight="bold", color="#0f172a", va="top")
    ax.text(0.05, 0.91, "Auditable AI for the cancer vaccine shortlist: standard + benchmark + first implementation.", fontsize=11, color="#475569")
    boxes = [
        (0.05, "Clinical need", "Scarce assay and\nmanufacturing slots.", "#fee2e2"),
        (0.285, "Evaluation gap", "Task mixing, overlap,\nsource bias, AUROC-only.", "#ffedd5"),
        (0.52, "NEO-PRIOR", "Minimum reporting\nand claim standard.", "#ecfeff"),
        (0.755, "NeoBench-Vax", "Public splits, audits,\nmodel cards, leaderboard.", "#f0fdf4"),
    ]
    for x, title, body, fc in boxes:
        box(ax, x, 0.66, 0.20, 0.18, title, body, fc)
    for x in [0.255, 0.49, 0.725]:
        ax.annotate("", xy=(x + 0.025, 0.75), xytext=(x, 0.75), arrowprops=dict(arrowstyle="->", lw=2, color="#64748b"))
    box(ax, 0.07, 0.38, 0.25, 0.16, "Paper 1", "Vaccine-ready AI for\npersonalized cancer\nimmunotherapy.", "#f8fafc")
    box(ax, 0.375, 0.38, 0.25, 0.16, "Paper 2", "Benchmarking vaccine-ready\nAI systems for neoantigen\nprioritization.", "#f8fafc")
    box(ax, 0.68, 0.38, 0.25, 0.16, "Product", "Auditable triage queue,\nevidence tags,\nabstention report.", "#f8fafc")
    box(ax, 0.18, 0.13, 0.64, 0.13, "Claim boundary", "Internal/locked retrospective only. No external validation. No quantum advantage. Clean comparator and product-assisted tracks remain separate.", "#fef9c3", "#ca8a04")
    fig.savefig(FIG / "fig_apex_one_slide_launch_story.png", dpi=260, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.axis("off")
    labels = ["datasets", "splits", "audits", "metrics", "models", "cards", "leaderboard", "adoption"]
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    for i, (lab, angle) in enumerate(zip(labels, angles)):
        x, y = 0.82 * np.cos(angle), 0.82 * np.sin(angle)
        ax.text(x, y, lab, ha="center", va="center", fontsize=12, weight="bold", color="#0f172a",
                bbox=dict(boxstyle="round,pad=0.45", fc="#f8fafc", ec="#475569"))
        x2, y2 = 0.82 * np.cos(angles[(i + 1) % len(labels)]), 0.82 * np.sin(angles[(i + 1) % len(labels)])
        ax.annotate("", xy=(x2 * 0.9, y2 * 0.9), xytext=(x * 0.9, y * 0.9),
                    arrowprops=dict(arrowstyle="->", color="#64748b", lw=1.4, connectionstyle="arc3,rad=0.2"))
    ax.text(0, 0.05, "NeoBench-Vax", ha="center", fontsize=21, weight="bold", color="#0f172a")
    ax.text(0, -0.10, "public benchmark flywheel", ha="center", fontsize=11, color="#475569")
    ax.text(0, 1.10, "Benchmark Infrastructure Flywheel", ha="center", fontsize=17, weight="bold", color="#0f172a")
    fig.savefig(FIG / "fig_apex_neobench_flywheel.png", dpi=260, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 5.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1)
    ax.axis("off")
    phases = [
        ("W1-2", "advisor asks\npresubmissions"),
        ("W3-4", "NEO-PRIOR v0.1\nmodel cards"),
        ("W5-6", "NeoBench-Vax\nleaderboard mockup"),
        ("W7-8", "Paper 1 submit\nproduct demo"),
        ("W9-12", "Paper 2 lock\npartner pilot"),
    ]
    ax.text(0.15, 0.92, "90-Day Apex Roadmap", fontsize=17, weight="bold", color="#0f172a", va="top")
    for i, (phase, body) in enumerate(phases):
        x = 0.35 + i * 1.85
        box(ax, x, 0.42, 1.45, 0.22, phase, body, "#ecfeff" if i % 2 == 0 else "#f0fdf4")
        if i < len(phases) - 1:
            ax.annotate("", xy=(x + 1.65, 0.53), xytext=(x + 1.47, 0.53), arrowprops=dict(arrowstyle="->", lw=1.7, color="#64748b"))
    ax.text(5, 0.18, "Decision: clinical advisor fast -> Nature Medicine/Cancer first; otherwise NBE/NMI benchmark route first.", ha="center", fontsize=11, weight="bold", color="#0f172a")
    fig.savefig(FIG / "fig_apex_90_day_roadmap.png", dpi=260, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.08, 0.94, "Claim Ladder", fontsize=17, weight="bold", color="#0f172a", va="top")
    ladder = [
        ("Now", "internal locked retrospective\nbenchmark + demo"),
        ("Next", "public benchmark release\nwith missing-data manifest"),
        ("Partner", "retrospective partner cohort\nunder locked protocol"),
        ("Clinical", "prospective assay/manufacturing\nshortlist study"),
        ("Not current", "clinical decision support\nregulatory claim"),
    ]
    for i, (title, body) in enumerate(ladder):
        box(ax, 0.14 + i * 0.035, 0.72 - i * 0.13, 0.70 - i * 0.07, 0.10, title, body, "#f8fafc")
    fig.savefig(FIG / "fig_apex_claim_ladder.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def landing_page() -> None:
    html = """
    <!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
    <title>NEO-PRIOR / NeoBench-Vax / CROSS-Neo</title>
    <style>
    body{margin:0;background:#09111f;color:#f8fafc;font-family:Inter,system-ui,sans-serif;line-height:1.55}
    header{min-height:88vh;display:grid;align-items:center;padding:72px 8vw;background:radial-gradient(circle at 30% 0,#12314f,#09111f 55%)}
    .k{color:#67e8f9;text-transform:uppercase;letter-spacing:.14em;font-size:12px;font-weight:800}
    h1{font-size:clamp(42px,6vw,84px);line-height:.96;max-width:980px;margin:18px 0}.lead{font-size:clamp(18px,2vw,24px);color:#a8b3c7;max-width:850px}
    main{padding:0 8vw 80px}.grid{display:grid;grid-template-columns:repeat(3,1fr);gap:16px}.card{background:#101b2d;border:1px solid rgba(148,163,184,.25);border-radius:8px;padding:22px}
    section{border-top:1px solid rgba(148,163,184,.25);padding:48px 0}h2{font-size:30px}.warn{color:#fde68a;font-weight:700}
    @media(max-width:820px){.grid{grid-template-columns:1fr}}
    </style></head><body><header><div><div class="k">NEO-PRIOR · NeoBench-Vax · CROSS-Neo</div>
    <h1>Auditable AI for the cancer vaccine shortlist.</h1>
    <p class="lead">A field-standard program for deciding whether neoantigen AI rankings are trustworthy enough to guide assay, manufacturing and clinical prioritization.</p>
    </div></header><main><section><h2>Three-layer contribution</h2><div class="grid">
    <div class="card"><h3>NEO-PRIOR</h3><p>Reporting standard: provenance, overlap, train-only pipelines, source/HLA stress, top-k, calibration, OOD and claim boundary.</p></div>
    <div class="card"><h3>NeoBench-Vax</h3><p>Public benchmark: split manifests, audits, model cards, public-score disclosure and clean/product-assisted tracks.</p></div>
    <div class="card"><h3>CROSS-Neo</h3><p>First implementation and stress-tested example. Failures become deployment boundaries, not hidden caveats.</p></div>
    </div></section><section><h2>Claim boundary</h2><p class="warn">Internal locked retrospective only. No external-validation claim. No quantum-advantage claim. Clean comparator and product-assisted tracks stay separate.</p></section></main></body></html>
    """
    text(OUT / "NEOBENCH_VAX_LANDING_PAGE_MOCKUP.html", html)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    text(
        OUT / "APEX_IMPACT_LAUNCH_MEMO.md",
        """
        # Apex Impact Launch Memo

        ## Final Form
        The high-impact story is a launch program, not a model report:

        1. **NEO-PRIOR**: reporting standard for vaccine-ready neoantigen AI.
        2. **NeoBench-Vax**: public benchmark and leaderboard infrastructure.
        3. **CROSS-Neo**: first implementation and stress-tested example.
        4. **Product-assisted triage**: business-facing audit layer with strict claim separation.

        ## Grand Claim
        Personalized cancer vaccines need auditable AI for the shortlist that decides what gets manufactured.

        ## Why This Raises Impact
        Standards create rules. Benchmarks create reusable infrastructure. CROSS-Neo becomes the first working example rather than the only reason the work matters.

        ## Immediate Moves
        Send presubmission inquiries, recruit a senior vaccine/immuno-oncology advisor, publish a benchmark mockup, and keep product-assisted demo claims separate from clean scientific validation.
        """,
    )
    text(
        OUT / "JOURNAL_SPECIFIC_ABSTRACTS.md",
        """
        # Journal-Specific Abstracts

        ## Nature Medicine / Nature Cancer
        Personalized cancer vaccines increasingly depend on computational prioritization of candidate neoantigens, yet the field lacks a shared standard for judging whether an AI-generated shortlist is trustworthy enough to guide immune testing, manufacturing or clinical prioritization. We propose NEO-PRIOR, a reporting and evaluation standard for vaccine-ready neoantigen AI, and outline NeoBench-Vax, a benchmark infrastructure with locked splits, overlap audits, model cards, top-k metrics, calibration and abstention reporting. The goal is not to claim that any current model is externally validated, but to define the evidence required before ranked shortlists can be interpreted clinically.

        ## Nature Biomedical Engineering / Nature Machine Intelligence
        Neoantigen prioritization is a small-n biomedical ranking problem where the output may determine assay allocation and vaccine manufacture. We introduce NeoBench-Vax, a benchmark infrastructure for vaccine-ready neoantigen AI, and CROSS-Neo as a first contamination-controlled implementation. The framework separates clean comparator and product-assisted tracks, requires exact/near/source/HLA/time holdouts, reports AUPRC and top-k precision as primary metrics, and publishes claim-boundary model cards.
        """,
    )
    text(
        OUT / "NEOBENCH_VAX_README_DRAFT.md",
        """
        # NeoBench-Vax

        NeoBench-Vax evaluates whether neoantigen AI systems can produce high-precision vaccine candidate shortlists while controlling public overlap, near-peptide memorization, source/study shift, HLA shortcuts, ambiguous negatives, calibration and OOD risk.

        ## Tracks
        - Clean Comparator: public predictor scores are not training features.
        - Product-Assisted: public predictor scores allowed if disclosed.
        - Source-Stress: source/study/prevalence shift.
        - OOD/Abstention: safe non-decision behavior.

        ## Primary Metrics
        AUPRC, top5/top10 precision, enrichment over prevalence, recall@k and abstention coverage versus precision.
        """,
    )
    text(
        OUT / "MODEL_CARD_TEMPLATE.md",
        """
        # NeoBench-Vax Model Card Template

        ## Model Name
        ## Intended Use
        ## Allowed Inputs
        ## Forbidden Inputs
        ## Training Data
        ## Split Contract
        ## Public Predictor Disclosure
        ## Primary Metrics
        ## Calibration and OOD
        ## Known Failure Modes
        ## Claim Boundary
        """,
    )
    text(
        OUT / "DATASET_CARD_TEMPLATE.md",
        """
        # NeoBench-Vax Dataset Card Template

        ## Dataset Name
        ## Source
        ## Cohort and Assay
        ## Labels and Negative-Label Ambiguity
        ## Candidate Generation
        ## HLA and Peptide Distribution
        ## Public Predictor Score Availability
        ## Overlap Audit
        ## Recommended Splits
        """,
    )
    text(
        OUT / "WEDNESDAY_DEMO_EXECUTIVE_SCRIPT.md",
        """
        # Wednesday Demo Executive Script

        We are not claiming a magic vaccine predictor. We are building the audit layer that tells a vaccine team which candidates to test first, why they were ranked, and when the model should not be trusted.

        The product output is a candidate queue with rank, score, HLA/peptide context, evidence tag, optional public comparator panel, OOD warning and NEO-PRIOR audit status.

        This demo is retrospective and internal. It is suitable for triage discussion and partner pilot design, not for clinical claims.
        """,
    )
    text(
        OUT / "PRESS_AND_FUNDING_PITCH.md",
        """
        # Press and Funding Pitch

        Personalized cancer vaccines are custom-built from a patient's tumor mutations. Before a vaccine can be made, scientists must choose a small number of candidate targets from a much larger list. That shortlist is increasingly generated by AI, but the field lacks a shared standard for deciding whether the AI shortlist is trustworthy.

        NEO-PRIOR is the proposed reporting standard. NeoBench-Vax is the benchmark framework. CROSS-Neo is the first implementation.
        """,
    )
    tsv(
        OUT / "apex_asset_inventory.tsv",
        [
            {"asset": "landing page mockup", "file": "NEOBENCH_VAX_LANDING_PAGE_MOCKUP.html", "audience": "editors/advisors/partners"},
            {"asset": "one-slide launch story", "file": "figures/fig_apex_one_slide_launch_story.png", "audience": "all"},
            {"asset": "benchmark flywheel", "file": "figures/fig_apex_neobench_flywheel.png", "audience": "editors/reviewers"},
            {"asset": "90-day roadmap", "file": "figures/fig_apex_90_day_roadmap.png", "audience": "internal/advisors"},
            {"asset": "claim ladder", "file": "figures/fig_apex_claim_ladder.png", "audience": "reviewers/partners"},
        ],
        ["asset", "file", "audience"],
    )
    landing_page()
    make_figures()
    summary = {
        "package": "cross_neo_apex_impact",
        "output_dir": str(OUT),
        "core_claim": "auditable AI for the cancer vaccine shortlist",
        "claim_boundary": "internal locked retrospective only; product-assisted demo separated from clean scientific validation",
    }
    (OUT / "apex_impact_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
