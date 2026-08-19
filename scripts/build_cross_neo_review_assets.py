#!/usr/bin/env python3
"""Build review paper figure/table/data provenance assets for CROSS-Neo.

This creates scaffolded review figures and precise data-use maps so the review
can be written as a figure-led, evidence-traceable manuscript.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "project/results/cross_neo_review_2026_05_10"
FIG = OUT / "figures"
PRODUCT = REPO / "project/results/cross_neo_product_demo_2026_05_10"
V0 = REPO / "project/results/cross_neo_v0"
V1 = REPO / "project/results/cross_neo_v1"
V2 = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
NEO = REPO / "project/results/p_neo_bayesian_2026_05_09"


PALETTE = {
    "clinical": "#2f6f9f",
    "model": "#4b8f73",
    "risk": "#b36b2c",
    "contract": "#5b5b8a",
    "grey": "#6b7280",
    "light": "#f4f7fa",
    "dark": "#17202a",
}


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def wrap(s: str, width: int = 22) -> str:
    return "\n".join(textwrap.wrap(s, width=width))


def add_box(ax, xy, w, h, text, color, fontsize=9, text_color="white", radius=0.06):
    box = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle=f"round,pad=0.015,rounding_size={radius}",
        facecolor=color,
        edgecolor="none",
        alpha=0.95,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=text_color,
        weight="bold",
    )
    return box


def add_arrow(ax, p1, p2, color="#333333", lw=1.8):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=13, linewidth=lw, color=color))


def fig1_funnel() -> Path:
    path = FIG / "fig1_vaccine_ready_prioritization_funnel.png"
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.03, 0.94, "Figure 1. Vaccine-ready neoantigen prioritization funnel", fontsize=16, weight="bold", color=PALETTE["dark"])
    ax.text(0.03, 0.89, "The bottleneck is no longer only finding HLA binders; it is ranking the few candidates worth manufacturing.", fontsize=10.5, color=PALETTE["grey"])

    steps = [
        ("Tumor/normal\nsequencing", "10^3-10^5\nvariants"),
        ("Candidate\npeptides", "10^3-10^4"),
        ("HLA presentation\nfilter", "10^2-10^3"),
        ("Immunogenicity\nranking", "10-100"),
        ("Vaccine/assay\nshortlist", "5-20"),
    ]
    xs = np.linspace(0.07, 0.83, len(steps))
    widths = [0.16, 0.145, 0.145, 0.145, 0.145]
    heights = [0.62, 0.50, 0.39, 0.29, 0.20]
    for i, ((label, n), x, w, h) in enumerate(zip(steps, xs, widths, heights)):
        y = 0.22 + (0.62 - h) / 2
        color = PALETTE["clinical"] if i < 3 else PALETTE["model"] if i == 3 else PALETTE["contract"]
        add_box(ax, (x, y), w, h, wrap(label, 16), color, fontsize=10)
        ax.text(x + w / 2, y - 0.045, n, ha="center", va="center", fontsize=10, color=PALETTE["dark"], weight="bold")
        if i < len(steps) - 1:
            add_arrow(ax, (x + w + 0.008, 0.53), (xs[i + 1] - 0.012, 0.53), color=PALETTE["grey"])

    add_box(ax, (0.78, 0.76), 0.18, 0.10, "Decision metric:\ntop-k utility", PALETTE["risk"], fontsize=10)
    ax.text(0.05, 0.10, "Review message: vaccine pipelines need prioritizers evaluated where decisions happen: the top 5-20 candidates.", fontsize=10.5, color=PALETTE["dark"])
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def fig2_taxonomy() -> Path:
    path = FIG / "fig2_predictor_taxonomy_by_biological_layer.png"
    fig, ax = plt.subplots(figsize=(13, 7.2))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.03, 0.95, "Figure 2. Predictor families mapped to biological decision layers", fontsize=16, weight="bold", color=PALETTE["dark"])
    layers = [
        ("Binding /\npresentation", "NetMHCpan\nMHCflurry\nMixMHCpred\nBigMHC-EL", PALETTE["clinical"]),
        ("Immunogenicity\npropensity", "PRIME\nDeepImmuno\nBigMHC-IM\nTransPHLA\nT-SCAPE", PALETTE["model"]),
        ("Patient / tumor\ncontext", "IMPROVE\nexpression\nclonality\nRNA support", PALETTE["contract"]),
        ("Structure / TCR\nrecognition", "pMHC geometry\nTCR-facing exposure\nTCR motif evidence", PALETTE["risk"]),
        ("Vaccine-ready\ntriage", "top-k rank\nenrichment\ncalibration\nabstention", "#3f7f7f"),
    ]
    xs = [0.05, 0.245, 0.44, 0.635, 0.83]
    for i, (title, examples, color) in enumerate(layers):
        add_box(ax, (xs[i], 0.54), 0.15, 0.22, title, color, fontsize=10)
        ax.text(xs[i] + 0.075, 0.40, examples, ha="center", va="center", fontsize=9, color=PALETTE["dark"])
        if i < len(layers) - 1:
            add_arrow(ax, (xs[i] + 0.155, 0.65), (xs[i + 1] - 0.01, 0.65), color=PALETTE["grey"])
    ax.text(0.06, 0.19, "Key point", fontsize=12, weight="bold", color=PALETTE["dark"])
    ax.text(0.06, 0.13, "A strong binder is not automatically a strong vaccine target. The review separates biological layers and then asks which metrics match the final decision.", fontsize=10.5, color=PALETTE["dark"])
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def fig3_failure_modes() -> Path:
    path = FIG / "fig3_benchmark_failure_modes_and_controls.png"
    fig, ax = plt.subplots(figsize=(13, 7.5))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.03, 0.95, "Figure 3. Benchmark failure modes that can mimic model progress", fontsize=16, weight="bold", color=PALETTE["dark"])
    failures = [
        ("Exact peptide-HLA\noverlap", "control:\nexact holdout"),
        ("Near-peptide /\nsource-protein overlap", "control:\ncluster + window holdout"),
        ("HLA allele /\nsupertype shortcut", "control:\nHLA/supertype holdout"),
        ("Source, study,\nassay shift", "control:\nstudy/source holdout"),
        ("Ambiguous\nnegative labels", "control:\nPU-aware reporting"),
        ("AUROC-only\nleaderboards", "control:\nAUPRC + top-k"),
    ]
    for i, (risk, control) in enumerate(failures):
        row = 0 if i < 3 else 1
        col = i % 3
        x = 0.07 + col * 0.30
        y = 0.62 - row * 0.32
        add_box(ax, (x, y), 0.22, 0.12, wrap(risk, 18), PALETTE["risk"], fontsize=10)
        add_arrow(ax, (x + 0.11, y - 0.01), (x + 0.11, y - 0.09), color=PALETTE["grey"])
        add_box(ax, (x, y - 0.22), 0.22, 0.10, wrap(control, 18), PALETTE["contract"], fontsize=9)
    ax.text(0.06, 0.07, "Review message: clean neoantigen benchmarking requires split families that remove both exact and near-reference information.", fontsize=10.5, color=PALETTE["dark"])
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def fig4_metrics_contract() -> Path:
    path = FIG / "fig4_metric_contract_and_competitor_gauntlet.png"
    gauntlet = PRODUCT / "strong_competitors/strong_competitor_gauntlet.tsv"
    fig, ax = plt.subplots(figsize=(13.5, 7.6))
    if gauntlet.exists():
        df = pd.read_csv(gauntlet, sep="\t")
        strict = df[df["benchmark_scope"] == "strict_product_demo_set_n89"].copy()
        strict = strict.sort_values(["top10_precision", "AUPRC"], ascending=False).head(10)
        labels = strict["method"].str.replace("_", " ", regex=False)
        y = np.arange(len(strict))
        colors = np.where(strict["method_family"].str.contains("CROSS", na=False), PALETTE["clinical"], PALETTE["risk"])
        ax.barh(y, strict["top10_precision"], color=colors, alpha=0.88, label="top10 precision")
        ax.scatter(strict["AUPRC"], y, color=PALETTE["dark"], s=42, label="AUPRC")
        ax.axvline(strict["prevalence"].iloc[0], color=PALETTE["grey"], linestyle="--", linewidth=1.2, label="prevalence")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlim(0, 1.02)
        ax.set_xlabel("Metric value")
        ax.set_title("Figure 4. Internal strict demo: top-k and AUPRC reveal practical prioritization", fontsize=15, weight="bold", loc="left")
        ax.legend(frameon=False, loc="lower right")
        ax.grid(axis="x", alpha=0.15)
    else:
        ax.set_axis_off()
        ax.text(0.1, 0.5, "Gauntlet file missing", fontsize=16)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def fig5_launch_flow() -> Path:
    path = FIG / "fig5_review_to_original_launch_flow.png"
    fig, ax = plt.subplots(figsize=(13, 6.8))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.03, 0.94, "Figure 5. Two-paper strategy: define the standard, then implement it", fontsize=16, weight="bold", color=PALETTE["dark"])
    add_box(ax, (0.08, 0.55), 0.32, 0.22, "Review / Perspective\nField standard:\nvaccine-ready prioritization", PALETTE["clinical"], fontsize=12)
    add_arrow(ax, (0.42, 0.66), (0.58, 0.66), color=PALETTE["grey"], lw=2.2)
    add_box(ax, (0.60, 0.55), 0.32, 0.22, "Original CROSS-Neo\nImplementation:\ncontamination-controlled\npan-allele triage", PALETTE["model"], fontsize=12)
    ax.text(0.10, 0.38, "Outputs:\n- evaluation contract\n- figure-led review\n- benchmark caveats\n- metric standards", fontsize=10.5, color=PALETTE["dark"])
    ax.text(0.62, 0.38, "Outputs:\n- locked split metrics\n- top-k candidate queue\n- OOD/abstention\n- product/demo separation", fontsize=10.5, color=PALETTE["dark"])
    ax.text(0.08, 0.16, "Category position: not another MHC-binding predictor; a framework for vaccine candidate prioritization.", fontsize=12, weight="bold", color=PALETTE["contract"])
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def write_data_use_map() -> None:
    rows = [
        {
            "data_or_source": "Recent personalized neoantigen vaccine clinical studies",
            "specific_examples": "mRNA-4157/V940 melanoma; autogene cevumeran pancreatic cancer; RCC personalized vaccine",
            "role_in_review": "clinical motivation that vaccine target ranking matters now",
            "used_in": "Fig 1; Table 1; Sections 1 and Why Now",
            "local_path": "review_reference_evidence_matrix.tsv",
            "claim_boundary": "clinical motivation, not predictor benchmark",
        },
        {
            "data_or_source": "Public predictor literature",
            "specific_examples": "NetMHCpan, MHCflurry, MixMHCpred, PRIME, BigMHC, DeepImmuno, TransPHLA, T-SCAPE, IMPROVE",
            "role_in_review": "taxonomy of biological layers and current computational landscape",
            "used_in": "Fig 2; Table 2; Section 3",
            "local_path": "review_reference_evidence_matrix.tsv; strong_competitor_gauntlet.tsv",
            "claim_boundary": "public predictors are comparators or literature context, not clean CROSS-Neo features",
        },
        {
            "data_or_source": "CROSS-Neo v0/v1 leakage and failure audits",
            "specific_examples": "retrieval leakage audit, branch errors, source shift, public overlap audit",
            "role_in_review": "internal examples motivating benchmark failure modes",
            "used_in": "Fig 3; Table 3; benchmarking failure modes section",
            "local_path": f"{V0}; {V1}",
            "claim_boundary": "internal audit examples only; do not claim external validation",
        },
        {
            "data_or_source": "Strong competitor gauntlet",
            "specific_examples": "MHCflurry, BigMHC, DeepImmuno, PRIME, TransPHLA, NetMHCpan, T-SCAPE against CROSS-Neo product components",
            "role_in_review": "example of why top-k/AUPRC can differ from AUROC-only public leaderboards",
            "used_in": "Fig 4; Table 4; optional Box 1",
            "local_path": f"{PRODUCT}/strong_competitors/strong_competitor_gauntlet.tsv",
            "claim_boundary": "internal retrospective demo benchmark; not external validation",
        },
        {
            "data_or_source": "CROSS-Neo product demo package",
            "specific_examples": "customer-safe priority queue; selected train/reference sources; product fixed pan-allele score",
            "role_in_review": "not central review evidence; supports original-paper bridge and product story",
            "used_in": "Fig 5; bridge section; original paper plan",
            "local_path": f"{PRODUCT}",
            "claim_boundary": "business/product artifact; separate from clean manuscript claim",
        },
        {
            "data_or_source": "CROSS-Neo v2 SOTA sprint",
            "specific_examples": "selective ensemble; ESM2/QK gate; source stress and overlap audit",
            "role_in_review": "optional internal implementation evidence for original paper, not review centerpiece",
            "used_in": "Original CROSS-Neo paper plan; not main review evidence",
            "local_path": f"{V2}",
            "claim_boundary": "internal locked results only; no quantum advantage claim",
        },
        {
            "data_or_source": "p_neo_bayesian wave outputs",
            "specific_examples": "wave3 public predictor scores; wave9 algorithm forest; wave11 expanded public benchmarks",
            "role_in_review": "competitor context and missing comparator inventory",
            "used_in": "Fig 2; Fig 4; Table 2; supplementary evidence matrix",
            "local_path": f"{NEO}",
            "claim_boundary": "public predictor scores are benchmark comparators; training overlap status must be caveated",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "review_data_use_map.tsv", sep="\t", index=False)


def write_catalogs(fig_paths: dict[str, Path]) -> None:
    figs = [
        {
            "figure": "Fig 1",
            "title": "Vaccine-ready neoantigen prioritization funnel",
            "main_message": "The clinical bottleneck is top-k candidate selection for manufacturing, not just HLA binding.",
            "panels": "sequencing -> candidate peptides -> presentation filter -> immunogenicity rank -> shortlist",
            "data_used": "clinical review literature; conceptual funnel",
            "output_path": str(fig_paths["fig1"]),
            "caption": "Personalized vaccine design is a constrained ranking problem in which thousands of candidates become a short manufacturable list.",
        },
        {
            "figure": "Fig 2",
            "title": "Predictor taxonomy by biological layer",
            "main_message": "Public predictors address different biological layers; binding, immunogenicity and vaccine triage are not identical tasks.",
            "panels": "binding/presentation; immunogenicity; patient context; structure/TCR; triage",
            "data_used": "review_reference_evidence_matrix.tsv",
            "output_path": str(fig_paths["fig2"]),
            "caption": "Neoantigen prioritization requires combining biological layers rather than treating all predictors as interchangeable binding scores.",
        },
        {
            "figure": "Fig 3",
            "title": "Benchmark failure modes and required controls",
            "main_message": "Exact/near overlap, HLA shortcuts, source shift and ambiguous negatives can mimic model progress.",
            "panels": "six failure modes paired with split/control recommendations",
            "data_used": "CROSS-Neo audits plus public benchmark logic",
            "output_path": str(fig_paths["fig3"]),
            "caption": "Clean evaluation requires split families that remove reference leakage and stress source/HLA generalization.",
        },
        {
            "figure": "Fig 4",
            "title": "Metric contract and competitor gauntlet",
            "main_message": "AUPRC and top-k precision align better with vaccine candidate decisions than AUROC alone.",
            "panels": "internal strict demo top10 precision and AUPRC across public and CROSS-Neo competitors",
            "data_used": "strong_competitor_gauntlet.tsv",
            "output_path": str(fig_paths["fig4"]),
            "caption": "Internal strict demo benchmark illustrates why top-k precision and AUPRC should be reported alongside prevalence.",
        },
        {
            "figure": "Fig 5",
            "title": "Review-to-original launch flow",
            "main_message": "The review defines the standard; the CROSS-Neo original paper implements it.",
            "panels": "review standard -> original implementation",
            "data_used": "launch strategy and CROSS-Neo result directories",
            "output_path": str(fig_paths["fig5"]),
            "caption": "Two-paper strategy separates field-level evaluation standards from the original algorithmic implementation.",
        },
    ]
    tables = [
        {
            "table": "Table 1",
            "title": "Clinical neoantigen vaccine signals",
            "main_message": "Clinical vaccine signals justify why prioritization now matters.",
            "data_used": "melanoma, pancreatic cancer, RCC and selected vaccine review references",
            "output_or_plan": "review_table_plan.tsv; review_reference_evidence_matrix.tsv",
        },
        {
            "table": "Table 2",
            "title": "Predictor families and biological layers",
            "main_message": "Different predictors solve different parts of the vaccine-design problem.",
            "data_used": "BigMHC, PRIME, IMPROVE, MHCflurry, NetMHCpan, DeepImmuno, TransPHLA, T-SCAPE, MixMHCpred literature",
            "output_or_plan": "review_reference_evidence_matrix.tsv",
        },
        {
            "table": "Table 3",
            "title": "Benchmark failure modes and controls",
            "main_message": "Each leakage/shift mode requires a matching split or audit.",
            "data_used": "CROSS-Neo v0/v1/v2 audit outputs and conceptual controls",
            "output_or_plan": "review_data_use_map.tsv; CROSS-Neo audit files",
        },
        {
            "table": "Table 4",
            "title": "Recommended metrics for vaccine-ready rankers",
            "main_message": "AUPRC, top-k, enrichment, calibration and OOD coverage should be primary reporting elements.",
            "data_used": "strong_competitor_gauntlet.tsv and metric definitions",
            "output_or_plan": "strong_competitor_gauntlet.tsv",
        },
        {
            "table": "Supplementary Table 1",
            "title": "Data provenance and claim boundary map",
            "main_message": "Every dataset is mapped to where it is used and what claims it can support.",
            "data_used": "all local CROSS-Neo and review package outputs",
            "output_or_plan": "review_data_use_map.tsv",
        },
    ]
    pd.DataFrame(figs).to_csv(OUT / "review_figure_catalog.tsv", sep="\t", index=False)
    pd.DataFrame(tables).to_csv(OUT / "review_table_catalog.tsv", sep="\t", index=False)


def write_storyboard(fig_paths: dict[str, Path]) -> None:
    text = f"""# CROSS-Neo Review Storyboard

## One-Line Story

Personalized neoantigen vaccines are clinically credible again, but the field needs vaccine-ready prioritization: contamination-controlled, source-aware, top-k ranking of the few candidates worth manufacturing.

## Whole Manuscript Flow

1. **Clinical urgency**: recent vaccine studies make candidate selection practically important.
2. **Biological layering**: binding/presentation is necessary but not sufficient for immunogenicity.
3. **Benchmark problem**: overlap, HLA shortcut, source shift and ambiguous negatives can inflate rankings.
4. **Metric contract**: AUPRC, top-k precision, enrichment, calibration and abstention match the vaccine decision better than AUROC alone.
5. **Two-paper bridge**: the review defines the standard; CROSS-Neo is the implementation paper.

## Representative Figures

### Figure 1

File: `{fig_paths['fig1']}`

Message: candidate neoantigens pass through a funnel, but the final bottleneck is only 5-20 practical vaccine/assay slots.

Where used: opening section and graphical abstract.

### Figure 2

File: `{fig_paths['fig2']}`

Message: NetMHCpan/MHCflurry/MixMHCpred/BigMHC-EL model presentation; PRIME/DeepImmuno/BigMHC-IM/TransPHLA/T-SCAPE/IMPROVE move toward immunogenicity or patient context; vaccine triage needs a separate evaluation layer.

Where used: predictor taxonomy section.

### Figure 3

File: `{fig_paths['fig3']}`

Message: exact overlap, near overlap, HLA shortcuts, source shift, ambiguous negatives and AUROC-only reporting each need a matching control.

Where used: benchmark failure-mode section.

### Figure 4

File: `{fig_paths['fig4']}`

Message: internal strict demo gauntlet shows how top-k/AUPRC reporting changes practical interpretation compared with public predictor-only baselines.

Where used: metrics section and optional supplementary CROSS-Neo bridge. Claim boundary: internal retrospective example, not external validation.

### Figure 5

File: `{fig_paths['fig5']}`

Message: review first defines the field standard; original CROSS-Neo paper implements it.

Where used: author strategy, cover pitch, final review outlook.

## Representative Tables

### Table 1. Clinical vaccine signals

Purpose: show why the topic is timely. Use melanoma, pancreatic cancer and renal cell carcinoma neoantigen vaccine studies.

### Table 2. Predictor families

Purpose: separate binding/presentation predictors from immunogenicity, patient-context and triage tools.

### Table 3. Benchmark failure modes

Purpose: make contamination and source shift a field-wide methodological issue.

### Table 4. Vaccine-ready metrics

Purpose: explain why AUPRC, top-k precision, enrichment and abstention should be primary.

### Supplementary Table 1. Data-use map

Purpose: show exactly what data is used where and what claim each data source can support.

## Data Provenance

The full data-use map is in `review_data_use_map.tsv`. Critical boundaries:

- Public predictor scores are used as literature context or benchmark comparators, not clean CROSS-Neo training features.
- CROSS-Neo product/demo metrics are internal retrospective examples only.
- Clinical vaccine papers motivate the problem; they do not validate CROSS-Neo.
- Public overlap/source-shift audits are used to motivate benchmark standards, not to attack public predictors.

## Editor Pitch Version

This is not another algorithm review. It is a translational evaluation standard for vaccine-ready neoantigen prioritization.
"""
    (OUT / "REVIEW_STORYBOARD_FIGURES_TABLES_DATA.md").write_text(text)


def main() -> None:
    ensure_dirs()
    fig_paths = {
        "fig1": fig1_funnel(),
        "fig2": fig2_taxonomy(),
        "fig3": fig3_failure_modes(),
        "fig4": fig4_metrics_contract(),
        "fig5": fig5_launch_flow(),
    }
    write_data_use_map()
    write_catalogs(fig_paths)
    write_storyboard(fig_paths)
    print(f"[review-assets] wrote {OUT}")


if __name__ == "__main__":
    main()
