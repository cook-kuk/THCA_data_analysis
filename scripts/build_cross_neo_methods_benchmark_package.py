#!/usr/bin/env python3
"""Build the CROSS-Neo methods/benchmarking manuscript package.

This is the second-paper package after the high-impact review/perspective.
It frames CROSS-Neo as a contamination-controlled benchmark/validation system
for vaccine-ready neoantigen prioritization, not as an external-valid product
claim or quantum-advantage paper.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "project/results/cross_neo_methods_benchmark_2026_05_10"
FIG = OUT / "figures"
REVIEW = REPO / "project/results/cross_neo_review_2026_05_10"
PRODUCT = REPO / "project/results/cross_neo_product_demo_2026_05_10"
V0 = REPO / "project/results/cross_neo_v0"
V1 = REPO / "project/results/cross_neo_v1"
V2 = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
NEO = REPO / "project/results/p_neo_bayesian_2026_05_09"


COLORS = {
    "blue": "#2f6f9f",
    "green": "#4b8f73",
    "amber": "#b36b2c",
    "purple": "#5b5b8a",
    "grey": "#6b7280",
    "dark": "#17202a",
    "red": "#9a3e3e",
}


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def add_box(ax, x, y, w, h, text, color, fs=9.5):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.015,rounding_size=0.035",
        facecolor=color,
        edgecolor="none",
        alpha=0.96,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color="white", fontsize=fs, weight="bold")


def add_arrow(ax, p1, p2):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=14, color=COLORS["grey"], linewidth=1.9))


def fig1_validation_system() -> Path:
    path = FIG / "fig1_cross_neo_validation_system.png"
    fig, ax = plt.subplots(figsize=(14, 7.5))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Figure 1. CROSS-Neo as a validation system for vaccine-ready prioritization", fontsize=16, weight="bold", color=COLORS["dark"])
    ax.text(0.04, 0.89, "The paper tests whether a prioritization stack can survive leakage, source shift and top-k decision constraints.", fontsize=10.5, color=COLORS["grey"])

    add_box(ax, 0.05, 0.66, 0.18, 0.12, "Inputs\nmutant peptide\nWT peptide\nHLA\nsource context", COLORS["blue"], 9)
    add_arrow(ax, (0.24, 0.72), (0.31, 0.72))
    add_box(ax, 0.32, 0.66, 0.18, 0.12, "Branches\ncounterfactual\nhard-decoy\nQK/sequence\nstructure/OOD", COLORS["green"], 9)
    add_arrow(ax, (0.51, 0.72), (0.58, 0.72))
    add_box(ax, 0.59, 0.66, 0.18, 0.12, "Fold-safe gate\ntrain-only scaling\ncalibration\nretrieval audit", COLORS["purple"], 9)
    add_arrow(ax, (0.78, 0.72), (0.85, 0.72))
    add_box(ax, 0.84, 0.66, 0.12, 0.12, "Top-k\ncandidate\nqueue", COLORS["amber"], 9)

    controls = [
        ("exact peptide-HLA\nholdout", 0.07, 0.38),
        ("near-peptide\ncluster holdout", 0.24, 0.38),
        ("HLA/supertype\nholdout", 0.41, 0.38),
        ("source/study\nholdout", 0.58, 0.38),
        ("public-overlap\naudit", 0.75, 0.38),
    ]
    for text, x, y in controls:
        add_box(ax, x, y, 0.13, 0.10, text, COLORS["red"], 8.5)

    ax.text(0.06, 0.16, "Primary outputs: AUPRC, top5/top10 precision, enrichment over prevalence, calibration and abstention.", fontsize=11, color=COLORS["dark"], weight="bold")
    ax.text(0.06, 0.10, "Claim boundary: internal locked validation and benchmark system; no external validation or quantum-advantage claim.", fontsize=10.5, color=COLORS["grey"])
    fig.tight_layout()
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def fig2_split_contract() -> Path:
    path = FIG / "fig2_leakage_safe_split_contract.png"
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Figure 2. Split-safe evaluation contract", fontsize=16, weight="bold", color=COLORS["dark"])

    add_box(ax, 0.08, 0.70, 0.24, 0.13, "Training fold only\nretrieval index\nscalers\ncalibrators\nfeature selection", COLORS["green"], 9)
    add_box(ax, 0.39, 0.70, 0.22, 0.13, "Locked test fold\nno index entry\nno weight tuning\nno threshold tuning", COLORS["red"], 9)
    add_box(ax, 0.68, 0.70, 0.22, 0.13, "Audit layer\nexact/near hits\nsource overlap\npublic corpus", COLORS["purple"], 9)
    add_arrow(ax, (0.32, 0.765), (0.38, 0.765))
    add_arrow(ax, (0.61, 0.765), (0.67, 0.765))

    rows = [
        ("allowed inside train fold", "retrieval evidence, scaling, nested fusion weights, calibration"),
        ("forbidden from test fold", "feature selection, gamma selection, threshold tuning, target prevalence, source label rate"),
        ("split families", "repeated, HLA-stratified, HLA-supertype, exact peptide-HLA, near peptide, source/study, time"),
        ("reporting", "public overlap flags, clean/reference-supported strata, OOD and abstention metrics"),
    ]
    y = 0.46
    for left, right in rows:
        ax.text(0.10, y, left, fontsize=10.5, weight="bold", color=COLORS["dark"])
        ax.text(0.36, y, right, fontsize=10.2, color=COLORS["dark"])
        y -= 0.09

    fig.tight_layout()
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def fig3_locked_performance() -> Path:
    path = FIG / "fig3_locked_performance_gauntlet.png"
    df = read_tsv(PRODUCT / "strong_competitors/strong_competitor_gauntlet.tsv")
    fig, ax = plt.subplots(figsize=(13.5, 8))
    if df.empty:
        ax.text(0.2, 0.5, "strong competitor gauntlet missing", fontsize=14)
    else:
        strict = df[df["benchmark_scope"] == "strict_product_demo_set_n89"].copy()
        strict = strict.sort_values(["top10_precision", "AUPRC"], ascending=False).head(14)
        labels = strict["method"].str.replace("_", " ", regex=False)
        y = np.arange(len(strict))
        colors = np.where(strict["method_family"].str.contains("CROSS", na=False), COLORS["blue"], COLORS["amber"])
        ax.barh(y, strict["top10_precision"], color=colors, alpha=0.9, label="top10 precision")
        ax.scatter(strict["AUPRC"], y, s=40, color=COLORS["dark"], label="AUPRC")
        ax.axvline(strict["prevalence"].iloc[0], color=COLORS["grey"], linestyle="--", linewidth=1.2, label="prevalence")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8.5)
        ax.invert_yaxis()
        ax.set_xlim(0, 1.02)
        ax.set_title("Figure 3. Internal locked product-demo gauntlet", loc="left", fontsize=15, weight="bold")
        ax.set_xlabel("Metric value")
        ax.legend(frameon=False, loc="lower right")
        ax.grid(axis="x", alpha=0.15)
    fig.tight_layout()
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def fig4_source_stress() -> Path:
    path = FIG / "fig4_source_stress_and_abstention_boundary.png"
    rows = []
    for file, label in [
        (V1 / "hard_decoy_focal_metrics.tsv", "hard_decoy"),
        (V1 / "source_balanced_metrics.tsv", "source_balanced"),
        (V2 / "metrics/source_heldout_metrics.tsv", "v2_source"),
        (V2 / "metrics/source_collapse_rescue_summary.tsv", "v2_rescue"),
    ]:
        df = read_tsv(file)
        if not df.empty:
            df["source_file_label"] = label
            rows.append(df)
    fig, ax = plt.subplots(figsize=(13, 7.5))
    if not rows:
        ax.text(0.2, 0.5, "source stress metrics missing", fontsize=14)
    else:
        allm = pd.concat(rows, ignore_index=True, sort=False)
        split_col = "split_name" if "split_name" in allm.columns else "split"
        model_col = "model" if "model" in allm.columns else "model_name" if "model_name" in allm.columns else None
        if model_col is None or split_col not in allm.columns:
            ax.text(0.1, 0.5, "source stress schema not recognized", fontsize=14)
        else:
            keep = allm[allm[split_col].astype(str).str.contains("source_heldout|TESLA|NEPdb|CEDAR", case=False, na=False)].copy()
            if keep.empty:
                keep = allm.copy().head(20)
            # Choose the best top10 per source-ish split.
            top_col = "top10_precision" if "top10_precision" in keep.columns else "top10"
            auprc_col = "AUPRC" if "AUPRC" in keep.columns else "auprc"
            keep[top_col] = pd.to_numeric(keep[top_col], errors="coerce")
            keep[auprc_col] = pd.to_numeric(keep[auprc_col], errors="coerce")
            best = keep.sort_values([split_col, top_col, auprc_col], ascending=[True, False, False]).groupby(split_col).head(1)
            best = best.head(12)
            labels = best[split_col].astype(str).str.replace("source_heldout_", "", regex=False)
            y = np.arange(len(best))
            ax.barh(y, best[top_col], color=COLORS["red"], alpha=0.78, label="top10 precision")
            ax.scatter(best[auprc_col], y, color=COLORS["dark"], s=42, label="AUPRC")
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=9)
            ax.invert_yaxis()
            ax.set_xlim(0, 1.02)
            ax.set_title("Figure 4. Source-heldout stress: where the model should abstain", loc="left", fontsize=15, weight="bold")
            ax.set_xlabel("Metric value")
            ax.legend(frameon=False, loc="lower right")
            ax.grid(axis="x", alpha=0.15)
    fig.tight_layout()
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def fig5_claim_boundary() -> Path:
    path = FIG / "fig5_claim_boundary_decision_tree.png"
    fig, ax = plt.subplots(figsize=(13.5, 7))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.04, 0.94, "Figure 5. Claim boundary decision tree", fontsize=16, weight="bold", color=COLORS["dark"])

    add_box(ax, 0.06, 0.72, 0.24, 0.12, "Was the split locked\nand train-only safe?", COLORS["purple"], 9)
    add_box(ax, 0.39, 0.72, 0.24, 0.12, "Does top-k improve\nwithout exact retrieval?", COLORS["purple"], 9)
    add_box(ax, 0.72, 0.72, 0.22, 0.12, "Does source/HLA stress\nremain above prevalence?", COLORS["purple"], 9)
    add_arrow(ax, (0.31, 0.78), (0.38, 0.78))
    add_arrow(ax, (0.64, 0.78), (0.71, 0.78))

    add_box(ax, 0.08, 0.42, 0.22, 0.11, "Yes\ninternal locked claim", COLORS["green"], 9)
    add_box(ax, 0.39, 0.42, 0.22, 0.11, "Partial\nsupplementary-ready", COLORS["amber"], 9)
    add_box(ax, 0.70, 0.42, 0.22, 0.11, "No\nhold / abstain", COLORS["red"], 9)
    ax.text(0.08, 0.22, "Forbidden statements:\nexternal validation\nquantum advantage\nuniversal superiority", fontsize=11, color=COLORS["dark"], weight="bold")
    ax.text(0.50, 0.22, "Allowed statements:\ninternal locked benchmark\npublic-comparator separation\nsource stress limitations\nvaccine triage utility", fontsize=11, color=COLORS["dark"], weight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def write_target_matrix() -> None:
    rows = [
        {
            "rank": 1,
            "journal": "Nature Biomedical Engineering",
            "format": "Article",
            "title": "Benchmarking vaccine-ready AI systems for neoantigen prioritization",
            "fit": "validation/deployment of computational systems for healthcare",
            "need_to_strengthen": "position as system validation, not just model leaderboard; keep main text to 3500 words and <=8 display items",
            "risk": "must show broad biomedical engineering significance beyond small-n immunology",
            "move": "presubmission enquiry with Article framing and full benchmark package",
        },
        {
            "rank": 2,
            "journal": "Nature Machine Intelligence",
            "format": "Analysis",
            "title": "Contamination-controlled benchmarking of neoantigen prioritization models",
            "fit": "comparative analysis of existing data with new conclusions for AI evaluation",
            "need_to_strengthen": "more public tool reproducibility, bootstrap/permutation, larger benchmark registry",
            "risk": "small n and domain-specific oncology could weaken broad AI interest",
            "move": "prepare Analysis synopsis; emphasize benchmark design and leakage",
        },
        {
            "rank": 3,
            "journal": "Cell Reports Medicine",
            "format": "Resource/Article",
            "title": "A contamination-controlled benchmark for personalized neoantigen vaccine prioritization",
            "fit": "translational clinical AI and immunotherapy decision support",
            "need_to_strengthen": "clinical framing, patient/source analysis, practical candidate queue",
            "risk": "needs stronger clinical utility and preferably prospective/independent data",
            "move": "secondary prestige route after NBE/NMI",
        },
        {
            "rank": 4,
            "journal": "Patterns",
            "format": "Resource/Article",
            "title": "A benchmark and model-card framework for neoantigen prioritization",
            "fit": "data-centric AI, reproducibility, benchmarks, model cards",
            "need_to_strengthen": "package registry, code release, reproducible score matrix, model cards",
            "risk": "less clinical prestige",
            "move": "strong realistic route",
        },
        {
            "rank": 5,
            "journal": "JITC",
            "format": "Original Research or Methods/Review hybrid",
            "title": "Contamination-controlled neoantigen prioritization for vaccine candidate triage",
            "fit": "immunotherapy audience and practical vaccine candidate ranking",
            "need_to_strengthen": "immunology interpretation and public predictor comparison",
            "risk": "lower AI/benchmark prestige",
            "move": "high probability fallback",
        },
        {
            "rank": 6,
            "journal": "npj Precision Oncology",
            "format": "Article",
            "title": "Pan-allele neoantigen prioritization under source shift and label uncertainty",
            "fit": "precision oncology computational method",
            "need_to_strengthen": "clinical precision-oncology framing",
            "risk": "less methods prestige",
            "move": "safe backup",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "methods_benchmark_target_matrix.tsv", sep="\t", index=False)


def write_synopses() -> None:
    nbe = """# Nature Biomedical Engineering Presubmission Synopsis

## Proposed Format
Article

## Proposed Title
Benchmarking vaccine-ready AI systems for neoantigen prioritization

## 175-Word Abstract Draft
Personalized neoantigen vaccines require computational systems that choose a small number of candidate peptides for manufacture or immune testing. Yet most neoantigen predictors are benchmarked as retrospective classifiers, often under public-corpus overlap, source shift, HLA shortcuts and incomplete negative labels. Here, we present CROSS-Neo, a contamination-controlled benchmarking and prioritization system for vaccine-ready neoantigen ranking. CROSS-Neo integrates mutant-WT counterfactual sequence features, hard-decoy pressure, fixed sequence/QK fallback signals, explicit retrieval-risk flags, source-aware calibration and OOD-aware abstention. We evaluate the system using locked split families that remove exact peptide-HLA overlap, near-peptide overlap, HLA/supertype shortcuts and source/study leakage. Primary outcomes are AUPRC, top-k precision, enrichment over prevalence and abstention behavior rather than AUROC alone. In internal retrospective benchmarks, the pan-allele product ranker improves top-k prioritization over public predictor-only baselines, while source-heldout stress tests reveal where abstention remains necessary. CROSS-Neo is therefore presented as a benchmarked validation framework, not an external-valid clinical product.

## Editorial Pitch
The article fits Nature Biomedical Engineering because it concerns the validation and deployment requirements of a computational system that may facilitate personalized therapy design. The central contribution is not a black-box model, but a leakage-controlled evaluation architecture and reporting standard for AI-guided vaccine candidate selection.

## Display Items
1. CROSS-Neo validation-system overview.
2. Leakage-safe split contract.
3. Internal locked competitor gauntlet.
4. Source-heldout stress and abstention boundary.
5. Claim-boundary decision tree.
6. NEO-PRIOR reporting checklist.
"""
    (OUT / "NATURE_BIOMEDICAL_ENGINEERING_SYNOPSIS.md").write_text(nbe)

    nmi = """# Nature Machine Intelligence Analysis Synopsis

## Proposed Format
Analysis

## Proposed Title
Contamination-controlled benchmarking of neoantigen prioritization models

## Core Claim
Neoantigen prediction leaderboards can change substantially when models are evaluated under exact/near overlap controls, source/HLA holdouts and top-k metrics that reflect vaccine candidate selection.

## Why Analysis
Nature Machine Intelligence defines Analysis as a new analysis of existing data or comparative analysis leading to novel conclusions for a broad audience. This manuscript would be framed as a comparative AI evaluation study, not a wet-lab validation paper.

## What Must Be Strengthened Before Submission
- Expand public competitor coverage where possible: IMPROVE, MixMHCpred standalone, DeepHLApan, GraphMHC if runnable locally.
- Add bootstrap/permutation tests for the final candidate.
- Add model-card and data-card package.
- Make the benchmark registry reproducible from raw or intermediate files.
- Keep biological claims modest and emphasize AI evaluation.
"""
    (OUT / "NATURE_MACHINE_INTELLIGENCE_ANALYSIS_SYNOPSIS.md").write_text(nmi)

    patterns = """# Patterns / Data-Centric AI Synopsis

## Proposed Title
A benchmark and model-card framework for neoantigen prioritization

## Core Fit
Patterns is the strongest realistic route if the paper is packaged as a reproducible benchmark, score matrix, model-card and data-card framework.

## Required Package
- canonical benchmark registry;
- public predictor score matrix;
- clean/comparator feature separation;
- split manifests;
- model cards;
- data provenance tables;
- failure-mode report;
- reproducible scripts and figures.
"""
    (OUT / "PATTERNS_SYNOPSIS.md").write_text(patterns)


def write_emails() -> None:
    text = """# Methods/Benchmarking Editor Emails

## Nature Biomedical Engineering

Subject: Presubmission enquiry: Article on benchmarking vaccine-ready AI systems for neoantigen prioritization

Dear Editors,

I am writing to ask whether Nature Biomedical Engineering would consider an Article provisionally titled "Benchmarking vaccine-ready AI systems for neoantigen prioritization."

Personalized neoantigen vaccines require computational systems that choose a small number of candidate peptides for manufacture or immune testing. However, most neoantigen predictors are evaluated as retrospective classifiers, often without explicit controls for public-corpus overlap, near-peptide leakage, source/study shift, HLA shortcuts or incomplete negative labels.

We present CROSS-Neo as a contamination-controlled benchmarking and prioritization system for vaccine-ready neoantigen ranking. The emphasis is on validation and deployment requirements for a computational system: locked split families, train-only retrieval and calibration, public-comparator separation, AUPRC/top-k/enrichment primary metrics, source-heldout stress tests and OOD-aware abstention. The manuscript would explicitly avoid external-validation and quantum-advantage claims.

I attach a short synopsis and figure plan and would be grateful for your advice on whether this would fit the journal.

Sincerely,

Seungho Cook

## Nature Machine Intelligence

Subject: Presubmission enquiry: Analysis on contamination-controlled neoantigen model benchmarking

Dear Editors,

I am writing to ask whether Nature Machine Intelligence would consider an Analysis provisionally titled "Contamination-controlled benchmarking of neoantigen prioritization models."

The manuscript asks a general AI-evaluation question in a clinically important domain: how do model rankings change when overlap-sensitive biomedical prediction benchmarks are evaluated using train-only retrieval controls, near-neighbor holdouts, source/HLA stress tests and top-k utility metrics rather than AUROC alone?

Using neoantigen vaccine prioritization as the case study, we compare public predictor baselines, internal sequence/structure/QK branches and leakage-aware fusion variants under locked split families. The aim is not to claim external clinical validation, but to show how benchmark design and reporting standards alter apparent model utility for personalized vaccine candidate selection.

Sincerely,

Seungho Cook
"""
    (OUT / "METHODS_BENCHMARK_EDITOR_EMAILS.md").write_text(text)


def write_skeleton() -> None:
    text = """# CROSS-Neo Methods/Benchmarking Manuscript Skeleton

## Working Title
Benchmarking vaccine-ready AI systems for neoantigen prioritization

## Claim Boundary
CROSS-Neo is an internal locked benchmark and prioritization framework. It does not claim external validation, clinical efficacy or quantum advantage.

## Abstract
Personalized neoantigen vaccines require computational systems that select a small number of candidates for manufacture or immune testing, but many predictors are evaluated as retrospective classifiers under overlap-sensitive public benchmarks. We introduce CROSS-Neo, a contamination-controlled benchmark and pan-allele prioritization stack for vaccine-ready neoantigen ranking. The system combines mutant-WT counterfactual features, hard-decoy pressure, fixed sequence/QK fallback signals, explicit retrieval-risk flags, source-aware calibration and OOD-aware abstention. Evaluation uses locked split families, including exact peptide-HLA, near-peptide, HLA/supertype, source/study and public-overlap audits. Primary outcomes are AUPRC, top-k precision, enrichment over prevalence and abstention behavior. Internal retrospective benchmarks show improved practical top-k prioritization over public predictor-only baselines, while source-heldout stress tests identify contexts where abstention remains necessary.

## Results

### A leakage-safe benchmark registry
- master table, split manifests, public overlap audit, retrieval risk flags.

### Public predictor baselines are strong but overlap-sensitive
- BigMHC, MHCflurry, PRIME, DeepImmuno, TransPHLA, NetMHCpan, T-SCAPE, MHCnuggets, NetMHCstabpan.
- Public predictor scores are comparators, not clean training features.

### CROSS-Neo branches have complementary error profiles
- counterfactual RF, QK fixed fallback, structure/geometry, retrieval evidence.
- naive concatenation fails; gated/fixed fusion is safer.

### Hard-decoy and pan-allele fusion improve top-k triage
- product fixed pan-allele score improves over public-only baselines in internal strict demo set.

### Source-heldout stress reveals abstention boundaries
- NEPdb improves, TESLA remains difficult.
- This is a limitation and deployment guardrail, not a failure to hide.

### Reporting standard and model-card outputs
- NEO-PRIOR checklist, data provenance, claim boundary, candidate queue.

## Methods
- data sources and inclusion criteria;
- feature construction;
- split design;
- retrieval index safety;
- calibration and abstention;
- public predictor score handling;
- metrics and statistics;
- software and reproducibility.

## Display Items
See `methods_figure_catalog.tsv`.

## Tables
See `methods_table_catalog.tsv`.
"""
    (OUT / "METHODS_BENCHMARK_MANUSCRIPT_SKELETON.md").write_text(text)


def write_catalogs(figs: dict[str, Path]) -> None:
    figure_rows = [
        {
            "figure": "Fig 1",
            "title": "CROSS-Neo validation system",
            "main_message": "CROSS-Neo is a leakage-controlled validation/prioritization stack, not just a classifier.",
            "data_used": "method architecture and result directory inventory",
            "path": str(figs["fig1"]),
        },
        {
            "figure": "Fig 2",
            "title": "Split-safe evaluation contract",
            "main_message": "All retrieval, calibration, feature selection and thresholding must be train-fold only.",
            "data_used": "CROSS-Neo split and audit design",
            "path": str(figs["fig2"]),
        },
        {
            "figure": "Fig 3",
            "title": "Internal locked competitor gauntlet",
            "main_message": "AUPRC/top-k metrics show practical prioritization gains over public-only baselines.",
            "data_used": "strong_competitor_gauntlet.tsv",
            "path": str(figs["fig3"]),
        },
        {
            "figure": "Fig 4",
            "title": "Source-heldout stress and abstention boundary",
            "main_message": "Source shift defines where the system should abstain rather than overclaim.",
            "data_used": "v1/v2 source-heldout metrics",
            "path": str(figs["fig4"]),
        },
        {
            "figure": "Fig 5",
            "title": "Claim boundary decision tree",
            "main_message": "Internal locked claims are allowed; external validation and quantum advantage are not.",
            "data_used": "decision/reporting logic",
            "path": str(figs["fig5"]),
        },
    ]
    table_rows = [
        {
            "table": "Table 1",
            "title": "Benchmark datasets and claim boundaries",
            "data_used": "master_table.tsv, registry, product demo, v2 registry",
        },
        {
            "table": "Table 2",
            "title": "Public predictor competitors and training-feature status",
            "data_used": "strong competitor gauntlet and public overlap audit",
        },
        {
            "table": "Table 3",
            "title": "Locked split performance",
            "data_used": "v1_model_comparison.tsv, v2 metrics, hard_decoy metrics",
        },
        {
            "table": "Table 4",
            "title": "Source-heldout failure and abstention modes",
            "data_used": "source_heldout metrics and source bias diagnostics",
        },
        {
            "table": "Supplementary Table 1",
            "title": "NEO-PRIOR compliance checklist",
            "data_used": "NEO_PRIOR_reporting_checklist.tsv",
        },
    ]
    pd.DataFrame(figure_rows).to_csv(OUT / "methods_figure_catalog.tsv", sep="\t", index=False)
    pd.DataFrame(table_rows).to_csv(OUT / "methods_table_catalog.tsv", sep="\t", index=False)


def write_data_map() -> None:
    rows = [
        {
            "source": "cross_neo_v0",
            "local_path": str(V0),
            "used_for": "master table, locked splits, retrieval leakage audit, baseline features, v0 metrics",
            "claim_boundary": "internal baseline and audit source",
        },
        {
            "source": "cross_neo_v1",
            "local_path": str(V1),
            "used_for": "branch error decomposition, gated MoE, source bias correction, PU ranking, hard-decoy repair",
            "claim_boundary": "internal method development and locked evaluation",
        },
        {
            "source": "cross_neo_v2_sota_sprint",
            "local_path": str(V2),
            "used_for": "selective ensemble, ESM2/QK gate, source stress, public overlap, model comparison",
            "claim_boundary": "internal locked/sprint results; no external validation",
        },
        {
            "source": "p_neo_bayesian public predictor waves",
            "local_path": str(NEO),
            "used_for": "public predictor competitor outputs and public benchmark context",
            "claim_boundary": "competitors only; public scores not clean training features",
        },
        {
            "source": "cross_neo_product_demo",
            "local_path": str(PRODUCT),
            "used_for": "product fixed pan-allele score, customer-safe queue, strong competitor gauntlet",
            "claim_boundary": "internal retrospective product demo; not external validation",
        },
        {
            "source": "cross_neo_review",
            "local_path": str(REVIEW),
            "used_for": "NEO-PRIOR standard and review-to-original bridge",
            "claim_boundary": "field standard and positioning; not original evidence",
        },
    ]
    pd.DataFrame(rows).to_csv(OUT / "methods_data_use_map.tsv", sep="\t", index=False)


def write_summary() -> None:
    summary = {
        "out_dir": str(OUT),
        "recommended_first_target": "Nature Biomedical Engineering Article",
        "recommended_second_target": "Nature Machine Intelligence Analysis",
        "paper_role": "methods/benchmarking companion after review/perspective",
        "forbidden_claims": ["external validation", "quantum advantage", "universal superiority"],
    }
    (OUT / "methods_benchmark_package_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main() -> None:
    ensure_dirs()
    figs = {
        "fig1": fig1_validation_system(),
        "fig2": fig2_split_contract(),
        "fig3": fig3_locked_performance(),
        "fig4": fig4_source_stress(),
        "fig5": fig5_claim_boundary(),
    }
    write_target_matrix()
    write_synopses()
    write_emails()
    write_skeleton()
    write_catalogs(figs)
    write_data_map()
    write_summary()
    print(f"[methods-benchmark] wrote {OUT}")


if __name__ == "__main__":
    main()
