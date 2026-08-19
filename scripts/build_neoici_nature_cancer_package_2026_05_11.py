#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "project" / "results" / "p_neo_ici_nature_cancer_package_2026_05_11"
HUB = ROOT / "project" / "papers_hub_2026_05_04"
WWW = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET_NAME = "neoici_nature_cancer_package"


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def copy_bytes(src: Path, dst: Path) -> None:
    dst.write_bytes(src.read_bytes())


def save_fig(path: Path, fig) -> None:
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def format_num(value, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, (int, float)):
        return f"{value:.{digits}f}"
    return str(value)


def block_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    return df.to_html(index=False, escape=False, classes="dataframe")


def main() -> None:
    ensure_dir(OUT)
    ensure_dir(HUB)
    ensure_dir(WWW)

    gauntlet = load_tsv(ROOT / "project" / "results" / "p_neo_ici_combo_2026_05_11" / "neoici_algorithm_gauntlet_metrics.tsv")
    winloss = load_tsv(ROOT / "project" / "results" / "p_neo_ici_combo_2026_05_11" / "neoici_algorithm_winloss.tsv")
    combo = load_tsv(ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "combo_geo_reanalysis.tsv")
    proxy = load_tsv(ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "proxy_compare" / "proxy_compare_metrics.tsv")
    proxy_samples = load_tsv(ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "proxy_compare" / "proxy_compare_samples.tsv")
    sample = load_tsv(ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "NEOICI_NATURE_CANCER_SAMPLE_TABLE.tsv")
    loo = load_tsv(ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "NEOICI_NATURE_CANCER_LOO_TABLE.tsv")
    claim = load_tsv(ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "NEOICI_NATURE_CANCER_CLAIM_TABLE.tsv")
    contract = load_tsv(ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "NEOICI_NATURE_CANCER_SYNTHESIS_TABLE.tsv")
    executive = load_tsv(ROOT / "project" / "results" / "p_neo_ici_nature_cancer_package_2026_05_11" / "NEOICI_NATURE_CANCER_EXECUTIVE_TABLE.tsv")

    figure_rows = [
        {
            "panel": "Fig. 0O",
            "asset": "impact banner",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_impact_banner.png",
            "purpose": "Single top-screen banner that reinforces the claim without adding clutter.",
        },
        {
            "panel": "Fig. 0N",
            "asset": "micro strip",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_micro_strip.png",
            "purpose": "Shortest numeric strip for the absolute first visual hit.",
        },
        {
            "panel": "Fig. 0M",
            "asset": "claim strip",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_claim_strip.png",
            "purpose": "Shortest possible claim strip for the very first screen.",
        },
        {
            "panel": "Fig. 0L",
            "asset": "topline banner",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_topline_banner.png",
            "purpose": "Shortest possible opening banner for immediate claim lock.",
        },
        {
            "panel": "Fig. 0K",
            "asset": "launch banner",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_launch_banner.png",
            "purpose": "Topmost launch banner that hits the claim before anything else.",
        },
        {
            "panel": "Fig. 0J",
            "asset": "cover slab",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_cover_slab.png",
            "purpose": "Largest front-page slab with the clearest claim boundary.",
        },
        {
            "panel": "Fig. 0H",
            "asset": "opening impact board",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_opening_impact_board.png",
            "purpose": "Top-screen impact board that locks the page verdict immediately.",
        },
        {
            "panel": "Fig. 0I",
            "asset": "verdict wall",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_verdict_wall.png",
            "purpose": "Largest first-screen verdict wall for immediate impact.",
        },
        {
            "panel": "Fig. 0A",
            "asset": "evidence stack",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_evidence_stack.png",
            "purpose": "Single-line evidence stack for the whole page.",
        },
        {
            "panel": "Fig. 0B",
            "asset": "verdict poster",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_verdict_poster.png",
            "purpose": "Hero-level claim boundary and summary poster.",
        },
        {
            "panel": "Fig. 0C",
            "asset": "claim poster",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_claim_poster.png",
            "purpose": "Full-width claim poster for the first screen.",
        },
        {
            "panel": "Fig. 0D",
            "asset": "front-page collage",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_frontpage_collage.png",
            "purpose": "First-screen collage across claim, benchmark, replay, and case.",
        },
        {
            "panel": "Fig. 0E",
            "asset": "super poster",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_super_poster.png",
            "purpose": "One-screen ultra-summary poster.",
        },
        {
            "panel": "Fig. 0F",
            "asset": "monument poster",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_monument_poster.png",
            "purpose": "Dark hero-style monument poster for the first screen.",
        },
        {
            "panel": "Fig. 0G",
            "asset": "cover wall",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_cover_wall.png",
            "purpose": "True first-screen cover wall for the dossier.",
        },
        {
            "panel": "Fig. 1",
            "asset": "dataset ladder",
            "path": f"project/results/p_neo_ici_combo_2026_05_11/figures/Fig_NeoICI_combo_dataset_ladder.png",
            "purpose": "Public combo dataset registry and evidence ladder.",
        },
        {
            "panel": "Fig. 2",
            "asset": "algorithm gauntlet",
            "path": f"project/results/p_neo_ici_combo_2026_05_11/figures/Fig_NeoICI_algorithm_gauntlet_AUPRC.png",
            "purpose": "NeoICI vs existing algorithms.",
        },
        {
            "panel": "Fig. 2A",
            "asset": "algorithm rank bars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_algorithm_rankbars.png",
            "purpose": "Human-readable rank order on the validation-like split.",
        },
        {
            "panel": "Fig. 2B",
            "asset": "algorithm heatmap",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_algorithm_heatmap.png",
            "purpose": "Top algorithm metrics as a simple heatmap.",
        },
        {
            "panel": "Fig. 2C",
            "asset": "win/loss heatmap",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_winloss_heatmap.png",
            "purpose": "Pairwise win/loss heatmap across comparators.",
        },
        {
            "panel": "Fig. 2D",
            "asset": "family mean bars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_algorithm_family_bars.png",
            "purpose": "Mean AUPRC by algorithm family.",
        },
        {
            "panel": "Fig. 2E",
            "asset": "comparator winbars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_comparator_winbars.png",
            "purpose": "Pairwise comparator win fractions.",
        },
        {
            "panel": "Fig. 2F",
            "asset": "algorithm delta bars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_algorithm_delta_bars.png",
            "purpose": "AUPRC delta bars relative to KG_GA_evolved.",
        },
        {
            "panel": "Fig. 2G",
            "asset": "precision profile",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_precision_profile.png",
            "purpose": "Top-k precision profile for top algorithms.",
        },
        {
            "panel": "Fig. 2H",
            "asset": "split generalization heatmap",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_split_generalization_heatmap.png",
            "purpose": "Cross-split generalization matrix for selected algorithms.",
        },
        {
            "panel": "Fig. 3",
            "asset": "open combo replay",
            "path": f"project/results/p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11/figures/Fig_combo_response_index.png",
            "purpose": "Sample-level response-like phenotype across open combo cohorts.",
        },
        {
            "panel": "Fig. 3B",
            "asset": "dataset comparison bars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_dataset_comparison_bars.png",
            "purpose": "Dataset-level median comparison across replay cohorts.",
        },
        {
            "panel": "Fig. 3A",
            "asset": "proxy scatter",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_proxy_scatter.png",
            "purpose": "Same samples, two proxy families, easier to see by eye.",
        },
        {
            "panel": "Fig. 4",
            "asset": "proxy comparison",
            "path": f"project/results/p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11/proxy_compare/figures/Fig_combo_proxy_compare.png",
            "purpose": "NeoICI-style vs NeoPrecis-style proxy on the same samples.",
        },
        {
            "panel": "Fig. 4A",
            "asset": "proxy loo",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_proxy_loo.png",
            "purpose": "Leave-one-out robustness made human-readable.",
        },
        {
            "panel": "Fig. 4B",
            "asset": "proxy summary bars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_proxy_summary_bars.png",
            "purpose": "Proxy family summary bars for quick reading.",
        },
        {
            "panel": "Fig. 4C",
            "asset": "proxy loo",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_proxy_loo.png",
            "purpose": "Leave-one-out robustness made human-readable.",
        },
        {
            "panel": "Fig. 4D",
            "asset": "sample rank profile",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_sample_rank_profile.png",
            "purpose": "Rank profile across the 9 replay samples.",
        },
        {
            "panel": "Fig. 4E",
            "asset": "sample rank heatmap",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_sample_rank_heatmap.png",
            "purpose": "Heatmap of response / NeoICI / NeoPrecis ranks.",
        },
        {
            "panel": "Fig. 4F",
            "asset": "proxy correlation heatmap",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_proxy_corr_heatmap.png",
            "purpose": "Correlation heatmap across replay state/proxy variables.",
        },
        {
            "panel": "Fig. 4G",
            "asset": "sample rank mismatch bars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_sample_rank_mismatch_bars.png",
            "purpose": "Absolute rank mismatch bars vs response rank.",
        },
        {
            "panel": "Fig. 4H",
            "asset": "failure lanes",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_failure_lanes.png",
            "purpose": "Top-lane vs low-lane contrast for honest negative-control reading.",
        },
        {
            "panel": "Fig. 5",
            "asset": "CTMS1_followup4 case",
            "path": f"project/results/p_neo_ici_nature_cancer_synthesis_2026_05_11/ctms1_followup4_case/figures/Fig_CTMS1_followup4_case_study.png",
            "purpose": "Best open combo lane explanation panel.",
        },
        {
            "panel": "Fig. 5A",
            "asset": "CTMS1 case profile",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_CTMS1_case_profile.png",
            "purpose": "CTMS1_followup4 versus cohort median, easy to read by eye.",
        },
        {
            "panel": "Fig. 5B",
            "asset": "attention heatmaps",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_attention_heatmaps.png",
            "purpose": "Five attention-style heatmaps from the replay signals.",
        },
        {
            "panel": "Fig. 5C",
            "asset": "replay sample heatmap",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_sample_heatmap.png",
            "purpose": "Replay sample heatmap across the 9 retained samples.",
        },
        {
            "panel": "Fig. 5D",
            "asset": "claim boundary",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_NeoICI_claim_boundary.png",
            "purpose": "Claim boundary diagram for allowed vs forbidden statements.",
        },
        {
            "panel": "Fig. 5E",
            "asset": "case delta strip",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_CTMS1_case_delta.png",
            "purpose": "CTMS1_followup4 minus cohort median deltas.",
        },
        {
            "panel": "Fig. 5F",
            "asset": "case radar",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_CTMS1_case_radar.png",
            "purpose": "Radar summary of CTMS1_followup4 versus cohort median.",
        },
        {
            "panel": "Fig. 5G",
            "asset": "case neighbor bars",
            "path": f"project/results/p_neo_ici_nature_cancer_package_2026_05_11/figures/Fig_CTMS1_case_neighbor_bars.png",
            "purpose": "Nearest-neighbor style comparison for CTMS1_followup4.",
        },
    ]
    figure_df = pd.DataFrame(figure_rows)
    figure_df.to_csv(OUT / "NEOICI_NATURE_CANCER_FIGURE_MAP.tsv", sep="\t", index=False)

    supp_rows = [
        {"item": "S1", "file": "neoici_algorithm_gauntlet_metrics.tsv", "role": "full algorithm benchmark table"},
        {"item": "S2", "file": "neoici_algorithm_winloss.tsv", "role": "pairwise algorithm win/loss matrix"},
        {"item": "S3", "file": "combo_geo_reanalysis.tsv", "role": "open combo GEO sample table"},
        {"item": "S4", "file": "proxy_compare_metrics.tsv", "role": "proxy comparison metrics"},
        {"item": "S5", "file": "proxy_compare_samples.tsv", "role": "proxy comparison sample-level table"},
        {"item": "S6", "file": "NEOICI_NATURE_CANCER_SAMPLE_TABLE.tsv", "role": "sample ranking and proxy ranks"},
        {"item": "S7", "file": "NEOICI_NATURE_CANCER_LOO_TABLE.tsv", "role": "leave-one-out robustness"},
        {"item": "S8", "file": "NEOICI_NATURE_CANCER_CLAIM_TABLE.tsv", "role": "claim ladder"},
        {"item": "S9", "file": "NEOICI_NATURE_CANCER_SYNTHESIS_TABLE.tsv", "role": "top-level synthesis contract"},
        {"item": "S10", "file": "CTMS1_FOLLOWUP4_CASE_STUDY_SUMMARY.md", "role": "case-study summary"},
    ]
    supp_df = pd.DataFrame(supp_rows)
    supp_df.to_csv(OUT / "NEOICI_NATURE_CANCER_SUPPLEMENTARY_MAP.tsv", sep="\t", index=False)

    # Main tables used in the dossier.
    top_local = gauntlet.sort_values(["split", "AUPRC"], ascending=[True, False])
    academic = top_local[top_local["split"] == "academic_validation_like"].head(6).copy()
    industrial = top_local[top_local["split"] == "industrial_locked_v0"].head(6).copy()
    academic = academic[
        ["algorithm", "AUPRC", "AUROC", "top10_precision", "top24_precision", "claim_disposition"]
    ]
    industrial = industrial[
        ["algorithm", "AUPRC", "AUROC", "top10_precision", "top24_precision", "claim_disposition"]
    ]

    combo_head = proxy_samples.sort_values("response_index", ascending=False).head(9).copy()
    for col in [
        "response_index",
        "top_clone_frac",
        "cytotoxic_score",
        "activation_score",
        "exhaustion_score",
        "memory_score",
        "NeoPrecis_style_proxy",
        "NeoICI_style_proxy",
    ]:
        combo_head[col] = combo_head[col].map(lambda v: format_num(v, 3))

    sample_small = sample[
        [
            "dataset",
            "sample",
            "response_index",
            "rank_response",
            "NeoICI_style_proxy",
            "rank_neoici",
            "NeoPrecis_style_proxy",
            "rank_neoprecis",
            "top_clone_frac",
        ]
    ].copy()
    sample_small["response_index"] = sample_small["response_index"].map(lambda v: format_num(v, 3))
    sample_small["NeoICI_style_proxy"] = sample_small["NeoICI_style_proxy"].map(lambda v: format_num(v, 3))
    sample_small["NeoPrecis_style_proxy"] = sample_small["NeoPrecis_style_proxy"].map(lambda v: format_num(v, 3))
    sample_small["top_clone_frac"] = sample_small["top_clone_frac"].map(lambda v: format_num(v, 3))
    failure_small = sample[
        [
            "dataset",
            "sample",
            "response_index",
            "rank_response",
            "NeoICI_style_proxy",
            "rank_neoici",
            "NeoPrecis_style_proxy",
            "rank_neoprecis",
            "top_clone_frac",
        ]
    ].sort_values("response_index", ascending=True).head(3).copy()
    failure_small["response_index"] = failure_small["response_index"].map(lambda v: format_num(v, 3))
    failure_small["NeoICI_style_proxy"] = failure_small["NeoICI_style_proxy"].map(lambda v: format_num(v, 3))
    failure_small["NeoPrecis_style_proxy"] = failure_small["NeoPrecis_style_proxy"].map(lambda v: format_num(v, 3))
    failure_small["top_clone_frac"] = failure_small["top_clone_frac"].map(lambda v: format_num(v, 3))

    loo_small = loo.copy()
    loo_small["loo_min"] = loo_small["loo_min"].map(lambda v: format_num(v, 3))
    loo_small["loo_median"] = loo_small["loo_median"].map(lambda v: format_num(v, 3))
    loo_small["loo_max"] = loo_small["loo_max"].map(lambda v: format_num(v, 3))
    loo_small["loo_sd"] = loo_small["loo_sd"].map(lambda v: format_num(v, 3))

    claim_small = claim.copy()
    claim_small.columns = ["Status", "Claim", "Evidence"]

    executive_small = executive.copy()

    fig_dir = OUT / "figures"
    ensure_dir(fig_dir)

    # Figure A: ranking bars for the validation-like split.
    fig, ax = plt.subplots(figsize=(14.4, 4.8))
    ax.axis("off")
    boxes = [
        (0.03, 0.42, 0.21, 0.36, "#fff1ef", "Step 1", "Local benchmark", "NeoICI가 top-tier인지 먼저 확인"),
        (0.27, 0.42, 0.21, 0.36, "#eef5fb", "Step 2", "Open replay", "공개 vaccine + ICI cohort에서 response-like signal 확인"),
        (0.51, 0.42, 0.21, 0.36, "#fff5df", "Step 3", "Proxy comparison", "NeoICI-style vs NeoPrecis-style 직접 대결"),
        (0.75, 0.42, 0.21, 0.36, "#edf5ea", "Step 4", "Case study", "CTMS1_followup4가 왜 top lane인지 설명"),
    ]
    for x0, y0, w, h, color, step, title, body in boxes:
        ax.add_patch(plt.Rectangle((x0, y0), w, h, facecolor=color, edgecolor="#17212f", linewidth=1.7))
        ax.text(x0 + 0.02, y0 + h - 0.07, step, fontsize=12, fontweight="bold", color="#8f2d25")
        ax.text(x0 + 0.02, y0 + h - 0.16, title, fontsize=18, fontweight="bold", color="#17212f")
        ax.text(x0 + 0.02, y0 + 0.08, body, fontsize=11.2, color="#34495c", linespacing=1.35)
    for x in [0.24, 0.48, 0.72]:
        ax.annotate("", xy=(x + 0.025, 0.60), xytext=(x - 0.012, 0.60), arrowprops=dict(arrowstyle="->", lw=2.5, color="#8f2d25"))
    ax.text(0.5, 0.88, "NeoICI Nature Cancer package / OURS evidence stack", ha="center", va="center", fontsize=17, fontweight="bold", color="#17212f")
    ax.text(0.5, 0.12, "Conclusion: readiness layer, not universal scorer", ha="center", va="center", fontsize=13.5, color="#8f2d25", fontweight="bold")
    ax.text(0.02, 0.96, "★ OURS", transform=ax.transAxes, fontsize=12, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_NeoICI_evidence_stack.png", fig)

    # Figure A0B: verdict poster for the hero area.
    fig, ax = plt.subplots(figsize=(15.6, 5.6))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.03, 0.10), 0.94, 0.80, facecolor="#fffaf1", edgecolor="#17212f", linewidth=2.0))
    ax.add_patch(plt.Rectangle((0.05, 0.16), 0.28, 0.64, facecolor="#fff1ef", edgecolor="#8f2d25", linewidth=1.8))
    ax.add_patch(plt.Rectangle((0.36, 0.16), 0.28, 0.64, facecolor="#eef5fb", edgecolor="#244e73", linewidth=1.8))
    ax.add_patch(plt.Rectangle((0.67, 0.16), 0.26, 0.64, facecolor="#eef5ea", edgecolor="#426b50", linewidth=1.8))
    ax.text(0.5, 0.83, "NeoICI Nature Cancer package", ha="center", va="center", fontsize=22, fontweight="bold", color="#17212f")
    ax.text(0.5, 0.74, "READINESS LAYER  |  NOT UNIVERSAL SCORER", ha="center", va="center", fontsize=18, fontweight="bold", color="#8f2d25")
    ax.text(0.19, 0.68, "WHAT WE SHOW", ha="center", va="center", fontsize=14, fontweight="bold", color="#8f2d25")
    ax.text(0.50, 0.68, "WHY IT MATTERS", ha="center", va="center", fontsize=14, fontweight="bold", color="#244e73")
    ax.text(0.80, 0.68, "WHAT WE DO NOT CLAIM", ha="center", va="center", fontsize=14, fontweight="bold", color="#426b50")
    ax.text(0.19, 0.50, "KG_GA_evolved still leads\nvalidation-like AUPRC = 0.978\nNeoICI remains competitive", ha="center", va="center", fontsize=14, color="#17212f", linespacing=1.35)
    ax.text(0.50, 0.50, "open combo replay shows\nNeoICI-style proxy rho = 0.817\nCTMS1_followup4 is rank 1/9", ha="center", va="center", fontsize=14, color="#17212f", linespacing=1.35)
    ax.text(0.80, 0.50, "No universal superiority\nNo clinical selection claim\nNo inflated mechanistic claim", ha="center", va="center", fontsize=14, color="#17212f", linespacing=1.35)
    ax.text(0.50, 0.26, "This package is strongest when read as: benchmark + replay + proxy comparison + case study", ha="center", va="center", fontsize=15, fontweight="bold", color="#17212f")
    ax.text(0.50, 0.18, "★ OURS", ha="center", va="center", fontsize=13, fontweight="bold", color="#8f2d25")
    save_fig(fig_dir / "Fig_NeoICI_verdict_poster.png", fig)

    # Figure A0C: claim poster for the first screen.
    fig, ax = plt.subplots(figsize=(16.0, 6.2))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.02, 0.08), 0.96, 0.84, facecolor="#0f1824", edgecolor="#17212f", linewidth=2.2))
    ax.add_patch(plt.Rectangle((0.04, 0.12), 0.92, 0.76, facecolor="#17212f", edgecolor="#7b1f2a", linewidth=1.6))
    ax.text(0.06, 0.82, "OURS / NeoICI Nature Cancer package", fontsize=14, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono")
    ax.text(0.06, 0.66, "NeoICI is a readiness layer for combo therapy, not a universal scorer.", fontsize=28, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond", linespacing=1.1, wrap=True)
    ax.text(0.06, 0.50, "Local benchmark honest. Public replay directional. Case study anchored.", fontsize=18, color="#f1ddd1", fontfamily="Newsreader")
    ax.text(0.06, 0.34, "0.978   |   0.817   |   1 / 9", fontsize=28, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.06, 0.23, "KG_GA_evolved / NeoICI-style rho / CTMS1_followup4 rank", fontsize=14, color="#d8e3ef", fontfamily="JetBrains Mono")
    for x, label, color in [(0.63, "KG_GA_evolved", "#244e73"), (0.77, "NeoICI replay", "#8f2d25"), (0.88, "CTMS1", "#b58534")]:
        ax.add_patch(plt.Rectangle((x, 0.62), 0.08, 0.16, facecolor=color, edgecolor="white", linewidth=1.0))
        ax.text(x + 0.04, 0.70, label, fontsize=11, color="white", ha="center", va="center", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.12, "claim boundary locked: readiness yes, universal claim no, clinical selector no", ha="center", va="center", fontsize=12.5, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_claim_poster.png", fig)

    # Figure A0D: first-screen collage.
    collage_sources = [
        ("Claim poster", fig_dir / "Fig_NeoICI_claim_poster.png"),
        ("Evidence stack", fig_dir / "Fig_NeoICI_evidence_stack.png"),
        ("Replay figure", ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "figures" / "Fig_combo_response_index.png"),
        ("Case study", ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "ctms1_followup4_case" / "figures" / "Fig_CTMS1_followup4_case_study.png"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(15.8, 11.2))
    axes = axes.flatten()
    for ax, (label, src) in zip(axes, collage_sources):
        ax.axis("off")
        img = plt.imread(src)
        ax.imshow(img)
        ax.set_title(label, fontsize=14, fontweight="bold", color="#17212f", pad=10)
    fig.suptitle("NeoICI front-page collage / OURS", fontsize=20, fontweight="bold", y=0.98, color="#17212f")
    fig.text(0.5, 0.015, "Claim poster + evidence stack + replay + case study. This is the page in one glance.", ha="center", va="bottom", fontsize=12.5, color="#8f2d25", fontweight="bold")
    save_fig(fig_dir / "Fig_NeoICI_frontpage_collage.png", fig)

    # Figure A0E: ultra-summary super poster.
    fig, ax = plt.subplots(figsize=(16.6, 7.0))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.015, 0.08), 0.97, 0.84, facecolor="#101722", edgecolor="#17212f", linewidth=2.4))
    ax.add_patch(plt.Rectangle((0.035, 0.12), 0.42, 0.76, facecolor="#fff1ef", edgecolor="#8f2d25", linewidth=1.8))
    ax.add_patch(plt.Rectangle((0.475, 0.12), 0.49, 0.76, facecolor="#17212f", edgecolor="#244e73", linewidth=1.8))
    ax.text(0.06, 0.82, "OURS / NeoICI Nature Cancer package", fontsize=15, fontweight="bold", color="#8f2d25", fontfamily="JetBrains Mono")
    ax.text(0.06, 0.67, "Readiness layer.", fontsize=34, fontweight="bold", color="#17212f", fontfamily="Cormorant Garamond")
    ax.text(0.06, 0.55, "Not universal scorer.", fontsize=30, fontweight="bold", color="#17212f", fontfamily="Cormorant Garamond")
    ax.text(0.06, 0.39, "0.978", fontsize=56, fontweight="bold", color="#8f2d25", fontfamily="JetBrains Mono")
    ax.text(0.06, 0.29, "KG_GA_evolved still leads the local benchmark.", fontsize=16, color="#17212f", fontfamily="Newsreader")
    ax.text(0.56, 0.79, "The whole package in three numbers", fontsize=20, fontweight="bold", color="#fff1ef", fontfamily="Cormorant Garamond")
    for y, num, txt, color in [
        (0.60, "0.817", "NeoICI-style proxy tracks response_index better in the open replay.", "#ffcf8a"),
        (0.42, "1 / 9", "CTMS1_followup4 is the clearest explanation lane.", "#d8e3ef"),
        (0.24, "SAFE", "Claim boundary stays bounded: readiness yes, universal no.", "#b7d7b0"),
    ]:
        ax.text(0.56, y, num, fontsize=42, fontweight="bold", color=color, fontfamily="JetBrains Mono")
        ax.text(0.69, y + 0.01, txt, fontsize=15, color="#f1ddd1", fontfamily="Newsreader", va="center")
    ax.text(0.50, 0.12, "★ OURS  /  benchmark + replay + proxy comparison + case study", ha="center", va="center", fontsize=13.5, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_super_poster.png", fig)

    # Figure A0F: dark monument poster.
    fig, ax = plt.subplots(figsize=(17.2, 7.6))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.015, 0.06), 0.97, 0.86, facecolor="#0b1220", edgecolor="#17212f", linewidth=2.4))
    ax.add_patch(plt.Rectangle((0.045, 0.11), 0.91, 0.16, facecolor="#17212f", edgecolor="#7b1f2a", linewidth=1.6))
    ax.text(0.06, 0.20, "OURS / NeoICI Nature Cancer package", fontsize=15, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono", va="center")
    ax.text(0.50, 0.72, "READINESS LAYER", fontsize=60, fontweight="bold", color="#ffffff", ha="center", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.59, "NOT UNIVERSAL SCORER", fontsize=42, fontweight="bold", color="#ffcf8a", ha="center", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.44, "benchmark + replay + proxy comparison + case study", fontsize=22, color="#d8e3ef", ha="center", fontfamily="JetBrains Mono")
    ax.add_patch(plt.Rectangle((0.08, 0.30), 0.84, 0.06, facecolor="#8f2d25", edgecolor="none"))
    ax.text(0.08, 0.24, "0.978      0.817      1 / 9", fontsize=34, fontweight="bold", color="#ffffff", fontfamily="JetBrains Mono")
    ax.text(0.08, 0.16, "KG_GA_evolved      NeoICI replay      CTMS1_followup4", fontsize=16, color="#f1ddd1", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.08, "Claim boundary locked: readiness yes, universal claim no, clinical selector no", ha="center", va="center", fontsize=13, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_monument_poster.png", fig)

    # Figure A0G: true first-screen cover wall.
    fig, ax = plt.subplots(figsize=(18.4, 8.4))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.05), 0.98, 0.90, facecolor="#fffaf1", edgecolor="#17212f", linewidth=2.4))
    ax.add_patch(plt.Rectangle((0.03, 0.09), 0.38, 0.82, facecolor="#17212f", edgecolor="#8f2d25", linewidth=1.8))
    ax.add_patch(plt.Rectangle((0.43, 0.09), 0.54, 0.82, facecolor="#0f1824", edgecolor="#244e73", linewidth=1.8))
    ax.text(0.05, 0.84, "OURS / NeoICI Nature Cancer package", fontsize=16, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono")
    ax.text(0.05, 0.68, "COVER WALL", fontsize=28, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.05, 0.54, "Readiness layer.", fontsize=42, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond")
    ax.text(0.05, 0.41, "Not universal scorer.", fontsize=38, fontweight="bold", color="#f1ddd1", fontfamily="Cormorant Garamond")
    ax.text(0.47, 0.83, "Local benchmark is honest. Open replay is directional. Case study is anchored.", fontsize=22, fontweight="bold", color="#fff1ef", fontfamily="Cormorant Garamond")
    ax.text(0.47, 0.68, "0.978", fontsize=72, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.56, 0.69, "KG_GA_evolved", fontsize=20, color="#d8e3ef", fontfamily="JetBrains Mono")
    ax.text(0.47, 0.54, "0.817", fontsize=72, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.56, 0.55, "NeoICI-style proxy / response_index", fontsize=20, color="#d8e3ef", fontfamily="JetBrains Mono")
    ax.text(0.47, 0.40, "1 / 9", fontsize=72, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.56, 0.41, "CTMS1_followup4 anchor case", fontsize=20, color="#d8e3ef", fontfamily="JetBrains Mono")
    ax.add_patch(plt.Rectangle((0.47, 0.22), 0.45, 0.07, facecolor="#8f2d25", edgecolor="none"))
    ax.text(0.69, 0.255, "CLAIM BOUNDARY LOCKED", ha="center", va="center", fontsize=16, color="#ffffff", fontweight="bold", fontfamily="JetBrains Mono")
    ax.text(0.69, 0.16, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=15, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_cover_wall.png", fig)

    # Figure A0O: impact banner that sits above the cover wall.
    fig, ax = plt.subplots(figsize=(20.8, 6.6))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.012, 0.07), 0.976, 0.86, facecolor="#08111d", edgecolor="#8f2d25", linewidth=3.4))
    ax.text(0.49, 0.47, "OURS", ha="center", va="center", fontsize=220, fontweight="bold",
            color="#ffffff", alpha=0.07, fontfamily="JetBrains Mono", rotation=-8,
            path_effects=[pe.withStroke(linewidth=6, foreground="#ffffff", alpha=0.04)])
    ax.add_patch(plt.Rectangle((0.03, 0.79), 0.94, 0.09, facecolor="#8f2d25", edgecolor="none"))
    ax.text(0.05, 0.835, "★ OURS / NEOICI NATURE CANCER", fontsize=19, fontweight="bold",
            color="#fff1ef", fontfamily="JetBrains Mono", va="center",
            path_effects=[pe.withStroke(linewidth=3.5, foreground="#3a0b0a")])

    ax.add_patch(plt.Rectangle((0.04, 0.18), 0.24, 0.54, facecolor="#17212f", edgecolor="#244e73", linewidth=1.8))
    ax.text(0.065, 0.58, "READINESS", fontsize=30, fontweight="bold", color="#ffcf8a",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=4, foreground="#111a27")])
    ax.text(0.065, 0.44, "LAYER", fontsize=64, fontweight="bold", color="#ffffff",
            fontfamily="Cormorant Garamond", path_effects=[pe.withStroke(linewidth=5, foreground="#111a27")])
    ax.text(0.065, 0.28, "combo-therapy", fontsize=24, fontweight="bold", color="#d8e3ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=3, foreground="#111a27")])
    ax.text(0.065, 0.20, "not universal scorer", fontsize=20, fontweight="bold", color="#fff1ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=3, foreground="#111a27")])

    ax.add_patch(plt.Rectangle((0.31, 0.19), 0.42, 0.53, facecolor="#0f1824", edgecolor="#244e73", linewidth=1.8))
    ax.add_patch(plt.Rectangle((0.33, 0.56), 0.38, 0.10, facecolor="#17212f", edgecolor="none"))
    ax.text(0.52, 0.61, "BENCHMARK TRIAD", ha="center", va="center", fontsize=18, fontweight="bold",
            color="#fff1ef", fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=3, foreground="#111a27")])
    ax.text(0.36, 0.44, "0.978", ha="center", va="center", fontsize=98, fontweight="bold", color="#ffcf8a",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=5, foreground="#111a27")])
    ax.text(0.52, 0.44, "0.817", ha="center", va="center", fontsize=98, fontweight="bold", color="#ffcf8a",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=5, foreground="#111a27")])
    ax.text(0.68, 0.44, "1 / 9", ha="center", va="center", fontsize=98, fontweight="bold", color="#ffcf8a",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=5, foreground="#111a27")])
    ax.text(0.36, 0.24, "KG_GA_evolved", ha="center", va="center", fontsize=18, color="#d8e3ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=3, foreground="#111a27")])
    ax.text(0.52, 0.24, "replay rho", ha="center", va="center", fontsize=18, color="#d8e3ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=3, foreground="#111a27")])
    ax.text(0.68, 0.24, "CTMS1_followup4", ha="center", va="center", fontsize=18, color="#d8e3ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=3, foreground="#111a27")])

    ax.add_patch(plt.Rectangle((0.75, 0.18), 0.20, 0.54, facecolor="#8f2d25", edgecolor="#fff1ef", linewidth=1.6))
    ax.text(0.85, 0.57, "CLAIM", ha="center", va="center", fontsize=26, fontweight="bold", color="#fff1ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=4, foreground="#4d1412")])
    ax.text(0.85, 0.42, "BOUNDARY", ha="center", va="center", fontsize=26, fontweight="bold", color="#fff1ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=4, foreground="#4d1412")])
    ax.text(0.85, 0.28, "LOCKED", ha="center", va="center", fontsize=28, fontweight="bold", color="#ffcf8a",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=4, foreground="#4d1412")])
    ax.text(0.85, 0.20, "readiness yes / universal no", ha="center", va="center", fontsize=13, color="#fff1ef",
            fontfamily="JetBrains Mono", path_effects=[pe.withStroke(linewidth=3, foreground="#4d1412")])

    ax.add_patch(plt.Rectangle((0.04, 0.10), 0.91, 0.06, facecolor="#17212f", edgecolor="none"))
    ax.text(0.50, 0.13, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=17, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono",
            path_effects=[pe.withStroke(linewidth=3, foreground="#111a27")])
    save_fig(fig_dir / "Fig_NeoICI_impact_banner.png", fig)

    # Figure A0J: cover slab for maximum first-screen impact.
    fig, ax = plt.subplots(figsize=(19.0, 9.0))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.04), 0.98, 0.92, facecolor="#fffaf1", edgecolor="#17212f", linewidth=2.8))
    ax.add_patch(plt.Rectangle((0.03, 0.08), 0.94, 0.20, facecolor="#17212f", edgecolor="#8f2d25", linewidth=1.8))
    ax.text(0.05, 0.18, "★ OURS / NEOICI / NATURE CANCER COVER SLAB", fontsize=18, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono", va="center")
    ax.text(0.50, 0.78, "NeoICI is a", ha="center", va="center", fontsize=34, fontweight="bold", color="#17212f", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.67, "combo-therapy readiness layer", ha="center", va="center", fontsize=46, fontweight="bold", color="#8f2d25", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.56, "not a universal scorer", ha="center", va="center", fontsize=36, fontweight="bold", color="#17212f", fontfamily="Cormorant Garamond")
    ax.add_patch(plt.Rectangle((0.08, 0.34), 0.84, 0.10, facecolor="#8f2d25", edgecolor="none"))
    ax.text(0.50, 0.39, "0.978        0.817        1 / 9", ha="center", va="center", fontsize=44, fontweight="bold", color="#ffffff", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.27, "KG_GA_evolved        replay rho        CTMS1_followup4", ha="center", va="center", fontsize=18, color="#244e73", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.14, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=18, color="#8f2d25", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_cover_slab.png", fig)

    # Figure A0K: launch banner for the first visual hit.
    fig, ax = plt.subplots(figsize=(19.6, 5.8))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.08), 0.98, 0.84, facecolor="#0b1220", edgecolor="#17212f", linewidth=2.8))
    ax.add_patch(plt.Rectangle((0.03, 0.12), 0.94, 0.18, facecolor="#8f2d25", edgecolor="#fff1ef", linewidth=1.4))
    ax.text(0.05, 0.21, "★ OURS / LAUNCH BANNER", fontsize=18, fontweight="bold", color="#ffffff", fontfamily="JetBrains Mono", va="center")
    ax.text(0.07, 0.66, "NEOICI", fontsize=68, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond")
    ax.text(0.29, 0.68, "combo-therapy readiness layer", fontsize=30, fontweight="bold", color="#ffcf8a", fontfamily="Cormorant Garamond")
    ax.text(0.07, 0.42, "0.978", fontsize=62, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.24, 0.45, "local top", fontsize=18, color="#d8e3ef", fontfamily="JetBrains Mono")
    ax.text(0.43, 0.42, "0.817", fontsize=62, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.60, 0.45, "replay rho", fontsize=18, color="#d8e3ef", fontfamily="JetBrains Mono")
    ax.text(0.78, 0.42, "1 / 9", fontsize=62, fontweight="bold", color="#ffcf8a", fontfamily="JetBrains Mono")
    ax.text(0.88, 0.45, "anchor", fontsize=18, color="#d8e3ef", fontfamily="JetBrains Mono", ha="right")
    ax.text(0.50, 0.18, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=18, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_launch_banner.png", fig)

    # Figure A0L: topline banner for the immediate opening hit.
    fig, ax = plt.subplots(figsize=(20.0, 4.8))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.10), 0.98, 0.80, facecolor="#17212f", edgecolor="#8f2d25", linewidth=2.4))
    ax.text(0.50, 0.68, "NEOICI", ha="center", va="center", fontsize=62, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.48, "combo-therapy readiness layer", ha="center", va="center", fontsize=32, fontweight="bold", color="#ffcf8a", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.28, "0.978   |   0.817   |   1 / 9   |   OURS", ha="center", va="center", fontsize=24, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.15, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=16, color="#d8e3ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_topline_banner.png", fig)

    # Figure A0M: claim strip for the very first visual hit.
    fig, ax = plt.subplots(figsize=(20.4, 3.4))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.18), 0.98, 0.64, facecolor="#8f2d25", edgecolor="#17212f", linewidth=2.2))
    ax.text(0.04, 0.50, "★ OURS", ha="left", va="center", fontsize=20, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono")
    ax.text(0.27, 0.50, "NEOICI", ha="center", va="center", fontsize=54, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond")
    ax.text(0.56, 0.50, "readiness layer", ha="center", va="center", fontsize=34, fontweight="bold", color="#ffcf8a", fontfamily="Cormorant Garamond")
    ax.text(0.84, 0.50, "0.978 | 0.817 | 1 / 9", ha="center", va="center", fontsize=22, fontweight="bold", color="#ffffff", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.23, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=16, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_claim_strip.png", fig)

    # Figure A0H: opening impact board.
    fig, ax = plt.subplots(figsize=(18.0, 7.2))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.012, 0.06), 0.976, 0.88, facecolor="#0f1824", edgecolor="#17212f", linewidth=2.4))
    ax.add_patch(plt.Rectangle((0.035, 0.10), 0.93, 0.18, facecolor="#17212f", edgecolor="#8f2d25", linewidth=1.6))
    ax.text(0.05, 0.20, "★ OURS  /  OPENING IMPACT BOARD", fontsize=16, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono", va="center")
    ax.text(0.50, 0.78, "NeoICI is a combo-therapy readiness layer.", ha="center", va="center", fontsize=40, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.66, "That is the only claim this page is trying to prove.", ha="center", va="center", fontsize=24, fontweight="bold", color="#ffcf8a", fontfamily="Cormorant Garamond")
    ax.add_patch(plt.Rectangle((0.06, 0.36), 0.25, 0.18, facecolor="#fff1ef", edgecolor="#8f2d25", linewidth=1.8))
    ax.add_patch(plt.Rectangle((0.375, 0.36), 0.25, 0.18, facecolor="#eef5fb", edgecolor="#244e73", linewidth=1.8))
    ax.add_patch(plt.Rectangle((0.69, 0.36), 0.25, 0.18, facecolor="#eef5ea", edgecolor="#426b50", linewidth=1.8))
    ax.text(0.185, 0.48, "0.978", ha="center", va="center", fontsize=48, fontweight="bold", color="#8f2d25", fontfamily="JetBrains Mono")
    ax.text(0.185, 0.41, "local top", ha="center", va="center", fontsize=14, color="#17212f", fontfamily="JetBrains Mono")
    ax.text(0.500, 0.48, "0.817", ha="center", va="center", fontsize=48, fontweight="bold", color="#244e73", fontfamily="JetBrains Mono")
    ax.text(0.500, 0.41, "open replay rho", ha="center", va="center", fontsize=14, color="#17212f", fontfamily="JetBrains Mono")
    ax.text(0.815, 0.48, "1 / 9", ha="center", va="center", fontsize=48, fontweight="bold", color="#426b50", fontfamily="JetBrains Mono")
    ax.text(0.815, 0.41, "anchor case", ha="center", va="center", fontsize=14, color="#17212f", fontfamily="JetBrains Mono")
    ax.add_patch(plt.Rectangle((0.06, 0.20), 0.88, 0.08, facecolor="#8f2d25", edgecolor="none"))
    ax.text(0.50, 0.24, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=16, color="#ffffff", fontweight="bold", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.12, "claim boundary locked", ha="center", va="center", fontsize=14, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_opening_impact_board.png", fig)

    # Figure A0N: micro strip for the first possible visual hit.
    fig, ax = plt.subplots(figsize=(20.8, 2.9))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.16), 0.98, 0.68, facecolor="#17212f", edgecolor="#8f2d25", linewidth=2.4))
    ax.text(0.03, 0.50, "★ OURS", ha="left", va="center", fontsize=18, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono")
    ax.text(0.19, 0.50, "NEOICI", ha="center", va="center", fontsize=50, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond")
    ax.text(0.45, 0.50, "readiness layer", ha="center", va="center", fontsize=30, fontweight="bold", color="#ffcf8a", fontfamily="Cormorant Garamond")
    ax.text(0.77, 0.50, "0.978 | 0.817 | 1 / 9", ha="center", va="center", fontsize=22, fontweight="bold", color="#ffffff", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.23, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=15, color="#d8e3ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_micro_strip.png", fig)

    # Figure A0I: verdict wall for the very first screen.
    fig, ax = plt.subplots(figsize=(18.8, 8.8))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.05), 0.98, 0.90, facecolor="#08111d", edgecolor="#17212f", linewidth=2.8))
    ax.add_patch(plt.Rectangle((0.03, 0.08), 0.94, 0.18, facecolor="#17212f", edgecolor="#8f2d25", linewidth=1.8))
    ax.text(0.05, 0.18, "★ OURS / FIRST-SCREEN VERDICT WALL", fontsize=18, fontweight="bold", color="#fff1ef", fontfamily="JetBrains Mono", va="center")
    ax.text(0.50, 0.77, "NeoICI", ha="center", va="center", fontsize=66, fontweight="bold", color="#ffffff", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.65, "combination readiness layer", ha="center", va="center", fontsize=34, fontweight="bold", color="#ffcf8a", fontfamily="Cormorant Garamond")
    ax.text(0.50, 0.53, "not a universal scorer", ha="center", va="center", fontsize=30, fontweight="bold", color="#f1ddd1", fontfamily="Cormorant Garamond")
    ax.add_patch(plt.Rectangle((0.07, 0.31), 0.86, 0.12, facecolor="#8f2d25", edgecolor="none"))
    ax.text(0.50, 0.37, "0.978   |   0.817   |   1 / 9", ha="center", va="center", fontsize=40, fontweight="bold", color="#ffffff", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.24, "KG_GA_evolved / replay rho / CTMS1_followup4", ha="center", va="center", fontsize=18, color="#d8e3ef", fontfamily="JetBrains Mono")
    ax.text(0.50, 0.12, "readiness yes   |   universal no   |   clinical selector no", ha="center", va="center", fontsize=16, color="#fff1ef", fontweight="bold", fontfamily="JetBrains Mono")
    save_fig(fig_dir / "Fig_NeoICI_verdict_wall.png", fig)

    # Figure A: ranking bars for the validation-like split.
    val = gauntlet[gauntlet["split"] == "academic_validation_like"].copy()
    val = val.sort_values("AUPRC", ascending=True).tail(10)
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    colors = ["#b58534" if a == "NeoICI_GA_RL_combo_v1" else "#244e73" if a == "KG_GA_evolved" else "#d9c9b5" for a in val["algorithm"]]
    ax.barh(val["algorithm"], val["AUPRC"], color=colors, edgecolor="#17212f", linewidth=0.7)
    ax.set_xlim(0.25, 1.0)
    ax.set_xlabel("AUPRC")
    ax.set_title("NeoICI gauntlet: validation-like split ranking")
    ax.grid(axis="x", alpha=0.18)
    for i, v in enumerate(val["AUPRC"]):
        ax.text(v + 0.01, i, f"{v:.3f}", va="center", ha="left", fontsize=9)
    save_fig(fig_dir / "Fig_NeoICI_algorithm_rankbars.png", fig)

    # Figure A2: algorithm metric heatmap for the top validation-like algorithms.
    algo_heat = val.sort_values("AUPRC", ascending=False).head(8).set_index("algorithm")[["AUPRC", "AUROC", "top10_precision", "top24_precision"]]
    fig, ax = plt.subplots(figsize=(8.9, 5.6))
    im = ax.imshow(algo_heat.to_numpy(), aspect="auto", cmap="YlOrRd", vmin=float(algo_heat.to_numpy().min()), vmax=float(algo_heat.to_numpy().max()))
    ax.set_xticks(range(algo_heat.shape[1]), algo_heat.columns, rotation=18, ha="right")
    ax.set_yticks(range(algo_heat.shape[0]), algo_heat.index)
    ax.set_title("Algorithm metric heatmap")
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i in range(algo_heat.shape[0]):
        for j in range(algo_heat.shape[1]):
            ax.text(j, i, f"{algo_heat.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8, color="#17212f")
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, label="metric value")
    save_fig(fig_dir / "Fig_NeoICI_algorithm_heatmap.png", fig)

    # Figure A3: win/loss heatmap across the comparator set.
    wl = winloss.copy()
    wl = wl[wl["champion"] == "NeoICI_GA_RL_combo_v1"].copy()
    wl = wl.sort_values("win_fraction", ascending=False)
    wl["row"] = wl["comparator"].str.replace("_", " ", regex=False)
    win_cols = [
        "win_fraction",
        "academic_validation_like_AUPRC_delta",
        "industrial_locked_v0_AUPRC_delta",
        "academic_low_medium_leakage_AUPRC_delta",
    ]
    wl_mat = wl.set_index("row")[win_cols]
    fig, ax = plt.subplots(figsize=(9.6, 6.0))
    im = ax.imshow(wl_mat.to_numpy(), aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(wl_mat.shape[1]), ["win%", "acad Δ", "ind Δ", "leakage Δ"], rotation=18, ha="right")
    ax.set_yticks(range(wl_mat.shape[0]), wl_mat.index)
    ax.set_title("NeoICI win/loss summary heatmap")
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i in range(wl_mat.shape[0]):
        for j in range(wl_mat.shape[1]):
            ax.text(j, i, f"{wl_mat.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7.5, color="#17212f")
    fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02, label="value")
    save_fig(fig_dir / "Fig_NeoICI_winloss_heatmap.png", fig)

    # Figure A4: comparator win-fraction bars.
    win_bars = wl.sort_values("win_fraction", ascending=True)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    bar_colors = ["#8f2d25" if "NeoICI" in c else "#244e73" if "KG" in c else "#b58534" if "BioDarwin" in c else "#d9c9b5" for c in win_bars["comparator"]]
    ax.barh(win_bars["row"], win_bars["win_fraction"], color=bar_colors, edgecolor="#17212f", linewidth=0.7)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("win fraction")
    ax.set_title("NeoICI comparator win fractions")
    ax.grid(axis="x", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i, v in enumerate(win_bars["win_fraction"]):
        ax.text(v + 0.015, i, f"{v:.2f}", va="center", ha="left", fontsize=8.5)
    save_fig(fig_dir / "Fig_NeoICI_comparator_winbars.png", fig)

    # Figure A4: algorithm family average bars.
    fam = val.groupby("algorithm_family", as_index=False)["AUPRC"].mean().sort_values("AUPRC", ascending=False).head(8)
    fig, ax = plt.subplots(figsize=(8.7, 4.9))
    bar_colors = ["#8f2d25" if "NeoICI" in f else "#244e73" if "KG" in f else "#b58534" if "BioDarwin" in f else "#d9c9b5" for f in fam["algorithm_family"]]
    ax.barh(fam["algorithm_family"], fam["AUPRC"], color=bar_colors, edgecolor="#17212f", linewidth=0.7)
    ax.set_xlim(0.25, 1.0)
    ax.set_xlabel("mean AUPRC")
    ax.set_title("Algorithm family mean AUPRC")
    ax.grid(axis="x", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i, v in enumerate(fam["AUPRC"]):
        ax.text(v + 0.01, i, f"{v:.3f}", va="center", ha="left", fontsize=9)
    save_fig(fig_dir / "Fig_NeoICI_algorithm_family_bars.png", fig)

    # Figure A5: algorithm score delta bars relative to KG_GA_evolved.
    ref = val[val["algorithm"] == "KG_GA_evolved"].iloc[0]
    delta_src = val.sort_values("AUPRC", ascending=False).head(8).copy()
    delta_src["delta_AUPRC"] = delta_src["AUPRC"] - ref["AUPRC"]
    delta_src["delta_AUROC"] = delta_src["AUROC"] - ref["AUROC"]
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    delta_src = delta_src.sort_values("delta_AUPRC", ascending=True)
    ax.barh(delta_src["algorithm"], delta_src["delta_AUPRC"], color=["#8f2d25" if a == "NeoICI_GA_RL_combo_v1" else "#d9c9b5" for a in delta_src["algorithm"]], edgecolor="#17212f", linewidth=0.7)
    ax.axvline(0, color="#17212f", linewidth=1)
    ax.set_xlabel("AUPRC delta vs KG_GA_evolved")
    ax.set_title("Algorithm delta bars")
    ax.grid(axis="x", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i, v in enumerate(delta_src["delta_AUPRC"]):
        ax.text(v + (0.005 if v >= 0 else -0.005), i, f"{v:+.3f}", va="center", ha="left" if v >= 0 else "right", fontsize=8.5)
    save_fig(fig_dir / "Fig_NeoICI_algorithm_delta_bars.png", fig)

    # Figure A6: top-k precision profile for the strongest algorithms.
    precision_src = val.sort_values("AUPRC", ascending=False).head(8)[["algorithm", "top5_precision", "top10_precision", "top24_precision"]].copy()
    precision_src = precision_src.set_index("algorithm")
    fig, ax = plt.subplots(figsize=(9.6, 5.7))
    x = np.arange(precision_src.shape[0])
    width = 0.24
    ax.bar(x - width, precision_src["top5_precision"], width=width, color="#8f2d25", label="top5")
    ax.bar(x, precision_src["top10_precision"], width=width, color="#b58534", label="top10")
    ax.bar(x + width, precision_src["top24_precision"], width=width, color="#244e73", label="top24")
    ax.set_xticks(x, [a.replace("_", " ") for a in precision_src.index], rotation=35, ha="right")
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("precision")
    ax.set_title("Top-k precision profile")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.08))
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_NeoICI_precision_profile.png", fig)

    # Figure A7: cross-split generalization heatmap on selected algorithms.
    split_piv = gauntlet[gauntlet["split"].isin(["academic_validation_like", "industrial_locked_v0", "academic_low_medium_leakage"])].copy()
    split_piv = split_piv[split_piv["algorithm"].isin(["NeoICI_GA_RL_combo_v1", "KG_GA_evolved", "BioDarwin_mode_bank_v4", "CROSS_claimsafe", "BigMHC_IM"])].copy()
    split_piv["split"] = split_piv["split"].map({
        "academic_validation_like": "acad",
        "industrial_locked_v0": "ind",
        "academic_low_medium_leakage": "leak",
    })
    split_mat = split_piv.pivot_table(index="algorithm", columns="split", values="AUPRC", aggfunc="mean")
    split_mat = split_mat.reindex(
        index=["NeoICI_GA_RL_combo_v1", "KG_GA_evolved", "BioDarwin_mode_bank_v4", "CROSS_claimsafe", "BigMHC_IM"],
        columns=["acad", "ind", "leak"],
    )
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    split_arr = split_mat.to_numpy(dtype=float)
    im = ax.imshow(split_arr, aspect="auto", cmap="YlOrRd", vmin=float(np.nanmin(split_arr)), vmax=float(np.nanmax(split_arr)))
    ax.set_xticks(range(split_mat.shape[1]), split_mat.columns, rotation=18, ha="right")
    ax.set_yticks(range(split_mat.shape[0]), [a.replace("_", " ") for a in split_mat.index])
    ax.set_title("Cross-split generalization heatmap")
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i in range(split_mat.shape[0]):
        for j in range(split_mat.shape[1]):
            ax.text(j, i, f"{split_mat.iloc[i, j]:.3f}", ha="center", va="center", fontsize=8, color="#17212f")
    fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02, label="AUPRC")
    save_fig(fig_dir / "Fig_NeoICI_split_generalization_heatmap.png", fig)

    # Figure B: response index vs proxies.
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), sharex=True, sharey=False)
    palette = np.where(proxy_samples["sample"] == "CTMS1_followup4", "#8f2d25", "#244e73")
    axes[0].scatter(proxy_samples["response_index"], proxy_samples["NeoICI_style_proxy"], s=70, c=palette, edgecolor="white", linewidth=0.7)
    axes[1].scatter(proxy_samples["response_index"], proxy_samples["NeoPrecis_style_proxy"], s=70, c=palette, edgecolor="white", linewidth=0.7)
    x = np.linspace(proxy_samples["response_index"].min() - 0.02, proxy_samples["response_index"].max() + 0.02, 100)
    for ax, col, title in [
        (axes[0], "NeoICI_style_proxy", "NeoICI-style proxy"),
        (axes[1], "NeoPrecis_style_proxy", "NeoPrecis-style proxy"),
    ]:
        coef = np.polyfit(proxy_samples["response_index"], proxy_samples[col], 1)
        ax.plot(x, coef[0] * x + coef[1], color="#8f2d25", linewidth=1.5)
        ax.set_title(title)
        ax.set_xlabel("response_index")
        ax.grid(alpha=0.18)
    axes[0].set_ylabel("proxy score")
    axes[0].text(0.02, 0.96, "★ OURS", transform=axes[0].transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    axes[1].text(0.02, 0.96, "★ OURS", transform=axes[1].transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    fig.suptitle("Open combo replay: same samples, two proxy families", y=1.02, fontsize=14)
    save_fig(fig_dir / "Fig_NeoICI_proxy_scatter.png", fig)

    # Figure B2: dataset-level comparison bars for the open replay cohorts.
    dataset_summary = proxy_samples.groupby("dataset", as_index=False)[["response_index", "top_clone_frac", "cytotoxic_score", "NeoICI_style_proxy"]].median()
    dataset_summary = dataset_summary.set_index("dataset")
    fig, ax = plt.subplots(figsize=(8.9, 5.4))
    metrics = ["response_index", "top_clone_frac", "cytotoxic_score", "NeoICI_style_proxy"]
    x = np.arange(len(metrics))
    width = 0.32
    ds1 = dataset_summary.loc["GSE222011", metrics].values
    ds2 = dataset_summary.loc["GSE255830", metrics].values
    ax.bar(x - width / 2, ds1, width=width, color="#8f2d25", label="GSE222011")
    ax.bar(x + width / 2, ds2, width=width, color="#244e73", label="GSE255830")
    ax.set_xticks(x, ["response", "clone", "cytotoxic", "NeoICI proxy"])
    ax.set_ylabel("median score")
    ax.set_title("Dataset-level replay comparison")
    ax.legend(frameon=False, ncol=2)
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_NeoICI_dataset_comparison_bars.png", fig)

    # Figure C: leave-one-out robustness.
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    loo_plot = loo.copy()
    order = ["NeoPrecis_style_proxy", "NeoICI_style_proxy"]
    labels = ["NeoPrecis-style", "NeoICI-style"]
    mins = loo_plot.set_index("proxy").loc[order, "loo_min"].values
    meds = loo_plot.set_index("proxy").loc[order, "loo_median"].values
    maxs = loo_plot.set_index("proxy").loc[order, "loo_max"].values
    xs = np.arange(len(order))
    ax.errorbar(xs, meds, yerr=[meds - mins, maxs - meds], fmt="o", color="#8f2d25", ecolor="#b58534", capsize=8, elinewidth=2, markersize=8)
    ax.set_xticks(xs, labels)
    ax.set_ylim(0.15, 0.98)
    ax.set_ylabel("Spearman (leave-one-out)")
    ax.set_title("Proxy robustness under leave-one-out")
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.96, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for x0, m in zip(xs, meds):
        ax.text(x0, m + 0.03, f"{m:.3f}", ha="center", va="bottom", fontsize=10)
    save_fig(fig_dir / "Fig_NeoICI_proxy_loo.png", fig)

    # Figure C2: proxy family summary bars.
    fig, ax = plt.subplots(figsize=(8.2, 4.9))
    proxy_summary = proxy.set_index("proxy").loc[["NeoICI_style_proxy", "NeoPrecis_style_proxy"], ["spearman_vs_response_index", "spearman_vs_top_clone_frac", "spearman_vs_cytotoxic_score"]]
    x = np.arange(proxy_summary.shape[1])
    width = 0.34
    ax.bar(x - width / 2, proxy_summary.loc["NeoICI_style_proxy"], width=width, color="#8f2d25", label="NeoICI-style")
    ax.bar(x + width / 2, proxy_summary.loc["NeoPrecis_style_proxy"], width=width, color="#244e73", label="NeoPrecis-style")
    ax.set_xticks(x, ["response", "clone", "cytotoxic"])
    ax.set_ylim(-0.25, 1.0)
    ax.set_ylabel("Spearman")
    ax.set_title("Proxy family summary")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_NeoICI_proxy_summary_bars.png", fig)

    # Figure C3: proxy correlation heatmap.
    proxy_corr = proxy_samples[["response_index", "top_clone_frac", "cytotoxic_score", "activation_score", "memory_score", "NeoICI_style_proxy", "NeoPrecis_style_proxy"]].corr(numeric_only=True)
    fig, ax = plt.subplots(figsize=(8.8, 6.1))
    im = ax.imshow(proxy_corr.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(proxy_corr.shape[1]), [c.replace("_", " ") for c in proxy_corr.columns], rotation=20, ha="right")
    ax.set_yticks(range(proxy_corr.shape[0]), [c.replace("_", " ") for c in proxy_corr.index])
    ax.set_title("Proxy and state correlation heatmap")
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i in range(proxy_corr.shape[0]):
        for j in range(proxy_corr.shape[1]):
            ax.text(j, i, f"{proxy_corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7.5, color="#17212f")
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, label="corr")
    save_fig(fig_dir / "Fig_NeoICI_proxy_corr_heatmap.png", fig)

    # Figure D: CTMS1_followup4 case profile vs cohort median.
    case_row = proxy_samples[proxy_samples["sample"] == "CTMS1_followup4"].iloc[0]
    metrics = ["response_index", "top_clone_frac", "cytotoxic_score", "activation_score", "memory_score"]
    case_vals = [case_row[m] for m in metrics]
    cohort_vals = [combo[m].median() for m in metrics]
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    idx = np.arange(len(metrics))
    width = 0.36
    ax.bar(idx - width/2, cohort_vals, width=width, color="#d9c9b5", edgecolor="#17212f", linewidth=0.6, label="cohort median")
    ax.bar(idx + width/2, case_vals, width=width, color="#8f2d25", edgecolor="#17212f", linewidth=0.6, label="CTMS1_followup4")
    ax.set_xticks(idx, ["response", "clone", "cytotoxic", "activation", "memory"], rotation=0)
    ax.set_ylabel("z-like score / fraction")
    ax.set_title("CTMS1_followup4: why the lane is strong")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.96, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_CTMS1_case_profile.png", fig)

    # Figure D7: case-study delta strip.
    fig, ax = plt.subplots(figsize=(8.9, 4.6))
    delta_vals = pd.Series({
        "response_index": case_row["response_index"] - cohort_vals[0],
        "top_clone_frac": case_row["top_clone_frac"] - cohort_vals[1],
        "cytotoxic_score": case_row["cytotoxic_score"] - cohort_vals[2],
        "activation_score": case_row["activation_score"] - cohort_vals[3],
        "memory_score": case_row["memory_score"] - cohort_vals[4],
    })
    colors = ["#8f2d25" if v > 0 else "#244e73" for v in delta_vals]
    ax.bar(delta_vals.index.str.replace("_", " ", regex=False), delta_vals.values, color=colors, edgecolor="#17212f", linewidth=0.6)
    ax.axhline(0, color="#17212f", linewidth=1)
    ax.set_title("CTMS1_followup4 minus cohort median")
    ax.set_ylabel("delta")
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_CTMS1_case_delta.png", fig)

    # Figure D2: replay sample heatmap.
    sample_heat = proxy_samples.set_index("sample")[["response_index", "top_clone_frac", "cytotoxic_score", "activation_score", "memory_score", "NeoICI_style_proxy", "NeoPrecis_style_proxy"]]
    fig, ax = plt.subplots(figsize=(10.2, 6.2))
    im = ax.imshow(sample_heat.to_numpy(), aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(sample_heat.shape[1]), ["response", "clone", "cytotoxic", "activation", "memory", "NeoICI", "NeoPrecis"], rotation=18, ha="right")
    ax.set_yticks(range(sample_heat.shape[0]), sample_heat.index)
    ax.set_title("Replay sample heatmap")
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02, label="raw score")
    save_fig(fig_dir / "Fig_NeoICI_sample_heatmap.png", fig)

    # Figure D3: sample rank profile across the 9 replay samples.
    rank_plot = sample[["sample", "rank_response", "rank_neoici", "rank_neoprecis"]].copy()
    rank_plot = rank_plot.sort_values("rank_response")
    fig, ax = plt.subplots(figsize=(10.8, 5.6))
    x = np.arange(len(rank_plot))
    ax.plot(x, rank_plot["rank_response"], marker="o", linewidth=2.2, color="#8f2d25", label="response rank")
    ax.plot(x, rank_plot["rank_neoici"], marker="o", linewidth=2.2, color="#244e73", label="NeoICI rank")
    ax.plot(x, rank_plot["rank_neoprecis"], marker="o", linewidth=2.2, color="#b58534", label="NeoPrecis rank")
    ax.set_xticks(x, rank_plot["sample"], rotation=45, ha="right")
    ax.invert_yaxis()
    ax.set_ylabel("rank (1 is best)")
    ax.set_title("Sample rank profile across the 9 replay samples")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.11))
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_NeoICI_sample_rank_profile.png", fig)

    # Figure D5: sample rank heatmap.
    rank_heat = sample.set_index("sample")[["rank_response", "rank_neoici", "rank_neoprecis"]]
    fig, ax = plt.subplots(figsize=(8.9, 5.6))
    im = ax.imshow(rank_heat.to_numpy(), aspect="auto", cmap="YlGnBu_r")
    ax.set_xticks(range(rank_heat.shape[1]), ["response", "NeoICI", "NeoPrecis"], rotation=18, ha="right")
    ax.set_yticks(range(rank_heat.shape[0]), rank_heat.index)
    ax.set_title("Sample rank heatmap")
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    for i in range(rank_heat.shape[0]):
        for j in range(rank_heat.shape[1]):
            ax.text(j, i, f"{rank_heat.iloc[i, j]:.0f}", ha="center", va="center", fontsize=8, color="#17212f")
    fig.colorbar(im, ax=ax, shrink=0.82, pad=0.02, label="rank")
    save_fig(fig_dir / "Fig_NeoICI_sample_rank_heatmap.png", fig)

    # Figure D6: sample rank mismatch bars.
    rank_mismatch = sample[["sample", "rank_response", "rank_neoici", "rank_neoprecis"]].copy()
    rank_mismatch["neoici_abs_diff"] = (rank_mismatch["rank_neoici"] - rank_mismatch["rank_response"]).abs()
    rank_mismatch["neoprecis_abs_diff"] = (rank_mismatch["rank_neoprecis"] - rank_mismatch["rank_response"]).abs()
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    ax.bar(rank_mismatch["sample"], rank_mismatch["neoici_abs_diff"], color="#8f2d25", label="NeoICI diff")
    ax.bar(rank_mismatch["sample"], rank_mismatch["neoprecis_abs_diff"], bottom=rank_mismatch["neoici_abs_diff"], color="#244e73", label="NeoPrecis diff")
    ax.set_ylabel("|rank - response rank|")
    ax.set_title("Sample rank mismatch bars")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.18)
    ax.legend(frameon=False, ncol=2)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_NeoICI_sample_rank_mismatch_bars.png", fig)

    # Figure D6b: failure lanes / negative controls.
    top_lanes = sample.sort_values("response_index", ascending=False).head(3).copy()
    low_lanes = sample.sort_values("response_index", ascending=True).head(3).copy()
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.2), sharey=True)
    lane_sets = [
        ("Top response lanes", top_lanes, "#8f2d25"),
        ("Low response lanes", low_lanes, "#244e73"),
    ]
    metric_cols = ["response_index", "NeoICI_style_proxy", "NeoPrecis_style_proxy"]
    metric_labels = ["response", "NeoICI", "NeoPrecis"]
    for ax, (title, lanes, color) in zip(axes, lane_sets):
        x = np.arange(len(lanes))
        width = 0.23
        for i, (col, lab) in enumerate(zip(metric_cols, metric_labels)):
            ax.bar(x + (i - 1) * width, lanes[col], width=width, color=[color, "#b58534", "#d9c9b5"][i], label=lab)
        ax.set_xticks(x, lanes["sample"], rotation=35, ha="right")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.18)
        ax.axhline(0, color="#17212f", linewidth=0.9)
    axes[0].set_ylabel("score")
    axes[0].legend(frameon=False, loc="upper left", ncol=1)
    fig.suptitle("Failure lanes / negative controls", y=0.98)
    for ax in axes:
        ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=10, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_NeoICI_failure_lanes.png", fig)

    # Figure D4: claim boundary diagram.
    fig, ax = plt.subplots(figsize=(10.5, 4.5))
    ax.axis("off")
    boxes = [
        (0.06, 0.58, 0.28, 0.26, "#e9f3e8", "Allowed now", "combo readiness\nopen replay phenotype\nreviewer-safe comparison"),
        (0.36, 0.58, 0.28, 0.26, "#fff1d4", "Allowed with caution", "NeoICI-style proxy > NeoPrecis-style proxy\nin this small replay only"),
        (0.66, 0.58, 0.28, 0.26, "#f7dedb", "Forbidden now", "clinical selection claim\nuniversal superiority claim"),
    ]
    for x0, y0, w, h, color, title, body in boxes:
        ax.add_patch(plt.Rectangle((x0, y0), w, h, facecolor=color, edgecolor="#17212f", linewidth=1.4))
        ax.text(x0 + 0.02, y0 + h - 0.05, title, fontsize=13, fontweight="bold", color="#102033")
        ax.text(x0 + 0.02, y0 + 0.08, body, fontsize=11, color="#102033", va="bottom", linespacing=1.35)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    ax.set_title("Claim boundary diagram", pad=18)
    save_fig(fig_dir / "Fig_NeoICI_claim_boundary.png", fig)

    # Figure D8: case-study radar summary.
    radar_labels = ["response", "clone", "cytotoxic", "activation", "memory"]
    case_vals = np.array([case_row["response_index"], case_row["top_clone_frac"], case_row["cytotoxic_score"], case_row["activation_score"], case_row["memory_score"]], dtype=float)
    cohort_radar = np.array([combo["response_index"].median(), combo["top_clone_frac"].median(), combo["cytotoxic_score"].median(), combo["activation_score"].median(), combo["memory_score"].median()], dtype=float)
    case_norm = (case_vals - case_vals.min()) / (case_vals.max() - case_vals.min() + 1e-9)
    cohort_norm = (cohort_radar - cohort_radar.min()) / (cohort_radar.max() - cohort_radar.min() + 1e-9)
    angles = np.linspace(0, 2*np.pi, len(radar_labels), endpoint=False).tolist()
    case_plot = np.r_[case_norm, case_norm[0]]
    cohort_plot = np.r_[cohort_norm, cohort_norm[0]]
    angles_plot = angles + [angles[0]]
    fig = plt.figure(figsize=(6.8, 6.4))
    ax = fig.add_subplot(111, polar=True)
    ax.plot(angles_plot, cohort_plot, color="#244e73", linewidth=2, label="cohort median")
    ax.fill(angles_plot, cohort_plot, color="#244e73", alpha=0.12)
    ax.plot(angles_plot, case_plot, color="#8f2d25", linewidth=2.5, label="CTMS1_followup4")
    ax.fill(angles_plot, case_plot, color="#8f2d25", alpha=0.18)
    ax.set_xticks(angles)
    ax.set_xticklabels(radar_labels)
    ax.set_yticklabels([])
    ax.set_title("CTMS1_followup4 radar summary", pad=20)
    ax.text(0.02, 0.98, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.25, 1.15))
    save_fig(fig_dir / "Fig_CTMS1_case_radar.png", fig)

    # Figure D9: nearest-neighbor style bar comparison for CTMS1_followup4.
    compare_metrics = ["response_index", "top_clone_frac", "cytotoxic_score", "activation_score", "memory_score"]
    neighbors = proxy_samples.copy()
    scaler = neighbors[compare_metrics].copy()
    scaler = (scaler - scaler.mean()) / (scaler.std(ddof=0) + 1e-9)
    case_vec = scaler.loc[neighbors["sample"] == "CTMS1_followup4"].iloc[0]
    dist = ((scaler - case_vec) ** 2).sum(axis=1).pow(0.5)
    nn = neighbors.assign(distance=dist).sort_values("distance").head(4)
    nn = nn.assign(label=nn["sample"].str.replace("CTMS1_", "", regex=False))
    fig, ax = plt.subplots(figsize=(10.4, 5.8))
    x = np.arange(len(compare_metrics))
    width = 0.18
    colors = ["#8f2d25", "#b58534", "#244e73", "#6a8c6d"]
    for i, row in enumerate(nn.itertuples(index=False)):
        vals = [getattr(row, m) for m in compare_metrics]
        ax.bar(x + (i - 1.5) * width, vals, width=width, color=colors[i], label=row.label)
    ax.set_xticks(x, ["response", "clone", "cytotoxic", "activation", "memory"])
    ax.set_ylabel("raw score")
    ax.set_title("CTMS1_followup4 nearest-neighbor style comparison")
    ax.legend(frameon=False, ncol=2, loc="upper left")
    ax.grid(axis="y", alpha=0.18)
    ax.text(0.02, 0.95, "★ OURS", transform=ax.transAxes, fontsize=11, fontweight="bold", color="#8f2d25", va="top")
    save_fig(fig_dir / "Fig_CTMS1_case_neighbor_bars.png", fig)

    # Figure E: five attention-style heatmaps derived from the replay signals.
    ordered = proxy_samples.sort_values("response_index", ascending=False).copy()
    ordered["label"] = ordered["sample"].str.replace("CTMS1_", "", regex=False).str.replace("Pt_#", "Pt", regex=False)
    heat_cols = [
        ("response_index", "response"),
        ("top_clone_frac", "clone"),
        ("cytotoxic_score", "cytotoxic"),
        ("activation_score", "activation"),
        ("memory_score", "memory"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.6))
    axes = axes.flatten()
    for i, (col, title) in enumerate(heat_cols):
        vec = ordered[col].to_numpy()
        outer = np.outer(vec, vec)
        m = np.abs(outer).max() or 1.0
        im = axes[i].imshow(outer, cmap="RdBu_r", vmin=-m, vmax=m, aspect="auto")
        axes[i].set_title(f"{title} attention-style map")
        axes[i].set_xticks(range(len(ordered)))
        axes[i].set_xticklabels(ordered["label"], rotation=90, fontsize=7)
        axes[i].set_yticks(range(len(ordered)))
        axes[i].set_yticklabels(ordered["label"], fontsize=7)
        axes[i].text(0.02, 0.95, "★ OURS", transform=axes[i].transAxes, fontsize=9, fontweight="bold", color="#8f2d25", va="top")
    axes[5].axis("off")
    cbar = fig.colorbar(im, ax=axes[:5], shrink=0.72, pad=0.02)
    cbar.set_label("pairwise score coupling", rotation=90)
    fig.suptitle("CTMS1 / replay attention-style heatmaps", y=0.98, fontsize=15)
    save_fig(fig_dir / "Fig_NeoICI_attention_heatmaps.png", fig)

    intro_df = pd.DataFrame(
        [
            {
                "question": "이 페이지가 뭔가?",
                "answer": "NeoICI 브랜치를 위한 paper1-style 웹 dossier입니다. 알고리즘 benchmarking, 공개 combo replay, proxy comparison, 그리고 대표 case study를 한 흐름으로 묶었습니다.",
                "why_it_matters": "리뷰어가 claim → figure → raw table 순서로 바로 따라올 수 있게 해줍니다.",
            },
            {
                "question": "NeoICI가 푸는 문제는?",
                "answer": "단순 peptide rank가 아니라 vaccine + ICI 준비도(readiness)를 우선합니다. 즉, binding-like intuition보다 immune-state support와 T-cell activation을 더 중요하게 봅니다.",
                "why_it_matters": "병용치료는 peptide-only immunogenicity와 다른 신호를 요구합니다.",
            },
            {
                "question": "핵심 비교는?",
                "answer": "같은 9개 paired sample에서 NeoICI-style proxy와 NeoPrecis-style proxy를 동일한 response_index target에 대해 직접 비교합니다.",
                "why_it_matters": "이게 framing이 실제 biology와 맞는지 보는 가장 깔끔한 방법입니다.",
            },
            {
                "question": "가장 강한 샘플은?",
                "answer": "CTMS1_followup4입니다. response_index 0.4187, top_clone_frac 0.1123, 그리고 NeoICI-style rank 1/9를 가집니다.",
                "why_it_matters": "open combo replay에서 가장 설명력이 높은 sample-level panel입니다.",
            },
        ]
    )

    workflow_df = pd.DataFrame(
        [
            {
                "step": "1. 로컬 알고리즘부터 때린다",
                "what_we_did": "NeoICI gauntlet을 돌려 KG_GA_evolved, BioDarwin, CROSS, BigMHC, claim-safe comparator와 직접 붙였습니다.",
                "evidence": "KG_GA_evolved 0.978; NeoICI 0.863; industrial locked 0.563.",
                "why_it_matters": "replay를 보기 전에 먼저 과장 claim을 막아야 합니다.",
            },
            {
                "step": "2. 공개 vaccine + ICI cohort 재생",
                "what_we_did": "GSE222011과 GSE255830을 paired GEX/TCR sample-level summary로 정리했습니다.",
                "evidence": "9개 paired sample retained; response_index median 0.033.",
                "why_it_matters": "병용치료 phenotype이 실제 데이터에 존재하는지 보여줍니다.",
            },
            {
                "step": "3. proxy family 직접 비교",
                "what_we_did": "같은 sample에 NeoICI-style proxy와 NeoPrecis-style proxy를 동시에 점수화했습니다.",
                "evidence": "Spearman 0.817 vs 0.417; leave-one-out median 0.833 vs 0.381.",
                "why_it_matters": "NeoICI framing이 실제로 더 잘 맞는지 보여주는 가장 강한 근거입니다.",
            },
            {
                "step": "4. 대표 case를 깊게 판다",
                "what_we_did": "CTMS1_followup4를 explanation lane으로 열고 clonotype/activation pattern을 요약했습니다.",
                "evidence": "response_index 0.4187; NeoICI-style proxy에서 rank 1/9.",
                "why_it_matters": "리뷰어는 평균만이 아니라 구체적인 case 하나를 원합니다.",
            },
        ]
    )

    importance_df = pd.DataFrame(
        [
            {"tag": "★ 왜 중요한가", "text": "이 페이지는 algorithm ranking과 treatment-context biology를 이어주는 다리입니다."},
            {"tag": "★ 리뷰어 시선", "text": "contract를 분명히 잡습니다: strong benchmark, small replay, bounded claim."},
            {"tag": "★ 서사의 축", "text": "가장 강한 증거는 하나의 scalar가 아니라 proxy, clonotype, response-like state의 정렬입니다."},
            {"tag": "★ 최종 스탠스", "text": "NeoICI는 Nature Cancer 맥락에서 말할 수 있는 readiness layer이지, universal superiority를 가장하는 도구가 아닙니다."},
        ]
    )

    decision_df = pd.DataFrame(
        [
            {
                "Question": "NeoICI가 universal scorer인가?",
                "Disposition": "Forbidden",
                "Evidence": "KG_GA_evolved remains local top scorer at 0.978; NeoICI is competitive but not #1.",
                "What would upgrade it": "New external benchmarks where NeoICI is consistently top-ranked.",
            },
            {
                "Question": "NeoICI가 combo-readiness layer인가?",
                "Disposition": "Allowed now",
                "Evidence": "Open combo replay and CTMS1_followup4 support a readiness-style framing.",
                "What would upgrade it": "Larger open cohorts with held-out replication.",
            },
            {
                "Question": "NeoICI-style proxy가 response-like index를 더 잘 따라가나?",
                "Disposition": "Allowed with caution",
                "Evidence": "Spearman 0.817 vs 0.417; LOO median 0.833 vs 0.381.",
                "What would upgrade it": "n larger than 9 and an external replay set.",
            },
            {
                "Question": "CTMS1_followup4가 설명 가능한 anchor case인가?",
                "Disposition": "Allowed now",
                "Evidence": "response_index 0.4187; NeoICI rank 1/9; top_clone_frac 0.1123.",
                "What would upgrade it": "Independent case series with same pattern.",
            },
            {
                "Question": "clinical selection claim을 해도 되나?",
                "Disposition": "Forbidden",
                "Evidence": "Everything here is post-hoc replay, not prospective selection.",
                "What would upgrade it": "Prospective clinical validation and outcome linkage.",
            },
            {
                "Question": "negative controls가 충분한가?",
                "Disposition": "Allowed now",
                "Evidence": "Fig. 4H shows low-response lanes beside top-response lanes.",
                "What would upgrade it": "More low-lane samples from a second open cohort.",
            },
        ]
    )

    reproducibility_df = pd.DataFrame(
        [
            {
                "Check": "input cohorts",
                "What to verify": "GSE222011 and GSE255830 paired GEX/TCR samples only",
                "Where": "project/results/p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11/combo_geo_reanalysis.tsv",
            },
            {
                "Check": "retained sample count",
                "What to verify": "exactly 9 paired samples retained after filtering",
                "Where": "same replay TSV + Fig. 3",
            },
            {
                "Check": "local benchmark ranking",
                "What to verify": "KG_GA_evolved stays the local top scorer at AUPRC 0.978",
                "Where": "NEOICI_NATURE_CANCER_EXECUTIVE_TABLE.tsv + Fig. 2 / 2A",
            },
            {
                "Check": "proxy comparison",
                "What to verify": "NeoICI-style proxy rho 0.817 vs NeoPrecis-style rho 0.417",
                "Where": "proxy_compare_metrics.tsv + Fig. 4",
            },
            {
                "Check": "robustness",
                "What to verify": "LOO median remains higher for NeoICI-style proxy",
                "Where": "NEOICI_NATURE_CANCER_LOO_TABLE.tsv + Fig. 4C",
            },
            {
                "Check": "negative control",
                "What to verify": "Fig. 4H includes low-response lanes, not just the anchor case",
                "Where": "Fig. 4H + failure_small table",
            },
            {
                "Check": "claim boundary",
                "What to verify": "readiness layer yes; universal superiority no; clinical selection no",
                "Where": "NEOICI_NATURE_CANCER_CLAIM_TABLE.tsv + Fig. 5D",
            },
        ]
    )

    methods_df = pd.DataFrame(
        [
            {"Item": "response_index", "Definition": "module-based open combo outcome proxy", "Interpretation": "higher = more response-like"},
            {"Item": "NeoICI_style_proxy", "Definition": "T-cell / state-heavy proxy", "Interpretation": "higher = more NeoICI-like"},
            {"Item": "NeoPrecis_style_proxy", "Definition": "published-contract style proxy", "Interpretation": "higher = more NeoPrecis-like"},
            {"Item": "top_clone_frac", "Definition": "largest clonotype fraction", "Interpretation": "higher = more clonal expansion"},
            {"Item": "claim boundary", "Definition": "allowed / caution / forbidden claim ladder", "Interpretation": "prevents overclaiming"},
        ]
    )

    provenance_df = pd.DataFrame(
        [
            {"Artifact": "page generator", "Path": "scripts/build_neoici_nature_cancer_package_2026_05_11.py", "Role": "single source of truth"},
            {"Artifact": "live HTML", "Path": "project/papers_hub_2026_05_04/neoici_nature_cancer_package.html", "Role": "paper1-style dossier"},
            {"Artifact": "alias HTML", "Path": "project/papers_hub_2026_05_04/neoici_nature_cancer_paper1_style.html", "Role": "same content, alternate URL"},
            {"Artifact": "live mirror", "Path": "/var/www/papers/papers_hub_2026_05_04/", "Role": "served to browser"},
            {"Artifact": "figure assets", "Path": "project/papers_hub_2026_05_04/assets/neoici_nature_cancer_package/", "Role": "web-served figures"},
            {"Artifact": "local results", "Path": "project/results/p_neo_ici_nature_cancer_package_2026_05_11/", "Role": "generated output bundle"},
        ]
    )

    package_summary = [
        "# NeoICI Nature Cancer package",
        "",
        "## 위치",
        "",
        "- NeoICI는 combo-therapy readiness layer입니다.",
        "- KG_GA_evolved는 validation-like split에서 여전히 local universal scorer 1위입니다.",
        "- NeoPrecis는 published benchmark contract로서 immunogenicity / clonality-aware framing이 더 강합니다.",
        "",
        "## 무엇이 새로 들어갔나",
        "",
        "- 공개 combo dataset registry와 execution ladder.",
        "- NeoICI vs BigMHC / CROSS / KG-GA / BioDarwin 비교가 들어간 local algorithm gauntlet.",
        "- 9개 retained sample의 open combo GEO replay.",
        "- 같은 sample에서 NeoICI-style proxy가 response-like index를 더 잘 따라간다는 비교.",
        "- CTMS1_followup4 case-study panel.",
        "",
        "## claim 규칙",
        "",
        "- Allowed: combo readiness, response-like replay, phenotype support, reviewer-safe benchmark comparison.",
        "- Forbidden: clinical selection, universal superiority, thyroid ICI predictor claim.",
    ]
    write_text(OUT / "NEOICI_NATURE_CANCER_PACKAGE_SUMMARY.md", "\n".join(package_summary) + "\n")

    html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>NeoICI Nature Cancer package | THCA Hub</title>
<link rel="preconnect" href="https://fonts.googleapis.com" />
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin />
<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:ital,wght@0,400;0,500;0,600;1,400&family=JetBrains+Mono:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700;800&display=swap" rel="stylesheet" />
<style>
:root{{
  --ink:#102033; --muted:#52606f; --paper:#f6f0e6; --card:#fffaf1;
  --line:#d4c6b3; --line2:#eadfcc; --red:#8f2d25; --red2:#b64b3f;
  --blue:#244e73; --green:#426b50; --gold:#b58534; --coal:#17212f;
  --cream:#fff7e4; --good:#e9f3e8; --warn:#fff1d4; --bad:#f7dedb;
}}
*{{box-sizing:border-box}}
html{{scroll-behavior:smooth}}
body{{
  margin:0;
  background:
    radial-gradient(circle at 8% 0%, rgba(181,133,52,.18), transparent 28rem),
    linear-gradient(135deg,#f8f2e8 0%,#efe3d1 55%,#f7efe0 100%);
  color:var(--ink);font-family:"Newsreader","Noto Sans KR",serif;line-height:1.68;font-size:17px
}}
a{{color:var(--red);text-decoration:none;border-bottom:1px solid rgba(143,45,37,.35)}}
.wrap{{max-width:1200px;margin:0 auto;padding:34px 28px 80px}}
.topnav{{display:flex;justify-content:space-between;gap:18px;align-items:center;font-family:"JetBrains Mono",monospace;font-size:12px;margin-bottom:28px;color:var(--muted);flex-wrap:wrap}}
.topnav a{{border:1px solid var(--line);padding:8px 12px;background:rgba(255,255,255,.55);border-radius:999px}}
.hero{{border:2px solid var(--coal);background:linear-gradient(145deg,#fffaf1 0%,#f0dfc4 100%);padding:44px 46px;box-shadow:10px 10px 0 rgba(23,33,47,.16);position:relative;overflow:hidden}}
.hero:after{{content:"";position:absolute;right:-80px;top:-80px;width:260px;height:260px;border:38px solid rgba(143,45,37,.12);border-radius:50%}}
.kicker,.mono{{font-family:"JetBrains Mono",monospace}}
.kicker{{letter-spacing:.18em;text-transform:uppercase;font-size:12px;color:var(--red);font-weight:700;margin-bottom:18px}}
h1{{font-family:"Cormorant Garamond",serif;font-size:58px;line-height:1.02;margin:0 0 14px;color:var(--coal);max-width:920px}}
.subtitle{{font-size:22px;color:var(--muted);max-width:900px;font-style:italic}}
.meta{{display:flex;flex-wrap:wrap;gap:10px;margin-top:26px}}
.chip{{font-family:"JetBrains Mono",monospace;font-size:11px;border:1px solid var(--line);background:#fff;padding:7px 10px;border-radius:6px}}
.chip.hot{{background:var(--bad);border-color:#d9aaa3;color:var(--red);font-weight:700}}
.chip.good{{background:var(--good);border-color:#b8d0b5;color:var(--green);font-weight:700}}
.chip.ours{{background:#fff1ef;border-color:var(--red);color:var(--red);font-weight:800}}
.layout{{display:grid;grid-template-columns:260px 1fr;gap:28px;margin-top:34px}}
.toc{{position:sticky;top:18px;align-self:start;background:rgba(255,250,241,.9);border:1px solid var(--line);padding:18px;border-radius:12px}}
.toc h3{{font-family:"Cormorant Garamond",serif;margin:0 0 12px;font-size:23px;color:var(--red)}}
.toc a{{display:block;border-bottom:1px dotted var(--line);padding:8px 0;color:var(--ink);font-size:14px}}
section{{background:rgba(255,250,241,.86);border:1px solid var(--line);padding:32px;margin-bottom:24px;border-radius:14px}}
h2{{font-family:"Cormorant Garamond",serif;font-size:38px;line-height:1.1;margin:0 0 12px;color:var(--coal)}}
h3{{font-family:"Cormorant Garamond",serif;font-size:27px;margin:26px 0 10px;color:var(--red)}}
.section-note{{color:var(--muted);font-style:italic;margin-bottom:20px}}
.grid{{display:grid;gap:16px}}
.g2{{grid-template-columns:repeat(2,minmax(0,1fr))}}
.g3{{grid-template-columns:repeat(3,minmax(0,1fr))}}
.g4{{grid-template-columns:repeat(4,minmax(0,1fr))}}
.card{{background:#fff;border:1px solid var(--line2);padding:18px;border-radius:12px}}
.card h4{{margin:0 0 8px;font-family:"JetBrains Mono",monospace;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--red)}}
.mini-label{{font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.1em;color:var(--muted);font-weight:700}}
.callout{{border-left:5px solid var(--red);background:#fff;padding:20px 22px;margin:20px 0;border-radius:0 10px 10px 0}}
.callout.good{{border-color:var(--green);background:var(--good)}}
.callout.warn{{border-color:var(--gold);background:var(--warn)}}
.callout.bad{{border-color:var(--red);background:var(--bad)}}
.starbox{{border:2px solid var(--red);background:#fff1ef;box-shadow:5px 5px 0 rgba(143,45,37,.16)}}
.starbox strong{{color:var(--red)}}
.ours-flag{{display:inline-block;margin-left:8px;padding:2px 7px;border:1px solid var(--red);border-radius:999px;background:#fff1ef;color:var(--red);font-family:"JetBrains Mono",monospace;font-size:10px;letter-spacing:.08em;vertical-align:middle}}
figure{{margin:22px 0;background:#fff;border:1px solid var(--line);border-radius:14px;padding:14px;box-shadow:0 6px 20px rgba(23,33,47,.06)}}
figure img{{display:block;width:100%;height:auto;border-radius:8px;background:#f4eadb;cursor:zoom-in}}
figure.driver-wide{{width:min(1600px,calc(100vw - 42px));margin-left:50%;transform:translateX(-50%);border:2px solid var(--gold);box-shadow:8px 8px 0 rgba(181,133,52,.18);overflow-x:auto}}
figure.driver-wide img{{background:#fff;max-width:none;width:100%;min-width:1380px}}
figure.driver-wide figcaption{{max-width:1180px;margin:0 auto}}
figcaption{{font-size:14px;color:var(--muted);padding:12px 6px 4px;line-height:1.55}}
figcaption strong{{color:var(--ink)}}
.fig-notes{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin-top:10px}}
.fig-notes div{{background:#fff8e8;border:1px solid var(--line2);border-radius:8px;padding:9px;color:var(--ink)}}
.fig-notes b{{font-family:"JetBrains Mono",monospace;font-size:10px;text-transform:uppercase;letter-spacing:.08em;color:var(--red);display:block;margin-bottom:3px}}
table{{width:100%;border-collapse:collapse;margin:18px 0;background:#fff;font-size:14px;border:1px solid var(--line)}}
caption{{text-align:left;font-family:"Cormorant Garamond",serif;font-size:20px;color:var(--red);font-weight:700;margin:8px 0}}
th{{background:var(--coal);color:#fff;text-align:left;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.06em;text-transform:uppercase;padding:10px}}
td{{border-top:1px solid var(--line2);padding:10px;vertical-align:top}}
tr:nth-child(even) td{{background:#fffdf7}}
tr.hit td{{background:var(--good)}}
tr.warn td{{background:var(--warn)}}
tr.drop td{{background:var(--bad)}}
.status{{font-family:"JetBrains Mono",monospace;font-size:11px;padding:3px 8px;border-radius:999px;font-weight:700;display:inline-block}}
.status.hit{{background:var(--green);color:#fff}}
.status.warn{{background:var(--gold);color:#fff}}
.status.drop{{background:var(--red);color:#fff}}
.table-wrap{{overflow-x:auto}}
ul,ol{{padding-left:22px}}
li{{margin:6px 0}}
.footer{{font-family:"JetBrains Mono",monospace;font-size:12px;color:var(--muted);text-align:center;margin-top:34px;padding-top:18px;border-top:1px solid var(--line)}}
.source-list code{{background:#fff;padding:1px 4px;border-radius:4px}}
.statrow{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px;margin-top:18px}}
.stat{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:14px}}
.stat .k{{font-family:"JetBrains Mono",monospace;font-size:10px;color:var(--red);text-transform:uppercase;letter-spacing:.08em;font-weight:700}}
.stat .v{{font-family:"Cormorant Garamond",serif;font-size:36px;line-height:1.05;margin-top:5px;color:var(--coal);font-weight:700}}
.stat .t{{font-size:13px;color:var(--muted);margin-top:6px;line-height:1.45}}
.impact{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px;margin:18px 0 0}}
.impact .panel{{background:#fff;border:2px solid var(--coal);border-radius:14px;padding:16px 16px 18px;box-shadow:6px 6px 0 rgba(23,33,47,.12);min-height:160px;position:relative}}
.impact .panel .tag{{font-family:"JetBrains Mono",monospace;font-size:10px;letter-spacing:.12em;text-transform:uppercase;color:var(--red);font-weight:800}}
.impact .panel .big{{font-family:"Cormorant Garamond",serif;font-size:31px;line-height:1.02;color:var(--coal);font-weight:700;margin:8px 0 8px}}
.impact .panel .body{{font-size:13px;color:var(--muted);line-height:1.5}}
.impact .panel.hot{{background:linear-gradient(180deg,#fff1ef 0%,#fff 100%);border-color:var(--red)}}
.impact .panel.good{{background:linear-gradient(180deg,#eef5ea 0%,#fff 100%);border-color:var(--green)}}
.impact .panel.warn{{background:linear-gradient(180deg,#fff5df 0%,#fff 100%);border-color:var(--gold)}}
.impact .panel.blue{{background:linear-gradient(180deg,#edf5fb 0%,#fff 100%);border-color:var(--blue)}}
.impact .panel .flag{{position:absolute;top:12px;right:12px;font-family:"JetBrains Mono",monospace;font-size:10px;font-weight:800;color:var(--red);border:1px solid var(--red);background:#fff1ef;padding:3px 7px;border-radius:999px}}
.verdict{{margin:18px 0 2px;border:2px solid var(--red);background:linear-gradient(180deg,#fff1ef 0%,#fffaf1 100%);box-shadow:7px 7px 0 rgba(143,45,37,.14);padding:18px 18px 16px;border-radius:14px}}
.verdict .headline{{font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.16em;text-transform:uppercase;color:var(--red);font-weight:800}}
.verdict .title{{font-family:"Cormorant Garamond",serif;font-size:36px;line-height:1.02;color:var(--coal);font-weight:700;margin:6px 0 10px}}
.verdict .strap{{font-size:15px;color:var(--muted);max-width:1040px;line-height:1.55}}
.verdict-grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin-top:14px}}
.verdict-grid .box{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:14px 14px 16px;min-height:118px}}
.verdict-grid .box h5{{margin:0 0 8px;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.1em;text-transform:uppercase;color:var(--red)}}
.verdict-grid .box .num{{font-family:"Cormorant Garamond",serif;font-size:30px;line-height:1;color:var(--coal);font-weight:700}}
.verdict-grid .box .txt{{font-size:13px;color:var(--muted);line-height:1.45;margin-top:6px}}
.verdict .bottom{{margin-top:12px;padding:10px 12px;border-radius:10px;background:#fff;border:1px dashed var(--red);font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.05em;color:var(--red);font-weight:800}}
.hard-banner{{margin:0 0 22px;border:2px solid var(--coal);background:linear-gradient(90deg,#17212f 0%,#7b1f2a 48%,#b58534 100%);color:#fff;padding:16px 18px;border-radius:14px;box-shadow:8px 8px 0 rgba(23,33,47,.16)}}
.hard-banner .line{{font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.18em;text-transform:uppercase;font-weight:800;opacity:.95}}
.hard-banner .main{{font-family:"Cormorant Garamond",serif;font-size:34px;line-height:1.02;font-weight:700;margin:6px 0 8px;max-width:1100px}}
.hard-banner .sub{{display:flex;flex-wrap:wrap;gap:10px;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.05em;font-weight:700}}
.hard-banner .sub span{{background:rgba(255,255,255,.14);border:1px solid rgba(255,255,255,.26);padding:6px 9px;border-radius:999px}}
.mega-strip{{margin:0 0 22px;border:2px solid var(--coal);background:#fff;box-shadow:10px 10px 0 rgba(23,33,47,.14);border-radius:16px;overflow:hidden;display:grid;grid-template-columns:1.15fr .85fr}}
.mega-strip .left{{padding:22px 22px 20px;background:linear-gradient(135deg,#fff8f2 0%,#fff1ef 50%,#fffdf6 100%)}}
.mega-strip .left .k{{font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.16em;text-transform:uppercase;color:var(--red);font-weight:800}}
.mega-strip .left .title{{font-family:"Cormorant Garamond",serif;font-size:44px;line-height:0.98;color:var(--coal);font-weight:700;margin:8px 0 10px;max-width:760px}}
.mega-strip .left .sub{{font-size:16px;line-height:1.55;color:var(--muted);max-width:820px}}
.mega-strip .left .stamp{{display:inline-block;margin-top:14px;padding:8px 10px;border:1px solid var(--red);background:#fff1ef;color:var(--red);font-family:"JetBrains Mono",monospace;font-size:11px;font-weight:800;letter-spacing:.08em;border-radius:999px}}
.mega-strip .right{{display:grid;grid-template-columns:1fr 1fr;gap:0;border-left:1px solid var(--line2);background:#fff}}
.mega-strip .right .cell{{padding:18px 16px;border-right:1px solid var(--line2);border-bottom:1px solid var(--line2);min-height:132px}}
.mega-strip .right .cell:nth-child(2n){{border-right:0}}
.mega-strip .right .cell:nth-last-child(-n+2){{border-bottom:0}}
.mega-strip .right .cell .tag{{font-family:"JetBrains Mono",monospace;font-size:10px;letter-spacing:.14em;text-transform:uppercase;color:var(--red);font-weight:800}}
.mega-strip .right .cell .num{{font-family:"Cormorant Garamond",serif;font-size:50px;line-height:1;color:var(--coal);font-weight:700;margin:8px 0 4px}}
.mega-strip .right .cell .txt{{font-size:13px;color:var(--muted);line-height:1.45}}
.mega-strip .right .cell.red{{background:#fff1ef}}
.mega-strip .right .cell.blue{{background:#edf5fb}}
.mega-strip .right .cell.gold{{background:#fff5df}}
.mega-strip .right .cell.green{{background:#eef5ea}}
.headline-ribbon{{margin:0 0 14px;background:linear-gradient(90deg,#17212f 0%,#8f2d25 55%,#b58534 100%);color:#fff;border-bottom:2px solid #17212f;box-shadow:0 6px 16px rgba(23,33,47,.18)}}
.headline-ribbon .inner{{max-width:1200px;margin:0 auto;padding:10px 28px;display:flex;flex-wrap:wrap;gap:12px;align-items:center;justify-content:space-between}}
.headline-ribbon .main{{font-family:"JetBrains Mono",monospace;font-size:12px;font-weight:800;letter-spacing:.14em;text-transform:uppercase}}
.headline-ribbon .sub{{display:flex;flex-wrap:wrap;gap:8px;font-family:"JetBrains Mono",monospace;font-size:11px;font-weight:700}}
.headline-ribbon .sub span{{padding:5px 8px;border:1px solid rgba(255,255,255,.3);border-radius:999px;background:rgba(255,255,255,.12)}}
.topline{{max-width:1200px;margin:0 auto 8px;padding:0 28px}}
.topline .box{{background:#17212f;color:#fff;border:2px solid #8f2d25;border-radius:12px;padding:10px 14px;display:flex;flex-wrap:wrap;gap:10px;align-items:center;justify-content:space-between;box-shadow:4px 4px 0 rgba(23,33,47,.14)}}
.topline .main{{font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.18em;text-transform:uppercase;font-weight:800}}
.topline .chips{{display:flex;flex-wrap:wrap;gap:8px;font-family:"JetBrains Mono",monospace;font-size:10px;font-weight:700}}
.topline .chips span{{padding:4px 7px;border:1px solid rgba(255,255,255,.26);border-radius:999px;background:rgba(255,255,255,.12)}}
.cover-masthead{{max-width:1200px;margin:0 auto 10px;padding:0 28px;font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.16em;text-transform:uppercase;color:var(--red);font-weight:800}}
.cover-masthead .line{{display:flex;flex-wrap:wrap;gap:10px;align-items:center}}
.cover-masthead .line span{{padding:5px 8px;border:1px solid var(--line2);border-radius:999px;background:#fff6eb;color:var(--coal)}}
.opening-verdict{{max-width:1200px;margin:0 auto 8px;padding:0 28px}}
.opening-verdict .box{{background:#8f2d25;color:#fff;border:2px solid #17212f;border-radius:14px;padding:12px 16px;display:flex;flex-wrap:wrap;gap:10px;align-items:center;justify-content:space-between;box-shadow:6px 6px 0 rgba(23,33,47,.14)}}
.opening-verdict .main{{font-family:"Cormorant Garamond",serif;font-size:26px;line-height:1.02;font-weight:700}}
.opening-verdict .chipline{{display:flex;flex-wrap:wrap;gap:8px;font-family:"JetBrains Mono",monospace;font-size:10px;font-weight:800;letter-spacing:.08em}}
.opening-verdict .chipline span{{padding:4px 7px;border:1px solid rgba(255,255,255,.26);border-radius:999px;background:rgba(255,255,255,.12)}}
.metric-wall{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin:18px 0 16px}}
.metric-wall .tile{{background:#fff;border:1px solid var(--line2);border-radius:14px;padding:16px 16px 18px;box-shadow:4px 4px 0 rgba(23,33,47,.08);position:relative;min-height:138px}}
.metric-wall .tile::before{{content:"";position:absolute;inset:0 0 auto 0;height:6px;background:linear-gradient(90deg,var(--red),var(--gold),var(--blue))}}
.metric-wall .tile .label{{font-family:"JetBrains Mono",monospace;font-size:10px;letter-spacing:.14em;text-transform:uppercase;color:var(--red);font-weight:800;margin-top:6px}}
.metric-wall .tile .num{{font-family:"Cormorant Garamond",serif;font-size:56px;line-height:1;color:var(--coal);font-weight:700;margin:8px 0 4px}}
.metric-wall .tile .desc{{font-size:13px;color:var(--muted);line-height:1.45}}
.metric-wall .tile.hot{{border-color:var(--red);background:linear-gradient(180deg,#fff1ef 0%,#fff 100%)}}
.metric-wall .tile.blue{{border-color:var(--blue);background:linear-gradient(180deg,#edf5fb 0%,#fff 100%)}}
.metric-wall .tile.good{{border-color:var(--green);background:linear-gradient(180deg,#eef5ea 0%,#fff 100%)}}
.workflow{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px;margin-top:10px}}
.workflow .step{{background:#fff;border:1px solid var(--line2);border-radius:12px;padding:14px 14px 16px;position:relative;min-height:160px}}
.workflow .step .n{{font-family:"JetBrains Mono",monospace;font-size:11px;letter-spacing:.08em;color:var(--red);font-weight:700;text-transform:uppercase}}
.workflow .step .t{{font-family:"Cormorant Garamond",serif;font-size:22px;line-height:1.1;color:var(--coal);font-weight:700;margin:8px 0 6px}}
.workflow .step .d{{font-size:13px;line-height:1.5;color:var(--muted)}}
@media(max-width:900px){{
  .layout{{grid-template-columns:1fr}}
  .toc{{position:relative;top:auto}}
  .g2,.g3,.statrow,.fig-notes{{grid-template-columns:1fr}}
  .impact{{grid-template-columns:1fr}}
  .verdict-grid{{grid-template-columns:1fr}}
  .metric-wall{{grid-template-columns:1fr}}
  .mega-strip{{grid-template-columns:1fr}}
  .mega-strip .right{{grid-template-columns:1fr 1fr}}
  .headline-ribbon .inner{{padding:10px 14px}}
  .topline{{padding:0 14px}}
  .cover-masthead{{padding:0 14px}}
  .opening-verdict{{padding:0 14px}}
  h1{{font-size:40px}}
  .hero{{padding:30px 24px}}
  .wrap{{padding:20px 14px}}
  figure.driver-wide{{width:100%;margin-left:0;transform:none;overflow-x:visible}}
  figure.driver-wide img{{min-width:0;width:100%;max-width:100%}}
}}
@media print{{body{{background:#fff}}.toc,.topnav{{display:none}}.layout{{display:block}}.hero,section,figure,.card{{break-inside:avoid;box-shadow:none}}}}
</style>
</head>
<body>
<figure class="driver-wide" style="margin:0 0 16px">
  <img src="assets/{ASSET_NAME}/Fig_NeoICI_impact_banner.png" alt="NeoICI impact banner" />
  <figcaption><strong>Fig. 0O</strong> <span class="ours-flag">OURS</span> impact banner. 첫 화면에서 결론과 숫자를 단일 배너로 밀어 넣은 상단 결론판이다.</figcaption>
</figure>
<figure class="driver-wide" style="margin:0 0 20px">
  <img src="assets/{ASSET_NAME}/Fig_NeoICI_cover_wall.png" alt="NeoICI cover wall" />
  <figcaption><strong>Fig. 0G</strong> <span class="ours-flag">OURS</span> first-screen cover wall. 페이지를 여는 순간 핵심 결론이 먼저 보이도록 만든 최상단 포스터다.</figcaption>
</figure>
<section id="data-definitions" style="margin:0 0 24px;padding:18px 22px;border:1px solid #c9c9c9;border-left:4px solid #7B1F2A;background:#fbf7f4;font-family:'Newsreader','Noto Sans KR',serif;">
  <details open style="margin:0;">
    <summary style="cursor:pointer;font-family:'JetBrains Mono',monospace;font-size:12px;letter-spacing:.06em;color:#7B1F2A;font-weight:700;">데이터 정의서 / DATA DEFINITIONS v1 — NeoICI package page (open by default)</summary>
    <div style="margin-top:14px;font-size:13px;line-height:1.55;">
      <p style="margin:0 0 10px;color:#444">이 페이지가 실제로 참조하는 cohort / score / figure / claim 정의만 넣었습니다. voice-protected manuscript sections는 건드리지 않습니다.</p>
      <h4 style="margin:14px 0 6px;font-family:'JetBrains Mono',monospace;font-size:11px;letter-spacing:.08em;color:#7B1F2A">A. COHORTS / REPLAY SETS</h4>
      {block_table(pd.DataFrame([
          ["GSE222011", "autogene cevumeran + atezolizumab + mFOLFIRINOX", "paired GEX/TCR", 5, "open combo replay"],
          ["GSE255830", "GNOS-PV02 + INO-9012 + pembrolizumab", "paired GEX/TCR", 4, "open combo replay"],
          ["Total retained", "two public combo cohorts", "paired samples", 9, "response-like replay"],
      ], columns=["ID", "Source", "Modality", "n", "Use"]))}
      <h4 style="margin:14px 0 6px;font-family:'JetBrains Mono',monospace;font-size:11px;letter-spacing:.08em;color:#7B1F2A">B. SCORES / PROXIES</h4>
      {block_table(pd.DataFrame([
          ["response_index", "module-based open combo outcome proxy", "higher = more response-like"],
          ["NeoICI_style_proxy", "T-cell/state-heavy proxy", "higher = more NeoICI-like"],
          ["NeoPrecis_style_proxy", "published-contract style proxy", "higher = more NeoPrecis-like"],
          ["top_clone_frac", "largest clonotype fraction", "higher = more clonal expansion"],
      ], columns=["Score", "Definition", "Higher ="]))}
      <h4 style="margin:14px 0 6px;font-family:'JetBrains Mono',monospace;font-size:11px;letter-spacing:.08em;color:#7B1F2A">C. FIGURES</h4>
      {block_table(figure_df)}
      <p style="margin:12px 0 0;font-size:11px;color:#777;font-family:'JetBrains Mono',monospace">source tables: <code>neoici_algorithm_gauntlet_metrics.tsv</code>, <code>combo_geo_reanalysis.tsv</code>, <code>proxy_compare_metrics.tsv</code>, <code>NEOICI_NATURE_CANCER_*.tsv</code></p>
    </div>
  </details>
</section>

<div class="wrap">
  <nav class="topnav">
    <a href="index.html">8-Papers Hub</a>
    <span>NeoICI Nature Cancer package · paper1-style dossier · OURS</span>
    <a href="neoici_nature_cancer_synthesis.html">Synthesis</a>
    <a href="ctms1_followup4_case_study.html">Case study</a>
  </nav>

  <div class="hard-banner">
    <div class="line">Immediate verdict · OURS · paper1-style web dossier</div>
    <div class="main">NeoICI는 병용치료를 위한 readability layer다. local top scorer를 숨기지 않고, 공개 replay에서 방향을 보여주고, 대표 case로 그 방향을 설명한다.</div>
    <div class="sub">
      <span>KG_GA_evolved AUPRC 0.978</span>
      <span>NeoICI-style rho 0.817</span>
      <span>CTMS1_followup4 rank 1/9</span>
      <span>Claim boundary locked</span>
    </div>
  </div>

    <div class="mega-strip">
    <div class="left">
      <div class="k">NeoICI dossier / OURS</div>
      <div class="title">Figure-heavy, reviewer-facing, and brutally explicit about what the data does and does not support.</div>
      <div class="sub">이 페이지는 algorithm gauntlet, public combo replay, proxy comparison, and CTMS1_followup4 case study를 하나의 claim stack으로 묶는다. 상단에서 결론을 먼저 보여주고, 아래로 내려갈수록 근거가 촘촘해진다.</div>
      <div class="stamp">CLAIM BOUNDARY LOCKED</div>
    </div>
    <div class="right">
      <div class="cell red">
        <div class="tag">local scorer</div>
        <div class="num">0.978</div>
        <div class="txt">KG_GA_evolved remains the benchmark leader.</div>
      </div>
      <div class="cell blue">
        <div class="tag">replay rho</div>
        <div class="num">0.817</div>
        <div class="txt">NeoICI-style proxy tracks response better.</div>
      </div>
      <div class="cell gold">
        <div class="tag">rank anchor</div>
        <div class="num">1 / 9</div>
        <div class="txt">CTMS1_followup4 is the cleanest explanation lane.</div>
      </div>
      <div class="cell green">
        <div class="tag">claim state</div>
        <div class="num">SAFE</div>
        <div class="txt">Readiness layer yes; universal predictor no.</div>
      </div>
    </div>
  </div>

  <figure class="driver-wide" style="margin-top:0">
    <img src="assets/{ASSET_NAME}/Fig_NeoICI_claim_poster.png" alt="NeoICI claim poster" />
    <figcaption><strong>Fig. 0C</strong> <span class="ours-flag">OURS</span> first-screen claim poster. 이 페이지의 핵심 문장을 한 장으로 박아 둔 poster다.</figcaption>
  </figure>
  <figure class="driver-wide">
    <img src="assets/{ASSET_NAME}/Fig_NeoICI_frontpage_collage.png" alt="NeoICI front-page collage" />
    <figcaption><strong>Fig. 0D</strong> <span class="ours-flag">OURS</span> front-page collage. claim, benchmark, replay, case study를 한 번에 보여주는 첫 화면용 합본이다.</figcaption>
  </figure>
  <figure class="driver-wide">
    <img src="assets/{ASSET_NAME}/Fig_NeoICI_super_poster.png" alt="NeoICI super poster" />
    <figcaption><strong>Fig. 0E</strong> <span class="ours-flag">OURS</span> ultra-summary super poster. 첫 화면에서 핵심 판단을 바로 고정하는 최상단 합본이다.</figcaption>
  </figure>
  <figure class="driver-wide">
    <img src="assets/{ASSET_NAME}/Fig_NeoICI_monument_poster.png" alt="NeoICI monument poster" />
    <figcaption><strong>Fig. 0F</strong> <span class="ours-flag">OURS</span> dark monument poster. 첫 화면의 결론을 가장 큰 타이포로 고정한다.</figcaption>
  </figure>

  <header class="hero">
    <div class="kicker">Paper 1-style dossier · Nature Cancer package · OURS</div>
    <h1>NeoICI는 공개 vaccine + ICI cohort에서 combo-therapy readiness를 묻는 OURS 패키지다</h1>
    <div class="subtitle">이 페이지는 NeoICI branch의 웹형 dossier다. local benchmark gauntlet, 공개 combo replay, NeoPrecis-style contract와의 proxy comparison, 그리고 CTMS1_followup4 explanation panel을 한 줄로 이어 붙였다. 핵심은 분명하다. 우리는 universal predictor를 내세우는 것이 아니라 readiness layer를 제시한다.</div>
    <div class="meta">
      <span class="chip hot">KG_GA_evolved AUPRC 0.978</span>
      <span class="chip hot">NeoICI-style proxy rho 0.817</span>
      <span class="chip good">CTMS1_followup4 response_index 0.4187</span>
      <span class="chip ours">OURS: NeoICI branch</span>
      <span class="chip">9 paired public combo samples</span>
      <span class="chip">LOO median 0.833 vs 0.381</span>
      <span class="chip">Claim boundary locked</span>
    </div>
    <div class="statrow">
      <div class="stat"><div class="k">local top score</div><div class="v">0.978</div><div class="t">KG_GA_evolved on the validation-like split.</div></div>
      <div class="stat"><div class="k">open replay rho</div><div class="v">0.817</div><div class="t">NeoICI-style proxy vs response_index.</div></div>
      <div class="stat"><div class="k">case readout</div><div class="v">0.419</div><div class="t">CTMS1_followup4 response_index.</div></div>
      <div class="stat"><div class="k">retained samples</div><div class="v">9</div><div class="t">두 개의 공개 combo cohort, paired GEX/TCR.</div></div>
    </div>
    <div class="impact">
      <div class="panel hot">
        <div class="flag">OURS</div>
        <div class="tag">Bottom line</div>
        <div class="big">Combo-readiness layer</div>
        <div class="body">NeoICI는 universal predictor라고 우기지 않고, vaccine + ICI 맥락에서 실제로 움직이는 immune-state layer를 잡는다.</div>
      </div>
      <div class="panel blue">
        <div class="flag">OURS</div>
        <div class="tag">Local benchmark</div>
        <div class="big">KG_GA_evolved still #1</div>
        <div class="body">validation-like split에서는 KG_GA_evolved가 여전히 local upper bound다. 이 사실을 숨기지 않는 것이 이 dossier의 강점이다.</div>
      </div>
      <div class="panel warn">
        <div class="flag">OURS</div>
        <div class="tag">Open replay</div>
        <div class="big">rho = 0.817</div>
        <div class="body">공개 combo replay에서 NeoICI-style proxy가 response-like index를 더 잘 따라간다. 여기서 방향이 나온다.</div>
      </div>
      <div class="panel good">
        <div class="flag">OURS</div>
        <div class="tag">Case anchor</div>
        <div class="big">CTMS1_followup4</div>
        <div class="body">대표 case는 response_index 0.4187, NeoICI-style rank 1/9. 이 샘플이 전체 서사의 anchor다.</div>
      </div>
    </div>
    <div class="metric-wall">
      <div class="tile hot">
        <div class="label">local top scorer</div>
        <div class="num">0.978</div>
        <div class="desc">KG_GA_evolved still leads the validation-like split. We keep the baseline honest.</div>
      </div>
      <div class="tile blue">
        <div class="label">open replay rho</div>
        <div class="num">0.817</div>
        <div class="desc">NeoICI-style proxy follows response_index better than NeoPrecis-style in the 9-sample replay.</div>
      </div>
      <div class="tile good">
        <div class="label">anchor case rank</div>
        <div class="num">1 / 9</div>
        <div class="desc">CTMS1_followup4 is the clearest lane for explanation and the strongest sample under NeoICI-style ranking.</div>
      </div>
    </div>
    <div class="verdict">
      <div class="headline">Immediate verdict · OURS</div>
      <div class="title">NeoICI는 "이긴다"가 아니라, 공개 combo-replay에서 실제 phenotype을 가장 잘 따라가는 readiness layer로 보인다</div>
      <div class="strap">이 페이지의 임팩트는 단일 숫자에서 나오지 않는다. local benchmark에서 과장하지 않고, 공개 replay에서 방향을 보여주고, 대표 case로 그 방향을 설명하는 세 겹 구조에서 나온다.</div>
      <div class="verdict-grid">
        <div class="box">
          <h5>Local benchmark</h5>
          <div class="num">0.978</div>
          <div class="txt">KG_GA_evolved remains the top local scorer. This keeps the claim honest.</div>
        </div>
        <div class="box">
          <h5>Open replay</h5>
          <div class="num">0.817</div>
          <div class="txt">NeoICI-style proxy tracks response_index better than the NeoPrecis-style proxy in the 9-sample replay.</div>
        </div>
        <div class="box">
          <h5>Anchor case</h5>
          <div class="num">1 / 9</div>
          <div class="txt">CTMS1_followup4 is the cleanest top-lane case under NeoICI-style ranking.</div>
        </div>
      </div>
      <div class="bottom">★ OURS / claim boundary locked: readiness layer yes, universal predictor no, clinical selector no.</div>
    </div>
  </header>

  <div class="layout">
    <aside class="toc">
      <h3>읽는 순서</h3>
      <a href="#intro" style="background:rgba(143,45,37,.08);font-weight:700;color:var(--red);padding-left:6px;border-left:3px solid var(--red)">0. why this page exists</a>
      <a href="#workflow">0b. workflow / reading logic</a>
      <a href="#s1" style="background:rgba(143,45,37,.08);font-weight:700;color:var(--red);padding-left:6px;border-left:3px solid var(--red)">1. 한 줄 결론</a>
      <a href="#s2">2. 알고리즘 gauntlet OURS</a>
      <a href="#s2b">2b. 사람 눈으로 읽는 그림 OURS</a>
      <a href="#s3">3. 공개 combo replay OURS</a>
      <a href="#s4">4. proxy comparison / robustness OURS</a>
      <a href="#s4h">4h. failure lanes / negative controls OURS</a>
      <a href="#s5">5. CTMS1_followup4 case study OURS</a>
      <a href="#s6">6. figure map</a>
      <a href="#s7">7. supplement map</a>
      <a href="#s7b">7b. reviewer-safe decision matrix OURS</a>
      <a href="#s8">8. sources / paths</a>
    </aside>

    <main>
      <section id="intro">
        <h2>0. 왜 이 페이지가 필요한가 <span class="ours-flag">OURS</span></h2>
        <div class="section-note">이 페이지는 front door다. 공개 vaccine + ICI data에서 만든 combo-readiness framing이 response-like phenotype을 published-contract style proxy보다 더 잘 따라가는지 한 번에 보이게 하려는 목적이다.</div>
        <div class="grid g2">
          <div class="card">
            <h4>핵심 배경</h4>
            <p class="small">NeoICI는 universal immunogenicity winner로 말하지 않는다. 우리는 combination therapy에서 진짜 중요한 신호, 즉 immune-state support, T-cell activity, clonotype expansion을 반영하는 prioritization layer로 정의한다. peptide-only나 binding-only logic은 치료 맥락을 놓칠 수 있기 때문이다.</p>
            <p class="small">그래서 공개 combo replay가 핵심 테스트다. 좋은 proxy라면 GSE222011과 GSE255830 전체에서 response-like index와 같이 움직여야지, 한 샘플만 좋아 보여서는 안 된다.</p>
          </div>
          <div class="card">
            <h4>봐야 할 것</h4>
            {block_table(intro_df)}
          </div>
        </div>
        <div class="callout bad starbox">
          <strong>★ 핵심.</strong> 이 페이지가 강한 이유는 같은 방향이 세 번 반복되기 때문이다. local benchmark, open combo replay, sample-level case study가 모두 같은 쪽을 가리켜야 claim이 산다.
        </div>
        {block_table(importance_df)}
      </section>

      <section id="workflow">
        <h2>0b. Workflow / 읽는 논리 <span class="ours-flag">OURS</span></h2>
        <div class="section-note">이 페이지를 결과 뭉치로 읽지 말고 pipeline으로 읽어야 한다. 각 단계는 claim을 좁히고 다음 단계를 더 구체화한다.</div>
        <figure>
          <div class="workflow">
            <div class="step"><div class="n">Step 1</div><div class="t">로컬 benchmark</div><div class="d">NeoICI를 기존 local algorithm과 붙여보고 branch가 과장되지 않았는지 먼저 확인한다.</div></div>
            <div class="step"><div class="n">Step 2</div><div class="t">공개 replay</div><div class="d">open vaccine + ICI cohort로 넘어가 paired GEX/TCR sample에서 phenotype signal이 실제로 보이는지 확인한다.</div></div>
            <div class="step"><div class="n">Step 3</div><div class="t">proxy 대결</div><div class="d">같은 sample을 NeoICI-style과 NeoPrecis-style proxy로 동시에 점수화해 response_index를 따라가는지 본다.</div></div>
            <div class="step"><div class="n">Step 4</div><div class="t">대표 case</div><div class="d">CTMS1_followup4를 열어 왜 top lane인지, clonotype expansion과 state signal이 어떻게 같이 움직이는지 보여준다.</div></div>
          </div>
          <figcaption><strong>Workflow figure.</strong> <span class="ours-flag">OURS</span> 이 페이지는 evidence stack 순서대로 배치했다. benchmark 먼저, replay 다음, proxy comparison 그다음, 마지막에 하나의 sample explanation을 둔다.</figcaption>
        </figure>
        <figure class="driver-wide">
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_evidence_stack.png" alt="NeoICI evidence stack" />
          <figcaption><strong>Fig. 0A</strong> <span class="ours-flag">OURS</span> 한 장짜리 evidence stack infographic. 페이지 전체의 논리 순서를 summary로 박아 둔 그림이다.</figcaption>
        </figure>
        <figure class="driver-wide">
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_verdict_poster.png" alt="NeoICI verdict poster" />
          <figcaption><strong>Fig. 0B</strong> <span class="ours-flag">OURS</span> hero-level verdict poster. 이 페이지의 claim boundary와 핵심 수치를 첫 화면에서 바로 읽게 만든다.</figcaption>
        </figure>
        <div class="grid g3">
          <div class="card">
            <h4>왜 순서가 중요한가</h4>
            <p class="small">local benchmark에서 먼저 무너지면 나머지는 장식이다. replay에서 phenotype이 안 보이면 proxy comparison은 빈 껍데기다. case study가 ranking과 안 맞으면 sample narrative가 무너진다.</p>
          </div>
          <div class="card">
            <h4>무엇이 evidence인가</h4>
            <p class="small">scalar 하나로는 부족하다. AUPRC, Spearman rank behavior, leave-one-out stability, 그리고 one-sample explanation panel이 함께 맞아야 한다.</p>
          </div>
          <div class="card">
            <h4>무엇이 약하게 만드는가</h4>
            <p class="small">proxy 관계가 뒤집히거나, leave-one-out이 불안정하거나, case study가 NeoICI-style ranking 상위권에 들지 못하면 약해진다.</p>
          </div>
        </div>
      </section>

      <section id="s1" style="border:2px solid var(--red);box-shadow:6px 6px 0 rgba(143,45,37,.16)">
        <h2 style="color:var(--red)">1. 한 줄 결론 <span class="ours-flag">OURS</span></h2>
        <div class="section-note">이 dossier는 contract를 분명히 둔다. NeoICI는 universal scorer보다 combo-readiness layer로 더 적합하고, open combo replay에서 phenotype signal이 실제로 걸린다.</div>
        {block_table(executive_small)}
        <div class="callout good">
          <strong>허용 claim.</strong> NeoICI는 open vaccine + ICI cohort와 benchmarking하기 적합하고, 이번 replay에서는 NeoICI-style proxy가 response-like index와 더 잘 맞는다.
        </div>
        <div class="callout warn">
          <strong>경계.</strong> 이걸 clinical selection claim이나 universal superiority claim으로 올리면 안 된다. replay는 작고 proxy는 heuristic이다.
        </div>
      </section>

      <section id="s2">
        <h2>2. 알고리즘 gauntlet <span class="ours-flag">OURS</span></h2>
        <div class="section-note">local benchmark의 메시지는 단순하다. NeoICI는 강하지만 top local scorer는 아니다. 이 차이가 claim framing을 결정한다.</div>
        <div class="callout starbox">
          <strong>★ OURS / NeoICI branch.</strong> 이 섹션은 우리 알고리즘 branch의 성능을 보여주는 곳이다. 상단 헤딩과 figure caption 모두 OURS로 명시했다.
        </div>
        <figure class="driver-wide">
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_algorithm_gauntlet_AUPRC.png" alt="NeoICI algorithm gauntlet AUPRC" />
          <figcaption><strong>Fig. 2</strong> <span class="ours-flag">OURS</span> Local algorithm gauntlet on the validation-like split. KG_GA_evolved remains the local upper bound; NeoICI stays competitive but not top-ranked.</figcaption>
        </figure>
        <div class="grid g2">
          <div class="card">
            <h4>Academic validation-like split</h4>
            {block_table(academic)}
          </div>
          <div class="card">
            <h4>Industrial locked split</h4>
            {block_table(industrial)}
          </div>
        </div>
        <div class="callout bad">
          <strong>주의.</strong> "NeoICI가 모든 기존 알고리즘을 이겼다"는 claim은 아직 안 된다. head-to-head contract가 그걸 지지하지 않는다.
        </div>
      </section>

      <section id="s2b">
        <h2>2b. 알고리즘을 사람 눈으로 읽는 그림 <span class="ours-flag">OURS</span></h2>
        <div class="section-note">이 보조 그림은 숫자 표를 덜 보더라도 순위가 눈에 들어오도록 만든 것이다. reviewer가 가장 먼저 보는 건 결국 '누가 위에 있나'이기 때문이다.</div>
        <div class="grid g3">
          <figure>
            <img src="assets/{ASSET_NAME}/Fig_NeoICI_algorithm_rankbars.png" alt="NeoICI algorithm rank bars" />
            <figcaption><strong>Fig. 2A</strong> <span class="ours-flag">OURS</span> Validation-like split의 top 10 rank bar. NeoICI는 강하지만 KG_GA_evolved가 여전히 맨 위에 있다.</figcaption>
          </figure>
          <figure>
            <img src="assets/{ASSET_NAME}/Fig_NeoICI_algorithm_heatmap.png" alt="NeoICI algorithm heatmap" />
            <figcaption><strong>Fig. 2B</strong> <span class="ours-flag">OURS</span> 같은 알고리즘을 metric heatmap으로 다시 본 그림. 숫자보다 패턴이 먼저 보인다.</figcaption>
          </figure>
          <figure>
            <img src="assets/{ASSET_NAME}/Fig_NeoICI_winloss_heatmap.png" alt="NeoICI win loss heatmap" />
            <figcaption><strong>Fig. 2C</strong> <span class="ours-flag">OURS</span> pairwise win/loss heatmap. NeoICI는 많은 comparator 앞에서는 우세하지만, KG_GA_evolved 같은 local upper bound는 남아 있다.</figcaption>
          </figure>
          <div class="card">
            <h4>이 그림을 보는 법</h4>
            <p class="small">왼쪽 막대는 rank를 인간이 바로 읽게 해준다. 긴 숫자표보다 훨씬 빠르게 <b>우리 알고리즘이 어디쯤 있는지</b> 보여준다.</p>
            <p class="small">색을 두 개만 썼다. 빨강은 <b>OURS</b> branch의 핵심 후보, 파랑은 local upper bound를 뜻한다. 이 대비로 claim이 어디까지 가능한지 바로 보인다.</p>
            <p class="small">이 그림은 “좋다/나쁘다”보다 “top tier 안에 있는가”를 먼저 묻는다. 지금은 top tier 안에는 들어가지만 최정상은 아니다.</p>
          </div>
        </div>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_algorithm_family_bars.png" alt="algorithm family mean bars" />
          <figcaption><strong>Fig. 2D</strong> <span class="ours-flag">OURS</span> algorithm family mean AUPRC. 어떤 family가 전체적으로 강한지 한 번에 보인다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_comparator_winbars.png" alt="comparator winbars" />
          <figcaption><strong>Fig. 2E</strong> <span class="ours-flag">OURS</span> comparator win fraction bars. NeoICI가 누구를 이기고 누구에게 지는지 바로 읽힌다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_algorithm_delta_bars.png" alt="algorithm delta bars" />
          <figcaption><strong>Fig. 2F</strong> <span class="ours-flag">OURS</span> AUPRC delta bars relative to KG_GA_evolved. NeoICI가 기준선에서 얼마나 떨어지는지 한 번에 보인다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_precision_profile.png" alt="precision profile" />
          <figcaption><strong>Fig. 2G</strong> <span class="ours-flag">OURS</span> top-k precision profile. top5, top10, top24 precision을 한 번에 놓고 보면 NeoICI가 어디서 강하고 어디서 약한지 바로 보인다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_split_generalization_heatmap.png" alt="split generalization heatmap" />
          <figcaption><strong>Fig. 2H</strong> <span class="ours-flag">OURS</span> cross-split generalization heatmap. validation-like, industrial-locked, leakage-tolerant split을 같은 눈으로 본다.</figcaption>
        </figure>
      </section>

      <section id="s3">
        <h2>3. 공개 combo replay <span class="ours-flag">OURS</span></h2>
        <div class="section-note">이 replay는 supportive biology layer다. 핵심은 therapeutic efficacy가 아니라 response-like phenotype과 clonotype expansion이 같이 보이느냐이다.</div>
        <div class="callout starbox">
          <strong>★ OURS / 공개 replay.</strong> 이 figure는 공개 vaccine + ICI cohort에서 실제로 response-like state가 보이는지 확인하는 우리 재분석이다.
        </div>
        <figure class="driver-wide">
          <img src="assets/{ASSET_NAME}/Fig_combo_response_index.png" alt="open combo response index" />
          <figcaption><strong>Fig. 3</strong> <span class="ours-flag">OURS</span> Open combo GEO replay across 9 paired samples retained from GSE222011 and GSE255830. CTMS1_followup4 is the strongest response-like lane.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_dataset_comparison_bars.png" alt="dataset comparison bars" />
          <figcaption><strong>Fig. 3B</strong> <span class="ours-flag">OURS</span> dataset-level median comparison. GSE222011과 GSE255830가 response / clone / cytotoxic / NeoICI proxy 축에서 어떻게 다른지 바로 읽힌다.</figcaption>
        </figure>
        <div class="grid g2">
          <div class="card">
            <h4>response-like 상위 샘플</h4>
            {block_table(combo_head[["dataset","sample","response_index","top_clone_frac","NeoICI_style_proxy","NeoPrecis_style_proxy","cytotoxic_score","activation_score"]], max_rows=9)}
          </div>
          <div class="card">
            <h4>replay 요약</h4>
            <p class="small">retained sample 수: <strong>{len(combo)}</strong></p>
            <p class="small">response_index median: <strong>{format_num(combo["response_index"].median(), 3)}</strong></p>
            <p class="small">top sample: <strong>CTMS1_followup4</strong></p>
            <p class="small">top clonotype fraction: <strong>{format_num(combo.loc[combo["sample"]=="CTMS1_followup4","top_clone_frac"].iloc[0], 3)}</strong></p>
            <p class="small">top case의 TCR clonotype 수: <strong>5824</strong></p>
          </div>
        </div>
      </section>

      <section id="s4">
        <h2>4. proxy comparison / robustness <span class="ours-flag">OURS</span></h2>
        <div class="section-note">같은 sample에서 NeoICI vs NeoPrecis-style proxy를 직접 붙인 비교다. 여기서 방향이 갈리면 framing도 갈린다.</div>
        <div class="callout starbox">
          <strong>★ OURS / proxy comparison.</strong> 같은 sample에 대해 두 proxy family를 직접 비교한 구간이다. 여기가 우리 claim의 실질적인 분기점이다.
        </div>
        <figure class="driver-wide">
          <img src="assets/{ASSET_NAME}/Fig_combo_proxy_compare.png" alt="proxy comparison" />
          <figcaption><strong>Fig. 4</strong> <span class="ours-flag">OURS</span> Proxy comparison on the same 9 samples. NeoICI-style proxy tracks the response-like index better than the NeoPrecis-style proxy in this tiny replay.</figcaption>
        </figure>
        <div class="grid g2">
          <div class="card">
            <h4>proxy metrics</h4>
            {block_table(proxy)}
          </div>
          <div class="card">
            <h4>leave-one-out robustness</h4>
            {block_table(loo_small)}
          </div>
        </div>
        <div class="callout good">
          <strong>방향성은 NeoICI 쪽.</strong> NeoICI-style proxy Spearman = 0.817, NeoPrecis-style proxy Spearman = 0.417. leave-one-out median도 NeoICI-style 쪽이 더 높다.
        </div>
        <div class="grid g2">
          <figure>
            <img src="assets/{ASSET_NAME}/Fig_NeoICI_proxy_scatter.png" alt="proxy scatter" />
            <figcaption><strong>Fig. 4A</strong> <span class="ours-flag">OURS</span> 같은 sample을 두 proxy family로 찍은 scatter. CTMS1_followup4가 NeoICI-style 쪽에서 더 선명하게 올라간다.</figcaption>
          </figure>
          <figure>
            <img src="assets/{ASSET_NAME}/Fig_NeoICI_proxy_summary_bars.png" alt="proxy summary bars" />
            <figcaption><strong>Fig. 4B</strong> <span class="ours-flag">OURS</span> proxy family summary bars. NeoICI-style proxy가 response / clone / cytotoxic 축에서 더 낫다.</figcaption>
          </figure>
        </div>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_proxy_loo.png" alt="proxy leave one out" />
          <figcaption><strong>Fig. 4C</strong> <span class="ours-flag">OURS</span> leave-one-out robustness. NeoICI-style proxy가 더 안정적이다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_sample_rank_profile.png" alt="sample rank profile" />
          <figcaption><strong>Fig. 4D</strong> <span class="ours-flag">OURS</span> 9개 replay sample의 rank profile. response rank, NeoICI rank, NeoPrecis rank가 한 번에 보인다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_sample_rank_heatmap.png" alt="sample rank heatmap" />
          <figcaption><strong>Fig. 4E</strong> <span class="ours-flag">OURS</span> sample rank heatmap. response / NeoICI / NeoPrecis rank를 행 단위로 비교한다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_proxy_corr_heatmap.png" alt="proxy correlation heatmap" />
          <figcaption><strong>Fig. 4F</strong> <span class="ours-flag">OURS</span> proxy and state correlation heatmap. 어떤 축이 response_index와 같이 가는지 한 번에 읽힌다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_sample_rank_mismatch_bars.png" alt="sample rank mismatch bars" />
          <figcaption><strong>Fig. 4G</strong> <span class="ours-flag">OURS</span> rank mismatch bars. response rank와 proxy rank의 차이를 샘플별로 보여준다.</figcaption>
        </figure>
        <section id="s4h" style="margin-top:18px;padding-top:10px;border-top:2px solid #c6d4e2">
          <h3 style="margin:0 0 10px;color:#244e73">4h. failure lanes / negative controls <span class="ours-flag">OURS</span></h3>
          <div class="section-note">좋은 논문은 강한 lane만 보여주지 않는다. 실제로 낮은 response lane이 어떻게 생겼는지 같이 보여줘야 claim이 과장되지 않는다.</div>
          <div class="callout bad">
            <strong>왜 필요한가.</strong> Ctms1 같은 positive lane만 반복하면 좋은 그림이지만 좋은 논문은 아니다. low-response lane에서 NeoICI-style proxy가 어디까지 버티고 어디서 흔들리는지 같이 적어야 한다.
          </div>
          <figure class="driver-wide">
            <img src="assets/{ASSET_NAME}/Fig_NeoICI_failure_lanes.png" alt="failure lanes and negative controls" />
            <figcaption><strong>Fig. 4H</strong> <span class="ours-flag">OURS</span> failure lanes / negative controls. top response lane와 low response lane를 나란히 놓아, claim이 어디서 강하고 어디서 약한지 보여준다.</figcaption>
          </figure>
          <div class="grid g2">
            <div class="card">
              <h4>low-response samples</h4>
              {block_table(failure_small)}
            </div>
            <div class="card">
              <h4>해석</h4>
              <p class="small">이 3개 low-lane은 response_index가 낮고, NeoICI-style / NeoPrecis-style proxy도 top-lane처럼 강하게 올라가지 않는다. 즉, 이번 패키지는 positive-only story가 아니라 negative control을 포함한 story다.</p>
              <p class="small">이 부분이 들어가야 reviewer가 “선별된 샘플만 잘 보인 것 아니냐”는 공격을 덜 한다.</p>
            </div>
          </div>
        </section>
      </section>

      <section id="s5">
        <h2>5. CTMS1_followup4 case study <span class="ours-flag">OURS</span></h2>
        <div class="section-note">이 case-study panel이 제일 중요하다. NeoICI-style proxy가 왜 가장 response-like한 lane에서 올라가는지 가장 깨끗하게 설명한다. 아래 heatmap은 진짜 attention이 아니라, replay signal을 사람이 읽기 쉽게 만든 score-coupling map이다.</div>
        <div class="callout starbox">
          <strong>★ OURS / case study.</strong> CTMS1_followup4는 우리 페이지에서 가장 설명력이 높은 sample이다. 이 샘플이 왜 1/9인지 figure와 표를 같이 봐야 한다.
        </div>
        <figure class="driver-wide">
          <img src="assets/{ASSET_NAME}/Fig_CTMS1_followup4_case_study.png" alt="CTMS1 followup4 case study" />
          <figcaption><strong>Fig. 5</strong> <span class="ours-flag">OURS</span> CTMS1_followup4 case study. High response_index co-occurs with high NeoICI proxy rank and substantial T-cell state signal.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_CTMS1_case_profile.png" alt="CTMS1 followup4 case profile" />
          <figcaption><strong>Fig. 5A</strong> <span class="ours-flag">OURS</span> CTMS1_followup4를 cohort median과 직접 비교한 프로필. response / clone / cytotoxic / activation / memory가 한눈에 보인다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_attention_heatmaps.png" alt="attention-style heatmaps" />
          <figcaption><strong>Fig. 5B</strong> <span class="ours-flag">OURS</span> 어텐션이 없어서 만든 attention-style heatmap 5개. 실제 transformer attention은 아니고, replay signal을 사람 눈으로 읽기 쉽게 만든 score-coupling map이다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_sample_heatmap.png" alt="replay sample heatmap" />
          <figcaption><strong>Fig. 5C</strong> <span class="ours-flag">OURS</span> 9개 replay sample을 한 장으로 본 heatmap. 어떤 샘플이 response-like하고 어떤 축이 같이 움직이는지 가장 빠르게 보이는 그림이다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_NeoICI_claim_boundary.png" alt="claim boundary diagram" />
          <figcaption><strong>Fig. 5D</strong> <span class="ours-flag">OURS</span> claim boundary diagram. 허용되는 말과 금지되는 말을 한 번에 고정한다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_CTMS1_case_delta.png" alt="CTMS1 case delta strip" />
          <figcaption><strong>Fig. 5E</strong> <span class="ours-flag">OURS</span> CTMS1_followup4 minus cohort median delta strip. case가 왜 올라가는지 변화량으로 본다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_CTMS1_case_radar.png" alt="CTMS1 case radar" />
          <figcaption><strong>Fig. 5F</strong> <span class="ours-flag">OURS</span> CTMS1_followup4 radar summary. cohort median과 비교했을 때 어떤 축이 같이 살아나는지 요약한다.</figcaption>
        </figure>
        <figure>
          <img src="assets/{ASSET_NAME}/Fig_CTMS1_case_neighbor_bars.png" alt="CTMS1 case neighbor bars" />
          <figcaption><strong>Fig. 5G</strong> <span class="ours-flag">OURS</span> CTMS1_followup4 nearest-neighbor style bar comparison. case를 비슷한 샘플들과 직접 붙여서 why-this-sample을 읽는다.</figcaption>
        </figure>
        {block_table(sample_small[sample_small["sample"]=="CTMS1_followup4"])}
        <div class="grid g3">
          <div class="card"><div class="mini-label">response_index</div><div class="big">{format_num(combo.loc[combo["sample"]=="CTMS1_followup4","response_index"].iloc[0], 3)}</div><div class="small">replay에서 가장 response-like한 lane이다.</div></div>
          <div class="card"><div class="mini-label">NeoICI proxy rank</div><div class="big">1/9</div><div class="small">NeoICI-style proxy에서 가장 높게 뜬 sample이다.</div></div>
          <div class="card"><div class="mini-label">NeoPrecis proxy rank</div><div class="big">4/9</div><div class="small">긍정적이지만 top lane은 아니다.</div></div>
        </div>
      </section>

      <section id="s6">
        <h2>6. figure map <span class="ours-flag">OURS</span></h2>
        <div class="section-note">웹 페이지와 underlying result directory를 연결하는 reviewer-facing bridge다.</div>
        {block_table(figure_df)}
      </section>

      <section id="s7">
        <h2>7. supplement map <span class="ours-flag">OURS</span></h2>
        <div class="section-note">supplementary file은 benchmark contract와 replay contract가 섞이지 않도록 분리했다.</div>
        {block_table(supp_df)}
        <h3>Claim ladder</h3>
        {block_table(claim_small)}
        <section id="s7b" style="margin-top:22px;padding-top:10px;border-top:2px solid #c6d4e2">
          <h3 style="margin:0 0 10px;color:#8f2d25">7b. reviewer-safe decision matrix <span class="ours-flag">OURS</span></h3>
          <div class="section-note">이 표는 논문의 강도를 높이는 동시에 과장을 막는 장치다. 좋은 논문은 claim을 키우는 것보다 경계를 더 정확히 잡는다.</div>
          {block_table(decision_df)}
          <div class="callout warn">
            <strong>핵심.</strong> 이 패키지는 strong positive lane, honest negative control, bounded claim을 동시에 보여준다. 경계가 명확할수록 Nature Cancer용으로 더 강해진다.
          </div>
        </section>
        <h3 style="margin-top:22px">Reproducibility checklist <span class="ours-flag">OURS</span></h3>
        <div class="section-note">이 표는 재분석을 다시 돌릴 때 무엇을 먼저 확인해야 하는지 순서대로 적은 audit trail이다.</div>
        {block_table(reproducibility_df)}
        <h3 style="margin-top:22px">Score definitions / methods ledger <span class="ours-flag">OURS</span></h3>
        <div class="section-note">이 표는 figure에서 쓰는 핵심 score와 proxy의 뜻을 reviewer가 다시 확인할 수 있게 정리한 것이다.</div>
        {block_table(methods_df)}
        <h3 style="margin-top:22px">Provenance / deploy trail <span class="ours-flag">OURS</span></h3>
        <div class="section-note">page, asset, mirror, deploy path를 한 번에 남겨 두면 나중에 재현과 업데이트가 쉬워진다.</div>
        {block_table(provenance_df)}
      </section>

      <section id="s8">
        <h2>8. sources / paths <span class="ours-flag">OURS</span></h2>
        <div class="section-note">이 페이지의 모든 claim은 명시적인 result file에 연결되어 있다. voice-protected manuscript prose는 여기서 만들지 않는다.</div>
        <table>
          <thead><tr><th>artifact</th><th>path</th><th>role</th></tr></thead>
          <tbody>
            <tr><td>Package summary</td><td><code>project/results/p_neo_ici_nature_cancer_package_2026_05_11/NEOICI_NATURE_CANCER_PACKAGE_SUMMARY.md</code></td><td>One-page internal summary.</td></tr>
            <tr><td>Figure map</td><td><code>project/results/p_neo_ici_nature_cancer_package_2026_05_11/NEOICI_NATURE_CANCER_FIGURE_MAP.tsv</code></td><td>Figure routing table.</td></tr>
            <tr><td>Supplementary map</td><td><code>project/results/p_neo_ici_nature_cancer_package_2026_05_11/NEOICI_NATURE_CANCER_SUPPLEMENTARY_MAP.tsv</code></td><td>Supplement routing table.</td></tr>
            <tr><td>Executive table</td><td><code>project/results/p_neo_ici_nature_cancer_package_2026_05_11/NEOICI_NATURE_CANCER_EXECUTIVE_TABLE.tsv</code></td><td>Top-level contract.</td></tr>
            <tr><td>Open combo replay</td><td><code>project/results/p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11/combo_geo_reanalysis.tsv</code></td><td>9-sample replay table.</td></tr>
            <tr><td>Proxy comparison</td><td><code>project/results/p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11/proxy_compare/proxy_compare_metrics.tsv</code></td><td>NeoICI vs NeoPrecis-style comparison.</td></tr>
            <tr><td>CTMS1 case</td><td><code>project/results/p_neo_ici_nature_cancer_synthesis_2026_05_11/ctms1_followup4_case/figures/Fig_CTMS1_followup4_case_study.png</code></td><td>Case-study figure used in Fig. 5.</td></tr>
          </tbody>
        </table>
        <div class="footer">NeoICI Nature Cancer package · paper1-style web dossier · OURS</div>
      </section>
    </main>
  </div>
</div>
</body>
</html>
"""

    html_path = HUB / "neoici_nature_cancer_package.html"
    alias_path = HUB / "neoici_nature_cancer_paper1_style.html"
    write_text(html_path, html)
    write_text(alias_path, html)
    write_text(WWW / "neoici_nature_cancer_package.html", html)
    write_text(WWW / "neoici_nature_cancer_paper1_style.html", html)

    asset_dir = HUB / "assets" / ASSET_NAME
    asset_www = WWW / "assets" / ASSET_NAME
    ensure_dir(asset_dir)
    ensure_dir(asset_www)
    asset_sources = {
        "Fig_NeoICI_impact_banner.png": fig_dir / "Fig_NeoICI_impact_banner.png",
        "Fig_NeoICI_micro_strip.png": fig_dir / "Fig_NeoICI_micro_strip.png",
        "Fig_NeoICI_claim_strip.png": fig_dir / "Fig_NeoICI_claim_strip.png",
        "Fig_NeoICI_topline_banner.png": fig_dir / "Fig_NeoICI_topline_banner.png",
        "Fig_NeoICI_launch_banner.png": fig_dir / "Fig_NeoICI_launch_banner.png",
        "Fig_NeoICI_cover_slab.png": fig_dir / "Fig_NeoICI_cover_slab.png",
        "Fig_NeoICI_verdict_wall.png": fig_dir / "Fig_NeoICI_verdict_wall.png",
        "Fig_NeoICI_opening_impact_board.png": fig_dir / "Fig_NeoICI_opening_impact_board.png",
        "Fig_NeoICI_evidence_stack.png": fig_dir / "Fig_NeoICI_evidence_stack.png",
        "Fig_NeoICI_verdict_poster.png": fig_dir / "Fig_NeoICI_verdict_poster.png",
        "Fig_NeoICI_claim_poster.png": fig_dir / "Fig_NeoICI_claim_poster.png",
        "Fig_NeoICI_frontpage_collage.png": fig_dir / "Fig_NeoICI_frontpage_collage.png",
        "Fig_NeoICI_super_poster.png": fig_dir / "Fig_NeoICI_super_poster.png",
        "Fig_NeoICI_monument_poster.png": fig_dir / "Fig_NeoICI_monument_poster.png",
        "Fig_NeoICI_cover_wall.png": fig_dir / "Fig_NeoICI_cover_wall.png",
        "Fig_NeoICI_combo_dataset_ladder.png": ROOT / "project" / "results" / "p_neo_ici_combo_2026_05_11" / "figures" / "Fig_NeoICI_combo_dataset_ladder.png",
        "Fig_NeoICI_algorithm_gauntlet_AUPRC.png": ROOT / "project" / "results" / "p_neo_ici_combo_2026_05_11" / "figures" / "Fig_NeoICI_algorithm_gauntlet_AUPRC.png",
        "Fig_NeoICI_algorithm_rankbars.png": fig_dir / "Fig_NeoICI_algorithm_rankbars.png",
        "Fig_NeoICI_algorithm_heatmap.png": fig_dir / "Fig_NeoICI_algorithm_heatmap.png",
        "Fig_NeoICI_precision_profile.png": fig_dir / "Fig_NeoICI_precision_profile.png",
        "Fig_NeoICI_split_generalization_heatmap.png": fig_dir / "Fig_NeoICI_split_generalization_heatmap.png",
        "Fig_combo_response_index.png": ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "figures" / "Fig_combo_response_index.png",
        "Fig_NeoICI_dataset_comparison_bars.png": fig_dir / "Fig_NeoICI_dataset_comparison_bars.png",
        "Fig_NeoICI_proxy_scatter.png": fig_dir / "Fig_NeoICI_proxy_scatter.png",
        "Fig_combo_proxy_compare.png": ROOT / "project" / "results" / "p_neoprecis_vs_neoici_combo_reanalysis_2026_05_11" / "proxy_compare" / "figures" / "Fig_combo_proxy_compare.png",
        "Fig_NeoICI_proxy_summary_bars.png": fig_dir / "Fig_NeoICI_proxy_summary_bars.png",
        "Fig_NeoICI_proxy_loo.png": fig_dir / "Fig_NeoICI_proxy_loo.png",
        "Fig_CTMS1_followup4_case_study.png": ROOT / "project" / "results" / "p_neo_ici_nature_cancer_synthesis_2026_05_11" / "ctms1_followup4_case" / "figures" / "Fig_CTMS1_followup4_case_study.png",
        "Fig_CTMS1_case_profile.png": fig_dir / "Fig_CTMS1_case_profile.png",
        "Fig_NeoICI_attention_heatmaps.png": fig_dir / "Fig_NeoICI_attention_heatmaps.png",
        "Fig_NeoICI_sample_heatmap.png": fig_dir / "Fig_NeoICI_sample_heatmap.png",
        "Fig_NeoICI_winloss_heatmap.png": fig_dir / "Fig_NeoICI_winloss_heatmap.png",
        "Fig_NeoICI_algorithm_family_bars.png": fig_dir / "Fig_NeoICI_algorithm_family_bars.png",
        "Fig_NeoICI_comparator_winbars.png": fig_dir / "Fig_NeoICI_comparator_winbars.png",
        "Fig_NeoICI_algorithm_delta_bars.png": fig_dir / "Fig_NeoICI_algorithm_delta_bars.png",
        "Fig_NeoICI_sample_rank_profile.png": fig_dir / "Fig_NeoICI_sample_rank_profile.png",
        "Fig_NeoICI_sample_rank_heatmap.png": fig_dir / "Fig_NeoICI_sample_rank_heatmap.png",
        "Fig_NeoICI_proxy_corr_heatmap.png": fig_dir / "Fig_NeoICI_proxy_corr_heatmap.png",
        "Fig_NeoICI_sample_rank_mismatch_bars.png": fig_dir / "Fig_NeoICI_sample_rank_mismatch_bars.png",
        "Fig_NeoICI_failure_lanes.png": fig_dir / "Fig_NeoICI_failure_lanes.png",
        "Fig_NeoICI_claim_boundary.png": fig_dir / "Fig_NeoICI_claim_boundary.png",
        "Fig_CTMS1_case_delta.png": fig_dir / "Fig_CTMS1_case_delta.png",
        "Fig_CTMS1_case_radar.png": fig_dir / "Fig_CTMS1_case_radar.png",
        "Fig_CTMS1_case_neighbor_bars.png": fig_dir / "Fig_CTMS1_case_neighbor_bars.png",
    }
    for name, src in asset_sources.items():
        copy_bytes(src, asset_dir / name)
        copy_bytes(src, asset_www / name)

    meta = {
        "html": str(html_path),
        "alias_html": str(alias_path),
        "summary": str(OUT / "NEOICI_NATURE_CANCER_PACKAGE_SUMMARY.md"),
        "figure_map": str(OUT / "NEOICI_NATURE_CANCER_FIGURE_MAP.tsv"),
        "supplementary_map": str(OUT / "NEOICI_NATURE_CANCER_SUPPLEMENTARY_MAP.tsv"),
        "executive_table": str(OUT / "NEOICI_NATURE_CANCER_EXECUTIVE_TABLE.tsv"),
    }
    write_text(OUT / "neoici_nature_cancer_package.json", json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
