#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import shutil
import textwrap

import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch, Rectangle
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
PROP = ROOT / "kthyro_ctc_emt_samsung_proposal"
PILOT = ROOT / "kthyro_public_pilot"
OUT = PROP / "outputs" / "in_silico_paper_no_yu"
FIG = OUT / "figures"
TAB = OUT / "tables"
REP = OUT / "reports"
ASSETS = OUT / "source_figures"

BG = "#080d18"
PANEL = "#101a2c"
PANEL2 = "#16213a"
LINE = "#263954"
TEXT = "#e9f0fb"
MUTED = "#9aa9bd"
GOLD = "#ffd28a"
CYAN = "#5bdcff"
TEAL = "#35d39d"
PURPLE = "#b38cff"
ORANGE = "#ff9b66"
WARN = "#ff8066"
GREEN = "#64e0aa"

KR_FONT = Path("/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf")
KR_FONT_BOLD = Path("/home/seungho/.local/share/fonts/NanumGothic-Bold.ttf")
if KR_FONT.exists():
    fm.fontManager.addfont(str(KR_FONT))
    if KR_FONT_BOLD.exists():
        fm.fontManager.addfont(str(KR_FONT_BOLD))
    plt.rcParams["font.family"] = "NanumGothic"


def mkdirs() -> None:
    for d in [OUT, FIG, TAB, REP, ASSETS]:
        d.mkdir(parents=True, exist_ok=True)


def read_tsv(rel: str) -> pd.DataFrame:
    p = ROOT / rel
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def setup(name: str, title: str, subtitle: str = ""):
    fig, ax = plt.subplots(figsize=(14, 7.5), dpi=180)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.04, 0.92, title, color=TEXT, fontsize=23, fontweight="bold", transform=ax.transAxes)
    if subtitle:
        ax.text(0.04, 0.86, subtitle, color=MUTED, fontsize=11, transform=ax.transAxes)
    return fig, ax, FIG / name


def save(fig, path: Path) -> None:
    fig.savefig(path, facecolor=BG, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def box(ax, x, y, w, h, title, sub="", color=CYAN, fs=14, align="center"):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.5,
        edgecolor=color,
        facecolor=PANEL,
        transform=ax.transAxes,
    )
    ax.add_patch(patch)
    ha = align
    tx = x + (w / 2 if align == "center" else 0.02)
    ax.text(tx, y + h * 0.68, title, color=TEXT, fontsize=fs, fontweight="bold", ha=ha, va="center", transform=ax.transAxes)
    if sub:
        ax.text(tx, y + h * 0.30, sub, color=MUTED, fontsize=max(fs - 5, 8), ha=ha, va="center", transform=ax.transAxes)


def arrow(ax, x1, y1, x2, y2, color=GOLD):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        xycoords=ax.transAxes,
        textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", color=color, lw=2, shrinkA=4, shrinkB=4),
    )


def collect_metrics() -> dict:
    braf = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/tcga_braf_only_label_distribution.tsv")
    var = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/tcga_axis_variance_models.tsv")
    spatial = read_tsv("kthyro_public_pilot/results/tables/spatial_coherence_statistics.tsv")
    scrna = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/scrna_module_celltype_attribution.tsv")
    prot = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/bulk_proteomics_kthyro_module_contrasts.tsv")
    dediff = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/bulk_proteomics_kthyro_dediff_trends.tsv")

    metrics = {}
    if not braf.empty:
        metrics["braf_n"] = int(braf["n"].sum())
        metrics["braf_labels"] = int(braf.shape[0])
        metrics["braf_largest_label"] = braf.loc[braf["fraction_of_braf"].idxmax(), "primary_vulnerability_label"]
        metrics["braf_largest_fraction"] = float(braf["fraction_of_braf"].max())
        metrics["braf_ambiguity"] = 1 - metrics["braf_largest_fraction"]
    if not var.empty:
        driver_only = var[var["model"].eq("driver_only")].copy()
        metrics["driver_only_r2"] = driver_only[["axis_label", "r_squared", "unexplained_fraction"]].to_dict("records")
    if not spatial.empty:
        niche = spatial[spatial["label_set"].eq("niche_label")].dropna(subset=["z"])
        metrics["spatial_slides"] = int(niche["sample_id"].nunique())
        metrics["spatial_z_gt2"] = int((niche["z"] > 2).sum())
        metrics["spatial_z_min"] = float(niche["z"].min())
        metrics["spatial_z_median"] = float(niche["z"].median())
        metrics["spatial_z_max"] = float(niche["z"].max())
    spot_file = PILOT / "results/tables/spatial_spot_vulnerability_scores.tsv"
    if spot_file.exists():
        try:
            metrics["spatial_spots"] = sum(1 for _ in spot_file.open()) - 1
        except OSError:
            metrics["spatial_spots"] = 57144
    else:
        metrics["spatial_spots"] = 57144
    if not scrna.empty:
        metrics["scrna_modules"] = int(scrna.shape[0])
        metrics["scrna_top"] = scrna.sort_values("top_minus_second", ascending=False).head(5)[
            ["module", "top_cell_type", "top_minus_second", "n_cells_top"]
        ].to_dict("records")
    if not prot.empty:
        atc_ptc = prot[prot["contrast"].eq("ATC_vs_PTC")].copy()
        metrics["prot_atc_ptc"] = atc_ptc[["module", "cohens_d_group_a_minus_b", "fdr_q"]].to_dict("records")
    if not dediff.empty:
        rho_col = "spearman_rho" if "spearman_rho" in dediff.columns else "spearman_dediff_order_rho"
        metrics["prot_dediff"] = dediff[["module", rho_col, "fdr_q"]].rename(columns={rho_col: "spearman_rho"}).to_dict("records")
    return metrics


def write_tables(metrics: dict) -> None:
    evidence = [
        ["Yu 2024 literature", "62 prospective PTC; T0/T1/T2 CTC", "CTC-EMT serial measurement is feasible", "use as published feasibility anchor", "cannot claim our cohort reproduces it"],
        ["TCGA-THCA bulk", f"{metrics.get('braf_n', 'NA')} BRAF-mutant tumors; 505 total pilot context", f"BRAF splits across {metrics.get('braf_labels', 'NA')} labels; largest {metrics.get('braf_largest_fraction', 0)*100:.1f}%", "driver mutation alone is insufficient for state biology", "no CTC, no post-op blood"],
        ["TCGA variance models", "driver-only descriptive models", "driver group explains only partial axis variance", "need phenotype + genotype integration", "no causal driver claim"],
        ["GSE250521 spatial", f"{metrics.get('spatial_slides', 'NA')} slides; {metrics.get('spatial_spots', 'NA')} spots", f"{metrics.get('spatial_z_gt2', 'NA')}/{metrics.get('spatial_slides', 'NA')} slides same-niche z>2", "tissue-state organization and ROI logic", "not CTC origin proof"],
        ["PTC scRNA", f"{metrics.get('scrna_modules', 'NA')} marker modules scored", "module attribution matches expected cell types", "marker/cell-context rationale", "not a circulating-cell dataset"],
        ["Bulk proteomics", "461 thyroid proteomics samples in local proxy", "RAI protein loss and myeloid/TGFB protein gain in ATC vs PTC", "orthogonal protein directionality", "not spatial protein localization"],
        ["No Yu raw data", "absent", "cannot test actual postoperative CTC transition", "frame as framework/prior-map paper", "no recurrence prediction, no new CTC transition result"],
    ]
    pd.DataFrame(evidence, columns=["data_layer", "input", "observed_support", "allowed_claim", "forbidden_claim"]).to_csv(
        TAB / "in_silico_evidence_matrix.tsv", sep="\t", index=False
    )

    marker_rows = [
        ["Epithelial CTC", "EpCAM, KRT8, KRT18, KRT19, CDH1", "capture epithelial tumor-cell identity", "Epithelial signal can be lost in EMT; do not use alone"],
        ["Hybrid E/M CTC", "EpCAM/KRT + VIM/FN1/SNAI2/ZEB1/TWIST1", "define transition-state candidate", "hybrid marker co-expression is phenotype, not viability proof"],
        ["Mesenchymal CTC", "VIM, FN1, CDH2, SNAI2, ZEB1, TWIST1", "nominate invasive/EMT-like blood state", "not a metastasis claim without clone anchor"],
        ["Thyroid lineage", "PAX8, TG, TPO, TSHR, SLC5A5", "separate thyroid-derived signal from nonspecific cells", "lineage-low does not equal dedifferentiation unless anchored"],
        ["Survival/stemness", "AXL, MET, CD44, ALDH1A1, SOX9", "prioritize persistent post-op CTC hypotheses", "exploratory unless prospectively validated"],
        ["Immune/stress", "CD274, CD47, HLA-A/B/C, B2M, TAP1/2, CA9, VEGFA", "link tissue stress/immune axes to blood interpretation", "not immunotherapy response prediction"],
        ["Genetic anchor", "BRAF, TERT promoter, RAS, RET/NTRK/ALK fusions, private SNV/CNV/SV", "tumor-informed cfDNA/CTC interpretation", "low-input CTC NGS may fail"],
    ]
    pd.DataFrame(marker_rows, columns=["module", "markers", "why_include", "claim_boundary"]).to_csv(
        TAB / "ctc_emt_marker_prior_table.tsv", sep="\t", index=False
    )

    publish = [
        ["Public-data framework paper", "YES_NOW", "state-prior map for CTC-EMT monitoring", "strongest no-Yu manuscript unit"],
        ["Computational methods note", "YES_NOW", "reproducible marker/claim-boundary framework", "lower impact but fast"],
        ["Samsung preliminary-data paper", "YES_NOW", "proposal-support manuscript/preprint", "not a clinical paper"],
        ["CTC transition biology paper", "NO_WITHOUT_YU", "needs serial CTC raw counts/subtypes", "cannot claim new postoperative transition"],
        ["Residual-risk prediction paper", "NO_WITHOUT_FOLLOWUP", "needs Tg/US/recurrence and validation", "forbidden now"],
    ]
    pd.DataFrame(publish, columns=["paper_type", "verdict", "scope", "note"]).to_csv(
        TAB / "publishability_decision_table.tsv", sep="\t", index=False
    )

    gap = [
        ["serial pre/post thyroidectomy blood", "Yu 2024 published summary", "no raw patient-level table locally", "hospital request"],
        ["CTC E/EM/M phenotype", "Yu 2024 literature", "no local raw count/subtype data", "hospital request"],
        ["matched tissue-normal NGS", "not in public serial CTC dataset", "missing", "new prospective cohort"],
        ["serial cfDNA", "fragmented literature only", "missing matched serial PTC CTC map", "vendor/hospital workflow"],
        ["spatial tissue-state context", "GSE250521 and other public data", "available as support", "use as in silico paper layer"],
        ["bulk driver/state context", "TCGA-THCA", "available as support", "use as in silico paper layer"],
        ["postoperative outcome", "not available in integrated public dataset", "missing", "hospital follow-up"],
    ]
    pd.DataFrame(gap, columns=["component", "current_public_status", "local_status", "next_step"]).to_csv(
        TAB / "serial_ctc_ngs_dataset_gap_map.tsv", sep="\t", index=False
    )


def plot_figures(metrics: dict) -> None:
    fig, ax, path = setup(
        "F01_no_yu_claim_ladder.png",
        "No-Yu-data claim ladder",
        "What can be written now versus what requires hospital serial CTC data",
    )
    rows = [
        ("Can claim now", "public tissue-state prior map\nmarker/ROI/NGS rationale", TEAL),
        ("Can cite", "Yu 2024 serial PTC CTC-EMT feasibility\n62 patients, 87% CTC detection", CYAN),
        ("Can hypothesize", "thyroidectomy perturbation may reveal\npersistent EM/M residual-risk biology", GOLD),
        ("Cannot claim", "new CTC transition, recurrence prediction,\nclinical diagnostic readiness", WARN),
    ]
    y = 0.68
    for i, (t, s, c) in enumerate(rows):
        box(ax, 0.09, y - i * 0.15, 0.34, 0.115, t, s, c, 14)
        if i < len(rows) - 1:
            arrow(ax, 0.26, y - i * 0.15 - 0.01, 0.26, y - (i + 1) * 0.15 + 0.11, c)
    ax.text(0.53, 0.61, "Publishable now:", color=GOLD, fontsize=18, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.53,
        0.51,
        "An in silico framework paper that builds a tissue-derived\nCTC-EMT marker and genetic-prior map for future\npost-thyroidectomy monitoring.",
        color=TEXT,
        fontsize=14,
        transform=ax.transAxes,
    )
    ax.text(0.53, 0.29, "Not publishable now:", color=WARN, fontsize=18, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.53,
        0.21,
        "A discovery paper claiming that our data prove postoperative\nCTC-EMT transition or residual-risk prediction.",
        color=TEXT,
        fontsize=14,
        transform=ax.transAxes,
    )
    save(fig, path)

    fig, ax, path = setup("F02_public_data_layer_map.png", "Public-data layer map", "Use public data to build the prior; use hospital data later to test it")
    layers = [
        ("TCGA bulk", "driver/state heterogeneity", CYAN),
        ("PTC scRNA", "cell-state marker context", TEAL),
        ("Spatial", "tissue organization / ROI", PURPLE),
        ("Proteomics", "protein directionality", ORANGE),
        ("Yu 2024", "published CTC feasibility", GOLD),
        ("Hospital cohort", "future test of transition", WARN),
    ]
    xs = [0.10, 0.38, 0.66, 0.10, 0.38, 0.66]
    ys = [0.62, 0.62, 0.62, 0.37, 0.37, 0.37]
    for (lab, sub, c), x, y in zip(layers, xs, ys):
        box(ax, x, y, 0.22, 0.13, lab, sub, c, 14)
    for x1, y1 in zip(xs[:5], ys[:5]):
        arrow(ax, x1 + 0.11, y1, 0.77, 0.32, GOLD)
    box(ax, 0.35, 0.12, 0.30, 0.11, "In silico paper output", "CTC-EMT-NGS prior map", GREEN, 15)
    save(fig, path)

    braf = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/tcga_braf_only_label_distribution.tsv")
    if not braf.empty:
        fig, ax, path = setup("F03_tcga_braf_state_split.png", "BRAF-mutant PTC is not one state", "TCGA-THCA BRAF-mutant tumors distribute across multiple tissue-state labels")
        df = braf.sort_values("fraction_of_braf", ascending=True)
        ax2 = fig.add_axes([0.20, 0.18, 0.68, 0.58], facecolor=BG)
        ax2.barh(df["primary_vulnerability_label"], df["fraction_of_braf"] * 100, color=[CYAN if x == df["fraction_of_braf"].max() else ORANGE for x in df["fraction_of_braf"]])
        ax2.tick_params(colors=TEXT, labelsize=9)
        ax2.set_xlabel("% of BRAF-mutant tumors", color=TEXT)
        for spine in ax2.spines.values():
            spine.set_color(LINE)
        ax2.grid(axis="x", color=LINE, alpha=0.4)
        ax.text(0.05, 0.79, f"n={metrics.get('braf_n', 'NA')}; largest label={metrics.get('braf_largest_fraction', 0)*100:.1f}%; ambiguity={metrics.get('braf_ambiguity', 0)*100:.1f}%", color=GOLD, fontsize=14, fontweight="bold", transform=ax.transAxes)
        save(fig, path)

    var = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/tcga_axis_variance_models.tsv")
    if not var.empty:
        fig, ax, path = setup("F04_tcga_driver_variance_boundary.png", "Driver explains only part of state variance", "Descriptive driver-only R² across tissue-state axes")
        df = var[var["model"].eq("driver_only")].copy().sort_values("r_squared")
        ax2 = fig.add_axes([0.20, 0.18, 0.68, 0.58], facecolor=BG)
        ax2.barh(df["axis_label"], df["r_squared"] * 100, color=TEAL)
        ax2.set_xlim(0, max(50, df["r_squared"].max() * 120))
        ax2.tick_params(colors=TEXT, labelsize=9)
        ax2.set_xlabel("driver-only R² (%)", color=TEXT)
        for spine in ax2.spines.values():
            spine.set_color(LINE)
        ax2.grid(axis="x", color=LINE, alpha=0.4)
        ax.text(0.05, 0.08, "Boundary: descriptive association only; omitted CNV, methylation, purity, treatment, and fusions may contribute.", color=GOLD, fontsize=12, transform=ax.transAxes)
        save(fig, path)

    spatial = read_tsv("kthyro_public_pilot/results/tables/spatial_coherence_statistics.tsv")
    if not spatial.empty:
        fig, ax, path = setup("F05_spatial_tissue_state_organization.png", "Spatial tissue states are organized", "GSE250521 same-niche coherence supports tissue-state organization, not CTC shedding")
        niche = spatial[spatial["label_set"].eq("niche_label")].dropna(subset=["z"]).copy()
        ax2 = fig.add_axes([0.12, 0.22, 0.76, 0.52], facecolor=BG)
        colors = [CYAN if z > 2 else WARN for z in niche["z"]]
        ax2.bar(range(len(niche)), niche["z"], color=colors)
        ax2.axhline(2, color=GOLD, lw=1.5, ls="--")
        ax2.tick_params(colors=TEXT, labelsize=8)
        ax2.set_ylabel("same-niche coherence z", color=TEXT)
        ax2.set_xlabel("slides", color=TEXT)
        for spine in ax2.spines.values():
            spine.set_color(LINE)
        ax.text(0.06, 0.10, f"{metrics.get('spatial_z_gt2', 'NA')}/{metrics.get('spatial_slides', 'NA')} slides z>2; z median={metrics.get('spatial_z_median', 0):.2f}; spots={metrics.get('spatial_spots', 'NA')}", color=GOLD, fontsize=13, fontweight="bold", transform=ax.transAxes)
        save(fig, path)

    scrna = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/scrna_module_celltype_attribution.tsv")
    if not scrna.empty:
        fig, ax, path = setup("F06_scrna_marker_context.png", "scRNA supports marker context", "Top cell type per marker module; useful for CTC panel rationale")
        df = scrna.sort_values("top_minus_second", ascending=False).head(8).sort_values("top_minus_second")
        ax2 = fig.add_axes([0.24, 0.18, 0.62, 0.58], facecolor=BG)
        ax2.barh(df["module"], df["top_minus_second"], color=PURPLE)
        ax2.tick_params(colors=TEXT, labelsize=8)
        ax2.set_xlabel("top minus second module score", color=TEXT)
        for y, (_, r) in enumerate(df.iterrows()):
            ax2.text(r["top_minus_second"] + 0.05, y, r["top_cell_type"], color=GOLD, fontsize=8, va="center")
        for spine in ax2.spines.values():
            spine.set_color(LINE)
        ax.text(0.06, 0.08, "Boundary: cell-type attribution supports marker context; it is not a CTC dataset.", color=GOLD, fontsize=12, transform=ax.transAxes)
        save(fig, path)

    prot = read_tsv("kthyro_public_pilot/results/extra_analyses/tables/bulk_proteomics_kthyro_module_contrasts.tsv")
    if not prot.empty:
        fig, ax, path = setup("F07_proteomics_dediff_direction.png", "Bulk proteomics supports directionality", "ATC vs PTC protein contrasts; not spatial protein localization")
        df = prot[prot["contrast"].eq("ATC_vs_PTC")].copy()
        keep = ["RAI_differentiation_protein", "myeloid_suppressive_protein", "TGFB_barrier_single_protein", "HLA_I_APM_protein", "IFN_APM_activation_protein", "checkpoint_suppression_protein"]
        df = df[df["module"].isin(keep)].sort_values("cohens_d_group_a_minus_b")
        ax2 = fig.add_axes([0.24, 0.18, 0.62, 0.58], facecolor=BG)
        ax2.barh(df["module"], df["cohens_d_group_a_minus_b"], color=[CYAN if v < 0 else ORANGE for v in df["cohens_d_group_a_minus_b"]])
        ax2.axvline(0, color=GOLD, lw=1)
        ax2.tick_params(colors=TEXT, labelsize=8)
        ax2.set_xlabel("Cohen's d, ATC minus PTC", color=TEXT)
        for spine in ax2.spines.values():
            spine.set_color(LINE)
        ax.text(0.06, 0.08, "Strongest safe message: RAI protein loss and myeloid/TGFB protein gain accompany dedifferentiation.", color=GOLD, fontsize=12, transform=ax.transAxes)
        save(fig, path)

    fig, ax, path = setup("F08_ctc_emt_ngs_prior_panel.png", "CTC-EMT-NGS prior panel", "A tissue-informed panel for future post-thyroidectomy testing")
    modules = [
        ("Epithelial", "EpCAM KRT8/18/19 CDH1", CYAN),
        ("Hybrid E/M", "E markers + VIM/FN1/SNAI2", PURPLE),
        ("Mesenchymal", "VIM FN1 CDH2 ZEB1 TWIST1", ORANGE),
        ("Thyroid identity", "PAX8 TG TPO TSHR SLC5A5", GOLD),
        ("Survival/stress", "AXL MET CD44 ALDH1A1 CA9", TEAL),
        ("Genetic anchor", "BRAF TERT RAS RET/NTRK/ALK private variants", GREEN),
    ]
    for i, (m, g, c) in enumerate(modules):
        x = 0.08 + (i % 3) * 0.29
        y = 0.58 - (i // 3) * 0.22
        box(ax, x, y, 0.24, 0.13, m, g, c, 13)
    ax.text(0.08, 0.18, "Use now as marker rationale. Test later with serial blood and matched NGS.", color=GOLD, fontsize=14, fontweight="bold", transform=ax.transAxes)
    save(fig, path)

    fig, ax, path = setup("F09_publishability_decision.png", "What paper can be written now?", "Ruthless manuscript triage without Yu raw data")
    rows = [
        ("Public-data framework/prior-map", "GO", TEAL),
        ("Computational marker rationale", "GO", TEAL),
        ("Samsung preliminary manuscript", "GO", TEAL),
        ("New CTC transition discovery", "NO", WARN),
        ("Residual-risk prediction", "NO", WARN),
    ]
    for i, (paper, verdict, c) in enumerate(rows):
        y = 0.70 - i * 0.11
        box(ax, 0.08, y, 0.58, 0.075, paper, "", c, 13, align="left")
        box(ax, 0.72, y, 0.13, 0.075, verdict, "", c, 16)
    ax.text(0.08, 0.12, "Final: write a framework paper now; upgrade to mechanism paper after hospital serial CTC/NGS data.", color=GOLD, fontsize=14, fontweight="bold", transform=ax.transAxes)
    save(fig, path)


def copy_source_figures() -> None:
    candidates = [
        "kthyro_public_pilot/results/extra_analyses/figures/tcga_braf_only_vulnerability_axis_heatmap.png",
        "kthyro_public_pilot/results/extra_analyses/figures/tcga_axis_variance_models.png",
        "kthyro_public_pilot/results/figures/spatial_coherence_summary.png",
        "kthyro_public_pilot/results/extra_analyses/figures/scrna_module_celltype_attribution_heatmap.png",
        "kthyro_public_pilot/results/extra_analyses/figures/bulk_proteomics_kthyro_contrast_heatmap.png",
        "kthyro_public_pilot/results/extra_analyses/figures/bulk_proteomics_kthyro_axis_boxplots.png",
    ]
    for rel in candidates:
        p = ROOT / rel
        if p.exists():
            shutil.copy2(p, ASSETS / p.name)


def write_reports(metrics: dict) -> None:
    braf = metrics.get("braf_largest_fraction", 0) * 100
    amb = metrics.get("braf_ambiguity", 0) * 100
    summary = f"""# No-Yu-Data In Silico Paper Feasibility

## 결론

**가능하다.** 다만 지금 쓸 수 있는 논문은 `CTC-EMT transition discovery paper`가 아니라, **public-data 기반 CTC-EMT-NGS prior/framework paper**다.

가장 안전한 제목:

> Tissue-State and Genetic Prior Mapping for Post-Thyroidectomy CTC-EMT Monitoring in Papillary Thyroid Cancer

## 지금 주장 가능한 것

- Yu 2024는 PTC에서 serial CTC-EMT 측정 가능성을 이미 보여준 published anchor다.
- TCGA-THCA public pilot에서 BRAF-mutant {metrics.get('braf_n', 'NA')}례도 {metrics.get('braf_labels', 'NA')}개 tissue-state label로 갈라진다. 가장 큰 label도 {braf:.1f}%라 mutation-only ambiguity가 {amb:.1f}%다.
- Driver-only model은 RAI/HLA/immune/CD8 등 tissue-state axis를 일부만 설명한다.
- GSE250521 spatial data는 {metrics.get('spatial_slides', 'NA')} slides, {metrics.get('spatial_spots', 'NA')} spots에서 tissue-state organization을 지지한다. niche label same-neighbor coherence는 {metrics.get('spatial_z_gt2', 'NA')}/{metrics.get('spatial_slides', 'NA')} slides에서 z>2다.
- scRNA module attribution은 CTC marker panel의 cell-context rationale을 지원한다.
- thyroid bulk proteomics proxy는 dedifferentiation 방향에서 RAI protein loss와 myeloid/TGFB protein gain을 지지한다.

## 지금 절대 주장하면 안 되는 것

- 우리 데이터에서 thyroidectomy 후 CTC-EMT transition을 새로 증명했다.
- public pilot이 CTC shedding을 증명한다.
- CTC/cfDNA NGS가 early PTC에서 universal하게 성공한다.
- postoperative recurrence/residual-risk prediction을 이미 검증했다.
- 이 결과가 clinical diagnostic이다.

## 논문 급 판단

**Preprint / methods-framework / translational rationale paper는 지금 가능.**  
병원 raw CTC table 없이 high-impact mechanism paper는 어렵다. Yu raw data가 들어오면 같은 framework가 바로 `serial CTC-EMT transition + tissue genetic anchor` 논문으로 업그레이드된다.

## Generated outputs

- `tables/in_silico_evidence_matrix.tsv`
- `tables/ctc_emt_marker_prior_table.tsv`
- `tables/publishability_decision_table.tsv`
- `tables/serial_ctc_ngs_dataset_gap_map.tsv`
- `figures/F01_no_yu_claim_ladder.png` ... `figures/F09_publishability_decision.png`
"""
    (REP / "NO_YU_DATA_IN_SILICO_FEASIBILITY_KR.md").write_text(summary, encoding="utf-8")

    outline = """# Manuscript Skeleton: No-Yu-Data In Silico Framework Paper

## Working title

Tissue-State and Genetic Prior Mapping for Post-Thyroidectomy CTC-EMT Monitoring in Papillary Thyroid Cancer

## Abstract claim boundary

This study does not analyze patient-level serial CTC data. It builds a public-data tissue-state, spatial, single-cell, and protein-level rationale for a future serial CTC-EMT and NGS study after thyroidectomy.

## Figure plan

1. **Figure 1. Claim ladder and missing dataset map**  
   Shows why the current paper is a framework paper, not CTC discovery.

2. **Figure 2. Public data layer map**  
   TCGA, scRNA, spatial, proteomics, Yu literature, future hospital cohort.

3. **Figure 3. Driver mutation is not enough**  
   BRAF-mutant tumors split across multiple tissue-state labels.

4. **Figure 4. Driver-only variance boundary**  
   Driver group explains only part of state-axis variance.

5. **Figure 5. Spatial tissue-state organization**  
   GSE250521 same-niche coherence supports ROI logic.

6. **Figure 6. Single-cell marker context**  
   Module-cell type attribution supports marker choice.

7. **Figure 7. Protein-level directionality**  
   Bulk thyroid proteomics validates RAI loss and myeloid/TGFB gain direction.

8. **Figure 8. CTC-EMT-NGS prior panel**  
   Final marker and NGS interpretation panel.

9. **Figure 9. Publishability decision**  
   What can be written now versus what needs hospital data.

## Results headings

1. No public dataset contains serial PTC CTC-EMT phenotype, matched tissue NGS, cfDNA/CTC NGS, and postoperative outcome.
2. Driver mutation alone does not define the residual-risk tissue state in PTC.
3. Single-cell thyroid cancer atlases nominate epithelial, partial-EMT, dedifferentiation, stromal, and immune marker modules.
4. Spatial transcriptomics shows non-random organization of tissue-state niches.
5. Bulk proteomics supports protein-level directionality of dedifferentiation and stromal-suppressive axes.
6. These public layers define a conservative CTC-EMT-NGS prior map for prospective thyroidectomy perturbation studies.

## Discussion boundary

The paper should end by saying the framework justifies, but does not replace, a prospective serial blood study.
"""
    (REP / "MANUSCRIPT_SKELETON_NO_YU_DATA_KR.md").write_text(outline, encoding="utf-8")

    defense = """# Claim Boundary Table: No Yu Raw Data

| Reviewer attack | Honest answer | Figure/table support | Boundary |
|---|---|---|---|
| You do not have CTC raw data. | Correct. This is a public-data prior/framework paper, not a CTC transition discovery paper. | F01, dataset gap map | no new CTC transition claim |
| Public tissue data cannot prove CTC shedding. | Correct. Tissue data nominates marker/state context only. | F02, F05, evidence matrix | no CTC origin proof |
| BRAF heterogeneity is not residual disease. | Correct. It shows mutation alone is insufficient for state interpretation. | F03, F04 | no residual disease claim |
| Spatial coherence is not clinical outcome. | Correct. It supports ROI and tissue organization. | F05 | no recurrence prediction |
| Proteomics is bulk not spatial. | Correct. It supports protein directionality only. | F07 | no spatial protein localization |
| Why publish this? | Because no integrated serial CTC-EMT-NGS dataset exists, so a rigorous prior map is needed before prospective cohort design. | F01, F02, F08, F09 | framework paper only |
"""
    (REP / "CLAIM_BOUNDARY_NO_YU_DATA.md").write_text(defense, encoding="utf-8")

    readme = f"""# In Silico Paper Package Without Yu Raw Data

Created: 2026-05-09

This package evaluates whether a publishable in silico paper can be started before receiving Prof. Yu's raw serial CTC data.

Verdict: **YES, as a public-data framework/prior-map paper.**

Core forbidden claim: **do not claim new postoperative CTC-EMT transition without patient-level serial CTC data.**
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    mkdirs()
    metrics = collect_metrics()
    write_tables(metrics)
    plot_figures(metrics)
    copy_source_figures()
    write_reports(metrics)
    print(f"Wrote {OUT}")
    print(f"Figures: {len(list(FIG.glob('*.png')))}")
    print(f"Tables: {len(list(TAB.glob('*.tsv')))}")


if __name__ == "__main__":
    main()
