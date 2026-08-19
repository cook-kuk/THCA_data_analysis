#!/usr/bin/env python3
"""8-gene panel vs TDS-16 vs eTDS-64 containment + parsimony analysis.

We confirmed via Paper 1 Round 8 audit that the canonical 8-gene panel is a
strict subset of TDS-16 (Landa 2016 / TCGA Cell 2014). Boucai 2023 CCR built
the enhanced TDS (eTDS, 64 genes) on top of TDS-16; the explicit eTDS list
sits in Supp Table S4 which is JS-gated on PMC. By construction
8-gene ⊂ TDS-16 ⊂ eTDS-64, so the containment argument is well-defined
without requiring the exact 64-gene roster.

This script produces:
  results/figures/figure_panel_vs_tds_etds_overlap.{png,pdf}
  results/tables/panel_vs_tds16_overlap.tsv
  results/reports/panel_vs_tds_etds_overlap.md
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_FIG = ROOT / "results" / "figures"
OUT_TAB = ROOT / "results" / "tables"
OUT_REP = ROOT / "results" / "reports"
for d in (OUT_FIG, OUT_TAB, OUT_REP):
    d.mkdir(parents=True, exist_ok=True)

# 8-gene panel (Paper 1 canonical)
PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]

# TDS-16 (TCGA-THCA Cell 2014 / Landa 2016) — sourced from
# project/metadata/tds16_genes.txt (TCGA-THCA Cell 2014 Table S5).
TDS_16 = ["DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
          "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR"]

# eTDS-64 (Boucai 2023 CCR Supp Table S4) — list is JS-gated on PMC; we assert
# the containment relation TDS_16 ⊂ eTDS_64 (declared in Boucai 2023 Methods)
# and treat the 48 extra genes as unknown until S4 is retrieved.
ETDS_EXTRA_KNOWN = []   # populate when we obtain the actual Supp S4 list

# Paper 1 Round 8 audit numbers
ROUND8 = {
    "panel_8_d": 1.783,  "panel_8_auc": 0.903,
    "tds_16_d":  1.973,  "tds_16_auc": 0.913,
    "spearman_panel_vs_tds16": 0.954,
}


# Colors
BLUE   = "#37618e"
BLUE_L = "#bccfe3"
GRAY_P = "#7e6e94"
GRAY_L = "#c2b9d0"
RED    = "#9c4742"
INK    = "#2a2a2a"
MUTED  = "#62656b"
IVORY  = "#fbf8f1"


def compute_overlap():
    panel_set = set(PANEL_8)
    tds_set = set(TDS_16)
    in_both = sorted(panel_set & tds_set)
    only_tds = sorted(tds_set - panel_set)
    only_panel = sorted(panel_set - tds_set)
    return in_both, only_tds, only_panel


def write_table(in_both, only_tds, only_panel):
    rows = []
    for g in PANEL_8:
        rows.append({"gene": g,
                     "in_panel_8": True,
                     "in_tds_16": g in TDS_16,
                     "category": "panel ∩ TDS-16" if g in TDS_16 else "panel only (unexpected)"})
    for g in only_tds:
        rows.append({"gene": g,
                     "in_panel_8": False,
                     "in_tds_16": True,
                     "category": "TDS-16 only (eTDS-extended iodide axis)"})
    import pandas as pd
    df = pd.DataFrame(rows).sort_values(["in_panel_8", "gene"], ascending=[False, True])
    df.to_csv(OUT_TAB / "panel_vs_tds16_overlap.tsv", sep="\t", index=False)
    print(f"wrote {OUT_TAB / 'panel_vs_tds16_overlap.tsv'}")


def write_report(in_both, only_tds, only_panel):
    lines = [
        "# 8-gene panel vs TDS-16 vs eTDS-64 — containment + parsimony",
        "",
        "## Source data",
        "- 8-gene panel (Paper 1 canonical): `config/eight_gene_panel.yaml`.",
        "- TDS-16 (Landa 2016, originally TCGA-THCA Cell 2014 Table S5): `project/metadata/tds16_genes.txt`.",
        "- eTDS-64 (Boucai 2023 CCR Supp Table S4, PMC10106408): the explicit gene roster is JS-gated on PMC and was not retrievable in this session. The containment claim relies on Boucai's own Methods statement that eTDS extends TDS-16.",
        "",
        "## Containment relations",
        f"- 8-gene ∩ TDS-16 = **{len(in_both)} / 8** genes — `{', '.join(in_both)}`.",
        f"- 8-gene \\ TDS-16 = **{len(only_panel)} / 8** — `{', '.join(only_panel) or 'none'}`.",
        f"- TDS-16 \\ 8-gene = **{len(only_tds)} / 16** — `{', '.join(only_tds)}`.",
        f"- TDS-16 ⊂ eTDS-64 (per Boucai 2023 Methods).",
        f"- Therefore 8-gene ⊂ TDS-16 ⊂ eTDS-64 (by transitivity).",
        "",
        "## Parsimony (Paper 1 R8 audit on TCGA-THCA)",
        f"- 8-gene Cohen d (DM1 vs DM2) = **{ROUND8['panel_8_d']:.3f}**, AUC = **{ROUND8['panel_8_auc']:.3f}**.",
        f"- TDS-16 Cohen d = **{ROUND8['tds_16_d']:.3f}**, AUC = **{ROUND8['tds_16_auc']:.3f}**.",
        f"- 8-gene captures **{ROUND8['panel_8_auc']/ROUND8['tds_16_auc']*100:.1f}%** of TDS-16's AUC discriminative power at **50% of the gene cost** (8 vs 16 genes).",
        f"- Within-sample Spearman ρ(8-gene panel z, TDS-16 score) = **{ROUND8['spearman_panel_vs_tds16']:.3f}** (Paper 1 memory `DM1 Round 8`).",
        "",
        "## Manuscript argument",
        "1. The 8-gene panel is a **parsimonious subset of the published TDS-16** that the field already uses, not an arbitrary cherry-pick.",
        "2. Boucai 2023 CCR built **eTDS-64** on top of TDS-16; our panel therefore lies inside the same iodide-handling signature space that Boucai's exceptional-responder analysis used.",
        "3. The 8 TDS-16 genes NOT in our panel (DIO2 · DUOX1 · DUOX2 · GLIS3 · SLC26A4 · SLC5A8 · THRA · THRB) form the **extended iodide axis** that can be added for sensitivity analysis without changing the headline score.",
        "4. The AUC parity (0.903 vs 0.913) supports presenting the 8-gene panel as the *operating* score with TDS-16 as the *robustness baseline* in the manuscript.",
        "",
        "## Caveats",
        f"- The explicit eTDS-64 gene roster (Boucai Supp Table S4) was NOT retrieved in this session. The PMC supplement is downloadable only via a JS-driven 'preparing to download' interstitial that defeats curl-style fetches. Either (a) manually download via browser, (b) request from corresponding author, or (c) check whether Landa 2016 / TCGA 2014 enumerates a comparable enhanced-iodide-axis set we can re-derive.",
        "- The Cohen d / AUC numbers above are from Paper 1 R8 audit on TCGA-THCA, a Tier-4 discovery cohort. Numbers will differ in Tier-1/2 anchored cohorts (Boucai, Mu, GSE151179).",
    ]
    (OUT_REP / "panel_vs_tds_etds_overlap.md").write_text("\n".join(lines))
    print(f"wrote {OUT_REP / 'panel_vs_tds_etds_overlap.md'}")


def make_figure(in_both, only_tds, only_panel):
    fig = plt.figure(figsize=(15, 7.8), facecolor=IVORY)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.18,
                          left=0.03, right=0.97, top=0.88, bottom=0.05)

    fig.suptitle(
        "8-gene panel ⊂ TDS-16 ⊂ eTDS-64 — parsimony with field-canonical containment",
        fontsize=14, fontweight="bold", color=INK, y=0.97
    )

    # ---------- Panel A: containment Venn (nested circles)
    ax = fig.add_subplot(gs[0, 0])
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.set_facecolor(IVORY)
    ax.text(50, 95, "Containment diagram", ha="center", fontsize=11.5, fontweight="bold", color=INK)

    # Outermost: eTDS-64 (Boucai)
    ax.add_patch(Circle((50, 50), 38, facecolor="#f0e8da", edgecolor="#7c6a45", lw=1.0,
                        alpha=0.85))
    ax.text(50, 86, "eTDS-64", fontsize=11.5, ha="center", color="#5b4d2e", fontweight="bold")
    ax.text(50, 82, "Boucai 2023 CCR — enhanced iodide axis", fontsize=8, ha="center",
            color=MUTED, style="italic")
    ax.text(50, 78, "(Supp S4; exact 48 extra genes JS-gated)", fontsize=7.2, ha="center",
            color=MUTED)

    # Middle: TDS-16 (Landa / TCGA)
    ax.add_patch(Circle((50, 45), 26, facecolor=GRAY_L, edgecolor=GRAY_P, lw=1.0, alpha=0.85))
    ax.text(50, 65, "TDS-16", fontsize=11.0, ha="center", color=GRAY_P, fontweight="bold")
    ax.text(50, 61, "Landa 2016 / TCGA Cell 2014", fontsize=7.8, ha="center", color=MUTED,
            style="italic")

    # Innermost: 8-gene panel (Paper 1)
    ax.add_patch(Circle((50, 40), 14, facecolor=BLUE_L, edgecolor=BLUE, lw=1.2, alpha=0.95))
    ax.text(50, 47, "8-gene panel", fontsize=10.5, ha="center", color=BLUE, fontweight="bold")
    # 8 gene names inside, compact
    inner_names = ["SLC5A5  ·  TPO  ·  TG  ·  TSHR", "PAX8  ·  NKX2-1  ·  FOXE1  ·  DIO1"]
    ax.text(50, 41, inner_names[0], fontsize=7.5, ha="center", color=INK)
    ax.text(50, 37.5, inner_names[1], fontsize=7.5, ha="center", color=INK)

    # 8 genes in TDS-16 but not in panel — listed in the gray ring
    extras_top = ", ".join(only_tds[:4])
    extras_bot = ", ".join(only_tds[4:])
    ax.text(50, 28, "TDS-16 extras: " + extras_top, fontsize=7.2, ha="center", color=MUTED)
    ax.text(50, 24.5, extras_bot, fontsize=7.2, ha="center", color=MUTED)

    ax.text(50, 10,
            "All 8 panel genes are within TDS-16 (8/8) — not a cherry-pick.",
            ha="center", fontsize=8.5, color=INK, fontweight="bold")
    ax.text(50, 5.5,
            "Boucai's eTDS-64 extends TDS-16 with 48 additional iodide-axis genes (sensitivity reserve).",
            ha="center", fontsize=7.8, color=MUTED, style="italic")

    # ---------- Panel B: parsimony bar chart
    ax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor(IVORY)
    ax.text(0.5, 1.06, "Parsimony (Paper 1 Round 8 audit on TCGA-THCA)",
            transform=ax.transAxes, ha="center", fontsize=11.5, fontweight="bold", color=INK)

    # Top row: AUC bar with 50% gene cost annotation
    panels = ["8-gene panel", "TDS-16"]
    aucs   = [ROUND8["panel_8_auc"], ROUND8["tds_16_auc"]]
    ds     = [ROUND8["panel_8_d"],   ROUND8["tds_16_d"]]
    ngenes = [8, 16]
    colors = [BLUE, GRAY_P]
    y = np.arange(len(panels))
    bars = ax.barh(y, aucs, color=colors, alpha=0.85, height=0.5, edgecolor="white")
    for yi, a, d, n, c in zip(y, aucs, ds, ngenes, colors):
        ax.text(a + 0.01, yi, f"AUC {a:.3f}   ·   Cohen d {d:.2f}   ·   {n} genes",
                va="center", fontsize=9.5, color=c, fontweight="bold")
    ax.set_yticks(y); ax.set_yticklabels(panels, fontsize=10)
    ax.set_xlim(0, 1.4)
    ax.set_xlabel("AUC (DM1 vs DM2, TCGA-THCA)", fontsize=9)
    ax.spines["right"].set_visible(False); ax.spines["top"].set_visible(False)
    ax.invert_yaxis()
    ax.axvline(0.913, color=GRAY_P, lw=0.6, ls="--", alpha=0.5)

    # Lower text block: key numbers + manuscript line
    ratio = ROUND8["panel_8_auc"] / ROUND8["tds_16_auc"] * 100
    ax.text(0.5, -0.30,
            f"8-gene captures {ratio:.1f}% of TDS-16's AUC at 50% gene cost  ·  Spearman ρ = {ROUND8['spearman_panel_vs_tds16']:.3f}",
            transform=ax.transAxes, ha="center", fontsize=10, color=INK, fontweight="bold")
    ax.text(0.5, -0.40,
            "→ Manuscript argument: parsimonious panel that lives inside the field-canonical iodide-handling signature space.",
            transform=ax.transAxes, ha="center", fontsize=8.8, color=MUTED, style="italic")
    ax.text(0.5, -0.50,
            "Caveat: Boucai eTDS-64 exact roster (Supp S4) is JS-gated on PMC; containment claim rests on Boucai's published Methods statement.",
            transform=ax.transAxes, ha="center", fontsize=7.8, color=MUTED, style="italic")

    out_png = OUT_FIG / "figure_panel_vs_tds_etds_overlap.png"
    out_pdf = OUT_FIG / "figure_panel_vs_tds_etds_overlap.pdf"
    fig.savefig(out_png, dpi=190, bbox_inches="tight", facecolor=IVORY)
    fig.savefig(out_pdf, bbox_inches="tight", facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


def main():
    in_both, only_tds, only_panel = compute_overlap()
    print(f"# 8-gene ∩ TDS-16 = {len(in_both)}/8 → {in_both}")
    print(f"# 8-gene \\ TDS-16 = {len(only_panel)} → {only_panel}")
    print(f"# TDS-16 \\ 8-gene = {len(only_tds)} → {only_tds}")
    write_table(in_both, only_tds, only_panel)
    write_report(in_both, only_tds, only_panel)
    make_figure(in_both, only_tds, only_panel)


if __name__ == "__main__":
    main()
