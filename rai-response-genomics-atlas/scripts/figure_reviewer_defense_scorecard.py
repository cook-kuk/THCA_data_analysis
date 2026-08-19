#!/usr/bin/env python3
"""Reviewer-defense scorecard — single composite of round-2 evidence.

Assembles four already-built sub-panels into a one-page summary:
  a) panel ⊂ TDS-16 ⊂ eTDS-64 containment (R2 cherry-pick defense)
  b) GSE151179 LOGO + random-panel null (R2)
  c) GSE112202 Tier 5 redifferentiation direction (mechanism plausibility)
  d) Figure-4-style Mu 2024 gray-zone bracket (R3 binary defense)

Pure assembly — no new computation. Reads PNG outputs and composes them.

Output: results/figures/reviewer_defense_scorecard.{png,pdf}
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "results" / "figures"

IVORY = "#fbf8f1"
INK = "#2a2a2a"
MUTED = "#62656b"
BLUE = "#37618e"

PANELS = [
    ("figure_panel_vs_tds_etds_overlap.png",
     "a · Panel ⊂ TDS-16 ⊂ eTDS-64 containment + parsimony",
     "8/8 panel genes are within TDS-16; 8-gene captures 98.9 % of TDS-16 AUC at 50 % gene cost.\n"
     "Spearman ρ = 0.954 (Paper 1 R8).  Reviewer R2 cherry-pick defense layer 1."),
    ("GSE151179_robustness.png",
     "b · GSE151179 robustness — LOGO + random-panel null + TDS-16 sensitivity",
     "LOGO d ∈ [−2.04, −1.49]; 1000-permutation empirical p = 0.014;\n"
     "TDS-16 d = −1.96, AUC = 0.96.  Reviewer R2 defense layer 2-4."),
    ("GSE112202_redifferentiation.png",
     "c · GSE112202 redifferentiation (Tier 5) — direction-of-effect",
     "6 / 8 panel genes upregulated by digoxin; median log2FC = +0.30.\n"
     "TSHR +0.95, SLC5A5 +0.73, NKX2-1 +0.52, TG +0.40.  Mechanism plausibility."),
    ("figure4_driver_grayzone.png",
     "d · Driver × zone × Mu 4-class gray zone",
     "Mu 2024 gray zone (G-RAIR + P-RAIR) = 29 / 214 = 14 %; drivers do not partition cleanly.\n"
     "BRAF-like and dark-matter zones share panel silencing but diverge on HT-overlap.  R3 defense."),
]


def main():
    fig = plt.figure(figsize=(17, 21), facecolor=IVORY)
    gs = fig.add_gridspec(4, 1, hspace=0.30, left=0.04, right=0.97, top=0.965, bottom=0.02)

    fig.suptitle("Reviewer-defense scorecard for the 8-gene RAI differentiation panel",
                 fontsize=15, fontweight="bold", color=INK, y=0.985)
    fig.text(0.5, 0.972,
             "Four reviewer-risk attack vectors with co-located evidence, one image",
             ha="center", fontsize=10.5, color=MUTED, style="italic")

    for i, (fname, title, caption) in enumerate(PANELS):
        ax = fig.add_subplot(gs[i, 0])
        path = FIG / fname
        if not path.exists():
            ax.text(0.5, 0.5, f"missing: {fname}", ha="center", color="red")
            ax.axis("off")
            continue
        img = mpimg.imread(path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(title, fontsize=12, fontweight="bold", color=INK, loc="left", pad=8)
        ax.text(0.5, -0.04, caption, transform=ax.transAxes, ha="center", va="top",
                fontsize=9.5, color=MUTED, style="italic")

    out_png = FIG / "reviewer_defense_scorecard.png"
    out_pdf = FIG / "reviewer_defense_scorecard.pdf"
    fig.savefig(out_png, dpi=140, bbox_inches="tight", facecolor=IVORY)
    fig.savefig(out_pdf, bbox_inches="tight", facecolor=IVORY)
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
