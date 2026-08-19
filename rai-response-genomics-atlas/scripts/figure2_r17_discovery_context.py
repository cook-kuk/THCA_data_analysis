#!/usr/bin/env python3
"""Figure 2 — RAI atlas discovery context: 4-pillar Paper 1 R17 cross-link.

Reuses Paper 1 R17 layer outputs (RNA bulk + sc + Korean + Mun protein + HM450 β).
Re-frames them as the RAI atlas discovery panel — the 8-gene differentiation-silencing
axis demonstrated across modalities BEFORE the Tier 1/2 label-anchored validation
covered in Figure 3.

Inputs (read-only):
  /home/seungho/personal/THCA_data_analysis/project/manuscript_biorxiv_2026_05_20/assets/
    fig_r17_two_axis_reconciliation.png
    fig_r17_sc_lu2023_two_axis.png
    fig_r17_mun2025_proteome.png
    fig_r17_zone_hm450_beta.png
    fig_r17_zone_gene_profile.png
    fig_r17_tert_zone_interaction.png

Output: results/figures/figure2_discovery_context.{png,pdf}
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

ROOT = Path(__file__).resolve().parent.parent
P1_ASSETS = Path("/home/seungho/personal/THCA_data_analysis/project/manuscript_biorxiv_2026_05_20/assets")
OUT_DIR = ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Sub-panel sources + captions
PANELS = [
    ("fig_r17_two_axis_reconciliation.png",
     "a · TCGA bulk RNA (n = 500)",
     "Two-axis decomposition · Tier 4 proxy · 14.8 % d4p2↔panel label concordance"),
    ("fig_r17_sc_lu2023_two_axis.png",
     "b · Lu 2023 single-cell (n = 14,624)",
     "ATC 38.3 % silenced + HT-overlap cells · cellular substrate"),
    ("fig_r17_mun2025_proteome.png",
     "c · Mun 2025 proteome (n = 336)",
     "ATC 58.4 % silenced + HT zone · Fisher OR = 8.54  p = 2.7×10⁻¹⁵"),
    ("fig_r17_zone_hm450_beta.png",
     "d · TCGA HM450 β (n = 518)",
     "Methylation-mediated silencing in 4/8 panel genes · DIO1·SLC5A5·TG·TPO"),
    ("fig_r17_zone_gene_profile.png",
     "e · Per-zone 8-gene profile",
     "Two silencing programs: BRAF-like → TPO/DIO1, dark-matter → TG/PAX8/TSHR"),
    ("fig_r17_tert_zone_interaction.png",
     "f · TERT × zone (n = 477)",
     "TERT⁺ uniquely enriched in dark-matter zone · OR = 2.34  p = 0.016"),
]


def main() -> None:
    fig = plt.figure(figsize=(15, 11.5), facecolor="white")
    gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.10,
                          left=0.04, right=0.97, top=0.92, bottom=0.04)

    fig.suptitle(
        "Discovery context — an 8-gene differentiation-silencing axis across 4 data pillars\n"
        "RNA bulk · single-cell · proteome · methylation  (Paper 1 R17 4-pillar reuse)",
        fontsize=13.5, fontweight="bold", color="#0e2a4a", y=0.985, family="DejaVu Sans",
    )

    positions = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
    for (r, c), (fname, title, caption) in zip(positions, PANELS):
        ax = fig.add_subplot(gs[r, c])
        img_path = P1_ASSETS / fname
        if not img_path.exists():
            ax.text(0.5, 0.5, f"missing: {fname}", ha="center", va="center", color="red")
            ax.axis("off")
            continue
        img = mpimg.imread(img_path)
        ax.imshow(img, aspect="equal")
        ax.axis("off")
        ax.set_title(title, fontsize=11, fontweight="bold", color="#0e2a4a", loc="left",
                     pad=6, family="DejaVu Sans")
        ax.text(0.5, -0.08, caption, transform=ax.transAxes, ha="center", va="top",
                fontsize=9, color="#444", style="italic", family="DejaVu Sans")

    out_png = OUT_DIR / "figure2_discovery_context.png"
    out_pdf = OUT_DIR / "figure2_discovery_context.pdf"
    fig.savefig(out_png, dpi=170, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
