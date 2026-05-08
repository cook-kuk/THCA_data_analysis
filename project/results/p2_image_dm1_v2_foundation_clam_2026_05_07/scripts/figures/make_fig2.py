#!/usr/bin/env python
"""Fig 2 — Phase 1 spatial validation grid (4x4).

Composites the 16 per-slide UNI x DM1 spatial overlay PNGs into a single
4x4 grid figure. Each panel is annotated with a short slide id, the PC1
rho_DM1 value, and the developmental stage (PT / PTC / LPTC / ATC).
Slides with |rho| > 0.3 are framed with a green border.

Inputs (LOCAL only):
  - phase1_gse250521/uni_dm1_correlation_per_slide.tsv
  - phase1_gse250521/uni_dm1_overlay_<slide>.png  (x16)

Outputs:
  - figures/fig2_phase1_spatial_grid.png
  - figures/fig2_phase1_spatial_grid.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from PIL import Image

# ---------------------------------------------------------------------------
ROOT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07"
)
PHASE1 = ROOT / "phase1_gse250521"
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

THRESH = 0.3

# Color-blind friendly: blue (positive), orange (negative), teal (pass-frame)
PASS_COLOR = "#117733"  # dark teal-green (Wong-safe; not pure red/green)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 10,
            "savefig.dpi": 300,
        }
    )

    df = pd.read_csv(PHASE1 / "uni_dm1_correlation_per_slide.tsv", sep="\t")
    pc1 = df[df["pc"] == 1].copy()

    # Stable order: by GSM id ascending => N-1..N-4, PTC-1..PTC-4, LPTC-1..LPTC-4, ATC-1..ATC-4
    pc1 = pc1.sort_values("slide").reset_index(drop=True)
    assert len(pc1) == 16, f"expected 16 slides, got {len(pc1)}"

    fig, axes = plt.subplots(4, 4, figsize=(13, 13))
    axes = axes.ravel()

    for ax, (_, row) in zip(axes, pc1.iterrows()):
        slide = str(row["slide"])
        rho = float(row["rho_DM1"])
        stage = str(row["stage"])

        # Short id: last token after the underscore (e.g. N-1, PTC-2, LPTC-4)
        short = slide.split("_", 1)[1] if "_" in slide else slide[:8]

        img_path = PHASE1 / f"uni_dm1_overlay_{slide}.png"
        if img_path.exists():
            img = Image.open(img_path)
            ax.imshow(np.asarray(img))
        else:
            ax.text(0.5, 0.5, "missing", ha="center", va="center")

        ax.set_xticks([])
        ax.set_yticks([])

        title = f"{short}  |  rho={rho:+.2f}  |  {stage}"
        ax.set_title(title, fontsize=10, pad=4)

        if abs(rho) > THRESH:
            for spine in ax.spines.values():
                spine.set_edgecolor(PASS_COLOR)
                spine.set_linewidth(3.0)
        else:
            for spine in ax.spines.values():
                spine.set_edgecolor("#999999")
                spine.set_linewidth(0.6)

    n_pass = int((pc1["rho_DM1"].abs() > THRESH).sum())
    fig.suptitle(
        f"Phase 1 spatial validation — UNI tile embeddings vs DM1 score "
        f"(GSE250521, n=16; |rho|>0.3 PASS = {n_pass}/16)",
        fontsize=12,
        y=0.995,
    )

    legend_handles = [
        mpatches.Patch(facecolor="white", edgecolor=PASS_COLOR, linewidth=3.0,
                       label=f"|rho| > {THRESH}  (PASS)"),
        mpatches.Patch(facecolor="white", edgecolor="#999999", linewidth=0.6,
                       label=f"|rho| <= {THRESH}"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.005),
        frameon=False,
        fontsize=10,
    )

    fig.tight_layout(rect=(0, 0.02, 1, 0.97))

    out_png = FIG_DIR / "fig2_phase1_spatial_grid.png"
    out_pdf = FIG_DIR / "fig2_phase1_spatial_grid.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")
    print(f"PASS = {n_pass}/16")


if __name__ == "__main__":
    main()
