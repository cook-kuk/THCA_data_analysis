#!/usr/bin/env python3
"""Fig 23 — HERO figure: perfect-separation strata."""
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

matplotlib.rcParams["font.family"] = ["NanumGothic", "DejaVu Sans"]

FIG = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
           "p2_image_dm1_v2_foundation_clam_2026_05_07/figures")
FIG.mkdir(parents=True, exist_ok=True)


def main():
    fig, ax = plt.subplots(figsize=(15, 7.2), dpi=220)
    ax.set_xlim(0, 15); ax.set_ylim(0, 7.2); ax.axis("off")

    title = ax.text(7.5, 6.7,
                    "Figure 23.  Perfect-separation strata — H&E only로 완벽 분리되는 환자 그룹",
                    ha="center", fontsize=16, weight="bold")

    # 4 hero cards — RAS_like, FVPTC (same 16 slides), BRAF_like, cPTC
    cards = [
        (0.4, 0.7, 3.45, 4.7, "RAS_like", "1.000", "n = 16 (3 DM1 / 13 DM2)",
         "molecular_subtype = RAS_like", "#7B1F2A", "#fff", "PERFECT"),
        (4.2, 0.7, 3.45, 4.7, "FVPTC", "1.000", "n = 16 (3 DM1 / 13 DM2)",
         "histology_subtype = FVPTC", "#7B1F2A", "#fff", "PERFECT"),
        (8.0, 0.7, 3.45, 4.7, "BRAF_like", "0.592", "n = 41 (26 DM1 / 15 DM2)",
         "molecular_subtype = BRAF_like", "#34547A", "#fff", "추가 신호 필요"),
        (11.8, 0.7, 2.85, 4.7, "cPTC", "0.592", "n = 41 (26 / 15)",
         "histology = cPTC", "#34547A", "#fff", "추가 신호 필요"),
    ]
    for (x, y, w, h, name, auc, n_str, lab, bg, fg, status) in cards:
        ax.add_patch(FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.18",
            facecolor=bg, edgecolor="black", lw=2.0, zorder=2))
        ax.text(x + w/2, y + h - 0.5, name,
                ha="center", va="center", fontsize=22, weight="bold", color=fg)
        ax.text(x + w/2, y + h - 1.15, lab,
                ha="center", va="center", fontsize=10,
                color=fg, alpha=0.85)
        ax.text(x + w/2, y + h/2 - 0.05, auc,
                ha="center", va="center", fontsize=58, weight="bold", color=fg,
                family="serif")
        ax.text(x + w/2, y + 1.3, n_str,
                ha="center", va="center", fontsize=11, color=fg, alpha=0.85)
        # status pill at bottom
        ax.add_patch(FancyBboxPatch(
            (x + 0.18, y + 0.18), w - 0.36, 0.7,
            boxstyle="round,pad=0.02,rounding_size=0.12",
            facecolor="white", edgecolor=bg, lw=1.4, zorder=3))
        ax.text(x + w/2, y + 0.53, status,
                ha="center", va="center",
                fontsize=11.5, weight="bold", color=bg)

    # bottom annotation banner
    ax.text(7.5, 0.32,
            "RAS_like (n=16, 3 DM1/13 DM2) 과 FVPTC (n=16, 3 DM1/13 DM2) 은 동일 16 슬라이드일 가능성 높음 "
            "— TCGA-THCA에서 FVPTC는 대부분 RAS-driven. 두 라벨이 같은 phenotype을 다른 축에서 잡은 것.",
            ha="center", fontsize=10.5, style="italic", color="#444")

    # left/right wing labels
    ax.text(2.05, 5.35, "H&E만으로 완벽 분리 (Δ = +0.41)",
            ha="center", fontsize=11.5, weight="bold", color="#7B1F2A")
    ax.text(11.6, 5.35, "추가 예측자 결합 권장 (분자 + image ensemble)",
            ha="center", fontsize=11, weight="bold", color="#34547A")

    fig.savefig(FIG / "fig23_hero_perfect_strata.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig23_hero_perfect_strata.pdf", bbox_inches="tight")
    plt.close(fig)
    print("WROTE", FIG / "fig23_hero_perfect_strata.png")


if __name__ == "__main__":
    main()
