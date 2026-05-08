"""Phase 3c — Build a Paper 1 supp candidate figure from Phase 1+2 results.

Output: project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase3_integration/supp_figure_image_dm1_v2.png

Layout (3-row composite):
  Row 1: Per-slide |ρ| bar plot (Phase 1)
  Row 2: Spatial overlay representative slides (PT/PTC/LPTC/ATC)
  Row 3 (if Phase 2 PASS/MARGINAL): CLAM 5-fold AUC + closure baseline comparison
"""
from __future__ import annotations
import argparse, os
from pathlib import Path
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

for _p in ["/home/seungho/.local/share/fonts/NanumGothic-Regular.ttf",
           "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"]:
    if os.path.exists(_p):
        try: fm.fontManager.addfont(_p)
        except Exception: pass
plt.rcParams.update({"font.size":10, "font.family":["NanumGothic","Noto Sans CJK JP","DejaVu Sans"],
                     "axes.unicode_minus": False})

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path,
                   default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"))
    args = p.parse_args()

    p1 = args.root/"phase1_gse250521"; p2 = args.root/"phase2_tcga_clam"
    out = args.root/"phase3_integration"; out.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(13, 13), gridspec_kw={"height_ratios":[1.0, 1.5, 1.2]})

    # Row 1 — per-slide |ρ| bar
    p1_corr = p1/"uni_dm1_correlation_per_slide.tsv"
    if p1_corr.exists():
        df = pd.read_csv(p1_corr, sep="\t")
        pc1 = df[df["pc"]==1].sort_values("rho_DM1")
        ax = axes[0]
        colors = ["#c0392b" if abs(r)>0.3 else "#7f8fa6" for r in pc1["rho_DM1"]]
        ax.barh(range(len(pc1)), pc1["rho_DM1"], color=colors, edgecolor="black")
        ax.set_yticks(range(len(pc1))); ax.set_yticklabels(pc1["slide"], fontsize=8)
        ax.axvline(0, color="black", lw=0.5)
        ax.axvline(0.3, color="grey", linestyle=":"); ax.axvline(-0.3, color="grey", linestyle=":")
        ax.set_xlabel("Spearman ρ (UNI PC1 ↔ DM1_like_score)")
        ax.set_title("Row 1 — Per-slide UNI ↔ DM1 correlation (kill-switch threshold |ρ|=0.3)", fontsize=11)
    else:
        axes[0].text(0.5, 0.5, "Phase 1 correlation TSV not found", ha="center", va="center")
        axes[0].axis("off")

    # Row 2 — placeholder for spatial overlay (will inset existing PNG when found)
    ax = axes[1]
    overlays = sorted(p1.glob("uni_dm1_overlay_*.png"))
    if overlays:
        # show one representative
        from matplotlib.image import imread
        img = imread(str(overlays[0]))
        ax.imshow(img); ax.axis("off")
        ax.set_title(f"Row 2 — Representative spatial overlay  ({overlays[0].name})", fontsize=11)
    else:
        ax.text(0.5, 0.5, "No spatial overlay PNG yet", ha="center", va="center")
        ax.axis("off")

    # Row 3 — Phase 2 CLAM AUC bar + closure baseline
    ax = axes[2]
    p2_report = p2/"PHASE2_REPORT.md"
    if p2_report.exists():
        import re
        text = p2_report.read_text()
        m = re.search(r"Pooled cross-fold AUC:\s*([\d.]+)", text)
        auc = float(m.group(1)) if m else None
        bars = {"closure ResNet50 baseline": 0.55, "v2 UNI+CLAM": auc} if auc else {"closure ResNet50 baseline":0.55}
        names = list(bars.keys()); vals = list(bars.values())
        colors = ["#7f8fa6", "#c0392b"]
        ax.bar(names, vals, color=colors[:len(names)], edgecolor="black")
        ax.axhline(0.70, color="black", linestyle="--", label="kill-switch threshold AUC=0.70")
        ax.set_ylim(0.4, 1.0); ax.set_ylabel("AUC (DM1 vs DM2)")
        for i, v in enumerate(vals):
            ax.text(i, v+0.01, f"{v:.3f}", ha="center", fontsize=10, fontweight="bold")
        ax.legend(); ax.set_title("Row 3 — Phase 2 CLAM AUC vs closure baseline", fontsize=11)
    else:
        ax.text(0.5, 0.5, "Phase 2 not run / report missing", ha="center", va="center")
        ax.axis("off")

    plt.suptitle("Paper 1 supp candidate — Image-DM1 v2 (foundation-model + CLAM)", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out/"supp_figure_image_dm1_v2.png", dpi=160, bbox_inches="tight")
    plt.close()
    print(f"[done] {out/'supp_figure_image_dm1_v2.png'}")

if __name__ == "__main__":
    main()
