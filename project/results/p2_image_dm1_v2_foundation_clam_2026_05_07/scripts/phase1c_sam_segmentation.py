"""
Phase 1c (Option A++) — SAM cell-level segmentation per spot tile + DM1 score density.

Workflow per slide:
  1. For each spot tile (size 224 or 448) with DM1_like score known:
  2. SAM ViT-H "automatic mask generator" → cell-level masks
  3. Per-cell area + count + density
  4. Spot-level → cell-density vs DM1_like correlation

Output:
  sam_per_spot_cell_density.tsv  (slide, spot_id, n_cells, mean_cell_area, density_per_mm2)
  sam_dm1_density_correlation.tsv  (per-slide ρ density ↔ DM1)
  sam_overlay_per_slide_*.png  (representative overlay)

Resource:
  ~30k spots × SAM ViT-H ~0.5-1 sec/spot on A100/L40S → ~5-10 hr full
  Restrict to 1 representative slide per stage (4 slides) for time.
"""
from __future__ import annotations
import argparse, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
import torch

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tiles_root", type=Path,
                   default=Path("project/data/processed/GSE250521/tiles"))
    p.add_argument("--scores_tsv", type=Path,
                   default=Path("project/results/01_spatial_score/all_spots_scored.tsv.gz"))
    p.add_argument("--sam_weights", type=Path,
                   default=Path("/workspace/sam_weights/sam_vit_h_4b8939.pth"))
    p.add_argument("--out_dir", type=Path,
                   default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/sam"))
    p.add_argument("--tile_size", type=int, default=224)
    p.add_argument("--device", default="cuda")
    p.add_argument("--rep_per_stage", type=int, default=1,
                   help="N representative slides per stage (PT/PTC/LPTC/ATC)")
    p.add_argument("--max_spots_per_slide", type=int, default=200,
                   help="cap to bound runtime")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # SAM
    print("[load] SAM ViT-H ...")
    from segment_anything import sam_model_registry, SamAutomaticMaskGenerator
    sam = sam_model_registry["vit_h"](checkpoint=str(args.sam_weights)).to(args.device).eval()
    mask_gen = SamAutomaticMaskGenerator(
        sam,
        points_per_side=32,
        pred_iou_thresh=0.88,
        stability_score_thresh=0.92,
        min_mask_region_area=20,
    )

    # spots metadata
    scores = pd.read_csv(args.scores_tsv, sep="\t")
    scores["spot_key"] = scores["sample_id"] + "_" + scores["spot_id"]

    # representative slides per stage
    stage_slides = {}
    for stage in ["PT", "PTC", "LPTC", "ATC"]:
        slides = sorted(scores[scores["stage"] == stage]["sample_id"].unique())
        stage_slides[stage] = slides[:args.rep_per_stage]

    rows = []
    t_start = time.time()
    for stage, slides in stage_slides.items():
        for slide in slides:
            print(f"\n[{stage}] {slide}")
            slide_dir = args.tiles_root / slide
            if not slide_dir.exists():
                print(f"  [skip] {slide} dir missing")
                continue
            slide_scores = scores[scores["sample_id"] == slide].copy()
            slide_scores = slide_scores[slide_scores["in_tissue"] == 1]
            slide_scores = slide_scores.head(args.max_spots_per_slide)

            for _, srow in slide_scores.iterrows():
                tile_path = slide_dir / f"{srow['spot_id']}_size{args.tile_size}.png"
                if not tile_path.exists(): continue
                img = np.array(Image.open(tile_path).convert("RGB"))
                with torch.no_grad():
                    masks = mask_gen.generate(img)
                n_cells = len(masks)
                mean_area = float(np.mean([m["area"] for m in masks])) if masks else 0.0
                rows.append({"slide": slide, "stage": stage, "spot_id": srow["spot_id"],
                             "n_cells": n_cells, "mean_cell_area": mean_area,
                             "DM1_like_score": srow["DM1_like_score"],
                             "RAI_8_score": srow["RAI_8_score"]})

            print(f"  spots done in {time.time()-t_start:.0f}s")

    df = pd.DataFrame(rows)
    df.to_csv(args.out_dir / "sam_per_spot_cell_density.tsv", sep="\t", index=False)

    # per-slide correlation
    from scipy.stats import spearmanr
    corrs = []
    for slide, g in df.groupby("slide"):
        if len(g) < 30: continue
        for col in ["n_cells", "mean_cell_area"]:
            rho, p = spearmanr(g[col], g["DM1_like_score"])
            corrs.append({"slide": slide, "stage": g["stage"].iloc[0],
                          "metric": col, "rho_DM1": float(rho), "p": float(p),
                          "n_spots": len(g)})
    pd.DataFrame(corrs).to_csv(args.out_dir / "sam_dm1_density_correlation.tsv", sep="\t", index=False)

    # representative overlay figure
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, (stage, slides) in zip(axes, stage_slides.items()):
        if not slides:
            ax.set_visible(False); continue
        slide = slides[0]
        sub = df[df["slide"] == slide]
        if len(sub) < 5:
            ax.set_visible(False); continue
        sc = ax.scatter(sub["DM1_like_score"], sub["n_cells"],
                        c=sub["mean_cell_area"], cmap="viridis", s=18)
        plt.colorbar(sc, ax=ax, fraction=0.04)
        ax.set_xlabel("DM1_like_score"); ax.set_ylabel("SAM n_cells per tile")
        ax.set_title(f"{stage} — {slide}")
    plt.suptitle("SAM cell density vs DM1_like score (representative slides)", fontsize=12)
    plt.tight_layout()
    plt.savefig(args.out_dir / "sam_dm1_overview.png", dpi=140, bbox_inches="tight")
    plt.close()

    print(f"\n[done] SAM Phase 1c results in {args.out_dir}")


if __name__ == "__main__":
    main()
