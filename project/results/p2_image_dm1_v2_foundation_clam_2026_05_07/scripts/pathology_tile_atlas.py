"""
Build a visual atlas of GSE250521 H&E tiles at the extremes of spatial DM1.

Panels:
  1. Highest transcriptomic DM1_like_score tiles
  2. Lowest transcriptomic DM1_like_score tiles
  3. Highest leave-slide-out H&E-predicted DM1 morphology tiles
  4. Lowest leave-slide-out H&E-predicted DM1 morphology tiles
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


REPO = Path(os.environ.get("THCA_REPO_ROOT", Path.cwd()))
P2_ROOT = Path(
    os.environ.get(
        "THCA_P2_ROOT",
        REPO / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07",
    )
)
TILES = REPO / "project/data/processed/GSE250521/tiles"
SCORES = REPO / "project/results/01_spatial_score/all_spots_scored.tsv.gz"
PRED = P2_ROOT / "analysis_supp/pathology_spatial_multiaxis_2026_05_09/spatial_multiaxis_spot_predictions.tsv"
OUT_DIR = P2_ROOT / "analysis_supp/pathology_tile_atlas_2026_05_09"


def load_rows() -> pd.DataFrame:
    scores = pd.read_csv(SCORES, sep="\t")
    pred = pd.read_csv(PRED, sep="\t")
    pred = pred[pred["target"] == "DM1_like_score"].copy()
    scores["slide"] = scores["sample_id"]
    df = scores.merge(
        pred[["slide", "spot_id", "pred_z"]],
        on=["slide", "spot_id"],
        how="inner",
    )
    df = df[df["in_tissue"] == 1].copy()
    df["tile_path"] = df.apply(
        lambda r: TILES / str(r["slide"]) / f"{r['spot_id']}_size224.png", axis=1
    )
    df = df[df["tile_path"].map(lambda p: Path(p).exists())].copy()
    return df


def pick_extremes(df: pd.DataFrame, value: str, n: int = 16):
    high = df.sort_values(value, ascending=False).head(n).copy()
    low = df.sort_values(value, ascending=True).head(n).copy()
    return high, low


def draw_panel(rows: pd.DataFrame, title: str, value_col: str, tile_size: int = 112, cols: int = 4) -> Image.Image:
    rows = rows.reset_index(drop=True)
    label_h = 34
    title_h = 34
    gap = 6
    panel_w = cols * tile_size + (cols + 1) * gap
    rows_n = int(np.ceil(len(rows) / cols))
    panel_h = title_h + rows_n * (tile_size + label_h) + (rows_n + 1) * gap
    canvas = Image.new("RGB", (panel_w, panel_h), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    draw.rectangle([0, 0, panel_w, title_h], fill=(22, 35, 49))
    draw.text((8, 10), title, fill=(255, 210, 138), font=font)
    for i, row in rows.iterrows():
        r = i // cols
        c = i % cols
        x = gap + c * (tile_size + gap)
        y = title_h + gap + r * (tile_size + label_h + gap)
        img = Image.open(row["tile_path"]).convert("RGB").resize((tile_size, tile_size))
        canvas.paste(img, (x, y))
        stage = str(row["stage"])
        val = float(row[value_col])
        text = f"{stage} {val:.2f}"
        draw.rectangle([x, y + tile_size, x + tile_size, y + tile_size + label_h], fill=(245, 245, 245))
        draw.text((x + 3, y + tile_size + 4), text, fill=(20, 20, 20), font=font)
        draw.text((x + 3, y + tile_size + 17), str(row["spot_id"])[:14], fill=(70, 70, 70), font=font)
    return canvas


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_rows()
    actual_high, actual_low = pick_extremes(df, "DM1_like_score", 16)
    pred_high, pred_low = pick_extremes(df, "pred_z", 16)
    actual_high.to_csv(OUT_DIR / "actual_dm1_high_tiles.tsv", sep="\t", index=False, na_rep="NA")
    actual_low.to_csv(OUT_DIR / "actual_dm1_low_tiles.tsv", sep="\t", index=False, na_rep="NA")
    pred_high.to_csv(OUT_DIR / "predicted_dm1_high_tiles.tsv", sep="\t", index=False, na_rep="NA")
    pred_low.to_csv(OUT_DIR / "predicted_dm1_low_tiles.tsv", sep="\t", index=False, na_rep="NA")

    panels = [
        draw_panel(actual_high, "Actual spatial DM1_like HIGH", "DM1_like_score"),
        draw_panel(actual_low, "Actual spatial DM1_like LOW", "DM1_like_score"),
        draw_panel(pred_high, "H&E-predicted DM1 morphology HIGH", "pred_z"),
        draw_panel(pred_low, "H&E-predicted DM1 morphology LOW", "pred_z"),
    ]
    gap = 18
    w = max(p.width for p in panels) * 2 + gap
    h = max(p.height for p in panels[:2]) + max(p.height for p in panels[2:]) + gap
    canvas = Image.new("RGB", (w, h), (235, 238, 240))
    canvas.paste(panels[0], (0, 0))
    canvas.paste(panels[1], (panels[0].width + gap, 0))
    canvas.paste(panels[2], (0, panels[0].height + gap))
    canvas.paste(panels[3], (panels[2].width + gap, panels[1].height + gap))
    out = OUT_DIR / "fig_spatial_dm1_tile_atlas.png"
    canvas.save(out)

    # Also make a compact PDF wrapper for easy insertion into dossiers.
    fig = plt.figure(figsize=(12, 12))
    plt.imshow(canvas)
    plt.axis("off")
    plt.title("GSE250521 H&E tile atlas: spatial DM1 transcriptomic and morphology-predicted extremes")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_spatial_dm1_tile_atlas.pdf")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
