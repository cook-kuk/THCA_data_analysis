#!/usr/bin/env python3
"""
Local cell detection on GSE250521 H&E tiles — Hovernext-compatible GeoJSON.

Why local:
  - Pod DM (RunPod) hasn't been triggered yet; we still need cell calls to drive
    the SPARK Analytical pipeline scaffold.
  - GSE250521 16-slide H&E tiles are already extracted (8521 tiles total).
  - Color-deconvolution (rgb2hed) + watershed nuclei segmentation is good enough
    to feed SPARK feature snippets without a pretrained foundation model.

Pipeline per tile:
  1. RGB → HED color deconvolution → hematoxylin channel
  2. Gaussian smoothing → peak_local_max → nuclei centroids
  3. Watershed → per-cell polygon (skimage regionprops bounding ellipses)
  4. Heuristic cell-class assignment based on nuclear morphology + eosin density:
     - lymphocyte:  small (area<60), round (eccentricity<0.6), high hematoxylin
     - thyrocyte:   medium (60-160), moderate eccentricity, surrounded by eosin
     - fibroblast:  elongated (eccentricity>0.7), low hematoxylin
     - other:       fallback

Output:
  project/results/03_pathology_poc/cell_geojson/{sample_id}.geojson
  project/results/03_pathology_poc/cell_features_per_tile.tsv.gz
    columns: tile_path, sample_id, x_hires, y_hires, n_thyrocyte, n_lymphocyte,
             n_fibroblast, n_other, mean_hematoxylin, area_total
"""
from __future__ import annotations
import gzip, json, os, sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from skimage.color import rgb2hed
from skimage.feature import peak_local_max
from skimage.filters import gaussian, threshold_otsu
from skimage.segmentation import watershed
from skimage.measure import label, regionprops
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parent.parent.parent.parent
TILES = ROOT / "project/data/processed/GSE250521/tiles"
META = ROOT / "project/results/03_pathology_poc/tile_metadata_resid.tsv.gz"
OUT_GEOJSON_DIR = ROOT / "project/results/03_pathology_poc/cell_geojson"
OUT_FEAT = ROOT / "project/results/03_pathology_poc/cell_features_per_tile.tsv.gz"
OUT_GEOJSON_DIR.mkdir(parents=True, exist_ok=True)


def detect_nuclei(img_rgb):
    """Return list of (cy, cx, area, eccentricity, hema_intensity, eosin_intensity)."""
    hed = rgb2hed(img_rgb / 255.0)
    h = hed[..., 0]
    e = hed[..., 1]
    h_smooth = gaussian(h, sigma=1.5)
    try:
        thr = threshold_otsu(h_smooth)
    except Exception:
        return []
    mask = h_smooth > thr
    if mask.sum() < 30:
        return []
    distance = ndi.distance_transform_edt(mask)
    coords = peak_local_max(distance, min_distance=4, threshold_abs=2, labels=mask)
    if len(coords) == 0:
        return []
    markers = np.zeros_like(distance, dtype=np.int32)
    for i, c in enumerate(coords, start=1):
        markers[c[0], c[1]] = i
    seg = watershed(-distance, markers, mask=mask)
    out = []
    for r in regionprops(seg, intensity_image=h_smooth):
        if r.area < 8 or r.area > 600:
            continue
        cy, cx = r.centroid
        eosin_box = e[max(0,int(cy)-3):int(cy)+4, max(0,int(cx)-3):int(cx)+4]
        out.append({
            "cy": float(cy), "cx": float(cx),
            "area": float(r.area),
            "eccentricity": float(r.eccentricity),
            "hema": float(r.mean_intensity),
            "eosin": float(eosin_box.mean()) if eosin_box.size else 0.0,
        })
    return out


def classify(cell):
    a = cell["area"]; ecc = cell["eccentricity"]; h = cell["hema"]; e = cell["eosin"]
    if a < 60 and ecc < 0.6 and h > 0.10:
        return "lymphocyte"
    if ecc > 0.72 and h < 0.10:
        return "fibroblast"
    if 60 <= a <= 200 and ecc < 0.75:
        return "thyrocyte"
    return "other"


def process_tile(tile_path):
    img = np.array(Image.open(tile_path).convert("RGB"))
    cells = detect_nuclei(img)
    counts = {"thyrocyte": 0, "lymphocyte": 0, "fibroblast": 0, "other": 0}
    feats = []
    for c in cells:
        cls = classify(c)
        counts[cls] += 1
        c["class"] = cls
        feats.append(c)
    return feats, counts, img.shape


def main():
    if not TILES.exists():
        print(f"missing {TILES}"); sys.exit(1)
    meta = pd.read_csv(META, sep="\t")
    # one row per (sample_id, spot_id, tile_size); use only the 224 tiles to keep runtime sane
    m = meta[meta.tile_size == 224].copy()
    print(f"processing {len(m)} tiles ({m.sample_id.nunique()} slides)")

    rows = []
    geo_per_slide = {}
    n_done = 0
    for _, r in m.iterrows():
        sid = r.sample_id
        tp = ROOT / r.tile_path
        if not tp.exists():
            continue
        feats, counts, shape = process_tile(tp)
        rows.append({
            "tile_path": str(tp.relative_to(ROOT)),
            "sample_id": sid,
            "spot_id": r.spot_id,
            "x_hires": r.x_hires, "y_hires": r.y_hires,
            "n_thyrocyte": counts["thyrocyte"],
            "n_lymphocyte": counts["lymphocyte"],
            "n_fibroblast": counts["fibroblast"],
            "n_other": counts["other"],
            "n_total": sum(counts.values()),
            "DM1_like_score_resid": r.DM1_like_score_resid,
            "RAI_8_score_resid": r.RAI_8_score_resid,
        })
        if sid not in geo_per_slide:
            geo_per_slide[sid] = []
        for c in feats:
            geo_per_slide[sid].append({
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [c["cx"], c["cy"]]},
                "properties": {
                    "classification": {"name": c["class"]},
                    "area": c["area"], "eccentricity": c["eccentricity"],
                    "hema": c["hema"], "eosin": c["eosin"],
                    "spot_id": r.spot_id, "tile_path": str(tp.name),
                    "tile_x_hires": float(r.x_hires), "tile_y_hires": float(r.y_hires),
                },
            })
        n_done += 1
        if n_done % 200 == 0:
            print(f"  ...{n_done}/{len(m)} tiles done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_FEAT, sep="\t", index=False, compression="gzip")
    print(f"\nwrote {OUT_FEAT.relative_to(ROOT)}  ({len(df)} tiles)")

    for sid, feats in geo_per_slide.items():
        out = OUT_GEOJSON_DIR / f"{sid}.geojson"
        out.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    print(f"wrote {len(geo_per_slide)} geojson files to {OUT_GEOJSON_DIR.relative_to(ROOT)}")

    # Summary by stage
    df["stage"] = df.sample_id.str.extract(r"_(N|PTC|LPTC|ATC)-")[0]
    g = df.groupby("stage").agg(
        n_tiles=("tile_path", "count"),
        thy_per_tile=("n_thyrocyte", "mean"),
        lym_per_tile=("n_lymphocyte", "mean"),
        fib_per_tile=("n_fibroblast", "mean"),
        oth_per_tile=("n_other", "mean"),
        total_per_tile=("n_total", "mean"),
    ).round(2)
    print("\n=== per-stage cell counts (per tile mean) ===")
    print(g)


if __name__ == "__main__":
    main()
