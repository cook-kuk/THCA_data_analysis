#!/usr/bin/env python3
"""
CCI spatial overlay rendering — per slide, 4-panel:
  H&E | CAF score | RAI_lineage score | CAF×RAI lag-product (avoidance map)

The 4th panel = CAF[i] × RAI_lag[i] (negative product → strong spatial avoidance).
This visualizes the headline finding (CAF × RAI = -0.27 spatial avoidance).

Output:
  project/results/03_pathology_poc/cci_overlays/{sid}.png
  project/results/03_pathology_poc/cci_overlays/index.html
"""
from __future__ import annotations
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent.parent.parent
GS = ROOT / "project/data/processed/GSE250521"
OUT = ROOT / "project/results/03_pathology_poc/cci_overlays"
OUT.mkdir(parents=True, exist_ok=True)

CAF = ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"]
RAI = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]


def get_scalefactor(spatial_uns):
    if not spatial_uns: return 1.0
    libs = list(spatial_uns.keys())
    return float(spatial_uns[libs[0]].get("scalefactors", {}).get("tissue_hires_scalef", 1.0))


def hires_path(sdir):
    cands = list(sdir.rglob("tissue_hires_image.png"))
    return cands[0] if cands else None


def module(X, sym, genes):
    cols = [sym[g] for g in genes if g in sym]
    if len(cols) < 2: return None
    sub = X[:, cols]
    return ((sub - sub.mean(0)) / (sub.std(0) + 1e-6)).mean(1)


def render(sid):
    sdir = GS / sid
    a = ad.read_h5ad(sdir / f"{sid}.scored.h5ad")
    sf = get_scalefactor(a.uns.get("spatial", {}))
    bg = hires_path(sdir)
    if bg is None: return None
    bg_img = np.array(Image.open(bg).convert("RGB"))
    coords = a.obsm["spatial"] * sf
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    sym = {s: i for i, s in enumerate(a.var.index.values)}
    caf = module(X, sym, CAF)
    rai = module(X, sym, RAI)
    if caf is None or rai is None: return None

    # CAF × RAI lag product on raw coords
    tree = cKDTree(a.obsm["spatial"])
    _, idx = tree.query(a.obsm["spatial"], k=7)
    rai_lag = rai[idx[:, 1:]].mean(axis=1)
    avoid_map = caf * (-rai_lag)  # high = CAF-rich + RAI-poor neighbor (avoidance hot zone)

    fig, axes = plt.subplots(1, 4, figsize=(20, 5.2), constrained_layout=True)
    panels = [
        ("H&E", None, None),
        ("CAF score (FAP/ACTA2/COL1A1…)", caf, "viridis"),
        ("RAI_lineage score (TPO/DIO1/PAX8…)", rai, "RdBu_r"),
        ("CAF × −RAI_lag (avoidance hot)", avoid_map, "magma"),
    ]
    for ax, (title, vals, cmap) in zip(axes, panels):
        ax.imshow(bg_img); ax.set_title(f"{sid} · {title}", fontsize=10); ax.axis("off")
        if vals is not None and np.std(vals) > 0:
            lo, hi = np.nanpercentile(vals, [5, 95])
            sc = ax.scatter(coords[:, 0], coords[:, 1], c=vals, cmap=cmap, s=4, alpha=0.75,
                            vmin=lo, vmax=hi, edgecolor="none")
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    out = OUT / f"{sid}.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    sids = sorted([d.name for d in GS.iterdir() if (d / f"{d.name}.scored.h5ad").exists()])
    out_paths = []
    for sid in sids:
        p = render(sid)
        print(f"  {sid}: {'ok' if p else 'skip'}")
        if p: out_paths.append(p)

    html = ['<!doctype html><html><head><meta charset="utf-8"><title>CCI · CAF–RAI avoidance maps</title>',
            '<style>body{background:#0a0d18;color:#e1e8f0;font-family:Inter,sans-serif;margin:24px}',
            'h1{margin:0 0 8px} .sub{color:#94a3b8;margin-bottom:24px}',
            'figure{margin:0 0 18px;border:1px solid #1f2a44;border-radius:10px;padding:8px;background:#0e1426}',
            'figcaption{font-size:12px;color:#94a3b8;padding:6px 4px}',
            'img{width:100%;display:block;border-radius:6px}</style></head><body>',
            '<h1>CCI spatial overlays — CAF × −RAI_lag avoidance maps (16 slides)</h1>',
            '<div class="sub">Right-most panel = CAF score × neighbor-mean(−RAI). Bright magma pixels = CAF-rich regions whose neighbors are RAI-poor (the −0.27 spatial avoidance signature).</div>']
    for p in out_paths:
        html.append(f'<figure><img src="{p.name}" /><figcaption>{p.stem}</figcaption></figure>')
    html.append("</body></html>")
    (OUT / "index.html").write_text("\n".join(html))
    print(f"\nwrote gallery {OUT / 'index.html'}  ({len(out_paths)} panels)")


if __name__ == "__main__":
    main()
