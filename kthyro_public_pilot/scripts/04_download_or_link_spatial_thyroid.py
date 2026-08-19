#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pilot_utils import append_missing_data, config, ensure_standard_dirs, maybe_symlink, pilot_root, setup_logging, stage_from_sample_id


def build_manifest(processed_dir: Path, raw_dir: Path) -> pd.DataFrame:
    rows = []
    if processed_dir.exists():
        for h5ad in sorted(processed_dir.glob("*/*.raw.h5ad")):
            sample_id = h5ad.parent.name
            image = h5ad.parent / "tissue_hires_image.png"
            link = pilot_root() / "data_raw" / "spatial" / "GSE250521" / sample_id / h5ad.name
            maybe_symlink(h5ad, link)
            if image.exists():
                maybe_symlink(image, link.parent / image.name)
            rows.append(
                {
                    "sample_id": sample_id,
                    "condition": stage_from_sample_id(sample_id),
                    "platform": "10x_Visium",
                    "has_image": image.exists(),
                    "has_spatial_coords": True,
                    "h5ad_path": str(link),
                    "image_path": str(link.parent / image.name) if image.exists() else "",
                    "raw_prefix": "",
                    "source": "local_processed_GSE250521",
                }
            )
    if not rows and raw_dir.exists():
        prefixes = sorted({p.name.split("_visium_", 1)[0] for p in raw_dir.glob("*_visium_matrix.mtx.gz")})
        for prefix in prefixes:
            rows.append(
                {
                    "sample_id": prefix,
                    "condition": stage_from_sample_id(prefix),
                    "platform": "10x_Visium_raw",
                    "has_image": (raw_dir / f"{prefix}_visium_tissue_hires_image.png.gz").exists(),
                    "has_spatial_coords": (raw_dir / f"{prefix}_visium_tissue_positions_list.csv.gz").exists(),
                    "h5ad_path": "",
                    "image_path": str(raw_dir / f"{prefix}_visium_tissue_hires_image.png.gz"),
                    "raw_prefix": str(raw_dir / f"{prefix}_visium"),
                    "source": "local_raw_GSE250521",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Link or document public thyroid spatial transcriptomics data.")
    parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("04_download_or_link_spatial_thyroid")
    cfg = config()
    local = cfg.get("local_sources", {})
    processed = Path(local.get("spatial_gse250521_processed", ""))
    raw = Path(local.get("spatial_gse250521_raw_extracted", ""))
    manifest = build_manifest(processed, raw)
    out = pilot_root() / "data_raw" / "spatial" / "spatial_sample_manifest.tsv"
    manifest.to_csv(out, sep="\t", index=False)
    if manifest.empty:
        append_missing_data(
            "GSE250521 spatial thyroid data",
            "No local processed h5ad or raw Visium files were detected. Manually download GSE250521 supplementary files from GEO, extract the Visium matrices/images, then place them under data_raw/spatial/GSE250521 or update config/config.yaml.",
        )
        logger.warning("No spatial thyroid data detected.")
    else:
        logger.info("Detected %d spatial samples.", len(manifest))
        logger.info("Wrote spatial manifest to %s", out)


if __name__ == "__main__":
    main()
