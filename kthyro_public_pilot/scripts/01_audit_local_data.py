#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

from pilot_utils import classify_file, config, ensure_standard_dirs, pilot_root, setup_logging


KEYWORDS = [
    "TCGA",
    "THCA",
    "GSE250521",
    "spatial",
    "Visium",
    "h5ad",
    "h5",
    "mtx",
    "matrix.mtx",
    "features.tsv",
    "barcodes.tsv",
    "tissue_positions",
    "scalefactors",
    "filtered_feature_bc_matrix",
    "scRNA",
    "Pu",
    "DepMap",
    "PRISM",
    "Achilles",
    "CCLE",
    "expression",
    "MAF",
    "mutation",
    "clinical",
    "cnv",
    "HLA",
    "DM1",
    "dark",
]


def within_depth(path: Path, root: Path, max_depth: int) -> bool:
    try:
        rel = path.relative_to(root)
    except ValueError:
        return False
    return len(rel.parts) <= max_depth


def search_files(paths: list[Path], max_depth: int) -> list[dict]:
    rows = []
    lower_keywords = [k.lower() for k in KEYWORDS]
    for root in paths:
        if not root.exists():
            continue
        for current, dirs, files in os.walk(root):
            cur = Path(current)
            if not within_depth(cur, root, max_depth):
                dirs[:] = []
                continue
            dirs[:] = [d for d in dirs if d not in {".git", "node_modules", ".venv", "__pycache__"}]
            for name in files:
                p = cur / name
                text = str(p).lower()
                if any(k in text for k in lower_keywords):
                    try:
                        size = p.stat().st_size
                    except OSError:
                        size = None
                    rows.append(
                        {
                            "path": str(p),
                            "name": p.name,
                            "size_bytes": size,
                            "data_type": classify_file(p),
                            "source_root": str(root),
                        }
                    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit local public-data assets.")
    parser.add_argument("--max-depth", type=int, default=5)
    args = parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("01_audit_local_data")
    cfg = config()
    paths = [Path(p) for p in cfg.get("search_paths", [])]
    rows = search_files(paths, args.max_depth)
    df = pd.DataFrame(rows).drop_duplicates("path") if rows else pd.DataFrame(columns=["path", "name", "size_bytes", "data_type", "source_root"])
    out = pilot_root() / "results" / "tables" / "local_data_inventory.tsv"
    df.to_csv(out, sep="\t", index=False)
    summary = df["data_type"].value_counts().to_dict() if not df.empty else {}
    logger.info("Discovered %d matching files.", len(df))
    for key, value in summary.items():
        logger.info("%s: %d", key, value)
    logger.info("Wrote local inventory to %s", out)


if __name__ == "__main__":
    main()
