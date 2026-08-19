#!/usr/bin/env python3
"""01 — print priority datasets from config/datasets.yaml with access state.

Usage: python3 scripts/01_search_geo_metadata.py
"""
from __future__ import annotations
from pathlib import Path
import sys
import yaml

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "datasets.yaml"


def main() -> None:
    if not CONFIG.exists():
        print(f"missing {CONFIG}", file=sys.stderr); sys.exit(1)
    cfg = yaml.safe_load(CONFIG.read_text())
    rows = []
    for key, ds in cfg.get("datasets", {}).items():
        if not isinstance(ds, dict): continue
        rows.append({
            "key": key,
            "accession": ds.get("accession", "?"),
            "cancer": ds.get("cancer_type", "?"),
            "omics": ds.get("omics", "?"),
            "tier": ds.get("label_tier", "?"),
            "access": ds.get("access", "?"),
            "priority": ds.get("priority", "?"),
            "next_action": ds.get("next_action", "?"),
        })
    fmt = "{key:30s} {accession:18s} {tier:32s} priority={priority} {access:24s} → {next_action}"
    print(f"# {len(rows)} datasets in inventory")
    for r in sorted(rows, key=lambda x: (str(x["priority"]), x["key"])):
        print(fmt.format(**{k: str(v) for k, v in r.items()}))


if __name__ == "__main__":
    main()
