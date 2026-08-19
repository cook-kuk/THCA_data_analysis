#!/usr/bin/env python3
"""07 — emit results/tables/dataset_inventory.csv from config/datasets.yaml.

Usage: python3 scripts/07_build_dataset_inventory.py
"""
from __future__ import annotations
from pathlib import Path
import yaml
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config" / "datasets.yaml"
OUT = ROOT / "results" / "tables" / "dataset_inventory.csv"


def main() -> None:
    cfg = yaml.safe_load(CONFIG.read_text())
    rows = []
    for key, ds in cfg.get("datasets", {}).items():
        if not isinstance(ds, dict):
            continue
        rows.append({
            "dataset_name": key,
            "accession": ds.get("accession", ""),
            "cancer_type": ds.get("cancer_type", ""),
            "sample_size": ds.get("sample_count") or ds.get("sample_count_hint", ""),
            "sample_type": ds.get("sample_type", ""),
            "omics_type": ds.get("omics", ""),
            "platform": ds.get("platform_hint", ""),
            "rai_label_type": ", ".join(str(x) for x in (ds.get("expected_labels", []) or [])),
            "label_tier": ds.get("label_tier", ""),
            "has_expression": "expression" in str(ds.get("omics", "")).lower() or "rna" in str(ds.get("omics", "")).lower(),
            "has_mutation": "ngs" in str(ds.get("omics", "")).lower() or "wes" in str(ds.get("omics", "")).lower() or bool(ds.get("driver_panel")),
            "has_clinical_response": "Tier_1" in str(ds.get("label_tier", "")) or "Tier_2" in str(ds.get("label_tier", "")),
            "public_download": ds.get("access", "") == "public",
            "controlled_access": "controlled" in str(ds.get("access", "")) or "request" in str(ds.get("access", "")),
            "useful_for": ds.get("paper_hint") or ds.get("paper") or "",
            "limitations": "",
            "priority": ds.get("priority", ""),
            "next_action": ds.get("next_action", ""),
        })
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"# wrote {OUT}  ({df.shape[0]} datasets)")


if __name__ == "__main__":
    main()
