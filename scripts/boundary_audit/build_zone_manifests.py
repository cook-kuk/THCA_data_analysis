#!/usr/bin/env python3
"""Generate per-zone manifests from data_registry/manifests/master_dataset_catalog.csv.

Each manifest has the simplified schema specified in PART 4 of the data-acquisition
task. Re-run this whenever master_dataset_catalog.csv changes.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List

REPO_ROOT = Path(__file__).resolve().parents[2]
MASTER = REPO_ROOT / "data_registry" / "manifests" / "master_dataset_catalog.csv"
MANIFEST_DIR = REPO_ROOT / "data_registry" / "manifests"

ZONE_TO_FILE = {
    "PAPER1_CANCER_CORE": "paper1_cancer_core_manifest.csv",
    "PAPER2_HLA_AUTOIMMUNE": "paper2_hla_autoimmune_manifest.csv",
    "BRIDGE_QUARANTINE_HT_PTC": "bridge_quarantine_manifest.csv",
    "REFERENCE_POPULATION_HLA": "reference_population_hla_manifest.csv",
    "RESTRICTED_OR_REJECTED": "restricted_or_rejected_manifest.csv",
}

MANIFEST_HEADER = [
    "dataset_id",
    "local_path",
    "source_url",
    "assay_type",
    "sample_count",
    "phenotype_labels",
    "allowed_analysis",
    "forbidden_analysis",
    "manual_review_required",
    "notes",
]


def needs_manual_review(row: Dict[str, str]) -> str:
    flags = []
    for col in ("publication", "year", "country_or_ancestry", "sample_count_reported", "phenotypes"):
        if "NEEDS_MANUAL_REVIEW" in row.get(col, ""):
            flags.append(col)
    if row.get("download_status", "").startswith("blocked"):
        flags.append("blocked_access")
    if row.get("title", "").startswith("DUPLICATE") or row.get("notes", "").startswith("Marked DUPLICATE"):
        flags.append("duplicate")
    return "YES" if flags else "NO"


def main() -> None:
    by_zone: Dict[str, List[Dict[str, str]]] = {z: [] for z in ZONE_TO_FILE}
    with MASTER.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            zone = row["zone"]
            if zone not in by_zone:
                continue
            by_zone[zone].append(row)

    for zone, rows in by_zone.items():
        out_path = MANIFEST_DIR / ZONE_TO_FILE[zone]
        with out_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(MANIFEST_HEADER)
            for row in rows:
                writer.writerow([
                    row["dataset_id"],
                    row["local_path"],
                    row["url"],
                    row["assay_type"],
                    row["sample_count_reported"],
                    row["phenotypes"],
                    row["allowed_use"],
                    row["forbidden_use"],
                    needs_manual_review(row),
                    row["notes"],
                ])
        print(f"wrote {out_path.relative_to(REPO_ROOT)} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
