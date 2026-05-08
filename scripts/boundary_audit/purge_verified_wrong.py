#!/usr/bin/env python3
"""Move VERIFIED-WRONG rows out of master_dataset_catalog.csv.

Reads master_dataset_catalog.csv. Rows whose `notes` field starts with
`[VERIFIED-WRONG` are appended to data_registry/manifests/_purged_verified_wrong.csv
(creating it on first run; appending on subsequent runs) and removed from master.
Re-runs build_zone_manifests.py at the end so per-zone manifests stay in sync.

Usage:
    python3 scripts/boundary_audit/purge_verified_wrong.py
    python3 scripts/boundary_audit/purge_verified_wrong.py --dry-run
    python3 scripts/boundary_audit/purge_verified_wrong.py --include-missing
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MASTER = REPO / "data_registry" / "manifests" / "master_dataset_catalog.csv"
PURGED = REPO / "data_registry" / "manifests" / "_purged_verified_wrong.csv"
BUILD_SCRIPT = REPO / "scripts" / "boundary_audit" / "build_zone_manifests.py"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    p.add_argument(
        "--include-missing",
        action="store_true",
        help="Also purge VERIFIED-MISSING rows (default: keep them, only purge WRONG)",
    )
    args = p.parse_args()

    rows = list(csv.DictReader(MASTER.open()))
    if not rows:
        print("master is empty; nothing to do", file=sys.stderr)
        return
    fieldnames = list(rows[0].keys())

    keep = []
    purge = []
    for row in rows:
        notes = (row.get("notes") or "").strip()
        if notes.startswith("[VERIFIED-WRONG"):
            purge.append(row)
        elif args.include_missing and notes.startswith("[VERIFIED-MISSING]"):
            purge.append(row)
        else:
            keep.append(row)

    print(f"  rows in master: {len(rows)}")
    print(f"  rows to purge:  {len(purge)}")
    print(f"  rows to keep:   {len(keep)}")
    if purge:
        print("  examples (first 5):")
        for row in purge[:5]:
            print(f"    - {row['dataset_id']:30s} {row['accession_or_id']:18s} {row['notes'][:80]}")

    if args.dry_run:
        print("dry-run: nothing written")
        return

    # Append to purged archive (write header on first run only)
    purged_exists = PURGED.exists() and PURGED.stat().st_size > 0
    with PURGED.open("a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, quoting=csv.QUOTE_MINIMAL)
        if not purged_exists:
            w.writeheader()
        w.writerows(purge)
    print(f"  appended {len(purge)} rows -> {PURGED.relative_to(REPO)}")

    # Rewrite master without purged rows
    with MASTER.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        w.writerows(keep)
    print(f"  master now: {len(keep)} rows -> {MASTER.relative_to(REPO)}")

    # Rebuild zone manifests
    subprocess.check_call(["python3", str(BUILD_SCRIPT)])


if __name__ == "__main__":
    main()
