#!/usr/bin/env python3
"""Apply accession verification findings back to master_dataset_catalog.csv.

Reads the on-disk verifier cache + master CSV, and:
  - For verified-but-non-thyroid entries: prepends [VERIFIED-WRONG: <real title>] to notes
  - For missing entries: sets download_status = VERIFIED_MISSING and prepends [VERIFIED-MISSING] to notes
  - For verified-thyroid entries: prepends [VERIFIED-OK] to notes (idempotent)

Re-runs build_zone_manifests.py at the end so per-zone manifests stay in sync.

Usage:
    python3 scripts/boundary_audit/apply_verification.py            # apply
    python3 scripts/boundary_audit/apply_verification.py --dry-run  # preview
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MASTER = REPO / "data_registry" / "manifests" / "master_dataset_catalog.csv"
CACHE = REPO / "data_registry" / "reports" / "_accession_verification_cache.json"
BUILD_SCRIPT = REPO / "scripts" / "boundary_audit" / "build_zone_manifests.py"

PREFIXES = (
    "[VERIFIED-OK]",
    "[VERIFIED-OK-XCANCER]",
    "[VERIFIED-WRONG",
    "[VERIFIED-MISSING]",
    "[VERIFIED-MANUAL-NEEDED]",
)


def strip_old(notes: str) -> str:
    """Remove any prior verification prefix so re-application is idempotent."""
    for p in PREFIXES:
        if notes.startswith(p):
            # cut up to and including the first '] '
            cut = notes.find("] ")
            if cut != -1:
                return notes[cut + 2 :]
    return notes


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    if not CACHE.exists():
        print(f"cache not found: {CACHE}. run verify_accessions.py first.", file=sys.stderr)
        sys.exit(2)

    cache = json.loads(CACHE.read_text())
    rows = list(csv.DictReader(MASTER.open()))
    fieldnames = list(rows[0].keys()) if rows else []

    n_ok = n_xcancer = n_wrong = n_missing = n_manual = n_skip = 0
    for row in rows:
        acc = (row.get("accession_or_id") or "").strip()
        source = (row.get("source_database") or "").strip()
        title = (row.get("title") or "").strip()

        # Cache key matches verifier's logic
        if source == "Publication":
            cache_key = f"PUB::{acc}::{title[:60]}"
        else:
            cache_key = acc

        if not acc or cache_key not in cache:
            n_skip += 1
            continue

        rec = cache[cache_key]
        exists = rec.get("exists")
        verified_title = (rec.get("title") or "").strip()
        thyroid = rec.get("thyroid_match")
        kind = rec.get("kind", "")

        old_notes = row.get("notes", "")
        clean_notes = strip_old(old_notes)

        # PubMed entries get gentler treatment because fuzzy-title search can miss
        # real papers when the catalog uses descriptive labels.
        if kind == "PubMed":
            if exists is True and thyroid:
                row["notes"] = f"[VERIFIED-OK] {clean_notes}".strip()
                n_ok += 1
            elif exists is True and not thyroid:
                # Cross-cancer ICI/TKI cohorts intentionally non-thyroid
                short = verified_title[:80]
                row["notes"] = f"[VERIFIED-OK-XCANCER: {short!r}] {clean_notes}".strip()
                n_xcancer += 1
            else:
                # exists False or error → manual lookup needed (NOT purgeable)
                row["notes"] = f"[VERIFIED-MANUAL-NEEDED] {clean_notes}".strip()
                n_manual += 1
            continue

        # Strict treatment for accession-based entries
        if exists is True and thyroid:
            row["notes"] = f"[VERIFIED-OK] {clean_notes}".strip()
            n_ok += 1
        elif exists is True and not thyroid:
            short = verified_title[:80]
            row["notes"] = f"[VERIFIED-WRONG: {short!r}] {clean_notes}".strip()
            row["download_status"] = "verified_wrong_subject"
            n_wrong += 1
        elif exists is False:
            row["notes"] = f"[VERIFIED-MISSING] {clean_notes}".strip()
            row["download_status"] = "verified_missing"
            n_missing += 1
        else:
            n_skip += 1

    print(f"  verified-ok:            {n_ok}")
    print(f"  verified-ok-xcancer:    {n_xcancer}")
    print(f"  verified-wrong:         {n_wrong}")
    print(f"  verified-missing:       {n_missing}")
    print(f"  verified-manual-needed: {n_manual}")
    print(f"  skipped (no cache):     {n_skip}")

    if args.dry_run:
        print("dry-run: master not modified")
        return

    with MASTER.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote updated master: {MASTER.relative_to(REPO)}")

    # Rebuild zone manifests
    subprocess.check_call(["python3", str(BUILD_SCRIPT)])


if __name__ == "__main__":
    main()
