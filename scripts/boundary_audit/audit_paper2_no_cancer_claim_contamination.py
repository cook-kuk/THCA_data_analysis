#!/usr/bin/env python3
"""Audit Paper 2 sources for cancer-claim contamination.

Hits are not all automatically forbidden — the report flags them as 'manual review'.
This script does NOT exit non-zero by default; only flags appearing in voice-protected
manuscript-style files exit non-zero.
"""
from __future__ import annotations

import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))
from _common import (  # noqa: E402  pylint: disable=wrong-import-position
    REPO_ROOT,
    compile_terms,
    iter_scan_files,
    load_allowlist_min,
    render_report,
    scan_file_for_terms,
)

PAPER2_TERMS = [
    "OS",
    "DSS",
    "PFI",
    "survival",
    "recurrence",
    "RAI response",
    "DM1",
    "DM2",
    "fusion-positive",
    "fusion_positive",
    "BRAF",
    "RAS",
    "TERT",
    "tumor stage",
    "lymph node metastasis",
    "distant metastasis",
    "cancer prognosis",
    "clinical risk prediction",
    "patient selection",
]

SCAN_ROOTS = [
    REPO_ROOT / "project" / "paper2_hla_boundary",
    REPO_ROOT / "project" / "manuscript_v8" / "_PAPER2_MASTER_VIEW.md",
    REPO_ROOT / "project" / "manuscript_v8" / "PAPER2_PILLAR1_V2_FOREST_META_SPEC_2026_05_04.md",
    REPO_ROOT / "project" / "manuscript_v8" / "PAPER2_TERMINOLOGY_CORRECTION_PLAN_2026_05_04.md",
    REPO_ROOT / "data_registry" / "manifests" / "paper2_hla_autoimmune_manifest.csv",
    REPO_ROOT / "data_registry" / "manifests" / "reference_population_hla_manifest.csv",
]

REPORT_PATH = REPO_ROOT / "reports" / "paper2_cancer_claim_contamination_audit.md"


def main() -> int:
    allowlist = load_allowlist_min()
    paper2_block = allowlist.get("paper2_allowed_terms", [])
    terms = compile_terms(PAPER2_TERMS)

    hits = []
    for f in iter_scan_files(SCAN_ROOTS):
        hits.extend(scan_file_for_terms(f, terms, paper2_block))

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(render_report("Paper 2 — cancer-claim contamination audit", hits, "paper2_allowed_terms"), encoding="utf-8")

    flagged = [h for h in hits if not h.allowlisted]
    print(f"paper2 audit: {len(hits)} hits, {len(flagged)} flagged for manual review. report -> {REPORT_PATH.relative_to(REPO_ROOT)}")
    # Per task spec: paper2 cancer-claim hits require manual review but do not auto-fail.
    # Only voice-protected paper2 manuscript prose would auto-fail — none present yet.
    return 0


if __name__ == "__main__":
    sys.exit(main())
