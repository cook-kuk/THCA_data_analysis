#!/usr/bin/env python3
"""Audit Paper 1 sources for HLA / Hashimoto / Graves / AITD contamination.

Exits non-zero if any unallowlisted hit is found.
Writes a markdown report to reports/paper1_hla_contamination_audit.md.
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

PAPER1_TERMS = [
    "HLA",
    "MHC",
    "DPB1",
    "DQB1",
    "DRB1",
    "DQA1",
    "DPA1",
    "AFND",
    "carrier frequency",
    "allele frequency",
    "Hashimoto",
    "Graves",
    "HT overlap",
    "GD",
    "autoimmune susceptibility",
    "AITD",
]

SCAN_ROOTS = [
    REPO_ROOT / "project" / "manuscript_v8",
    REPO_ROOT / "data_registry" / "manifests" / "paper1_cancer_core_manifest.csv",
    REPO_ROOT / "reports",
    REPO_ROOT / "submission",
]

EXTRA_SCAN_FILES = [
    REPO_ROOT / "project" / "manuscript_v8" / "10_full_manuscript_compiled.md",
    REPO_ROOT / "project" / "manuscript_v8" / "05_figure_captions.md",
    REPO_ROOT / "project" / "manuscript_v8" / "08_cover_letter.md",
    REPO_ROOT / "project" / "manuscript_v8" / "09_reviewer_qa.md",
    REPO_ROOT / "project" / "manuscript_v8" / "paper1.html",
]

REPORT_PATH = REPO_ROOT / "reports" / "paper1_hla_contamination_audit.md"


def main() -> int:
    allowlist = load_allowlist_min()
    paper1_block = allowlist.get("paper1_allowed_terms", [])
    terms = compile_terms(PAPER1_TERMS)

    roots = list(SCAN_ROOTS)
    roots.extend(p for p in EXTRA_SCAN_FILES if p.exists())

    hits = []
    for f in iter_scan_files(roots):
        hits.extend(scan_file_for_terms(f, terms, paper1_block))

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(render_report("Paper 1 — HLA contamination audit", hits, "paper1_allowed_terms"), encoding="utf-8")

    flagged = [h for h in hits if not h.allowlisted]
    print(f"paper1 audit: {len(hits)} hits, {len(flagged)} flagged. report -> {REPORT_PATH.relative_to(REPO_ROOT)}")
    return 1 if flagged else 0


if __name__ == "__main__":
    sys.exit(main())
