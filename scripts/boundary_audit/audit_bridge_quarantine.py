#!/usr/bin/env python3
"""Audit bridge-zone outputs.

Enforces:
  - bridge files must NOT export HLA-allele association tables
  - bridge files must NOT export cancer-survival × HLA association tables
  - HLA gene-expression context is allowed only as transcriptomic context
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

# Forbidden patterns in bridge zone:
#   - allele-typed HLA references (e.g., "HLA-DRB1*", "DPB1*05:01")
#   - cancer outcome × HLA join language ("survival", "OS", "DSS", "recurrence")
#     when co-occurring with HLA terms in the same file (we flag the outcome term;
#     reviewer reads context to decide).
BRIDGE_FORBIDDEN_TERMS = [
    "HLA-DRB1*",
    "HLA-DPB1*",
    "HLA-DQB1*",
    "HLA-DQA1*",
    "HLA-DPA1*",
    "DPB1*0",
    "DQB1*0",
    "DRB1*0",
    "allele association",
    "carrier frequency",
    "OS",
    "DSS",
    "PFI",
    "survival",
    "recurrence",
    "patient selection",
    "clinical risk prediction",
]

SCAN_ROOTS = [
    REPO_ROOT / "data_registry" / "bridge_quarantine_ht_ptc",
    REPO_ROOT / "data_registry" / "manifests" / "bridge_quarantine_manifest.csv",
    REPO_ROOT / "project" / "data" / "external" / "gse286332",  # cross-zone bridge
]

REPORT_PATH = REPO_ROOT / "reports" / "bridge_quarantine_audit.md"


def main() -> int:
    allowlist = load_allowlist_min()
    bridge_block = allowlist.get("bridge_allowed_terms", [])
    terms = compile_terms(BRIDGE_FORBIDDEN_TERMS)

    hits = []
    for f in iter_scan_files(SCAN_ROOTS):
        hits.extend(scan_file_for_terms(f, terms, bridge_block))

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(render_report("Bridge quarantine audit", hits, "bridge_allowed_terms"), encoding="utf-8")

    flagged = [h for h in hits if not h.allowlisted]
    print(f"bridge audit: {len(hits)} hits, {len(flagged)} flagged. report -> {REPORT_PATH.relative_to(REPO_ROOT)}")
    return 1 if flagged else 0


if __name__ == "__main__":
    sys.exit(main())
