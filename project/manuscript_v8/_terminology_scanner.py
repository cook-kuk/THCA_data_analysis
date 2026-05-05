#!/usr/bin/env python3
"""
Paper 2 terminology scanner — dry-run only, NO file edits.

Per PAPER2_TERMINOLOGY_CORRECTION_PLAN_2026_05_04.md (§4 Dry-run scan plan).

Output: _terminology_scan_hits.tsv with columns:
  file_path | line_no | matching_text | pattern | category | territory | suggested_action

Usage:
  python3 _terminology_scanner.py
"""

from __future__ import annotations
from pathlib import Path
import re
import sys

# -----------------------------------------------------------------------------
# Patterns (case-insensitive, word-boundary where applicable)
# -----------------------------------------------------------------------------

# §1 Replacement candidates: (pattern, suggested_replacement, category, territory, action)
REPLACE_CANDIDATES = [
    (r"\bautoimmune-PTC\b",                "Hashimoto-overlap PTC",            "Replace_§1", "Paper2_terminology", "REPLACE"),
    (r"\bautoimmune PTC\b",                "Hashimoto-overlap PTC",            "Replace_§1", "Paper2_terminology", "REPLACE"),
    (r"\bautoimmune-thyroid continuum\b",  "Hashimoto-overlap context",        "Replace_§1", "Paper2_terminology", "REPLACE_OR_REMOVE"),
    (r"\bautoimmune thyroid continuum\b",  "Hashimoto-overlap context",        "Replace_§1", "Paper2_terminology", "REPLACE_OR_REMOVE"),
    (r"\bPan-Asian autoimmune-thyroid\b",  "REMOVE or Paper 4 reserve",        "Replace_§1", "Paper2_terminology", "REMOVE_OR_MOVE"),
    (r"\bPan-Asian autoimmune\b",          "manual review",                    "Replace_§1", "Paper2_terminology", "MANUAL_REVIEW"),
]

# §2.1 Forbidden — Paper 4 (GD) territory
FORBIDDEN_GD = [
    r"\bGraves\b",
    r"\bGD\b",
    r"\bTSAb\b",
    r"\bTSI\b",
    r"\bthyroid-stimulating antibody\b",
    r"\bhyperthyroidism\b",
    r"\bthyrotoxicosis\b",
    r"\bexophthalmos\b",
    r"Graves['']?\s+ophthalmopathy",
    r"\bTED\b",
    r"\bthyroid eye disease\b",
    r"\bBundang Graves'?\s*cohort\b",
]

# §2.2 Forbidden — Paper 1 (cancer-driver-deep) territory
FORBIDDEN_PAPER1 = [
    r"\bTROP2\b",
    r"\bsacituzumab\b",
    r"\bH&E-DM1\b",
    r"\bimage-DM1\b",
    r"\bTCGA WSI\b",
    r"\bDM1 cluster mechanism deep\b",
    r"\bBRAF/RAS driver mechanism deep\b",
    r"\bTERT promoter kinetics\b",
    r"\bkinase fusion mechanism\b",
    r"\bdark matter subtype mechanism\b",
]

# §2.3 Forbidden — Paper 3 (ICI/immune) territory
FORBIDDEN_PAPER3 = [
    r"\bICI response prediction\b",
    r"\bHLA LOH\b",
    r"\bneoantigen\b",
    r"\bDIAL audit\b",
    r"\bpan-cancer ICI\b",
    r"\bC5AR1 anti-PD-1 synergy\b",
]

FORBIDDEN_GROUPS = [
    ("Forbidden_§2.1", "Paper4_GD",      "REVIEW_OR_REMOVE_OR_MOVE_PAPER4", FORBIDDEN_GD),
    ("Forbidden_§2.2", "Paper1_cancer",  "REVIEW_OR_KEEP_IN_PAPER1_FILE",   FORBIDDEN_PAPER1),
    ("Forbidden_§2.3", "Paper3_ICI",     "REVIEW_OR_MOVE_PAPER3",           FORBIDDEN_PAPER3),
]

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent  # manuscript_v8/
MEMORY_ROOT = Path.home() / ".claude" / "projects" / "-home-seungho-personal-THCA-data-analysis" / "memory"

# Whitelist (review-only, scan but tag as WHITELIST)
WHITELIST_PATHS = {
    "_audit_HT_vs_GD_2026_05_02.md",                      # audit doc itself — terms appear deliberately
    "_PAPER2_MASTER_VIEW.md",                             # meta-doc
    "_session_origin_paper2_isolated.md",                 # meta-doc
    "_paper2_isolated_resume_response.md",                # meta-doc
    "PAPER2_PILLAR1_V2_FOREST_META_SPEC_2026_05_04.md",   # spec — terms appear in interpretation guard
    "PAPER2_TERMINOLOGY_CORRECTION_PLAN_2026_05_04.md",   # plan itself
    "_terminology_scanner.py",                             # this script
    "_terminology_scan_hits.tsv",                          # output (will not exist on first run)
    "v18_paper2_HT_isolated.md",                          # already clean memory
    "v17_landa2016_cite_save.md",                         # Paper 1 cite, terms by design
    "v19_paper3_ici_track_a.md",                          # Paper 3 territory
    "v19_paper4_GD_backlog.md",                           # Paper 4 territory by design
    "v19_paper3_GD_gating.md",                            # Paper 4 (legacy) territory
    "paper_numbering_2026_05_04.md",                      # canonical numbering, terms by design
    "v17_graves_pivot.md",                                # Paper 4 backburner memory
    "v17_marathon_mode_post_pillar1.md",                  # marathon agreement
    "v17_sprint_vs_marathon_violation.md",                # marathon agreement
}

# Exclusions (do not scan)
EXCLUDE_PATHS = {
    "_archive_2026_05_03_sprint_draft.md",                # frozen archive
    "04_intro_1_1_hook.md",                               # empty W1 voice-first
}

def collect_targets() -> list[Path]:
    targets: list[Path] = []
    # manuscript_v8 markdown files
    for md in sorted(PROJECT_ROOT.glob("*.md")):
        if md.name in EXCLUDE_PATHS:
            continue
        targets.append(md)
    # memory markdown files (Paper 2 relevant + global index)
    if MEMORY_ROOT.exists():
        for md in sorted(MEMORY_ROOT.glob("*.md")):
            if md.name in EXCLUDE_PATHS:
                continue
            targets.append(md)
    return targets


# -----------------------------------------------------------------------------
# Scanning
# -----------------------------------------------------------------------------
def scan_file(path: Path) -> list[dict]:
    hits: list[dict] = []
    is_whitelist = path.name in WHITELIST_PATHS
    # Display path: short form depending on root
    p = str(path)
    if "manuscript_v8" in p:
        display_path = "manuscript_v8/" + path.name
    elif ".claude" in p and "memory" in p:
        display_path = "memory/" + path.name
    else:
        display_path = str(path)
    try:
        text = path.read_text(encoding="utf-8")
    except Exception as e:
        return []
    for line_no, line in enumerate(text.splitlines(), start=1):
        # Replace candidates
        for pattern, replacement, category, territory, action in REPLACE_CANDIDATES:
            for m in re.finditer(pattern, line, flags=re.IGNORECASE):
                act = "WHITELIST_REVIEW" if is_whitelist else action
                hits.append({
                    "file_path": display_path,
                    "line_no": line_no,
                    "matching_text": m.group(0),
                    "pattern": pattern,
                    "category": category,
                    "territory": territory,
                    "suggested_action": act,
                    "replacement_or_note": replacement,
                    "line_excerpt": line.strip()[:120],
                })
        # Forbidden groups
        for category, territory, action, patterns in FORBIDDEN_GROUPS:
            for pattern in patterns:
                for m in re.finditer(pattern, line, flags=re.IGNORECASE):
                    act = "WHITELIST_REVIEW" if is_whitelist else action
                    hits.append({
                        "file_path": display_path,
                        "line_no": line_no,
                        "matching_text": m.group(0),
                        "pattern": pattern,
                        "category": category,
                        "territory": territory,
                        "suggested_action": act,
                        "replacement_or_note": "",
                        "line_excerpt": line.strip()[:120],
                    })
    return hits


def main() -> int:
    targets = collect_targets()
    all_hits: list[dict] = []
    for path in targets:
        all_hits.extend(scan_file(path))

    # Output TSV
    out = PROJECT_ROOT / "_terminology_scan_hits.tsv"
    cols = ["file_path", "line_no", "matching_text", "pattern", "category",
            "territory", "suggested_action", "replacement_or_note", "line_excerpt"]
    with out.open("w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for h in all_hits:
            f.write("\t".join(str(h[c]).replace("\t", " ").replace("\n", " ") for c in cols) + "\n")

    # Summary by category
    from collections import Counter
    cat_counts = Counter(h["category"] for h in all_hits)
    territory_counts = Counter(h["territory"] for h in all_hits)
    action_counts = Counter(h["suggested_action"] for h in all_hits)
    file_counts = Counter(h["file_path"] for h in all_hits)

    print(f"Scan complete. {len(all_hits)} total hits across {len(targets)} files.")
    print(f"Output: manuscript_v8/_terminology_scan_hits.tsv")
    print()
    print("Hits by category:")
    for k, v in sorted(cat_counts.items()):
        print(f"  {k}: {v}")
    print()
    print("Hits by territory:")
    for k, v in sorted(territory_counts.items()):
        print(f"  {k}: {v}")
    print()
    print("Hits by suggested_action:")
    for k, v in sorted(action_counts.items()):
        print(f"  {k}: {v}")
    print()
    print(f"Top 10 files by hit count (out of {len(file_counts)} files with hits):")
    for path, count in file_counts.most_common(10):
        print(f"  {count:4d}  {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
