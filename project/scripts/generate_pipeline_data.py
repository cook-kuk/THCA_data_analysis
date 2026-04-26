#!/usr/bin/env python3
"""Generate THYRAI pipeline.json from drug_discovery_targets.tsv.

Applies conservative stage rules:
  - pChEMBL >= 8 AND literature_count >= 3  -> hit_to_lead
  - otherwise                                -> target_id

All 8 targets remain at target_id or hit_to_lead. No fabricated clinical
progression. Output schema is documented in the rebrand spec.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

PROJECT = Path("/opt/thyroid-dash/project")
TARGETS_TSV = PROJECT / "results/tables/drug_discovery_targets.tsv"
LIT_TSV = PROJECT / "results/tables/drug_discovery_literature.tsv"
REPORT_MD = PROJECT / "reports/biomarker_to_drug_report.md"
OUT_JSON = PROJECT / "reports/html/assets/data/pipeline.json"

INDICATION_MAP = {
    "TACSTD2": "Epithelial-origin thyroid cancer",
    "TMPRSS4": "Invasive papillary thyroid cancer",
    "PLEKHA6": "Dedifferentiated thyroid cancer",
    "CYP1B1":  "BRAF-like thyroid cancer",
    "LDLR":    "Metabolic-axis thyroid cancer",
    "GABRB2":  "Neural-marker thyroid cancer",
    "B3GNT3":  "Glycosylation-shift thyroid cancer",
    "PTPRE":   "Signaling-rewired thyroid cancer",
}


def to_float(v: str) -> float:
    try:
        return float(v) if v not in (None, "", "NA") else 0.0
    except ValueError:
        return 0.0


def to_int(v: str) -> int:
    try:
        return int(float(v)) if v not in (None, "", "NA") else 0
    except ValueError:
        return 0


def count_literature() -> dict:
    counts: dict = {}
    if not LIT_TSV.exists():
        return counts
    with LIT_TSV.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("gene") or row.get("target") or ""
            if not gene:
                continue
            counts[gene] = counts.get(gene, 0) + 1
    return counts


def load_rationales() -> dict:
    """Extract simple per-gene bullet sentences from the report md."""
    rationales: dict = {}
    if not REPORT_MD.exists():
        return rationales
    txt = REPORT_MD.read_text()
    # Lines of the form `- **GENE** — ...` or `### GENE`
    current = None
    for line in txt.splitlines():
        m = re.match(r"^#{1,4}\s+([A-Z][A-Z0-9]{1,12})\b", line)
        if m:
            current = m.group(1)
            continue
        m2 = re.match(r"^-\s+\*\*([A-Z][A-Z0-9]{1,12})\*\*\s*[—\-:]\s*(.+)$", line)
        if m2:
            rationales[m2.group(1)] = m2.group(2).strip()
            continue
        if current and line.strip() and not line.startswith("#"):
            rationales.setdefault(current, line.strip()[:280])
    return rationales


def stage_rule(pchembl: float, literature: int, compounds: int) -> str:
    if pchembl >= 8.0 and literature >= 3:
        return "hit_to_lead"
    return "target_id"


def main() -> int:
    if not TARGETS_TSV.exists():
        raise SystemExit(f"Missing: {TARGETS_TSV}")

    lit_counts = count_literature()
    rationales = load_rationales()

    rows = []
    with TARGETS_TSV.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            gene = r["gene"].strip()
            pchembl = to_float(r.get("best_pchembl", ""))
            compounds = to_int(r.get("n_compounds", "0"))
            classification = r.get("classification", "target").strip()
            is_novel = r.get("is_novel", "False").strip().lower() == "true"
            top_compound = (r.get("top_compound") or "").strip() or None
            lit = lit_counts.get(gene, 0)
            stage = stage_rule(pchembl, lit, compounds)
            rows.append({
                "target": gene,
                "novel": bool(is_novel),
                "classification": classification,
                "indication": INDICATION_MAP.get(gene, "Thyroid cancer subtype"),
                "compounds": compounds,
                "best_pchembl": round(pchembl, 2) if pchembl else None,
                "literature_count": lit,
                "stage": stage,
                "licensing": True,
                "top_compound": top_compound,
                "rationale": rationales.get(gene, ""),
            })

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(rows, indent=2))
    print(f"[pipeline] wrote {len(rows)} targets -> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
