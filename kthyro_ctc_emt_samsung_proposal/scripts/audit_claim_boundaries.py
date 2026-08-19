#!/usr/bin/env python3
"""Simple conservative claim-boundary audit for generated reports."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "outputs" / "reports"
LOG = ROOT / "logs" / "claim_boundary_audit.log"

FORBIDDEN_PATTERNS = [
    "already predicts recurrence in our cohort",
    "proves CTC biology",
    "will work in all early PTC",
    "ready clinical diagnostic",
    "vendor technology is the novelty",
    "clinical deployment evidence",
    "proven recurrence predictor",
]

ALLOWED_CONTEXTS = [
    "Forbidden",
    "금지",
    "does not",
    "none of this",
    "not ",
    "아니다",
    "주장하지",
    "claim 없음",
]


def line_is_boundary(line: str) -> bool:
    lower = line.lower()
    return any(ctx.lower() in lower for ctx in ALLOWED_CONTEXTS)


def main() -> int:
    findings = []
    for path in sorted(REPORTS.glob("*.md")):
        in_forbidden_section = False
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip().lower()
            if stripped.startswith("##") and "forbidden" in stripped:
                in_forbidden_section = True
            elif stripped.startswith("##"):
                in_forbidden_section = False
            if in_forbidden_section:
                continue
            low = line.lower()
            for pat in FORBIDDEN_PATTERNS:
                if pat.lower() in low and not line_is_boundary(line):
                    findings.append((path.name, n, pat, line.strip()))

    LOG.parent.mkdir(parents=True, exist_ok=True)
    if findings:
        rows = ["FAIL: possible overclaim patterns"]
        for item in findings:
            rows.append("\t".join(map(str, item)))
        LOG.write_text("\n".join(rows) + "\n", encoding="utf-8")
        print(LOG)
        return 1

    LOG.write_text("PASS: no unbounded forbidden claim patterns detected.\n", encoding="utf-8")
    print(LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
