"""Pattern 5 -- Honest Failure Catalog.

Every research sprint accumulates failed approaches. The instinct is to
delete them. The discipline is to log them with enough context that future
agents (or future-you) don't repeat the mistake.

Used in v3 -> v17p35 to catalog things like:
  * "PRE-2 ComBat under LODO produced DIAL=0.000 -- abandon non-LODO claim"
  * "v1 TERT recovery missed thca_tcga_pub -- check ALL cBio thyroid studies"
  * "A3 5-cohort meta-analysis returned empty -- panels not aligned, not a bug"
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class FailureRecord:
    failure_id: str
    title: str
    what_was_tried: str
    why_it_failed: str
    next_action: str = ""
    severity: str = "MEDIUM"  # HIGH | MEDIUM | LOW
    artifacts: list[str] = field(default_factory=list)
    timestamp: float = field(default_factory=time.time)


class FailureCatalog:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def log(self, failure: FailureRecord) -> None:
        with self.path.open("a") as f:
            f.write(json.dumps(failure.__dict__) + "\n")

    def list(self, severity: str | None = None) -> list[FailureRecord]:
        if not self.path.exists():
            return []
        out: list[FailureRecord] = []
        for line in self.path.read_text().splitlines():
            if not line.strip():
                continue
            d = json.loads(line)
            if severity and d.get("severity") != severity:
                continue
            out.append(FailureRecord(**d))
        return out

    def render_markdown(self) -> str:
        records = self.list()
        if not records:
            return "_No failures logged._"
        lines = ["# Failure catalog", ""]
        for r in records:
            lines.append(f"## {r.failure_id}: {r.title}  _(sev: {r.severity})_")
            lines.append(f"**Tried:** {r.what_was_tried}")
            lines.append(f"**Why it failed:** {r.why_it_failed}")
            if r.next_action:
                lines.append(f"**Next action:** {r.next_action}")
            if r.artifacts:
                lines.append(f"**Artifacts:** " + ", ".join(r.artifacts))
            lines.append("")
        return "\n".join(lines)
