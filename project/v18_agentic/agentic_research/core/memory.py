"""ProjectMemory — file-backed audit + structured fact store.

Two layers:
  * `record()` writes per-iteration artifacts to disk for the audit trail.
  * `fact()` / `recall()` store interpretable claims (orthogonality, AUC,
    survival p-value) keyed by topic for future agents to pull.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class Fact:
    topic: str
    claim: str
    evidence: dict[str, Any] = field(default_factory=dict)
    confidence: str = "MEDIUM"  # HIGH | MEDIUM | LOW | FLAGGED
    timestamp: float = field(default_factory=time.time)


class ProjectMemory:
    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.audit_dir = self.root / "audit"
        self.audit_dir.mkdir(exist_ok=True)
        self.facts_file = self.root / "facts.jsonl"
        self._facts_cache: list[Fact] | None = None

    def record(self, name: str, artifacts: dict[str, Any]) -> Path:
        path = self.audit_dir / f"{int(time.time() * 1000)}_{name}.json"
        path.write_text(json.dumps(artifacts, indent=2, default=str))
        return path

    def fact(self, fact: Fact) -> None:
        with self.facts_file.open("a") as f:
            f.write(json.dumps(fact.__dict__) + "\n")
        self._facts_cache = None

    def recall(self, topic: str, min_confidence: str = "MEDIUM") -> list[Fact]:
        order = {"HIGH": 3, "MEDIUM": 2, "LOW": 1, "FLAGGED": 0}
        threshold = order.get(min_confidence, 2)
        if self._facts_cache is None:
            self._facts_cache = []
            if self.facts_file.exists():
                for line in self.facts_file.read_text().splitlines():
                    if line.strip():
                        d = json.loads(line)
                        self._facts_cache.append(Fact(**d))
        return [f for f in self._facts_cache if f.topic == topic and order.get(f.confidence, 0) >= threshold]
