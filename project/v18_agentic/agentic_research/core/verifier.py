"""Verifier — judge an artifact against the goal's success criteria.

Default implementation runs deterministic checks (presence of keys, numeric
ranges, file existence). Subclasses can plug in LLM-as-judge.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable


@dataclass
class VerificationResult:
    passed: bool
    failures: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        return {"passed": self.passed, "failures": self.failures, "notes": self.notes}


Check = Callable[[dict[str, Any]], tuple[bool, str]]


class Verifier:
    def __init__(self, checks: list[Check] | None = None):
        self.checks = checks or []

    async def verify(self, artifacts: dict[str, Any], goal) -> VerificationResult:
        result = VerificationResult(passed=True)
        for check in self.checks:
            ok, msg = check(artifacts)
            if not ok:
                result.passed = False
                result.failures.append(msg)
        # Cross-check goal.success_criteria — these are free-text; record them
        for criterion in getattr(goal, "success_criteria", []):
            result.notes.append(f"criterion: {criterion}")
        return result
