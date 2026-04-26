#!/usr/bin/env python3
"""v5 DIAL Phase 4 — verify theory .tex exists (hand-written, not generated).

The Lemma and corollary text are hand-written in
reports/v5/v5_dial_theory.tex and committed as part of the sprint. This
phase just verifies presence and minimum length so the orchestrator can
continue.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from v5_dial_common import REPORTS_V5, LOGS, log_line

LOGFILE = LOGS / "v5_dial_run.log"


def main():
    p = REPORTS_V5 / "v5_dial_theory.tex"
    if not p.exists():
        raise FileNotFoundError(f"theory file missing: {p}")
    n_chars = len(p.read_text())
    n_words = len(p.read_text().split())
    log_line(LOGFILE, f"PHASE4 theory OK words={n_words} chars={n_chars}")
    if n_words < 200:
        raise ValueError(f"theory too short: {n_words} words")


if __name__ == "__main__":
    main()
