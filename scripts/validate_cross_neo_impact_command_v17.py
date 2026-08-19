#!/usr/bin/env python3
"""Validate CROSS-Neo impact command v17 outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/impact_command_v17_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
ASSETS = HUB / "assets/cross_neo_impact_command_v17"
HTML = HUB / "cross_neo_impact_command_v17.html"


def main() -> None:
    summary = json.loads((BASE / "impact_command_v17_summary.json").read_text())
    assert summary["task_count"] == 120
    assert summary["control_count"] == 20
    assert summary["top96_hits"] == 91
    assert summary["top10_hits"] == 10
    assert summary["risk_sum"] < 0.15
    assert summary["execution_wells"] == 96
    assert summary["claim_lock_count"] == 6
    assert summary["live_deploy_ok"] is True

    checks = {
        "executive_scoreboard_v17.tsv": 7,
        "claim_lock_matrix_v17.tsv": 6,
        "rebuttal_register_v17.tsv": 8,
        "action_queue_v17.tsv": 6,
        "evidence_spine_v17.tsv": 6,
    }
    for filename, expected_rows in checks.items():
        df = pd.read_csv(BASE / filename, sep="\t")
        assert len(df) >= expected_rows, (filename, len(df), expected_rows)

    assert HTML.exists() and HTML.stat().st_size > 1000
    pngs = sorted(ASSETS.glob("*.png"))
    assert len(pngs) >= 6, len(pngs)
    for path in pngs:
        image = Image.open(path).convert("RGB")
        nonflat = sum(hi - lo for lo, hi in image.getextrema())
        assert nonflat > 10, path

    print("validation_ok", summary["task_count"], summary["top96_hits"], summary["risk_sum"], len(pngs), summary["live_deploy_ok"])


if __name__ == "__main__":
    main()
