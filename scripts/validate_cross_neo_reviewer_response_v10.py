#!/usr/bin/env python3
"""Validate CROSS-Neo reviewer response v10 outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/reviewer_response_v10_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
ASSETS = HUB / "assets/cross_neo_reviewer_response_v10"
HTML = HUB / "cross_neo_reviewer_response_v10.html"


def main() -> None:
    summary = json.loads((BASE / "reviewer_response_v10_summary.json").read_text())
    assert summary["top96_hits"] == 91
    assert summary["top10_hits"] == 10
    assert summary["risk_sum"] < 0.15
    assert summary["execution_wells"] == 96
    assert summary["order_lines"] == 114
    assert summary["live_deploy_ok"] is True

    checks = {
        "reviewer_task_board_v10.tsv": 20,
        "control_manifest_v10.tsv": 8,
        "failure_mode_registry_v10.tsv": 8,
        "reviewer_objection_matrix_v10.tsv": 5,
        "five_slide_scaffold_v10.tsv": 5,
    }
    for filename, min_rows in checks.items():
        df = pd.read_csv(BASE / filename, sep="\t")
        assert len(df) >= min_rows, (filename, len(df), min_rows)

    assert HTML.exists() and HTML.stat().st_size > 1000
    pngs = sorted(ASSETS.glob("*.png"))
    assert len(pngs) >= 7, len(pngs)
    for path in pngs:
        image = Image.open(path).convert("RGB")
        nonflat = sum(hi - lo for lo, hi in image.getextrema())
        assert nonflat > 10, path

    print("validation_ok", summary["top96_hits"], summary["top10_hits"], summary["risk_sum"], len(pngs), summary["live_deploy_ok"])


if __name__ == "__main__":
    main()
