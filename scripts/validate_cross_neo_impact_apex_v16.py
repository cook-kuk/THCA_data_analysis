#!/usr/bin/env python3
"""Validate CROSS-Neo impact apex v16 outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/impact_apex_v16_2026_05_11"
HUB = ROOT / "project/papers_hub_2026_05_04"
ASSETS = HUB / "assets/cross_neo_impact_apex_v16"
HTML = HUB / "cross_neo_impact_apex_v16.html"


def main() -> None:
    summary = json.loads((BASE / "impact_apex_v16_summary.json").read_text())
    assert summary["task_count"] == 120
    assert summary["control_count"] == 20
    assert summary["top96_hits"] == 91
    assert summary["top10_hits"] == 10
    assert summary["risk_sum"] < 0.15
    assert summary["execution_wells"] == 96
    assert summary["caption_count"] == 8
    assert summary["live_deploy_ok"] is True

    checks = {
        "claim_ladder_v16.tsv": 6,
        "rebuttal_moat_v16.tsv": 8,
        "decision_grid_v16.tsv": 6,
        "evidence_bridge_v16.tsv": 6,
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
