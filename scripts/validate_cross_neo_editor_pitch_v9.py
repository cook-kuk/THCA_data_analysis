#!/usr/bin/env python3
"""Validate CROSS-Neo editor pitch v9 outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/editor_pitch_v9_2026_05_10"
HUB = ROOT / "project/papers_hub_2026_05_04"
ASSETS = HUB / "assets/cross_neo_editor_pitch_v9"
HTML = HUB / "cross_neo_editor_pitch_v9.html"


def main() -> None:
    summary = json.loads((BASE / "editor_pitch_v9_summary.json").read_text())
    assert summary["top96_hits"] == 91
    assert summary["top10_hits"] == 10
    assert summary["risk_sum"] < 0.15
    assert summary["execution_wells"] == 96
    assert summary["order_lines"] == 114

    expected_rows = {
        "editor_pitch_deck_v9.tsv": 5,
        "reviewer_objection_response_map_v9.tsv": 5,
        "submission_claim_boundary_v9.tsv": 5,
        "figure_bundle_manifest_v9.tsv": 9,
    }
    for name, min_rows in expected_rows.items():
        df = pd.read_csv(BASE / name, sep="\t")
        assert len(df) >= min_rows, (name, len(df), min_rows)

    assert HTML.exists() and HTML.stat().st_size > 1000
    pngs = sorted(ASSETS.glob("*.png"))
    assert len(pngs) >= 9, len(pngs)
    for path in pngs:
        image = Image.open(path).convert("RGB")
        nonflat = sum(hi - lo for lo, hi in image.getextrema())
        assert nonflat > 10, path

    print(
        "validation_ok",
        summary["top96_hits"],
        summary["top10_hits"],
        summary["risk_sum"],
        len(pngs),
        "live",
        summary["live_deploy_ok"],
    )


if __name__ == "__main__":
    main()
