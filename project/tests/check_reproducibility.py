#!/usr/bin/env python3
"""F10 reproducibility check.

Compares the freshly generated v5p1_dial_all_cancers.tsv against the
checked-in reference (tests/expected_dial.tsv).

Pass criteria (per row, joined on (cancer, classifier)):
  - |dial_observed - dial_expected| < 0.01
  - interpretation_observed == interpretation_expected

Exit code 0 on success, 1 on any mismatch (with a per-row diff report).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
OBSERVED = PROJECT / "results" / "v5" / "v5p1_dial_all_cancers.tsv"
EXPECTED = PROJECT / "tests" / "expected_dial.tsv"

DIAL_TOL = 0.01
KEY = ["cancer", "classifier"]


def _load(p: Path) -> pd.DataFrame:
    if not p.exists():
        sys.stderr.write(f"[FAIL] missing file: {p}\n")
        sys.exit(2)
    df = pd.read_csv(p, sep="\t")
    needed = {"cancer", "classifier", "dial", "interpretation"}
    missing = needed - set(df.columns)
    if missing:
        sys.stderr.write(f"[FAIL] {p} missing columns: {missing}\n")
        sys.exit(2)
    return df


def main() -> int:
    obs = _load(OBSERVED)
    exp = _load(EXPECTED)

    if len(obs) != len(exp):
        print(f"[WARN] row-count differs: observed={len(obs)} expected={len(exp)}",
              file=sys.stderr)

    merged = exp.merge(obs, on=KEY, how="outer", suffixes=("_exp", "_obs"),
                       indicator=True)

    failures: list[str] = []
    only_exp = merged[merged["_merge"] == "left_only"]
    only_obs = merged[merged["_merge"] == "right_only"]
    for _, r in only_exp.iterrows():
        failures.append(f"missing in observed: {r['cancer']}/{r['classifier']}")
    for _, r in only_obs.iterrows():
        failures.append(f"unexpected in observed: {r['cancer']}/{r['classifier']}")

    both = merged[merged["_merge"] == "both"].copy()
    both["dial_diff"] = (both["dial_obs"] - both["dial_exp"]).abs()

    n_pass = 0
    for _, r in both.iterrows():
        diff = float(r["dial_diff"])
        same_interp = str(r["interpretation_obs"]) == str(r["interpretation_exp"])
        if diff < DIAL_TOL and same_interp:
            n_pass += 1
        else:
            failures.append(
                f"{r['cancer']}/{r['classifier']}: "
                f"dial_diff={diff:.4f} (tol={DIAL_TOL}), "
                f"interp obs={r['interpretation_obs']!r} exp={r['interpretation_exp']!r}"
            )

    print(f"matched rows: {len(both)}, passing: {n_pass}, failing: {len(failures)}")

    if failures:
        print("\nFAILURES:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1

    print("REPRODUCIBILITY CHECK: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
