#!/usr/bin/env python3
"""TCGA-LUAD pathology-AI audit entry point.

At present this writes an input-status report. Once LUAD UNI features,
manifest, metadata, and OOF predictions exist, this wrapper should be
extended to call the same C1-C5/E1-E8/S1-S5 routines used for THCA.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent


def main() -> None:
    subprocess.run(
        [sys.executable, str(HERE / "audit_tcga_cohort_scaffold.py"), "LUAD"],
        check=True,
    )


if __name__ == "__main__":
    main()
