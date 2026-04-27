#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from v17p3_common import *  # noqa: F401,F403

RES35 = ROOT / "results" / "v17p35"
TAB35 = RES35 / "tables"
FIG35 = RES35 / "figs"
RPT35 = ROOT / "reports" / "v17p35"
LOG35 = RES35 / "orchestrator.log"

for d in (TAB35, FIG35, RPT35):
    d.mkdir(parents=True, exist_ok=True)
