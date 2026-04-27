#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from v17p35_common import LOG35, ROOT, log_line

TIER1 = [
    "v17p35_FIX4_gsea_compare.py",
]


def run_one(script: str) -> tuple[str, int]:
    t0 = time.time()
    path = ROOT / "notebooks_or_scripts" / script
    log_line(LOG35, f"start {script}")
    p = subprocess.run(
        [str(ROOT / ".venv" / "bin" / "python"), str(path)],
        cwd=str(ROOT.parent),
        capture_output=True,
        text=True,
    )
    if p.stdout:
        log_line(LOG35, f"{script} stdout {p.stdout[-1200:]}")
    if p.stderr:
        log_line(LOG35, f"{script} stderr {p.stderr[-1200:]}")
    log_line(LOG35, f"end {script} code={p.returncode} sec={time.time()-t0:.1f}")
    return script, p.returncode


def main() -> None:
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(run_one, s) for s in TIER1]
        for fut in as_completed(futs):
            fut.result()


if __name__ == "__main__":
    main()
