#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from v17p2_common import LOG, ROOT, log_line

SCRIPTS = {
    "layer1": "v17p2_robustness.py",
    "layer2": "v17p2_clinical.py",
    "layer3": "v17p2_biology.py",
    "layer4": "v17p2_external.py",
    "layer5": "v17p2_dial_method.py",
    "report": "v17p2_build_report.py",
}


def run_script(name: str) -> tuple[str, int, float]:
    path = ROOT / "notebooks_or_scripts" / SCRIPTS[name]
    t0 = time.time()
    log_line(LOG, f"{name} start script={path.name}")
    proc = subprocess.run([str(ROOT / ".venv" / "bin" / "python"), str(path)], cwd=str(ROOT.parent), capture_output=True, text=True)
    dt = time.time() - t0
    if proc.stdout:
        log_line(LOG, f"{name} stdout {proc.stdout[-1000:]}")
    if proc.stderr:
        log_line(LOG, f"{name} stderr {proc.stderr[-1000:]}")
    log_line(LOG, f"{name} end code={proc.returncode} sec={dt:.1f}")
    return name, proc.returncode, dt


def main() -> None:
    results = {}
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs = {ex.submit(run_script, name): name for name in ["layer1", "layer2", "layer3"]}
        for fut in as_completed(futs):
            name, code, dt = fut.result()
            results[name] = code
    if results.get("layer1", 1) == 0:
        results["layer4"] = run_script("layer4")[1]
    if results.get("layer4", 1) == 0:
        results["layer5"] = run_script("layer5")[1]
    if all(v == 0 for v in results.values()):
        results["report"] = run_script("report")[1]
    failed = {k: v for k, v in results.items() if v != 0}
    if failed:
        log_line(LOG, f"orchestrator failed {failed}")
        sys.exit(1)
    log_line(LOG, f"orchestrator complete {results}")


if __name__ == "__main__":
    main()
