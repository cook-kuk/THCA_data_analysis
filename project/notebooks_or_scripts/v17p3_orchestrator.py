#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from v17p3_common import LOG, ROOT, log_line

TASKS = [
    "v17p3_F1_external_recovery.py",
    "v17p3_F2_msigdb_proper.py",
    "v17p3_F3_clinical_reframe.py",
    "v17p3_A1_scrna_heterogeneity.py",
    "v17p3_A2_dm_in_braf_ras.py",
    "v17p3_A3_drug_response.py",
    "v17p3_A4_pancancer.py",
    "v17p3_A5_genomic_immune.py",
    "v17p3_A6_tf_circuit.py",
]


def run_one(script: str) -> tuple[str, int]:
    t0 = time.time()
    path = ROOT / "notebooks_or_scripts" / script
    log_line(LOG, f"start {script}")
    p = subprocess.run([str(ROOT / ".venv" / "bin" / "python"), str(path)], cwd=str(ROOT.parent), capture_output=True, text=True)
    if p.stdout:
        log_line(LOG, f"{script} stdout {p.stdout[-1200:]}")
    if p.stderr:
        log_line(LOG, f"{script} stderr {p.stderr[-1200:]}")
    log_line(LOG, f"end {script} code={p.returncode} sec={time.time()-t0:.1f}")
    return script, p.returncode


def main() -> None:
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(run_one, s) for s in TASKS[:3]]
        for fut in as_completed(futs):
            fut.result()
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs = [ex.submit(run_one, s) for s in TASKS[3:6]]
        for fut in as_completed(futs):
            fut.result()
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs = [ex.submit(run_one, s) for s in TASKS[6:]]
        for fut in as_completed(futs):
            fut.result()
    run_one("v17p3_build_report.py")
    run_one("v17p3_make_docs.py")


if __name__ == "__main__":
    main()
