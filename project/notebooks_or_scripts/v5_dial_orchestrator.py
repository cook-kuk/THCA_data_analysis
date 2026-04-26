#!/usr/bin/env python3
"""v5 DIAL Cross-Cancer Orchestrator.

Runs phases 1-7 sequentially with timestamped logging. Continues on phase
failure. Health-check at the end.
"""
from __future__ import annotations

import os
# Limit BLAS thread oversubscription when 5 workers run in parallel.
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(k, "2")

import sys
import time
import importlib
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from v5_dial_common import LOGS, PROJECT, RESULTS_V5, REPORTS_V5, PAGES, log_line, peak_rss_gb

LOGFILE = LOGS / "v5_dial_run.log"


def run_phase(name: str, module: str):
    log_line(LOGFILE, f"================ START {name} ({module}) ================")
    t0 = time.time()
    try:
        if module in sys.modules:
            importlib.reload(sys.modules[module])
        else:
            importlib.import_module(module)
        m = sys.modules[module]
        m.main()
        dt = time.time() - t0
        log_line(LOGFILE, f"================ OK    {name} in {dt:.1f}s RSS={peak_rss_gb():.2f}GB ================")
        return True
    except Exception as e:
        dt = time.time() - t0
        tb = traceback.format_exc(limit=10)
        log_line(LOGFILE, f"================ FAIL  {name} after {dt:.1f}s: {type(e).__name__}: {e}")
        log_line(LOGFILE, tb)
        log_line(LOGFILE, f"================ CONTINUE with next phase ================")
        return False


def health_check() -> dict:
    """Verify outputs; return dict with pass/fail summary."""
    import pandas as pd
    out = {}
    # aggregated DIAL
    tsv = RESULTS_V5 / "v5_dial_all_cancers.tsv"
    out["v5_dial_all_cancers_exists"] = tsv.exists()
    if tsv.exists():
        df = pd.read_csv(tsv, sep="\t")
        out["v5_dial_rows"] = len(df)
        out["v5_dial_rows_ok"] = len(df) >= 10
    else:
        out["v5_dial_rows"] = 0
        out["v5_dial_rows_ok"] = False

    # paper
    tex = REPORTS_V5 / "v5_dial_paper.tex"
    out["v5_dial_paper_exists"] = tex.exists()
    if tex.exists():
        wc = len(tex.read_text().split())
        out["v5_dial_paper_words"] = wc
        out["v5_dial_paper_words_ok"] = wc > 3000
    else:
        out["v5_dial_paper_words"] = 0
        out["v5_dial_paper_words_ok"] = False

    # dashboard page
    page = PAGES / "v4a_dial_cross_cancer.html"
    out["v5_dial_page_exists"] = page.exists()
    out["v5_dial_page_bytes"] = page.stat().st_size if page.exists() else 0

    # figures
    fig_dir = PROJECT / "reports" / "html" / "figs_interactive" / "v5"
    out["v5_dial_fig_count"] = len(list(fig_dir.glob("v5_dial_fig*.html")))

    # theory
    out["v5_dial_theory_exists"] = (REPORTS_V5 / "v5_dial_theory.tex").exists()

    # try curl
    try:
        import subprocess
        r = subprocess.run(
            ["curl", "-sI", "--max-time", "5",
             "http://127.0.0.1:8012/reports/html/pages/v4a_dial_cross_cancer.html"],
            capture_output=True, timeout=10,
        )
        first_line = r.stdout.decode().splitlines()[0] if r.stdout else ""
        out["curl_status"] = first_line
    except Exception as e:
        out["curl_status"] = f"error: {e}"

    return out


def main():
    log_line(LOGFILE, "############# v5 DIAL ORCHESTRATOR START #############")
    t0 = time.time()

    phases = [
        ("PHASE 1 download",   "v5_dial_phase1_download"),
        ("PHASE 2 harmonize",  "v5_dial_phase2_harmonize"),
        ("PHASE 3 DIAL compute", "v5_dial_phase3_compute"),
        # Phase 4 (theory) is hand-written; verify it exists
        ("PHASE 4 theory check", "v5_dial_phase4_theory_check"),
        ("PHASE 5 figures",    "v5_dial_phase5_figures"),
        ("PHASE 6 paper",      "v5_dial_phase6_paper"),
        ("PHASE 7 dashboard",  "v5_dial_phase7_dashboard"),
    ]

    verdicts = []
    for name, mod in phases:
        ok = run_phase(name, mod)
        verdicts.append((name, ok))

    dt = time.time() - t0
    log_line(LOGFILE, f"############# ALL PHASES FINISHED in {dt/60:.1f} min #############")

    hc = health_check()
    log_line(LOGFILE, f"HEALTH CHECK: {hc}")

    any_fail = any(not ok for _, ok in verdicts)
    verdict_str = "PARTIAL" if any_fail else "COMPLETE"
    log_line(LOGFILE, f"ORCHESTRATOR VERDICT = {verdict_str}")
    # final summary
    log_line(LOGFILE, "Phase verdicts: " + ", ".join(f"{n}={'ok' if o else 'FAIL'}" for n, o in verdicts))


if __name__ == "__main__":
    main()
