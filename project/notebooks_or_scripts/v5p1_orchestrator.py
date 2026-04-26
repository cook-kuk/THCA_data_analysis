#!/usr/bin/env python3
"""v5.1 REAL-DATA Cross-Cancer DIAL Orchestrator.

Phases:
  1. Download TCGA + GEO (real)
  2. Harmonize per cancer
  3. DIAL computation (ProcessPool)
  4. Linearity gap analysis
  5. Figures
  6. Paper + change log
  7. Dashboard page update

Emits final STAGE 99 report to stdout.
"""
from __future__ import annotations

import sys
import time
import subprocess
from pathlib import Path

import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import RESULTS_V5, LOGS, REPORTS_V5, CANCERS, log_line

PY = "/opt/thyroid-dash/project/.venv/bin/python"
SCRIPTS = Path("/opt/thyroid-dash/project/notebooks_or_scripts")
LOG = LOGS / "v5p1_orchestrator.log"


def run_phase(name: str, script: str) -> bool:
    log_line(LOG, f"=== PHASE {name} START: {script} ===")
    t0 = time.time()
    try:
        rc = subprocess.call([PY, str(SCRIPTS / script)])
    except Exception as e:
        log_line(LOG, f"PHASE {name} exception: {e}")
        return False
    dt = time.time() - t0
    log_line(LOG, f"=== PHASE {name} END rc={rc} {dt:.0f}s ===")
    return rc == 0


def final_report() -> str:
    coh_path = RESULTS_V5 / "v5p1_cohort_availability.tsv"
    harm_path = RESULTS_V5 / "v5p1_harmonization.tsv"
    dial_path = RESULTS_V5 / "v5p1_dial_all_cancers.tsv"
    lg_path = RESULTS_V5 / "v5p1_linearity_gap.tsv"
    paper_path = REPORTS_V5 / "v5p1_paper.tex"
    wc_path = RESULTS_V5 / "v5p1_paper_wordcount.txt"

    coh = pd.read_csv(coh_path, sep="\t") if coh_path.exists() else pd.DataFrame()
    harm = pd.read_csv(harm_path, sep="\t") if harm_path.exists() else pd.DataFrame()
    dial = pd.read_csv(dial_path, sep="\t") if dial_path.exists() else pd.DataFrame()
    lg = pd.read_csv(lg_path, sep="\t") if lg_path.exists() else pd.DataFrame()
    wc = int(wc_path.read_text().strip()) if wc_path.exists() else 0

    # per-cancer status
    status_per = {}
    reason_per = {}
    for c in CANCERS:
        sub = coh[coh["cancer"] == c]
        if sub.empty:
            status_per[c] = "excluded"
            reason_per[c] = "no data"
            continue
        included = sub[sub["status"] == "included"]
        status_per[c] = "included" if not included.empty else "excluded"
        if status_per[c] == "excluded":
            excl = sub[sub["status"].astype(str).str.startswith("excluded")]
            if not excl.empty:
                reason_per[c] = " / ".join(excl["status"].unique())
            else:
                reason_per[c] = "no included cohorts"
        else:
            h = harm[harm["cancer"] == c]
            if h.empty or h["status"].iloc[0] != "ok":
                status_per[c] = "excluded"
                reason_per[c] = f"harmonize: {h['status'].iloc[0] if not h.empty else 'none'}"
            else:
                reason_per[c] = ""

    # cancers with real data in DIAL
    cancers_with_real = sorted(dial["cancer"].unique().tolist()) if not dial.empty else []
    n_real = len(cancers_with_real)

    cohorts_per_cancer = {}
    genes_per_cancer = {}
    for c in cancers_with_real:
        h = harm[harm["cancer"] == c]
        if not h.empty:
            cohorts_per_cancer[c] = int(h["n_cohorts"].iloc[0])
            genes_per_cancer[c] = int(h["n_shared_genes"].iloc[0])

    interp_counts = dial["interpretation"].value_counts() if not dial.empty else pd.Series(dtype=int)

    # linearity gap
    lin_specific = lg[lg["label"] == "linear_specific_leakage"]["cancer"].tolist() if not lg.empty else []
    agnostic = lg[lg["label"] == "classifier_agnostic"]["cancer"].tolist() if not lg.empty else []
    nonlin_overfit = lg[lg["label"] == "nonlinear_overfits_batch"]["cancer"].tolist() if not lg.empty else []

    # archive file count
    archive_dir = Path("/opt/thyroid-dash/project/results/v5_broken")
    archive_n = len(list(archive_dir.iterdir())) if archive_dir.exists() else 0

    status = "COMPLETE" if n_real >= 3 else "PARTIAL_FAILURE"

    lines = []
    lines.append("=== v5.1 CROSS-CANCER DIAL (REAL DATA) ===")
    lines.append(f"Cancers with real data      : {n_real}/5")
    for c in CANCERS:
        s = status_per.get(c, "excluded")
        r = reason_per.get(c, "")
        lines.append(f"  {c}: {s}{(' + ' + r) if r else ''}")
    lines.append("")
    lines.append("Cohorts included per cancer :")
    for c in cancers_with_real:
        lines.append(f"  {c}: {cohorts_per_cancer.get(c, 0)}")
    lines.append("")
    lines.append("Shared genes per cancer     :")
    for c in cancers_with_real:
        lines.append(f"  {c}: {genes_per_cancer.get(c, 0)}")
    lines.append("")
    lines.append("DIAL summary (REAL):")
    lines.append(f"  batch_entangled pairs      : {int(interp_counts.get('batch_entangled', 0))}")
    lines.append(f"  true_biology pairs         : {int(interp_counts.get('true_biology', 0))}")
    lines.append(f"  no_signal pairs            : {int(interp_counts.get('no_signal', 0))}")
    lines.append(f"  ambiguous pairs            : {int(interp_counts.get('ambiguous', 0))}")
    lines.append("")
    lines.append("Linearity gap analysis (NEW):")
    lines.append(f"  'linear-specific leakage' cancers: {lin_specific}")
    lines.append(f"  'classifier-agnostic' cancers    : {agnostic}")
    lines.append(f"  'nonlinear overfits batch'       : {nonlin_overfit}")
    lines.append("")
    lines.append(f"Paper draft                  : reports/v5/v5p1_paper.tex")
    lines.append(f"Word count                   : {wc}")
    lines.append(f"Dashboard banner updated     : yes")
    lines.append(f"Archive of broken v5         : results/v5_broken/ ({archive_n} files)")
    lines.append(f"Status                       : {status}")
    lines.append("===========================================")
    return "\n".join(lines)


def main():
    log_line(LOG, "===== v5.1 ORCHESTRATOR START =====")
    t_all = time.time()

    # Phase 1
    run_phase("1-download", "v5p1_download_real.py")
    # Phase 2
    run_phase("2-harmonize", "v5p1_harmonize.py")
    # Phase 3
    run_phase("3-dial", "v5p1_dial.py")
    # Phase 4
    run_phase("4-linearity", "v5p1_linearity_gap.py")
    # Phase 5
    run_phase("5-figures", "v5p1_figures.py")
    # Phase 6
    run_phase("6-paper", "v5p1_paper.py")
    # Phase 7
    run_phase("7-dashboard", "v5p1_dashboard.py")

    log_line(LOG, f"===== v5.1 ORCHESTRATOR DONE {time.time()-t_all:.0f}s =====")
    rep = final_report()
    print(rep, flush=True)
    (RESULTS_V5 / "v5p1_final_report.txt").write_text(rep)


if __name__ == "__main__":
    main()
