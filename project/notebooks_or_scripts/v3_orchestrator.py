"""v3 STEP 10 — Orchestrator.

P0 (blocking)         : v3_honesty_audit.py — all 6 tasks (run internally in
                        parallel via ProcessPoolExecutor)
P1 (parallel)         : v3_panel_size_curve, v3_three_class, v3_survival,
                        v3_bethesda_sim
P2 (parallel, optional): v3_multimodal
P3 (sequential)       : v3_build_pages, v3_write_reports

STAGE 99 block emitted at the end with real numbers.
"""
from __future__ import annotations

import json
import logging
import multiprocessing as mp
import os
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
NOTEBOOKS = ROOT / "notebooks_or_scripts"
LOGS = ROOT / "logs"
LOGS.mkdir(parents=True, exist_ok=True)

PYTHON = "/opt/thyroid-dash/project/.venv/bin/python"

LOG_FILE = LOGS / "v3_run.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, mode="w"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("v3")


def _run_script(script: str, extra: list | None = None) -> dict:
    """Spawn a python process running the script; capture stdout+stderr."""
    cmd = [PYTHON, str(NOTEBOOKS / script)] + (extra or [])
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=2400)
        return dict(script=script, ok=(proc.returncode == 0),
                     elapsed=time.time() - t0,
                     stdout=proc.stdout[-4000:], stderr=proc.stderr[-4000:],
                     returncode=proc.returncode)
    except Exception as e:
        return dict(script=script, ok=False,
                     elapsed=time.time() - t0, error=str(e),
                     tb=traceback.format_exc())


def _run_honesty_task(key: str) -> dict:
    """Run a single honesty audit task in-process via subprocess."""
    return _run_script("v3_honesty_audit.py", ["--task", key])


def run_p0_honesty():
    log.info("P0: honesty audit (6 tasks in parallel)")
    keys = ["A", "B", "C", "D", "E", "F"]
    max_workers = min(mp.cpu_count(), 4)
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(_run_honesty_task, k): k for k in keys}
        for fut in as_completed(futures):
            k = futures[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = dict(script="v3_honesty_audit.py", key=k, ok=False,
                         error=str(e))
            r["key"] = k
            log.info(f"  task {k}: ok={r.get('ok')} "
                      f"elapsed={r.get('elapsed', 0):.1f}s")
            if not r.get("ok"):
                log.warning(f"  task {k} stderr: {r.get('stderr', '')[-500:]}")
            results.append(r)
    return results


def run_p1_parallel():
    log.info("P1: parallel scripts (panel, three_class, survival, bethesda)")
    scripts = ["v3_panel_size_curve.py", "v3_three_class.py",
                "v3_survival.py", "v3_bethesda_sim.py"]
    max_workers = min(mp.cpu_count(), 4)
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(_run_script, s): s for s in scripts}
        for fut in as_completed(futures):
            s = futures[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = dict(script=s, ok=False, error=str(e))
            log.info(f"  {s}: ok={r.get('ok')} "
                      f"elapsed={r.get('elapsed', 0):.1f}s")
            if not r.get("ok"):
                log.warning(f"  {s} stderr: {r.get('stderr', '')[-500:]}")
            results.append(r)
    return results


def run_p2_multimodal():
    log.info("P2: multimodal (best-effort)")
    r = _run_script("v3_multimodal.py")
    log.info(f"  v3_multimodal.py: ok={r.get('ok')} "
              f"elapsed={r.get('elapsed', 0):.1f}s")
    if not r.get("ok"):
        log.warning(f"  stderr: {r.get('stderr', '')[-500:]}")
    return [r]


def run_p3_sequential():
    log.info("P3: build pages + write reports")
    results = []
    for s in ["v3_build_pages.py", "v3_write_reports.py"]:
        r = _run_script(s)
        log.info(f"  {s}: ok={r.get('ok')} elapsed={r.get('elapsed', 0):.1f}s")
        if not r.get("ok"):
            log.warning(f"  {s} stderr: {r.get('stderr', '')[-500:]}")
        results.append(r)
    return results


# ---------- STAGE 99 ----------
def _safe_tsv(p):
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def _safe_json(p):
    try:
        return json.loads(Path(p).read_text())
    except Exception:
        return {}


def emit_stage99(timings):
    RT = ROOT / "results" / "tables"
    lev = _safe_tsv(RT / "v3_leakage_curve.tsv")
    mapk = _safe_tsv(RT / "v3_mapk_ablation.tsv")
    perm = _safe_tsv(RT / "v3_permutation_null.tsv")
    ident = _safe_tsv(RT / "v3_dataset_identifiability.tsv")
    lodo = _safe_tsv(RT / "v3_lodo.tsv")
    panel = _safe_tsv(RT / "v3_panel_size_curve.tsv")
    best = _safe_json(RT / "v3_panel_best.json")
    tri = _safe_tsv(RT / "v3_three_class_metrics.tsv")
    cox = _safe_tsv(RT / "v3_survival_cox.tsv")
    beth = _safe_tsv(RT / "v3_bethesda_surgery_reduction.tsv")
    meth = _safe_tsv(RT / "v3_methylation_results.tsv")

    def _pick_lev(k):
        r = lev[lev.k_removed == k]
        return float(r.iloc[0]["mean_auc"]) if not r.empty else float("nan")

    auc_k0 = _pick_lev(0)
    auc_k10 = _pick_lev(10)
    auc_k30 = _pick_lev(30)

    mapk_full = (float(mapk[mapk.variant == "full_TierA67_clean"].iloc[0][
        "mean_cv_auc"]) if not mapk.empty else float("nan"))
    mapk_minus = (float(mapk[mapk.variant == "minus_MAPK_OUTPUT"].iloc[0][
        "mean_cv_auc"]) if not mapk.empty else float("nan"))
    perm_p = float(perm.iloc[0]["empirical_p"]) if not perm.empty else float("nan")
    real_auc = float(perm.iloc[0]["real_auc"]) if not perm.empty else float("nan")
    ident_mean = float(ident["ovr_auc_mean"].mean()) if not ident.empty else float("nan")

    logreg_macro = float(tri[(tri.model == "logreg") &
                              (tri.class_name == "MACRO")].iloc[0]["ovr_auc"]) \
        if not tri.empty else float("nan")

    beth_rows = beth.to_dict(orient="records") if not beth.empty else []

    lines = [
        "",
        "=" * 78,
        "STAGE 99 — v3 RESEARCH PHASE COMPLETE",
        "=" * 78,
        f"build_time_utc         : {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "HONESTY AUDIT (page 19)",
        f"  TierA67_clean TCGA AUC @ k=0       : {auc_k0:.3f}",
        f"  TierA67_clean TCGA AUC @ k=10      : {auc_k10:.3f}",
        f"  TierA67_clean TCGA AUC @ k=30      : {auc_k30:.3f}",
        f"  full TierA67_clean CV AUC          : {mapk_full:.3f}",
        f"  minus MAPK_OUTPUT CV AUC           : {mapk_minus:.3f}",
        f"  permutation null p (1000 shuffles) : {perm_p:.4f}  (real AUC {real_auc:.3f})",
        f"  dataset identifiability mean OvR   : {ident_mean:.3f}",
        "",
        "LODO (page 19)",
    ]
    if not lodo.empty:
        for _, r in lodo.iterrows():
            lines.append(
                f"  {r['external']:10s} [{r['mode']:25s}] "
                f"AUC={r['auc']:.3f}  PR={r['pr_auc']:.3f}  "
                f"bACC={r['bacc']:.3f}  NPV@Se95={r['npv_at_se95']:.3f}")
    else:
        lines.append("  (empty)")

    lines += [
        "",
        "PANEL SIZE CURVE (page 20)",
        f"  best strategy : {best.get('strategy', 'n/a')}",
        f"  best k        : {best.get('k', 'n/a')}",
        f"  TCGA 5-fold AUC           : {best.get('tcga_auc', float('nan')):.3f}",
        f"  GSE27155 isotonic AUC     : {best.get('gse27155_auc', float('nan')):.3f}",
        "",
        "THREE-CLASS (page 21)",
        f"  LogReg macro OvR AUC      : {logreg_macro:.3f}",
        "",
        "SURVIVAL (page 22)",
        f"  Cox rows                  : {len(cox)}",
        "",
        "BETHESDA SIM (page 23)",
    ]
    for r in beth_rows:
        lines.append(
            f"  prev={r['prev']:.0%}  thr={r['chosen_threshold']:.2f}  "
            f"Se={r['sens']:.3f}  Sp={r['spec']:.3f}  "
            f"NPV={r['npv']:.3f}  surgery_reduction={r['reduction_vs_treat_all'] * 100:.0f}%")
    if not beth_rows:
        lines.append("  (empty)")

    lines += [
        "",
        "MULTIMODAL (page 24)",
        f"  methylation GSE97466 AUC  : "
        f"{(float(meth.iloc[0]['mean_cv_auc']) if not meth.empty else float('nan')):.3f}",
        "  fusion                    : skipped_no_raw",
        "  SCNA                      : skipped_no_raw",
        "",
        "TIMINGS",
    ]
    for stage, t in timings:
        lines.append(f"  {stage:10s}  {t:.1f}s")
    lines += [
        "",
        "=" * 78,
    ]
    out = "\n".join(lines)
    print(out)
    LOG_FILE_99 = LOGS / "v3_stage99.txt"
    LOG_FILE_99.write_text(out + "\n", encoding="utf-8")
    return out


def main():
    timings = []
    t0 = time.time()
    try:
        p0 = run_p0_honesty()
    finally:
        timings.append(("P0", time.time() - t0))

    t1 = time.time()
    try:
        p1 = run_p1_parallel()
    finally:
        timings.append(("P1", time.time() - t1))

    t2 = time.time()
    try:
        p2 = run_p2_multimodal()
    finally:
        timings.append(("P2", time.time() - t2))

    t3 = time.time()
    try:
        p3 = run_p3_sequential()
    finally:
        timings.append(("P3", time.time() - t3))

    (LOGS / "v3_orchestrator_results.json").write_text(
        json.dumps(dict(p0=p0, p1=p1, p2=p2, p3=p3), indent=2, default=str))

    emit_stage99(timings)


if __name__ == "__main__":
    main()
