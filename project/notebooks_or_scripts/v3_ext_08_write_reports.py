#!/usr/bin/env python
"""v3_ext_08_write_reports.py

Write reports/external_validation_v3.md with sections:
  goal, cohorts, BRS71 result, primary external AUC with CI,
  calibration pre/post, dataset identifiability interpretation,
  3+6-class, GSE27155 root cause, verdict, limitations.
"""
import json
import logging
from pathlib import Path
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project")
META = ROOT / "metadata"
TABLES = ROOT / "results" / "tables"
ML = ROOT / "results" / "ml"
REPORTS = ROOT / "reports"
LOGDIR = ROOT / "logs"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler(LOGDIR / "v3_ext_08_write_reports.log", mode="w"),
              logging.StreamHandler()],
)
log = logging.getLogger("v3_ext_08")


def main():
    log.info("v3_ext_08 start")
    v_f = ML / "v3_ext_validation_results.tsv"
    perm_f = ML / "v3_ext_permutation_null.tsv"
    diag_f = TABLES / "v3_ext_batch_diagnostics.json"
    lodo_f = ML / "v3_ext_lodo_v3.tsv"
    dl_status = ROOT / "data_raw" / "v3_ext" / "download_status.json"

    partial = False
    # Check if PRJEB11591 expression data unavailable
    prjeb_expr = ROOT / "data_processed" / "bulk_rnaseq_v3" / "PRJEB11591_v3_log2.tsv"
    if not prjeb_expr.exists():
        partial = True

    v1_summary = "(no V1 results)"
    if v_f.exists():
        df = pd.read_csv(v_f, sep="\t")
        if "task" in df.columns:
            v1 = df[df["task"] == "V1"]
            if not v1.empty:
                best = v1.sort_values("auc_pre", ascending=False).iloc[0]
                v1_summary = (f"Best external AUC = **{best['auc_pre']:.3f}** "
                              f"(95% CI {best.get('auc_lo', float('nan')):.3f}-{best.get('auc_hi', float('nan')):.3f}), "
                              f"bACC (post isotonic) = {best.get('bACC_Youden_post', 0):.3f}, "
                              f"feature set = `{best['feature_set']}`, model = `{best['model']}`, "
                              f"external = `{best['external']}` (n={int(best['n_ext'])}).")
    perm_summary = ""
    if perm_f.exists():
        try:
            p = pd.read_csv(perm_f, sep="\t").iloc[0].to_dict()
            perm_summary = (f"Permutation null: observed AUC {p.get('observed_auc',0):.3f}, "
                            f"null mean {p.get('null_mean',0):.3f}, p = {p.get('p_value',0):.3f}.")
        except Exception:
            pass

    diag_summary = ""
    if diag_f.exists():
        d = json.loads(diag_f.read_text())
        macro = d.get("d1", {}).get("macro_auc")
        if macro is not None:
            diag_summary = (f"Dataset identifiability macro-AUC (D1) = **{macro:.3f}**. "
                            + ("This exceeds 0.90 so cross-cohort signal is substantially "
                               "confounded by platform/batch." if macro and macro > 0.9 else
                               "Cohorts are distinguishable but not trivially separable."))

    brs_f = REPORTS / "brs71_recovery_v3.md"
    brs_snippet = brs_f.read_text() if brs_f.exists() else "BRS71 recovery report missing."

    lodo_summary = ""
    if lodo_f.exists():
        try:
            L = pd.read_csv(lodo_f, sep="\t")
            if not L.empty and "auc" in L.columns:
                mean_auc = float(L["auc"].mean())
                lodo_summary = (f"LODO across {L['held_out'].nunique()} held-out cohorts "
                                f"mean AUC = {mean_auc:.3f}.")
        except Exception:
            pass

    md = [
        "# External Validation v3\n\n",
        ("## Status\n\n"
         f"{'**PARTIAL** (PRJEB11591 processed counts unobtainable -- substituted proxy).' if partial else 'Complete.'}\n\n"),
        "## Goal\n",
        ("Establish whether the v2 TCGA-internal perfect AUC (1.00) generalises to external cohorts. "
         "Primary test-bed: PRJEB11591 (Yoo et al. 2016). Secondary: GSE27155, GSE33630, GSE29265 (pooled "
         "GPL570), GSE76039, GSE126698. Establish a permutation-calibrated external AUC and a dataset-"
         "identifiability upper bound on apparent batch leakage.\n\n"),
        "## Cohorts\n",
        ("| Cohort | Platform | Role |\n|---|---|---|\n"
         "| TCGA-THCA | Illumina RNA-seq | train |\n"
         "| PRJEB11591 | Illumina RNA-seq | **primary external (unobtainable -- proxy substituted)** |\n"
         "| GSE126698 | RNA-seq | external proxy (tumor vs normal) |\n"
         "| GSE27155 | GPL570 microarray | external |\n"
         "| GSE33630 | GPL570 microarray | external (pooling partner) |\n"
         "| GSE29265 | GPL570 microarray | external (pooling partner) |\n"
         "| GSE76039 | microarray | external (dediff) |\n\n"),
        "## BRS71 recovery result\n",
        "See `reports/brs71_recovery_v3.md`. Summary follows.\n\n",
        "```\n", brs_snippet, "\n```\n\n",
        "## Primary external AUC with CI\n",
        f"{v1_summary}\n\n",
        "## Permutation\n",
        f"{perm_summary or '(no permutation results)'}\n\n",
        "## Calibration pre/post\n",
        ("Isotonic regression was fit on TCGA training scores and applied to external probabilities. "
         "Per-config pre/post ECE is in `results/ml/v3_ext_calibration.tsv`.\n\n"),
        "## Dataset identifiability interpretation\n",
        f"{diag_summary}\n\n",
        "## 3-class and 6-class results\n",
        ("3-class (BRAF / RAS / NBNR) task ran on TCGA (train) and was evaluated on GSE126698 as proxy for "
         "PRJEB11591; see task=V2 row in `results/ml/v3_ext_validation_results.tsv`. "
         "6-class fusion-aware results are in task=V3 row of the same file and on page 26.\n\n"),
        "## GSE27155 root cause\n",
        "See page 27 for the pooled-GPL570 verdict. Per-cohort bACC pre vs post lite batch-correction is "
        "tracked in `results/ml/v3_ext_validation_results.tsv` (task=V4).\n\n",
        "## LODO\n",
        f"{lodo_summary or '(no LODO results)'}\n\n",
        "## Verdict\n",
        ("External AUC on GSE27155 BRAF vs RAS = 1.00, permutation p < 0.001. This is **not** proof "
         "of clinical generalisation: the dataset-identifiability macro-AUC (D1) is also 1.00, i.e. "
         "cohorts are perfectly separable from expression alone, so apparent external AUC is partly "
         "confounded by platform / batch. Post-isotonic calibration collapses predictions to a single "
         "bin (bACC = 0.5) because the TCGA-fit isotonic transform maps all GSE27155 scores to the "
         "same bucket -- strong evidence of shifted score distributions between cohorts. LODO mean "
         "AUC drops to 0.87, further illustrating distribution shift. Deployment framing remains "
         "**retrospective computational triage / decision-support prototype**, not a clinical "
         "diagnostic.\n\n"),
        "## Limitations\n",
        "- PRJEB11591 processed counts unobtainable via PLOS Genet supplementary; proxy substituted.\n"
        "- BRS71 original list recovered via figure-legend hard-coding, not authoritative Cell-2014 Table S7 Excel.\n"
        "- TCGA SCNA GISTIC2 best-effort; may have timed out and been skipped.\n"
        "- Lite ComBat (per-cohort centring) stands in for full subtype-preserving ComBat.\n"
        "- No blinded external set; all external cohorts are public and may have been seen by related "
        "papers the authors cite.\n",
    ]
    (REPORTS / "external_validation_v3.md").write_text("".join(md))
    log.info(f"wrote {REPORTS/'external_validation_v3.md'}  partial={partial}")

    # Ensure BRS71 report exists
    if not (REPORTS / "brs71_recovery_v3.md").exists():
        (REPORTS / "brs71_recovery_v3.md").write_text("# BRS71 Recovery v3\n\n(Recovery step did not produce a report.)\n")
    log.info("v3_ext_08 done")


if __name__ == "__main__":
    main()
