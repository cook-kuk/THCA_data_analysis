#!/usr/bin/env python3
"""AUDIT 01 — TCGA RAI cohort construction and index-course sensitivity.

Why this exists
---------------
The 2026-08-06 first-pass analysis aggregated `treatment_best_response` to one value per
patient by taking the BEST response across that patient's radioiodine courses. Auditing the
raw table shows that choice is not innocuous:

  TCGA-FK-A3S3  day 162, 106.4 mCi -> Stable Disease
                day 581, 197.8 mCi -> Complete Response
  TCGA-FE-A3PA  day  91, primary   -> Partial Response
                day 378, distant   -> Complete Response
  TCGA-H2-A3RI  day  84, primary   -> Radiographic Progressive Disease

Because the value varies WITHIN a patient and tracks course timing and treatment site, the
field is course-specific, not a patient-level best response. That resolves the provenance
question empirically (a patient-level field would be constant across a patient's rows), and
it also means "best across courses" biases the cohort toward complete response.

This script rebuilds the cohort under five explicit index-course rules and reports the
initial-response test under each, so the reader can see how much the conclusion depends on
a decision that was previously implicit.

Outputs
  audit/rai_integration_20260806/03_TCGA_COHORT_CONSTRUCTION_AUDIT.md
  results/tables/audit01_cohort_ledger_2026_08_06.tsv
  results/tables/audit01_index_course_sensitivity_2026_08_06.tsv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
REPO = Path("/data/thca/repo_results")
PANEL = REPO / "r17_tcga_panel_d4p2_reconciliation" / "r17_per_sample_merged.tsv"
RADIO = Path("/data/rai_atlas/raw/TCGA_THCA_biotab/"
             "nationwidechildrens.org_clinical_radiation_thca.txt")
TAB = ROOT / "results" / "tables"
AUDIT = ROOT.parent / "audit" / "rai_integration_20260806"
TAB.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)

STAMP = "2026_08_06"
SEED = 20260806
RESP = ["Complete Response", "Partial Response", "Stable Disease",
        "Radiographic Progressive Disease"]
RANK = {r: i for i, r in enumerate(RESP)}


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                 / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp if sp > 0 else np.nan


def boot_ci(a, b, n=5000, seed=SEED):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    v = [cohens_d(rng.choice(a, len(a), True), rng.choice(b, len(b), True)) for _ in range(n)]
    v = np.asarray([x for x in v if np.isfinite(x)])
    return (np.percentile(v, 2.5), np.percentile(v, 97.5)) if len(v) else (np.nan, np.nan)


def main():
    rad = pd.read_csv(RADIO, sep="\t", skiprows=[1, 2], low_memory=False)
    rad["unit"] = rad["radiation_adjuvant_units"].astype(str).str.strip().str.lower()
    rad["dose"] = pd.to_numeric(rad["radiation_total_dose"], errors="coerce")
    rad["start"] = pd.to_numeric(rad["radiation_therapy_started_days_to"], errors="coerce")
    rai = rad[rad["unit"] == "mci"].copy()
    ebrt = set(rad.loc[rad["unit"].isin(["gy", "cgy"]), "bcr_patient_barcode"])

    panel = pd.read_csv(PANEL, sep="\t")
    panel["patient"] = panel["sampleId"].str.slice(0, 12)
    panel = panel.drop_duplicates("patient")

    # ---------------- ledger ----------------
    led = [
        ("radiation course rows (all units)", len(rad)),
        ("unique patients with any radiation record", rad.bcr_patient_barcode.nunique()),
        ("mCi courses (= radioiodine)", len(rai)),
        ("Gy / cGy courses (= external beam)", int(rad["unit"].isin(["gy", "cgy"]).sum())),
        ("courses with unit not recorded", int((rad["unit"] == "[not available]").sum())),
        ("unique patients with >=1 mCi course", rai.bcr_patient_barcode.nunique()),
        ("patients with >1 mCi course", int((rai.groupby('bcr_patient_barcode').size() > 1).sum())),
        ("max mCi courses per patient", int(rai.groupby('bcr_patient_barcode').size().max())),
        ("mCi patients who also received Gy/cGy", len(set(rai.bcr_patient_barcode) & ebrt)),
        ("mCi courses with an evaluable response", int(rai.treatment_best_response.isin(RESP).sum())),
        ("mCi courses with start day recorded", int(rai["start"].notna().sum())),
        ("patients whose courses DISAGREE on response",
         int((rai[rai.treatment_best_response.isin(RESP)]
              .groupby('bcr_patient_barcode').treatment_best_response.nunique() > 1).sum())),
        ("patients in R17 panel table", len(panel)),
        ("mCi patients matched to a panel score",
         len(set(rai.bcr_patient_barcode) & set(panel.patient))),
    ]
    ledger = pd.DataFrame(led, columns=["quantity", "value"])
    ledger.to_csv(TAB / f"audit01_cohort_ledger_{STAMP}.tsv", sep="\t", index=False)
    print(ledger.to_string(index=False))

    # ---------------- index-course rules ----------------
    ev = rai[rai.treatment_best_response.isin(RESP)].copy()
    ev["rank"] = ev.treatment_best_response.map(RANK)

    rules = {}
    rules["first course (earliest start day)"] = (
        ev.sort_values("start").groupby("bcr_patient_barcode").first()["rank"])
    rules["last course"] = (
        ev.sort_values("start").groupby("bcr_patient_barcode").last()["rank"])
    rules["highest-dose course"] = (
        ev.sort_values("dose").groupby("bcr_patient_barcode").last()["rank"])
    rules["best response across courses"] = ev.groupby("bcr_patient_barcode")["rank"].min()
    rules["worst response across courses"] = ev.groupby("bcr_patient_barcode")["rank"].max()

    rows = []
    for name, s in rules.items():
        df = pd.DataFrame({"rank": s})
        df["patient"] = df.index
        df = df.merge(panel[["patient", "RAI_8"]], on="patient", how="inner").dropna()
        df["nonCR"] = (df["rank"] > 0).astype(int)
        a = df.loc[df.nonCR == 1, "RAI_8"].values
        b = df.loc[df.nonCR == 0, "RAI_8"].values
        if len(a) < 2 or len(b) < 2:
            continue
        d = cohens_d(a, b)
        lo, hi = boot_ci(a, b)
        p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        rho, prho = stats.spearmanr(df["RAI_8"], df["rank"])
        rows.append(dict(index_rule=name, n_total=len(df), n_nonCR=len(a), n_CR=len(b),
                         cohens_d=round(float(d), 3), ci_low=round(float(lo), 3),
                         ci_high=round(float(hi), 3), p_mwu=float(p),
                         spearman_rho_ordinal=round(float(rho), 4), p_ordinal=float(prho)))
        print(f"\n{name}: n={len(df)} (non-CR {len(a)} / CR {len(b)})  "
              f"d={d:+.3f} [{lo:+.3f},{hi:+.3f}]  P={p:.3g}  rho={rho:+.3f}")

    sens = pd.DataFrame(rows)
    sens.to_csv(TAB / f"audit01_index_course_sensitivity_{STAMP}.tsv", sep="\t", index=False)

    # ---------------- equivalence framing ----------------
    # A null is only informative against a pre-specified margin. None was pre-specified in the
    # first pass, so report the smallest effect the CI excludes rather than claiming absence.
    prim = sens[sens.index_rule.str.startswith("first course")]
    if len(prim):
        r = prim.iloc[0]
        print("\n--- equivalence framing (first-course rule) ---")
        print(f"observed d = {r.cohens_d}; 95% CI [{r.ci_low}, {r.ci_high}]")
        print(f"the data exclude |d| > {max(abs(r.ci_low), abs(r.ci_high)):.2f} but no "
              f"equivalence margin was pre-specified, so this is NOT evidence of no effect")

    with open(AUDIT / "03_TCGA_COHORT_CONSTRUCTION_AUDIT.md", "w") as fh:
        fh.write("# AUDIT 03 — TCGA radioiodine cohort construction\n\n")
        fh.write("Generated 2026-08-06 by "
                 "`rai-response-genomics-atlas/scripts/audit/01_tcga_cohort_and_index_course.py`\n\n")
        fh.write("## Provenance of `treatment_best_response`\n\n")
        fh.write("The field sits in the radiation table keyed by `bcr_radiation_barcode`, its "
                 "human label is `measure_of_response` (CDE 2857291), and every mCi course "
                 "carries a start day. Decisively, the value **varies within a patient across "
                 "courses** and tracks course timing and treatment site:\n\n")
        fh.write("| patient | course | dose (mCi) | site | response |\n|---|---|---|---|---|\n")
        fh.write("| TCGA-FK-A3S3 | day 162 | 106.4 | primary | Stable Disease |\n")
        fh.write("| TCGA-FK-A3S3 | day 581 | 197.8 | primary | Complete Response |\n")
        fh.write("| TCGA-FE-A3PA | day 91 | 217 | primary | Partial Response |\n")
        fh.write("| TCGA-FE-A3PA | day 378 | 217 | distant | Complete Response |\n\n")
        fh.write("A patient-level best-response field would be constant across a patient's "
                 "rows. It is not. The field is therefore course-specific.\n\n")
        fh.write("**Consequence for the first-pass analysis.** That analysis aggregated to "
                 "the *best* response across courses, which discards the earlier, worse "
                 "response (Stable Disease for TCGA-FK-A3S3) and biases the cohort toward "
                 "complete response. For the question 'does the score relate to the initial "
                 "response to radioiodine', the **first course** is the correct index.\n\n")
        fh.write("## Cohort ledger\n\n| quantity | value |\n|---|---|\n")
        for q, v in led:
            fh.write(f"| {q} | {v} |\n")
        fh.write("\n## Index-course sensitivity\n\n")
        fh.write(sens.to_markdown(index=False) if len(sens) else "_no rule produced a test_")
        fh.write("\n\n## Statistical framing correction\n\n")
        fh.write("The first pass wrote that the cohort 'detects d >= 0.63 at 80% power, so "
                 "clinically meaningful effects are excluded'. That conflates a "
                 "minimum-detectable-effect calculation with an equivalence test. No "
                 "equivalence margin was pre-specified. The defensible statement is that the "
                 "confidence interval excludes effects larger than its own bounds, and that "
                 "no association was detected — not that no effect exists.\n")
    print(f"\nwrote {AUDIT / '03_TCGA_COHORT_CONSTRUCTION_AUDIT.md'}")


if __name__ == "__main__":
    main()
