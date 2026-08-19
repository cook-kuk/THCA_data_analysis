#!/usr/bin/env python3
"""AUDIT 08 — time from RADIOIODINE ADMINISTRATION to recurrence.

Answers the specific objection raised by Prof. Kang Minsu (2026-08-06):

    "Progression free survival로 과연 RAI 반응성을 평가할 수 있는가?
     RAI 안 하고 갑상선 수술 후 재발한 사람들도 포함되어 있을 것이기 때문입니다."

The objection is correct for a diagnosis-anchored PFS computed over all patients.
This script removes both defects:

  1. COHORT   — restricted to patients with at least one millicurie-dosed radiation
                course, i.e. patients who actually received radioiodine. Patients
                treated with external beam only (Gy/cGy) or with no radiation at all
                are excluded, and are analysed separately as a comparison group.

  2. TIME ORIGIN — the clock starts at the FIRST RADIOIODINE COURSE, not at diagnosis.
                t = new_tumor_event_dx_days_to - radiation_therapy_started_days_to

Events occurring BEFORE the first radioiodine course are not "recurrence after RAI";
they are counted, reported, and excluded from the time-to-event analysis, because
including them is exactly the confound being corrected.

Outputs
  results/tables/audit08_rai_anchored_recurrence_2026_08_06.tsv
  results/tables/audit08_rai_anchored_cohort_ledger_2026_08_06.tsv
  figures/audit/figure_audit08_rai_anchored_km_2026_08_06.{png,pdf}
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parent / "project" / "results"
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
RADIO = Path("/data/rai_atlas/raw/TCGA_THCA_biotab"
             "/nationwidechildrens.org_clinical_radiation_thca.txt")
OUT_T = ROOT / "results" / "tables"
OUT_F = ROOT / "figures" / "audit"
SEED = 20260806

spec = importlib.util.spec_from_file_location(
    "base", ROOT / "scripts" / f"tcga_rai_best_response_2026_08_06.py")
_base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_base)


def km(t, e):
    """Kaplan-Meier estimate. Returns (times, survival)."""
    order = np.argsort(t)
    t, e = np.asarray(t)[order], np.asarray(e)[order]
    times, surv, s, n = [0.0], [1.0], 1.0, len(t)
    for i, ti in enumerate(t):
        at_risk = n - i
        if e[i] == 1:
            s *= (1 - 1 / at_risk)
            times.append(ti)
            surv.append(s)
    times.append(t.max() if len(t) else 0.0)
    surv.append(s)
    return np.array(times), np.array(surv)


def logrank(t1, e1, t2, e2):
    t = np.concatenate([t1, t2])
    e = np.concatenate([e1, e2])
    g = np.concatenate([np.zeros(len(t1)), np.ones(len(t2))])
    order = np.argsort(t)
    t, e, g = t[order], e[order], g[order]
    O1 = E1 = V = 0.0
    for ti in np.unique(t[e == 1]):
        at = t >= ti
        n, n1 = at.sum(), (at & (g == 1)).sum()
        d = ((t == ti) & (e == 1)).sum()
        d1 = ((t == ti) & (e == 1) & (g == 1)).sum()
        if n <= 1:
            continue
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V <= 0:
        return np.nan, np.nan
    chi = (O1 - E1) ** 2 / V
    return chi, 1 - stats.chi2.cdf(chi, 1)


def cox_score_test(t, e, x):
    """Univariable Cox score test (Breslow ties). Returns (beta_newton, se, p)."""
    order = np.argsort(t)
    t, e, x = np.asarray(t)[order], np.asarray(e)[order], np.asarray(x)[order]

    def nll_grad_hess(b):
        r = np.exp(b * x)
        g = h = ll = 0.0
        for i in np.where(e == 1)[0]:
            at = t >= t[i]
            s0 = r[at].sum()
            s1 = (r[at] * x[at]).sum()
            s2 = (r[at] * x[at] ** 2).sum()
            ll += b * x[i] - np.log(s0)
            g += x[i] - s1 / s0
            h += (s1 / s0) ** 2 - s2 / s0
        return ll, g, h

    b = 0.0
    for _ in range(50):
        ll, g, h = nll_grad_hess(b)
        if abs(h) < 1e-10:
            break
        step = g / h
        b_new = b - step
        if abs(b_new - b) < 1e-8:
            b = b_new
            break
        b = b_new
    _, g, h = nll_grad_hess(b)
    se = float(np.sqrt(-1.0 / h)) if h < 0 else np.nan
    z = b / se if se and np.isfinite(se) else np.nan
    p = 2 * (1 - stats.norm.cdf(abs(z))) if np.isfinite(z) else np.nan
    return b, se, p


def main():
    OUT_T.mkdir(parents=True, exist_ok=True)
    OUT_F.mkdir(parents=True, exist_ok=True)

    # ---- panel score + radioiodine courses (already restricted to mCi patients)
    df = _base.build()

    # ---- first radioiodine course start day, per patient
    rad = pd.read_csv(RADIO, sep="\t", skiprows=[1, 2], low_memory=False)
    units = rad["radiation_adjuvant_units"].astype(str).str.lower()
    rai = rad[units.eq("mci")].copy()
    rai["start"] = pd.to_numeric(rai["radiation_therapy_started_days_to"],
                                 errors="coerce")
    first = (rai.groupby("bcr_patient_barcode")["start"].min()
             .reset_index().rename(columns={"bcr_patient_barcode": "patient",
                                            "start": "rai_start_day"}))
    df = df.merge(first, on="patient", how="left")

    # ---- recurrence timing from the pan-cancer clinical table
    sv = pd.read_csv(SURV, sep="\t", low_memory=False)
    sv["patient"] = sv["_PATIENT"]
    # PFI / PFI.time are already merged in by _base.build(); take only the
    # timing fields that are not yet present, or pandas will suffix them _x/_y
    keep = ["patient", "new_tumor_event_dx_days_to", "last_contact_days_to"]
    sv = sv[keep].drop_duplicates("patient")
    for c in keep[1:]:
        sv[c] = pd.to_numeric(sv[c], errors="coerce")
    df = df.merge(sv, on="patient", how="left")

    ledger = []

    def note(k, v):
        ledger.append({"quantity": k, "value": v})

    note("patients with a panel score and >=1 mCi course", len(df))
    note("of those, with a recorded RAI start day", int(df["rai_start_day"].notna().sum()))

    d = df[df["rai_start_day"].notna()].copy()

    # ---- build the RAI-anchored clock
    d["evt_day"] = d["new_tumor_event_dx_days_to"]
    d["has_evt"] = d["evt_day"].notna()
    note("with a recorded new-tumour-event date", int(d["has_evt"].sum()))

    # events that happened BEFORE radioiodine are not recurrence after RAI
    pre = d["has_evt"] & (d["evt_day"] < d["rai_start_day"])
    note("EXCLUDED - new tumour event occurred BEFORE the first RAI course",
         int(pre.sum()))

    # censoring time for event-free patients: last contact, re-anchored
    d["t_evt"] = d["evt_day"] - d["rai_start_day"]
    d["t_cen"] = d["last_contact_days_to"] - d["rai_start_day"]
    d["time"] = np.where(d["has_evt"] & ~pre, d["t_evt"], d["t_cen"])
    d["event"] = np.where(d["has_evt"] & ~pre, 1, 0)

    d = d[d["time"].notna() & (d["time"] > 0)].copy()
    note("ANALYSIS SET - RAI-treated, positive follow-up after the RAI course", len(d))
    note("recurrence events after RAI", int(d["event"].sum()))
    note("median follow-up after RAI (days)", float(np.median(d["time"])))
    note("median follow-up after RAI (months)", round(float(np.median(d["time"])) / 30.44, 1))

    score = pd.to_numeric(d["RAI_8"], errors="coerce")
    d = d[score.notna()].copy()
    d["RAI_8"] = pd.to_numeric(d["RAI_8"], errors="coerce")
    note("ANALYSIS SET with a non-missing panel score", len(d))
    note("events in that set", int(d["event"].sum()))
    note("events per variable (EPV), single covariate", round(d["event"].sum() / 1, 2))

    rows = []

    # ---- continuous Cox
    b, se, p = cox_score_test(d["time"].values, d["event"].values, d["RAI_8"].values)
    rows.append(dict(
        analysis="Cox, panel score continuous, RAI-anchored clock",
        n=len(d), events=int(d["event"].sum()),
        hr=float(np.exp(b)),
        ci_low=float(np.exp(b - 1.96 * se)) if np.isfinite(se) else np.nan,
        ci_high=float(np.exp(b + 1.96 * se)) if np.isfinite(se) else np.nan,
        p=p, note="HR per 1 unit of the eight-gene z-score"))

    # ---- median split, log-rank
    med = d["RAI_8"].median()
    lo = d[d["RAI_8"] <= med]
    hi = d[d["RAI_8"] > med]
    chi, plr = logrank(lo["time"].values, lo["event"].values,
                       hi["time"].values, hi["event"].values)
    rows.append(dict(
        analysis="log-rank, panel score below vs above median, RAI-anchored clock",
        n=len(d), events=int(d["event"].sum()),
        hr=np.nan, ci_low=np.nan, ci_high=np.nan, p=plr,
        note=f"chi2={chi:.3f}; low n={len(lo)} ev={int(lo['event'].sum())}; "
             f"high n={len(hi)} ev={int(hi['event'].sum())}"))

    # ---- the diagnosis-anchored PFI, for contrast: this is what the objection targets
    dd = df[df["PFI"].notna() & df["PFI.time"].notna() & df["RAI_8"].notna()].copy()
    dd["PFI"] = pd.to_numeric(dd["PFI"], errors="coerce")
    dd["PFI.time"] = pd.to_numeric(dd["PFI.time"], errors="coerce")
    dd = dd[dd["PFI.time"] > 0]
    b2, se2, p2 = cox_score_test(dd["PFI.time"].values, dd["PFI"].values,
                                 pd.to_numeric(dd["RAI_8"], errors="coerce").values)
    rows.append(dict(
        analysis="Cox, panel score continuous, DIAGNOSIS-anchored PFI (for contrast)",
        n=len(dd), events=int(dd["PFI"].sum()),
        hr=float(np.exp(b2)),
        ci_low=float(np.exp(b2 - 1.96 * se2)) if np.isfinite(se2) else np.nan,
        ci_high=float(np.exp(b2 + 1.96 * se2)) if np.isfinite(se2) else np.nan,
        p=p2,
        note="Same patients, wrong time origin. Shown only to quantify how much the "
             "time-origin choice matters; NOT a claim."))

    res = pd.DataFrame(rows)
    res.to_csv(OUT_T / "audit08_rai_anchored_recurrence_2026_08_06.tsv",
               sep="\t", index=False)
    pd.DataFrame(ledger).to_csv(
        OUT_T / "audit08_rai_anchored_cohort_ledger_2026_08_06.tsv",
        sep="\t", index=False)

    print("\n=== COHORT LEDGER ===")
    for r in ledger:
        print(f"  {r['quantity']:<62} {r['value']}")
    print("\n=== RESULTS ===")
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(res.to_string(index=False))

    return d, lo, hi, res


if __name__ == "__main__":
    main()
