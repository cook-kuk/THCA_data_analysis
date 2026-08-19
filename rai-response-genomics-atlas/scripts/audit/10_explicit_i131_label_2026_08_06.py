#!/usr/bin/env python3
"""AUDIT 10 — the explicit I-131 field, and the treatment x score interaction.

The full TCGA-THCA BCR Biotab clinical supplement (downloaded 2026-08-06) carries a
field the radiation file does not:

    i_131_radiation_first_tx_method  in {Thyroxine withdrawal, rhTSH,
                                         "Patient did not receive I-131 treatment"}

That gives, for the first time, an EXPLICIT non-radioiodine group rather than one
inferred from the absence of a radiation record - which is what Prof. Kang's objection
actually requires, and what makes a treatment-by-score interaction estimable at all.

It also exposes a conflict inside TCGA that has to be reported, not smoothed over:
some patients are explicitly recorded as not having received I-131 while carrying a
millicurie-dosed course in the radiation file.

Analyses
  1  concordance of the two sources, and the size of the conflict
  2  cohort re-definition under three rules (permissive / strict / explicit-only)
  3  the RAI-anchored recurrence result under each rule
  4  TREATMENT x SCORE INTERACTION using the explicit labels
  5  TSH preparation method (withdrawal vs rhTSH) as an efficacy covariate

Outputs
  results/tables/audit10_explicit_i131_concordance_2026_08_06.tsv
  results/tables/audit10_explicit_i131_results_2026_08_06.tsv
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT.parent / "external_data" / "tcga_thca_biotab_full" / "source"
RADIO = Path("/data/rai_atlas/raw/TCGA_THCA_biotab"
             "/nationwidechildrens.org_clinical_radiation_thca.txt")
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
OUT = ROOT / "results" / "tables"
SEED = 20260806

spec = importlib.util.spec_from_file_location(
    "base", ROOT / "scripts" / "tcga_rai_best_response_2026_08_06.py")
_base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_base)

a9spec = importlib.util.spec_from_file_location(
    "a9", ROOT / "scripts" / "audit" / "09_rai_anchored_robustness_2026_08_06.py")
_a9 = importlib.util.module_from_spec(a9spec)
a9spec.loader.exec_module(_a9)
cox_fit = _a9.cox_fit

NO_I131 = "Patient did not receive I-131 treatment"
PREP = ["Thyroxine withdrawal", "rhTSH"]


def build():
    panel = pd.read_csv(_base.PANEL, sep="\t")
    panel["patient"] = panel["sampleId"].str.slice(0, 12)
    panel = panel.drop_duplicates("patient")
    keep = [c for c in ["patient", "RAI_8", "AJCC_PATHOLOGIC_TUMOR_STAGE"]
            if c in panel.columns]
    panel = panel[keep]

    pt = pd.read_csv(EXT / "nationwidechildrens.org_clinical_patient_thca.txt",
                     sep="\t", skiprows=[1, 2], low_memory=False)
    pt = pt.rename(columns={"bcr_patient_barcode": "patient"})
    pt = pt[["patient", "i_131_radiation_first_tx_method",
             "i_131_radiation_tx_cumulative_dose", "radiation_treatment_adjuvant",
             "age_at_diagnosis", "ajcc_pathologic_tumor_stage"]]

    rad = pd.read_csv(RADIO, sep="\t", skiprows=[1, 2], low_memory=False)
    u = rad["radiation_adjuvant_units"].astype(str).str.lower()
    rai = rad[u.eq("mci")].copy()
    rai["start"] = pd.to_numeric(rai["radiation_therapy_started_days_to"], errors="coerce")
    first = (rai.groupby("bcr_patient_barcode")["start"].min().reset_index()
             .rename(columns={"bcr_patient_barcode": "patient", "start": "rai_start_day"}))

    sv = pd.read_csv(SURV, sep="\t", low_memory=False)
    sv["patient"] = sv["_PATIENT"]
    sv = sv[["patient", "new_tumor_event_dx_days_to", "last_contact_days_to",
             "PFI", "PFI.time"]].drop_duplicates("patient")
    for c in sv.columns[1:]:
        sv[c] = pd.to_numeric(sv[c], errors="coerce")

    d = (panel.merge(pt, on="patient", how="left")
              .merge(first, on="patient", how="left")
              .merge(sv, on="patient", how="left"))
    d["RAI_8"] = pd.to_numeric(d["RAI_8"], errors="coerce")
    d["age"] = pd.to_numeric(d["age_at_diagnosis"], errors="coerce")
    d["stage_high"] = d["ajcc_pathologic_tumor_stage"].astype(str).str.contains(
        "III|IV", regex=True).astype(int)
    meth = d["i_131_radiation_first_tx_method"].astype(str)
    d["explicit_no"] = meth.eq(NO_I131)
    d["explicit_yes"] = meth.isin(PREP)
    d["prep"] = np.where(meth.eq("rhTSH"), 1, np.where(meth.eq("Thyroxine withdrawal"), 0, np.nan))
    d["radfile_mci"] = d["rai_start_day"].notna()
    return d


def anchored(d):
    """RAI-anchored clock for patients with a radiation start day."""
    x = d.copy()
    has = x["new_tumor_event_dx_days_to"].notna()
    pre = has & (x["new_tumor_event_dx_days_to"] < x["rai_start_day"])
    x["time"] = np.where(has & ~pre,
                         x["new_tumor_event_dx_days_to"] - x["rai_start_day"],
                         x["last_contact_days_to"] - x["rai_start_day"])
    x["event"] = np.where(has & ~pre, 1, 0)
    return x[x["time"].notna() & (x["time"] > 0) & x["RAI_8"].notna()]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = build()
    conc, res = [], []

    def C(k, v):
        conc.append({"quantity": k, "value": v})

    # ------------------------------------------------------------- 1 concordance
    ct = pd.crosstab(
        np.where(d["explicit_no"], "explicit NO I-131",
                 np.where(d["explicit_yes"], "explicit YES I-131", "not recorded")),
        d["radfile_mci"])
    print("=== concordance: explicit field vs radiation-file mCi course ===")
    print(ct.to_string(), "\n")
    C("panel patients", len(d))
    C("explicit YES I-131 (withdrawal or rhTSH)", int(d["explicit_yes"].sum()))
    C("explicit NO I-131", int(d["explicit_no"].sum()))
    C("I-131 status not recorded", int((~d["explicit_yes"] & ~d["explicit_no"]).sum()))
    C("has an mCi course in the radiation file", int(d["radfile_mci"].sum()))
    n_conf = int((d["explicit_no"] & d["radfile_mci"]).sum())
    n_miss = int((d["explicit_yes"] & ~d["radfile_mci"]).sum())
    C("CONFLICT - explicit NO but mCi course present", n_conf)
    C("MISSED - explicit YES but no mCi course", n_miss)
    C("conflict share of the mCi cohort (%)",
      round(100 * n_conf / max(int(d["radfile_mci"].sum()), 1), 1))

    # ------------------------------------------------------------- 2/3 cohort rules
    rules = {
        "permissive (AUDIT 08: any mCi course)": d["radfile_mci"],
        "strict (mCi course AND not explicitly denied)": d["radfile_mci"] & ~d["explicit_no"],
        "explicit-only (explicit YES I-131, mCi course present)": d["explicit_yes"] & d["radfile_mci"],
    }
    for lab, mask in rules.items():
        a = anchored(d[mask])
        if len(a) < 20 or a["event"].sum() < 5:
            continue
        b, se, pv = cox_fit(a["time"].values, a["event"].values,
                            a["RAI_8"].values.reshape(-1, 1))
        res.append(dict(analysis=f"cohort rule · {lab}", n=len(a),
                        events=int(a["event"].sum()),
                        hr=float(np.exp(b[0])),
                        ci_low=float(np.exp(b[0] - 1.96 * se[0])),
                        ci_high=float(np.exp(b[0] + 1.96 * se[0])),
                        p=float(pv[0]),
                        note="Cox, panel score continuous, RAI-anchored clock"))
        print(f"cohort {lab:52} n={len(a):3d} ev={int(a['event'].sum()):2d} "
              f"HR={np.exp(b[0]):.3f} P={pv[0]:.4f}")

    # ------------------------------------------------- 4 treatment x score interaction
    print("\n=== treatment x score interaction, using the EXPLICIT labels ===")
    it = d[(d["explicit_yes"] | d["explicit_no"]) & d["RAI_8"].notna()
           & d["PFI"].notna() & d["PFI.time"].notna()].copy()
    it["trt"] = it["explicit_yes"].astype(int)
    it = it[it["PFI.time"] > 0]
    for g, lab in [(1, "explicit I-131"), (0, "explicit NO I-131")]:
        s = it[it["trt"] == g]
        print(f"  {lab:20} n={len(s):3d}  PFI events={int(s['PFI'].sum()):2d}")
        C(f"{lab} — n with score and PFI", len(s))
        C(f"{lab} — PFI events", int(s["PFI"].sum()))

    ev_no = int(it.loc[it["trt"] == 0, "PFI"].sum())
    ev_yes = int(it.loc[it["trt"] == 1, "PFI"].sum())
    if ev_no >= 5 and ev_yes >= 5:
        X = np.column_stack([it["RAI_8"].values, it["trt"].values,
                             it["RAI_8"].values * it["trt"].values])
        b, se, pv = cox_fit(it["PFI.time"].values, it["PFI"].values, X)
        res.append(dict(analysis="INTERACTION · score x explicit I-131 treatment (Cox, PFI)",
                        n=len(it), events=int(it["PFI"].sum()),
                        hr=float(np.exp(b[2])),
                        ci_low=float(np.exp(b[2] - 1.96 * se[2])),
                        ci_high=float(np.exp(b[2] + 1.96 * se[2])),
                        p=float(pv[2]),
                        note=(f"interaction term; events {ev_yes} treated / {ev_no} untreated; "
                              f"main score HR={np.exp(b[0]):.3f} P={pv[0]:.4f}; "
                              f"treatment HR={np.exp(b[1]):.3f} P={pv[1]:.4f}; "
                              "ESTIMABLE but low-event — treat as exploratory")))
        print(f"  INTERACTION HR={np.exp(b[2]):.3f} "
              f"[{np.exp(b[2]-1.96*se[2]):.3f},{np.exp(b[2]+1.96*se[2]):.3f}] P={pv[2]:.4f}")
        # stratum-specific
        for g, lab in [(1, "within explicit I-131"), (0, "within explicit NO I-131")]:
            s = it[it["trt"] == g]
            if s["PFI"].sum() < 3:
                continue
            bb, ss, pp = cox_fit(s["PFI.time"].values, s["PFI"].values,
                                 s["RAI_8"].values.reshape(-1, 1))
            res.append(dict(analysis=f"stratum · {lab}", n=len(s),
                            events=int(s["PFI"].sum()), hr=float(np.exp(bb[0])),
                            ci_low=float(np.exp(bb[0] - 1.96 * ss[0])),
                            ci_high=float(np.exp(bb[0] + 1.96 * ss[0])),
                            p=float(pp[0]), note="Cox, PFI, diagnosis-anchored"))
            print(f"  {lab:26} n={len(s):3d} ev={int(s['PFI'].sum()):2d} "
                  f"HR={np.exp(bb[0]):.3f} P={pp[0]:.4f}")
    else:
        res.append(dict(analysis="INTERACTION · score x explicit I-131 treatment",
                        n=len(it), events=int(it["PFI"].sum()), hr=np.nan,
                        ci_low=np.nan, ci_high=np.nan, p=np.nan,
                        note=(f"NOT ESTIMABLE — events {ev_yes} treated / {ev_no} untreated; "
                              "fewer than 5 events in a stratum")))
        print(f"  NOT ESTIMABLE — {ev_yes} treated / {ev_no} untreated events")

    # ------------------------------------------------------ 5 TSH preparation method
    print("\n=== TSH preparation method as an efficacy covariate ===")
    pr = anchored(d[d["radfile_mci"] & d["prep"].notna()])
    if len(pr) > 20 and pr["event"].sum() >= 5:
        for g, lab in [(0, "thyroxine withdrawal"), (1, "rhTSH")]:
            s = pr[pr["prep"] == g]
            print(f"  {lab:22} n={len(s):3d} events={int(s['event'].sum()):2d}")
        b, se, pv = cox_fit(pr["time"].values, pr["event"].values,
                            pr[["RAI_8", "prep"]].values)
        res.append(dict(analysis="Cox · score + TSH preparation method (RAI-anchored)",
                        n=len(pr), events=int(pr["event"].sum()),
                        hr=float(np.exp(b[0])),
                        ci_low=float(np.exp(b[0] - 1.96 * se[0])),
                        ci_high=float(np.exp(b[0] + 1.96 * se[0])),
                        p=float(pv[0]),
                        note=(f"rhTSH vs withdrawal HR={np.exp(b[1]):.3f} P={pv[1]:.4f}; "
                              "preparation method affects iodine uptake and is never modelled")))
        print(f"  score HR={np.exp(b[0]):.3f} P={pv[0]:.4f} | "
              f"rhTSH vs withdrawal HR={np.exp(b[1]):.3f} P={pv[1]:.4f}")
        rho, prho = stats.mannwhitneyu(pr.loc[pr["prep"] == 1, "RAI_8"],
                                       pr.loc[pr["prep"] == 0, "RAI_8"])
        C("panel score differs by preparation method (Mann-Whitney P)", round(float(prho), 4))

    pd.DataFrame(conc).to_csv(
        OUT / "audit10_explicit_i131_concordance_2026_08_06.tsv", sep="\t", index=False)
    pd.DataFrame(res).to_csv(
        OUT / "audit10_explicit_i131_results_2026_08_06.tsv", sep="\t", index=False)
    print(f"\nwrote 2 tables to {OUT}")


if __name__ == "__main__":
    main()
