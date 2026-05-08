"""
Moffitt 2015 Nat Genet classical/basal classifier on TCGA-PAAD.

  - Score = mean z-score across panel genes
  - Class = sign(z_basal - z_classical), with delta magnitude as confidence
  - KM survival on OS by class
  - Cross-walk to KRAS / TP53 / SMAD4 mutation status
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

RAW = Path("/data/pdac_poc/raw")
OUT = Path("/data/pdac_poc/processed")
RES = Path("/data/pdac_poc/results")
RES.mkdir(parents=True, exist_ok=True)

MOFFITT_BASAL = ["VGLL1", "UCA1", "S100A2", "LY6D", "SPRR3", "SPRR1B", "LEMD1",
                 "LYPD3", "KRT15", "DHRS9", "AREG", "CST6", "SERPINB3",
                 "KRT6A", "SERPINB4", "FAM83A", "SCEL", "FGFBP1", "KRT7",
                 "KRT17", "GPR87", "TNS4", "SLC2A1"]
MOFFITT_CLASSICAL = ["BTNL8", "FAM3D", "AGR3", "CTSE", "LYZ", "TFF2", "TFF1",
                     "ANXA10", "LGALS4", "ECT2", "CLRN3", "MYO1A", "CLDN18",
                     "LRRC31", "TFF3", "CDX2", "SERPINA10", "VSIG2", "TSPAN8",
                     "ST6GALNAC1", "AGR2", "TOX3"]


def main():
    z = pd.read_csv(RAW / "mrna_z_Zscores.tsv", sep="\t", index_col=0)
    clin = pd.read_csv(RAW / "clinical.tsv", sep="\t")
    mut = pd.read_csv(RAW / "mutations.tsv", sep="\t")

    basal_genes = [g for g in MOFFITT_BASAL if g in z.columns]
    class_genes = [g for g in MOFFITT_CLASSICAL if g in z.columns]
    print(f"[03] basal genes available: {len(basal_genes)}/{len(MOFFITT_BASAL)}")
    print(f"[03] classical genes available: {len(class_genes)}/{len(MOFFITT_CLASSICAL)}")

    z["score_basal"] = z[basal_genes].mean(axis=1)
    z["score_classical"] = z[class_genes].mean(axis=1)
    z["delta_basal_minus_classical"] = z["score_basal"] - z["score_classical"]
    z["moffitt_call"] = np.where(z["delta_basal_minus_classical"] > 0,
                                  "basal-like", "classical")
    z["confidence"] = z["delta_basal_minus_classical"].abs()

    # Driver gene cross-walk
    driver_status = (mut.groupby(["sampleId", "hugo"])
                     .size().unstack(fill_value=0).clip(upper=1))
    for g in ["KRAS", "TP53", "SMAD4", "CDKN2A"]:
        z[f"mut_{g}"] = z.index.map(lambda s: int(driver_status.get(g, {}).get(s, 0))
                                    if g in driver_status.columns else 0)

    out_cols = ["score_basal", "score_classical", "delta_basal_minus_classical",
                "moffitt_call", "confidence", "mut_KRAS", "mut_TP53",
                "mut_SMAD4", "mut_CDKN2A"]
    out = z[out_cols].copy()
    out.index.name = "sampleId"
    out.to_csv(RES / "moffitt_calls.tsv", sep="\t")

    summary = {
        "n_total": int(len(out)),
        "n_basal": int((out["moffitt_call"] == "basal-like").sum()),
        "n_classical": int((out["moffitt_call"] == "classical").sum()),
        "fraction_basal": float((out["moffitt_call"] == "basal-like").mean()),
        "kras_in_basal": float(out.query("moffitt_call=='basal-like'")["mut_KRAS"].mean()),
        "kras_in_classical": float(out.query("moffitt_call=='classical'")["mut_KRAS"].mean()),
        "tp53_in_basal": float(out.query("moffitt_call=='basal-like'")["mut_TP53"].mean()),
        "tp53_in_classical": float(out.query("moffitt_call=='classical'")["mut_TP53"].mean()),
        "smad4_in_basal": float(out.query("moffitt_call=='basal-like'")["mut_SMAD4"].mean()),
        "smad4_in_classical": float(out.query("moffitt_call=='classical'")["mut_SMAD4"].mean()),
    }

    # Survival
    clin_idx = clin.set_index("sampleId")
    os_months = pd.to_numeric(clin_idx.get("OS_MONTHS"), errors="coerce")
    os_status = clin_idx.get("OS_STATUS")
    if os_months is not None and os_status is not None:
        ev = os_status.fillna("").str.startswith("1").astype(int)
        merged = (out.join(pd.DataFrame({"os_months": os_months, "event": ev}))
                  .dropna(subset=["os_months"]))
        merged = merged[merged["os_months"] > 0]
        b = merged.query("moffitt_call=='basal-like'")
        c = merged.query("moffitt_call=='classical'")
        if len(b) >= 5 and len(c) >= 5:
            lr = logrank_test(b["os_months"], c["os_months"],
                              b["event"], c["event"])
            kmf = KaplanMeierFitter()
            kmf.fit(b["os_months"], b["event"])
            med_b = kmf.median_survival_time_
            kmf.fit(c["os_months"], c["event"])
            med_c = kmf.median_survival_time_
            summary.update({
                "survival_n_basal": int(len(b)),
                "survival_n_classical": int(len(c)),
                "median_OS_months_basal": float(med_b) if not np.isnan(med_b) else None,
                "median_OS_months_classical": float(med_c) if not np.isnan(med_c) else None,
                "logrank_p_basal_vs_classical": float(lr.p_value),
                "logrank_test_statistic": float(lr.test_statistic),
            })
            print(f"[03] median OS basal={med_b:.1f}m  classical={med_c:.1f}m  p={lr.p_value:.4g}")

    with open(OUT / "moffitt_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[03] basal={summary['n_basal']}  classical={summary['n_classical']}  basal-fraction={summary['fraction_basal']:.2%}")
    print("[03] done")


if __name__ == "__main__":
    main()
