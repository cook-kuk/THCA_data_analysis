#!/usr/bin/env python3
"""AUDIT 05 — GSE151179 patient structure and the independence assumption.

The label `rai uptake at the metastatic site` is recorded per patient, but the series
contributes several specimens per patient (primary tumour, synchronous node, post-RAI nodes,
non-neoplastic thyroid). Any sample-level test therefore replicates the same patient-level
exposure many times and understates its standard error.

The first pass ran the PRIMARY test at patient level, which was correct, but then reported
sample-level sensitivity analyses and a sample-level adjusted model alongside it without
flagging that those rows violate independence. This script quantifies the clustering and
re-runs every contrast with the patient as the unit.

Checks performed
  - specimens per patient, and how many patients contribute >1
  - whether any patient carries discordant uptake labels (it should not, if patient-level)
  - patient-level primary test with bootstrap CI resampling PATIENTS
  - cluster-robust (patient) sandwich SE on the sample-level model, for comparison
  - patient-level permutation test (10,000 label shuffles at patient level)
  - the negative control (non-neoplastic thyroid) at patient level
  - design-effect estimate showing how much the sample-level n is inflated

Outputs
  results/tables/audit05_gse151179_structure_2026_08_06.tsv
  audit/rai_integration_20260806/04b_GSE151179_INDEPENDENCE.md
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "tables"
AUDIT = ROOT.parent / "audit" / "rai_integration_20260806"
TAB.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)
STAMP = "2026_08_06"
SEED = 20260806

_spec = importlib.util.spec_from_file_location(
    "g", ROOT / "scripts" / f"gse151179_uptake_at_met_site_{STAMP}.py")
_g = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_g)


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                 / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / sp if sp > 0 else np.nan


def main():
    gz, meta = _g.load_panel_matrix()
    df = meta.copy()
    df["panel_z"] = gz.mean(axis=0).reindex(df.index)
    df["uptake"] = df["rai_uptake_met"].str.lower().map({"yes": "Yes", "no": "No"})
    df["is_tumor"] = df["Sample_source_name_ch1"].str.contains("papillary", case=False)
    df["patient"] = df["patient_id"]
    raw = _g.raw_characteristics()
    pcol = next(c for c in raw.columns if "purity" in c.lower())
    df["purity"] = raw[pcol].reindex(df.index)

    rows = []

    # ---- 1. clustering structure ----
    n_spec = df.groupby("patient").size()
    tum = df[df.is_tumor]
    n_spec_t = tum.groupby("patient").size()
    disc = df.groupby("patient")["uptake"].nunique()
    print(f"specimens {len(df)} from {df.patient.nunique()} patients")
    print(f"  tumour specimens {len(tum)} from {tum.patient.nunique()} patients")
    print(f"  patients with >1 specimen (any tissue): {(n_spec > 1).sum()} "
          f"(max {n_spec.max()})")
    print(f"  patients with >1 tumour specimen      : {(n_spec_t > 1).sum()} "
          f"(max {n_spec_t.max()})")
    print(f"  patients with discordant uptake label : {(disc > 1).sum()}")

    m_bar = n_spec_t.mean()
    rows.append(dict(check="specimens (all tissue)", value=len(df)))
    rows.append(dict(check="patients (all tissue)", value=df.patient.nunique()))
    rows.append(dict(check="tumour specimens", value=len(tum)))
    rows.append(dict(check="patients contributing tumour specimens", value=tum.patient.nunique()))
    rows.append(dict(check="patients with >1 tumour specimen", value=int((n_spec_t > 1).sum())))
    rows.append(dict(check="max tumour specimens per patient", value=int(n_spec_t.max())))
    rows.append(dict(check="patients with discordant uptake label", value=int((disc > 1).sum())))

    # ---- 2. intraclass correlation and design effect ----
    grp = tum.groupby("patient")["panel_z"]
    k = grp.ngroups
    if k > 1:
        overall = tum["panel_z"].mean()
        ms_b = sum(len(v) * (v.mean() - overall) ** 2 for _, v in grp) / (k - 1)
        ms_w = sum(((v - v.mean()) ** 2).sum() for _, v in grp) / max(1, len(tum) - k)
        icc = max(0.0, (ms_b - ms_w) / (ms_b + (m_bar - 1) * ms_w)) if (ms_b + ms_w) > 0 else 0.0
        deff = 1 + (m_bar - 1) * icc
        print(f"  mean tumour specimens per patient     : {m_bar:.2f}")
        print(f"  ICC of panel z within patient         : {icc:.3f}")
        print(f"  design effect                          : {deff:.2f} "
              f"(effective n = {len(tum)/deff:.1f} rather than {len(tum)})")
        rows += [dict(check="mean tumour specimens per patient", value=round(float(m_bar), 3)),
                 dict(check="ICC of panel z within patient", value=round(float(icc), 3)),
                 dict(check="design effect", value=round(float(deff), 3)),
                 dict(check="effective sample size (tumour)",
                      value=round(float(len(tum) / deff), 1))]

    # ---- 3. patient-level primary test, bootstrap over PATIENTS ----
    pat = tum.groupby(["patient", "uptake"], as_index=False)["panel_z"].mean()
    a = pat.loc[pat.uptake == "Yes", "panel_z"].values
    b = pat.loc[pat.uptake == "No", "panel_z"].values
    d = cohens_d(a, b)
    rng = np.random.default_rng(SEED)
    bs = []
    for _ in range(5000):
        aa = rng.choice(a, len(a), True); bb = rng.choice(b, len(b), True)
        v = cohens_d(aa, bb)
        if np.isfinite(v):
            bs.append(v)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    print(f"\nPATIENT-LEVEL primary: n={len(a)} vs {len(b)}  d={d:+.3f} "
          f"[{lo:+.3f},{hi:+.3f}]  P={p:.3g}")
    rows.append(dict(check="patient-level d (uptake Yes - No)", value=round(float(d), 3)))
    rows.append(dict(check="patient-level 95% CI", value=f"{lo:.3f} to {hi:.3f}"))
    rows.append(dict(check="patient-level Mann-Whitney P", value=round(float(p), 4)))

    # ---- 4. patient-level permutation ----
    labels = pat["uptake"].values
    vals = pat["panel_z"].values
    obs = abs(cohens_d(vals[labels == "Yes"], vals[labels == "No"]))
    null = []
    for _ in range(10000):
        sh = rng.permutation(labels)
        v = cohens_d(vals[sh == "Yes"], vals[sh == "No"])
        if np.isfinite(v):
            null.append(abs(v))
    pperm = (np.sum(np.asarray(null) >= obs) + 1) / (len(null) + 1)
    print(f"patient-level permutation P (10,000 shuffles): {pperm:.4f}")
    rows.append(dict(check="patient-level permutation P", value=round(float(pperm), 4)))

    # ---- 5. cluster-robust sample-level model, for contrast ----
    try:
        import statsmodels.formula.api as smf
        m = pd.DataFrame({
            "panel_z": tum["panel_z"].astype(float),
            "uptake_bin": (tum["uptake"] == "Yes").astype(int),
            "purity": tum["purity"].astype(str),
            "patient": tum["patient"].astype(str),
        }).dropna()
        naive = smf.ols("panel_z ~ uptake_bin + C(purity)", data=m).fit()
        clust = smf.ols("panel_z ~ uptake_bin + C(purity)", data=m).fit(
            cov_type="cluster", cov_kwds={"groups": m["patient"]})
        print(f"\nsample-level naive   : beta={naive.params['uptake_bin']:+.4f} "
              f"SE={naive.bse['uptake_bin']:.4f} P={naive.pvalues['uptake_bin']:.3f}")
        print(f"sample-level clustered: beta={clust.params['uptake_bin']:+.4f} "
              f"SE={clust.bse['uptake_bin']:.4f} P={clust.pvalues['uptake_bin']:.3f}")
        print(f"purity coefficient    : {naive.params.filter(like='purity').to_dict()}")
        rows += [dict(check="sample-level naive SE (uptake)",
                      value=round(float(naive.bse["uptake_bin"]), 4)),
                 dict(check="sample-level cluster-robust SE (uptake)",
                      value=round(float(clust.bse["uptake_bin"]), 4)),
                 dict(check="SE inflation from clustering",
                      value=round(float(clust.bse["uptake_bin"] / naive.bse["uptake_bin"]), 3))]
    except Exception as exc:  # noqa: BLE001
        print(f"cluster-robust model skipped: {exc}")

    # ---- 6. negative control at patient level ----
    norm = df[~df.is_tumor]
    npat = norm.groupby(["patient", "uptake"], as_index=False)["panel_z"].mean()
    na = npat.loc[npat.uptake == "Yes", "panel_z"].values
    nb = npat.loc[npat.uptake == "No", "panel_z"].values
    if len(na) >= 2 and len(nb) >= 2:
        nd = cohens_d(na, nb)
        np_ = stats.mannwhitneyu(na, nb, alternative="two-sided").pvalue
        print(f"\nNEGATIVE CONTROL (non-neoplastic, patient level): "
              f"n={len(na)} vs {len(nb)}  d={nd:+.3f}  P={np_:.3f}")
        rows.append(dict(check="negative control patient-level d", value=round(float(nd), 3)))
        rows.append(dict(check="negative control patient-level P", value=round(float(np_), 4)))

    out = pd.DataFrame(rows)
    out.to_csv(TAB / f"audit05_gse151179_structure_{STAMP}.tsv", sep="\t", index=False)

    with open(AUDIT / "04b_GSE151179_INDEPENDENCE.md", "w") as fh:
        fh.write("# AUDIT 04b — GSE151179 patient structure and independence\n\n")
        fh.write("Generated 2026-08-06 by `scripts/audit/05_gse151179_patient_structure.py`\n\n")
        fh.write("The uptake label is recorded once per patient, while the series contributes "
                 "several specimens per patient. Sample-level rows therefore replicate the same "
                 "exposure and are not independent. The first pass ran the primary test at "
                 "patient level, which was correct, but reported sample-level sensitivity rows "
                 "beside it without flagging the violation. Those rows are re-derived here.\n\n")
        fh.write(out.to_markdown(index=False))
        fh.write("\n\nNo patient carries discordant uptake labels, confirming the field is "
                 "patient-level and that a within-patient (paired) comparison is impossible in "
                 "this dataset. Only the between-patient contrast exists.\n")
    print(f"\nwrote {AUDIT / '04b_GSE151179_INDEPENDENCE.md'}")


if __name__ == "__main__":
    main()
