#!/usr/bin/env python3
"""Pooled test: does driver mutation class identify radioiodine refractoriness?

Three independent cohorts, all obtained today from open supplementary material with no
data-access application:

  Siraj 2022   KFSHRC, Saudi Arabia   n = 158   Cancers 14:1584, Table S1 + S3
  Zhang 2026   Fudan, China           n = 113   Cell Rep Med 7:102661, Table S1 + S2
  Boucai 2023  MSKCC, USA             n =  24   Clin Cancer Res 29:1620, Table S2

Endpoint differs slightly by cohort and that is stated rather than hidden: Siraj uses a
seven-criterion refractoriness definition, Zhang uses RAI-avid versus refractory on
post-therapy scans, Boucai uses RECIST v1.1 exceptional response. All three are
"did radioiodine work", which is the question.

The test is deliberately narrow: BRAF V600E versus RAS hotspot versus BRAF/RAS-negative,
against refractory yes/no. A random-effects meta-analysis pools the log odds ratio for the
BRAF/RAS-negative compartment, which is the compartment the DM1/DM2 manuscript is about.

Outputs:
  results/tables/pooled_driver_vs_rai_2026_08_06.tsv
  results/figures/figure_pooled_driver_vs_rai_2026_08_06.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
BOUCAI = Path("/data/rai_atlas/external/boucai2023/boucai_table_S2.xlsx")
STAMP = "2026_08_06"
RAS = {"NRAS", "HRAS", "KRAS"}


def boucai_counts():
    """Per-patient driver class and response from the ER/NR mutation sheet."""
    d = pd.read_excel(BOUCAI, "Somatic mutations in ER and NR", header=1)
    cols = list(d.columns)
    er = d.iloc[:, [0, 1, 2]].copy()
    er.columns = ["sample", "gene", "mutation"]
    er["group"] = "Exceptional responder"
    nr = d.iloc[:, [7, 8, 9]].copy()
    nr.columns = ["sample", "gene", "mutation"]
    nr["group"] = "Non-responder"
    m = pd.concat([er, nr], ignore_index=True).dropna(subset=["sample", "gene"])
    m = m[~m["sample"].astype(str).str.contains("Sample ID", na=False)]
    m["gene"] = m["gene"].astype(str).str.upper().str.strip()
    m["mutation"] = m["mutation"].astype(str).str.upper().str.strip()

    rows = []
    for (s, g), sub in m.groupby(["sample", "group"]):
        genes = set(sub["gene"])
        muts = " ".join(sub["mutation"])
        is_braf = any("BRAF" in x for x in genes) and "V600E" in muts
        is_ras = bool(genes & RAS) and any(k in muts for k in ("Q61", "G12", "G13"))
        cls = "BRAF V600E" if is_braf else ("RAS hotspot" if is_ras else "BRAF/RAS-negative")
        rows.append(dict(cohort="Boucai 2023", patient=s, group=g, driver=cls,
                         refractory=int(g == "Non-responder")))
    return pd.DataFrame(rows)


def load_cohort(fname, cohort, driver_col, resp_col, refractory_value):
    p = TAB / fname
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    d = d[d[resp_col].notna() & d[driver_col].notna()]
    out = pd.DataFrame({
        "cohort": cohort,
        "patient": d.iloc[:, 0].astype(str),
        "driver": d[driver_col],
        "refractory": (d[resp_col].astype(str).str.strip() == refractory_value).astype(int),
    })
    return out


def main():
    parts = [
        load_cohort(f"siraj2022_patient_level_{STAMP}.tsv", "Siraj 2022", "driver", "rai", "Refractory"),
        load_cohort(f"zhang2026_patient_table_{STAMP}.tsv", "Zhang 2026", "driver",
                    "rai_sensitivity", "Refractory"),
        boucai_counts(),
    ]
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    print("pooled cohorts:", df.groupby("cohort").size().to_dict())
    print("total patients:", len(df), "| refractory:", int(df.refractory.sum()))

    rows = []
    for coh, sub in df.groupby("cohort"):
        ct = pd.crosstab(sub["driver"], sub["refractory"])
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            continue
        chi = stats.chi2_contingency(ct.values)
        rows.append(dict(cohort=coh, test="driver class x refractory (chi-square)",
                         n=len(sub), statistic=float(chi.statistic), p=float(chi.pvalue),
                         detail=str(ct.to_dict())))
        print(f"\n{coh} (n={len(sub)}): chi-square P = {chi.pvalue:.3g}")
        for d_ in ct.index:
            a = int(ct.loc[d_, 1]); b = int(ct.loc[d_, 0])
            print(f"   {d_:20s} {a}/{a+b} refractory ({100*a/max(a+b,1):.0f}%)")

    # random-effects meta-analysis of the BRAF/RAS-negative compartment
    meta = []
    for coh, sub in df.groupby("cohort"):
        neg = sub["driver"] == "BRAF/RAS-negative"
        a = int(((neg) & (sub.refractory == 1)).sum()); b = int(((neg) & (sub.refractory == 0)).sum())
        c = int(((~neg) & (sub.refractory == 1)).sum()); d_ = int(((~neg) & (sub.refractory == 0)).sum())
        # Haldane-Anscombe correction for any zero cell
        aa, bb, cc, dd = (a + .5, b + .5, c + .5, d_ + .5) if 0 in (a, b, c, d_) else (a, b, c, d_)
        lor = np.log((aa * dd) / (bb * cc))
        se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
        meta.append(dict(cohort=coh, n=len(sub), a=a, b=b, c=c, d=d_, lor=lor, se=se,
                         orr=np.exp(lor), lo=np.exp(lor - 1.96 * se), hi=np.exp(lor + 1.96 * se)))
    md = pd.DataFrame(meta)

    w = 1 / md.se ** 2
    fe = float((w * md.lor).sum() / w.sum())
    q = float((w * (md.lor - fe) ** 2).sum())
    dfree = len(md) - 1
    tau2 = max(0.0, (q - dfree) / (w.sum() - (w ** 2).sum() / w.sum())) if dfree > 0 else 0.0
    wr = 1 / (md.se ** 2 + tau2)
    re = float((wr * md.lor).sum() / wr.sum())
    se_re = float(np.sqrt(1 / wr.sum()))
    z = re / se_re
    p_re = 2 * (1 - stats.norm.cdf(abs(z)))
    i2 = max(0.0, (q - dfree) / q * 100) if q > 0 else 0.0
    print(f"\nrandom-effects OR (BRAF/RAS-negative vs rest) = {np.exp(re):.2f} "
          f"[{np.exp(re-1.96*se_re):.2f}, {np.exp(re+1.96*se_re):.2f}], P = {p_re:.3g}, "
          f"I2 = {i2:.0f}%, tau2 = {tau2:.3f}")

    rows.append(dict(cohort="POOLED (random effects)",
                     test="BRAF/RAS-negative vs rest — refractory",
                     n=int(md.n.sum()), statistic=float(np.exp(re)), p=float(p_re),
                     detail=f"95% CI {np.exp(re-1.96*se_re):.2f}-{np.exp(re+1.96*se_re):.2f}; "
                            f"I2={i2:.0f}%"))
    pd.DataFrame(rows).to_csv(TAB / f"pooled_driver_vs_rai_{STAMP}.tsv", sep="\t", index=False)
    df.to_csv(TAB / f"pooled_driver_patient_level_{STAMP}.tsv", sep="\t", index=False)

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.0), facecolor="white", layout="constrained",
                             gridspec_kw={"width_ratios": [1.25, 1.0]})

    ax = axes[0]
    order = ["BRAF V600E", "RAS hotspot", "BRAF/RAS-negative"]
    cohorts = list(df.cohort.unique())
    colors = {"Siraj 2022": "#37618e", "Zhang 2026": "#9c4742", "Boucai 2023": "#b5761f"}
    width = 0.26
    for i, coh in enumerate(cohorts):
        sub = df[df.cohort == coh]
        fr, ns = [], []
        for o in order:
            g = sub[sub.driver == o]
            fr.append(100 * g.refractory.mean() if len(g) else np.nan)
            ns.append(len(g))
        xs = np.arange(len(order)) + (i - (len(cohorts) - 1) / 2) * width
        ax.bar(xs, fr, width=width * 0.92, color=colors.get(coh, "#666"), alpha=0.88,
               edgecolor="white", label=f"{coh} (n={len(sub)})")
        for x_, f_, n_ in zip(xs, fr, ns):
            if np.isfinite(f_):
                ax.text(x_, f_ + 1.5, f"{n_}", ha="center", fontsize=7.5, color="#555")
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([o.replace(" ", "\n") for o in order], fontsize=9.5)
    ax.set_ylabel("% radioiodine-refractory")
    ax.legend(fontsize=8.5, frameon=False, loc="upper left")
    per = {r["cohort"]: r["p"] for r in rows if "chi-square" in r["test"]}
    lab = " · ".join(f"{c.split()[0]} P={per[c]:.3g}" for c in cohorts if c in per)
    ax.set_title("a · Refractoriness by driver class in three independent cohorts\n" + lab,
                 fontsize=10.5, loc="left", fontweight="bold")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    ax = axes[1]
    y = np.arange(len(md))
    ax.errorbar(md.orr, y, xerr=[md.orr - md.lo, md.hi - md.orr], fmt="s",
                color="#37618e", ecolor="#8fa8c4", capsize=3, ms=7)
    ax.errorbar([np.exp(re)], [len(md)], xerr=[[np.exp(re) - np.exp(re - 1.96 * se_re)],
                                               [np.exp(re + 1.96 * se_re) - np.exp(re)]],
                fmt="D", color="#9c2b25", ecolor="#9c2b25", capsize=4, ms=10)
    ax.axvline(1, color="#666", lw=0.8, ls="--")
    ax.set_yticks(list(y) + [len(md)])
    ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(md.cohort, md.n)]
                       + [f"Pooled (n={int(md.n.sum())})"], fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio for refractoriness, BRAF/RAS-negative vs rest")
    ax.set_title(f"b · Random-effects meta-analysis\nOR = {np.exp(re):.2f} "
                 f"[{np.exp(re-1.96*se_re):.2f}, {np.exp(re+1.96*se_re):.2f}], "
                 f"P = {p_re:.2f}, I² = {i2:.0f}%",
                 fontsize=10.5, loc="left", fontweight="bold")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)

    fig.suptitle("Driver mutation class does not identify radioiodine-refractory disease — "
                 "three cohorts, %d patients" % len(df), fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / f"figure_pooled_driver_vs_rai_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_pooled_driver_vs_rai_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"\nwrote {FIG / f'figure_pooled_driver_vs_rai_{STAMP}.png'}")


if __name__ == "__main__":
    main()
